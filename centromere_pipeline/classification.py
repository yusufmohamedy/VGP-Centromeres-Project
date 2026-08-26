"""
Candidate classification and YAML serialization module.

Implements:
1. Copy-number and chromosome position calculation
2. Noise pre-filtering and quality gating
3. Global dominant unit detection
4. Primary and secondary centromere candidate selection with identity tiebreaking
5. Uncertainty flag assignment (short, low_score, other_region, other_unit)
6. Chromosome shape classification (metacentric, submetacentric, acrocentric, telocentric)
7. Structured YAML summary serialization (including consensus FASTA sequences)
"""

from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import polars as pl

from .config import PipelineConfig


def classify_shape(midpoint: float, chrom_length: float) -> str:
    """Classify chromosome morphology based on centromere relative position."""
    if chrom_length <= 0:
        return "metacentric"
    pos = midpoint / float(chrom_length)
    if pos < 0.05 or pos > 0.95:
        return "telocentric"
    elif pos < 0.25 or pos > 0.75:
        return "acrocentric"
    elif (0.25 <= pos <= 0.40) or (0.60 <= pos <= 0.75):
        return "submetacentric"
    elif 0.40 < pos < 0.60:
        return "metacentric"
    return ""


def classify_clusters(
    raw_clusters_df: pl.DataFrame,
    chrom_lengths: Dict[str, int],
    config: PipelineConfig,
    is_scaffolds: bool = False,
    has_identity: bool = False,
) -> pl.DataFrame:
    """Classify raw Kadane clusters for chromosomes or scaffolds.

    For scaffolds: outputs simplified table sorted by total_chunk_bp descending.
    For chromosomes: executes 9-phase classification logic with primary/secondary
    calls, uncertainty flags, shape assignment, and sort by total_chunk_bp descending.
    """
    if len(raw_clusters_df) == 0:
        if is_scaffolds:
            cols = [
                ("scaffold", pl.Utf8),
                ("cluster_start", pl.Int64),
                ("cluster_end", pl.Int64),
                ("unit_id", pl.Utf8),
                ("density_score", pl.Float64),
            ]
            if has_identity:
                cols.append(("identity_score", pl.Float64))
            cols.extend([
                ("unit_size", pl.Int64),
                ("cluster_size", pl.Float64),
                ("total_chunk_bp", pl.Float64),
                ("copy_number", pl.Float64),
            ])
            return pl.DataFrame(schema=dict(cols))
        else:
            cols = [
                ("chromosome", pl.Utf8),
                ("cluster_start", pl.Int64),
                ("cluster_end", pl.Int64),
                ("unit_id", pl.Utf8),
                ("density_score", pl.Float64),
            ]
            if has_identity:
                cols.append(("identity_score", pl.Float64))
            cols.extend([
                ("unit_size", pl.Int64),
                ("cluster_size", pl.Float64),
                ("total_chunk_bp", pl.Float64),
                ("copy_number", pl.Float64),
                ("chromosome_position", pl.Float64),
                ("classification", pl.Utf8),
                ("uncertainty", pl.Utf8),
                ("shape", pl.Utf8),
            ])
            return pl.DataFrame(schema=dict(cols))

    # Basic computed columns
    df = raw_clusters_df.with_columns([
        pl.col("chromosome").cast(pl.Utf8),
        pl.col("unit_id").cast(pl.Utf8),
        pl.col("unit_size").cast(pl.Int64),
        pl.col("cluster_start").cast(pl.Int64),
        pl.col("cluster_end").cast(pl.Int64),
        (pl.col("cluster_end") - pl.col("cluster_start")).cast(pl.Float64).alias("cluster_size"),
        pl.col("total_chunk_bp").cast(pl.Float64),
        pl.col("density_score").cast(pl.Float64),
        (pl.col("total_chunk_bp") / pl.col("unit_size")).cast(pl.Float64).alias("copy_number"),
    ])

    if has_identity and "identity_score" in df.columns:
        df = df.with_columns(pl.col("identity_score").cast(pl.Float64))

    # For scaffolds: Return immediately sorted descending by total_chunk_bp
    if is_scaffolds:
        out_cols = ["scaffold", "cluster_start", "cluster_end", "unit_id", "density_score"]
        if has_identity and "identity_score" in df.columns:
            out_cols.append("identity_score")
        out_cols.extend(["unit_size", "cluster_size", "total_chunk_bp", "copy_number"])

        scaf_df = (
            df.rename({"chromosome": "scaffold"})
            .select(out_cols)
            .sort("total_chunk_bp", descending=True)
        )
        return scaf_df

    # Add chromosome length and position
    # Map chrom_length using Python dict
    chrom_lens_col = [chrom_lengths.get(c, 0) for c in df["chromosome"].to_list()]
    df = df.with_columns([
        pl.Series("chrom_length", chrom_lens_col, dtype=pl.Int64),
        pl.int_range(0, len(df), dtype=pl.Int64).alias("row_id"),
    ])

    midpoints = (df["cluster_start"] + df["cluster_end"]) / 2.0
    lens_arr = df["chrom_length"].to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        pos_raw = np.where(lens_arr > 0, midpoints.to_numpy() / lens_arr, 0.0)
        positions = np.round(pos_raw, 4)

    df = df.with_columns(pl.Series("chromosome_position", positions, dtype=pl.Float64))

    # Phase 1: Noise Pre-Filter & Quality Gate
    # Disqualify clusters below min_copy_multiplier * unit_size
    valid_clusters_df = df.filter(
        pl.col("total_chunk_bp") >= (config.min_copy_multiplier * pl.col("unit_size"))
    )

    quality_df = valid_clusters_df.filter(
        pl.col("density_score") >= config.min_density
    )
    if len(quality_df) == 0 and len(valid_clusters_df) > 0:
        quality_df = valid_clusters_df

    primary_eligible_df = quality_df

    # Phase 2: Find biggest cluster per chromosome (among eligible)
    biggest_per_chr: List[Dict[str, Any]] = []
    if len(primary_eligible_df) > 0:
        for chr_key, group in primary_eligible_df.partition_by("chromosome", as_dict=True).items():
            chr_name = chr_key[0] if isinstance(chr_key, tuple) else str(chr_key)
            top = group.sort("total_chunk_bp", descending=True).row(0, named=True)
            biggest_per_chr.append({
                "chromosome": chr_name,
                "biggest_unit_id": top["unit_id"],
                "biggest_bp": top["total_chunk_bp"],
                "biggest_row_id": top["row_id"],
            })

    # Phase 3: Find best potential unit (unit_id biggest on most chromosomes)
    unit_stats: Dict[str, Dict[str, Any]] = defaultdict(lambda: {"n_chrs": 0, "tot_bp": 0.0})
    for b in biggest_per_chr:
        u = b["biggest_unit_id"]
        unit_stats[u]["n_chrs"] += 1
        unit_stats[u]["tot_bp"] += b["biggest_bp"]

    sorted_units = sorted(
        unit_stats.items(), key=lambda x: (x[1]["n_chrs"], x[1]["tot_bp"]), reverse=True
    )
    best_potential_unit = sorted_units[0][0] if sorted_units else ""

    # Phase 4 & 5: Calculate avg primary total_chunk_bp from Group A
    group_a_biggest = [b for b in biggest_per_chr if b["biggest_unit_id"] == best_potential_unit]
    if group_a_biggest:
        avg_primary_bp = sum(b["biggest_bp"] for b in group_a_biggest) / float(len(group_a_biggest))
    elif biggest_per_chr:
        avg_primary_bp = sum(b["biggest_bp"] for b in biggest_per_chr) / float(len(biggest_per_chr))
    else:
        avg_primary_bp = 0.0

    # Phase 6: Primary & Secondary Candidate Selection
    primary_row_ids: set[int] = set()
    secondary_row_ids: set[int] = set()

    all_chrs = df["chromosome"].unique(maintain_order=True).to_list()

    for chr_name in all_chrs:
        chr_clusters = primary_eligible_df.filter(pl.col("chromosome") == chr_name).sort(
            "total_chunk_bp", descending=True
        )
        if len(chr_clusters) == 0:
            chr_clusters = quality_df.filter(pl.col("chromosome") == chr_name).sort(
                "total_chunk_bp", descending=True
            )
        if len(chr_clusters) == 0:
            continue

        biggest_cluster = chr_clusters.row(0, named=True)

        # Same-unit tiebreaker using identity score (if available)
        if has_identity and "identity_score" in chr_clusters.columns:
            same_unit = chr_clusters.filter(pl.col("unit_id") == biggest_cluster["unit_id"])
            if len(same_unit) > 1:
                biggest_same_bp = same_unit.row(0, named=True)["total_chunk_bp"]
                similar = same_unit.filter(
                    pl.col("total_chunk_bp") >= (config.identity_tiebreak_pct * biggest_same_bp)
                )
                if len(similar) > 1:
                    biggest_cluster = similar.sort(
                        by=["identity_score", "total_chunk_bp"],
                        descending=[True, True],
                        nulls_last=True,
                    ).row(0, named=True)

        # Check if best_potential_unit exists on this chromosome
        gwc_clusters = chr_clusters.filter(pl.col("unit_id") == best_potential_unit).sort(
            "total_chunk_bp", descending=True
        )

        selected_primary = biggest_cluster
        if len(gwc_clusters) > 0:
            if has_identity and "identity_score" in gwc_clusters.columns:
                biggest_gwc_bp = gwc_clusters.row(0, named=True)["total_chunk_bp"]
                similar_gwc = gwc_clusters.filter(
                    pl.col("total_chunk_bp") >= (config.identity_tiebreak_pct * biggest_gwc_bp)
                )
                if len(similar_gwc) > 1:
                    top_gwc = similar_gwc.sort(
                        by=["identity_score", "total_chunk_bp"],
                        descending=[True, True],
                        nulls_last=True,
                    ).row(0, named=True)
                else:
                    top_gwc = gwc_clusters.row(0, named=True)
            else:
                top_gwc = gwc_clusters.row(0, named=True)

            if top_gwc["total_chunk_bp"] >= (config.min_size_pct * avg_primary_bp):
                selected_primary = top_gwc

        sel_row_id = int(selected_primary["row_id"])
        sel_total_chunk_bp = float(selected_primary["total_chunk_bp"])
        primary_row_ids.add(sel_row_id)

        # Secondary candidates:
        chr_quality = quality_df.filter(
            (pl.col("chromosome") == chr_name)
            & (pl.col("row_id") != sel_row_id)
            & (pl.col("total_chunk_bp") >= (config.other_region_pct * sel_total_chunk_bp))
        )
        for r_id in chr_quality["row_id"].to_list():
            secondary_row_ids.add(int(r_id))

    # Phase 7: Assign Classifications
    classifications = []
    for r_id in df["row_id"].to_list():
        if r_id in primary_row_ids:
            classifications.append("primary")
        elif r_id in secondary_row_ids:
            classifications.append("secondary")
        else:
            classifications.append("")

    df = df.with_columns(pl.Series("classification", classifications, dtype=pl.Utf8))

    # Phase 8: Uncertainty Flags (Primary candidates only)
    other_region_row_ids: set[int] = set()
    other_unit_row_ids: set[int] = set()

    for row in df.iter_rows(named=True):
        if row["classification"] == "primary":
            r_id = row["row_id"]
            chr_name = row["chromosome"]
            u_id = row["unit_id"]
            tot_bp = row["total_chunk_bp"]

            other_clusters = df.filter(
                (pl.col("chromosome") == chr_name)
                & (pl.col("row_id") != r_id)
                & (pl.col("density_score") >= config.min_density)
            )

            # Same-unit competition
            same_unit_comp = other_clusters.filter(
                (pl.col("unit_id") == u_id)
                & (pl.col("total_chunk_bp") >= (config.other_region_pct * tot_bp))
            )
            if len(same_unit_comp) > 0:
                other_region_row_ids.add(r_id)

            # Different-unit competition
            diff_unit_comp = other_clusters.filter(
                (pl.col("unit_id") != u_id)
                & (pl.col("total_chunk_bp") >= (config.other_region_pct * tot_bp))
            )
            if len(diff_unit_comp) > 0:
                other_unit_row_ids.add(r_id)

    uncertainties = []
    for row in df.iter_rows(named=True):
        if row["classification"] == "primary":
            r_id = row["row_id"]
            flags = []
            if row["total_chunk_bp"] < (config.min_size_pct * avg_primary_bp):
                flags.append("short")
            if row["density_score"] < config.low_score_primary:
                flags.append("low_score")
            if r_id in other_region_row_ids:
                flags.append("other_region")
            if r_id in other_unit_row_ids:
                flags.append("other_unit")
            uncertainties.append(";".join(flags))
        else:
            uncertainties.append("")

    df = df.with_columns(pl.Series("uncertainty", uncertainties, dtype=pl.Utf8))

    # Phase 9: Shape Classification
    shapes = []
    for row in df.iter_rows(named=True):
        if not row["classification"]:
            shapes.append("")
        else:
            pos = row["chromosome_position"]
            if pos < 0.05 or pos > 0.95:
                shapes.append("telocentric")
            elif pos < 0.25 or pos > 0.75:
                shapes.append("acrocentric")
            elif (0.25 <= pos <= 0.40) or (0.60 <= pos <= 0.75):
                shapes.append("submetacentric")
            elif 0.40 < pos < 0.60:
                shapes.append("metacentric")
            else:
                shapes.append("")

    df = df.with_columns(pl.Series("shape", shapes, dtype=pl.Utf8))

    # Select final columns and sort by total_chunk_bp descending
    out_cols = ["chromosome", "cluster_start", "cluster_end", "unit_id", "density_score"]
    if has_identity and "identity_score" in df.columns:
        out_cols.append("identity_score")
    out_cols.extend([
        "unit_size",
        "cluster_size",
        "total_chunk_bp",
        "copy_number",
        "chromosome_position",
        "classification",
        "uncertainty",
        "shape",
    ])

    final_df = df.select(out_cols).sort("total_chunk_bp", descending=True)
    return final_df


def dump_yaml_clean(data: Any, indent: int = 0, is_list_with_spacing: bool = False) -> str:
    """Custom YAML formatter matching original pipeline output with clean indentation and spacing."""
    lines: List[str] = []
    ind = "  " * indent

    if isinstance(data, dict):
        for k, v in data.items():
            if v is None:
                continue
            if k in ("best_candidate", "alternative_candidates", "chromosomes", "alternate_centromeres") and lines and lines[-1] != "":
                lines.append("")
            if isinstance(v, (dict, list)):
                lines.append(f"{ind}{k}:")
                is_spaced_list = k in ("chromosomes", "alternative_candidates", "alternate_centromeres")
                lines.append(dump_yaml_clean(v, indent + 1, is_list_with_spacing=is_spaced_list))
            else:
                if isinstance(v, str) and (":" in v or "\n" in v or v.startswith("<")):
                    val_str = f'"{v}"'
                else:
                    val_str = str(v)
                lines.append(f"{ind}{k}: {val_str}")
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            if is_list_with_spacing and idx > 0:
                lines.append("")  # Blank line between entries for visual readability

            if isinstance(item, dict):
                first = True
                for k, v in item.items():
                    if v is None:
                        continue
                    if k in ("alternate_centromeres",) and lines and lines[-1] != "":
                        lines.append("")
                    is_spaced_list = k in ("chromosomes", "alternative_candidates", "alternate_centromeres")
                    if isinstance(v, str) and (":" in v or "\n" in v or v.startswith("<")):
                        val_str = f'"{v}"'
                    else:
                        val_str = str(v)
                    if first:
                        if isinstance(v, (dict, list)):
                            lines.append(f"{ind}- {k}:")
                            lines.append(dump_yaml_clean(v, indent + 2, is_list_with_spacing=is_spaced_list))
                        else:
                            lines.append(f"{ind}- {k}: {val_str}")
                        first = False
                    else:
                        if isinstance(v, (dict, list)):
                            lines.append(f"{ind}  {k}:")
                            lines.append(dump_yaml_clean(v, indent + 2, is_list_with_spacing=is_spaced_list))
                        else:
                            lines.append(f"{ind}  {k}: {val_str}")
            else:
                lines.append(f"{ind}- {item}")

    return "\n".join(lines)


def generate_centromere_yaml(
    chr_clusters_df: pl.DataFrame,
    scaffold_clusters_df: pl.DataFrame,
    chrom_lengths: Dict[str, int],
    scaffold_lengths: Dict[str, int],
    species_name: str,
    fasta_sequences: Dict[str, str],
    output_path: Path,
    has_identity: bool = False,
) -> None:
    """Generate structured centromere summary YAML file matching original specifications."""
    chr_clusters = chr_clusters_df.to_dicts() if len(chr_clusters_df) > 0 else []
    scaffold_clusters = scaffold_clusters_df.to_dicts() if len(scaffold_clusters_df) > 0 else []

    # Estimate missing chrom lengths if necessary
    for cl in chr_clusters:
        c = cl["chromosome"]
        pos = cl.get("chromosome_position", 0.0)
        if c not in chrom_lengths and pos and pos > 0:
            mid = (cl["cluster_start"] + cl["cluster_end"]) / 2.0
            chrom_lengths[c] = int(round(mid / pos))

    # Group clusters by chromosome
    clusters_by_chr: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for cl in chr_clusters:
        clusters_by_chr[cl["chromosome"]].append(cl)

    # Unit-level repeat totals and unit lengths
    unit_chr_bp: Dict[str, int] = defaultdict(int)
    unit_scaf_bp: Dict[str, int] = defaultdict(int)
    unit_lengths: Dict[str, int] = {}

    for cl in chr_clusters:
        u = str(cl["unit_id"])
        unit_chr_bp[u] += int(cl["total_chunk_bp"])
        if u not in unit_lengths and int(cl.get("unit_size", 0)) > 0:
            unit_lengths[u] = int(cl["unit_size"])

    for cl in scaffold_clusters:
        u = str(cl["unit_id"])
        unit_scaf_bp[u] += int(cl["total_chunk_bp"])
        if u not in unit_lengths and int(cl.get("unit_size", 0)) > 0:
            unit_lengths[u] = int(cl["unit_size"])

    all_detected_units = set(unit_chr_bp.keys()).union(set(unit_scaf_bp.keys()))
    for u in all_detected_units:
        if u not in unit_lengths:
            parts = u.split("-")
            try:
                unit_lengths[u] = int(parts[-1])
            except ValueError:
                unit_lengths[u] = 100

    unit_pct_scaf: Dict[str, float] = {}
    for u in all_detected_units:
        tot_unit = unit_chr_bp[u] + unit_scaf_bp[u]
        unit_pct_scaf[u] = round((unit_scaf_bp[u] / float(tot_unit)) * 100.0, 5) if tot_unit > 0 else 0.0

    primary_clusters_by_unit: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    primary_chrs_by_unit: Dict[str, set[str]] = defaultdict(set)
    present_chrs_by_unit: Dict[str, set[str]] = defaultdict(set)

    yaml_chromosomes: List[Dict[str, Any]] = []

    for chr_id, chrom_len in chrom_lengths.items():
        chr_list = clusters_by_chr.get(chr_id, [])

        primary_cluster = next((cl for cl in chr_list if cl.get("classification") == "primary"), None)

        shape: Optional[str] = None
        centromere_dict: Optional[Dict[str, Any]] = None

        if primary_cluster:
            sp_start = int(primary_cluster["cluster_start"])
            sp_end = int(primary_cluster["cluster_end"])
            sp_score = int(primary_cluster["density_score"])
            sp_ident = primary_cluster.get("identity_score")

            u = str(primary_cluster["unit_id"])
            primary_clusters_by_unit[u].append(primary_cluster)
            primary_chrs_by_unit[u].add(chr_id)
            present_chrs_by_unit[u].add(chr_id)

            shape = primary_cluster.get("shape")
            if not shape:
                midpoint = (sp_start + sp_end) / 2.0
                shape = classify_shape(midpoint, chrom_len)

            # Uncertainty flags
            raw_uncertainty = primary_cluster.get("uncertainty", "")
            centromere_uncertain = None
            if raw_uncertainty:
                formatted_flags = [flag.replace("_", " ") for flag in raw_uncertainty.split(";") if flag.strip()]
                if formatted_flags:
                    centromere_uncertain = ", ".join(formatted_flags)

            centromere_dict = {
                "start": sp_start,
                "end": sp_end,
                "unit_name": primary_cluster["unit_id"],
                "unit_length": int(primary_cluster["unit_size"]),
                "density_score": sp_score,
            }
            if has_identity and sp_ident is not None:
                centromere_dict["identity_score"] = int(round(float(sp_ident)))

            if shape:
                centromere_dict["centromere_position"] = shape

            if centromere_uncertain:
                centromere_dict["centromere_uncertain"] = centromere_uncertain

        # Secondary clusters (alternate centromeres)
        secondary_clusters = [cl for cl in chr_list if cl.get("classification") == "secondary"]
        secondary_clusters.sort(key=lambda x: x["total_chunk_bp"], reverse=True)
        alternates: List[Dict[str, Any]] = []

        for cl in secondary_clusters:
            u_sec = str(cl["unit_id"])
            present_chrs_by_unit[u_sec].add(chr_id)

            alt_shape = cl.get("shape")
            if not alt_shape:
                alt_midpoint = (cl["cluster_start"] + cl["cluster_end"]) / 2.0
                alt_shape = classify_shape(alt_midpoint, chrom_len)

            alt_entry = {
                "start": int(cl["cluster_start"]),
                "end": int(cl["cluster_end"]),
                "unit_name": cl["unit_id"],
                "unit_length": int(cl["unit_size"]),
                "density_score": int(cl["density_score"]),
            }
            if has_identity and cl.get("identity_score") is not None:
                alt_entry["identity_score"] = int(round(float(cl["identity_score"])))
            if alt_shape:
                alt_entry["centromere_position"] = alt_shape

            alternates.append(alt_entry)

        chr_entry: Dict[str, Any] = {
            "accession": chr_id,
            "chromosome_length": int(chrom_len),
        }
        if centromere_dict:
            chr_entry["centromere"] = centromere_dict
        if alternates:
            chr_entry["alternate_centromeres"] = alternates

        yaml_chromosomes.append(chr_entry)

    # Candidate summaries: best_candidate & alternative_candidates
    candidate_units = [u for u in all_detected_units if len(present_chrs_by_unit[u]) > 0]

    best_unit: Optional[str] = None
    if candidate_units:
        units_with_primary = [u for u in candidate_units if len(primary_chrs_by_unit[u]) > 0]
        if units_with_primary:
            best_unit = max(
                units_with_primary,
                key=lambda u: (
                    len(primary_chrs_by_unit[u]),
                    sum(cl["total_chunk_bp"] for cl in primary_clusters_by_unit[u]),
                    len(present_chrs_by_unit[u]),
                    unit_chr_bp[u] + unit_scaf_bp[u],
                ),
            )
        else:
            best_unit = max(
                candidate_units,
                key=lambda u: (
                    len(present_chrs_by_unit[u]),
                    unit_chr_bp[u] + unit_scaf_bp[u],
                ),
            )

    best_candidate_dict: Optional[Dict[str, Any]] = None
    if best_unit:
        best_candidate_dict = {
            "unit_name": best_unit,
            "unit_length": unit_lengths.get(best_unit, 0),
        }
        if best_unit in fasta_sequences and fasta_sequences[best_unit]:
            best_candidate_dict["consensus_sequence"] = fasta_sequences[best_unit]

        pri_list = primary_clusters_by_unit[best_unit]
        if pri_list:
            avg_len = int(round(sum((cl["cluster_end"] - cl["cluster_start"]) for cl in pri_list) / float(len(pri_list))))
            avg_dens = int(round(sum(cl["density_score"] for cl in pri_list) / float(len(pri_list))))
            best_candidate_dict["average_centromere_candidate_length"] = avg_len
            best_candidate_dict["average_centromere_candidate_density"] = avg_dens

            if has_identity:
                ident_scores = [cl["identity_score"] for cl in pri_list if cl.get("identity_score") is not None]
                if ident_scores:
                    best_candidate_dict["average_centromere_candidate_identity"] = int(round(sum(ident_scores) / float(len(ident_scores))))

        best_candidate_dict["number_of_chromosomes_where_it_is_a_centromeric_candidate"] = len(primary_chrs_by_unit[best_unit])
        best_candidate_dict["number_of_chromosomes_where_it_is_present"] = len(present_chrs_by_unit[best_unit])
        best_candidate_dict["percentage_of_best_candidate_in_scaffolds"] = unit_pct_scaf.get(best_unit, 0.0)

    alternative_candidates_list: List[Dict[str, Any]] = []
    alt_units = [u for u in candidate_units if u != best_unit]
    alt_units.sort(
        key=lambda u: (
            len(primary_chrs_by_unit[u]),
            len(present_chrs_by_unit[u]),
            unit_chr_bp[u] + unit_scaf_bp[u],
        ),
        reverse=True,
    )

    for u in alt_units:
        alt_entry = {
            "unit_name": u,
            "unit_length": unit_lengths.get(u, 0),
        }
        if u in fasta_sequences and fasta_sequences[u]:
            alt_entry["consensus_sequence"] = fasta_sequences[u]

        pri_list = primary_clusters_by_unit[u]
        if pri_list:
            avg_len = int(round(sum((cl["cluster_end"] - cl["cluster_start"]) for cl in pri_list) / float(len(pri_list))))
            avg_dens = int(round(sum(cl["density_score"] for cl in pri_list) / float(len(pri_list))))
            alt_entry["average_centromere_candidate_length"] = avg_len
            alt_entry["average_centromere_candidate_density"] = avg_dens

            if has_identity:
                ident_scores = [cl["identity_score"] for cl in pri_list if cl.get("identity_score") is not None]
                if ident_scores:
                    alt_entry["average_centromere_candidate_identity"] = int(round(sum(ident_scores) / float(len(ident_scores))))

        alt_entry["number_of_chromosomes_where_it_is_a_centromeric_candidate"] = len(primary_chrs_by_unit[u])
        alt_entry["number_of_chromosomes_where_it_is_present"] = len(present_chrs_by_unit[u])
        alt_entry["percentage_of_alt_candidate_in_scaffolds"] = unit_pct_scaf.get(u, 0.0)

        alternative_candidates_list.append(alt_entry)

    # Genome-wide summary stats
    total_chrom_length = sum(chrom_lengths.values())
    total_scaffold_length = sum(scaffold_lengths.values())
    genome_size = total_chrom_length + total_scaffold_length

    total_chr_repeats = sum(cl["total_chunk_bp"] for cl in chr_clusters)
    total_scaffold_repeats = sum(cl["total_chunk_bp"] for cl in scaffold_clusters)
    total_genome_repeats = total_chr_repeats + total_scaffold_repeats

    pct_genome_in_scaffolds = round((total_scaffold_length / float(genome_size)) * 100.0, 5) if genome_size > 0 else 0.0
    pct_all_repeats_in_scaffolds = round((total_scaffold_repeats / float(total_genome_repeats)) * 100.0, 5) if total_genome_repeats > 0 else 0.0

    # Assemble YAML structure
    species_dict: Dict[str, Any] = {
        "name": species_name,
        "total_chromosomes": len(chrom_lengths),
        "genome_size": genome_size,
        "total_scaffold_length": total_scaffold_length,
        "percentage_of_genome_in_scaffolds": pct_genome_in_scaffolds,
        "percentage_of_all_repeats_in_scaffolds": pct_all_repeats_in_scaffolds,
    }

    if best_candidate_dict:
        species_dict["best_candidate"] = best_candidate_dict

    if alternative_candidates_list:
        species_dict["alternative_candidates"] = alternative_candidates_list

    species_dict["chromosomes"] = yaml_chromosomes

    yaml_doc = {"species": species_dict}
    yaml_output = dump_yaml_clean(yaml_doc)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(yaml_output + "\n")
