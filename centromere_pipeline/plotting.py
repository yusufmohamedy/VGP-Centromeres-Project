"""
Karyotype visualization module for chromosomes and concatenated scaffolds.

Generates publication-ready visualizations:
1. Combined All-Clusters Karyotype Plot (matching plot_combined_all_clusters.R)
2. Combined YAML Summary Karyotype Plot (matching plot_combined_yaml.R)
"""

from collections import defaultdict
import math
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import polars as pl
import yaml

from .config import PipelineConfig


def get_chromosome_sort_key(name: str) -> Tuple[int, Any]:
    """Generate sort key for chromosome identifiers.

    Supports 'chr1', 'chr2', ..., 'chrX', 'chrY', 'chrM' and accession styles.
    """
    clean = re.sub(r"^chr", "", name, flags=re.IGNORECASE).strip()
    if clean.isdigit():
        return (0, int(clean))
    elif clean.upper() == "X":
        return (1, 1000)
    elif clean.upper() == "Y":
        return (1, 1001)
    elif clean.upper() in ("M", "MT"):
        return (1, 1002)
    else:
        # Natural alphanumeric split
        parts = [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", name)]
        return (2, parts)


def sort_chromosome_records(
    chrom_lengths: Dict[str, int]
) -> List[Tuple[str, int]]:
    """Sort chromosomes by biological/natural sorting convention."""
    names = list(chrom_lengths.keys())
    is_named_style = any(re.match(r"^chr", n, flags=re.IGNORECASE) for n in names)

    if is_named_style:
        sorted_names = sorted(names, key=get_chromosome_sort_key)
    else:
        sorted_names = sorted(
            names,
            key=lambda x: [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", x)],
        )

    return [(name, chrom_lengths[name]) for name in sorted_names]


def build_scaffold_offsets(
    scaffold_lengths: Dict[str, int]
) -> Tuple[List[Dict[str, Any]], int]:
    """Sort scaffolds descending by length and compute cumulative offsets.

    Returns:
        (scaffold_offset_records, total_scaffold_length)
    """
    if not scaffold_lengths:
        return [], 0

    sorted_scafs = sorted(scaffold_lengths.items(), key=lambda x: x[1], reverse=True)
    records = []
    current_offset = 0

    for name, length in sorted_scafs:
        records.append({
            "name": name,
            "length": length,
            "offset": current_offset,
            "end_offset": current_offset + length,
        })
        current_offset += length

    return records, current_offset


def get_color_palette(unit_ids: List[str]) -> Tuple[Dict[str, Any], float]:
    """Generate color mapping for distinct satellite family units matching RColorBrewer palettes.

    <=8 units: Set2
    <=12 units: Paired
    >12 units: Interpolated tab20
    """
    sorted_units = sorted(set(unit_ids))
    n_units = len(sorted_units)

    if n_units <= 5:
        alpha = 0.88
    elif n_units <= 12:
        alpha = 0.80
    else:
        alpha = 0.70

    if n_units <= 8:
        cmap = matplotlib.colormaps.get_cmap("Set2")
        colors = [cmap(i / max(1, n_units - 1)) if n_units > 1 else cmap(0) for i in range(n_units)]
    elif n_units <= 12:
        cmap = matplotlib.colormaps.get_cmap("Paired")
        colors = [cmap(i / max(1, n_units - 1)) for i in range(n_units)]
    else:
        cmap = matplotlib.colormaps.get_cmap("tab20")
        colors = [cmap(i / max(1, n_units - 1)) for i in range(n_units)]

    color_map = {unit: colors[i] for i, unit in enumerate(sorted_units)}
    return color_map, alpha


def get_legend_layout(n_items: int, n_rows: int = 15) -> Tuple[int, float, float]:
    """Calculate legend columns, right subplot margin, and width adjustment based on item count.

    Wraps to multiple columns only if the legend item count exceeds the vertical span of all tracks.

    Returns:
        (ncols, right_margin, width_bonus)
    """
    # Max items that fit vertically in a single column without exceeding the tracks
    max_items_per_col = max(1, n_rows)
    ncols = max(1, math.ceil(n_items / max_items_per_col))
    ncols = min(4, ncols)

    if ncols == 1:
        right_margin = 0.89
        width_bonus = 0.0
    elif ncols == 2:
        right_margin = 0.82
        width_bonus = 1.6
    elif ncols == 3:
        right_margin = 0.75
        width_bonus = 3.2
    else:
        right_margin = 0.68
        width_bonus = 4.8
    return ncols, right_margin, width_bonus


def plot_combined_all_clusters(
    chr_clusters_df: pl.DataFrame,
    scaf_clusters_df: pl.DataFrame,
    chrom_lengths: Dict[str, int],
    scaffold_lengths: Dict[str, int],
    species_name: str,
    output_png: Path,
    config: PipelineConfig,
) -> None:
    """Generate the Combined All-Clusters Karyotype Plot.

    Ported with 100% visual fidelity from plot_combined_all_clusters.R.
    """
    # 1. Backbones & Scaffolds Offset
    sorted_chroms = sort_chromosome_records(chrom_lengths)
    scaf_records, total_scaf_len = build_scaffold_offsets(scaffold_lengths)

    backbones = list(sorted_chroms)
    has_scaffolds = total_scaf_len > 0
    if has_scaffolds:
        backbones.append(("Scaffolds", total_scaf_len))

    if not backbones:
        return

    # 2. Filter Clusters by Copy Multiplier
    chr_filt = chr_clusters_df.filter(
        pl.col("total_chunk_bp") >= (config.min_copy_multiplier * pl.col("unit_size"))
    ) if len(chr_clusters_df) > 0 else chr_clusters_df

    scaf_filt = scaf_clusters_df.filter(
        pl.col("total_chunk_bp") >= (config.min_copy_multiplier * pl.col("unit_size"))
    ) if len(scaf_clusters_df) > 0 else scaf_clusters_df

    # Scaffold offset lookup
    scaf_offset_map = {r["name"]: r["offset"] for r in scaf_records}

    # Prepare plot items
    plot_items = []
    all_units = set()
    has_valid_ident = False

    score_col = config.plot_score_type

    for row in chr_filt.iter_rows(named=True):
        u = str(row["unit_id"])
        all_units.add(u)
        score = row.get("density_score", 0.0)
        ident = row.get("identity_score")
        if ident is not None and ident > 0:
            has_valid_ident = True

        chosen_score = ident if (score_col == "identity" and ident is not None) else score
        plot_items.append({
            "track": str(row["chromosome"]),
            "start": float(row["cluster_start"]),
            "end": float(row["cluster_end"]),
            "unit": u,
            "score": float(chosen_score or 0.0),
        })

    for row in scaf_filt.iter_rows(named=True):
        scaf_name = str(row.get("scaffold", row.get("chromosome", "")))
        if scaf_name in scaf_offset_map:
            off = scaf_offset_map[scaf_name]
            u = str(row["unit_id"])
            all_units.add(u)
            score = row.get("density_score", 0.0)
            ident = row.get("identity_score")
            if ident is not None and ident > 0:
                has_valid_ident = True

            chosen_score = ident if (score_col == "identity" and ident is not None) else score
            plot_items.append({
                "track": "Scaffolds",
                "start": float(row["cluster_start"]) + off,
                "end": float(row["cluster_end"]) + off,
                "unit": u,
                "score": float(chosen_score or 0.0),
            })

    backbone_names = set(name for name, _ in backbones)
    valid_plot_items = [c for c in plot_items if c["track"] in backbone_names]
    all_units = sorted(set(c["unit"] for c in valid_plot_items))
    color_map, alpha_val = get_color_palette(all_units)
    n_rows = len(backbones)

    ncols, right_margin, width_bonus = get_legend_layout(len(all_units), n_rows=n_rows)
    # Responsive plot sizing: default 16.0 inches plus extra width for multi-column legends
    plot_width = float(config.custom_plot_width) if config.custom_plot_width is not None else (16.0 + width_bonus)
    plot_height = max(5.5, n_rows * 0.48 + 1.8)

    y_axis_label = "Sequence Identity Score" if (score_col == "identity" and has_valid_ident) else "Density Score"

    fig, axes = plt.subplots(
        nrows=n_rows,
        ncols=1,
        figsize=(plot_width, plot_height),
        sharex=True,
        gridspec_kw={"hspace": 0.35},
    )
    if n_rows == 1:
        axes = [axes]

    max_len_mb = max(length for _, length in backbones) / 1e6

    for idx, (track_name, track_length) in enumerate(backbones):
        ax = axes[idx]
        track_len_mb = track_length / 1e6

        # Backbone segment (scaffold track slightly darker)
        bb_color = "#555555" if track_name == "Scaffolds" else "grey"
        ax.plot([0, track_len_mb], [0, 0], color=bb_color, linewidth=1.0 if track_name == "Scaffolds" else 0.9, zorder=1)

        # Scaffold tick boundaries if Scaffolds track (darker hairline ticks)
        if track_name == "Scaffolds" and scaf_records:
            b_mbs = [scaf["end_offset"] / 1e6 for scaf in scaf_records]
            ax.vlines(b_mbs, ymin=0, ymax=1000, color="#111111", linewidth=0.08, alpha=0.85, zorder=2)
            ax.vlines(b_mbs, ymin=-30, ymax=60, color="#111111", linewidth=0.08, alpha=0.85, zorder=2)

        # Plot cluster rectangles on this track
        track_clusters = [c for c in plot_items if c["track"] == track_name]
        for cl in track_clusters:
            x_min = cl["start"] / 1e6
            x_max = cl["end"] / 1e6
            width = max(x_max - x_min, 0.001)
            y_val = cl["score"]
            c_color = color_map.get(cl["unit"], "blue")

            rect = patches.Rectangle(
                (x_min, 0),
                width,
                y_val,
                facecolor=c_color,
                edgecolor="none",
                alpha=alpha_val,
                zorder=3,
            )
            ax.add_patch(rect)

        # Styling
        ax.set_ylim(-35, 1050)
        ax.set_yticks([0, 500, 1000])
        ax.tick_params(axis="y", labelsize=7.5)
        ax.grid(axis="y", color="grey", linestyle="-", linewidth=0.3, alpha=0.3)
        ax.grid(axis="x", color="grey", linestyle="--", linewidth=0.35, alpha=0.3)

        # Label on left (consistently aligned with breathing room)
        ax.set_ylabel(track_name, rotation=0, ha="right", va="center", fontsize=8.5, fontweight="bold")
        ax.yaxis.set_label_coords(-0.04, 0.5)
        for spine in ax.spines.values():
            spine.set_visible(False)

    axes[-1].tick_params(axis="x", colors="#444444", labelsize=8.5, bottom=True)
    axes[-1].set_xlim(0, max_len_mb * 1.04)
    axes[-1].set_xlabel("Genomic Position (Mb)", fontsize=11, fontweight="bold", labelpad=10)

    # Title & Subtitle on top (aligned left with ample padding)
    title = f"{species_name} — Centromeric Satellite Karyotype (All Clusters: Chromosomes & Scaffolds)"
    subtitle = f"Filtered clusters (total_chunk_bp >= {config.min_copy_multiplier} * unit_size) across {len(all_units)} satellite families [Metric: {y_axis_label}]"
    caption = "Bottom row = Concatenated scaffolds (black tick marks = scaffold boundaries)"

    fig.text(0.10, 0.985, title, fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(0.10, 0.955, subtitle, fontsize=10, color="#555555", ha="left", va="top")
    fig.text(0.10, 0.015, caption, ha="left", va="bottom", fontsize=8.5, color="#666666", style="italic")

    # Legend
    legend_patches = [
        patches.Patch(facecolor=color_map[u], edgecolor="none", label=u)
        for u in sorted(all_units)
    ]
    if legend_patches:
        fig.legend(
            handles=legend_patches,
            title="Satellite Unit",
            loc="center left",
            bbox_to_anchor=(right_margin + 0.012, 0.5),
            ncols=ncols,
            columnspacing=0.8,
            handlelength=0.9,
            handletextpad=0.3,
            frameon=False,
            fontsize=8,
            title_fontsize=9,
        )

    plt.subplots_adjust(left=0.10, right=right_margin, top=0.93, bottom=0.07)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_combined_yaml(
    yaml_file: Path,
    scaf_clusters_df: pl.DataFrame,
    chrom_lengths: Dict[str, int],
    scaffold_lengths: Dict[str, int],
    species_name: str,
    output_png: Path,
    config: PipelineConfig,
) -> None:
    """Generate the Combined YAML Summary Karyotype Plot.

    Ported with 100% visual fidelity from plot_combined_yaml.R.
    """
    if not yaml_file.exists():
        return

    with open(yaml_file, "r", encoding="utf-8") as f:
        yaml_data = yaml.safe_load(f)

    if not yaml_data or "species" not in yaml_data:
        return

    sp_info = yaml_data["species"]
    chrom_list = sp_info.get("chromosomes", [])

    # Extract chromosome clusters from YAML
    chr_rows = []
    yaml_chrom_lengths: Dict[str, int] = {}

    for c_obj in chrom_list:
        chr_name = str(c_obj.get("accession", ""))
        c_len = int(c_obj.get("chromosome_length", chrom_lengths.get(chr_name, 0)))
        yaml_chrom_lengths[chr_name] = c_len

        # Primary candidate
        cent = c_obj.get("centromere")
        if cent:
            u_name = str(cent.get("unit_name", cent.get("unit_length", "")))
            chr_rows.append({
                "chr": chr_name,
                "start": float(cent.get("start", 0)),
                "end": float(cent.get("end", 0)),
                "unit": u_name,
                "score": float(cent.get("density_score", 0)),
                "identity_score": float(cent["identity_score"]) if cent.get("identity_score") is not None else None,
                "candidate_type": "Primary Candidate",
            })

        # Alternate candidates
        alts = c_obj.get("alternate_centromeres", [])
        for alt in alts:
            u_name = str(alt.get("unit_name", alt.get("unit_length", "")))
            chr_rows.append({
                "chr": chr_name,
                "start": float(alt.get("start", 0)),
                "end": float(alt.get("end", 0)),
                "unit": u_name,
                "score": float(alt.get("density_score", 0)),
                "identity_score": float(alt["identity_score"]) if alt.get("identity_score") is not None else None,
                "candidate_type": "Alternate Candidate",
            })

    # Scaffold clusters
    sorted_chroms = sort_chromosome_records(yaml_chrom_lengths if yaml_chrom_lengths else chrom_lengths)
    scaf_records, total_scaf_len = build_scaffold_offsets(scaffold_lengths)
    scaf_offset_map = {r["name"]: r["offset"] for r in scaf_records}

    # Filter scaffold clusters to only primary/secondary candidate clusters
    if "classification" in scaf_clusters_df.columns:
        scaf_candidates_df = scaf_clusters_df.filter(
            pl.col("classification").is_in(["primary", "secondary"])
        )
    else:
        scaf_candidates_df = scaf_clusters_df

    scaf_rows = []
    for row in scaf_candidates_df.iter_rows(named=True):
        scaf_name = str(row.get("scaffold", row.get("chromosome", "")))
        if scaf_name in scaf_offset_map:
            off = scaf_offset_map[scaf_name]
            scaf_rows.append({
                "chr": "Scaffolds",
                "start": float(row["cluster_start"]) + off,
                "end": float(row["cluster_end"]) + off,
                "unit": str(row["unit_id"]),
                "score": float(row.get("density_score", 0)),
                "identity_score": float(row["identity_score"]) if row.get("identity_score") is not None else None,
                "candidate_type": "Alternate Candidate",
            })

    all_clusters = chr_rows + scaf_rows
    if not all_clusters:
        return

    backbones = list(sorted_chroms)
    has_scaffolds = total_scaf_len > 0
    if has_scaffolds:
        backbones.append(("Scaffolds", total_scaf_len))

    backbone_names = set(name for name, _ in backbones)
    # Only keep units actually present and drawn on the visible tracks
    valid_clusters = [c for c in all_clusters if c["chr"] in backbone_names]
    all_units = sorted(set(c["unit"] for c in valid_clusters))
    color_map, alpha_val = get_color_palette(all_units)
    n_rows = len(backbones)

    ncols, right_margin, width_bonus = get_legend_layout(len(all_units) + 1, n_rows=n_rows)
    # Responsive plot sizing: default 16.0 inches plus extra width for multi-column legends
    plot_width = float(config.custom_plot_width) if config.custom_plot_width is not None else (16.0 + width_bonus)
    plot_height = max(5.5, n_rows * 0.48 + 1.8)

    score_col = config.plot_score_type
    has_valid_ident = any(c.get("identity_score") is not None for c in all_clusters)
    y_axis_label = "Sequence Identity Score" if (score_col == "identity" and has_valid_ident) else "Density Score"

    fig, axes = plt.subplots(
        nrows=n_rows,
        ncols=1,
        figsize=(plot_width, plot_height),
        sharex=True,
        gridspec_kw={"hspace": 0.35},
    )
    if n_rows == 1:
        axes = [axes]

    max_len_mb = max(length for _, length in backbones) / 1e6

    # Unit labels per chromosome for superprimaries
    best_candidates = [c for c in all_clusters if c["candidate_type"] == "Primary Candidate" and c["chr"] != "Scaffolds"]
    super_labels: Dict[str, List[str]] = defaultdict(list)
    for b in best_candidates:
        if b["unit"] not in super_labels[b["chr"]]:
            super_labels[b["chr"]].append(b["unit"])

    for idx, (track_name, track_length) in enumerate(backbones):
        ax = axes[idx]
        track_len_mb = track_length / 1e6

        # Backbone segment (scaffold track slightly darker)
        bb_color = "#555555" if track_name == "Scaffolds" else "grey"
        ax.plot([0, track_len_mb], [0, 0], color=bb_color, linewidth=1.0 if track_name == "Scaffolds" else 0.9, zorder=1)

        # Scaffold tick boundaries (darker hairline ticks)
        if track_name == "Scaffolds" and scaf_records:
            b_mbs = [scaf["end_offset"] / 1e6 for scaf in scaf_records]
            ax.vlines(b_mbs, ymin=0, ymax=1000, color="#111111", linewidth=0.08, alpha=0.85, zorder=2)
            ax.vlines(b_mbs, ymin=-30, ymax=60, color="#111111", linewidth=0.08, alpha=0.85, zorder=2)

        # Clusters on this track (draw Alternate first, Best Candidate on top)
        track_clusters = [c for c in all_clusters if c["chr"] == track_name]
        track_clusters.sort(key=lambda c: 1 if c.get("candidate_type") == "Primary Candidate" else 0)
        for cl in track_clusters:
            x_min = cl["start"] / 1e6
            x_max = cl["end"] / 1e6
            width = max(x_max - x_min, 0.001)

            ident = cl.get("identity_score")
            chosen_score = ident if (score_col == "identity" and ident is not None) else cl["score"]
            c_color = color_map.get(cl["unit"], "blue")
            is_best = cl.get("candidate_type") == "Primary Candidate"

            rect = patches.Rectangle(
                (x_min, 0),
                width,
                chosen_score,
                facecolor=c_color,
                edgecolor="none",
                alpha=alpha_val,
                zorder=4 if is_best else 3,
            )
            ax.add_patch(rect)

            # Dot marker for Best Candidate
            if is_best:
                mid_mb = (x_min + x_max) / 2.0
                ax.plot(mid_mb, 0, marker="o", markersize=4.5, color="black", zorder=5)

        # Right-side text label for superprimaries
        if track_name in super_labels:
            lbl_str = ", ".join(super_labels[track_name])
            ax.text(
                track_len_mb + (max_len_mb * 0.008),
                600,
                lbl_str,
                fontsize=8.5,
                fontweight="bold",
                color="#222222",
                va="center",
                ha="left",
            )

        # Styling
        ax.set_ylim(-35, 1050)
        ax.set_yticks([0, 500, 1000])
        ax.tick_params(axis="y", labelsize=7.5)
        ax.grid(axis="y", color="grey", linestyle="-", linewidth=0.3, alpha=0.3)
        ax.grid(axis="x", color="grey", linestyle="--", linewidth=0.35, alpha=0.3)

        ax.set_ylabel(track_name, rotation=0, ha="right", va="center", fontsize=8.5, fontweight="bold")
        ax.yaxis.set_label_coords(-0.04, 0.5)
        for spine in ax.spines.values():
            spine.set_visible(False)

    axes[-1].tick_params(axis="x", colors="#444444", labelsize=8.5, bottom=True)
    axes[-1].set_xlim(0, max_len_mb * 1.08)
    axes[-1].set_xlabel("Genomic Position (Mb)", fontsize=11, fontweight="bold", labelpad=10)

    title = f"{species_name} — Centromeric Satellite Summary Karyotype (Chromosomes & Scaffolds)"
    subtitle = f"Primary candidates (dot marker) & alternate candidates across {len(all_units)} satellite families [Metric: {y_axis_label}]"
    caption = "Bottom row = Concatenated scaffolds (black tick marks = scaffold boundaries)"

    fig.text(0.10, 0.985, title, fontsize=13, fontweight="bold", ha="left", va="top")
    fig.text(0.10, 0.955, subtitle, fontsize=10, color="#555555", ha="left", va="top")
    fig.text(0.10, 0.015, caption, ha="left", va="bottom", fontsize=8.5, color="#666666", style="italic")

    # Legend
    legend_patches = [
        patches.Patch(facecolor=color_map[u], edgecolor="none", label=u)
        for u in all_units
    ]
    # Add Dot marker legend item
    marker_line = matplotlib.lines.Line2D(
        [], [], color="black", marker="o", linestyle="None", markersize=6, label="Primary Candidate"
    )
    legend_handles = legend_patches + [marker_line]

    fig.legend(
        handles=legend_handles,
        title="Unit Size / Family",
        loc="center left",
        bbox_to_anchor=(right_margin + 0.012, 0.5),
        ncols=ncols,
        columnspacing=0.8,
        handlelength=0.9,
        handletextpad=0.3,
        frameon=False,
        fontsize=8,
        title_fontsize=9,
    )

    plt.subplots_adjust(left=0.10, right=right_margin, top=0.93, bottom=0.07)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
