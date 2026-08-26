"""
I/O utilities module for reading and writing genomic data formats.

Supports BED intervals, chromosome length files, consensus FASTA files,
and TSV / BED table exports.
"""

from pathlib import Path
import re
from typing import Dict, List, Optional, Set, Tuple
import polars as pl


def derive_species_name(bed_path: Path) -> str:
    """Auto-derive species/sample identifier from BED filename.

    Mirrors the regex used in the original shell script:
    sed -E 's/(-consensus_sat|_consensus|-consensus|-sat_clusters|-sat|-tan)?\\.(bed|tsv|txt|fa|fasta)$//'
    """
    stem = bed_path.name
    # Strip known extensions and suffixes
    pattern = r"(-consensus_sat|_consensus|-consensus|-sat_clusters|-sat|-tan)?\.(bed|tsv|txt|fa|fasta)$"
    species = re.sub(pattern, "", stem, flags=re.IGNORECASE)
    return species


def auto_detect_fasta(
    bed_path: Path, species: str, script_dir: Optional[Path] = None
) -> Optional[Path]:
    """Auto-detect consensus sequences FASTA file from standard location patterns."""
    bed_dir = bed_path.resolve().parent
    cwd = Path.cwd()
    search_dirs = [
        bed_dir,
        cwd,
        cwd / "consensus_outputs",
        cwd / "consensus_outputs" / species,
        bed_dir.parent,
        bed_dir.parent / "consensus_outputs",
        bed_dir.parent / "consensus_outputs" / species,
        bed_dir / "consensus_outputs",
        bed_dir / "consensus_outputs" / species,
    ]
    if script_dir:
        search_dirs.extend([
            script_dir,
            script_dir / "consensus_outputs",
            script_dir / "consensus_outputs" / species,
        ])

    candidate_names = [
        f"{species}-sat_clusters.fa",
        f"{species}_sat_clusters.fa",
        "sat_clusters.fa",
        f"{species}-tancons.fa",
        f"{species}_tancons.fa",
    ]

    for d in search_dirs:
        for name in candidate_names:
            cand = d / name
            if cand.exists() and cand.is_file() and cand.stat().st_size > 0:
                return cand.resolve()

    return None


def detect_bed_identity_column(bed_path: Path) -> bool:
    """Detect whether BED file has 5+ columns with a numeric 5th column (identity score)."""
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return False

    with open(bed_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 5:
                val = parts[4].strip()
                # Check if it matches a float / integer regex
                if re.match(r"^[0-9]+(\.[0-9]+)?$", val):
                    return True
            break
    return False


def load_sequence_lengths(lengths_path: Optional[Path]) -> Dict[str, int]:
    """Load sequence lengths from a two-column tab-separated or whitespace-separated file.

    Format:
        <sequence_id> <length_in_bp>

    Returns:
        Dict mapping sequence identifier to length in bp.
    """
    lengths: Dict[str, int] = {}
    if lengths_path is None or not Path(lengths_path).exists():
        return lengths

    path = Path(lengths_path)
    if path.stat().st_size == 0:
        return lengths

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) >= 2:
                seq_id = parts[0].strip()
                try:
                    lengths[seq_id] = int(parts[1].strip())
                except ValueError:
                    continue
    return lengths


def load_fasta_sequences(fasta_path: Optional[Path]) -> Dict[str, str]:
    """Load consensus sequences from a FASTA file.

    Returns:
        Dict mapping sequence header ID (first whitespace-delimited word) to sequence string.
    """
    sequences: Dict[str, str] = {}
    if fasta_path is None or not Path(fasta_path).exists():
        return sequences

    path = Path(fasta_path)
    if path.stat().st_size == 0:
        return sequences

    current_id: Optional[str] = None
    current_seq_parts: List[str] = []

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq_parts)
                header = line[1:].strip()
                current_id = header.split()[0]
                current_seq_parts = []
            else:
                current_seq_parts.append(line)

        if current_id is not None:
            sequences[current_id] = "".join(current_seq_parts)

    return sequences


def read_bed_file(bed_path: Path) -> Tuple[pl.DataFrame, bool]:
    """Read a satellite annotation BED file into a Polars DataFrame.

    Returns:
        (df, has_identity): DataFrame with columns [chromosome, start, end, satellite_name, (optional identity_score)]
        and boolean flag indicating whether identity score column was present.
    """
    has_identity = detect_bed_identity_column(bed_path)

    if not bed_path.exists() or bed_path.stat().st_size == 0:
        schema = {
            "chromosome": pl.Utf8,
            "start": pl.Int64,
            "end": pl.Int64,
            "satellite_name": pl.Utf8,
        }
        if has_identity:
            schema["identity_score"] = pl.Float64
        return pl.DataFrame(schema=schema), has_identity

    # Read tab-delimited BED
    if has_identity:
        schema_overrides = {
            "column_1": pl.Utf8,
            "column_2": pl.Int64,
            "column_3": pl.Int64,
            "column_4": pl.Utf8,
            "column_5": pl.Float64,
        }
        df = pl.read_csv(
            bed_path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            schema_overrides=schema_overrides,
            truncate_ragged_lines=True,
        )
        df = df.select([
            pl.col("column_1").alias("chromosome"),
            pl.col("column_2").alias("start"),
            pl.col("column_3").alias("end"),
            pl.col("column_4").alias("satellite_name"),
            pl.col("column_5").alias("identity_score"),
        ])
    else:
        schema_overrides = {
            "column_1": pl.Utf8,
            "column_2": pl.Int64,
            "column_3": pl.Int64,
            "column_4": pl.Utf8,
        }
        df = pl.read_csv(
            bed_path,
            separator="\t",
            has_header=False,
            comment_prefix="#",
            schema_overrides=schema_overrides,
            truncate_ragged_lines=True,
        )
        df = df.select([
            pl.col("column_1").alias("chromosome"),
            pl.col("column_2").alias("start"),
            pl.col("column_3").alias("end"),
            pl.col("column_4").alias("satellite_name"),
        ])

    return df, has_identity


def filter_bed_by_sequences(
    bed_df: pl.DataFrame, valid_sequences: Set[str]
) -> pl.DataFrame:
    """Filter BED DataFrame to records belonging to specified sequence IDs."""
    if len(bed_df) == 0 or not valid_sequences:
        return bed_df.filter(pl.lit(False))
    return bed_df.filter(pl.col("chromosome").is_in(list(valid_sequences)))
