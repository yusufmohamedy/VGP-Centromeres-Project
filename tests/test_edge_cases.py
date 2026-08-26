"""
Unit tests for edge cases, missing data, empty files, and synthetic 5-column BED.
"""

from pathlib import Path
import polars as pl
import pytest

from centromere_pipeline.config import PipelineConfig
from centromere_pipeline.io_utils import read_bed_file
from centromere_pipeline.density_clustering import process_dataset
from centromere_pipeline.classification import classify_clusters, generate_centromere_yaml


def test_empty_bed_file(tmp_path: Path):
    empty_bed = tmp_path / "empty.bed"
    empty_bed.write_text("")

    chrom_lengths = {"chr1": 100000}
    scaf_lengths = {}
    config = PipelineConfig(
        bed_file=empty_bed,
        chroms_file=tmp_path / "chroms.tsv",
        scaffolds_file=tmp_path / "scafs.tsv",
    )

    bed_df, has_id = read_bed_file(empty_bed)
    assert len(bed_df) == 0

    raw_chr = process_dataset(bed_df, set(chrom_lengths.keys()), config, has_id)
    assert len(raw_chr) == 0

    classified = classify_clusters(raw_chr, chrom_lengths, config, is_scaffolds=False)
    assert len(classified) == 0
    assert "chromosome" in classified.columns
    assert "classification" in classified.columns

    yaml_out = tmp_path / "empty_summary.yaml"
    generate_centromere_yaml(
        classified, pl.DataFrame(), chrom_lengths, scaf_lengths,
        "empty_species", {}, yaml_out, has_identity=has_id
    )
    assert yaml_out.exists()


def test_no_matching_chromosomes(tmp_path: Path):
    bed_file = tmp_path / "mismatch.bed"
    bed_file.write_text("chrX\t100\t500\tsat-1-100\n")

    chrom_lengths = {"chr1": 100000, "chr2": 200000}  # chrX not present
    config = PipelineConfig(
        bed_file=bed_file,
        chroms_file=tmp_path / "chroms.tsv",
        scaffolds_file=tmp_path / "scafs.tsv",
    )

    bed_df, has_id = read_bed_file(bed_file)
    raw_chr = process_dataset(bed_df, set(chrom_lengths.keys()), config, has_id)
    assert len(raw_chr) == 0


def test_synthetic_5col_identity_dataset(tmp_path: Path):
    """Test full pipeline execution on a synthetic 5-column BED dataset with sequence identity."""
    bed_file = tmp_path / "synthetic_5col.bed"
    chrom_file = tmp_path / "chroms.tsv"
    scaf_file = tmp_path / "scaffolds.tsv"

    chrom_file.write_text("chr1\t1000000\nchr2\t800000\n")
    scaf_file.write_text("scaf_1\t50000\n")

    # Generate records: chr1 has two competing sat-1 arrays, one bigger but lower identity
    bed_content = []
    # Cluster A: start 100k-150k, ident 850
    for s in range(100_000, 150_000, 200):
        bed_content.append(f"chr1\t{s}\t{s+180}\tsat-1-170\t850")
    # Cluster B: start 400k-440k, ident 980 (higher identity, close size)
    for s in range(400_000, 440_000, 200):
        bed_content.append(f"chr1\t{s}\t{s+180}\tsat-1-170\t980")
    # Scaffold record
    bed_content.append("scaf_1\t1000\t2500\tsat-1-170\t920")

    bed_file.write_text("\n".join(bed_content) + "\n")

    config = PipelineConfig(
        bed_file=bed_file,
        chroms_file=chrom_file,
        scaffolds_file=scaf_file,
        species="synthetic_test",
        output_dir=tmp_path / "out",
        identity_tiebreak_pct=0.70,
    )

    bed_df, has_id = read_bed_file(bed_file)
    assert has_id is True
    assert "identity_score" in bed_df.columns

    chrom_lengths = {"chr1": 1000000, "chr2": 800000}
    scaf_lengths = {"scaf_1": 50000}

    raw_chr = process_dataset(bed_df, set(chrom_lengths.keys()), config, has_identity=True)
    assert len(raw_chr) >= 2
    assert "identity_score" in raw_chr.columns

    classified_chr = classify_clusters(raw_chr, chrom_lengths, config, is_scaffolds=False, has_identity=True)
    assert "identity_score" in classified_chr.columns

    # Check identity tiebreaker: Cluster B has identity 980, Cluster A has 850
    # Cluster B total_chunk_bp is ~36,000, Cluster A is ~45,000 (36k / 45k = 0.80 >= 0.70 id_pct)
    primary_row = classified_chr.filter(pl.col("classification") == "primary").row(0, named=True)
    assert primary_row["cluster_start"] == 400_000
    assert primary_row["identity_score"] == 980


def test_cli_consensus_fa_argument(tmp_path: Path):
    from centromere_pipeline.cli import build_arg_parser

    parser = build_arg_parser()
    args = parser.parse_args([
        "sample.bed",
        "-C", "chroms.tsv",
        "-Cs", "scafs.tsv",
        "-consensus_fa", str(tmp_path / "clusters.fa")
    ])
    assert args.consensus_fa == str(tmp_path / "clusters.fa")

    # Test config backwards-compatible alias
    cfg1 = PipelineConfig(
        bed_file=Path("sample.bed"),
        chroms_file=Path("chroms.tsv"),
        scaffolds_file=Path("scafs.tsv"),
        consensus_fa=tmp_path / "clusters.fa",
    )
    assert cfg1.consensus_fa == tmp_path / "clusters.fa"
    assert cfg1.fasta_file == tmp_path / "clusters.fa"

    cfg2 = PipelineConfig(
        bed_file=Path("sample.bed"),
        chroms_file=Path("chroms.tsv"),
        scaffolds_file=Path("scafs.tsv"),
        fasta_file=tmp_path / "legacy.fa",
    )
    assert cfg2.consensus_fa == tmp_path / "legacy.fa"
    assert cfg2.fasta_file == tmp_path / "legacy.fa"


def test_4col_bed_with_plotting(tmp_path: Path):
    """Ensure 4-column BED datasets without identity scores plot without KeyError."""
    from centromere_pipeline.cli import run_pipeline

    bed_file = tmp_path / "test_4col.bed"
    bed_file.write_text(
        "chr1\t100000\t150000\tsat_170\n"
        "chr1\t150050\t200000\tsat_170\n"
        "scaf_1\t10000\t25000\tsat_170\n"
    )

    chr_file = tmp_path / "chroms.tsv"
    chr_file.write_text("chr1\t1000000\n")

    scaf_file = tmp_path / "scafs.tsv"
    scaf_file.write_text("scaf_1\t100000\n")

    out_dir = tmp_path / "output"

    config = PipelineConfig(
        bed_file=bed_file,
        chroms_file=chr_file,
        scaffolds_file=scaf_file,
        species="test_4col_species",
        output_dir=out_dir,
        auto_plot=True,
    )

    results = run_pipeline(config)
    assert results["yaml_file"].exists()
    assert results["summary_plot"] is not None and results["summary_plot"].exists()
    assert results["all_clusters_plot"] is not None and results["all_clusters_plot"].exists()
