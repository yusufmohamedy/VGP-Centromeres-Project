"""
Unit tests for classification module.
"""

from pathlib import Path
import polars as pl
import pytest
import yaml

from centromere_pipeline.config import PipelineConfig
from centromere_pipeline.classification import (
    classify_clusters,
    classify_shape,
    dump_yaml_clean,
    generate_centromere_yaml,
)


def test_classify_shape():
    chrom_len = 100_000_000
    # < 5% -> telocentric
    assert classify_shape(2_000_000, chrom_len) == "telocentric"
    # > 95% -> telocentric
    assert classify_shape(98_000_000, chrom_len) == "telocentric"
    # 10% (< 25%) -> acrocentric
    assert classify_shape(10_000_000, chrom_len) == "acrocentric"
    # 85% (> 75%) -> acrocentric
    assert classify_shape(85_000_000, chrom_len) == "acrocentric"
    # 30% (25-40%) -> submetacentric
    assert classify_shape(30_000_000, chrom_len) == "submetacentric"
    # 70% (60-75%) -> submetacentric
    assert classify_shape(70_000_000, chrom_len) == "submetacentric"
    # 50% (40-60%) -> metacentric
    assert classify_shape(50_000_000, chrom_len) == "metacentric"


def test_classify_clusters_primary_and_secondary():
    raw_df = pl.DataFrame({
        "unit_id": ["sat-1-100", "sat-1-100", "sat-2-200"],
        "unit_size": [100, 100, 200],
        "chromosome": ["chr1", "chr1", "chr1"],
        "cluster_start": [10_000, 50_000, 80_000],
        "cluster_end": [20_000, 55_000, 82_000],
        "gap_before": ["NA", "30000", "25000"],
        "gap_after": ["30000", "25000", "NA"],
        "cluster_id": ["cluster_1", "cluster_2", "cluster_3"],
        "n_chunks": [10, 5, 2],
        "total_chunk_bp": [8000.0, 3500.0, 1000.0],
        "max_internal_gap": [100, 50, 20],
        "max_sum": [7500.0, 3000.0, 900.0],
        "density": [0.8, 0.7, 0.5],
        "density_score": [800.0, 700.0, 500.0],
    })

    chrom_lengths = {"chr1": 100_000}
    config = PipelineConfig(
        bed_file=Path("dummy.bed"),
        chroms_file=Path("dummy.tsv"),
        scaffolds_file=Path("dummy.tsv"),
        min_copy_multiplier=5,  # min_bp = 5 * 100 = 500 for sat-1, 1000 for sat-2
        min_density=300,
        other_region_pct=0.25,
    )

    classified = classify_clusters(raw_df, chrom_lengths, config, is_scaffolds=False)
    assert len(classified) == 3
    # Top cluster should be primary
    top_row = classified.row(0, named=True)
    assert top_row["classification"] == "primary"
    assert top_row["total_chunk_bp"] == 8000.0

    # Second cluster (3500 bp >= 0.25 * 8000 = 2000 bp) should be secondary
    sec_row = classified.row(1, named=True)
    assert sec_row["classification"] == "secondary"


def test_generate_centromere_yaml_with_fasta(tmp_path: Path):
    chr_df = pl.DataFrame({
        "chromosome": ["chr1"],
        "cluster_start": [1000],
        "cluster_end": [5000],
        "unit_id": ["sat-1-100"],
        "density_score": [950.0],
        "unit_size": [100],
        "cluster_size": [4000.0],
        "total_chunk_bp": [3800.0],
        "copy_number": [38.0],
        "chromosome_position": [0.5],
        "classification": ["primary"],
        "uncertainty": [""],
        "shape": ["metacentric"],
    })
    scaf_df = pl.DataFrame()
    chrom_lens = {"chr1": 10000}
    scaf_lens = {}
    fasta_seqs = {"sat-1-100": "ATGCATGCATGC"}

    out_yaml = tmp_path / "test_fasta_summary.yaml"
    generate_centromere_yaml(
        chr_df, scaf_df, chrom_lens, scaf_lens, "test_species", fasta_seqs, out_yaml
    )

    assert out_yaml.exists()
    with open(out_yaml) as f:
        data = yaml.safe_load(f)

    assert data["species"]["best_candidate"]["consensus_sequence"] == "ATGCATGCATGC"

