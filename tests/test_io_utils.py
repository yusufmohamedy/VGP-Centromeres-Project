"""
Unit tests for io_utils module.
"""

from pathlib import Path
import polars as pl
import pytest

from centromere_pipeline.io_utils import (
    auto_detect_fasta,
    derive_species_name,
    detect_bed_identity_column,
    filter_bed_by_sequences,
    load_fasta_sequences,
    load_sequence_lengths,
    read_bed_file,
)


def test_derive_species_name():
    assert derive_species_name(Path("aSpeBom-consensus_sat.bed")) == "aSpeBom"
    assert derive_species_name(Path("bTaeGut_consensus.tsv")) == "bTaeGut"
    assert derive_species_name(Path("mDelDel-sat_clusters.fa")) == "mDelDel"
    assert derive_species_name(Path("sample123.bed")) == "sample123"


def test_load_sequence_lengths(tmp_path: Path):
    len_file = tmp_path / "test_lengths.tsv"
    len_file.write_text("chr1\t1000000\nchr2\t500000\n# comment\nchr3  250000\n")

    lengths = load_sequence_lengths(len_file)
    assert lengths == {"chr1": 1000000, "chr2": 500000, "chr3": 250000}

    # Empty / non-existent
    assert load_sequence_lengths(tmp_path / "non_existent.tsv") == {}
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("")
    assert load_sequence_lengths(empty_file) == {}


def test_load_fasta_sequences(tmp_path: Path):
    fa_file = tmp_path / "test.fa"
    fa_file.write_text(">sat-1-170 monomer repeat\nACGTACGT\nACGT\n>sat-2-350\nGGCC\n")

    seqs = load_fasta_sequences(fa_file)
    assert seqs == {
        "sat-1-170": "ACGTACGTACGT",
        "sat-2-350": "GGCC",
    }
    assert load_fasta_sequences(None) == {}


def test_read_bed_file_4col(tmp_path: Path):
    bed_file = tmp_path / "test_4col.bed"
    bed_file.write_text("chr1\t100\t200\tsat-1-100\nchr1\t300\t400\tsat-1-100\n")

    assert not detect_bed_identity_column(bed_file)
    df, has_id = read_bed_file(bed_file)
    assert not has_id
    assert len(df) == 2
    assert df.columns == ["chromosome", "start", "end", "satellite_name"]
    assert df["chromosome"].to_list() == ["chr1", "chr1"]
    assert df["start"].to_list() == [100, 300]


def test_read_bed_file_5col_identity(tmp_path: Path):
    bed_file = tmp_path / "test_5col.bed"
    bed_file.write_text("chr1\t100\t200\tsat-1-100\t95.5\nchr1\t300\t400\tsat-1-100\t88.0\n")

    assert detect_bed_identity_column(bed_file)
    df, has_id = read_bed_file(bed_file)
    assert has_id
    assert len(df) == 2
    assert df.columns == ["chromosome", "start", "end", "satellite_name", "identity_score"]
    assert df["identity_score"].to_list() == [95.5, 88.0]


def test_filter_bed_by_sequences():
    df = pl.DataFrame({
        "chromosome": ["chr1", "chr2", "NW_123"],
        "start": [100, 200, 300],
        "end": [150, 250, 350],
        "satellite_name": ["sat-1-50", "sat-1-50", "sat-1-50"],
    })
    filtered = filter_bed_by_sequences(df, {"chr1", "chr2"})
    assert len(filtered) == 2
    assert "NW_123" not in filtered["chromosome"].to_list()


def test_auto_detect_fasta(tmp_path: Path):
    bed_dir = tmp_path / "inputs"
    bed_dir.mkdir()
    bed_file = bed_dir / "sample_test.bed"
    bed_file.write_text("chr1\t10\t20\tsat-1-100\n")

    # 1. Test detection in consensus_outputs
    cons_dir = tmp_path / "consensus_outputs"
    cons_dir.mkdir()
    fa_file = cons_dir / "sample_test-sat_clusters.fa"
    fa_file.write_text(">sat-1-100\nACGT\n")

    detected = auto_detect_fasta(bed_file, "sample_test", script_dir=tmp_path)
    assert detected == fa_file.resolve()

    # 2. Test scaffolds dir is NOT detected
    scaf_dir = tmp_path / "scaffolds"
    scaf_dir.mkdir()
    scaf_fa = scaf_dir / "other_species-sat_clusters.fa"
    scaf_fa.write_text(">sat-1-100\nACGT\n")
    assert auto_detect_fasta(bed_file, "other_species", script_dir=tmp_path) is None
