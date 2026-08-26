"""
Test suite validating the end-to-end execution of the self-contained toy example dataset.
"""

from pathlib import Path
import yaml
import pytest

from centromere_pipeline.config import PipelineConfig
from centromere_pipeline.cli import run_pipeline


def test_toy_example_end_to_end(tmp_path: Path):
    """Run pipeline on the rich Unicorn example dataset and assert all features."""
    example_dir = Path(__file__).resolve().parent.parent / "example"
    bed_file = example_dir / "consensus_outputs" / "Unicorn-consensus_sat.bed"
    chroms_file = example_dir / "Unicorn_chromosome_lengths.tsv"
    scaffolds_file = example_dir / "Unicorn_scaffold_lengths.tsv"
    fasta_file = example_dir / "consensus_outputs" / "Unicorn-sat_clusters.fa"

    assert bed_file.exists(), f"Missing Unicorn BED: {bed_file}"
    assert chroms_file.exists(), f"Missing Unicorn chroms: {chroms_file}"
    assert scaffolds_file.exists(), f"Missing Unicorn scaffolds: {scaffolds_file}"
    assert fasta_file.exists(), f"Missing Unicorn fasta: {fasta_file}"

    out_dir = tmp_path / "unicorn_outputs"
    config = PipelineConfig(
        bed_file=bed_file,
        chroms_file=chroms_file,
        scaffolds_file=scaffolds_file,
        species="Unicorn",
        output_dir=out_dir,
        auto_plot=True,
    )

    results = run_pipeline(config)

    # 1. Output files exist with correct species prefix
    assert results["chr_clusters"].exists()
    assert results["scaf_clusters"].exists()
    assert results["yaml_file"].exists()
    assert results["summary_plot"].exists()
    assert results["all_clusters_plot"].exists()

    # 2. Validate YAML content
    with open(results["yaml_file"], "r") as f:
        data = yaml.safe_load(f)

    species_info = data["species"]
    assert species_info["name"] == "Unicorn"
    assert species_info["total_chromosomes"] == 8
    assert species_info["total_scaffold_length"] == 35_500_000

    # Best candidate validation
    best_cand = species_info["best_candidate"]
    assert best_cand["unit_name"] == "sat-1-171"
    assert best_cand["unit_length"] == 171
    assert "consensus_sequence" in best_cand
    assert best_cand["average_centromere_candidate_identity"] > 900  # 0-1000 scale
    # sat-1-171 wins on 6 chromosomes; chr7 is sat-2-340 and chr8 is sat-9-240
    assert best_cand["number_of_chromosomes_where_it_is_a_centromeric_candidate"] == 6

    # Alternative candidates validation
    alt_cands = species_info.get("alternative_candidates", [])
    alt_unit_names = [a["unit_name"] for a in alt_cands]
    assert "sat-2-340" in alt_unit_names   # Wins primary on chr7
    assert "sat-9-240" in alt_unit_names   # Wins primary on chr8
    assert "sat-5-171" in alt_unit_names   # Same unit size as sat-1-171 (tiebreaker loser on chr3)

    # chr7: sat-2-340 wins primary (different satellite unit size wins entire chromosome)
    chroms = {c["accession"]: c for c in species_info["chromosomes"]}
    # sat-1-171 IS present (~107kb) but too small to override avg_primary_bp threshold (~102kb)
    # -> shows as alternate with other_unit flag on chr7
    assert chroms["chr7"]["centromere"]["unit_name"] == "sat-2-340"
    assert chroms["chr7"]["centromere"]["unit_length"] == 340
    assert chroms["chr7"]["centromere"]["centromere_position"] == "metacentric"
    assert chroms["chr7"]["centromere"].get("centromere_uncertain") == "other unit"
    chr7_alts = chroms["chr7"].get("alternate_centromeres", [])
    chr7_alt_units = [a["unit_name"] for a in chr7_alts]
    assert "sat-1-171" in chr7_alt_units  # Best candidate present but too small to win

    # chr8: sat-9-240 wins primary (different satellite unit size wins entire chromosome)
    assert chroms["chr8"]["centromere"]["unit_name"] == "sat-9-240"
    assert chroms["chr8"]["centromere"]["unit_length"] == 240
    assert chroms["chr8"]["centromere"]["centromere_position"] == "acrocentric"

    # chr3: SAME-UNIT identity tiebreaker — two sat-1-171 clusters at similar sizes
    # Cluster A (4.2 Mb, 340kb, identity 984) wins over Cluster B (18 Mb, 315kb, identity 921)
    # because 315/340 = 92.6% >= 70% tiebreak threshold AND identity 984 > 921
    assert chroms["chr3"]["centromere"]["unit_name"] == "sat-1-171"
    assert chroms["chr3"]["centromere"]["identity_score"] > 970  # High-identity winner
    # The lower-identity sat-1-171 cluster appears as 'other region' alternate (same unit, diff locus)
    assert chroms["chr3"]["centromere"].get("centromere_uncertain") is not None
    assert "other region" in chroms["chr3"]["centromere"].get("centromere_uncertain", "")
    chr3_alts = chroms["chr3"].get("alternate_centromeres", [])
    chr3_alt_units = [a["unit_name"] for a in chr3_alts]
    assert "sat-1-171" in chr3_alt_units  # Same family, different locus (tiebreak loser)
    sat1_alt = next(a for a in chr3_alts if a["unit_name"] == "sat-1-171")
    assert sat1_alt["identity_score"] < chroms["chr3"]["centromere"]["identity_score"]  # Lower identity loser

    # Validate all 4 chromosome shapes
    assert chroms["chr1"]["centromere"]["centromere_position"] == "metacentric"
    assert chroms["chr2"]["centromere"]["centromere_position"] == "submetacentric"
    assert chroms["chr3"]["centromere"]["centromere_position"] == "acrocentric"
    assert chroms["chr4"]["centromere"]["centromere_position"] == "telocentric"

    # Validate all uncertainty flags
    # 1. other_unit flag on chr1 (sat-2-340 secondary array >= 25% of primary)
    assert chroms["chr1"]["centromere"].get("centromere_uncertain") == "other unit"
    # 2. clean on chr2
    assert "centromere_uncertain" not in chroms["chr2"]["centromere"]
    # 3. other_region on chr3 — same-unit tiebreaker: lower-identity sat-1-171 cluster at 18 Mb
    #    (same unit_id, different locus, within 70% size ratio -> other_region flag)
    assert "other region" in chroms["chr3"]["centromere"].get("centromere_uncertain", "")
    # 4. short flag on chr5 (< 25% of avg candidate length)
    assert chroms["chr5"]["centromere"].get("centromere_uncertain") == "short"
    # 5. low_score flag on chr6 (density score < 500)
    assert chroms["chr6"]["centromere"].get("centromere_uncertain") == "low score"

    # Validate identity score on 0-1000 scale throughout
    for chr_name, chr_data in chroms.items():
        ident = chr_data["centromere"].get("identity_score")
        assert ident is not None
        assert 500 <= ident <= 1000, f"Identity score {ident} not in 0-1000 range for {chr_name}"

    # Validate All-Clusters BED contains wide low-density cluster and minor families
    # that do NOT appear in the summary YAML
    with open(results["chr_clusters"], "r") as f:
        chr_bed_lines = [line.strip().split("\t") for line in f if line.strip()]

    # Check wide low-density cluster on chr4 (~2.5 Mb span, density ~100.0)
    wide_clusters = [
        row for row in chr_bed_lines
        if row[0] == "chr4" and row[3] == "sat-8-85" and float(row[4]) <= 150.0
    ]
    assert len(wide_clusters) >= 1, "Wide low-density cluster on chr4 not found in all-clusters BED"
    wide_cl = wide_clusters[0]
    span = float(wide_cl[2]) - float(wide_cl[1])
    assert span > 2_000_000, f"Wide cluster span {span} bp is not > 2 Mb"
    # Verify it has no primary/secondary classification
    assert wide_cl[11] in ("", '""'), f"Expected empty classification, got {wide_cl[11]}"
    assert wide_cl[11] not in ("primary", "secondary")

    # Verify novel unit sat-10-450 is present in all-clusters BED
    sat10_clusters = [row for row in chr_bed_lines if row[3] == "sat-10-450"]
    assert len(sat10_clusters) >= 1, "sat-10-450 should be present in all-clusters BED"
    # But NOT in summary YAML
    yaml_all_units = [best_cand["unit_name"]] + [a["unit_name"] for a in alt_cands]
    assert "sat-10-450" not in yaml_all_units, "Minor family sat-10-450 should not be in summary YAML"
