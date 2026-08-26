"""
Command-line interface and main orchestration entry point for the centromere pipeline.

Replaces scaffold_pipeline.sh with a unified, high-performance Python application.
"""

import argparse
import logging
from pathlib import Path
import sys
import time
from typing import Any, Dict, List, Optional

from .config import PipelineConfig
from .io_utils import (
    auto_detect_fasta,
    derive_species_name,
    load_fasta_sequences,
    load_sequence_lengths,
    read_bed_file,
)
from .density_clustering import process_dataset
from .classification import classify_clusters, generate_centromere_yaml
from .plotting import plot_combined_all_clusters, plot_combined_yaml


def setup_logger(verbose: bool = False) -> logging.Logger:
    """Configure structured console logging."""
    logger = logging.getLogger("centromere_pipeline")
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    logger.handlers.clear()

    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter(
        "[%(asctime)s] %(levelname)s: %(message)s", datefmt="%H:%M:%S"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


def build_arg_parser() -> argparse.ArgumentParser:
    """Build argument parser matching scaffold_pipeline.sh flags and defaults."""
    parser = argparse.ArgumentParser(
        description="Chromosome & Scaffold Centromeric Satellite Clustering Pipeline (Satmatch Mode)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional or flagged BED file
    parser.add_argument("bed_positional", nargs="?", default=None, help="Path to satellite BED file")
    parser.add_argument("-b", "--bed", dest="bed_flag", default=None, help="Path to satellite BED file")

    # Required sequence length files
    parser.add_argument("-C", "--chroms", dest="chroms", required=True, help="Path to chromosome lengths file (tsv: chr\\tlength)")
    parser.add_argument("-Cs", "--scaffolds", dest="scaffolds", required=True, help="Path to scaffold lengths file (tsv: scaffold\\tlength)")

    # Sample & Output Metadata
    parser.add_argument("-s", "-S", "--species", "--prefix", dest="species", default=None, help="Species/Sample identifier (default: auto-derived from BED)")
    parser.add_argument("-o", "--output", dest="output_dir", default=None, help="Custom output directory")
    parser.add_argument(
        "-consensus_fa", "--consensus-fa", "--consensus_fa", "-cfa", "-fa", "--fasta",
        "-consensus_sequences_fa", "--consensus-sequences-fa", "--consensus_sequences_fa",
        dest="consensus_fa", default=None, help="Path to sat_clusters.fa containing consensus sequences"
    )

    # Filtering & Copy-number Parameters
    parser.add_argument("-su", "-min_unit", "-min_u", dest="min_unit_size", type=int, default=0, help="Minimum unit size in bp")
    parser.add_argument("-cm", "-copy_mult", dest="min_copy_multiplier", type=int, default=20, help="Minimum copy number multiplier")

    # Classification Thresholds
    parser.add_argument("-min_density", "-min_score", dest="min_density", type=float, default=300.0, help="Minimum density score for primary qualification")
    parser.add_argument("-low_score", dest="low_score_primary", type=float, default=500.0, help="Score threshold below which primary is flagged low_score")
    parser.add_argument("-min_size_pct", dest="min_size_pct", type=float, default=0.25, help="Fraction of average primary bp for short flag")
    parser.add_argument("-other_region_pct", dest="other_region_pct", type=float, default=0.25, help="Fraction of primary bp for other_region flag")
    parser.add_argument("-id_pct", "-identity_tiebreak_pct", dest="identity_tiebreak_pct", type=float, default=0.70, help="Identity score tiebreak size ratio threshold")

    # Kadane Tuning
    parser.add_argument("-g", dest="max_gap", type=float, default=1e7, help="Hard gap limit (bp)")
    parser.add_argument("-l", dest="min_bp", type=float, default=0.0, help="Min total chunk bp to retain cluster")
    parser.add_argument("-sg", dest="start_gap", type=float, default=1e3, help="Zero-penalty gap threshold (bp)")
    parser.add_argument("-tg", dest="target_gap", type=float, default=1e6, help="Gap size where penalty reaches target_coef")
    parser.add_argument("-w", dest="window_bp", type=int, default=100_000, help="Density window +/- bp")
    parser.add_argument("-min_coef", dest="min_coef", type=float, default=0.1, help="Lower floor on auto-calculated target_coef")
    parser.add_argument("-max_coef", dest="max_coef", type=float, default=1.0, help="Upper cap on auto-calculated target_coef")
    parser.add_argument("-coef", dest="target_coef_override", type=float, default=None, help="Manual override for target_coef")

    # Visualization
    parser.add_argument("-P", dest="auto_plot", action="store_true", help="Generate combined karyotype plots (YAML & All-Clusters)")
    parser.add_argument("-plot_score", "--score-type", "-score_type", dest="plot_score_type", choices=["density", "identity"], default="density", help="Plot score metric")
    parser.add_argument("-plot_width", "--plot-width", dest="custom_plot_width", type=float, default=None, help="Custom plot width in inches (default: auto-fill screen ~16.0)")

    return parser


def run_pipeline(config: PipelineConfig, logger: Optional[logging.Logger] = None) -> Dict[str, Any]:
    """Execute the full centromere pipeline from a PipelineConfig."""
    if logger is None:
        logger = setup_logger()

    start_time = time.time()

    # Resolve BED file
    bed_path = Path(config.bed_file).resolve()
    if not bed_path.exists():
        raise FileNotFoundError(f"Satellite BED file not found: {bed_path}")

    chroms_path = Path(config.chroms_file).resolve()
    if not chroms_path.exists():
        raise FileNotFoundError(f"Chromosome lengths file not found: {chroms_path}")

    scaffolds_path = Path(config.scaffolds_file).resolve()
    if not scaffolds_path.exists():
        raise FileNotFoundError(f"Scaffold lengths file not found: {scaffolds_path}")

    # Derive species identifier
    species = config.species if config.species else derive_species_name(bed_path)
    config.species = species

    # Output directory
    output_dir = config.resolve_output_dir(Path.cwd())
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect / load FASTA consensus sequences
    fa_source = config.consensus_fa or config.fasta_file
    fasta_path = Path(fa_source).resolve() if fa_source else auto_detect_fasta(bed_path, species)
    fasta_sequences = load_fasta_sequences(fasta_path) if fasta_path else {}

    # Load sequence lengths
    chrom_lengths = load_sequence_lengths(chroms_path)
    scaffold_lengths = load_sequence_lengths(scaffolds_path)

    # Read BED dataset
    bed_df, has_identity = read_bed_file(bed_path)

    logger.info("===================================================================")
    logger.info(f"  Scaffold & Chromosome Centromeric Pipeline (Python): {species}")
    logger.info("===================================================================")
    logger.info(f"  Mode:            Satmatch (vectorized Polars/NumPy)")
    logger.info(f"  Identity Col:    {'Present (Column 5)' if has_identity else 'None'}")
    logger.info(f"  Consensus FASTA: {fasta_path if fasta_path else 'None'}")
    logger.info(f"  BED file:        {bed_path} ({len(bed_df):,} intervals)")
    logger.info(f"  Chroms file:     {chroms_path} ({len(chrom_lengths)} chromosomes)")
    logger.info(f"  Scaffolds file:  {scaffolds_path} ({len(scaffold_lengths)} scaffolds)")
    logger.info(f"  Min Unit Size:   {config.min_unit_size} bp")
    logger.info(f"  Output Dir:      {output_dir}")
    logger.info(f"  Plots enabled:   {config.auto_plot}")
    logger.info("-------------------------------------------------------------------")

    prefix_str = f"{species}_" if species else ""

    # Step 1: Process Chromosomes
    logger.info("  1. Processing Chromosome Arrays...")
    raw_chr_df = process_dataset(bed_df, set(chrom_lengths.keys()), config, has_identity)
    classified_chr_df = classify_clusters(
        raw_chr_df, chrom_lengths, config, is_scaffolds=False, has_identity=has_identity
    )
    chr_bed_out = output_dir / f"{prefix_str}all_clusters_chr.bed"
    classified_chr_df.write_csv(chr_bed_out, separator="\t", null_value="")
    logger.info(f"     Successfully wrote {len(classified_chr_df)} chromosome clusters to {chr_bed_out}")

    # Step 2: Process Scaffolds
    logger.info("  2. Processing Scaffold Arrays...")
    raw_scaf_df = process_dataset(bed_df, set(scaffold_lengths.keys()), config, has_identity)
    classified_scaf_df = classify_clusters(
        raw_scaf_df, scaffold_lengths, config, is_scaffolds=True, has_identity=has_identity
    )
    scaf_bed_out = output_dir / f"{prefix_str}all_clusters_scaffolds.bed"
    classified_scaf_df.write_csv(scaf_bed_out, separator="\t", null_value="")
    logger.info(f"     Successfully wrote {len(classified_scaf_df)} scaffold clusters to {scaf_bed_out}")

    # Step 3: Generate Centromere Summary YAML
    logger.info("  3. Generating Centromere Summary YAML (Chromosomes Only)...")
    yaml_out = output_dir / (f"{species}_centromere_summary.yaml" if species else "centromere_summary.yaml")
    generate_centromere_yaml(
        classified_chr_df,
        classified_scaf_df,
        chrom_lengths,
        scaffold_lengths,
        species,
        fasta_sequences,
        yaml_out,
        has_identity=has_identity,
    )
    logger.info(f"     Successfully generated YAML: {yaml_out}")

    # Step 4: Combined Karyotype Visualizations (-P flag)
    summary_plot_out = None
    all_clusters_plot_out = None

    if config.auto_plot:
        logger.info("  4. Generating Combined Karyotype Plots...")

        # Plot 1: Centromere Summary / YAML Plot
        summary_plot_out = output_dir / (f"{species}_centromere_summary_karyotype.png" if species else "centromere_summary_karyotype.png")
        plot_combined_yaml(
            yaml_out,
            classified_scaf_df,
            chrom_lengths,
            scaffold_lengths,
            species,
            summary_plot_out,
            config,
        )
        logger.info(f"     Generated YAML Summary Plot: {summary_plot_out}")

        # Plot 2: Combined All-Clusters Plot
        all_clusters_plot_out = output_dir / f"{prefix_str}combined_all_clusters_karyotype.png"
        plot_combined_all_clusters(
            classified_chr_df,
            classified_scaf_df,
            chrom_lengths,
            scaffold_lengths,
            species,
            all_clusters_plot_out,
            config,
        )
        logger.info(f"     Generated All-Clusters Plot: {all_clusters_plot_out}")

    elapsed = time.time() - start_time
    minutes = int(elapsed // 60)
    seconds = int(elapsed % 60)

    logger.info("===================================================================")
    logger.info(f"  Pipeline execution completed successfully for {species}!")
    logger.info(f"  Outputs saved to: {output_dir}")
    logger.info(f"  Chromosome clusters: {chr_bed_out}")
    logger.info(f"  Scaffold clusters:   {scaf_bed_out}")
    logger.info(f"  Centromere YAML:     {yaml_out}")
    if config.auto_plot:
        logger.info(f"  Centromere Summary Plot:   {summary_plot_out}")
        logger.info(f"  Combined All-Clusters Plot:{all_clusters_plot_out}")
    logger.info(f"  Total Time: {minutes}m {seconds}s ({elapsed:.2f} seconds)")
    logger.info("===================================================================")

    return {
        "species": species,
        "output_dir": output_dir,
        "chr_clusters": chr_bed_out,
        "scaf_clusters": scaf_bed_out,
        "yaml_file": yaml_out,
        "summary_plot": summary_plot_out,
        "all_clusters_plot": all_clusters_plot_out,
        "elapsed_seconds": elapsed,
    }


def main(argv: Optional[List[str]] = None) -> None:
    """Main CLI entry point."""
    if argv is None:
        argv = sys.argv[1:]

    parser = build_arg_parser()
    args = parser.parse_args(argv)

    # Determine input BED path (from positional or flagged option)
    bed_path_str = args.bed_flag or args.bed_positional
    if not bed_path_str:
        parser.print_help()
        sys.exit(1)

    config = PipelineConfig(
        bed_file=Path(bed_path_str),
        chroms_file=Path(args.chroms),
        scaffolds_file=Path(args.scaffolds),
        species=args.species or "",
        output_dir=Path(args.output_dir) if args.output_dir else None,
        consensus_fa=Path(args.consensus_fa) if args.consensus_fa else None,
        min_unit_size=args.min_unit_size,
        min_copy_multiplier=args.min_copy_multiplier,
        min_density=args.min_density,
        low_score_primary=args.low_score_primary,
        min_size_pct=args.min_size_pct,
        other_region_pct=args.other_region_pct,
        identity_tiebreak_pct=args.identity_tiebreak_pct,
        max_gap=args.max_gap,
        min_bp=args.min_bp,
        start_gap=args.start_gap,
        target_gap=args.target_gap,
        window_bp=args.window_bp,
        min_coef=args.min_coef,
        max_coef=args.max_coef,
        target_coef_override=args.target_coef_override,
        auto_plot=args.auto_plot,
        plot_score_type=args.plot_score_type,
    )

    logger = setup_logger()
    try:
        run_pipeline(config, logger)
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
