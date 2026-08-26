"""
Pipeline configuration module.

Defines the PipelineConfig dataclass holding all configurable thresholds,
tuning hyperparameters, and input/output path specifications.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class PipelineConfig:
    """Configuration container for the centromere satellite pipeline."""

    # Required file paths
    bed_file: Path
    chroms_file: Path
    scaffolds_file: Path

    # Sample & Output Metadata
    species: str = ""
    output_dir: Optional[Path] = None
    consensus_fa: Optional[Path] = None
    fasta_file: Optional[Path] = None

    def __post_init__(self):
        if self.consensus_fa is None and self.fasta_file is not None:
            self.consensus_fa = self.fasta_file
        elif self.consensus_fa is not None and self.fasta_file is None:
            self.fasta_file = self.consensus_fa

    # Filtering & Copy-number Parameters
    min_unit_size: int = 0
    min_copy_multiplier: int = 20

    # Classification Thresholds (matching all_clusters.R defaults passed from scaffold_pipeline.sh)
    min_density: float = 300.0
    low_score_primary: float = 500.0
    min_size_pct: float = 0.25
    other_region_pct: float = 0.25
    identity_tiebreak_pct: float = 0.70

    # Kadane Algorithm Tuning Parameters
    max_gap: float = 1e7
    min_bp: float = 0.0
    start_gap: float = 1e3
    target_gap: float = 1e6
    window_bp: int = 100_000
    min_coef: float = 0.1
    max_coef: float = 1.0
    target_coef_override: Optional[float] = None

    # Visualization Options
    auto_plot: bool = False
    plot_score_type: str = "density"  # "density" or "identity"
    custom_plot_width: Optional[float] = None


    def resolve_output_dir(self, base_dir: Optional[Path] = None) -> Path:
        """Resolve output directory, falling back to ./output/<species> if not specified."""
        if self.output_dir is not None:
            return Path(self.output_dir)
        species_name = self.species if self.species else "unknown_species"
        if base_dir is not None:
            return base_dir / "output" / species_name
        return Path("output") / species_name
