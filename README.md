# VGP Satellite Centromere Identification Pipeline

A pipeline for automated identification, density-scoring, maximal subarray clustering (using Kadane's algorithm), and centromeric satellite classification across chromosomes and scaffolds.

*Created by Yusuf Mohamedy, Chenxi Zhou and Richard Durbin*

---

## Architecture Overview

```
VGP_Centromere_Identification_Pipeline/
├── pyproject.toml              # Build & dependency configuration
├── README.md                   # Documentation & CLI reference
├── example/                    # Self-contained demonstration dataset
│   ├── Unicorn_chromosome_lengths.tsv
│   ├── Unicorn_scaffold_lengths.tsv
│   ├── consensus_outputs/
│   │   ├── Unicorn-consensus_sat.bed
│   │   └── Unicorn-sat_clusters.fa
│   └── run_example.sh          # 1-click demonstration runner
├── centromere_pipeline/        # Production Python package
│   ├── __init__.py
│   ├── cli.py                  # CLI entry point, structured logging, orchestration
│   ├── config.py               # Single source of truth for pipeline parameters (dataclass)
│   ├── io_utils.py             # High-speed I/O (Polars BED reading, FASTA parser, lengths loader)
│   ├── density_clustering.py   # Vectorized interval merging, auto-coef, Kadane peel-off algorithm
│   ├── classification.py       # Candidate classification, uncertainty flags, shape, clean YAML export
│   └── plotting.py             # Combined Karyotype visualizer (All-Clusters & YAML Summary)
└── tests/
    ├── test_classification.py   # Unit tests for classification logic & YAML
    ├── test_density_clustering.py # Unit tests for Kadane algorithm & auto-coef
    ├── test_edge_cases.py       # Edge cases (empty files, missing chroms, 5-col BED)
    ├── test_example.py          # End-to-end Unicorn dataset integration test
    └── test_io_utils.py         # Unit tests for parsing & loading
```

---

## Key Features

- **High-Speed Vectorized Performance**: Core numerical loops leverage NumPy arrays and Polars DataFrames for ultra-fast interval operations and zero-copy slicing.
- **Kadane's Maximal Subarray Algorithm**: Dynamic array clustering with automated gap penalties and peel-off decomposition for robust cluster boundary identification.
- **Automated Centromere Classification**: Classifies primary and alternate candidate arrays, detects chromosome morphology (metacentric, submetacentric, acrocentric, telocentric), and flags uncertainty (short arrays, low scores, competing loci/units).
- **Dual Karyotype Visualizations**: Automatically generates publication-ready karyotype plots for both summary candidate arrays and genome-wide satellite clusters.
- **Clean CLI & Python API**: Full command-line interface (`centromere-pipeline`) and direct Python library interface (`from centromere_pipeline import run_pipeline, PipelineConfig`).
- **Self-Contained Demonstration Dataset**: Instant test dataset in `example/` demonstrating all pipeline features in under 1 second.
- **Automated Test Suite**: 21 unit and integration tests covering algorithmic edge cases, coordinate transformations, and end-to-end execution.

---

## Quickstart (Toy Example)

```bash
# 1. Install package in editable mode
pip install -e .

# 2. Run the included toy demonstration
./example/run_example.sh
```

---

## Installation

```bash
git clone https://github.com/<your-username>/VGP_Centromere_Identification_Pipeline.git
cd VGP_Centromere_Identification_Pipeline
pip install -e .
```

### Dependencies
- Python `>= 3.10`
- `numpy >= 1.24`
- `polars >= 1.0`
- `matplotlib >= 3.7`
- `pyyaml >= 6.0`
- `pytest >= 7.0` (for testing)

---

## Command-Line Usage

### Basic Execution
```bash
centromere-pipeline <sat.bed> -C <chroms_file> -Cs <scaffolds_file> [options]
```
Or via Python module:
```bash
python -m centromere_pipeline.cli <sat.bed> -C <chroms_file> -Cs <scaffolds_file> [options]
```

### Example: Running a Full Pipeline with Visualizations
```bash
centromere-pipeline scaffolds/aSpeBom-consensus_sat.bed \
    -C scaffolds/output/aSpeBom/aSpeBom_chromosome_lengths.tsv \
    -Cs scaffolds/output/aSpeBom/aSpeBom_scaffold_lengths.tsv \
    -s aSpeBom \
    -o output/aSpeBom \
    -P
```

### Command-Line Arguments Reference

| Option | Type | Default | Description |
|---|---|---|---|
| `<sat.bed>` or `-b` | Path | *Required* | Path to satellite BED annotation file (4 or 5 columns) |
| `-C`, `--chroms` | Path | *Required* | Path to chromosome lengths file (`chr\tlength`) |
| `-Cs`, `--scaffolds` | Path | *Required* | Path to scaffold lengths file (`scaffold\tlength`) |
| `-s`, `--species` | String | Auto-derived | Species / Sample identifier (derived from BED name if omitted) |
| `-o`, `--output` | Path | `./output/<species>` | Output directory |
| `-consensus_fa`, `--consensus-fa`, `-fa` | Path | Auto-detected | Path to consensus FASTA file (`consensus_outputs/<species>-sat_clusters.fa`) |
| `-P` | Flag | `False` | Generate combined karyotype visualizations (PNG) |
| `-plot_score` | String | `density` | Metric for plotting (`density` or `identity`) |
| `-su` | Int | `0` | Minimum monomer unit size in bp (filters smaller repeats) |
| `-cm` | Int | `20` | Minimum copy number multiplier (noise filter) |
| `-min_density` | Float | `300.0` | Minimum Kadane density score for primary qualification |
| `-low_score` | Float | `500.0` | Score threshold below which primary is flagged `low_score` |
| `-min_size_pct` | Float | `0.25` | Fraction of average primary bp for `short` flag |
| `-other_region_pct` | Float | `0.25` | Fraction of primary bp for `other_region` / `other_unit` flags |
| `-id_pct` | Float | `0.70` | Identity score tiebreak size ratio threshold |
| `-g` | Float | `1e7` | Hard maximum gap limit in bp |
| `-l` | Float | `0.0` | Minimum cumulative repeat chunk length (bp) to keep cluster |
| `-sg` | Float | `1e3` | Zero-penalty gap threshold (bp) |
| `-tg` | Float | `1e6` | Gap distance where logarithmic penalty reaches target coefficient |
| `-w` | Int | `100,000` | Sliding window (+/- bp) for local density estimation |
| `-min_coef` | Float | `0.1` | Lower floor on auto-calculated penalty coefficient |
| `-max_coef` | Float | `1.0` | Upper cap on auto-calculated penalty coefficient |
| `-coef` | Float | `None` | Manual override for penalty coefficient |

---

## Python API Usage

```python
from pathlib import Path
from centromere_pipeline import PipelineConfig, run_pipeline

config = PipelineConfig(
    bed_file=Path("example/consensus_outputs/Unicorn-consensus_sat.bed"),
    chroms_file=Path("example/Unicorn_chromosome_lengths.tsv"),
    scaffolds_file=Path("example/Unicorn_scaffold_lengths.tsv"),
    species="Unicorn",
    auto_plot=True,
)

results = run_pipeline(config)
print(f"Chromosome BED saved to: {results['chr_clusters']}")
print(f"YAML Summary saved to: {results['yaml_file']}")
print(f"Elapsed Time: {results['elapsed_seconds']:.2f}s")
```

---

## Output Files

For a species `SPECIES`, the pipeline produces:
1. `<SPECIES>_all_clusters_chr.bed`: Tab-separated BED file of all classified chromosome clusters (sorted descending by `total_chunk_bp`).
2. `<SPECIES>_all_clusters_scaffolds.bed`: Tab-separated BED file of all scaffold repeat clusters (sorted descending by `total_chunk_bp`).
3. `<SPECIES>_centromere_summary.yaml`: Structured YAML specification including genome-wide summary statistics, best candidate, alternative candidates, and per-chromosome centromere calls.
4. (Optional with `-P`) `<SPECIES>_centromere_summary_karyotype.png`: Publication-quality karyotype diagram showing primary (dot marker) and alternate centromere arrays with concatenated scaffolds.
5. (Optional with `-P`) `<SPECIES>_combined_all_clusters_karyotype.png`: Karyotype diagram displaying all satellite arrays across chromosomes and scaffolds.

---

## Running the Automated Test Suite

```bash
pytest -v
```

The test suite validates algorithmic correctness, Kadane cluster boundaries, edge cases (empty files, missing scaffolds, coordinate mismatches), 5-column identity scoring, and end-to-end execution on the included `Unicorn` dataset.
