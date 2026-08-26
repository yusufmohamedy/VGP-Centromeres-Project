# Comprehensive Example Dataset (`Unicorn`)

This directory contains a complete, self-contained demonstration dataset designed to showcase and explain **all selection criteria, shape classifications, uncertainty flags, shared unit sizes, same-unit identity tiebreaking, multi-satellite family primary winners, wide low-density arrays, fragmented scaffolds with multiple clusters, and 0–1000 scale sequence identity scores** handled by the `centromere-pipeline`.

---

## 📁 Directory Structure

```
example/
├── Unicorn_chromosome_lengths.tsv         # 8 chromosomes (60Mb down to 12Mb, Total = 265 Mb)
├── Unicorn_scaffold_lengths.tsv           # 50 unplaced scaffolds (1.92Mb down to 100kb, Total = 35.5 Mb)
├── consensus_outputs/
│   ├── Unicorn-consensus_sat.bed          # 5-column BED (3,091 intervals, 0-1000 identity scores)
│   └── Unicorn-sat_clusters.fa            # Consensus FASTA sequences for 10 repeat families
├── run_example.sh                         # 1-click demonstration execution script
├── README.md                              # Detailed dataset and selection criteria guide
└── pipeline_outputs/                      # Generated results (BED, YAML, PNG)
```

---

## 🧬 Full Breakdown of Showcased Pipeline Features

### 1. Realistic Fragmented Scaffolds (50 Scaffolds, Multi-Cluster Architecture)
- **8 Chromosomes** (`chr1`–`chr8`, 265 Mb total) ranging in length from 60 Mb down to 12 Mb.
- **50 Unplaced Scaffolds** (`scaf_1`–`scaf_50`, 35.5 Mb total, ~11.7% of the genome) ranging in size from **1.92 Mb down to 100 kb**, accurately modeling a realistic draft genome assembly N50 distribution.
- **Multiple Clusters per Scaffold**:
  - `scaf_1` (1.92 Mb): `sat-1-171` (300 kb) + `sat-5-171` (120 kb, distinct family sharing 171 bp unit size).
  - `scaf_3` (1.77 Mb): `sat-1-171` (250 kb) + `sat-8-85` (80 kb).
  - `scaf_6` (1.55 Mb): `sat-4-1200` (180 kb) + `sat-6-340` (100 kb).
  - `scaf_8` (1.42 Mb): `sat-1-171` (200 kb) + `sat-3-55` (40 kb).
  - `scaf_12` (1.18 Mb): `sat-1-171` (160 kb) + `sat-9-240` (70 kb).
  - `scaf_15` (1.01 Mb): `sat-1-171` (150 kb) + `sat-7-55` (35 kb).
  - `scaf_20` (770 kb): `sat-1-171` (100 kb) + `sat-10-450` (30 kb).
- **25 Empty Scaffolds**: Half of the scaffolds (including `scaf_11`, `scaf_16`, `scaf_18`, `scaf_21`, `scaf_23`, `scaf_25`, `scaf_27`, and the tail `scaf_43`–`scaf_50`) contain **zero satellite repeats**, demonstrating realistic unplaced non-centromeric contigs.
- **Visualized as Concatenated Bottom Track**: All 50 scaffolds are concatenated end-to-end with hairline tick lines separating adjacent scaffolds.
- **Scaffold Summary Metrics Generated in YAML**:
  - `percentage_of_genome_in_scaffolds: 11.73554%`
  - `percentage_of_all_repeats_in_scaffolds: 46.90223%`
  - `percentage_of_best_candidate_in_scaffolds: 50.28494%`

---

### 2. Multi-Satellite Families & Distinct Primary Winners
Demonstrates how the pipeline handles multiple distinct satellite families that share identical monomer lengths, and how non-dominant families can win primary calls on individual chromosomes:

| Satellite Name | Unit Size (bp) | Role & Genomic Distribution |
| :--- | :--- | :--- |
| **`sat-1-171`** | **171 bp** | **Global Best Candidate**: Dominant repeat winning primary centromeres on `chr1`–`chr6` (6 chromosomes) and abundant across scaffolds. |
| **`sat-2-340`** | **340 bp** | **Alternative Candidate / Primary Winner on `chr7`**: Higher-order dimer; forms secondary array on `chr1` and outcompetes `sat-1-171` to win primary centromere on `chr7`. |
| **`sat-9-240`** | **240 bp** | **Alternative Candidate / Primary Winner on `chr8`**: Distinct repeat family winning primary centromere on `chr8` and present on scaffolds. |
| **`sat-5-171`** | **171 bp** | **Alternative Candidate**: Distinct sequence sharing the *same 171 bp unit size* with `sat-1-171`; present at 30 Mb on `chr3` and on `scaf_1` and `scaf_35`. |
| **`sat-6-340`** | **340 bp** | **Alternative Candidate**: Distinct sequence sharing the *same 340 bp unit size* with `sat-2-340`; present on `scaf_6`, `scaf_13`, `chr3`, and `chr5`. |
| **`sat-3-55`** | **55 bp** | Short pentanucleotide repeat cluster on `chr1`, `chr2`, `chr6`, and `scaf_8`. |
| **`sat-7-55`** | **55 bp** | Distinct sequence sharing the *same 55 bp unit size* with `sat-3-55`; arm repeat on `chr1`, `chr4`, `chr5`, `chr8`, `scaf_15`, and `scaf_24`. |
| **`sat-4-1200`** | **1200 bp** | Mega-satellite enriched predominantly on unplaced scaffolds (`scaf_4`, `scaf_6`). |
| **`sat-8-85`** | **85 bp** | Minor satellite repeat forming a wide low-density ribbon on `chr4`, and clusters on `chr2`, `chr7`, `scaf_3`, and `scaf_29`. |
| **`sat-10-450`** | **450 bp** | High-order repeat unit present as minor interstitial clusters on `chr1`, `chr3`, `chr7`, and `scaf_20`. |

---

### 3. Really Wide, Very Low-Density Cluster & Non-Summary Arrays
Highlights the distinction between the **All-Clusters Karyotype Plot** and the **YAML Summary Karyotype Plot**:

- **Wide Low-Density Cluster (`chr4` at 14.0–16.5 Mb)**:
  - Formed by `sat-8-85` across a **2.5 Mb continuous span** with sparse monomer spacing ($10\%$ coverage).
  - Calculated density score: **`100.0`** (exactly $100$ on the 0–1000 scale).
  - **All-Clusters Plot**: Rendered as a wide, short yellow rectangle (height = 100).
  - **Summary Plot**: **Completely excluded** because its density falls below the default quality threshold (`min_density: 500.0`) and non-candidate classification (`classification: ""`).
- **Minor Interstitial Arrays**:
  - Small arm clusters of `sat-3-55`, `sat-7-55`, `sat-10-450`, `sat-6-340`, etc., are correctly filtered from the centromere summary YAML ($< 25\%$ of primary candidate size) but remain visible on the All-Clusters plot.
- **Dynamic Multi-Column Legend**:
  - The All-Clusters plot tracks 10 satellite families across 9 tracks (8 chromosomes + 1 scaffold row). Because $10 > 9$, the legend automatically splits into **2 columns** to preserve vertical layout alignment.

---

### 4. Sequence Identity on Standard 0–1000 Scale
- In Column 5 of the BED input, sequence match quality is formatted on the standard bioinformatics **0–1000 scale** (e.g. `985` = 98.5% identity, `920` = 92.0% identity).
- The pipeline calculates weighted cluster identities and outputs them strictly on the 0–1000 scale throughout all YAML and table outputs:
  - Global metric: `average_centromere_candidate_identity: 979`
  - Per-chromosome: `identity_score: 980` (chr1), `985` (chr2), `985` (chr3), `988` (chr4), `976` (chr5), `960` (chr6), `969` (chr7), `975` (chr8).

---

### 5. Same-Family Sequence Identity Tiebreaking
Showcases how sequence identity resolves ambiguity between multiple arrays of the **same satellite family** on a chromosome:

- On **`chr3`**, two distinct arrays of **`sat-1-171`** are present:
  - **Cluster A (at 4.2 Mb)**: ~272 kb, density ~802, identity **`985`** (thinner and less dense).
  - **Cluster B (at 18 Mb)**: ~340 kb, density ~947, identity **`920`** (wider and denser).
- **Tiebreak Mechanism**: Because `272 / 340 = 80.0%` (exceeding the default 70% `identity_tiebreak_pct` threshold), size and density alone cannot decide the call. The pipeline uses sequence identity to break the tie, selecting **Cluster A at 4.2 Mb** (`identity: 985`) as the primary centromere.
- **Visual Distinction**: Cluster A is visibly skinnier and lower density in the karyotype plot, clearly demonstrating that the pipeline prioritizes sequence fidelity over raw array size when arrays are of comparable length.
- Cluster B triggers the diagnostic **`other_region`** uncertainty flag.

---

### 6. Best Candidate Present but Too Small to Override
- On **`chr7`**, `sat-2-340` (300 kb) forms the largest array and is selected as the primary centromere.
- The global best candidate (`sat-1-171`) is also present at 3.5 Mb (~100 kb).
- Because `100 kb < 25%` of the genome-wide average candidate size (~508 kb), it is **too small to override** `sat-2-340`.
- However, because `100 kb ≥ 25%` of the primary array on `chr7` (300 kb), it is retained as an alternate candidate and triggers the **`other_unit`** uncertainty flag.

---

### 7. All 4 Chromosome Morphology / Shape Classifications

| Chromosome | Length | Primary Centromere Unit | Centromere Center | Relative Position | Shape Classification |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`chr1`** | 60 Mb | `sat-1-171` (171 bp) | 30.0 Mb | **50.0%** (40%–60%) | **`metacentric`** |
| **`chr2`** | 50 Mb | `sat-1-171` (171 bp) | 16.0 Mb | **32.0%** (25%–40%) | **`submetacentric`** |
| **`chr3`** | 42 Mb | `sat-1-171` (171 bp) | 4.2 Mb | **10.0%** (5%–25%) | **`acrocentric`** |
| **`chr4`** | 35 Mb | `sat-1-171` (171 bp) | 0.7 Mb | **2.0%** (< 5%) | **`telocentric`** |
| **`chr5`** | 28 Mb | `sat-1-171` (171 bp) | 14.0 Mb | **50.0%** (40%–60%) | **`metacentric`** |
| **`chr6`** | 22 Mb | `sat-1-171` (171 bp) | 6.6 Mb | **30.0%** (25%–40%) | **`submetacentric`** |
| **`chr7`** | 18 Mb | **`sat-2-340` (340 bp)** | 9.0 Mb | **50.0%** (40%–60%) | **`metacentric`** |
| **`chr8`** | 12 Mb | **`sat-9-240` (240 bp)** | 1.5 Mb | **12.5%** (5%–25%) | **`acrocentric`** |

---

### 8. All Uncertainty Flags Fully Showcased

The pipeline assigns automated diagnostic uncertainty flags in the YAML output:

| Chromosome | Uncertainty Flag | Mechanism / Trigger Condition |
| :--- | :--- | :--- |
| **`chr1`** | **`other_unit`** | A large secondary array of a *different* satellite family (`sat-2-340`, 200 kb) exists with $\ge 25\%$ the size of the primary centromere. |
| **`chr2`** | *(none / clean)* | Clear, unambiguous single primary centromere with no competing arrays. |
| **`chr3`** | **`other_region, other_unit`** | **`other_region`**: same-family tiebreak loser (`sat-1-171`, 340 kb at 18 Mb). **`other_unit`**: competing array of `sat-5-171` at 30 Mb. |
| **`chr4`** | *(none / clean)* | Clear telocentric centromere with no competing secondary arrays. |
| **`chr5`** | **`short`** | Primary centromere array is only ~69 kb, which is $< 25\%$ of the genome-wide average candidate length (~508 kb). |
| **`chr6`** | **`low_score`** | Primary array density score is `418`, which falls below the default `low_score_primary` quality threshold of `500.0`. |
| **`chr7`** | **`other_unit`** | `sat-2-340` is primary winner, with `sat-1-171` secondary array at 3.5 Mb ($\ge 25\%$ of primary size). |
| **`chr8`** | *(none / clean)* | Clear acrocentric centromere won by `sat-9-240`. |

---

### 9. FASTA Consensus Sequence Integration
- Consensus nucleotide sequences are automatically loaded from `consensus_outputs/Unicorn-sat_clusters.fa` (or `-consensus_fa`).
- Full sequences are embedded in the generated YAML under `best_candidate` and each entry in `alternative_candidates`.

---

## 🚀 How to Run

### Command Line:
```bash
centromere-pipeline example/consensus_outputs/Unicorn-consensus_sat.bed \
    -C example/Unicorn_chromosome_lengths.tsv \
    -Cs example/Unicorn_scaffold_lengths.tsv \
    -s Unicorn \
    -o example/pipeline_outputs \
    -P
```

### Or using the script:
```bash
./example/run_example.sh
```

### Generated Output Files in `example/pipeline_outputs/`:
1. `Unicorn_all_clusters_chr.bed`: Classified chromosome clusters (31 total clusters).
2. `Unicorn_all_clusters_scaffolds.bed`: Scaffold repeat clusters (32 total clusters across 50 scaffolds).
3. `Unicorn_centromere_summary.yaml`: Complete YAML report.
4. `Unicorn_centromere_summary_karyotype.png`: Summary karyotype plot with primary dots and present-only legend.
5. `Unicorn_combined_all_clusters_karyotype.png`: All-clusters karyotype plot showing all 10 repeat families, the wide low-density array on `chr4`, and 50 scaffolds.
