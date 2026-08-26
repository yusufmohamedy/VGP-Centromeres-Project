#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# run_example.sh
# Demonstrates running the centromere-pipeline on the toy example dataset.
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$SCRIPT_DIR"
PREFIX="Unicorn"

echo ">>> Running centromere-pipeline on toy example dataset..."

centromere-pipeline "$BASE_DIR/consensus_outputs/${PREFIX}-consensus_sat.bed" \
    -C "$BASE_DIR/${PREFIX}_chromosome_lengths.tsv" \
    -Cs "$BASE_DIR/${PREFIX}_scaffold_lengths.tsv" \
    -s "$PREFIX" \
    -o "$BASE_DIR/pipeline_outputs" \
    -consensus_fa "$BASE_DIR/consensus_outputs/${PREFIX}-sat_clusters.fa" \
    -P

echo ">>> Done! Check $BASE_DIR/pipeline_outputs for results."
