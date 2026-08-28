#!/bin/bash

#SBATCH --job-name=pipeline
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=~/pipeline_logs/pipeline_%A_%a.out
#SBATCH --error=~/pipeline_logs/pipeline_%A_%a.err

set -eo pipefail

module load python/3.11.0-icl


# make sure to update the array, basedir etc


TASK=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ~/example_slurm_array_task_list.tsv)
group=$(echo "$TASK" | awk '{print $1}')
species=$(echo "$TASK" | awk '{print $2}')
prefix=$(echo "$TASK" | awk '{print $3}')

base_dir=~/$group/$species

format_duration() {
    local s=$1
    local mins=$((s / 60))
    local secs=$((s % 60))
    echo "${mins}m ${secs}s (${s}s)"
}

START_TIME=$(date +%s)

echo "============================================"
echo " Job ID      : ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo " Species dir : $group/$species"
echo " Prefix      : $prefix"
echo " Output dir  : $base_dir/pipeline_outputs"
echo " Started at  : $(date)"
echo "============================================"

# Input file validation
if [[ ! -d "$base_dir" ]]; then
    echo "ERROR: Species directory does not exist: $base_dir" >&2
    exit 1
fi

if [[ ! -s "$base_dir/consensus_outputs/${prefix}-consensus_sat.bed" ]]; then
    echo "ERROR: Missing or empty consensus_sat.bed file: $base_dir/consensus_outputs/${prefix}-consensus_sat.bed" >&2
    exit 1
fi

if [[ ! -s "$base_dir/${prefix}_chromosome_lengths.tsv" ]]; then
    echo "ERROR: Missing or empty chromosome_lengths.tsv file: $base_dir/${prefix}_chromosome_lengths.tsv" >&2
    exit 1
fi

# Scaffold lengths file is optional for T2T / chromosome-only assemblies (can be empty or absent)
scaf_file="$base_dir/${prefix}_scaffold_lengths.tsv"
if [[ ! -f "$scaf_file" ]]; then
    scaf_file="/dev/null"
fi

if [[ ! -s "$base_dir/consensus_outputs/${prefix}-sat_clusters.fa" ]]; then
    echo "ERROR: Missing or empty consensus fasta file: $base_dir/consensus_outputs/${prefix}-sat_clusters.fa" >&2
    exit 1
fi

mkdir -p "$base_dir/pipeline_outputs"


echo ">>> running pipeline..."

centromere-pipeline $base_dir/consensus_outputs/${prefix}-consensus_sat.bed \
    -C $base_dir/${prefix}_chromosome_lengths.tsv \
    -Cs "$scaf_file" \
    -s $prefix \
    -o $base_dir/pipeline_outputs \
    -consensus_fa $base_dir/consensus_outputs/${prefix}-sat_clusters.fa \
    -P

if [[ ! -s "$base_dir/pipeline_outputs/${prefix}_all_clusters_chr.bed" ]]; then
    echo "ERROR: ${prefix}_all_clusters_chr.bed file not produced by pipeline: $base_dir/pipeline_outputs/${prefix}_all_clusters_chr.bed" >&2
    exit 1
fi

if [[ ! -f "$base_dir/pipeline_outputs/${prefix}_all_clusters_scaffolds.bed" ]]; then
    echo "ERROR: ${prefix}_all_clusters_scaffolds.bed file not produced by pipeline: $base_dir/pipeline_outputs/${prefix}_all_clusters_scaffolds.bed" >&2
    exit 1
fi

if [[ ! -s "$base_dir/pipeline_outputs/${prefix}_centromere_summary.yaml" ]]; then
    echo "ERROR: ${prefix}_centromere_summary.yaml file not produced by pipeline: $base_dir/pipeline_outputs/${prefix}_centromere_summary.yaml" >&2
    exit 1
fi

if [[ ! -s "$base_dir/pipeline_outputs/${prefix}_centromere_summary_karyotype.png" ]]; then
    echo "ERROR: ${prefix}_centromere_summary_karyotype.png file not produced by pipeline: $base_dir/pipeline_outputs/${prefix}_centromere_summary_karyotype.png" >&2
    exit 1
fi

if [[ ! -s "$base_dir/pipeline_outputs/${prefix}_combined_all_clusters_karyotype.png" ]]; then
    echo "ERROR: ${prefix}_combined_all_clusters_karyotype.png file not produced by pipeline: $base_dir/pipeline_outputs/${prefix}_combined_all_clusters_karyotype.png" >&2
    exit 1
fi

echo ">>> pipeline completed!!!"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))



echo "============================================"
echo " Finished at     : $(date)"
echo " Total Duration  : $(format_duration $RUNTIME)"
echo " Output directory: ${base_dir}/pipeline_outputs"
echo "============================================"
