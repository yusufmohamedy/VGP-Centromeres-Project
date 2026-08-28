#!/bin/bash

#SBATCH --job-name=consensus
#SBATCH --array=1-10
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=~/consensus_logs/consensus_%A_%a.out
#SBATCH --error=~/consensus_logs/consensus_%A_%a.err

set -eo pipefail



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
echo " Output dir  : $base_dir/consensus_outputs"
echo " Started at  : $(date)"
echo "============================================"

# Input file validation
if [[ ! -d "$base_dir" ]]; then
    echo "ERROR: Species directory does not exist: $base_dir" >&2
    exit 1
fi

if [[ ! -s "$base_dir/${prefix}.1ano" ]]; then
    echo "ERROR: Missing or empty 1ano file: $base_dir/${prefix}.1ano" >&2
    exit 1
fi

if [[ ! -s "$base_dir/${prefix}.fa.gz" ]]; then
    echo "ERROR: Missing or empty fa.gz file: $base_dir/${prefix}.fa.gz" >&2
    exit 1
fi

if [[ ! -s "$base_dir/${prefix}.fa" ]]; then
    echo "ERROR: Missing uncompressed genome FASTA file: $base_dir/${prefix}.fa" >&2
    exit 1
fi

mkdir -p "$base_dir/consensus_outputs"

# Step 1: Run tancons on 1ano file to obtain consensus of each tandem repeat found
echo ">>> Step 1: Running tancons..."
T_STEP1_START=$(date +%s)
tancons -o "$base_dir/consensus_outputs/${prefix}-tancons.fa" \
        -s "$base_dir/${prefix}.fa.gz" \
        "$base_dir/${prefix}.1ano" \
        > "$base_dir/consensus_outputs/${prefix}-tancons.log" 2>&1 || {
    echo "ERROR: tancons failed for $prefix! Check log: $base_dir/consensus_outputs/${prefix}-tancons.log" >&2
    exit 1
}

if [[ ! -s "$base_dir/consensus_outputs/${prefix}-tancons.fa" ]]; then
    echo "ERROR: tancons produced an empty or missing FASTA file: $base_dir/consensus_outputs/${prefix}-tancons.fa" >&2
    exit 1
fi

seqstat "$base_dir/consensus_outputs/${prefix}-tancons.fa"
T_STEP1_END=$(date +%s)
DUR_TANCONS=$((T_STEP1_END - T_STEP1_START))
echo ">>> Step 1 (tancons) completed in $(format_duration $DUR_TANCONS)"

# Step 2: Use satmatch to cluster consensus sequences from tancons
echo ">>> Step 2: Running satmatch..."
T_STEP2_START=$(date +%s)
satmatch -fa "$base_dir/consensus_outputs/${prefix}-sat_clusters.fa" \
         -loop "$base_dir/consensus_outputs/${prefix}-tancons.fa" \
         > "$base_dir/consensus_outputs/${prefix}-satmatch.log" 2>&1 || {
    echo "ERROR: satmatch failed for $prefix! Check log: $base_dir/consensus_outputs/${prefix}-satmatch.log" >&2
    exit 1
}

if [[ ! -s "$base_dir/consensus_outputs/${prefix}-sat_clusters.fa" ]]; then
    echo "ERROR: satmatch produced an empty or missing FASTA file: $base_dir/consensus_outputs/${prefix}-sat_clusters.fa" >&2
    exit 1
fi

seqstat "$base_dir/consensus_outputs/${prefix}-sat_clusters.fa"
T_STEP2_END=$(date +%s)
DUR_SATMATCH=$((T_STEP2_END - T_STEP2_START))
echo ">>> Step 2 (satmatch) completed in $(format_duration $DUR_SATMATCH)"

# Step 3: Run nhmmer to map consensus repeats across genome
echo ">>> Step 3: Running nhmmer..."
T_STEP3_START=$(date +%s)
nhmmer --cpu=8 \
       --dna \
       --qformat fasta \
       --tblout "$base_dir/consensus_outputs/${prefix}-sat_positions.tbl" \
       -o "$base_dir/consensus_outputs/${prefix}-nhmmer.out" \
       "$base_dir/consensus_outputs/${prefix}-sat_clusters.fa" \
       "$base_dir/${prefix}.fa" || {
    echo "ERROR: nhmmer mapping failed for $prefix!" >&2
    exit 1
}

if [[ ! -s "$base_dir/consensus_outputs/${prefix}-sat_positions.tbl" ]]; then
    echo "ERROR: nhmmer produced an empty or missing table file: $base_dir/consensus_outputs/${prefix}-sat_positions.tbl" >&2
    exit 1
fi

T_STEP3_END=$(date +%s)
DUR_NHMMER=$((T_STEP3_END - T_STEP3_START))
echo ">>> Step 3 (nhmmer) completed in $(format_duration $DUR_NHMMER)"

# Step 4: Generate merged satellite BED directly from nhmmer table
echo ">>> Step 4: Generating merged satellite BED..."
T_STEP4_START=$(date +%s)
awk -v OFS='\t' '!/^#/ && $1 != "" {
    start = ($7 < $8) ? $7 : $8
    end   = ($7 < $8) ? $8 : $7
    print $1, start - 1, end, $3
}' "$base_dir/consensus_outputs/${prefix}-sat_positions.tbl" \
| sort -k1,1V -k4,4 -k2,2n \
| awk -v OFS='\t' '
{
    split($4, a, "-")
    unit_size = a[length(a)] + 0
    if (unit_size <= 0) unit_size = 100

    if (NR == 1) {
        chr = $1; start = $2; end = $3; sat = $4; max_gap = unit_size
        next
    }

    # Merge if same chromosome, same satellite, and gap <= unit_size
    if ($1 == chr && $4 == sat && ($2 - end) <= max_gap) {
        if ($3 > end) end = $3
    } else {
        print chr, start, end, sat
        chr = $1; start = $2; end = $3; sat = $4; max_gap = unit_size
    }
}
END {
    if (NR > 0) print chr, start, end, sat
}' \
| sort -k1,1V -k2,2n \
> "$base_dir/consensus_outputs/${prefix}-consensus_sat.bed"

T_STEP4_END=$(date +%s)
DUR_BED=$((T_STEP4_END - T_STEP4_START))
echo ">>> Step 4 (BED merge) completed in $(format_duration $DUR_BED)"

END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))

echo "============================================"
echo " Finished at     : $(date)"
echo " Total Duration  : $(format_duration $RUNTIME)"
echo " Step Breakdown  :"
echo "   - 1. tancons  : $(format_duration $DUR_TANCONS)"
echo "   - 2. satmatch : $(format_duration $DUR_SATMATCH)"
echo "   - 3. nhmmer   : $(format_duration $DUR_NHMMER)"
echo "   - 4. BED merge: $(format_duration $DUR_BED)"
echo " Output BED      : ${base_dir}/consensus_outputs/${prefix}-consensus_sat.bed"
echo "============================================"
