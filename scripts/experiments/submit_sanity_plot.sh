#!/bin/bash
#SBATCH --job-name=sot_sanity_plot
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --qos=standard

# ============================================================================
# submit_sanity_plot.sh -- one SLURM job per long-shedder scenario, running
# plot_sot_sanity_regressions.py (heavy per-lineage hamming_iw cost belongs
# on SLURM, not inline after impact_long_shedders_exp.py).
#
# Usage: sbatch submit_sanity_plot.sh <exp_num> <scenario> <target_osr_std> <target_osr_long>
#
# No 'set -u': conda's own activate.d hooks (e.g. referencing $CONDA_BUILD
# with no default) are not nounset-safe and would abort here before
# 'conda activate' even finishes. Same reason submit_plots.sh has no set -u.
# ============================================================================
set -eo pipefail

EXP_NUM="$1"
SCENARIO="$2"
TARGET_OSR_STD="$3"
TARGET_OSR_LONG="$4"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

LOG_DIR="slurm_logs/sot_sanity"
mkdir -p "$LOG_DIR"

echo "Sanity plot: exp_num=$EXP_NUM scenario=$SCENARIO"
date

python scripts/plot_sot_sanity_regressions.py \
    --exp-num "$EXP_NUM" \
    --scenario "$SCENARIO" \
    --target-osr-std "$TARGET_OSR_STD" \
    --target-osr-long "$TARGET_OSR_LONG" \
    > "${LOG_DIR}/${SCENARIO}_${SLURM_JOB_ID}.log" \
    2> "${LOG_DIR}/${SCENARIO}_${SLURM_JOB_ID}.err"

echo "Sanity plot for scenario=$SCENARIO exp_num=$EXP_NUM completed."
date
