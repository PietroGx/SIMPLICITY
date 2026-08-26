#!/bin/bash
#SBATCH --job-name=sot_sanity_plot
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --qos=standard

# ============================================================================
# submit_sanity_plot.sh -- ONE SLURM job producing the combined sanity
# grid (one row per long-shedder scenario) via plot_sot_sanity_regressions.py
# (heavy per-lineage hamming_iw cost belongs on SLURM, not inline after
# impact_long_shedders_exp.py). Previously one job per scenario in parallel;
# now one job runs every scenario sequentially in-process (needed to compute
# shared x/y limits across the whole grid before drawing), hence the time
# budget bumped 30min -> 2h (roughly 4x the old single-scenario budget).
#
# Usage: sbatch submit_sanity_plot.sh <exp_num> <target_osr_std> <target_osr_long>
#
# No 'set -u': conda's own activate.d hooks (e.g. referencing $CONDA_BUILD
# with no default) are not nounset-safe and would abort here before
# 'conda activate' even finishes. Same reason submit_plots.sh has no set -u.
# ============================================================================
set -eo pipefail

EXP_NUM="$1"
TARGET_OSR_STD="$2"
TARGET_OSR_LONG="$3"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

LOG_DIR="slurm_logs/sot_sanity"
mkdir -p "$LOG_DIR"

echo "Sanity plot: exp_num=$EXP_NUM (all long-shedder scenarios, combined grid)"
date

python scripts/experiments/plot_sot_sanity_regressions.py \
    --exp-num "$EXP_NUM" \
    --target-osr-std "$TARGET_OSR_STD" \
    --target-osr-long "$TARGET_OSR_LONG" \
    > "${LOG_DIR}/combined_${SLURM_JOB_ID}.log" \
    2> "${LOG_DIR}/combined_${SLURM_JOB_ID}.err"

echo "Combined sanity plot for exp_num=$EXP_NUM completed."
date
