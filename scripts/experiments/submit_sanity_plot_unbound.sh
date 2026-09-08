#!/bin/bash
#SBATCH --job-name=unbound_sanity_plot
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --qos=standard

# ============================================================================
# submit_sanity_plot_unbound.sh -- ONE SLURM job producing the combined sanity
# grid for the UNBOUND pipeline (one row per long-shedder scenario) via
# plot_sot_sanity_regressions.py.
#
# Separate from submit_sanity_plot.sh so the bound pipeline stays untouched.
# The only differences are --exp-name and the scenario list, which now
# includes edge_case; the plotting script itself is shared, so both pipelines
# measure their clocks with the same code.
#
# Time budget 3h rather than 2: four scenario rows instead of three, and
# edge_case's 350-day infections carry many more intra-host lineages through
# the per-lineage hamming_iw cost.
#
# Usage: sbatch submit_sanity_plot_unbound.sh <exp_num> <target_osr_std> <target_osr_long>
#
# No 'set -u': conda's activate.d hooks are not nounset-safe (same reason
# submit_sanity_plot.sh omits it).
# ============================================================================
set -eo pipefail

EXP_NUM="$1"
TARGET_OSR_STD="$2"
TARGET_OSR_LONG="$3"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

LOG_DIR="slurm_logs/unbound_sanity"
mkdir -p "$LOG_DIR"

echo "Unbound sanity plot: exp_num=$EXP_NUM (all long-shedder scenarios)"
date

python scripts/experiments/plot_sot_sanity_regressions.py \
    --exp-num "$EXP_NUM" \
    --exp-name impact_long_shedders_unbound \
    --scenarios SOT,HIV_low,HIV_high,edge_case \
    --target-osr-std "$TARGET_OSR_STD" \
    --target-osr-long "$TARGET_OSR_LONG" \
    > "${LOG_DIR}/combined_${SLURM_JOB_ID}.log" \
    2> "${LOG_DIR}/combined_${SLURM_JOB_ID}.err"

echo "Unbound sanity plot for exp_num=$EXP_NUM completed."
date
