#!/bin/bash
#SBATCH --job-name=sanity_plot
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --qos=standard

# ============================================================================
# submit_sanity_plot.sh -- ONE SLURM job producing the combined sanity
# grid (one row per long-shedder scenario) via
# plot_sot_sanity_regressions.py.
#
# control is included: the standard rate is
# calibrated once in a clean population and the population clock is a free
# output, so control IS the baseline every elevation is measured against --
# unbound run #1 had to infer it from SOT/HIV_low agreeing. Its right-hand
# panel is empty (no long shedders); plot_sot_sanity_regressions.py handles
# that.
#
# Time budget 3h rather than 2: four scenario rows instead of three, and
# edge_case's 350-day infections carry many more intra-host lineages through
# the per-lineage hamming_iw cost.
#
# Usage: sbatch submit_sanity_plot.sh <exp_num> <target_osr_std> \
#            <target_osr_long> [exp_name]
#
# exp_name selects the consensus pipeline and defaults to the argmax one. The
# plotting script derives both its inputs (<exp_name>_<scenario>_#<num>) and its
# output directory (<exp_name>_sanity_#<num>) from it, so the two pipelines
# cannot collide.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe (same reason
# submit_sanity_plot.sh omits it).
# ============================================================================
set -eo pipefail

EXP_NUM="$1"
TARGET_OSR_STD="$2"
TARGET_OSR_LONG="$3"
EXP_NAME="${4:-impact_long_shedders}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

LOG_DIR="slurm_logs/sanity"
mkdir -p "$LOG_DIR"

echo "Sanity plot: exp_num=$EXP_NUM exp_name=$EXP_NAME (all long-shedder scenarios)"
date

python scripts/experiments/plot_sot_sanity_regressions.py \
    --exp-num "$EXP_NUM" \
    --exp-name "$EXP_NAME" \
    --scenarios control,SOT,HIV_low,HIV_high,edge_case \
    --target-osr-std "$TARGET_OSR_STD" \
    --target-osr-long "$TARGET_OSR_LONG" \
    > "${LOG_DIR}/${EXP_NAME}_${SLURM_JOB_ID}.log" \
    2> "${LOG_DIR}/${EXP_NAME}_${SLURM_JOB_ID}.err"

echo "Sanity plot for exp_num=$EXP_NUM exp_name=$EXP_NAME completed."
date
