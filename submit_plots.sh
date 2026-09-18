#!/bin/bash
#SBATCH --job-name=simplicity_plots
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=main
#SBATCH --qos=standard

# ============================================================================
# submit_plots.sh -- generate the paper figures for ONE pipeline arm.
#
# Usage:
#   sbatch submit_plots.sh <exp_num> [exp_name] [figure] [seed] [format]
#
#   sbatch submit_plots.sh 7                                   # bound arm, all figures
#   sbatch submit_plots.sh 1 impact_long_shedders_unbound      # unbound arm, all figures
#   sbatch submit_plots.sh 7 impact_long_shedders 3            # bound arm, figure 3 only
#
# exp_num is REQUIRED and has no default. The figure scripts dropped their
# --exp-num default of 4 in v2.4.33 because #4 was one of the runs that
# produced invalid science, and a plotting job that silently falls back to it
# is worse than one that refuses to start.
#
# Output lands in Data/figures/Figure_<n>_<arm>_#<exp_num>.<format>, so the
# two arms cannot overwrite each other.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe (same reason
# submit_sanity_plot.sh omits it).
# ============================================================================
set -eo pipefail

EXP_NUM="$1"
EXP_NAME=${2:-"impact_long_shedders"}
FIG_CHOICE=${3:-"all"}
SEED=${4:-"1"}
FORMAT=${5:-"png"}

if [[ -z "$EXP_NUM" ]]; then
    echo "Error: exp_num is required."
    echo "Usage: sbatch submit_plots.sh <exp_num> [exp_name] [figure] [seed] [format]"
    exit 1
fi

# Portable conda activation -- the old hardcoded /home/gerlep93/miniconda3
# path broke as soon as the job ran under a different HPC account. Every other
# submit script in this repo already uses this form.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

mkdir -p slurm_logs/fig1 slurm_logs/fig2 slurm_logs/fig3 slurm_logs/fig4

echo "SIMPLICITY figures: exp_num=$EXP_NUM  exp_name=$EXP_NAME"
echo "  figure=$FIG_CHOICE  seed=$SEED  format=$FORMAT"
date

run_fig() {
    local n="$1"
    shift
    echo "  -> Figure $n"
    python scripts/nature_plots/fig${n}_save.py \
        --exp-num "$EXP_NUM" \
        --exp-name "$EXP_NAME" \
        --format "$FORMAT" \
        "$@" \
        > "slurm_logs/fig${n}/plots_${SLURM_JOB_ID}.log" \
        2> "slurm_logs/fig${n}/plots_${SLURM_JOB_ID}.err"
}

# Figure 1 takes no --seed; figures 2-4 do.
case "$FIG_CHOICE" in
    1)   run_fig 1 ;;
    2)   run_fig 2 --seed "$SEED" ;;
    3)   run_fig 3 --seed "$SEED" ;;
    4)   run_fig 4 --seed "$SEED" ;;
    all) run_fig 1
         run_fig 2 --seed "$SEED"
         run_fig 3 --seed "$SEED"
         run_fig 4 --seed "$SEED" ;;
    *)   echo "Error: invalid figure '$FIG_CHOICE' (use 1, 2, 3, 4 or all)."
         exit 1 ;;
esac

echo "Job completed. Figures in Data/figures/"
date
