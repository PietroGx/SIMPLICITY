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
# Every scenario of the run's settings file, control included: the standard rate is
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
# Usage: sbatch submit_sanity_plot.sh <stamp> <target_osr_std> \
#            <target_osr_long> [exp_name]
#
# exp_name selects the consensus mode and defaults to the distribution one.
# Since 8e a scenario is a GROUP of one experiment, so the plotting script
# reads <exp_name>_#<num> and scopes by exp__<scenario>. Where the figure goes
# is NOT decided here: it comes from impact_long_shedders_config.sanity_plot_path,
# the one place that path is spelled, which nests it inside the run's own
# directory and keeps the two consensus modes apart.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

STAMP="$1"
TARGET_OSR_STD="$2"
TARGET_OSR_LONG="$3"
EXP_NAME="${4:-impact_long_shedders}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

LOG_DIR="slurm_logs/sanity"
mkdir -p "$LOG_DIR"

# the run's scenarios, from the settings file the pipeline exported
# (SIMPLICITY_EXPERIMENT_SETTINGS reaches this job through sbatch's environment)
SCENARIOS=$(python - <<'PY'
import sys
sys.path.insert(0, "scripts/experiments_pipeline")
import impact_long_shedders_config as cfg
print(",".join(s["name"] for s in cfg.SCENARIOS))
PY
)

echo "Sanity plot: stamp=$STAMP exp_name=$EXP_NAME scenarios=$SCENARIOS"
echo "settings: ${SIMPLICITY_EXPERIMENT_SETTINGS:-default file of the config}"
date

python scripts/experiments_pipeline/plot_sot_sanity_regressions.py \
    --stamp "$STAMP" \
    --exp-name "$EXP_NAME" \
    --scenarios "$SCENARIOS" \
    --target-osr-std "$TARGET_OSR_STD" \
    --target-osr-long "$TARGET_OSR_LONG" \
    > "${LOG_DIR}/${EXP_NAME}_${SLURM_JOB_ID}.log" \
    2> "${LOG_DIR}/${EXP_NAME}_${SLURM_JOB_ID}.err"

echo "Sanity plot for stamp=$STAMP exp_name=$EXP_NAME completed."
date
