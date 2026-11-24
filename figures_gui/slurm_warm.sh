#!/bin/bash
#SBATCH --job-name=fig2_warm
#SBATCH --array=0-31
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/warm_%A_%a.log
#SBATCH --error=slurm_logs/fig2_precompute/warm_%A_%a.err

# ============================================================================
# slurm_warm.sh -- the panel warm-up (figures_gui/warm_panels.py) as an array:
# task I of N warms every N-th job from I. run_precompute.sh submits it after
# the grid and before the pack job; the pack job's own warm pass then answers
# from disk and computes only what a failed shard left.
#
#     STAMP=2 sbatch figures_gui/slurm_warm.sh
#     STAMP=2 sbatch --array=0-63 figures_gui/slurm_warm.sh    # 64 shards
#
# N is the array's size, so it must be one contiguous 0-(N-1) range.
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"

STAMP="${STAMP:-1}"
EXP_NAME="${EXP_NAME:-}"
STORE="${STORE:-$ROOT/Data/figures_studio}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
N="${SLURM_ARRAY_TASK_COUNT:-1}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I figures_gui/sweep.py --print-paths --store "$STORE" \
        --stamp "$STAMP" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME STAMP _REL EXP_STORE _CACHE <<<"$PATHS"

echo "fig2 warm: shard $IDX/$N on $(hostname -s), $EXP_NAME #$STAMP"
date
python -I figures_gui/warm_panels.py --store "$EXP_STORE" \
    --stamp "$STAMP" --exp-name "$EXP_NAME" \
    --procs "${SLURM_CPUS_PER_TASK:-4}" --shard "$IDX/$N"
