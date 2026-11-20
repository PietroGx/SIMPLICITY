#!/bin/bash
#SBATCH --job-name=fig2_warm
#SBATCH --array=0-15
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/warm_%A_%a.log
#SBATCH --error=slurm_logs/fig2_precompute/warm_%A_%a.err

# ============================================================================
# slurm_warm.sh -- the panel warm-up (tests/studio/warm_panels.py) as an array:
# task I of N warms every N-th job from I. run_precompute.sh submits it after
# the grid and before the pack job; the pack job's own warm pass then answers
# from disk and computes only what a failed shard left.
#
#     EXP_NUM=2 sbatch tests/studio/slurm_warm.sh
#     EXP_NUM=2 sbatch --array=0-31 tests/studio/slurm_warm.sh    # 32 shards
#
# N is the array's size, so it must be one contiguous 0-(N-1) range.
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT"

EXP_NUM="${EXP_NUM:-1}"
EXP_NAME="${EXP_NAME:-}"
STORE="${STORE:-$ROOT/fig2_studio}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
N="${SLURM_ARRAY_TASK_COUNT:-1}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I tests/studio/sweep.py --print-paths --store "$STORE" \
        --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME EXP_NUM _REL EXP_STORE _CACHE <<<"$PATHS"

echo "fig2 warm: shard $IDX/$N on $(hostname -s), $EXP_NAME #$EXP_NUM"
date
python -I tests/studio/warm_panels.py --store "$EXP_STORE" \
    --exp-num "$EXP_NUM" --exp-name "$EXP_NAME" \
    --procs "${SLURM_CPUS_PER_TASK:-4}" --shard "$IDX/$N"
