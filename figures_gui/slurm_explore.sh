#!/bin/bash
#SBATCH --job-name=fig2_explore
#SBATCH --array=0-9
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/explore_%A_%a.log
#SBATCH --error=slurm_logs/fig2_precompute/explore_%A_%a.err

# ============================================================================
# slurm_explore.sh -- figures_gui/warm_explore.py as an array: each task fills
# the daily records and hidden-time anchors of every N-th run, at the studio's
# operating point (K, MT). A one-task follow-up (--array=0-0, afterany) finds
# them all on disk and warms the Figure 5 panel at Figure 5's rules.
#
#     STAMP=2 sbatch figures_gui/slurm_explore.sh
#     STAMP=2 sbatch --dependency=afterany:<id> --array=0-0 figures_gui/slurm_explore.sh
#
# run_precompute.sh submits both alongside the grid; for a stamped experiment,
# whose grid is not built yet, scripts/launchers/run_pipeline_and_grid.sh does.
# Needs only the run output, not a clade build.
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"

STAMP="${STAMP:-latest}"
EXP_NAME="${EXP_NAME:-impact_long_shedders}"
STORE="${STORE:-$ROOT/Data/figures_studio}"
K="${K:-7}"
MT="${MT:-0.01}"
SEEDS="${SEEDS:-50}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"
N="${SLURM_ARRAY_TASK_COUNT:-1}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I figures_gui/sweep.py --print-paths --store "$STORE" \
        --stamp "$STAMP" --exp-name "$EXP_NAME")
read -r EXP_NAME STAMP _REL EXP_STORE _CACHE <<<"$PATHS"

echo "explore warm: shard $IDX/$N on $(hostname -s), $EXP_NAME #$STAMP, k=$K m_T=$MT"
date
python -I figures_gui/warm_explore.py --store "$EXP_STORE" \
    --stamp "$STAMP" --exp-name "$EXP_NAME" --k "$K" --mt "$MT" --seeds "$SEEDS" \
    --procs "${SLURM_CPUS_PER_TASK:-4}" --shard "$IDX/$N" \
    $([ "$N" -gt 1 ] && echo --no-panel)
