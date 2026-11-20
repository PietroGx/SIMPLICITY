#!/bin/bash
#SBATCH --job-name=fig2_fig1
#SBATCH --time=00:45:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/fig1_%j.log
#SBATCH --error=slurm_logs/fig2_precompute/fig1_%j.err

# ============================================================================
# slurm_fig1.sh -- the Figure 1 data (tests/studio/fig1.py build), written into
# the experiment's store as fig1/. It reads the raw run output and the
# real-world data, not the grid, so run_precompute.sh starts it alongside the
# grid; the pack job waits for it and the archive carries fig1/.
# One process, ~5 min, ~1 GB peak, measured on a developer box.
#
#     EXP_NUM=2 sbatch tests/studio/slurm_fig1.sh
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT"
[ -f tests/studio/fig1.py ] || { echo "[fig1] no tests/studio/fig1.py in $ROOT" >&2; exit 1; }

EXP_NUM="${EXP_NUM:-1}"
EXP_NAME="${EXP_NAME:-}"
STORE="${STORE:-$ROOT/Data/figures_studio}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I tests/studio/sweep.py --print-paths --store "$STORE" \
        --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME EXP_NUM _REL EXP_STORE _CACHE <<<"$PATHS"

echo "fig2 fig1: $EXP_NAME #$EXP_NUM -> $EXP_STORE/fig1 on $(hostname -s)"
date
T0=$SECONDS
python -I tests/studio/fig1.py build --exp-num "$EXP_NUM" \
    --exp-name "$EXP_NAME" --store "$EXP_STORE"
echo "fig2 fig1: done in $(( SECONDS - T0 ))s"
date
