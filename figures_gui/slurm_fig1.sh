#!/bin/bash
#SBATCH --job-name=fig2_fig1
#SBATCH --time=00:45:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=Data/logs/precompute/fig1_%j.log
#SBATCH --error=Data/logs/precompute/fig1_%j.err

# ============================================================================
# slurm_fig1.sh -- the Figure 1 data (figures_gui/fig1.py build), written into
# the experiment's store as fig1/. It reads the raw run output and the
# real-world data, not the grid, so run_precompute.sh starts it alongside the
# grid; the pack job waits for it and the archive carries fig1/.
# One process, ~5 min, ~1 GB peak, measured on a developer box.
#
#     STAMP=2 sbatch figures_gui/slurm_fig1.sh
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"
[ -f figures_gui/fig1.py ] || { echo "[fig1] no figures_gui/fig1.py in $ROOT" >&2; exit 1; }

STAMP="${STAMP:-1}"
EXP_NAME="${EXP_NAME:-}"
STORE="${STORE:-$ROOT/Data/figures}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I figures_gui/sweep.py --print-paths --store "$STORE" \
        --stamp "$STAMP" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME STAMP _REL EXP_STORE _CACHE EXP_ROOT <<<"$PATHS"

echo "fig2 fig1: $EXP_NAME #$STAMP -> $EXP_ROOT/fig1 on $(hostname -s)"
date
T0=$SECONDS
python -I figures_gui/fig1.py build --stamp "$STAMP" \
    --exp-name "$EXP_NAME" --store "$EXP_ROOT"
echo "fig2 fig1: done in $(( SECONDS - T0 ))s"
date
