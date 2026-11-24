#!/bin/bash
#SBATCH --job-name=fig2_tx
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=Data/logs/precompute/tx_%j.log
#SBATCH --error=Data/logs/precompute/tx_%j.err

# ============================================================================
# slurm_tx.sh -- the transmission table panel 2 reads, written into the
# experiment's store as tx/ (figures_gui/export_transmissions.py). It reads only the
# run output, not the grid, so run_precompute.sh starts it alongside the grid;
# the warm array and the pack job wait for it, and the archive carries tx/.
#
#     STAMP=2 sbatch figures_gui/slurm_tx.sh
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"

STAMP="${STAMP:-1}"
EXP_NAME="${EXP_NAME:-}"
SEEDS="${SEEDS:-50}"
STORE="${STORE:-$ROOT/Data/figures}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I figures_gui/sweep.py --print-paths --store "$STORE" \
        --stamp "$STAMP" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME STAMP _REL EXP_STORE _CACHE EXP_ROOT <<<"$PATHS"

echo "fig2 tx: $EXP_NAME #$STAMP -> $EXP_STORE/tx"
python -I figures_gui/export_transmissions.py --exp-num "$STAMP" \
    --exp-name "$EXP_NAME" --seeds "$SEEDS" \
    --procs "${SLURM_CPUS_PER_TASK:-4}" --out "$EXP_STORE/tx"
