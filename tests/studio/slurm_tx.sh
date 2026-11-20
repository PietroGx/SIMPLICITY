#!/bin/bash
#SBATCH --job-name=fig2_tx
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/tx_%j.log
#SBATCH --error=slurm_logs/fig2_precompute/tx_%j.err

# ============================================================================
# slurm_tx.sh -- the transmission table panel 2 reads, written into the
# experiment's store as tx/ (tests/export_transmissions.py). It reads only the
# run output, not the grid, so run_precompute.sh starts it alongside the grid;
# the warm array and the pack job wait for it, and the archive carries tx/.
#
#     EXP_NUM=2 sbatch tests/studio/slurm_tx.sh
# ============================================================================
set -eo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT"

EXP_NUM="${EXP_NUM:-1}"
EXP_NAME="${EXP_NAME:-}"
SEEDS="${SEEDS:-50}"
STORE="${STORE:-$ROOT/Data/fig2_studio}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

PATHS=$(python -I tests/studio/sweep.py --print-paths --store "$STORE" \
        --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME EXP_NUM _REL EXP_STORE _CACHE <<<"$PATHS"

echo "fig2 tx: $EXP_NAME #$EXP_NUM -> $EXP_STORE/tx"
python -I tests/export_transmissions.py --exp-num "$EXP_NUM" \
    --exp-name "$EXP_NAME" --seeds "$SEEDS" \
    --procs "${SLURM_CPUS_PER_TASK:-4}" --out "$EXP_STORE/tx"
