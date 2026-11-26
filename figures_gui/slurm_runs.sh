#!/bin/bash
#SBATCH --job-name=fig2_runs
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=Data/logs/precompute/runs_%j.log
#SBATCH --error=Data/logs/precompute/runs_%j.err

# ============================================================================
# slurm_runs.sh -- the per-seed export runs/ (figures_gui/export_clade_explorer.py):
# the phylogeny and sparse lineage frequencies Figure 2, the Clustering tab
# and the episode grid draw. Reads only the run output, so run_precompute.sh
# starts it alongside the grid; the pack job waits for it.
#
#     STAMP=2 sbatch figures_gui/slurm_runs.sh
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

echo "fig2 runs: $EXP_NAME #$STAMP -> $EXP_STORE/runs"
python -I figures_gui/export_clade_explorer.py --stamp "$STAMP" \
    --exp-name "$EXP_NAME" --seeds "$SEEDS" \
    --procs "${SLURM_CPUS_PER_TASK:-4}" --out "$EXP_STORE/runs"
