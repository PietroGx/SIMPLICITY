#!/bin/bash
#SBATCH --job-name=fig2_precompute
#SBATCH --array=0-383
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/%A_%a.log
#SBATCH --error=slurm_logs/fig2_precompute/%A_%a.err

# ============================================================================
# slurm_precompute.sh -- the figure-2 clade-build grid, one build per array
# task.
#
# Why this is an HPC job and not a local one. Rule B designation makes m_T part
# of the clade definition: a lineage founds a clade only if its subtree carries
# at least m_T of the run's own transmission events, so moving m_T moves the
# PARTITION and cannot be applied over an existing build. The grid therefore
# gains m_T as an axis:
#
#     3 k  x  16 m_T  x  4 smoothing windows  x  2 prominences  =  384 builds
#
# At roughly 6 minutes each that is about 38 hours sequentially, which is a day
# and a half of a box that also serves production vhosts. As an array it is one
# task per build.
#
# The grid is NOT restated here. `sweep.py --list` prints it and
# `sweep.py --only N` runs entry N, so the axes live in one place and the array
# index is the contract between them. Check before submitting:
#
#     python -I tests/studio/sweep.py --list | wc -l      # must equal 384
#
# Usage, from the REPO ROOT:
#
#     mkdir -p slurm_logs/fig2_precompute
#     sbatch tests/studio/slurm_precompute.sh                      # all 384
#     sbatch --array=0-15 tests/studio/slurm_precompute.sh         # one k block
#     sbatch tests/studio/slurm_precompute.sh /scratch/$USER/cache # another out dir
#
# A task whose build already carries the current code hash and the full run
# count exits immediately, so a resubmission after a partial failure is free.
#
# What it needs: the run output for the experiment (Data/impact_long_shedders_#1
# by default) readable from the compute node, and the `simplicity` conda env.
#
# Disk, measured per build extrapolated to 250 runs: m_T=0 is 52 MB (that rung
# is divergence alone), 0.005 is 6.9 MB, 0.05 is 3.9 MB, 1.0 is 1.7 MB. The
# whole grid is about 2.9 GB -- the low rungs dominate.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

OUT="${1:-$HOME/fig2_studio/cache}"
SEEDS="${SEEDS:-50}"
EXP_NUM="${EXP_NUM:-1}"
PROCS="${SLURM_CPUS_PER_TASK:-2}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

mkdir -p "$OUT"

echo "fig2 precompute: grid entry $IDX -> $OUT"
echo "  seeds=$SEEDS exp_num=$EXP_NUM procs=$PROCS"
python -I tests/studio/sweep.py --list | sed -n "$((IDX + 1))p" \
    | awk '{print "  k="$1" m_T="$2" win="$3" prom="$4" band "$5"-"$6" "$7}'
date

python -I tests/studio/sweep.py \
    --only "$IDX" \
    --out "$OUT" \
    --procs "$PROCS" \
    --seeds "$SEEDS"

echo "fig2 precompute: grid entry $IDX done"
date
