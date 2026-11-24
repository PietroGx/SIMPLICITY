#!/bin/bash
#SBATCH --job-name=fig2_precompute
#SBATCH --array=0-167
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=Data/logs/precompute/%A_%a.log
#SBATCH --error=Data/logs/precompute/%A_%a.err

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
#     3 k  x  7 m_T  x  4 smoothing windows  x  2 prominences  =  168 builds
#
# At roughly 6 minutes each that is about 38 hours sequentially, which is a day
# and a half of a box that also serves production vhosts. As an array it is one
# task per build.
#
# The grid is NOT restated here. `sweep.py --list` prints it and
# `sweep.py --only N` runs entry N, so the axes live in one place and the array
# index is the contract between them. Check before submitting:
#
#     python -I figures_gui/sweep.py --list | wc -l      # must equal 168
#
# Usage, from the REPO ROOT:
#
#     mkdir -p Data/logs/precompute
#     sbatch figures_gui/slurm_precompute.sh                      # all 168
#     sbatch --array=0-55 figures_gui/slurm_precompute.sh         # one k block
#     sbatch figures_gui/slurm_precompute.sh /scratch/$USER/cache # another out dir
#     STAMP=3 sbatch figures_gui/slurm_precompute.sh            # another experiment
#
# A task whose build already carries the current code hash and the full run
# count exits immediately, so a resubmission after a partial failure is free.
# `run_precompute.sh --resubmit` submits an array of just the unfinished
# indices. A task that dies writes Data/logs/precompute/failed/<IDX>.txt
# naming itself, its parameters and its log, so the cause is findable without
# reading 168 logs.
#
# What it needs: the run output for the experiment (Data/impact_long_shedders_#1
# by default) readable from the compute node, and the `simplicity` conda env.
#
# STAMP/EXP_NAME are FORWARDED to sweep.py, which forwards them to
# compute.py. They used to be read and printed here and nowhere else: sweep.py
# had no --stamp, so compute.py took its own default of 1 and the line above
# was a claim about a build that read experiment 1. Because
# compute.params_key holds no experiment identity, each experiment also gets
# its own store ($STORE/exp_<name>_<num>/cache); $1 still overrides, and
# sweep.py refuses an out directory belonging to a different experiment.
#
# Disk, measured per build extrapolated to 250 runs: m_T=0 is 52 MB (that value
# is divergence alone), 0.005 is 6.9 MB, 0.05 is 3.9 MB, 1.0 is 1.7 MB. The
# whole grid is about 2.9 GB -- the low values dominate.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

# sbatch runs this from the node spool dir, so SLURM_SUBMIT_DIR is the repo
# root; BASH_SOURCE is the right answer for a hand run.
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

OUT="${1:-}"
SEEDS="${SEEDS:-50}"
STAMP="${STAMP:-1}"
EXP_NAME="${EXP_NAME:-}"
STORE="${STORE:-$ROOT/Data/figures}"
PROCS="${SLURM_CPUS_PER_TASK:-2}"
IDX="${SLURM_ARRAY_TASK_ID:-0}"

LOGS="Data/logs/precompute"
TIMES="$LOGS/timings"
FAILED="$LOGS/failed"
STEM="${SLURM_ARRAY_JOB_ID:-local}_${IDX}"
mkdir -p "$TIMES" "$FAILED"

# One file per task, never a shared append: 168 tasks appending to one file is
# a race. A task that dies leaves $FAILED/<IDX>.txt naming itself and its log,
# so the cause is found without reading 168 logs. A rerun clears it first.
rm -f "$FAILED/${IDX}.txt"
trap 'rc=$?; [ "$rc" -eq 0 ] || {
        echo "idx     $IDX"
        echo "rc      $rc"
        echo "job     $STEM"
        echo "host    $(hostname -s)"
        echo "when    $(date -Is)"
        echo "params  ${PARAMS:-unread}"
        echo "out     $OUT"
        st=$(python -I figures_gui/sweep.py --status --out "$OUT" \
             --seeds "$SEEDS" --stamp "$STAMP" \
             ${EXP_NAME:+--exp-name "$EXP_NAME"} 2>/dev/null \
             | awk -v i="$IDX" "\$1==i {print \$2}" || true)
        echo "state   ${st:-unknown}"
        echo "log     $LOGS/${STEM}.log"
        echo "err     $LOGS/${STEM}.err"
      } > "$FAILED/${IDX}.txt"' EXIT

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

# The out directory, resolved the same way run_precompute.sh and pack_cache.sh
# resolve it -- `sweep.py --print-paths`, needing the env above. $1 wins.
PATHS=$(python -I figures_gui/sweep.py --print-paths --store "$STORE" \
        --stamp "$STAMP" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME STAMP _REL _STORE CACHE _ROOT <<<"$PATHS"
OUT="${OUT:-$CACHE}"

mkdir -p "$OUT"

PARAMS=$(python -I figures_gui/sweep.py --list | sed -n "$((IDX + 1))p" \
    | awk '{print "k="$1" m_T="$2" win="$3" prom="$4" band "$5"-"$6" "$7}')

echo "fig2 precompute: grid entry $IDX -> $OUT"
echo "  seeds=$SEEDS procs=$PROCS"
echo "  experiment $EXP_NAME #$STAMP  (forwarded to compute.py)"
echo "  $PARAMS"
date

T0=$SECONDS
python -I figures_gui/sweep.py \
    --only "$IDX" \
    --out "$OUT" \
    --procs "$PROCS" \
    --seeds "$SEEDS" \
    --stamp "$STAMP" \
    --exp-name "$EXP_NAME"
ELAPSED=$(( SECONDS - T0 ))

# The point of this record is to replace a guessed rate with a measured one.
# run_precompute.sh averages whatever files exist.
printf '%s %s %s %s\n' "$IDX" "$ELAPSED" "$PROCS" "$(hostname -s)" \
    > "$TIMES/${IDX}.txt"

echo "fig2 precompute: grid entry $IDX done in ${ELAPSED}s on $PROCS cpus"
date
