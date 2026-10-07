#!/usr/bin/env bash
# Run every calibrated-grid cell one after another, resumably, ON A COMPUTE NODE.
#
#   bash tests/run_calibrated_grid.sh                      # all 18, #1..#18
#   bash tests/run_calibrated_grid.sh --cells 14           # just that one, #14
#   bash tests/run_calibrated_grid.sh --cells 17 --rerun   # redo a stuck cell
#   bash tests/run_calibrated_grid.sh --exp-base 20        # if #1..#18 are taken
#
# NOTHING RUNS ON THE LOGIN NODE. Invoked from the login node this submits
# itself with sbatch and returns immediately; the orchestrator then runs as a
# Slurm job and submits the per-stage arrays from there. That matters because
# the pipeline does not only wait between stages: run_isolated_calibration
# calls plot_and_fit_long_nsr_calibration inline and cal_2 fits inline too, so
# a login-node orchestrator is doing real compute.
#
# SIMPLICITY_GRID_LOCAL=1 runs it here in the foreground instead, for debugging.
#
# One cell at a time, each getting the whole Slurm release cap, so a cell runs
# as fast as the cluster allows and only one array is ever queued. Expect this
# to take a day or more: 18 pipelines of three dependent stages each.
#
# Resumable. --skip-completed passes over any cell that already has a
# calibration table AND its full production run, so re-submitting this exact
# command after an interruption picks up where it stopped. A cell that died
# part-way also needs --rerun, which clears its Data/*_#<n> trees first --
# without it that cell stops with "You already run an experiment with the same
# name!".
#
# Cell N is experiment #N. --cells runs a subset without renumbering the rest,
# so a trial cell keeps the number it will have in the full grid.
#
# Anything passed here goes through to tests/test_calibrated_grid.py:
#   bash tests/run_calibrated_grid.sh --populations 5000 --r-values 1.06
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1   # the repo root; running from tests/ would
                                    # scatter a stray Data/ the .gitignore misses
REPO="$PWD"

# The system python has no numpy. Prefer an explicit interpreter, then the
# active conda env, then conda's own resolution.
PY="${SIMPLICITY_PYTHON:-}"
if [ -z "$PY" ]; then
    if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/python" ]; then
        PY="$CONDA_PREFIX/bin/python"
    elif command -v conda >/dev/null 2>&1; then
        PY="$(conda info --base)/envs/simplicity/bin/python"
    else
        PY="python"
    fi
fi
if ! "$PY" -c 'import numpy' >/dev/null 2>&1; then
    echo "error: $PY cannot import numpy." >&2
    echo "       activate the simplicity env, or set SIMPLICITY_PYTHON." >&2
    exit 1
fi

mkdir -p Data/pipeline_logs

# --------------------------------------------------------------- submit path
# Not already inside Slurm, and not explicitly told to run here: hand the whole
# orchestrator to sbatch and get off the login node.
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${SIMPLICITY_GRID_LOCAL:-0}" != "1" ]; then
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "error: no sbatch here. Set SIMPLICITY_GRID_LOCAL=1 to run in the" >&2
        echo "       foreground if you really mean to." >&2
        exit 1
    fi
    STAMP="$(date +%Y%m%d_%H%M%S)"
    OUT="$REPO/Data/pipeline_logs/calibrated_grid_${STAMP}_%j.log"
    # The orchestrator waits and fits between stages rather than computing hard,
    # so one core is right; the walltime has to cover the WHOLE grid, since
    # every cell runs inside this one job.
    JOB=$(sbatch --parsable \
        --job-name=calibrated_grid \
        --time="${SIMPLICITY_GRID_TIME:-3-00:00:00}" \
        --mem="${SIMPLICITY_GRID_MEM:-8G}" \
        --cpus-per-task=1 \
        --output="$OUT" --error="$OUT" \
        --wrap "cd '$REPO' && SIMPLICITY_PYTHON='$PY' SIMPLICITY_GRID_LOCAL=1 \
                bash tests/run_calibrated_grid.sh $(printf '%q ' "$@")")
    rc=$?
    if [ "$rc" -ne 0 ] || [ -z "$JOB" ]; then
        echo "error: sbatch refused the orchestrator job (exit $rc)." >&2
        exit 1
    fi
    echo "submitted orchestrator as job $JOB"
    echo "log    : ${OUT/\%j/$JOB}"
    echo "watch  : tail -f ${OUT/\%j/$JOB}"
    echo "cancel : scancel $JOB"
    echo
    echo "Per-cell pipeline logs still go to Data/pipeline_logs/."
    exit 0
fi

# ------------------------------------------------------------------ run path
echo "interpreter : $PY"
echo "node        : $(hostname)"
echo "slurm job   : ${SLURM_JOB_ID:-none (local)}"
echo "started     : $(date -Is)"
echo

# -u so the log fills as it goes rather than at the end -- this runs for hours
# and a buffered log tells you nothing while it matters. Slurm already captures
# stdout to --output, so there is no tee here: a second copy would only drift.
"$PY" -u tests/test_calibrated_grid.py \
    --runner slurm \
    --concurrency 1 \
    --skip-completed \
    "$@"

status="$?"
echo
echo "finished    : $(date -Is)   exit $status"
if [ "$status" -ne 0 ]; then
    echo "re-submit the same command to resume; completed cells are skipped."
    echo "a cell that died part-way also needs --rerun to clear its Data trees."
fi
exit "$status"
