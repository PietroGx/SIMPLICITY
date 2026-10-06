#!/usr/bin/env bash
# Run every calibrated-grid cell one after another, resumably.
#
#   bash tests/run_calibrated_grid.sh                 # all 18 cells, #1..#18
#   bash tests/run_calibrated_grid.sh --cells 14      # just that one, as #14
#   bash tests/run_calibrated_grid.sh --exp-base 20   # if #1..#18 are taken
#
# Run it from scripts/start_simplicity_session.sh, which attaches the tmux
# session with the conda env and SBATCH_QOS=standard already set -- this takes
# hours and has to outlive the login. It logs to Data/pipeline_logs/ either way,
# so detaching loses nothing.
#
# One cell at a time, each getting the whole Slurm release cap, so a cell runs
# as fast as the cluster allows and only one array is ever queued. Expect this
# to take a day or more: 18 pipelines of three dependent stages each.
#
# Resumable. --skip-completed passes over any cell that already has a
# calibration table AND its full production run, so re-running this exact
# command after an interruption picks up where it stopped.
#
# Cell N is experiment #N. --cells runs a subset without renumbering the rest,
# so a trial cell keeps the number it will have in the full grid. If those
# numbers are already used on the cluster, pass --exp-base: a collision stops
# the stage with "You already run an experiment with the same name!" rather
# than overwriting anything.
#
# Anything passed here goes through to tests/test_calibrated_grid.py, so the
# grid can be narrowed the usual way:
#   bash tests/run_calibrated_grid.sh --populations 5000 --r-values 1.06
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1   # the repo root; running from tests/ would
                                    # scatter a stray Data/ the .gitignore misses

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
LOG="Data/pipeline_logs/calibrated_grid_sequential_$(date +%Y%m%d_%H%M%S).log"

echo "interpreter : $PY"
echo "log         : $LOG"
echo "started     : $(date -Is)"
echo

# -u so the log fills as it goes rather than at the end -- this runs for hours
# and a buffered log tells you nothing while it matters.
"$PY" -u tests/test_calibrated_grid.py \
    --runner slurm \
    --concurrency 1 \
    --skip-completed \
    "$@" 2>&1 | tee -a "$LOG"

status="${PIPESTATUS[0]}"
echo
echo "finished    : $(date -Is)   exit $status"
if [ "$status" -ne 0 ]; then
    echo "re-run the same command to resume; completed cells are skipped."
fi
exit "$status"
