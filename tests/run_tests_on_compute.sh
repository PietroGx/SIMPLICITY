#!/usr/bin/env bash
# Run the test suite ON A COMPUTE NODE. Nothing here touches the login node.
#
#   bash tests/run_tests_on_compute.sh                  # the local suite
#   bash tests/run_tests_on_compute.sh --figures 2      # + figures against #2
#   bash tests/run_tests_on_compute.sh --slurm 90       # + resume/reconcile
#   bash tests/run_tests_on_compute.sh --figures 2 --slurm 90
#   bash tests/run_tests_on_compute.sh --no-watch       # submit and return
#
# Most of these tests build a real experiment and run real simulations -- that
# is the point of them, and it is also why they must not run where you log in.
# Invoked from the login node this submits itself with sbatch and returns; the
# tests then run as a Slurm job. SIMPLICITY_TESTS_LOCAL=1 runs them here
# instead, for a laptop.
#
# --slurm submits jobs from inside a Slurm job. That is the same thing
# run_calibrated_grid.sh does with its per-stage arrays, so it works here.
set -uo pipefail

FIGURES=""
SLURM_NUM=""
WATCH=0
[ -t 1 ] && WATCH=1
ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --figures)  FIGURES="$2"; shift 2 ;;
        --slurm)    SLURM_NUM="$2"; shift 2 ;;
        --watch)    WATCH=1; shift ;;
        --no-watch) WATCH=0; shift ;;
        *)          ARGS+=("$1"); shift ;;
    esac
done

cd "$(dirname "$0")/.." || exit 1      # the repo root; running from tests/
REPO="$PWD"                            # would scatter a stray Data/

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
    echo "error: $PY cannot import numpy. Activate the simplicity env, or set" >&2
    echo "       SIMPLICITY_PYTHON." >&2
    exit 1
fi

mkdir -p Data/test_logs

# --------------------------------------------------------------- submit path
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${SIMPLICITY_TESTS_LOCAL:-0}" != "1" ]; then
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "error: no sbatch here. Set SIMPLICITY_TESTS_LOCAL=1 to run in the" >&2
        echo "       foreground if you really mean to." >&2
        exit 1
    fi
    STAMP="$(date +%Y%m%d_%H%M%S)"
    OUT="$REPO/Data/test_logs/tests_${STAMP}_%j.log"
    FORWARD=""
    [ -n "$FIGURES" ]   && FORWARD="$FORWARD --figures $FIGURES"
    [ -n "$SLURM_NUM" ] && FORWARD="$FORWARD --slurm $SLURM_NUM"
    # 2 cpus because test_multiprocessing_runner uses a pool of 2
    JOB=$(sbatch --parsable \
        --job-name=simplicity_tests \
        --time="${SIMPLICITY_TESTS_TIME:-02:00:00}" \
        --mem="${SIMPLICITY_TESTS_MEM:-8G}" \
        --cpus-per-task=2 \
        --output="$OUT" --error="$OUT" \
        --wrap "cd '$REPO' && SIMPLICITY_PYTHON='$PY' SIMPLICITY_TESTS_LOCAL=1 \
                bash tests/run_tests_on_compute.sh$FORWARD")
    rc=$?
    if [ "$rc" -ne 0 ] || [ -z "$JOB" ]; then
        echo "error: sbatch refused the job (exit $rc)." >&2
        exit 1
    fi
    LOGFILE="${OUT/\%j/$JOB}"
    echo "submitted as job $JOB"
    echo "log    : $LOGFILE"
    echo "watch  : tail -f $LOGFILE"
    echo "cancel : scancel $JOB"
    if [ "$WATCH" = "1" ]; then
        echo
        echo "waiting for it to start (Ctrl-C stops watching, not the job)..."
        while [ ! -f "$LOGFILE" ]; do
            state=$(squeue -j "$JOB" -h -o %T 2>/dev/null)
            if [ -z "$state" ]; then
                echo "job $JOB left the queue without writing a log -- sacct -j $JOB"
                exit 1
            fi
            sleep 5
        done
        exec tail -f "$LOGFILE"
    fi
    exit 0
fi

# ------------------------------------------------------------------ run path
echo "interpreter : $PY"
echo "node        : $(hostname)"
echo "slurm job   : ${SLURM_JOB_ID:-none (local)}"
echo "started     : $(date -Is)"
echo

PASSED=(); FAILED=()

run_test () {
    local name="$1"; shift
    echo
    echo "================================================================"
    echo "  $name"
    echo "================================================================"
    local start; start=$(date +%s)
    if "$PY" -u "$@"; then
        PASSED+=("$name ($(( $(date +%s) - start ))s)")
    else
        FAILED+=("$name ($(( $(date +%s) - start ))s)")
    fi
}

# self-contained: these build a real experiment and run real simulations
run_test "settings round trip"      tests/test_settings_roundtrip.py
run_test "stale state guard"        tests/test_stale_state.py
run_test "reconcile completed"      tests/test_reconcile_completed.py
run_test "resume"                   tests/test_resume.py
run_test "slurm lifecycle (fake)"   tests/test_slurm_lifecycle.py
run_test "analysis reads layout"    tests/test_analysis_reads_new_layout.py
run_test "trees/clustering/archive" tests/test_trees_clustering_archive.py
run_test "multiprocessing runner"   tests/test_multiprocessing_runner.py

# needs real production output
if [ -n "$FIGURES" ]; then
    run_test "paper figures (#$FIGURES)" \
        tests/test_figures_on_real_output.py --exp-num "$FIGURES"
fi

# submits its own Slurm jobs, from inside this one
if [ -n "$SLURM_NUM" ]; then
    run_test "slurm resume + reconcile (#$SLURM_NUM)" \
        tests/check_slurm_resume_and_reconcile.py --exp-num "$SLURM_NUM"
fi

echo
echo "================================================================"
echo "  summary"
echo "================================================================"
for t in "${PASSED[@]}"; do echo "  PASS  $t"; done
for t in "${FAILED[@]+"${FAILED[@]}"}"; do echo "  FAIL  $t"; done
echo
echo "finished    : $(date -Is)"
echo "${#PASSED[@]} passed, ${#FAILED[@]} failed"
exit "${#FAILED[@]}"
