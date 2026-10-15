#!/usr/bin/env bash
# The HPC tier: the tests that need something only the cluster has -- a real
# Slurm controller, or the output of a real production run. Submits itself with
# sbatch and runs on a COMPUTE node. Nothing here touches the login node.
#
#   bash tests/run_tests_on_compute.sh                  # slurm interface, ~2 min
#   bash tests/run_tests_on_compute.sh --calib 2        # re-fit #2, diff it
#   bash tests/run_tests_on_compute.sh --figures 2      # figures against #2
#   bash tests/run_tests_on_compute.sh --slurm 90       # resume/reconcile, real jobs
#   bash tests/run_tests_on_compute.sh --no-interface --figures 2
#   bash tests/run_tests_on_compute.sh --with-local     # + the local tier here
#   bash tests/run_tests_on_compute.sh --no-watch       # submit and return
#
# The local tier -- serial runner, multiprocessing runner, faked Slurm, built-
# from-scratch experiments -- is NOT run here. It needs no cluster, so it runs
# on a development machine: bash tests/run_tests_local.sh. --with-local forces
# it anyway, for confirming the cluster's own python and filesystem.
#
# --slurm and the interface probe submit jobs from inside a Slurm job. That is
# the same thing run_calibrated_grid.sh does with its per-stage arrays.
set -uo pipefail

FIGURES=""
CALIB=""
CALIB_PIPELINE="bound"
SLURM_NUM=""
INTERFACE=1
WITH_LOCAL=0
WATCH=0
[ -t 1 ] && WATCH=1
while [ $# -gt 0 ]; do
    case "$1" in
        --figures)        FIGURES="$2"; shift 2 ;;
        --calib)          CALIB="$2"; shift 2 ;;
        --calib-pipeline) CALIB_PIPELINE="$2"; shift 2 ;;
        --slurm)          SLURM_NUM="$2"; shift 2 ;;
        --interface)      INTERFACE=1; shift ;;
        --no-interface)   INTERFACE=0; shift ;;
        --with-local)     WITH_LOCAL=1; shift ;;
        --watch)          WATCH=1; shift ;;
        --no-watch)       WATCH=0; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

cd "$(dirname "$0")/.." || exit 1      # the repo root; running from tests/
REPO="$PWD"                            # would scatter a stray Data/
. tests/_suite.sh

if [ "$INTERFACE" = "0" ] && [ "$WITH_LOCAL" = "0" ] \
   && [ -z "$CALIB$FIGURES$SLURM_NUM" ]; then
    echo "nothing to run: --no-interface with no --calib/--figures/--slurm." >&2
    exit 2
fi

resolve_python
mkdir -p Data/test_logs

# --------------------------------------------------------------- submit path
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${SIMPLICITY_TESTS_LOCAL:-0}" != "1" ]; then
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "error: no sbatch here, so there is no HPC tier to run. The tests" >&2
        echo "       that need no cluster are tests/run_tests_local.sh." >&2
        exit 1
    fi
    STAMP="$(date +%Y%m%d_%H%M%S)"
    OUT="$REPO/Data/test_logs/tests_${STAMP}_%j.log"
    FORWARD=""
    [ "$INTERFACE" = "0" ] && FORWARD="$FORWARD --no-interface"
    [ "$WITH_LOCAL" = "1" ] && FORWARD="$FORWARD --with-local"
    [ -n "$FIGURES" ]   && FORWARD="$FORWARD --figures $FIGURES"
    [ -n "$CALIB" ]     && FORWARD="$FORWARD --calib $CALIB --calib-pipeline $CALIB_PIPELINE"
    [ -n "$SLURM_NUM" ] && FORWARD="$FORWARD --slurm $SLURM_NUM"
    # the local tier's multiprocessing test uses a pool of 2; the HPC tier
    # itself is single-threaded and waits on Slurm
    CPUS=1
    [ "$WITH_LOCAL" = "1" ] && CPUS=2
    JOB=$(sbatch --parsable \
        --job-name=simplicity_tests \
        --time="${SIMPLICITY_TESTS_TIME:-02:00:00}" \
        --mem="${SIMPLICITY_TESTS_MEM:-8G}" \
        --cpus-per-task="$CPUS" \
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
echo "slurm job   : ${SLURM_JOB_ID:-none (forced local)}"
echo "tier        : HPC$([ "$WITH_LOCAL" = 1 ] && echo ' + local')"
echo "started     : $(date -Is)"

# does this cluster's Slurm still reply the way the parsers assume
if [ "$INTERFACE" = "1" ]; then
    run_test "slurm interface" tests/check_slurm_interface.py
fi

# re-fits a finished calibration and diffs it against what that run froze
if [ -n "$CALIB" ]; then
    run_test "calibration unchanged (#$CALIB, $CALIB_PIPELINE)" \
        tests/check_calibration_unchanged.py --exp-num "$CALIB" \
        --pipeline "$CALIB_PIPELINE"
fi

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

if [ "$WITH_LOCAL" = "1" ]; then
    run_local_tier 0
fi

summary
exit "${#FAILED[@]}"
