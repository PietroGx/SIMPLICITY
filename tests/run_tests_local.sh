#!/usr/bin/env bash
# The LOCAL tier: run these on your own machine, not on the cluster.
#
#   bash tests/run_tests_local.sh            # the whole tier, ~5 min
#   bash tests/run_tests_local.sh --quick    # skip the two slow ones, ~40 s
#   bash tests/run_tests_local.sh --list     # what would run
#
# Every test here builds its own experiment and runs its own simulations --
# serial, multiprocessing, and a faked Slurm. None of them touches a real
# controller or a real run's Data/, so a compute node buys nothing: the cluster
# is for tests/run_tests_on_compute.sh.
set -uo pipefail

cd "$(dirname "$0")/.." || exit 1      # the repo root; running from tests/
. tests/_suite.sh

QUICK=0; LIST=0
while [ $# -gt 0 ]; do
    case "$1" in
        --quick) QUICK=1; shift ;;
        --list)  LIST=1; shift ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

if [ "$LIST" = "1" ]; then
    echo "local tier:"
    for entry in "${LOCAL_TESTS[@]}"; do
        name="${entry%%|*}"
        if [[ "|$LOCAL_SLOW|" == *"|$name|"* ]]; then
            echo "  $name  (slow -- dropped by --quick)"
        else
            echo "  $name"
        fi
    done
    echo
    echo "the HPC tier lives in tests/run_tests_on_compute.sh"
    exit 0
fi

# A cluster login node has sbatch and no job id. Running minutes of simulation
# there is the one thing this project never does.
if [ -z "${SLURM_JOB_ID:-}" ] && command -v sbatch >/dev/null 2>&1 \
   && [ "${SIMPLICITY_ALLOW_LOGIN_NODE:-0}" != "1" ]; then
    echo "error: this looks like a cluster login node (sbatch is on PATH and" >&2
    echo "       there is no SLURM_JOB_ID). The local tier runs real" >&2
    echo "       simulations for minutes -- not here." >&2
    echo >&2
    echo "       on your own machine : bash tests/run_tests_local.sh" >&2
    echo "       on a compute node   : bash tests/run_tests_on_compute.sh --with-local" >&2
    exit 1
fi

resolve_python
mkdir -p Data/test_logs

echo "interpreter : $PY"
echo "host        : $(hostname)"
echo "tier        : local$([ "$QUICK" = 1 ] && echo ' (--quick)')"
echo "started     : $(date -Is)"

run_local_tier "$QUICK"
summary
exit "${#FAILED[@]}"
