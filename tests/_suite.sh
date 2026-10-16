# Shared by run_tests_local.sh and run_tests_on_compute.sh. Source it from the
# repo root; it is not executable on its own.
#
# THE TWO TIERS
#
# local  self-contained. Builds its own experiments, runs its own simulations
#        through the serial and multiprocessing runners, fakes Slurm. Needs a
#        python env and nothing else. A cluster makes none of it truer, so it
#        runs on a development machine.
#
# HPC    needs something only the cluster has: a real Slurm controller, or the
#        output of a real finished run. Runs on a COMPUTE node, never where
#        you log in.
#
# A test that drifts between tiers is a test running somewhere it proves less.

LOCAL_TESTS=(
    "settings round trip|tests/test_settings_roundtrip.py"
    "settings append (8a)|tests/test_settings_append.py"
    "stale state guard|tests/test_stale_state.py"
    "reconcile completed|tests/test_reconcile_completed.py"
    "resume|tests/test_resume.py"
    "slurm lifecycle (fake)|tests/test_slurm_lifecycle.py"
    "group partition (7c)|tests/test_group_partition.py"
    "stage dispatch|tests/test_stage_dispatch.py"
    "merged experiment (8e)|tests/test_merged_experiment.py"
    "consensus modes|tests/test_consensus_pipelines.py"
    "runtime guard|tests/test_runtime_guard.py"
    "figure pipelines|tests/test_figure_pipelines.py"
    "analysis reads layout|tests/test_analysis_reads_new_layout.py"
    "trees/clustering/archive|tests/test_trees_clustering_archive.py"
    "multiprocessing runner|tests/test_multiprocessing_runner.py"
)
# the two that simulate end to end; minutes, not seconds. --quick drops them.
LOCAL_SLOW="analysis reads layout|multiprocessing runner"

PASSED=(); FAILED=()

resolve_python () {
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
        echo "error: $PY cannot import numpy. Activate the simplicity env, or" >&2
        echo "       set SIMPLICITY_PYTHON." >&2
        exit 1
    fi
}

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

run_local_tier () {
    local quick="${1:-0}" entry name path
    for entry in "${LOCAL_TESTS[@]}"; do
        name="${entry%%|*}"; path="${entry##*|}"
        if [ "$quick" = "1" ] && [[ "|$LOCAL_SLOW|" == *"|$name|"* ]]; then
            echo; echo "  SKIP  $name (--quick)"
            continue
        fi
        run_test "$name" "$path"
    done
}

summary () {
    echo
    echo "================================================================"
    echo "  summary"
    echo "================================================================"
    for t in "${PASSED[@]+"${PASSED[@]}"}"; do echo "  PASS  $t"; done
    for t in "${FAILED[@]+"${FAILED[@]}"}"; do echo "  FAIL  $t"; done
    echo
    echo "finished    : $(date -Is)"
    echo "${#PASSED[@]} passed, ${#FAILED[@]} failed"
}
