#!/usr/bin/env bash
# ============================================================================
# run_precompute.sh -- one command for the whole clade-build grid: check, ask,
# submit, pack.
#
# Runs on the LOGIN node. It submits two jobs and returns:
#
#   1. the array, one task per build (tests/studio/slurm_precompute.sh)
#   2. a pack job with --dependency=afterok on the array, which tars the grid
#      into one downloadable file (tests/studio/pack_cache.sh)
#
# The dependency is the point: the archive is only written if every build
# succeeded, so a partial grid cannot quietly become a file you download.
#
# Before submitting it runs the checks that are cheap on a login node and
# expensive to get wrong 384 tasks in, then shows what it found and asks. The
# clade-rule check is the one that matters: if the clade definition is broken,
# every build in the grid is wrong in the same way and nothing downstream will
# say so.
#
# Usage, from the repo root:
#     bash tests/studio/run_precompute.sh                # check, ask, submit
#     bash tests/studio/run_precompute.sh --yes           # no prompt
#     bash tests/studio/run_precompute.sh --limit 24      # at most 24 at once
#     bash tests/studio/run_precompute.sh --check-only     # checks, submit nothing
#     EXP_NUM=3 SEEDS=50 bash tests/studio/run_precompute.sh
#
# The env: this activates `simplicity` itself. scripts/start_simplicity_session.sh
# is the interactive way in (tmux + env + SBATCH_QOS) and cannot be called from
# here -- it ends in `exec tmux attach`, which would replace this process.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

YES=0
LIMIT=""
CHECK_ONLY=0
EXP_NUM="${EXP_NUM:-1}"
SEEDS="${SEEDS:-50}"
STORE="${STORE:-$HOME/fig2_studio}"

while [ $# -gt 0 ]; do
    case "$1" in
        --yes|-y)     YES=1 ;;
        --limit)      LIMIT="$2"; shift ;;
        --check-only) CHECK_ONLY=1 ;;
        --exp-num)    EXP_NUM="$2"; shift ;;
        --seeds)      SEEDS="$2"; shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

HASHED=(tests/studio/compute.py simplicity/clustering.py
        scripts/nature_plots/clade_preprocess.py)

if [ "${CONDA_DEFAULT_ENV:-}" != "simplicity" ]; then
    echo "[env] activating the simplicity env"
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate simplicity
fi

echo
echo "================================================================"
echo " fig2 clade-build grid -- checks"
echo "================================================================"

fail=0
note() { printf "  %-54s %s\n" "$1" "$2"; }
bad()  { printf "  %-54s %s\n" "$1" "$2"; fail=1; }

# --- the grid, and the hash the builds will carry ------------------------
WANT=$(python -I tests/studio/sweep.py --list | wc -l)
[ "$WANT" -gt 0 ] && note "grid size" "$WANT builds" \
                  || bad  "grid size" "sweep.py --list printed nothing"

HASH=$(python - <<'PY'
import sys
sys.path.insert(0, "tests/studio")
import compute
print(compute.code_hash())
PY
)
note "build hash" "$HASH"

# A build records the hash of the three files that produced it, so the studio
# can refuse a build computed by different code. Building from an uncommitted
# tree gives a hash nobody can reproduce from the repo.
# "unmodified" is not enough: an UNTRACKED hashed file is not in the repo at
# all, so the hash cannot be reproduced from a clone and the cluster cannot
# even get the file by pulling. git status says nothing about ignored files,
# which is exactly how that would pass unnoticed.
UNTRACKED=()
for f in "${HASHED[@]}"; do
    git ls-files --error-unmatch "$f" >/dev/null 2>&1 || UNTRACKED+=("$f")
done
DIRTY=$(git status --porcelain -- "${HASHED[@]}" | wc -l)
if [ "${#UNTRACKED[@]}" -gt 0 ]; then
    bad "the 3 hashed files" "${#UNTRACKED[@]} NOT IN THE REPO"
    for f in "${UNTRACKED[@]}"; do echo "      untracked: $f"; done
    echo "      A build's hash covers these files. One that is not tracked"
    echo "      cannot be pulled onto the cluster, and the hash it produces"
    echo "      cannot be reproduced from the repo."
elif [ "$DIRTY" -ne 0 ]; then
    bad "the 3 hashed files" "$DIRTY MODIFIED -- this hash is not reproducible"
    git status --short -- "${HASHED[@]}" | sed 's/^/      /'
else
    note "the 3 hashed files" "tracked and clean at $(git rev-parse --short HEAD)"
fi

# --- is the run output there, and complete? -----------------------------
python - "$EXP_NUM" "$SEEDS" <<'PY' || fail=1
import sys, os
sys.path.insert(0, "tests/studio")
import compute
import simplicity.dir_manager as dm
exp_num, seeds = int(sys.argv[1]), int(sys.argv[2])
try:
    present, _ = compute.resolve_scenarios(compute.EXP_NAME, exp_num,
                                           required=compute.SCEN_ORDER)
except Exception as exc:
    print(f"  {'run output':54} UNREADABLE: {exc}")
    raise SystemExit(1)
missing = [s for s in compute.SCEN_ORDER if s not in present]
if missing:
    print(f"  {'run output':54} MISSING {missing}")
    raise SystemExit(1)
short = []
for sc in compute.SCEN_ORDER:
    n = len(dm.get_seeded_simulation_output_dirs(
        compute.one_sod(compute.EXP_NAME, sc, exp_num)))
    if n < seeds:
        short.append(f"{sc}={n}")
print(f"  {'run output':54} 5/5 scenarios, "
      + ("all >= %d seeds" % seeds if not short else "SHORT: " + ", ".join(short)))
raise SystemExit(1 if short else 0)
PY

# --- the clade rules themselves -----------------------------------------
if python -I tests/check_clustering_module.py > /tmp/fig2_clade_check.$$ 2>&1; then
    note "clade rules (check_clustering_module.py)" "$(tail -1 /tmp/fig2_clade_check.$$)"
else
    bad "clade rules (check_clustering_module.py)" "$(tail -1 /tmp/fig2_clade_check.$$)"
    sed 's/^/      /' /tmp/fig2_clade_check.$$ | grep FAIL | head -5
fi
rm -f /tmp/fig2_clade_check.$$

# --- what is left to build, and what it costs ---------------------------
# Captured whole, then read -- NOT piped into `head`. `head` closes the pipe
# after one line, sweep.py dies on SIGPIPE, pipefail turns that into a failed
# pipeline and `set -e` exits the script silently, AFTER every check has
# printed green. It won that race on one box and lost it on the cluster.
DRY=$(python -I tests/studio/sweep.py --dry-run --seeds "$SEEDS" 2>/dev/null) || true
TODO=$(printf '%s\n' "$DRY" | sed -n '1s/.*already built, \([0-9]*\) to run.*/\1/p')
if [ -z "$TODO" ]; then
    bad "builds still to run" "could not read sweep.py --dry-run"
    TODO=0
else
    note "builds still to run" "$TODO of $WANT  (a finished build is skipped)"
    note "rough cost" "$(( TODO * 6 / 60 ))h of task time, about $(( TODO * 8 / 1000 )).$(( TODO * 8 % 1000 / 100 )) GB"
fi
note "archive afterwards" "~350 MB as one tar.gz (measured 8.3x on JSON)"

echo
if [ "$fail" -ne 0 ]; then
    echo "  One or more checks FAILED. Nothing submitted."
    echo "  Fix the above, or re-run with --check-only to re-check."
    exit 1
fi
echo "  All checks passed."

[ "$CHECK_ONLY" -eq 1 ] && { echo "  --check-only: submitting nothing."; exit 0; }

if [ "$YES" -ne 1 ]; then
    echo
    read -r -p "Submit $TODO builds plus the packing job? [y/N] " reply
    case "$reply" in
        y|Y|yes|YES) ;;
        *) echo "  Nothing submitted."; exit 0 ;;
    esac
fi

mkdir -p slurm_logs/fig2_precompute

ARRAY_SPEC="0-$(( WANT - 1 ))"
[ -n "$LIMIT" ] && ARRAY_SPEC="${ARRAY_SPEC}%${LIMIT}"

JOB=$(EXP_NUM="$EXP_NUM" SEEDS="$SEEDS" \
      sbatch --parsable --array="$ARRAY_SPEC" tests/studio/slurm_precompute.sh)
echo
echo "  array      $JOB   ($ARRAY_SPEC)"

PACK=$(STORE="$STORE" sbatch --parsable --dependency=afterok:"$JOB" \
       tests/studio/pack_cache.sh)
echo "  pack       $PACK  (runs only if every build succeeds)"
echo
echo "  watch:     squeue -j $JOB,$PACK"
echo "  logs:      slurm_logs/fig2_precompute/"
echo "  result:    $STORE/fig2_cache_${HASH}_<date>.tar.gz + .manifest.txt"
