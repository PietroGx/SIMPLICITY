#!/usr/bin/env bash
# ============================================================================
# run_precompute.sh -- one command for the whole clade-build grid: check, ask,
# submit, pack.
#
# Runs on the LOGIN node. It submits five jobs and returns:
#
#   1. the array, one task per build (tests/studio/slurm_precompute.sh)
#   2. the transmission table panel 2 reads, alongside 1 -- it needs only the
#      run output (tests/studio/slurm_tx.sh)
#   2b. the Figure 1 data, alongside 1 -- raw run output and real-world data
#      only (tests/studio/slurm_fig1.sh); --no-fig1 skips it
#   3. the warm array, WARM_SHARDS tasks (default 32, --warm-shards N) that
#      memoise every panel answer, afterany on 1 and 2
#      (tests/studio/slurm_warm.sh)
#   4. a pack job, afterany on 3 and 2b, which tars the grid, tx/ and fig1/ into one
#      downloadable file (tests/studio/pack_cache.sh). Its own warm pass then
#      reads from disk and computes only what a failed shard left.
#      --warm-shards 0 skips 3 and leaves the warming to the pack job alone.
#
# pack_cache.sh still refuses a partial grid, so a partial grid cannot quietly
# become a file you download. The dependency is afterany, not afterok, because
# afterok CANCELLED the pack job whenever a task failed and nothing said so:
# the only symptom was an archive that never appeared. Under afterany the pack
# job runs, refuses, and prints which grid indices are missing.
#
# Before submitting it runs the checks that are cheap on a login node and
# expensive to get wrong 168 tasks in, then shows what it found and asks. The
# clade-rule check is the one that matters: if the clade definition is broken,
# every build in the grid is wrong in the same way and nothing downstream will
# say so.
#
# Usage, from the repo root:
#     bash tests/studio/run_precompute.sh                # check, ask, submit
#     bash tests/studio/run_precompute.sh --yes           # no prompt
#     bash tests/studio/run_precompute.sh --limit 24      # at most 24 at once
#     bash tests/studio/run_precompute.sh --check-only     # checks, submit nothing
#     bash tests/studio/run_precompute.sh --calibrate      # time ONE build, then stop
#     bash tests/studio/run_precompute.sh --resubmit        # only unfinished builds
#     bash tests/studio/run_precompute.sh --warm-shards 32  # warm on 32 tasks
#     bash tests/studio/run_precompute.sh --no-fig1        # no Figure 1 job
#     bash tests/studio/run_precompute.sh --dry-run        # print the sbatch chain
#     EXP_NUM=3 bash tests/studio/run_precompute.sh         # another experiment
#
# EXP_NUM/EXP_NAME (or --exp-num/--exp-name) pick the run output the builds
# read, and reach compute.py through sweep.py. They used to be printed and
# gated on but never forwarded, so EXP_NUM=3 checked experiment 3 and built
# experiment 1. Since compute.params_key holds no experiment identity, each
# experiment gets its OWN STORE -- experiment 1 of impact_long_shedders keeps
# $STORE/cache, anything else builds into $STORE/exp_<name>_<num>/cache -- and
# sweep.py refuses an --out belonging to another experiment. The resolution is
# `sweep.py --print-paths`, in one place, and the gate/build agreement is
# printed rather than assumed.
#
# --resubmit submits an array covering just the grid indices that sweep.py
# reports as not done, instead of all 168. Slurm takes an explicit index list
# but not an unbounded one, so the list is cut into arrays of at most
# RESUBMIT_CHUNK (default 100, or --chunk N) indices and the pack job depends
# on all of them with afterany. Note --limit throttles EACH array, so with 3
# chunks and --limit 24 the ceiling is 72 tasks at once, not 24.
#
# --dry-run runs the checks, then prints each sbatch with its resources and
# dependencies instead of submitting; `sbatch` is shadowed so nothing can
# reach Slurm. A failed check does not stop it.
#
# --calibrate submits a single task, waits for it, and reports what a build
# actually costs on THIS cluster. Until it has been run the cost line is an
# extrapolation from a developer box, which is the wrong machine with the wrong
# core count: it read 6 min a build at --procs 2, and a task here gets
# --cpus-per-task 4. Run it once before committing to 168.
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
CALIBRATE=0
CAL_IDX=0
RESUBMIT=0
CHUNK="${RESUBMIT_CHUNK:-100}"
TIMES="slurm_logs/fig2_precompute/timings"
EXP_NUM="${EXP_NUM:-1}"
EXP_NAME="${EXP_NAME:-}"
SEEDS="${SEEDS:-50}"
STORE="${STORE:-$ROOT/Data/figures_studio}"
WARM_SHARDS="${WARM_SHARDS:-32}"
FIG1=1
DRY_RUN=0

while [ $# -gt 0 ]; do
    case "$1" in
        --yes|-y)     YES=1 ;;
        --limit)      LIMIT="$2"; shift ;;
        --check-only) CHECK_ONLY=1 ;;
        --calibrate)  CALIBRATE=1 ;;
        --calibrate-index) CAL_IDX="$2"; shift ;;
        --resubmit)   RESUBMIT=1 ;;
        --chunk)      CHUNK="$2"; shift ;;
        --exp-num)    EXP_NUM="$2"; shift ;;
        --exp-name)   EXP_NAME="$2"; shift ;;
        --seeds)      SEEDS="$2"; shift ;;
        --warm-shards) WARM_SHARDS="$2"; shift ;;
        --no-fig1)    FIG1=0 ;;
        --dry-run)    DRY_RUN=1; YES=1 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done
case "$WARM_SHARDS" in
    ''|*[!0-9]*) echo "--warm-shards wants a whole number, got '$WARM_SHARDS'" >&2; exit 2 ;;
esac

# shadowed, not just skipped: no path below can submit under --dry-run
[ "$DRY_RUN" -eq 1 ] && sbatch() { echo "[dry-run] sbatch blocked" >&2; exit 9; }

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

# --- which experiment, and where its grid goes ---------------------------
# One implementation of the cache location, shared with slurm_precompute.sh
# and pack_cache.sh. EXP_NAME unset means sweep.py's default.
PATHS=$(python -I tests/studio/sweep.py --print-paths --store "$STORE" \
        --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME EXP_NUM EXP_REL EXP_STORE CACHE <<<"$PATHS"
note "experiment to build" "$EXP_NAME #$EXP_NUM"
note "cache directory" "$CACHE"

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
# It must be the experiment the array will BUILD, not merely the one printed:
# the gate reports the identity it resolved and the note below sets it beside
# the build's. It prints that line FIRST, so a gate that then dies on
# unreadable output still says what it was looking at.
GATE=$(python - "$EXP_NUM" "$SEEDS" "$EXP_NAME" <<'PY'
import sys, os
sys.path.insert(0, "tests/studio")
import compute
import simplicity.dir_manager as dm
exp_num, seeds, exp_name = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
print(f"gate-exp {exp_name} {exp_num}")
try:
    present, _ = compute.resolve_scenarios(exp_name, exp_num,
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
        compute.one_sod(exp_name, sc, exp_num)))
    if n < seeds:
        short.append(f"{sc}={n}")
print(f"  {'run output':54} 5/5 scenarios, "
      + ("all >= %d seeds" % seeds if not short else "SHORT: " + ", ".join(short)))
raise SystemExit(1 if short else 0)
PY
) || fail=1
printf '%s\n' "$GATE" | grep -v '^gate-exp ' || true
read -r GATE_NAME GATE_NUM \
    <<<"$(printf '%s\n' "$GATE" | sed -n 's/^gate-exp //p')" || true
if [ "$GATE_NAME" = "$EXP_NAME" ] && [ "$GATE_NUM" = "$EXP_NUM" ]; then
    note "experiment, gate == build" "$GATE_NAME #$GATE_NUM  ->  $CACHE"
else
    bad "experiment, gate vs build" \
        "gate ${GATE_NAME:-?} #${GATE_NUM:-?} != build $EXP_NAME #$EXP_NUM"
fi

# --- changes staged for the next rebuild ---------------------------------
# a pending_*.patch still applying cleanly has not been applied: building
# without it would cost a whole grid to redo
for pp in tests/studio/pending_*.patch; do
    [ -e "$pp" ] || continue
    if git apply --check -p1 "$pp" 2>/dev/null; then
        bad "staged change $(basename "$pp")" "NOT applied: git apply -p1 $pp"
    else
        note "staged change $(basename "$pp")" "applied (or no longer applies: check, then delete it)"
    fi
done

# --- the clade rules themselves -----------------------------------------
if python -I tests/check_clustering_module.py --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"} > /tmp/fig2_clade_check.$$ 2>&1; then
    note "clade rules (check_clustering_module.py)" "$(tail -1 /tmp/fig2_clade_check.$$)"
else
    bad "clade rules (check_clustering_module.py)" "$(tail -1 /tmp/fig2_clade_check.$$)"
    sed 's/^/      /' /tmp/fig2_clade_check.$$ | grep FAIL | head -5
fi
rm -f /tmp/fig2_clade_check.$$

# --- how long a build actually takes, if anything has measured it -------
RATE_MIN=6
RATE_SRC="assumed, from a developer box at --procs 2 -- run --calibrate"
if [ -d "$TIMES" ] && [ -n "$(ls -A "$TIMES" 2>/dev/null)" ]; then
    read -r RATE_MIN NTIME <<<"$(cat "$TIMES"/*.txt 2>/dev/null | awk \
        '{s+=$2; n++} END {if (n) printf "%.1f %d", s/n/60, n; else printf "6 0"}')"
    [ "${NTIME:-0}" -gt 0 ] && RATE_SRC="measured on $NTIME task(s) here"
fi
note "minutes per build" "$RATE_MIN  ($RATE_SRC)"

# --- what is left to build, and what it costs ---------------------------
# sweep.py --missing is the same done() test the array and pack_cache.sh use,
# so all three agree on what a finished build is. Captured whole, then read --
# NOT piped into `head`. `head` closes the pipe after one line, sweep.py dies
# on SIGPIPE, pipefail turns that into a failed pipeline and `set -e` exits the
# script silently, AFTER every check has printed green. It won that race on one
# box and lost it on the cluster.
MISS=$(python -I tests/studio/sweep.py --missing \
       --out "$CACHE" --seeds "$SEEDS" \
       --exp-num "$EXP_NUM" --exp-name "$EXP_NAME") || MISS="__err__"
if [ "$MISS" = "__err__" ]; then
    bad "builds still to run" "could not read sweep.py --missing"
    TODO=0
    MISS=""
else
    TODO=$(printf '%s' "$MISS" | grep -c . || true)
    note "builds still to run" "$TODO of $WANT  (a finished build is skipped)"
    # `|| true`: a new experiment's cache does not exist yet, find exits 1 and
    # pipefail would kill the script after every check had printed green.
    DIRS=$(find "$CACHE" -maxdepth 1 -type d -name "*__${HASH}" 2>/dev/null | wc -l || true)
    note "directories at this hash" "$DIRS  (not the same question as done())"
    COST=$(awk -v t="$TODO" -v r="$RATE_MIN" 'BEGIN {printf "%.0f", t*r/60}')
    note "rough cost" "${COST}h of task time, about $(( TODO * 8 / 1000 )).$(( TODO * 8 % 1000 / 100 )) GB"
fi
# --- Figure 1 -------------------------------------------------------------
if [ "$FIG1" -eq 0 ]; then
    note "Figure 1 (slurm_fig1.sh)" "skipped, --no-fig1"
elif [ ! -f tests/studio/fig1.py ]; then
    bad "Figure 1 (slurm_fig1.sh)" "tests/studio/fig1.py MISSING -- or --no-fig1"
elif ! git ls-files --error-unmatch tests/studio/fig1.py >/dev/null 2>&1; then
    note "Figure 1 (slurm_fig1.sh)" "fig1.py UNTRACKED -- a pull will not bring it"
else
    note "Figure 1 (slurm_fig1.sh)" "-> $EXP_STORE/fig1"
fi
note "archive afterwards" "~350 MB as one tar.gz (measured 8.3x on JSON)"

echo
if [ "$fail" -ne 0 ] && [ "$DRY_RUN" -eq 1 ]; then
    echo "  One or more checks FAILED. --dry-run: printing the chain anyway."
elif [ "$fail" -ne 0 ]; then
    echo "  One or more checks FAILED. Nothing submitted."
    echo "  Fix the above, or re-run with --check-only to re-check."
    exit 1
else
    echo "  All checks passed."
fi

[ "$CHECK_ONLY" -eq 1 ] && { echo "  --check-only: submitting nothing."; exit 0; }

if [ "$CALIBRATE" -eq 1 ] && [ "$DRY_RUN" -eq 1 ]; then
    echo "  --dry-run: would run: sbatch --wait --array=${CAL_IDX}-${CAL_IDX}" \
         "tests/studio/slurm_precompute.sh $CACHE"
    exit 0
fi
if [ "$CALIBRATE" -eq 1 ]; then
    mkdir -p slurm_logs/fig2_precompute
    echo
    echo "  timing ONE build (grid entry $CAL_IDX). sbatch --wait blocks until"
    echo "  it finishes, including any queue time -- expect minutes."
    rm -f "$TIMES/${CAL_IDX}.txt"
    EXP_NUM="$EXP_NUM" EXP_NAME="$EXP_NAME" SEEDS="$SEEDS" \
        sbatch --wait --array="${CAL_IDX}-${CAL_IDX}" \
        tests/studio/slurm_precompute.sh "$CACHE" || true
    if [ ! -f "$TIMES/${CAL_IDX}.txt" ]; then
        echo "  no timing was written. Either the task failed -- check"
        echo "  slurm_logs/fig2_precompute/ -- or that build already existed,"
        echo "  in which case pass --calibrate-index N for one not yet built."
        exit 1
    fi
    read -r _i SECS CPUS HOST < "$TIMES/${CAL_IDX}.txt"
    echo
    awk -v s="$SECS" -v c="$CPUS" -v h="$HOST" -v n="$WANT" 'BEGIN {
        printf "  one build took %d s (%.1f min) on %s cpus, node %s\n", s, s/60, c, h
        printf "  the whole grid is therefore about %.0f h of task time\n", n*s/3600
        printf "  at 24 tasks at once that is roughly %.1f h of wall clock\n", n*s/3600/24
    }'
    echo
    echo "  Re-run without --calibrate to submit; the cost line now uses this."
    exit 0
fi

# --- which indices go in, and in how many arrays -------------------------
# One spec per sbatch. The full grid is a range; --resubmit is an explicit list
# cut into chunks, because Slurm takes `--array=12,45,97` but not a list of any
# length. Each chunk is its own array job and the pack job waits on all of them.
SPECS=()
if [ "$RESUBMIT" -eq 1 ]; then
    n=0; cur=""
    for i in $MISS; do
        cur="${cur:+$cur,}$i"; n=$(( n + 1 ))
        [ "$n" -ge "$CHUNK" ] && { SPECS+=("$cur"); cur=""; n=0; }
    done
    [ -n "$cur" ] && SPECS+=("$cur")
    echo
    if [ "${#SPECS[@]}" -eq 0 ]; then
        echo "  --resubmit: nothing is missing. Submitting tx, Figure 1, warm and pack only."
    else
        echo "  --resubmit: $TODO missing index(es) in ${#SPECS[@]} array(s)" \
             "of at most $CHUNK"
        IDS=$(printf '%s' "$MISS" | tr '\n' ' ')
        if [ "$TODO" -le 40 ]; then
            printf '  missing:   %s\n' "$IDS"
        else
            printf '  missing:   %s...\n' "${IDS:0:200}"
            echo "             (all $TODO: sweep.py --missing --out $CACHE" \
                 "--exp-num $EXP_NUM --exp-name $EXP_NAME)"
        fi
    fi
else
    SPECS=("0-$(( WANT - 1 ))")
fi

F1TXT=""
[ "$FIG1" -eq 1 ] && F1TXT=", the Figure 1 data"
if [ "$YES" -ne 1 ]; then
    echo
    read -r -p "Submit ${#SPECS[@]} array job(s) covering $TODO builds, the transmission table$F1TXT, $WARM_SHARDS warm task(s) and the packing job? [y/N] " reply
    case "$reply" in
        y|Y|yes|YES) ;;
        *) echo "  Nothing submitted."; exit 0 ;;
    esac
fi

export STORE SEEDS EXP_NUM EXP_NAME
# sub NAME ARGS...: sbatch --parsable ARGS, printing the job id. Under --dry-run
# it prints the command and the script's #SBATCH resources and returns <NAME>.
sub() {
    local name="$1" a sh=""
    shift
    if [ "$DRY_RUN" -eq 1 ]; then
        for a in "$@"; do case "$a" in *.sh) sh="$a" ;; esac; done
        echo "  [dry-run] <$name> = sbatch --parsable $*" >&2
        echo "  [dry-run]     $(grep -oE '^#SBATCH --(time|mem|cpus-per-task)=[^ ]+' "$sh" \
             | sed 's/^#SBATCH --//' | paste -sd' ' -)" >&2
        echo "<$name>"
        return 0
    fi
    sbatch --parsable "$@"
}

if [ "$DRY_RUN" -eq 1 ]; then
    echo
    echo "  [dry-run] env STORE=$STORE SEEDS=$SEEDS EXP_NUM=$EXP_NUM EXP_NAME=$EXP_NAME"
else
    mkdir -p slurm_logs/fig2_precompute
fi

JOBS=()
echo
n=0
for spec in "${SPECS[@]}"; do
    [ -n "$LIMIT" ] && spec="${spec}%${LIMIT}"
    JOB=$(sub "array$n" --array="$spec" tests/studio/slurm_precompute.sh "$CACHE")
    n=$(( n + 1 ))
    JOBS+=("$JOB")
    SHOWN="$spec"
    [ "${#spec}" -gt 72 ] && SHOWN="${spec:0:72} ..."
    echo "  array      $JOB   ($SHOWN)"
done

TXJ=$(sub tx tests/studio/slurm_tx.sh)
echo "  tx         $TXJ   (transmission table, alongside the grid)"
JOBS+=("$TXJ")
FIG1J=""
if [ "$FIG1" -eq 1 ]; then
    # not fatal: a refused submission must not leave the grid without a pack job
    if FIG1J=$(sub fig1 tests/studio/slurm_fig1.sh); then
        echo "  fig1       $FIG1J   (Figure 1 data, alongside the grid)"
    else
        FIG1J=""
        echo "  fig1       NOT SUBMITTED -- the archive will lack fig1/" >&2
    fi
fi
DEP="afterany:$(IFS=:; echo "${JOBS[*]}")"
WARMJ=""
if [ "$WARM_SHARDS" -gt 0 ]; then
    WARMJ=$(sub warm --dependency="$DEP" --array="0-$(( WARM_SHARDS - 1 ))" \
            tests/studio/slurm_warm.sh)
    echo "  warm       $WARMJ   ($WARM_SHARDS shards, afterany on the grid and tx)"
    DEP="afterany:$WARMJ:$TXJ"
    JOBS+=("$WARMJ")
fi
[ -n "$FIG1J" ] && { DEP="$DEP:$FIG1J"; JOBS+=("$FIG1J"); }
PACK=$(sub pack --dependency="$DEP" tests/studio/pack_cache.sh)
echo "  pack       $PACK  (afterany: runs even if tasks fail, and says which)"
echo
if [ "$DRY_RUN" -eq 1 ]; then
    echo "  --dry-run: nothing submitted."
    exit 0
fi
WATCH="$PACK"
[ "${#JOBS[@]}" -gt 0 ] && WATCH="$(IFS=,; echo "${JOBS[*]}"),$PACK"
echo "  watch:     squeue -j $WATCH"
echo "  logs:      slurm_logs/fig2_precompute/"
echo "  failures:  slurm_logs/fig2_precompute/failed/<IDX>.txt"
echo "  result:    $EXP_STORE/fig2_cache_${HASH}_<date>.tar.gz + .manifest.txt"
echo "  if the grid is short the pack job exits 3 and names the indices;"
echo "  then: bash tests/studio/run_precompute.sh --resubmit"
