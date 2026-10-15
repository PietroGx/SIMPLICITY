#!/usr/bin/env bash
# Every NSR range probe a run needs, in one job: the three cal_1 cohorts and
# the single cal_2 sweep.
#
#   bash tests/run_nsr_range_finder.sh                  # all four
#   bash tests/run_nsr_range_finder.sh --base 10        # number them #10..#13
#   bash tests/run_nsr_range_finder.sh --only stage1
#   bash tests/run_nsr_range_finder.sh --steps 8 --seeds 20
#   bash tests/run_nsr_range_finder.sh --no-watch       # submit and return
#   bash tests/run_nsr_range_finder.sh --dry-run        # print the probes only
#
# Anything after `--` goes to test_nsr_range_finder.py unchanged, e.g. a
# cheap shakedown:  --runner serial --steps 2 --seeds 1 -- --population-size 60
#
# WHY THESE FOUR
#
# cal_1 probes ONE (tau_3_long, R_long) cohort per run and there are three.
# They sweep the SAME grid -- NSR_RANGES['cal_long'] is one range for all of
# them -- over DIFFERENT windows, so one probe cannot speak for the others:
#
#                        tau_3_long   infection   cal_1 sim window
#     SOT                     48.23       63 d             365 d
#     HIV_low/HIV_high        94.23      109 d             365 d
#     edge_case              350.23      365 d            1095 d
#
# Stage 1 measures the INTRA-HOST clock, so its saturation window is one
# infection (tau_3_long + phase_offset), not the simulated window: d =
# mean_osr * infection / 365.25. edge_case's infection is 5.8x SOT's, so it
# reaches the d>0.1 ceiling at roughly a SIXTH of SOT's NSR -- probing only
# SOT would recommend a max already saturated for edge_case.
#
# HIV_low and HIV_high are ONE cohort (same tau, same R_long), so only
# HIV_low is probed. `control` has no long side, so cal_1 has no group for it.
#
# cal_2 is the 0%-long-shedder sweep: one group, one standard rate shared by
# every scenario, nothing scenario-specific to probe.
#
# The sweep bounds are READ FROM THE CONFIG, not typed here, so each probe
# spans exactly the range a run would currently use and the recap answers
# "is NSR_RANGES right?" rather than "is some range right?".
#
# Each probe BLOCKS -- the Slurm runner polls its own array to completion --
# so they run one after another inside a single job. This script submits
# itself with sbatch; nothing runs on the login node. A failing probe does not
# stop the rest.
set -uo pipefail

BASE=1
STEPS=6
SEEDS=10
ONLY=""
RUNNER="slurm"
DRY=0
EXTRA=()
WATCH=0
[ -t 1 ] && WATCH=1
while [ $# -gt 0 ]; do
    case "$1" in
        --base)     BASE="$2"; shift 2 ;;
        --steps)    STEPS="$2"; shift 2 ;;
        --seeds)    SEEDS="$2"; shift 2 ;;
        --only)     ONLY="$2"; shift 2 ;;
        --runner)   RUNNER="$2"; shift 2 ;;
        --dry-run)  DRY=1; shift ;;
        --watch)    WATCH=1; shift ;;
        --no-watch) WATCH=0; shift ;;
        --)         shift; EXTRA=("$@"); break ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
case "$ONLY" in ""|stage1|stage2) ;;
    *) echo "--only takes stage1 or stage2, got '$ONLY'" >&2; exit 2 ;;
esac

cd "$(dirname "$0")/.." || exit 1      # the repo root; running from tests/
REPO="$PWD"                            # would scatter a stray Data/
. tests/_suite.sh
resolve_python
mkdir -p Data/test_logs

# --------------------------------------------------------------- submit path
if [ -z "${SLURM_JOB_ID:-}" ] && [ "$RUNNER" = "slurm" ] && [ "$DRY" = "0" ]; then
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "error: no sbatch here. This probe needs a real Slurm controller;" >&2
        echo "       pass --runner serial to run it on this machine instead." >&2
        exit 1
    fi
    STAMP="$(date +%Y%m%d_%H%M%S)"
    OUT="$REPO/Data/test_logs/nsr_range_${STAMP}_%j.log"
    FORWARD="--base $BASE --steps $STEPS --seeds $SEEDS"
    [ -n "$ONLY" ] && FORWARD="$FORWARD --only $ONLY"
    [ "${#EXTRA[@]}" -gt 0 ] && FORWARD="$FORWARD -- ${EXTRA[*]}"
    JOB=$(sbatch --parsable \
        --job-name=nsr_range_finder \
        --time="${SIMPLICITY_RANGE_TIME:-12:00:00}" \
        --mem="${SIMPLICITY_RANGE_MEM:-8G}" \
        --cpus-per-task=1 \
        --output="$OUT" --error="$OUT" \
        --wrap "cd '$REPO' && SIMPLICITY_PYTHON='$PY' \
                bash tests/run_nsr_range_finder.sh $FORWARD")
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
        # An empty squeue reply is not proof the job died: squeue fails
        # transiently, and the shared filesystem can lag the log a compute
        # node just wrote. Four tries before believing it, then say what
        # sacct says.
        gone=0
        while [ ! -f "$LOGFILE" ]; do
            if [ -n "$(squeue -j "$JOB" -h -o %T 2>/dev/null)" ]; then
                gone=0
            else
                gone=$((gone + 1))
                if [ "$gone" -ge 4 ]; then
                    echo
                    echo "job $JOB is no longer queued and no log appeared at"
                    echo "  $LOGFILE"
                    command -v sacct >/dev/null 2>&1 && sacct -j "$JOB" \
                        --format=JobID,JobName%20,State,ExitCode,Reason%30,Start,End
                    exit 1
                fi
            fi
            sleep 5
        done
        exec tail -f "$LOGFILE"
    fi
    exit 0
fi

# ------------------------------------------------------------------ run path
# The bounds a run would actually use, straight from the config.
BOUNDS="$("$PY" - <<'PY'
import sys
sys.path[:0] = ['.', 'scripts', 'scripts/experiments']
from impact_long_shedders_config import NSR_RANGES as R
print(R['cal_long']['min'], R['cal_long']['max'],
      R['cal_std']['min'], R['cal_std']['max'])
PY
)" || { echo "error: could not read NSR_RANGES from the config." >&2; exit 1; }
read -r S1_MIN S1_MAX S2_MIN S2_MAX <<<"$BOUNDS"

# stage|scenario|min|max   -- scenario is unused by stage 2 and not passed
PROBES=(
    "1|SOT|$S1_MIN|$S1_MAX"
    "1|HIV_low|$S1_MIN|$S1_MAX"
    "1|edge_case|$S1_MIN|$S1_MAX"
    "2||$S2_MIN|$S2_MAX"
)

echo "interpreter : $PY"
echo "node        : $(hostname)"
echo "slurm job   : ${SLURM_JOB_ID:-none}"
echo "runner      : $RUNNER"
echo "grid        : $STEPS nodes, $SEEDS seeds/node"
echo "cal_long    : $S1_MIN .. $S1_MAX   (NSR_RANGES['cal_long'])"
echo "cal_std     : $S2_MIN .. $S2_MAX   (NSR_RANGES['cal_std'])"
echo "started     : $(date -Is)"

OK=(); BAD=(); RECAPS=()
N=$BASE
for probe in "${PROBES[@]}"; do
    IFS='|' read -r stage scenario lo hi <<<"$probe"
    [ "$ONLY" = "stage1" ] && [ "$stage" != "1" ] && continue
    [ "$ONLY" = "stage2" ] && [ "$stage" != "2" ] && continue

    label="stage $stage${scenario:+ / $scenario}"
    args=(tests/test_nsr_range_finder.py "$N" --stage "$stage"
          --runner "$RUNNER" --nsr-min "$lo" --nsr-max "$hi"
          --steps "$STEPS" --n-seeds "$SEEDS")
    [ -n "$scenario" ] && args+=(--scenario "$scenario")
    [ "${#EXTRA[@]}" -gt 0 ] && args+=("${EXTRA[@]}")

    echo
    echo "================================================================"
    echo "  #$N  $label   $(date -Is)"
    echo "================================================================"
    if [ "$DRY" = "1" ]; then
        echo "  would run: $PY ${args[*]}"
        OK+=("#$N  $label  (dry run)")
        RECAPS+=("Data/tests/test_nsr_range_out_#$N.txt")
        N=$((N + 1))
        continue
    fi
    start=$SECONDS
    if "$PY" "${args[@]}"; then
        OK+=("#$N  $label  ($((SECONDS - start))s)")
    else
        BAD+=("#$N  $label  (exit $?, $((SECONDS - start))s)")
    fi
    RECAPS+=("Data/tests/test_nsr_range_out_#$N.txt")
    N=$((N + 1))
done

echo
echo "================================================================"
echo "  summary"
echo "================================================================"
for line in "${OK[@]}";  do echo "  PASS  $line"; done
for line in "${BAD[@]}"; do echo "  FAIL  $line"; done
echo
echo "recaps:"
for f in "${RECAPS[@]}"; do
    [ -f "$f" ] && echo "  $f" || echo "  $f   (missing)"
done
echo
echo "finished    : $(date -Is)"
echo "${#OK[@]} ok, ${#BAD[@]} failed"

# The point of the run: every recommendation side by side.
[ "$DRY" = "1" ] && exit 0
for f in "${RECAPS[@]}"; do
    [ -f "$f" ] || continue
    echo
    echo "---------------------------------------------------------------- $f"
    cat "$f"
done

[ "${#BAD[@]}" -eq 0 ]
