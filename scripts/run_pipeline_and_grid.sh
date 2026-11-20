#!/usr/bin/env bash
# ============================================================================
# run_pipeline_and_grid.sh -- one overnight command: the impact_long_shedders
# production pipeline, then the studio's clade-build grid over the run it just
# produced.
#
# Runs on the LOGIN node, INSIDE tmux, from the repo root:
#
#     bash scripts/run_pipeline_and_grid.sh --check-only   # every check, no job
#     bash scripts/run_pipeline_and_grid.sh                # check, ask, run
#     bash scripts/run_pipeline_and_grid.sh --yes           # no prompt
#
# stage 1  scripts/experiments/run_impact_long_shedders_pipeline.py
#          --runner slurm: submits each stage's array and BLOCKS until done,
#          so this process is the run. Hence tmux: a login session that dies
#          takes the orchestrator with it, and nohup is not the answer here --
#          it leaves no session to reattach to and no terminal to confirm in.
# stage 2  tests/studio/run_precompute.sh --yes, with EXP_NUM/EXP_NAME/SEEDS
#          exported so the grid is built from THIS run. ONLY if stage 1 exited
#          0. It submits the 384-task array plus the afterany pack job and
#          returns without blocking.
#
# Every check runs before anything is submitted; one failure submits nothing.
# The two checks that also guard the grid -- the three hashed files and
# compute.code_hash() -- are spelled as tests/studio/run_precompute.sh spells
# them, not a second time in a second way. run_precompute.sh re-runs its own
# gate at stage 2, against run output that only exists by then.
#
# Memory: --slurm-mem reaches cal_1, where edge_case lives, and the experiment
# phase (run_impact_long_shedders_pipeline.py:377,395,428). cal_2 takes its
# own --slurm-mem-cal2 (cal2_memory, :75-87) and stays at the 1G default: it
# profiles at 301 MB and is the cheapest stage in the pipeline.
#
# Nothing is destructive by default: the pipeline keeps each experiment's
# slurm/ logs and still packs Data_Export. --drop-slurm-logs and --no-compress
# are opt-in and are named on the plan line.
#
# One log, Data/pipeline_logs/<name>_pipeline_and_grid_#N_<stamp>.log, with
# every stage boundary timestamped. --check-only writes no log.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PIPELINE="scripts/experiments/run_impact_long_shedders_pipeline.py"
PRECOMPUTE="tests/studio/run_precompute.sh"
HASHED=(tests/studio/compute.py simplicity/clustering.py
        scripts/nature_plots/clade_preprocess.py)

EXP_NUM=2
CAL_SEEDS=30
EXP_SEEDS=50
CONSENSUS=distribution
SLURM_MEM=3G
SLURM_MEM_CAL2=1G
SLURM_TIME=""
CAL_WORKERS=""
GRID_SEEDS=""
STORE="${STORE:-$ROOT/Data/figures_studio}"
LIMIT=""
LOG=""
YES=0
CHECK_ONLY=0
NO_COMPRESS=0
DROP_SLURM_LOGS=0

usage() {
    sed -n '2,42p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'
Flags (each default is printed in the checks table):
  --exp-num N           experiment number, must be free       [2]
  --cal-seeds N         seeds per calibration grid point       [30]
  --exp-seeds N         seeds per production scenario          [50]
  --consensus MODE      distribution | argmax, one mode        [distribution]
  --slurm-mem SPEC      per-task memory, cal_1 + experiment    [3G]
  --slurm-mem-cal2 SPEC per-task memory, cal_2 only            [1G]
  --slurm-time SPEC     per-task walltime                      [pipeline default]
  --cal-workers N       SIMPLICITY_CALIBRATION_WORKERS for the
                        login-node per-seed fit, ~1.6 GB each  [auto, cap 4]
  --grid-seeds N        seeds run_precompute.sh's gate demands [--exp-seeds]
  --store PATH          studio store for the grid              [<repo>/Data/figures_studio]
  --limit N             throttle the grid array                [none]
  --log PATH            the one log                            [Data/pipeline_logs/...]
  --no-compress         skip the pipeline's Data_Export packing
  --drop-slurm-logs     let the packager DELETE each experiment's slurm/ logs
  --yes, -y             skip the confirmation prompt
  --check-only          run every check, launch nothing
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --exp-num)          EXP_NUM="$2"; shift ;;
        --cal-seeds)        CAL_SEEDS="$2"; shift ;;
        --exp-seeds)        EXP_SEEDS="$2"; shift ;;
        --consensus)        CONSENSUS="$2"; shift ;;
        --slurm-mem)        SLURM_MEM="$2"; shift ;;
        --slurm-mem-cal2)   SLURM_MEM_CAL2="$2"; shift ;;
        --slurm-time)       SLURM_TIME="$2"; shift ;;
        --cal-workers)      CAL_WORKERS="$2"; shift ;;
        --grid-seeds)       GRID_SEEDS="$2"; shift ;;
        --store)            STORE="$2"; shift ;;
        --limit)            LIMIT="$2"; shift ;;
        --log)              LOG="$2"; shift ;;
        --no-compress)      NO_COMPRESS=1 ;;
        --drop-slurm-logs)  DROP_SLURM_LOGS=1 ;;
        --yes|-y)           YES=1 ;;
        --check-only)       CHECK_ONLY=1 ;;
        -h|--help)          usage; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

GRID_SEEDS="${GRID_SEEDS:-$EXP_SEEDS}"
STAMP="$(date +%Y%m%d_%H%M%S)"

stamp() { printf '[%s] %s\n' "$(date -Is)" "$*"; }

if [ "$CHECK_ONLY" -eq 0 ]; then
    LOG="${LOG:-Data/pipeline_logs/impact_long_shedders_pipeline_and_grid_#${EXP_NUM}_${STAMP}.log}"
    mkdir -p "$(dirname "$LOG")"
    # everything below -- checks, prompt, both stages -- goes to the one log
    exec > >(tee -a "$LOG") 2>&1
fi

echo
echo "================================================================"
echo " pipeline + clade-build grid -- checks"
echo "================================================================"

fail=0
note() { printf "  %-54s %s\n" "$1" "$2"; }
bad()  { printf "  %-54s %s\n" "$1" "$2"; fail=1; }

# --- the session this must survive in -----------------------------------
if [ -n "${TMUX:-}" ]; then
    note "tmux" "$(tmux display-message -p '#S' 2>/dev/null || echo yes)"
else
    bad "tmux" "NOT in tmux (\$TMUX unset)"
    echo "      Stage 1 blocks for hours on the login node: a session that"
    echo "      dies takes it with it. nohup is not the answer -- there would"
    echo "      be nothing to reattach to. Use:"
    echo "        bash scripts/start_simplicity_session.sh"
fi

# --- the env and the tools ----------------------------------------------
if [ "${CONDA_DEFAULT_ENV:-}" = "simplicity" ]; then
    note "conda env" "simplicity"
else
    bad "conda env" "'${CONDA_DEFAULT_ENV:-none}', not simplicity"
    echo "      conda activate simplicity   (the system python has no numpy)"
fi
PY="$(command -v python || true)"
[ -n "$PY" ] && note "python" "$PY" || bad "python" "not on PATH"

SBATCH="$(command -v sbatch || true)"
[ -n "$SBATCH" ] && note "sbatch" "$SBATCH" \
                 || bad  "sbatch" "not on PATH -- this is not a login node"

# --- the tree that will produce the output ------------------------------
# Tracked modifications are refused: an overnight run whose code is not in the
# repo cannot be attributed afterwards. Untracked files are counted, not
# refused -- they are not on any read path here, and the three files whose
# content IS baked into a build are checked for tracking separately below.
HEAD_SHA="$(git rev-parse --short HEAD)"
MODIFIED=$(git status --porcelain --untracked-files=no | wc -l)
if [ "$MODIFIED" -eq 0 ]; then
    note "tracked files" "clean at $HEAD_SHA"
else
    bad "tracked files" "$MODIFIED MODIFIED -- the run would not be attributable"
    git status --short --untracked-files=no | sed 's/^/      /'
fi
UNTRACKED=$(git status --porcelain --untracked-files=normal | grep -c '^??' || true)
note "untracked files" "$UNTRACKED (not part of the run's provenance)"

if UPSTREAM="$(git rev-parse --abbrev-ref '@{u}' 2>/dev/null)"; then
    if [ "$(git rev-parse HEAD)" = "$(git rev-parse '@{u}')" ]; then
        note "HEAD == pushed commit" "$UPSTREAM at $HEAD_SHA"
    else
        bad "HEAD vs pushed commit" \
            "$HEAD_SHA != $(git rev-parse --short '@{u}') ($UPSTREAM)"
        echo "      (the comparison is against the last fetch; git pull first)"
    fi
else
    bad "HEAD vs pushed commit" "this branch has no upstream"
fi

# --- the three hashed files, spelled as run_precompute.sh spells them ----
MISSING_HASHED=()
for f in "${HASHED[@]}"; do
    git ls-files --error-unmatch "$f" >/dev/null 2>&1 || MISSING_HASHED+=("$f")
done
DIRTY=$(git status --porcelain -- "${HASHED[@]}" | wc -l)
if [ "${#MISSING_HASHED[@]}" -gt 0 ]; then
    bad "the 3 hashed files" "${#MISSING_HASHED[@]} NOT IN THE REPO"
    for f in "${MISSING_HASHED[@]}"; do echo "      untracked: $f"; done
elif [ "$DIRTY" -ne 0 ]; then
    bad "the 3 hashed files" "$DIRTY MODIFIED -- this hash is not reproducible"
    git status --short -- "${HASHED[@]}" | sed 's/^/      /'
else
    note "the 3 hashed files" "tracked and clean at $HEAD_SHA"
fi

HASH=$(python - <<'PY' || true
import sys
sys.path.insert(0, "tests/studio")
import compute
print(compute.code_hash())
PY
)
[ -n "$HASH" ] && note "build hash" "$HASH" \
               || bad  "build hash" "compute.code_hash() printed nothing"

# --- is the experiment number free? -------------------------------------
# The names come from the config's own helpers (run_exp_name/EXP_NAME), not
# from a second spelling here: the pipeline writes exactly these.
IDENT=$(python - "$EXP_NUM" "$CONSENSUS" <<'PY'
import os, sys
sys.path.insert(0, "scripts/experiments")
import simplicity.dir_manager as dm
import impact_long_shedders_config as cfg
num, mode = int(sys.argv[1]), sys.argv[2]
writes = list(dict.fromkeys([f"{cfg.EXP_NAME}_#{num}",
                             cfg.run_exp_name(num, mode)]))
taken = [n for n in writes
         if os.path.exists(os.path.join(dm.get_data_dir(), n))]
print("ident", cfg.EXP_NAME, cfg.exp_name(mode), len(cfg.SCENARIOS))
print("writes", " ".join(writes))
print("taken", " ".join(taken))
raise SystemExit(1 if taken else 0)
PY
) || ident_fail=1
read -r RUN_NAME GRID_EXP_NAME N_SCEN \
    <<<"$(printf '%s\n' "$IDENT" | sed -n 's/^ident //p')"
WRITES=$(printf '%s\n' "$IDENT" | sed -n 's/^writes //p')
TAKEN=$(printf '%s\n' "$IDENT" | sed -n 's/^taken //p')
if [ -z "$RUN_NAME" ]; then
    bad "experiment #$EXP_NUM free" "could not resolve the names to write"
    printf '%s\n' "$IDENT" | sed 's/^/      /'
elif [ -n "${ident_fail:-}" ] || [ -n "$TAKEN" ]; then
    bad "experiment #$EXP_NUM free" "TAKEN: $TAKEN"
    echo "      Pick another --exp-num, or move those directories aside."
else
    note "experiment #$EXP_NUM free" "will write: $WRITES"
fi

# --- the grid that follows ----------------------------------------------
WANT=$(python -I tests/studio/sweep.py --list | wc -l)
[ "$WANT" -gt 0 ] && note "grid size" "$WANT builds" \
                  || bad  "grid size" "sweep.py --list printed nothing"

# --- what was resolved ---------------------------------------------------
echo
echo "  resolved settings"
note "experiment" "${RUN_NAME:-?}_#$EXP_NUM  (consensus: $CONSENSUS)"
note "--cal-seeds" "$CAL_SEEDS per calibration grid point"
note "--exp-seeds" "$EXP_SEEDS per scenario, ${N_SCEN:-?} scenarios"
note "--slurm-mem (cal_1 + experiment)" \
     "$SLURM_MEM  (edge_case peaks ~1.57 GB; the default is 1G)"
note "--slurm-mem-cal2" "$SLURM_MEM_CAL2  (cal_2 profiles at 301 MB)"
note "--slurm-time" "${SLURM_TIME:-pipeline default (SIMPLICITY_SLURM_TIME)}"
note "login-node fit workers" \
     "${CAL_WORKERS:-auto, cap 4}  (~1.6 GB each, edge_case seeds)"
note "grid experiment" "${GRID_EXP_NAME:-?} #$EXP_NUM, gate at $GRID_SEEDS seeds"
note "grid store" "$STORE${LIMIT:+  (array throttled to $LIMIT)}"
note "packing" "$([ "$NO_COMPRESS" -eq 1 ] \
     && echo 'SKIPPED (--no-compress)' || echo 'Data_Export tar.xz (kept)')"
note "slurm logs" "$([ "$DROP_SLURM_LOGS" -eq 1 ] \
     && echo 'DELETED by the packager (--drop-slurm-logs)' || echo 'kept')"
note "log" "${LOG:---check-only: none}"

# --- the plan -------------------------------------------------------------
STAGE1=(python "$PIPELINE" --exp-num "$EXP_NUM" --runner slurm
        --cal-seeds "$CAL_SEEDS" --exp-seeds "$EXP_SEEDS"
        --consensus "$CONSENSUS"
        --slurm-mem "$SLURM_MEM" --slurm-mem-cal2 "$SLURM_MEM_CAL2")
if [ -n "$SLURM_TIME" ]; then STAGE1+=(--slurm-time "$SLURM_TIME"); fi
if [ "$NO_COMPRESS" -eq 1 ]; then STAGE1+=(--no-compress); fi
if [ "$DROP_SLURM_LOGS" -eq 1 ]; then STAGE1+=(--drop-slurm-logs); fi

STAGE2=(bash "$PRECOMPUTE" --yes)
if [ -n "$LIMIT" ]; then STAGE2+=(--limit "$LIMIT"); fi

echo
echo "  plan"
echo "    stage 1  ${STAGE1[*]}"
echo "             blocks: ${N_SCEN:-?} scenarios x $EXP_SEEDS seeds plus both"
echo "             calibration sweeps at $CAL_SEEDS seeds per grid point"
echo "    stage 2  EXP_NUM=$EXP_NUM EXP_NAME=$GRID_EXP_NAME SEEDS=$GRID_SEEDS" \
     "STORE=$STORE ${STAGE2[*]}"
echo "             only if stage 1 exits 0; $WANT builds + the afterany pack"
echo "             job, and it gates again on its own before submitting"

echo
if [ "$fail" -ne 0 ]; then
    echo "  One or more checks FAILED. Nothing launched, nothing submitted."
    exit 1
fi
echo "  All checks passed."

if [ "$CHECK_ONLY" -eq 1 ]; then
    echo "  --check-only: launching nothing."
    exit 0
fi

if [ "$YES" -ne 1 ]; then
    echo
    read -r -p "Run stage 1 (blocks for hours) and then stage 2? [y/N] " reply
    case "$reply" in
        y|Y|yes|YES) ;;
        *) echo "  Nothing launched."; exit 0 ;;
    esac
fi

if [ -n "$CAL_WORKERS" ]; then export SIMPLICITY_CALIBRATION_WORKERS="$CAL_WORKERS"; fi

echo
stamp "stage 1 START  ${STAGE1[*]}"
T0=$SECONDS
set +e
"${STAGE1[@]}"
RC1=$?
set -e
stamp "stage 1 END    exit $RC1 after $(( (SECONDS - T0) / 60 )) min"

if [ "$RC1" -ne 0 ]; then
    stamp "stage 2 SKIPPED -- stage 1 exited $RC1, so there is no finished run" \
          "to build the grid from. Nothing was submitted for the grid."
    echo "  log: $LOG"
    exit "$RC1"
fi

echo
stamp "stage 2 START  EXP_NUM=$EXP_NUM EXP_NAME=$GRID_EXP_NAME" \
      "SEEDS=$GRID_SEEDS STORE=$STORE ${STAGE2[*]}"
T1=$SECONDS
set +e
EXP_NUM="$EXP_NUM" EXP_NAME="$GRID_EXP_NAME" SEEDS="$GRID_SEEDS" \
    STORE="$STORE" "${STAGE2[@]}"
RC2=$?
set -e
stamp "stage 2 END    exit $RC2 after $(( (SECONDS - T1) / 60 )) min"

echo
stamp "done. stage 1 exit $RC1, stage 2 exit $RC2"
echo "  log:   $LOG"
echo "  watch: squeue -u $USER      (stage 2's array + pack job)"
exit "$RC2"
