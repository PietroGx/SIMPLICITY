#!/usr/bin/env bash
# ============================================================================
# run_pipeline_and_grid.sh -- one overnight command: the impact_long_shedders
# production pipeline, then the studio's clade-build grid over the run it just
# produced.
#
# Runs on the LOGIN node, INSIDE tmux, from the repo root:
#
#     bash scripts/launchers/run_pipeline_and_grid.sh --settings F --check-only
#     bash scripts/launchers/run_pipeline_and_grid.sh --settings F   # check, ask, run
#     bash scripts/launchers/run_pipeline_and_grid.sh --settings F --yes
#
# F is an experiment settings file in scripts/experiments_settings/, written by
# hand or by scripts/experiment_settings_gui/. It defines the run: population,
# scenarios, windows, NSR ranges, targets, seeds, consensus and Slurm resources.
# The flags below override single values; the pipeline records what it actually
# ran with in the experiment's setup_data/. The experiment is named by its stamp,
# the minute this launcher started (GLOSSARY.md: stamp): no number to pick.
#
# stage 1  scripts/experiments_pipeline/run_impact_long_shedders_pipeline.py
#          --settings F --runner slurm: submits each stage's array and BLOCKS until done,
#          so this process is the run. Hence tmux: a login session that dies
#          takes the orchestrator with it, and nohup is not the answer here --
#          it leaves no session to reattach to and no terminal to confirm in.
# stage 2  figures_gui/prepare.py --slurm --yes for THIS run's stamp and
#          experiment: run_precompute.sh's checks, then the grid array, tx,
#          runs, explore, warm and (unless --no-compress) the afterany pack job.
#          ONLY if stage 1 exited 0. Returns without blocking. Figure 1 is not
#          built on the cluster: the studio builds it locally.
#          (The pipeline's own --prepare-figures does the same for a run started
#          without this launcher.)
#
# Every check runs before anything is submitted; one failure submits nothing.
# The two checks that also guard the grid -- the three hashed files and
# compute.code_hash() -- are spelled as figures_gui/run_precompute.sh spells
# them, not a second time in a second way. run_precompute.sh re-runs its own
# gate at stage 2, against run output that only exists by then.
#
# Memory: --slurm-mem reaches cal_1, where edge_case lives, and the experiment
# phase (run_impact_long_shedders_pipeline.py:441,460,493). cal_2 takes its
# own --slurm-mem-cal2 (cal2_memory, :84-96), 1G in the default settings: it
# profiles at 301 MB and is the cheapest stage in the pipeline.
#
# Nothing is destructive by default: the pipeline keeps each experiment's
# slurm/ logs and this launcher still packs Data_Export (it passes the
# pipeline's --pack-results, and --pack to stage 2). --drop-slurm-logs and
# --no-compress are opt-in and are named on the plan line.
#
# One log, Data/logs/pipeline/<name>_pipeline_and_grid_#<stamp>_<time>.log, with
# every stage boundary timestamped. --check-only writes no log.
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

PIPELINE="scripts/experiments_pipeline/run_impact_long_shedders_pipeline.py"
HASHED=(figures_gui/compute.py simplicity/clustering.py
        figures_gui/prep/clade_preprocess.py)

# empty = take it from the settings file
SETTINGS=""
STAMP=""
CAL_SEEDS=""
EXP_SEEDS=""
CONSENSUS=""
SLURM_MEM=""
SLURM_MEM_CAL2=""
SLURM_TIME=""
CAL_WORKERS=""
GRID_SEEDS=""
STORE="${STORE:-$ROOT/Data/figures}"
LIMIT=""
LOG=""
YES=0
CHECK_ONLY=0
NO_COMPRESS=0
DROP_SLURM_LOGS=0

usage() {
    sed -n '2,47p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'
Flags (each resolved value is printed in the checks table):
  --settings FILE       experiment settings file    [the config's default file]
  --stamp S             names the experiment, must be free [now, YYYYMMDD-HHMM]
  --cal-seeds N         seeds per calibration grid point       [settings]
  --exp-seeds N         seeds per production scenario          [settings]
  --consensus MODE      distribution | argmax, one mode        [settings]
  --slurm-mem SPEC      per-task memory, cal_1 + experiment    [settings]
  --slurm-mem-cal2 SPEC per-task memory, cal_2 only            [settings]
  --slurm-time SPEC     per-task walltime        [settings, else pipeline default]
  --cal-workers N       SIMPLICITY_CALIBRATION_WORKERS for the
                        login-node per-seed fit, ~1.6 GB each  [auto, cap 4]
  --grid-seeds N        seeds run_precompute.sh's gate demands [--exp-seeds]
  --store PATH          studio store for the grid              [<repo>/Data/figures]
  --limit N             throttle the grid array                [none]
  --log PATH            the one log                            [Data/logs/pipeline/...]
  --no-compress         skip the Data_Export packing (run data and figure data)
  --drop-slurm-logs     let the packager DELETE each experiment's slurm/ logs
  --yes, -y             skip the confirmation prompt
  --check-only          run every check, launch nothing
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --settings)         SETTINGS="$2"; shift ;;
        --stamp)          STAMP="$2"; shift ;;
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

# The settings file, through the pipeline's own loader: a file it would
# refuse is refused here, before anything else runs.
RESOLVED=$(python - "$SETTINGS" "$STAMP" <<'PY'
import os, sys
sys.path.insert(0, "scripts/experiments_pipeline")
import experiment_settings as es
import impact_long_shedders_config as cfg
s = es.load(sys.argv[1] or es.DEFAULT_SETTINGS)
blank = lambda v: "" if v is None else v
for v in (os.path.relpath(os.path.abspath(sys.argv[1] or es.DEFAULT_SETTINGS)),
          cfg.check_stamp(sys.argv[2]) if sys.argv[2] else cfg.new_stamp(),
          s["seeds"]["cal"], s["seeds"]["exp"], s["consensus"],
          s["slurm"]["mem"], blank(s["slurm"]["mem_cal2"]), blank(s["slurm"]["time"])):
    print(v)
PY
) || { echo "error: the settings file or --stamp was refused (above). Nothing launched." >&2; exit 2; }
mapfile -t S <<<"$RESOLVED"
SETTINGS="${S[0]}"
STAMP="${S[1]}"
CAL_SEEDS="${CAL_SEEDS:-${S[2]}}"
EXP_SEEDS="${EXP_SEEDS:-${S[3]}}"
CONSENSUS="${CONSENSUS:-${S[4]}}"
SLURM_MEM="${SLURM_MEM:-${S[5]}}"
SLURM_MEM_CAL2="${SLURM_MEM_CAL2:-${S[6]}}"
SLURM_TIME="${SLURM_TIME:-${S[7]}}"
# every python below, and both stages, read the config through this file
export SIMPLICITY_EXPERIMENT_SETTINGS="$ROOT/$SETTINGS"

GRID_SEEDS="${GRID_SEEDS:-$EXP_SEEDS}"
LOG_TIME="$(date +%Y%m%d_%H%M%S)"

log_line() { printf '[%s] %s\n' "$(date -Is)" "$*"; }

if [ "$CHECK_ONLY" -eq 0 ]; then
    LOG="${LOG:-Data/logs/pipeline/impact_long_shedders_pipeline_and_grid_#${STAMP}_${LOG_TIME}.log}"
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
    echo "        bash scripts/launchers/start_simplicity_session.sh"
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
sys.path.insert(0, "figures_gui")
import compute
print(compute.code_hash())
PY
)
[ -n "$HASH" ] && note "build hash" "$HASH" \
               || bad  "build hash" "compute.code_hash() printed nothing"

# --- is the stamp free? --------------------------------------------------
# The names come from the config's own helpers (run_exp_name/EXP_NAME), not
# from a second spelling here: the pipeline writes exactly these.
IDENT=$(python - "$STAMP" "$CONSENSUS" <<'PY'
import os, sys
sys.path.insert(0, "scripts/experiments_pipeline")
import simplicity.dir_manager as dm
import impact_long_shedders_config as cfg
num, mode = sys.argv[1], sys.argv[2]
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
    bad "stamp $STAMP free" "could not resolve the names to write"
    printf '%s\n' "$IDENT" | sed 's/^/      /'
elif [ -n "${ident_fail:-}" ] || [ -n "$TAKEN" ]; then
    bad "stamp $STAMP free" "TAKEN: $TAKEN"
    echo "      Launched twice this minute? Wait a minute, or pass --stamp."
else
    note "stamp $STAMP free" "will write: $WRITES"
fi

# --- the settings file ---------------------------------------------------
# Untracked is reported, not refused: the pipeline records the file it ran
# with in the run's setup_data/, which is what attributes the run.
if ! git ls-files --error-unmatch "$SETTINGS" >/dev/null 2>&1; then
    note "settings file" "$SETTINGS  (untracked; recorded in the run's setup_data/)"
elif [ -n "$(git status --porcelain -- "$SETTINGS")" ]; then
    note "settings file" "$SETTINGS  (modified; recorded in the run's setup_data/)"
else
    note "settings file" "$SETTINGS  (tracked, clean)"
fi

# --- the grid that follows ----------------------------------------------
WANT=$(python -I figures_gui/sweep.py --list | wc -l)
[ "$WANT" -gt 0 ] && note "grid size" "$WANT builds" \
                  || bad  "grid size" "sweep.py --list printed nothing"

# --- what was resolved ---------------------------------------------------
echo
echo "  resolved settings"
note "settings" "$SETTINGS"
note "experiment" "${RUN_NAME:-?}_#$STAMP  (consensus: $CONSENSUS)"
note "--cal-seeds" "$CAL_SEEDS per calibration grid point"
note "--exp-seeds" "$EXP_SEEDS per scenario, ${N_SCEN:-?} scenarios"
note "--slurm-mem (cal_1 + experiment)" \
     "$SLURM_MEM  (edge_case peaks ~1.57 GB; the default is 1G)"
note "--slurm-mem-cal2" "${SLURM_MEM_CAL2:---slurm-mem}  (cal_2 profiles at 301 MB)"
note "--slurm-time" "${SLURM_TIME:-pipeline default (SIMPLICITY_SLURM_TIME)}"
note "login-node fit workers" \
     "${CAL_WORKERS:-auto, cap 4}  (~1.6 GB each, edge_case seeds)"
note "grid experiment" "${GRID_EXP_NAME:-?} #$STAMP, gate at $GRID_SEEDS seeds"
note "grid store" "$STORE${LIMIT:+  (array throttled to $LIMIT)}"
note "packing" "$([ "$NO_COMPRESS" -eq 1 ] \
     && echo 'SKIPPED (--no-compress)' || echo 'Data_Export: run tar.xz + figure tar.gz')"
note "slurm logs" "$([ "$DROP_SLURM_LOGS" -eq 1 ] \
     && echo 'DELETED by the packager (--drop-slurm-logs)' || echo 'kept')"
note "log" "${LOG:---check-only: none}"

# --- the plan -------------------------------------------------------------
STAGE1=(python "$PIPELINE" --settings "$SETTINGS" --stamp "$STAMP"
        --runner slurm --cal-seeds "$CAL_SEEDS" --exp-seeds "$EXP_SEEDS"
        --consensus "$CONSENSUS" --slurm-mem "$SLURM_MEM")
if [ -n "$SLURM_MEM_CAL2" ]; then STAGE1+=(--slurm-mem-cal2 "$SLURM_MEM_CAL2"); fi
if [ -n "$SLURM_TIME" ]; then STAGE1+=(--slurm-time "$SLURM_TIME"); fi
if [ "$NO_COMPRESS" -eq 0 ]; then STAGE1+=(--pack-results); fi
if [ "$DROP_SLURM_LOGS" -eq 1 ]; then STAGE1+=(--drop-slurm-logs); fi

STAGE2=(python -I figures_gui/prepare.py --slurm --yes --stamp "$STAMP"
        --exp-name "$GRID_EXP_NAME" --seeds "$GRID_SEEDS" --store "$STORE")
if [ "$NO_COMPRESS" -eq 0 ]; then STAGE2+=(--pack); fi
if [ -n "$LIMIT" ]; then STAGE2+=(-- --limit "$LIMIT"); fi

echo
echo "  plan"
echo "    stage 1  ${STAGE1[*]}"
echo "             blocks: ${N_SCEN:-?} scenarios x $EXP_SEEDS seeds plus both"
echo "             calibration sweeps at $CAL_SEEDS seeds per grid point"
echo "    stage 2  STAMP=$STAMP EXP_NAME=$GRID_EXP_NAME SEEDS=$GRID_SEEDS" \
     "STORE=$STORE ${STAGE2[*]}"
echo "             only if stage 1 exits 0; $WANT builds, the transmission table,"
echo "             the per-seed export, the explore warm and the afterany pack job;"
echo "             it gates again on its own before submitting"

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
log_line "stage 1 START  ${STAGE1[*]}"
T0=$SECONDS
set +e
"${STAGE1[@]}"
RC1=$?
set -e
log_line "stage 1 END    exit $RC1 after $(( (SECONDS - T0) / 60 )) min"

if [ "$RC1" -ne 0 ]; then
    log_line "stage 2 SKIPPED -- stage 1 exited $RC1, so there is no finished run" \
          "to build the grid from. Nothing was submitted for the grid."
    echo "  log: $LOG"
    exit "$RC1"
fi

echo
log_line "stage 2 START  STAMP=$STAMP EXP_NAME=$GRID_EXP_NAME" \
      "SEEDS=$GRID_SEEDS STORE=$STORE ${STAGE2[*]}"
T1=$SECONDS
set +e
STAMP="$STAMP" EXP_NAME="$GRID_EXP_NAME" SEEDS="$GRID_SEEDS" \
    STORE="$STORE" "${STAGE2[@]}"
RC2=$?
set -e
log_line "stage 2 END    exit $RC2 after $(( (SECONDS - T1) / 60 )) min"

echo
log_line "done. stage 1 exit $RC1, stage 2 exit $RC2"
echo "  log:   $LOG"
echo "  watch: squeue -u $USER      (stage 2's array + pack job)"
exit "$RC2"
