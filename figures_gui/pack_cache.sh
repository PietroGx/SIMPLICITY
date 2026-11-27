#!/bin/bash
#SBATCH --job-name=fig2_pack
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=Data/logs/precompute/pack_%j.log
#SBATCH --error=Data/logs/precompute/pack_%j.err

# ============================================================================
# pack_cache.sh -- one experiment's studio data, archived into
# Data_Export/figures/<name>_<stamp>/: one tar.gz per section (exploration,
# fig3..fig5; never fig1, which holds patient-derived values), split into parts of at most 1 GB that each
# unpack on their own (figures_gui/layout.py pack).
#
# `run_precompute.sh` submits this with --dependency=afterany on the array
# (through the warm array) and on the tx and runs jobs, so
# it runs even when tasks failed and SAYS which builds are missing. It used to
# be afterok, which cancelled this job on any failure: the only sign of a
# broken grid was an archive that never appeared. It can also be run by hand.
#
# It REFUSES to pack a partial grid: an archive holding 150 of 168 builds looks
# the same from the outside and would be found out later, on the wrong machine.
#
# Completeness comes from `sweep.py --missing`, which is the same `done()` test
# that makes a resubmitted array skip finished work. It is NOT a directory
# count: a build that died halfway leaves a directory with a missing or partial
# index.json, and a directory count calls that present.
#
# Exit codes, so a refusal is not just "FAILED" in sacct:
#     0  packed  (or, on a self-submit, submitted -- the pack job's own
#                 exit code is then in sacct, not here)
#     3  refused, the grid is incomplete
#     1  something else broke
#
# It packs ONE experiment's store, named by STAMP/EXP_NAME and resolved by
# `sweep.py --print-paths` -- the same resolution run_precompute.sh and
# slurm_precompute.sh use. The manifest records which experiment the runs came
# from: nothing inside a build does, because compute.params_key and index.json
# hold the clade parameters and the code hash and no experiment identity.
#
# Packing is minutes of CPU and a few hundred MB of I/O, so it belongs on a
# compute node. A hand run on a login node SUBMITS ITSELF and prints the job
# id; --here does the work in the calling shell instead.
#
# Usage, from the repo root:
#     bash   figures_gui/pack_cache.sh            # submits itself to SLURM
#     bash   figures_gui/pack_cache.sh --here     # packs in this shell
#     bash   figures_gui/pack_cache.sh --partial  # pack what exists anyway
#     bash   figures_gui/pack_cache.sh --no-warm  # skip the panel warm-up
#     STAMP=3 bash figures_gui/pack_cache.sh    # another experiment's store
#     sbatch figures_gui/pack_cache.sh            # submit explicitly
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

# sbatch copies this script to the node's spool dir, so BASH_SOURCE points
# there, not at the repo. SLURM_SUBMIT_DIR is the repo root under sbatch;
# BASH_SOURCE is the right answer when the script is run by hand.
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$ROOT"
[ -f figures_gui/compute.py ] || { echo "[pack] not a SIMPLICITY checkout: $ROOT" >&2; exit 1; }

STORE="${STORE:-$ROOT/Data/figures}"
SEEDS="${SEEDS:-50}"
STAMP="${STAMP:-1}"
EXP_NAME="${EXP_NAME:-}"
PARTIAL=0
HERE=0
WARM=1
for a in "$@"; do
    case "$a" in
        --partial) PARTIAL=1 ;;
        --here)    HERE=1 ;;
        --no-warm) WARM=0 ;;
        *) echo "[pack] unknown argument: $a" >&2; exit 1 ;;
    esac
done

# Not already inside a job, and sbatch is reachable: submit rather than pack
# here. Done before conda, which costs seconds and is only needed for the work.
if [ -z "$SLURM_JOB_ID" ] && [ "$HERE" -eq 0 ] && command -v sbatch >/dev/null 2>&1; then
    mkdir -p Data/logs/precompute
    JOB=$(STORE="$STORE" SEEDS="$SEEDS" STAMP="$STAMP" EXP_NAME="$EXP_NAME" \
          sbatch --parsable figures_gui/pack_cache.sh "$@")
    echo "[pack] submitted job $JOB"
    echo "[pack]   watch   squeue -j $JOB"
    echo "[pack]   log     Data/logs/precompute/pack_$JOB.log"
    echo "[pack]   errors  Data/logs/precompute/pack_$JOB.err"
    echo "[pack]   --here would have packed in this shell instead"
    exit 0
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

HASH=$(python - <<'PY'
import sys
sys.path.insert(0, "figures_gui")
import compute
print(compute.code_hash())
PY
)
PATHS=$(python -I figures_gui/sweep.py --print-paths --store "$STORE" \
        --stamp "$STAMP" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME STAMP EXP_REL EXP_STORE CACHE EXP_ROOT <<<"$PATHS"

WANT=$(python -I figures_gui/sweep.py --list | wc -l)
MISS=$(python -I figures_gui/sweep.py --missing \
       --out "$CACHE" --seeds "$SEEDS" \
       --stamp "$STAMP" --exp-name "$EXP_NAME")
NMISS=$(printf '%s' "$MISS" | grep -c . || true)
HAVE=$(( WANT - NMISS ))
DIRS=$(find "$CACHE" -maxdepth 1 -type d -name "*__${HASH}" 2>/dev/null | wc -l || true)

echo "[pack] build hash   $HASH"
echo "[pack] experiment   $EXP_NAME #$STAMP"
echo "[pack] cache        $CACHE"
echo "[pack] grid expects $WANT builds, $HAVE finished ($DIRS directories)"

if [ "$NMISS" -ne 0 ]; then
    IDS=$(printf '%s' "$MISS" | tr '\n' ' ')
    if [ "$PARTIAL" -eq 0 ]; then
        {
        echo "[pack] REFUSING: the grid is incomplete ($HAVE of $WANT)."
        python -I figures_gui/sweep.py --status \
            --out "$CACHE" --seeds "$SEEDS" \
            --stamp "$STAMP" --exp-name "$EXP_NAME" \
            2>&1 | sed 's/^/[pack] /' || true
        echo "[pack] $NMISS missing: $IDS"
        echo "[pack] why each task died, one file per task:"
        echo "[pack]     Data/logs/precompute/failed/<IDX>.txt"
        echo "[pack] resubmit ONLY those builds, then this packs itself:"
        echo "[pack]     bash figures_gui/run_precompute.sh --resubmit"
        if [ "$NMISS" -le 64 ]; then
            echo "[pack] or by hand:"
            echo "[pack]     STAMP=$STAMP EXP_NAME=$EXP_NAME sbatch --array=$(printf '%s' "$MISS" | paste -sd, -) figures_gui/slurm_precompute.sh $CACHE"
        fi
        echo "[pack] or pass --partial to pack anyway."
        } >&2
        exit 3
    fi
    echo "[pack] --partial: packing anyway, $NMISS missing: $IDS"
fi

# Warm the panels BEFORE packing, so the archive carries them. The grid is the
# expensive precomputation but not the whole of it: panel 4 walks every build
# and memoises a row per (build, m_T, burn-in), ~200 s per m_T value on a
# 168-build grid. Warming here, where the data already is, costs one core-hour
# of a compute node and saves that wait on every machine the archive reaches.
# Under run_precompute.sh the slurm_warm.sh array has already done it, so this
# pass reads every answer from disk and computes only what a failed shard left.
WARMED=no
if [ "$WARM" -eq 1 ]; then
    echo "[pack] warming panels before packing"
    # --stamp matters: the daily kind reads the RAW run output, and its memo
    # key does not name the experiment, so warm_panels' default of 1 would
    # read experiment 1's runs into experiment 2's store, or fail where there
    # is no experiment 1
    if python -I figures_gui/warm_panels.py --store "$EXP_STORE" \
             --stamp "$STAMP" ${EXP_NAME:+--exp-name "$EXP_NAME"} \
             --procs "${SLURM_CPUS_PER_TASK:-4}" 2>&1 | sed -u 's/^/[pack] /'; then
        WARMED=yes
    else
        echo "[pack] WARNING: the warm-up failed; packing the grid anyway" >&2
        WARMED=partial
    fi
fi

# Only builds at the current hash are packed; stale builds and caches stay
# behind (layout.py migrate --prune deletes them).
[ -f "$EXP_STORE/tx/index.json" ] \
    || echo "[pack] WARNING: no $EXP_STORE/tx -- Figures 3 and 4 will lack transmissions" >&2
[ -n "$(ls "$EXP_STORE/runs" 2>/dev/null)" ] \
    || echo "[pack] WARNING: no $EXP_STORE/runs -- Figure 2 and the Clustering tab will be empty" >&2
python -I figures_gui/layout.py pack --base "$STORE" --stamp "$STAMP" --exp-name "$EXP_NAME" \
    --note "builds=$HAVE of $WANT finished (sweep.py done(), not a directory count)" \
    --note "panels=$WARMED" \
    --note "dirty=$(git status --porcelain -- figures_gui/compute.py simplicity/clustering.py figures_gui/prep/clade_preprocess.py | wc -l) of the 3 hashed files"
