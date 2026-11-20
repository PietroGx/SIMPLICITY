#!/bin/bash
#SBATCH --job-name=fig2_pack
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=main
#SBATCH --qos=standard
#SBATCH --output=slurm_logs/fig2_precompute/pack_%j.log
#SBATCH --error=slurm_logs/fig2_precompute/pack_%j.err

# ============================================================================
# pack_cache.sh -- the clade-build grid as ONE downloadable file.
#
# `run_precompute.sh` submits this with --dependency=afterany on the array, so
# it runs even when tasks failed and SAYS which builds are missing. It used to
# be afterok, which cancelled this job on any failure: the only sign of a
# broken grid was an archive that never appeared. It can also be run by hand.
#
# It REFUSES to pack a partial grid. The whole point of one file is that what
# arrives is the whole grid; a tarball holding 150 of 168 builds looks the same
# from the outside and would be found out later, on the wrong machine.
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
# It packs ONE experiment's store, named by EXP_NUM/EXP_NAME and resolved by
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
#     bash   tests/studio/pack_cache.sh            # submits itself to SLURM
#     bash   tests/studio/pack_cache.sh --here     # packs in this shell
#     bash   tests/studio/pack_cache.sh --partial  # pack what exists anyway
#     bash   tests/studio/pack_cache.sh --no-warm  # skip the panel warm-up
#     EXP_NUM=3 bash tests/studio/pack_cache.sh    # another experiment's store
#     sbatch tests/studio/pack_cache.sh            # submit explicitly
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

# sbatch copies this script to the node's spool dir, so BASH_SOURCE points
# there, not at the repo. SLURM_SUBMIT_DIR is the repo root under sbatch;
# BASH_SOURCE is the right answer when the script is run by hand.
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "$ROOT"
[ -f tests/studio/compute.py ] || { echo "[pack] not a SIMPLICITY checkout: $ROOT" >&2; exit 1; }

STORE="${STORE:-$ROOT/Data/fig2_studio}"
SEEDS="${SEEDS:-50}"
EXP_NUM="${EXP_NUM:-1}"
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
    mkdir -p slurm_logs/fig2_precompute
    JOB=$(STORE="$STORE" SEEDS="$SEEDS" EXP_NUM="$EXP_NUM" EXP_NAME="$EXP_NAME" \
          sbatch --parsable tests/studio/pack_cache.sh "$@")
    echo "[pack] submitted job $JOB"
    echo "[pack]   watch   squeue -j $JOB"
    echo "[pack]   log     slurm_logs/fig2_precompute/pack_$JOB.log"
    echo "[pack]   errors  slurm_logs/fig2_precompute/pack_$JOB.err"
    echo "[pack]   --here would have packed in this shell instead"
    exit 0
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate simplicity

HASH=$(python - <<'PY'
import sys
sys.path.insert(0, "tests/studio")
import compute
print(compute.code_hash())
PY
)
PATHS=$(python -I tests/studio/sweep.py --print-paths --store "$STORE" \
        --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"})
read -r EXP_NAME EXP_NUM EXP_REL EXP_STORE CACHE <<<"$PATHS"

WANT=$(python -I tests/studio/sweep.py --list | wc -l)
MISS=$(python -I tests/studio/sweep.py --missing \
       --out "$CACHE" --seeds "$SEEDS" \
       --exp-num "$EXP_NUM" --exp-name "$EXP_NAME")
NMISS=$(printf '%s' "$MISS" | grep -c . || true)
HAVE=$(( WANT - NMISS ))
DIRS=$(find "$CACHE" -maxdepth 1 -type d -name "*__${HASH}" 2>/dev/null | wc -l || true)

echo "[pack] build hash   $HASH"
echo "[pack] experiment   $EXP_NAME #$EXP_NUM"
echo "[pack] cache        $CACHE"
echo "[pack] grid expects $WANT builds, $HAVE finished ($DIRS directories)"

if [ "$NMISS" -ne 0 ]; then
    IDS=$(printf '%s' "$MISS" | tr '\n' ' ')
    if [ "$PARTIAL" -eq 0 ]; then
        {
        echo "[pack] REFUSING: the grid is incomplete ($HAVE of $WANT)."
        python -I tests/studio/sweep.py --status \
            --out "$CACHE" --seeds "$SEEDS" \
            --exp-num "$EXP_NUM" --exp-name "$EXP_NAME" \
            2>&1 | sed 's/^/[pack] /' || true
        echo "[pack] $NMISS missing: $IDS"
        echo "[pack] why each task died, one file per task:"
        echo "[pack]     slurm_logs/fig2_precompute/failed/<IDX>.txt"
        echo "[pack] resubmit ONLY those builds, then this packs itself:"
        echo "[pack]     bash tests/studio/run_precompute.sh --resubmit"
        if [ "$NMISS" -le 64 ]; then
            echo "[pack] or by hand:"
            echo "[pack]     EXP_NUM=$EXP_NUM EXP_NAME=$EXP_NAME sbatch --array=$(printf '%s' "$MISS" | paste -sd, -) tests/studio/slurm_precompute.sh $CACHE"
        fi
        echo "[pack] or pass --partial to pack anyway."
        } >&2
        exit 3
    fi
    echo "[pack] --partial: packing anyway, $NMISS missing: $IDS"
fi

# EXP_REL is "." for the default experiment, which extracts into the store
# itself; any other experiment extracts into its own store under it.
EXTRACT="<repo>/Data/fig2_studio"
[ "$EXP_REL" != "." ] && EXTRACT="<repo>/Data/fig2_studio/$EXP_REL"

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
    # --exp-num matters: the daily kind reads the RAW run output, and its memo
    # key does not name the experiment, so warm_panels' default of 1 would
    # read experiment 1's runs into experiment 2's store, or fail where there
    # is no experiment 1
    if python -I tests/studio/warm_panels.py --store "$EXP_STORE" \
             --exp-num "$EXP_NUM" ${EXP_NAME:+--exp-name "$EXP_NAME"} \
             --procs "${SLURM_CPUS_PER_TASK:-4}" 2>&1 | sed -u 's/^/[pack] /'; then
        WARMED=yes
    else
        echo "[pack] WARNING: the warm-up failed; packing the grid anyway" >&2
        WARMED=partial
    fi
fi

STAMP=$(date +%Y%m%d)
OUT="$EXP_STORE/fig2_cache_${HASH}_${STAMP}.tar.gz"
MAN="${OUT%.tar.gz}.manifest.txt"

echo "[pack] writing $OUT"
# cache/_panels IS packed. It used to be excluded as "rebuilt in milliseconds",
# which is true of panels 1 to 3 and wrong about panel 4: that one reads every
# build in the grid, and a cold value measured 197 s on the 168-build run #2.
# It is ~13 MB against a 416 MB archive.
# tx/ is panel 2's transmission table (slurm_tx.sh); without it panel 2 has
# nothing to read. pigz compresses on every core where gzip uses one.
PARTS=(cache)
[ -f "$EXP_STORE/tx/index.json" ] && PARTS+=(tx) \
    || echo "[pack] WARNING: no $EXP_STORE/tx -- panel 2 will be empty" >&2
if command -v pigz >/dev/null 2>&1; then
    tar -cf - -C "$EXP_STORE" "${PARTS[@]}" | pigz -p "${SLURM_CPUS_PER_TASK:-4}" > "$OUT"
else
    tar -czf "$OUT" -C "$EXP_STORE" "${PARTS[@]}"
fi

SIZE=$(du -h "$OUT" | cut -f1)
SUM=$(sha256sum "$OUT" | cut -d' ' -f1)
NFILES=$(find "$CACHE" -maxdepth 2 -type f -path "*__${HASH}/*" 2>/dev/null | wc -l || true)

{
    echo "fig2 clade-build grid"
    echo "built        $(date -Is)"
    echo "build hash   $HASH"
    echo "experiment   $EXP_NAME #$EXP_NUM  (the run output these builds read)"
    echo "builds       $HAVE of $WANT finished (sweep.py done(), not dir count)"
    echo "directories  $DIRS at this hash"
    echo "panels       $WARMED  (cache/_panels is packed with the grid)"
    echo "contents     ${PARTS[*]}"
    echo "files        $NFILES"
    echo "archive      $(basename "$OUT")  $SIZE"
    echo "sha256       $SUM"
    echo "commit       $(git rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "dirty        $(git status --porcelain -- tests/studio/compute.py simplicity/clustering.py scripts/nature_plots/clade_preprocess.py | wc -l) of the 3 hashed files"
    echo
    echo "extract on the studio box with:"
    echo "    tar -xzf $(basename "$OUT") -C $EXTRACT"
} > "$MAN"

echo "[pack] $SIZE  sha256 ${SUM:0:16}..."
echo "[pack] manifest $MAN"
cat "$MAN"
