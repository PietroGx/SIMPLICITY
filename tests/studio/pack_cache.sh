#!/bin/bash
#SBATCH --job-name=fig2_pack
#SBATCH --time=00:40:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
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
# arrives is the whole grid; a tarball holding 300 of 384 builds looks the same
# from the outside and would be found out later, on the wrong machine.
#
# Completeness comes from `sweep.py --missing`, which is the same `done()` test
# that makes a resubmitted array skip finished work. It is NOT a directory
# count: a build that died halfway leaves a directory with a missing or partial
# index.json, and a directory count calls that present.
#
# Exit codes, so a refusal is not just "FAILED" in sacct:
#     0  packed
#     3  refused, the grid is incomplete
#     1  something else broke
#
# It packs ONE experiment's store, named by EXP_NUM/EXP_NAME and resolved by
# `sweep.py --print-paths` -- the same resolution run_precompute.sh and
# slurm_precompute.sh use. The manifest records which experiment the runs came
# from: nothing inside a build does, because compute.params_key and index.json
# hold the clade parameters and the code hash and no experiment identity.
#
# Usage, from the repo root:
#     sbatch tests/studio/pack_cache.sh            # the whole grid
#     bash   tests/studio/pack_cache.sh --partial  # pack what exists anyway
#     EXP_NUM=3 bash tests/studio/pack_cache.sh    # another experiment's store
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

STORE="${STORE:-$HOME/fig2_studio}"
SEEDS="${SEEDS:-50}"
EXP_NUM="${EXP_NUM:-1}"
EXP_NAME="${EXP_NAME:-}"
PARTIAL=0
[ "${1:-}" = "--partial" ] && PARTIAL=1

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
EXTRACT="~/fig2_studio"
[ "$EXP_REL" != "." ] && EXTRACT="~/fig2_studio/$EXP_REL"

STAMP=$(date +%Y%m%d)
OUT="$EXP_STORE/fig2_cache_${HASH}_${STAMP}.tar.gz"
MAN="${OUT%.tar.gz}.manifest.txt"

echo "[pack] writing $OUT"
# _panels is the local panel cache: derived, keyed on the same hashes, and
# rebuilt in milliseconds. It is not part of the grid.
tar -czf "$OUT" -C "$EXP_STORE" --exclude="cache/_panels" cache

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
