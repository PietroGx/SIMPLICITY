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
# `run_precompute.sh` submits this with --dependency=afterok on the array, so
# it runs only if every build succeeded. It can also be run by hand.
#
# It REFUSES to pack a partial grid. The whole point of one file is that what
# arrives is the whole grid; a tarball holding 300 of 384 builds looks the same
# from the outside and would be found out later, on the wrong machine.
#
# Builds are counted by the current code hash, not by directory name: a
# directory left over from an earlier checkout is not part of this grid.
#
# Usage, from the repo root:
#     sbatch tests/studio/pack_cache.sh            # the whole grid
#     bash   tests/studio/pack_cache.sh --partial  # pack what exists anyway
#
# No 'set -u': conda's activate.d hooks are not nounset-safe.
# ============================================================================
set -eo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

STORE="${STORE:-$HOME/fig2_studio}"
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
WANT=$(python -I tests/studio/sweep.py --list | wc -l)
HAVE=$(find "$STORE/cache" -maxdepth 1 -type d -name "*__${HASH}" 2>/dev/null | wc -l)

echo "[pack] build hash   $HASH"
echo "[pack] grid expects $WANT builds, found $HAVE at this hash"

if [ "$HAVE" -ne "$WANT" ] && [ "$PARTIAL" -eq 0 ]; then
    echo "[pack] REFUSING: the grid is incomplete ($HAVE of $WANT)." >&2
    echo "[pack] Resubmit the array -- a finished build is skipped, so it only" >&2
    echo "[pack] rebuilds what is missing. Or pass --partial to pack anyway." >&2
    exit 1
fi

STAMP=$(date +%Y%m%d)
OUT="$STORE/fig2_cache_${HASH}_${STAMP}.tar.gz"
MAN="${OUT%.tar.gz}.manifest.txt"

echo "[pack] writing $OUT"
# _panels is the local panel cache: derived, keyed on the same hashes, and
# rebuilt in milliseconds. It is not part of the grid.
tar -czf "$OUT" -C "$STORE" --exclude="cache/_panels" cache

SIZE=$(du -h "$OUT" | cut -f1)
SUM=$(sha256sum "$OUT" | cut -d' ' -f1)
NFILES=$(find "$STORE/cache" -maxdepth 2 -type f -path "*__${HASH}/*" | wc -l)

{
    echo "fig2 clade-build grid"
    echo "built        $(date -Is)"
    echo "build hash   $HASH"
    echo "builds       $HAVE of $WANT expected"
    echo "files        $NFILES"
    echo "archive      $(basename "$OUT")  $SIZE"
    echo "sha256       $SUM"
    echo "commit       $(git rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "dirty        $(git status --porcelain -- tests/studio/compute.py simplicity/clustering.py scripts/nature_plots/clade_preprocess.py | wc -l) of the 3 hashed files"
    echo
    echo "extract on the studio box with:"
    echo "    tar -xzf $(basename "$OUT") -C ~/fig2_studio"
} > "$MAN"

echo "[pack] $SIZE  sha256 ${SUM:0:16}..."
echo "[pack] manifest $MAN"
cat "$MAN"
