#!/bin/bash
# ==============================================================================
# archive_paper_runs.sh -- one-shot bundle of everything behind the paper.
#
# Collects, for the two runs the figures are built from:
#
#   bound pipeline    impact_long_shedders_*_#<BOUND_NUM>            (default 7)
#   unbound pipeline  impact_long_shedders_unbound_*_#<UNBOUND_NUM>  (default 1)
#
# ...plus the calibration setup tables, pipeline logs, artifact zips and the
# generated figures, and hands the lot to scripts/package_simplicity_data.sh
# (pixz, multi-core, ~1 GB parts, integrity-verified) so it can be downloaded.
#
# Usage:
#   ./scripts/archive_paper_runs.sh                 # bound #7 + unbound #1
#   ./scripts/archive_paper_runs.sh 7 1             # explicit run numbers
#   ./scripts/archive_paper_runs.sh 7 1 --slurm     # submit the packing to SLURM
#
# Run from the repo root. Writes to Data_Export/ (gitignored).
#
# Nothing is deleted or moved: this only reads Data/ and writes archives.
# ==============================================================================
set -eo pipefail

BOUND_NUM=${1:-7}
UNBOUND_NUM=${2:-1}
SLURM_FLAG=${3:-""}

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

DATA_DIR="Data"
STAGING="${DATA_DIR}/paper_archive_bound${BOUND_NUM}_unbound${UNBOUND_NUM}"
EXPORT_DIR="Data_Export"

if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: no $DATA_DIR directory here. Run from the repo root on the HPC."
    exit 1
fi

echo "=============================================================="
echo " Paper archive: bound #${BOUND_NUM} + unbound #${UNBOUND_NUM}"
echo "=============================================================="

# --- STEP 1: collect, by HARD LINK so nothing is duplicated on disk --------
# Not symlinks: package_simplicity_data.sh runs `tar -I pixz -cf` without
# --dereference, so a symlinked staging tree would archive the links and none
# of the data. Hard links cost no extra space on the same filesystem and tar
# sees them as ordinary files. `cp -r` is the fallback if Data/ and the
# staging directory somehow straddle filesystems.
rm -rf "$STAGING"
mkdir -p "$STAGING"

collected=0
collect() {
    local pattern="$1"
    local label="$2"
    local found=0
    for path in $pattern; do
        [ -e "$path" ] || continue
        if ! cp -al "$path" "$STAGING/" 2>/dev/null; then
            echo "    (hard link failed for $(basename "$path") -- copying)"
            cp -r "$path" "$STAGING/"
        fi
        found=$((found + 1))
        collected=$((collected + 1))
    done
    printf '  %-34s %s\n' "$label" "$found"
}

echo ""
echo "--- STEP 1: COLLECTING ---"
collect "${DATA_DIR}/impact_long_shedders_*_#${BOUND_NUM}"            "bound run #${BOUND_NUM}"
collect "${DATA_DIR}/impact_long_shedders_unbound_*_#${UNBOUND_NUM}"  "unbound run #${UNBOUND_NUM}"
collect "${DATA_DIR}/pipeline_logs/*_#${BOUND_NUM}.log"               "bound pipeline log"
collect "${DATA_DIR}/pipeline_logs/*_#${UNBOUND_NUM}.log"             "unbound pipeline log"
collect "${DATA_DIR}/pipeline_artifacts/*"                            "artifact zips"
collect "${DATA_DIR}/figures"                                         "figures"
collect "${DATA_DIR}/RealWorldData"                                   "real-world inputs"

if [[ "$collected" -eq 0 ]]; then
    echo ""
    echo "Error: nothing matched. Check the run numbers and that you are on the"
    echo "machine holding Data/ (the HPC, not a laptop checkout)."
    rm -rf "$STAGING"
    exit 1
fi

echo ""
echo "  total items: $collected"
echo "  uncompressed: $(du -sh "$STAGING" 2>/dev/null | cut -f1)"

# --- STEP 2: a manifest, so the archive says what it is --------------------
MANIFEST="$STAGING/MANIFEST.txt"
{
    echo "SIMPLICITY paper archive"
    echo "created:  $(date -Iseconds)"
    echo "host:     $(hostname)"
    echo "commit:   $(git rev-parse --short HEAD 2>/dev/null || echo 'n/a')"
    echo "branch:   $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'n/a')"
    echo ""
    echo "bound pipeline   : impact_long_shedders_*_#${BOUND_NUM}"
    echo "unbound pipeline : impact_long_shedders_unbound_*_#${UNBOUND_NUM}"
    echo ""
    echo "contents:"
    for item in "$STAGING"/*; do
        [ "$(basename "$item")" = "MANIFEST.txt" ] && continue
        printf '  %s\n' "$(basename "$item")"
    done
} > "$MANIFEST"
echo "  manifest: $MANIFEST"

# --- STEP 3: hand to the existing pixz packager ----------------------------
echo ""
echo "--- STEP 2: COMPRESSING (pixz, via package_simplicity_data.sh) ---"
SUFFIX="paper_archive_bound${BOUND_NUM}_unbound${UNBOUND_NUM}"
./scripts/package_simplicity_data.sh "$DATA_DIR" "$SUFFIX" $SLURM_FLAG

if [[ "$SLURM_FLAG" == "--slurm" ]]; then
    echo ""
    echo "Packing submitted to SLURM. Archives will appear in ${EXPORT_DIR}/"
    echo "Staging kept at $STAGING until it finishes."
    exit 0
fi

# --- STEP 4: say plainly what to download ----------------------------------
echo ""
echo "=============================================================="
echo " WHAT TO DOWNLOAD"
echo "=============================================================="
shopt -s nullglob
archives=("${EXPORT_DIR}"/*"${SUFFIX}"*.xz "${EXPORT_DIR}"/simplicity_export_*.xz)
if [[ ${#archives[@]} -eq 0 ]]; then
    echo "  [warn] no archive found in ${EXPORT_DIR}/ -- check the output above."
else
    for a in "${archives[@]}"; do
        printf '  %s  (%s)\n' "$(cd "$(dirname "$a")" && pwd)/$(basename "$a")" \
                              "$(du -h "$a" | cut -f1)"
    done
    echo ""
    echo "  From your laptop:"
    echo "    scp '<user>@<hpc>:$(pwd)/${EXPORT_DIR}/*.xz' ."
fi
echo ""
echo "  Staging (hard links, safe to delete):"
echo "    $STAGING"
echo "=============================================================="
date
