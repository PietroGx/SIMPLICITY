#!/bin/bash
# ==============================================================================
# SIMPLICITY PACKAGING (SLURM + MULTI-CORE + AUTO-VERIFY)
# ==============================================================================
# Two ways to say what to pack:
#
#   package_simplicity_data.sh SOURCE_DIR SUFFIX [--slurm]
#       Every directory under SOURCE_DIR matching *SUFFIX. Used by
#       archive_paper_runs.sh and nsr_sensitivity_analysis.py, which build a
#       staging directory or own their whole name -- the match is safe there.
#
#   package_simplicity_data.sh --dir D [--dir D ...] [--tag TAG] [--slurm]
#       Exactly these directories. The pipeline uses this: a run's output
#       directories are knowable (impact_long_shedders_config.run_export_dirs),
#       and `Data/*_#2` also matched anything else ending in the same number --
#       a range-finder probe, a scratch experiment -- which then rode along in
#       the export.
#
# TAG goes into the archive name. Without it two runs packed on the same day
# both write simplicity_export_<date>_part1.tar.xz and the second overwrites
# the first.
EXPORT_DIR="Data_Export"
DATE_STAMP=$(date +%Y%m%d)

USAGE="Usage: ./scripts/package_simplicity_data.sh [SOURCE_DIR] [SUFFIX] [--slurm]
       ./scripts/package_simplicity_data.sh --dir D [--dir D ...] [--tag TAG] [--slurm]"

ORIGINAL_ARGS=("$@")
SOURCE_DIR=""
SUFFIX=""
TAG=""
USE_SLURM=0
EXPLICIT_DIRS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dir)   [[ -n "$2" ]] || { echo "--dir needs a directory"; exit 1; }
                 EXPLICIT_DIRS+=("$2"); shift 2 ;;
        --tag)   [[ -n "$2" ]] || { echo "--tag needs a value"; exit 1; }
                 TAG="$2"; shift 2 ;;
        --slurm) USE_SLURM=1; shift ;;
        -*)      echo "Unknown option: $1"; echo "$USAGE"; exit 1 ;;
        *)       if [[ -z "$SOURCE_DIR" ]]; then SOURCE_DIR="$1"
                 elif [[ -z "$SUFFIX" ]]; then SUFFIX="$1"
                 else echo "Unexpected argument: $1"; echo "$USAGE"; exit 1; fi
                 shift ;;
    esac
done

if [[ ${#EXPLICIT_DIRS[@]} -gt 0 ]]; then
    if [[ -n "$SOURCE_DIR" ]]; then
        echo "--dir and SOURCE_DIR/SUFFIX are different modes; pick one."
        echo "$USAGE"; exit 1
    fi
elif [[ -z "$SOURCE_DIR" || -z "$SUFFIX" ]]; then
    echo "$USAGE"; exit 1
fi

if [[ "$USE_SLURM" == "1" ]]; then
    echo "Submitting packaging job to SLURM..."
    # Rebuilt from the real arguments rather than re-spelled, so a new flag
    # cannot be dropped on the way to the compute node.
    WRAPPED="./scripts/package_simplicity_data.sh"
    for arg in "${ORIGINAL_ARGS[@]}"; do
        [[ "$arg" == "--slurm" ]] && continue
        WRAPPED+=" $(printf '%q' "$arg")"
    done
    sbatch --job-name="simp_pkg" --cpus-per-task=16 --mem=32G --time=04:00:00 \
           --mail-type=END --export=ALL,SIMPLICITY_KEEP_SLURM_LOGS="${SIMPLICITY_KEEP_SLURM_LOGS:-0}" \
           --wrap="$WRAPPED"
    exit 0
fi

if ! command -v pixz &> /dev/null; then echo "Error: pixz not found."; exit 1; fi

mkdir -p "$EXPORT_DIR"

if [[ ${#EXPLICIT_DIRS[@]} -gt 0 ]]; then
    TARGET_FOLDERS=()
    for folder in "${EXPLICIT_DIRS[@]}"; do
        if [[ -d "$folder" ]]; then
            TARGET_FOLDERS+=("$folder")
        else
            echo " -> [skip] not on disk: $folder"
        fi
    done
    if [[ ${#TARGET_FOLDERS[@]} -eq 0 ]]; then
        echo "Error: none of the requested directories exist; nothing to pack."
        exit 1
    fi
    echo "--- PACKING ${#TARGET_FOLDERS[@]} named director(ies) ---"
    printf '  %s\n' "${TARGET_FOLDERS[@]}"
else
    mapfile -t TARGET_FOLDERS < <(find "$SOURCE_DIR" -maxdepth 1 -type d -name "*$SUFFIX" | sort)
fi

if [[ -n "$TAG" ]]; then
    ARCHIVE_PREFIX="$EXPORT_DIR/simplicity_export_${TAG}_${DATE_STAMP}"
else
    ARCHIVE_PREFIX="$EXPORT_DIR/simplicity_export_${DATE_STAMP}"
fi

# Set SIMPLICITY_KEEP_SLURM_LOGS=1 to archive the slurm/ directories instead of
# deleting them. They are the only record of why a task failed -- a truncated
# pickle, a numpy ValueError -- and this runs in place on the real Data/, not
# on a copy, so the default is destructive and irreversible.
if [[ "${SIMPLICITY_KEEP_SLURM_LOGS:-0}" == "1" ]]; then
    echo "--- STEP 1: KEEPING LOGS (SIMPLICITY_KEEP_SLURM_LOGS=1) ---"
else
    echo "--- STEP 1: CLEANING LOGS ---"
    for folder in "${TARGET_FOLDERS[@]}"; do [ -d "$folder/slurm" ] && rm -rf "$folder/slurm"; done
fi

echo "--- STEP 2: MULTI-CORE COMPRESSION ---"
current_part=1
current_batch=()
current_batch_size=0
MAX_SIZE_BYTES=$((1 * 1024 * 1024 * 1024))

for folder in "${TARGET_FOLDERS[@]}"; do
    folder_size=$(du -sb "$folder" | cut -f1)
    if (( current_batch_size + folder_size > MAX_SIZE_BYTES )) && [ ${#current_batch[@]} -gt 0 ]; then
        archive_name="${ARCHIVE_PREFIX}_part${current_part}.tar.xz"
        echo " -> Creating $archive_name..."
        tar -I pixz -cf "$archive_name" "${current_batch[@]}"
        ((current_part++))
        current_batch=()
        current_batch_size=0
    fi
    current_batch+=("$folder")
    current_batch_size=$((current_batch_size + folder_size))
done
[ ${#current_batch[@]} -gt 0 ] && tar -I pixz -cf "${ARCHIVE_PREFIX}_part${current_part}.tar.xz" "${current_batch[@]}"

echo "--- STEP 3: INTEGRITY VERIFICATION ---"
# `pixz -l` lists the tarball index: a read-only check that the archive is
# both valid xz AND a readable tar.
#
# This used to be `pixz -t`, which is NOT a test flag -- in pixz, -t means
# "don't assume input is in tar format", i.e. it COMPRESSES. Every archive was
# therefore silently re-compressed to <name>.tar.xz.xz and the original
# deleted, while the script reported SUCCESS. The result needed `pixz -d`
# twice to open and could not be read by plain `tar -xf` or `xz -d` at all.
ERRORS=0
shopt -s nullglob
ARCHIVES=("${ARCHIVE_PREFIX}"_*.tar.xz)
if [ ${#ARCHIVES[@]} -eq 0 ]; then
    echo " -> No archives produced; nothing to verify."
    ERRORS=1
fi
for archive in "${ARCHIVES[@]}"; do
    echo -n " -> Verifying $(basename "$archive")... "
    if pixz -l "$archive" >/dev/null 2>&1; then
        echo "OK"
    else
        echo "FAILED"
        ERRORS=$((ERRORS+1))
    fi
done

echo "=========================================================="
if [ "$ERRORS" -eq 0 ]; then
    echo "SUCCESS: Packaging and verification complete."
    echo "Extract with:  tar -I pixz -xf <archive>"
else
    echo "WARNING: $ERRORS archive(s) failed the integrity check!"
fi
echo "=========================================================="
