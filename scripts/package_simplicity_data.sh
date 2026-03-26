#!/bin/bash
# ==============================================================================
# SIMPLICITY PACKAGING (SLURM + MULTI-CORE + AUTO-VERIFY)
# ==============================================================================
SOURCE_DIR="$1"
SUFFIX="$2"
FLAG="$3"
EXPORT_DIR="Data_Export"
DATE_STAMP=$(date +%Y%m%d)

if [[ $# -lt 2 ]]; then
    echo "Usage: ./scripts/package_simplicity_data.sh [SOURCE_DIR] [SUFFIX] [--slurm]"
    exit 1
fi

if [[ "$FLAG" == "--slurm" ]]; then
    echo "Submitting packaging job to SLURM..."
    sbatch --job-name="simp_pkg" --cpus-per-task=16 --mem=32G --time=04:00:00 \
           --mail-type=END --wrap="./scripts/package_simplicity_data.sh $SOURCE_DIR '$SUFFIX'"
    exit 0
fi

if ! command -v pixz &> /dev/null; then echo "Error: pixz not found."; exit 1; fi

mkdir -p "$EXPORT_DIR"
mapfile -t TARGET_FOLDERS < <(find "$SOURCE_DIR" -maxdepth 1 -type d -name "*$SUFFIX" | sort)

echo "--- STEP 1: CLEANING LOGS ---"
for folder in "${TARGET_FOLDERS[@]}"; do [ -d "$folder/slurm" ] && rm -rf "$folder/slurm"; done

echo "--- STEP 2: MULTI-CORE COMPRESSION ---"
current_part=1
current_batch=()
current_batch_size=0
MAX_SIZE_BYTES=$((1 * 1024 * 1024 * 1024))

for folder in "${TARGET_FOLDERS[@]}"; do
    folder_size=$(du -sb "$folder" | cut -f1)
    if (( current_batch_size + folder_size > MAX_SIZE_BYTES )) && [ ${#current_batch[@]} -gt 0 ]; then
        archive_name="$EXPORT_DIR/simplicity_export_${DATE_STAMP}_part${current_part}.tar.xz"
        echo " -> Creating $archive_name..."
        tar -I pixz -cf "$archive_name" "${current_batch[@]}"
        ((current_part++))
        current_batch=()
        current_batch_size=0
    fi
    current_batch+=("$folder")
    current_batch_size=$((current_batch_size + folder_size))
done
[ ${#current_batch[@]} -gt 0 ] && tar -I pixz -cf "$EXPORT_DIR/simplicity_export_${DATE_STAMP}_part${current_part}.tar.xz" "${current_batch[@]}"

echo "--- STEP 3: INTEGRITY VERIFICATION ---"
ERRORS=0
for archive in "$EXPORT_DIR"/simplicity_export_${DATE_STAMP}_*.tar.xz; do
    echo -n " -> Verifying $(basename "$archive")... "
    if pixz -t "$archive" 2>/dev/null; then
        echo "OK"
    else
        echo "FAILED"
        ERRORS=$((ERRORS+1))
    fi
done

echo "=========================================================="
if [ "$ERRORS" -eq 0 ]; then
    echo "SUCCESS: Packaging and verification complete."
else
    echo "WARNING: $ERRORS archive(s) failed the integrity check!"
fi
echo "=========================================================="
