#!/bin/bash
# ==============================================================================
# SIMPLICITY RECONSTRUCTION (SLURM + HIGH-SPEED + AUTO-VERIFY)
# ==============================================================================
EXPORT_DIR="Data_Export"
TARGET_DIR="Data"
FLAG="$1"

if [[ "$FLAG" == "--slurm" ]]; then
    echo "Submitting reconstruction job to SLURM..."
    sbatch --job-name="simp_rec" --cpus-per-task=16 --mem=64G --time=02:00:00 \
           --mail-type=END --wrap="./scripts/reconstruct_simplicity_data.sh --internal"
    exit 0
fi

ARCHIVES=($(ls "$EXPORT_DIR"/simplicity_export_*.tar.xz 2>/dev/null))
if [ ${#ARCHIVES[@]} -eq 0 ]; then echo "No archives found."; exit 1; fi

mkdir -p "$TARGET_DIR"

echo "--- STEP 1: MULTI-CORE EXTRACTION ---"
for archive in "${ARCHIVES[@]}"; do
    echo " -> Extracting $(basename "$archive")"
    pixz -d -i "$archive" | tar -x -C .
done

echo "--- STEP 2: DATA VERIFICATION ---"
ERRORS=0
mapfile -t EXTRACTED < <(find "$TARGET_DIR" -maxdepth 1 -mindepth 1 -type d)
echo " -> Found ${#EXTRACTED[@]} directories in $TARGET_DIR/"

if [ ${#EXTRACTED[@]} -eq 0 ]; then
    echo " -> ERROR: No directories were extracted."
    ERRORS=$((ERRORS+1))
else
    for folder in "${EXTRACTED[@]}"; do
        SIZE=$(du -sk "$folder" | cut -f1)
        if [ "$SIZE" -lt 100 ]; then
            echo " -> WARNING: $folder is suspiciously empty (${SIZE}KB)."
            ERRORS=$((ERRORS+1))
        fi
    done
fi

echo "=========================================================="
if [ "$ERRORS" -eq 0 ]; then
    echo "SUCCESS: Reconstruction and verification complete."
else
    echo "WARNING: $ERRORS anomaly/anomalies detected in extracted data!"
fi
echo "=========================================================="
