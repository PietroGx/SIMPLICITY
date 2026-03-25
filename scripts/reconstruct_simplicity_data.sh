#!/bin/bash
# SIMPLICITY DATA RECONSTRUCTION (SLURM + PIXZ + PROGRESS)
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
for archive in "${ARCHIVES[@]}"; do
    echo "Processing $(basename "$archive")..."
    if [[ "$FLAG" == "--internal" ]]; then
        pixz -d < "$archive" | tar -x -C .
    else
        if command -v pv &> /dev/null; then
            pv "$archive" | pixz -d | tar -x -C .
        else
            pixz -d < "$archive" | tar -x -C .
        fi
    fi
done
echo "Reconstruction complete in $TARGET_DIR/"
