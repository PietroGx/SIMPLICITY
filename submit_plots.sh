#!/bin/bash
#SBATCH --job-name=simplicity_plots
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=main
#SBATCH --qos=standard

# Initialize Conda
eval "$(/home/gerlep93/miniconda3/bin/conda shell.bash hook)"
conda activate simplicity

# Capture optional arguments for Figure 2 customization
FIG_CHOICE=${1:-"all"}
SEED=${2:-"1"}
FORMAT=${3:-"png"}

# Create dedicated log directory
mkdir -p slurm_logs/fig${FIG_CHOICE}

echo "Starting SIMPLICITY plotting job for Figure $FIG_CHOICE (Seed: $SEED, Format: $FORMAT)..."
date

# Helper to run python script
run_fig2() {
    python scripts/nature_plots/fig2_save.py \
        --seed "$SEED" \
        --format "$FORMAT" \
        > slurm_logs/fig2/plots_${SLURM_JOB_ID}.log 2> slurm_logs/fig2/plots_${SLURM_JOB_ID}.err
}

if [ "$FIG_CHOICE" == "1" ]; then
    python scripts/nature_plots/fig1_save.py > slurm_logs/fig1/plots_${SLURM_JOB_ID}.log 2> slurm_logs/fig1/plots_${SLURM_JOB_ID}.err
elif [ "$FIG_CHOICE" == "2" ]; then
    run_fig2
elif [ "$FIG_CHOICE" == "all" ]; then
    python scripts/nature_plots/fig1_save.py > slurm_logs/fig1/plots_${SLURM_JOB_ID}.log 2> slurm_logs/fig1/plots_${SLURM_JOB_ID}.err
    run_fig2
else
    echo "Error: Invalid argument '$FIG_CHOICE'."
    exit 1
fi

echo "Job completed."
date