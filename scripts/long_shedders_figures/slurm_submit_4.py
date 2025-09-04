#!/usr/bin/env python3
# ASCII-only, UTF-8 safe

import os
import stat
import argparse
import subprocess
import shlex

BASH_TEMPLATE = """#!/usr/bin/env bash
set -euo pipefail

# Optional: activate env or modules (uncomment and edit one of these lines)
# source ~/.bashrc && conda activate simplicity
# module load python/3.10

# Arrays of parameter values (bash arrays)
CLUSTERS=({clusters})
MINDAYS=({mindays})

# Derived from array task id
IDX="${{SLURM_ARRAY_TASK_ID:-0}}"
N_MD="${{#MINDAYS[@]}}"
CLUSTER="${{CLUSTERS[$((IDX / N_MD))]}}"
MD="${{MINDAYS[$((IDX % N_MD))]}}"

echo "[INFO] JobID=${{SLURM_JOB_ID:-NA}} Task=${{SLURM_ARRAY_TASK_ID:-0}} -> CLUSTER=$CLUSTER MD=$MD"

# Ensure output directory exists
mkdir -p "{outdir}"
mkdir -p "{logdir}"

# Run the figure script
python {script_path} --cluster-threshold "$CLUSTER" --min-days "$MD" --index-csv {index_csv} --outdir "{outdir}"
"""

def main():
    ap = argparse.ArgumentParser(description="Submit a single Slurm job array to run 4try.py across parameter grids.")
    ap.add_argument("--script-path", default="scripts/long_shedders_figures/4try.py",
                    help="Path to the figure-building Python script (default: scripts/long_shedders_figures/figure_4_long_grid.py)")
    ap.add_argument("--index-csv", default="scripts/long_shedders_experiments/INDEX_exp_scripts.csv",
                    help="Path to INDEX_exp_scripts.csv")
    ap.add_argument("--outdir", default="Data/Fig4_Results",
                    help="Output directory for the two figures written by 4try.py")
    ap.add_argument("--clusters", default="5,7,9,11,13",
                    help="Comma-separated cluster thresholds (default: 5,7,9,11,13)")
    ap.add_argument("--mindays", default="600,900",
                    help="Comma-separated min-days values (default: 600,900)")
    ap.add_argument("--job-name", default="fig4_grid",
                    help="Slurm job name (default: fig4_grid)")
    ap.add_argument("--logdir", default="logs",
                    help="Directory for Slurm logs (default: logs)")
    ap.add_argument("--time", default="02:00:00",
                    help="Slurm time limit (default: 02:00:00)")
    ap.add_argument("--cpus", default="4",
                    help="Slurm --cpus-per-task (default: 4)")
    ap.add_argument("--mem", default="16G",
                    help="Slurm --mem (default: 16G)")
    args = ap.parse_args()

    # Prepare arrays
    clusters_list = [x.strip() for x in args.clusters.split(",") if x.strip()]
    mindays_list = [x.strip() for x in args.mindays.split(",") if x.strip()]
    if not clusters_list or not mindays_list:
        raise SystemExit("Both --clusters and --mindays must contain at least one value.")

    n_jobs = len(clusters_list) * len(mindays_list)

    # Make sure logs dir exists
    os.makedirs(args.logdir, exist_ok=True)

    # Build bash script content
    clusters_str = " ".join(clusters_list)
    mindays_str = " ".join(mindays_list)
    bash_body = BASH_TEMPLATE.format(
        clusters=clusters_str,
        mindays=mindays_str,
        outdir=args.outdir,
        logdir=args.logdir,
        script_path=shlex.quote(args.script_path),
        index_csv=shlex.quote(args.index_csv),
    )

    # Write the bash runner
    runner_path = os.path.join("scripts", "long_shedders_figures", "run_fig4_grid.sh")
    os.makedirs(os.path.dirname(runner_path), exist_ok=True)
    with open(runner_path, "w", encoding="utf-8") as f:
        f.write(bash_body)

    # Make it executable
    st = os.stat(runner_path)
    os.chmod(runner_path, st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    # Submit array job
    sbatch_cmd = [
        "sbatch",
        f"--job-name={args.job_name}",
        f"--output={os.path.join(args.logdir, args.job_name)}_%A_%a.out",
        f"--error={os.path.join(args.logdir, args.job_name)}_%A_%a.err",
        f"--time={args.time}",
        f"--cpus-per-task={args.cpus}",
        f"--mem={args.mem}",
        f"--array=0-{n_jobs-1}",
        runner_path,
    ]

    print("[SUBMIT]", " ".join(sbatch_cmd))
    out = subprocess.check_output(sbatch_cmd, text=True).strip()
    print(out)
    print(f"[INFO] Array size: {n_jobs} tasks "
          f"({len(clusters_list)} clusters x {len(mindays_list)} min-days)")

if __name__ == "__main__":
    main()

