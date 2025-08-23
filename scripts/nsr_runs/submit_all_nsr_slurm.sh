#!/usr/bin/env bash
set -euo pipefail

# -------- light resources for the DRIVER jobs (not the seeds) --------
TIME_PER_DRIVER="1-00:00:00"   # long enough for a driver to manage its array
CPUS_DRIVER=2
MEM_DRIVER="2G"

# -------- paths --------
DRIVERS_GLOB="scripts/nsr_runs/[0-9][0-9]_generate_data_*_NSR_run.py"
LOGDIR="scripts/nsr_runs/logs"
WATCHDOG="scripts/fit_OSR/run_driver_with_watchdog.sh"   # reuse the same watchdog
mkdir -p "$LOGDIR"

# -------- collect drivers --------
readarray -t FILES < <(ls -1 $DRIVERS_GLOB | sort)
if [ "${#FILES[@]}" -eq 0 ]; then
  echo "No NSR-run scripts found under $DRIVERS_GLOB" >&2
  exit 1
fi

# -------- submit drivers chained with dependencies --------
prev_job=""
for py in "${FILES[@]}"; do
  base=$(basename "$py")
  out=$(sbatch \
    --time="$TIME_PER_DRIVER" \
    --cpus-per-task="$CPUS_DRIVER" \
    --mem="$MEM_DRIVER" \
    -o "$LOGDIR/${base%.py}_%j.out" \
    -e "$LOGDIR/${base%.py}_%j.err" \
    --export=ALL,SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM=500,QT_QPA_PLATFORM=offscreen,PYTHONPATH=$PWD \
    ${prev_job:+--dependency=afterok:$prev_job} \
    --wrap "$WATCHDOG \"$py\"")
  jid=$(awk '{print $4}' <<<"$out")
  echo "Submitted $base as job $jid (afterok: ${prev_job:-none})"
  prev_job="$jid"
done

echo
echo "Monitor with:"
echo "  squeue -u \$USER -o '%18i %20j %8T %10M %10l %5D %R' | grep -E 'NSR_run|osr|batchfit'"
