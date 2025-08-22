#!/usr/bin/env bash
set -euo pipefail

# -------- light resources for the DRIVER jobs (not the seeds) --------
TIME_PER_DRIVER="1-00:00:00"   # keep driver alive while it manages its array
CPUS_DRIVER=2
MEM_DRIVER="2G"

# -------- paths --------
DRIVERS_GLOB="scripts/fit_OSR/[0-9][0-9]_generate_data_*.py"
BATCH_FIT="scripts/fit_OSR/batch_fit_long.py"
LOGDIR="scripts/fit_OSR/logs"
mkdir -p "$LOGDIR"

# -------- collect the 32 drivers --------
readarray -t FILES < <(ls -1 $DRIVERS_GLOB | sort)
if [ "${#FILES[@]}" -ne 32 ]; then
  echo "Expected 32 driver scripts, found ${#FILES[@]}." >&2
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
    --export=ALL,SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM=500,QT_QPA_PLATFORM=offscreen \
    ${prev_job:+--dependency=afterok:$prev_job} \
    --wrap "scripts/fit_OSR/run_driver_with_watchdog.sh \"$py\"")
  jid=$(awk '{print $4}' <<<"$out")
  echo "Submitted $base as job $jid (afterok: ${prev_job:-none})"
  prev_job="$jid"
done

# -------- submit the batch fit after the last driver completes --------
batch_out=$(sbatch \
  --time="00:30:00" \
  --cpus-per-task=2 \
  --mem="4G" \
  -o "$LOGDIR/batchfit_%j.out" \
  -e "$LOGDIR/batchfit_%j.err" \
  --dependency=afterok:${prev_job} \
  --export=ALL,QT_QPA_PLATFORM=offscreen \
  --wrap "python \"$BATCH_FIT\"")
batch_jid=$(awk '{print $4}' <<<"$batch_out")
echo "Submitted batch_fit_long.py as job $batch_jid (afterok: $prev_job)"

echo
echo "Monitor with:"
echo "  squeue -u $USER -o '%18i %20j %8T %10M %10l %5D %R'"
