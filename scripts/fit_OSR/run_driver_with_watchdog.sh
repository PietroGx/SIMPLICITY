#!/usr/bin/env bash
set -euo pipefail

PY="$1"                          # path to the driver .py (e.g., scripts/fit_OSR/01_generate_data_...py)

# Derive experiment name from filename: NN_generate_data_<EXP>.py
BN="$(basename "$PY")"
EXP="${BN#*_generate_data_}"
EXP="${EXP%.py}"

echo "[watchdog] starting driver for $EXP"
python "$PY" slurm 1 &
PID=$!

# Poll every 30s; cancel ONLY tasks held by admin and mark them .failed so the runner advances.
while kill -0 "$PID" 2>/dev/null; do
  # Get the array JobID for this experiment (if submitted)
  JID="$(squeue --name "$EXP" -h -o %A | head -n1 || true)"
  if [[ -n "${JID:-}" ]]; then
    # Find PENDING tasks whose reason is JobHeldAdmin, return as <jobid>_<taskid>
    mapfile -t HELD < <(squeue -h -j "$JID" -t PENDING -o "%A_%a %R" | awk '$2=="JobHeldAdmin"{print $1}')
    if (( ${#HELD[@]} )); then
      echo "[watchdog] $EXP: found ${#HELD[@]} admin-held tasks; cancelling those and marking .failed"
      for jid_tid in "${HELD[@]}"; do
        scancel "$jid_tid" || true
        tid="${jid_tid#*_}"            # SLURM_ARRAY_TASK_ID (1-based)

        # Touch the .failed signal for that seed so the runner's accounting is consistent
        EXP_NAME="$EXP" TID="$tid" PYTHONPATH="$PWD" python - <<'PY'
import os, pathlib
import simplicity.settings_manager as sm
exp = os.environ["EXP_NAME"]
tid = int(os.environ["TID"])  # 1-based index into seeded parameters
paths = sm.get_seeded_simulation_parameters_paths(exp)
p = paths[tid-1]
pathlib.Path(p + ".failed").touch()
PY
      done
    fi
  fi
  sleep 30
done

wait "$PID"
exit $?
