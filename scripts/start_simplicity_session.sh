#!/usr/bin/env bash
# ============================================================================
# start_simplicity_session.sh
# Launch (or reattach to) a tmux session with the `simplicity` conda env
# active and SBATCH_QOS=standard exported inside the session.
# ============================================================================

set -euo pipefail

SESSION="simplicity"
CONDA_ENV="simplicity"

# If the session already exists, just attach to it (avoid stacking activations).
if tmux has-session -t "$SESSION" 2>/dev/null; then
    echo "[start] Session '$SESSION' already exists -- attaching."
    exec tmux attach -t "$SESSION"
fi

# Create a new detached session.
echo "[start] Creating tmux session '$SESSION'..."
tmux new-session -d -s "$SESSION"

# Inside the session: source conda's shell hook, activate the env, export QOS.
tmux send-keys -t "$SESSION" \
    'source "$(conda info --base)/etc/profile.d/conda.sh"' C-m
tmux send-keys -t "$SESSION" "conda activate ${CONDA_ENV}" C-m
tmux send-keys -t "$SESSION" 'export SBATCH_QOS=standard' C-m
tmux send-keys -t "$SESSION" 'echo "[ready] env=${CONDA_DEFAULT_ENV}  SBATCH_QOS=${SBATCH_QOS}"' C-m

# Attach.
exec tmux attach -t "$SESSION"
