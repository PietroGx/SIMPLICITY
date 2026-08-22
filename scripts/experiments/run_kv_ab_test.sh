#!/usr/bin/env bash
# ============================================================================
# run_kv_ab_test.sh
# ----------------------------------------------------------------------------
# Launches two detached tmux sessions (same conda-activation pattern as
# scripts/start_simplicity_session.sh), each running a full
# run_impact_long_shedders_pipeline.py end to end -- one with
# --ih-virus-emergence-rate 0, one with 0.1 -- to test directly whether k_v
# is what's driving cal_2's inflated standard-NSR floor (see BACKLOG).
#
# Both runs use --cal-seeds 5 --exp-seeds 5 (cheap A/B pass, not a final
# calibration) and distinct --exp-num values so they write to separate
# Data/ trees and can run concurrently without clobbering each other. Each
# pipeline invocation is unchanged otherwise, so it still produces its own
# nsr_calibration_table.csv and, at the end, its own
# Data/pipeline_artifacts/impact_long_shedders_#{exp_num}_artifacts.zip
# (log + calibration-fit plots + sanity grid) via the pipeline's existing
# write_artifacts_archive step -- nothing new needed here for that.
#
# Sessions are left detached (not attached) since two run concurrently;
# reattach manually with `tmux attach -t <session>`.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONDA_ENV="simplicity"

EXP_NUM_KV0=2
EXP_NUM_KV1=3
CAL_SEEDS=5
EXP_SEEDS=5

SESSION_KV0="simplicity_kv0"
SESSION_KV1="simplicity_kv1"

launch() {
    local session="$1"
    local exp_num="$2"
    local kv="$3"

    if tmux has-session -t "$session" 2>/dev/null; then
        echo "[start] Session '$session' already exists -- leaving it alone."
        return
    fi

    echo "[start] Creating tmux session '$session' (exp_num=${exp_num}, k_v=${kv})..."
    tmux new-session -d -s "$session"

    tmux send-keys -t "$session" \
        'source "$(conda info --base)/etc/profile.d/conda.sh"' C-m
    tmux send-keys -t "$session" "conda activate ${CONDA_ENV}" C-m
    tmux send-keys -t "$session" 'export SBATCH_QOS=standard' C-m
    tmux send-keys -t "$session" "cd '${REPO_ROOT}'" C-m
    tmux send-keys -t "$session" \
        "python scripts/experiments/run_impact_long_shedders_pipeline.py \
--exp-num ${exp_num} \
--cal-seeds ${CAL_SEEDS} \
--exp-seeds ${EXP_SEEDS} \
--ih-virus-emergence-rate ${kv}" C-m
}

launch "$SESSION_KV0" "$EXP_NUM_KV0" 0
launch "$SESSION_KV1" "$EXP_NUM_KV1" 0.1

echo
echo "[done] Both pipelines launched detached."
echo "  k_v=0   -> tmux attach -t ${SESSION_KV0}   (exp_num=${EXP_NUM_KV0})"
echo "  k_v=0.1 -> tmux attach -t ${SESSION_KV1}   (exp_num=${EXP_NUM_KV1})"
