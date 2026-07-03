# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

- [x] **Temp tau hack in cal_2** — fixed: `tau_3_long` is now derived as
      `inf_duration_long - (tau_1 + tau_2 + tau_4)` (shared
      `impact_long_shedders_config.derive_scenario_params`). cal_1's long
      calibration sweep now uses the corrected values too
      (`unique_long_taus`), so both sides agree. **Recalibration still
      required**: run the full pipeline under a fresh `--exp-num` — the old
      `calibrate_long_nsr_#2` was keyed by the raw (uncorrected) values and is
      no longer valid.
- [x] cal_1 / cal_2 import consistency — `impact_long_shedders_cal_2.py` was
      the one using the broken `from experiments.experiment_script_runner
      import ...`; fixed to the sibling-import form already used by cal_1 and
      the exp runner.

------------------------------------------------------------------------

## Active Refactor: Single Pipeline

Goal: collapse the current three-script flow (cal_1 → cal_2 → exp) into one
coherent, fully sequential pipeline. No manual multi-script invocation.

- [x] Single shared parameter set for cal_1 and cal_2 (one source of truth;
      removes the USER_FIXED_PARAMS duplication between cal_2 and the exp
      runner) — new `impact_long_shedders_config.py`.
- [x] Remove the per-scenario Python loop in cal_2; submit all jobs at once
      using the native slurm reader (as done in the long-calibration step) —
      via the new `_scenario_groups` mechanism in
      `settings_manager.generate_experiment_settings`.
- [x] Unify into one fully sequential pipeline so the user does not run 3
      scripts (or bash script sequentially calling them) — new
      `run_impact_long_shedders_pipeline.py`.
- [x] `hamming_iw` per-lineage cost in the sanity/intra-host plotting is heavy;
      run plotting on slurm for production-scale runs — new
      `submit_sanity_plot.sh`, dispatched once per long-shedder scenario by
      the pipeline script.


------------------------------------------------------------------------

## Backlog


------------------------------------------------------------------------

## Bug Fixes

- [ ] Slurm runner gets stuck on abrupt job terminations (e.g. OOM or OOT): the
      job is not reported as failed and the runner does not submit a new job.

------------------------------------------------------------------------

## Technical Debt

- [ ] Memory and runtime requests to slurm are hardcoded in the runner.
- [ ] Env variables are set in dir_manager (should be set elsewhere, or rename
      dir_manager).
- [ ] Tree builder (lines 43, 78): a tree node stores only the first of the IH
      lineages (for correct tree lineage coloring). If multiple lineages are
      stored, decide on a coloring strategy.
- [ ] Double substitutions at the same genome position are not explicitly handled.



