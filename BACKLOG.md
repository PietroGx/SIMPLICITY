# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

- [ ] **Temp tau hack in cal_2** — `impact_long_shedders_cal_2.py` currently sets
      `tau_3_long = D` with NO subtraction of (tau_1 + tau_2 + tau_4), to match the
      keys stored by the existing long calibration (`calibrate_long_nsr_#2`, which
      itself stored raw 63/109/365). This is a known inconsistency with the intended
      derivation `tau_3_long = D - (tau_1 + tau_2 + tau_4)`. Fix the derivation AND
      recalibrate before rerunning anything for production.
- [ ] cal_1 / cal_2 import consistency: if they still use
      `from experiments.experiment_script_runner import ...` and are launched by
      direct path, they hit the import error already fixed in the exp runner.

------------------------------------------------------------------------

## Active Refactor: Single Pipeline

Goal: collapse the current three-script flow (cal_1 → cal_2 → exp) into one
coherent, fully sequential pipeline. No manual multi-script invocation.

- [ ] Single shared parameter set for cal_1 and cal_2 (one source of truth;
      removes the USER_FIXED_PARAMS duplication between cal_2 and the exp runner).
- [ ] Remove the per-scenario Python loop in cal_2; submit all jobs at once using
      the native slurm reader (as done in the long-calibration step).
- [ ] Remove the per-scenario loop in the experiment runner (same de-looping).
- [ ] Unify into one fully sequential pipeline so the user does not run 3 scripts 
      (or bash script sequentially calling them)
- [ ] `hamming_iw` per-lineage cost in the sanity/intra-host plotting is heavy;
      run plotting on slurm for production-scale runs.


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



