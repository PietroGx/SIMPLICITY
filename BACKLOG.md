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
- [x] **cal_1 R_long inconsistency** — fixed: cal_1's isolated calibration
      was feeding the same raw `R_long=1.5` (a weekly per-long-shedder
      infection rate) directly into every `tau_3_long` grid point as the
      final R_long, unscaled by duration, while production derives
      `R_long = R_long_per_week * tau_3_long / 7.0`
      (`derive_scenario_params`). Most tau groups were calibrated at the
      wrong R_long. Fixed via a shared
      `derive_r_long(R_long_per_week, tau_3_long)` helper and a per-tau
      `_scenario_groups` submission in cal_1.py (same mechanism cal_2
      already uses) -- each context scales its own weekly baseline
      (`CAL1_ISOLATED_FIXED_PARAMS["R_long"]`=1.5 for cal_1, `sp["R_long"]`
      for cal_2/exp) by that point's own `tau_3_long/7`.
      **Recalibration required**: any experiment run before this fix used
      the wrong R_long for cal_1's grid and must be redone under a fresh
      `--exp-num`.

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
- [x] cal_1's isolated-calibration `fixed_params` moved out of the script body
      into `impact_long_shedders_config.CAL1_ISOLATED_FIXED_PARAMS` — every
      hardcoded parameter set the pipeline uses now lives in that one config
      file (still intentionally distinct values from `USER_FIXED_PARAMS`,
      just no longer duplicated inline in a script).


------------------------------------------------------------------------

## Figures 3 & 4 (long-shedders paper)

- [x] Figure 3: pies (one selected scenario) + metrics PCA (Peak/Burden/
      Survival/Growth, all scenarios, one point per seed) + sequence-space
      PCA (selected scenario) + per-seed excess-Hamming-divergence
      consistency check — new
      `scripts/nature_plots/fig3_{preprocess_data,plots,save}.py`.
- [x] Extended `get_clade_winners` into `get_clade_metrics`
      (`long_shedders_preprocess.py`): same clustering pass now also returns
      continuous long-origin values per metric, not just the categorical
      winner. Added `get_lineage_labels` (per-lineage, not per-clade).
- [x] Figure 4: 5-scenario phylogenetic tree comparison (same seed), colored
      by the same long/standard/mixed/founder labels as Figure 3 — new
      `scripts/nature_plots/fig4_{preprocess_data,plots,save}.py`, via
      `baltic.loadNewick` + an opt-in `label_internal_nodes` fix to
      `simplicity/tree/newick.py` (gated to phylogenetic trees only —
      infection-tree Newick export is unaffected).
- [x] `submit_plots.sh`: added `"3"` and `"4"` branches (and into `"all"`).

------------------------------------------------------------------------

## Backlog

- [ ] Figure 3 Panel B: try the long-vs-standard contrast version of the 4
      metrics (instead of/alongside the current raw long-origin values) —
      deferred during design, worth comparing once real data exists.
- [ ] Figure 3 Panel C: test whether `max_long_per_individual` subsampling is
      needed for the sequence-space PCA (full long-shedder data vs.
      subsampled) once real experiment output exists.
- [ ] Figure 3 Panel C: verify the per-seed pairwise-Hamming cap (200 sampled
      pairs) in the consistency check is generous enough once real genome
      counts per seed are known; tune if not.

------------------------------------------------------------------------

## Bug Fixes

- [ ] Slurm runner gets stuck on abrupt job terminations (e.g. OOM or OOT): the
      job is not reported as failed and the runner does not submit a new job.

------------------------------------------------------------------------

## Monitor visibility

- [x] SLURM monitor had no visibility into a running simulation's internal
      state (only completed/running/failed counts) — fixed: `extrande.py`'s
      `ProgressReporter` optionally writes a throttled JSON snapshot
      (`time`/`final_time`/`infected`) to
      `<seeded_simulation_parameters_path>.progress` (same signal-file
      convention as `.started`/`.completed`/`.failed`), threaded through
      `simulation.py` and `unit_run.py`. `simplicity/runners/slurm.py`'s
      polling loop gained a separate hourly timer that reports the internal
      state of any simulation running longer than an hour, reading that
      snapshot — independent of the existing ~17s `SimulationsStatus` line.

------------------------------------------------------------------------

## Technical Debt

- [x] Memory and runtime requests to slurm are hardcoded in the runner —
      fixed: `simplicity/runners/slurm.py` now reads `SIMPLICITY_SLURM_MEM`/
      `SIMPLICITY_SLURM_TIME` env vars (same pattern already used for
      `SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM`, so the shared
      `run_seeded_simulations(...)` interface across serial/multiprocessing/
      slurm didn't need to change). `impact_long_shedders_config.py` exposes
      `add_slurm_resource_args`/`set_slurm_resource_env` so cal_1/cal_2/exp/
      the pipeline orchestrator all expose identical `--slurm-mem`/
      `--slurm-time` flags. Memory default lowered to `2G` based on real
      `sacct` data on the isolated long-NSR calibration grid (MaxRSS
      0.29–0.45G across 36 tasks) — down from the old blanket `8G` every job
      used to request regardless of actual need. Time limit kept at the
      existing `1-00:00:00` default (observed elapsed time up to ~88 min
      already has plenty of headroom there, and unlike memory there's no
      cost to leaving it generous).
- [ ] Env variables are set in dir_manager (should be set elsewhere, or rename
      dir_manager).
- [ ] Tree builder (lines 43, 78): a tree node stores only the first of the IH
      lineages (for correct tree lineage coloring). If multiple lineages are
      stored, decide on a coloring strategy.
- [ ] Double substitutions at the same genome position are not explicitly handled.



