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
      `derive_r_long(r_long_weekly_rate, tau_3_long)` helper and a per-tau
      `_scenario_groups` submission in cal_1.py (same mechanism cal_2
      already uses) -- each context scales its own weekly baseline
      (`CAL1_R_LONG_WEEKLY_RATE`=1.5 for cal_1, `sp["R_long"]` for cal_2/exp)
      by that point's own `tau_3_long/7`.
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
- [x] **Config/pipeline separation** — settings-construction
      (`user_set_experiment_settings`, `build_scenario_groups`,
      `build_scenario_settings`) was still defined inline in cal_1.py/
      cal_2.py/exp.py, not in the config file, so "what parameters get
      built" and "how the pipeline runs" were mixed together in all three
      scripts. Moved all of it into `impact_long_shedders_config.py`
      (`build_cal1_settings`, `build_cal2_scenario_groups`/
      `build_cal2_settings`, `build_exp_scenario_settings`, `lookup_long_nsr`)
      — the three stage scripts now only parse args, call these builders,
      and dispatch/handle I/O. No stage script constructs a parameter dict
      itself anymore.


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
- [ ] **cal_1's isolated calibration: edge_case sequencing shortfall vs.
      SOT depletion, sharing one parameter** — edge_case
      (`tau_3_long=350.23`) reaches the full 365-day `final_time` but only
      accumulates 3-8 sequenced diagnoses per seed (needs `min_seq=30`) --
      `tau_3_long` is almost as long as `final_time`, so individuals barely
      complete a full cycle. Raising `R_long` doesn't help (confirmed).
      Raising `CAL1_ISOLATED_FIXED_PARAMS['infected_individuals_at_start']`
      (10->60) does (3-8 -> 20-39 seqs), but at 60 it also pushes SOT
      (`R_long=1.3`) into early depletion (most seeds end via "No
      susceptibles left" well before 365 days) -- a single shared parameter
      with competing pressures across scenarios. Needs either an
      intermediate value that satisfies both (untested), or a per-scenario
      starting-cohort size (would need a new SCENARIOS field, out of
      tonight's tuning scope). `tests/test_long_shedder_r_long_grid.py`
      (new) can help explore this systematically on the HPC. Also open:
      SOT/HIV_low's calibration fits are currently poor (R² near 0 / weak)
      at only 5 diagnostic seeds/NSR-point -- not yet determined whether
      that's a genuine `NSR_RANGES['cal1_long_nsr']` mismatch or just
      statistical noise from the small diagnostic seed count.

------------------------------------------------------------------------

## Bug Fixes

- [x] **Slurm runner gets stuck on abrupt job terminations (e.g. OOM or OOT)**
      — fixed: root-caused via a real HPC run (`impact_long_shedders_calibration_lng_nsr_#1`,
      `edge_case` tau group) using the new `[long-running]` monitor report
      itself, which showed a task's progress snapshot frozen unchanged for
      30+ hours. A Slurm walltime/memory kill terminates job()'s process
      directly, so its `except` block never runs and `.completed`/`.failed`
      never gets touched — the task stays "started" forever,
      `poll_simulations_status`'s `left` count never reaches 0 for it, and
      the polling loop could only exit via the unrelated, previously
      silently-swallowed `release_simulations` "no held task" exception (see
      below). Fixed with `reconcile_terminated_tasks` (new): periodically
      (`RECONCILE_INTERVAL_S`, 15 min) cross-checks `.started`-but-unresolved
      tasks against `sacct` (which retains terminal-state history after a
      job leaves `squeue`, unlike `squeue` itself) via the job-ID mapping
      files `job()` already wrote but nothing previously read, and touches
      `.failed` on Slurm's behalf for any task Slurm reports as
      `TIMEOUT`/`OUT_OF_MEMORY`/`CANCELLED`/etc.
- [x] **`release_simulations` false "no held task" exception, silently
      swallowed by the pipeline** — fixed as two related bugs found analyzing
      the same HPC run: (1) once every task has been released at least once,
      zombie/stuck tasks (see above) inflate `status.left` without inflating
      `status.pending`, so the polling loop keeps computing `n > 0` and
      calling `release_simulations` long after there's nothing left to
      release — its own on-disk scan correctly finds nothing, but it still
      unconditionally queried `squeue` and raised once the (finished/aged-out)
      array job no longer appeared there. Fixed: `release_simulations`
      returns immediately when its own scan finds nothing to release, before
      ever touching `squeue`. (2) `experiment_script_runner.run_experiment_script`
      caught *any* exception from a stage (this one included) and printed
      "COMPLETED" regardless, so real failures — this one, and separately a
      genuine `sbatch: Unable to contact slurm controller` outage observed in
      the same run — were invisible until an unrelated downstream stage
      crashed with a confusing secondary error several steps later. Fixed:
      it now re-raises after printing, aborting the pipeline immediately at
      the actual point of failure.
- [x] **Tasks that fail to launch (Slurm auto-requeues them held, e.g.
      squeue reason "launch failed requeued held") were never detected or
      resolved** — fixed: observed live on `--exp-num 2` (4 array tasks
      stuck `PD` with that reason; worked around manually at the time by
      releasing/failing the specific stuck task IDs by hand). Distinct from
      the OOT/OOM case above -- `reconcile_terminated_tasks` only looks at
      tasks that reached `job()` (have a `.started` signal); these never do,
      so there's no job-ID mapping file either (`job()` only writes that
      after `.started`). New `reconcile_launch_failures`, on its own
      2-minute timer (`LAUNCH_FAILURE_RECONCILE_INTERVAL_S`, shorter than
      the OOT/OOM reconciler since a re-release is cheap and prompt recovery
      matters for what's usually a transient node issue): finds
      `.released`-but-never-`.started` tasks, maps array task ID directly
      (1-based) to `sm.get_seeded_simulation_parameters_paths`' order (same
      mapping `submit_simulations`/`job()` already rely on -- no job-ID file
      needed here), matches `squeue`'s free-text `Reason` column against
      "launch failed", and either re-`scontrol release`s it (tracked via an
      on-disk `.launch_retries` counter, up to `LAUNCH_FAILURE_MAX_RETRIES`
      = 3 attempts) or marks `.failed` once retries are exhausted.

- [x] **R_long had no single source of truth, and didn't actually hold
      beta_long == beta_standard** — fixed: `sp['R_long']=1` (a "weekly rate"
      proxy in `standard_values.json`) and `CAL1_R_LONG_WEEKLY_RATE` (a
      second, separate weekly-rate constant in
      `impact_long_shedders_config.py`) were each fed through
      `derive_r_long(rate, tau) = rate * tau / 7.0` -- neither was tied to
      the R the standard population actually runs with, so beta_long ended
      up ~1.4-1.5x beta_standard instead of matching it. New
      `simplicity.settings_manager.compute_r_long(R, tau_2, tau_3,
      tau_3_long, multiplier=1.0)` is the single shared formula:
      `R_long = R * multiplier * (tau_2+tau_3_long) / (tau_2+tau_3)`, so
      beta_long == beta_standard * multiplier exactly. `multiplier` is a new
      per-scenario field on `SCENARIOS` (`beta_multiplier`, defaults 1.0 --
      a deliberate knob for making long shedders more/less infectious per
      day than standard individuals, replacing the old *accidental* gap).
      `R_long` is no longer a stored default anywhere: removed from
      `standard_values.json` and `write_standard_parameters_values()`;
      `read_settings_and_write_simulation_parameters` now computes it via
      `compute_r_long` whenever an experiment doesn't explicitly override
      it. This required also fixing `check_parameters_names`, which
      validated parameter names against `standard_values.json`'s keys (so
      removing `R_long` from there would have broken every experiment that
      explicitly sets it) -- now validates against `parameter_specs.json`
      instead, correctly decoupling "valid parameter name" from "has a
      static default". **Recalibration required**: any experiment run
      before this fix used the old, uncorrected R_long.
      **Superseded same day** — see the entry below: this beta-matching
      design (`R_long` derived from `R`) turned out to make the isolated
      cal_1 context deplete its population before producing usable data
      (run #3/#4), and was replaced with directly-specified per-scenario
      `R_long` values instead of any formula.
- [x] **R_long redesigned again: directly-specified per scenario, no
      formula at all** — after 2.4.14's population/duration bump for cal_1
      still wasn't enough (run #4: only the mildest tau group produced any
      data, with a garbage R²=0.033 fit), traced the root cause to the
      *beta-matching* design itself: holding beta_long == beta_standard
      necessarily makes R_long (total secondary infections) grow with
      duration, so cal_1's 100%-long-shedder isolated population had an
      effective R of 4.6-31 and burned through any population size fast.
      Decided to hold total R_long constant instead (same lifetime
      secondary infections regardless of duration, daily rate drops for
      longer-shedding scenarios) -- and to make it a plain, directly-set
      value per scenario rather than derived from R by any formula at all.
      `SCENARIOS` gained an `R_long` field per entry (replacing
      `beta_multiplier`), defaulting to `1.0` for all four long-shedder
      scenarios. `simplicity.settings_manager.compute_r_long` removed
      entirely (reverted to the pre-2.4.12 shape): `R_long` is a plain
      flat default in `standard_values.json` again (like
      `long_shedders_ratio`), `check_parameters_names` reverted to
      validating against `standard_values.json`, no special-casing left in
      `read_settings_and_write_simulation_parameters`. Consequence: two
      scenarios can now share a duration but have different R_long
      (HIV_low/HIV_high), so cal_1's grid (`unique_long_tau_r_long_pairs`),
      `long_nsr_calibration_plot.py`, `cal_2.py::compute_long_nsr_per_tau`,
      and `lookup_long_nsr` all group/key by the `(tau_3_long, R_long)`
      pair instead of tau alone. Also removed NSR_RANGES' JSON reference
      file (`impact_long_shedders_nsr_ranges.json`) the same way -- both
      NSR ranges and R_long are now plain code in
      `impact_long_shedders_config.py`, no external Data/ file, matching
      that file's own "this is the only file to touch" principle with no
      exceptions left. **Recalibration required.**
- [x] **`infect_long_shedder` recruitment probability decoupled from
      `long_shedders_ratio`** — fixed: `population_model.py:139-146`'s
      per-new-infection recruitment check used a hardcoded `0.01` literal
      instead of `long_shedders_ratio`, so every scenario recruited long
      shedders at the same fixed 1% rate regardless of its configured
      ratio -- coincided for `edge_case` (ratio=0.01) but made `HIV_high`'s
      12% cap largely aspirational. Now
      `population.rng4.uniform() < population.long_shedders_ratio`.
      Verified empirically (20k draws): recruitment rate now tracks ratio
      for both 0.01 and 0.12. Side effect: `CAL1_ISOLATED_FIXED_PARAMS`'s
      `long_shedders_ratio=1.0` ("100% LS" isolated context) now actually
      behaves as advertised, instead of still drawing at 1% per infection.
      **Recalibration required.**

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

- [x] `scripts/` mixed impact_long_shedders-specific files with generic
      utilities — fixed: moved `check_calibration.py` and
      `plot_sot_sanity_regressions.py` into `scripts/experiments/` (both are
      only meaningful to this pipeline); left generic tools
      (`archive_experiment.py`, `check_completed_simulations.py`,
      `fit_OSR.py`, `generate_user_set_parameters_file.py`,
      `get_NSR_for_model_from_OSR.py`, `nextstrain.py`, the `.sh` scripts) in
      `scripts/`. Deleted `scripts/plot_calibration_step_1.py`: its own
      header referenced `calibrate_impact_long_shedders.py`
      (`imp_ls_cal_{scenario}_#{exp_num}`), a script/naming convention that
      no longer exists anywhere in the repo (superseded by the cal_1/cal_2/
      exp pipeline) — genuinely dead code, nothing else referenced it.
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



