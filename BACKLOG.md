# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

- [x] **`long_shedders_ratio` was an incidence rate, not a population
      share** — fixed in v2.4.28. Status was drawn per new infection in
      `infect_long_shedder`, so the parameter set the FLOW of long-shedder
      infections; prevalence is flow x residence time, and long shedders
      stay infected 3-5x longer, so the realized share of infected
      individuals came out 2.3-4.1x the nominal ratio (measured
      2.32/4.08/36.3% against a duration-weighted prediction of
      2.78/4.71/40.0%, i.e. 83-91% of prediction across all three
      scenarios). `long_shedders_ratio=0.12` therefore produced a scenario
      in which ~36% of infected individuals were long shedders, peaking at
      58%. The same number was also serving as a concurrent-count cap
      (`population_size * ratio`), which is why incidence read slightly
      BELOW nominal. Now a host trait fixed at creation across all
      reservoir records, reusing the existing `type` field.
      **Note the inflation does not disappear and should not**: with x% of
      the population carrying the trait, more than x% of prevalent
      infections are in those hosts, because they stay infected longer.
      That is now a genuine model output rather than a mis-specified
      parameter. **Recalibration required** (exp #6).
- [x] **Long shedders were censused, not sampled** — fixed in v2.4.28.
      `sequence_long_shedders` sequenced every recovered long shedder, every
      lineage, with no diagnosis requirement and no `sequencing_rate` gate,
      alongside a ~0.5% sample of standard individuals: long shedders were
      1.1% of hosts and 85% of sequences (SOT), 12% and 98% (HIV_high), with
      the census supplying 99-100% of their rows. Both cohorts now sampled
      only on diagnosis at the same rate. Kept on for cal_1 only (100%
      long-shedder population, so no asymmetry possible, and its fit reads
      lineage trajectories not sequences). **Recalibration required.**
- [x] **47% of sequence rows were duplicate genomes** — fixed in v2.4.28.
      `IH_lineages` is a multiset and a duplicated slot carries a
      byte-identical genome until it mutates, so a host sequenced once
      emitted the same sequence several times and was weighted by lineage
      count rather than diversity. This affected standard individuals too.
      Both sequencing paths now emit one row per distinct genome.
      **Recalibration required.**
- [x] **cal_2 keeps only ~12-15% of its seeds, and it is the `min_seq`
      gate on STANDARD sequences** — fixed in v2.4.29 by sequencing every
      diagnosed individual in cal_2 only (`CAL2_SEQUENCING_RATE = 1.0`);
      production and cal_1 keep 0.05. Measured 4/1/3 standard sequences per
      seed at 0.05 (gate passed 0/3) against 50/41/72 at 1.0 (3/3).
      Unbiased: `rng6` drives nothing but the sequencing coin flip and
      recorded sequences are never read back, so the OSR expectation is
      unchanged and only its variance drops; it also removes the
      survivorship bias whereby the gate admitted only seeds with unusually
      large epidemics. Original diagnosis below. Standard individuals are
      sequenced only on diagnosis and only at `sequencing_rate`, giving
      ~0.007 sequences per host, so a seed needs thousands of standard
      infections to clear `min_seq=30`. This is what limited run #5's cal_2
      (24/30/26/**7** points per scenario) and fit HIV_high on 7 points with
      A uncertain to 22.7%. It also makes cal_2 unprobeable at any locally
      runnable scale. Candidate fixes, none applied: raise
      `sequencing_rate` for cal_2 only (observation-side, unbiased, does not
      touch dynamics since `rng6` is dedicated and sequences are never read
      back); raise `CAL2_FINAL_TIME` 365 -> 1095, which would also close the
      calibration/production window gap below; or lower `min_seq`, which
      just fits noisier regressions on the same thin data.
- [ ] **cal_2 calibrates at 365 d while production runs 1095 d** —
      unresolved, introduced in v2.4.25 and first exposed by run #5. Each
      scenario's Stage 2 fit reproduces the 0.0013 target exactly at its
      calibrated NSR *at 365 days*, but production then measured
      0.00209/0.00196/0.00156 (1.61x/1.51x/1.20x). A root-to-tip regression
      through the origin has a window-dependent slope whenever divergence
      accumulates non-linearly, which inherited divergence through
      transmission chains makes it do. Free test, no new simulations:
      refit run #5's sanity data restricted to t <= 1 year; if the standard
      slope drops to ~0.0013 the mismatch is confirmed and the fix is
      `CAL2_FINAL_TIME = 1095`.

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
- [x] **cal_1/cal_2 `IH_virus_emergence_rate` (k_v) mismatch — hypothesis
      tested and ruled out; real cause found and fixed** — real run #1's
      Stage 2 aborted via `NSR_SANITY_MAX` because SOT's standard-NSR fit
      was garbage (R²≈0.0003). Initially traced to a structural difference
      between the two stages (`USER_FIXED_PARAMS['IH_virus_emergence_rate']`
      = 0.1 used by cal_2/exp vs `CAL1_ISOLATED_FIXED_PARAMS` silently
      defaulting to 0) and fixed as a consistency issue regardless
      (`--ih-virus-emergence-rate` now threads through cal_1 too), but the
      dedicated A/B test (`scripts/experiments/run_kv_ab_test.sh`, exp #2
      k_v=0 vs exp #3 k_v=0.1) did NOT confirm k_v as the mechanism:
      control/SOT calibrated near-identically in both arms, and HIV_low
      actually fit *worse* at k_v=0 -- the opposite of what the hypothesis
      predicted. What tracked cleanly with the failure, in both arms, was
      `tau_3_long`, independent of k_v.
      **Real cause**: Stage 1's calibration fit (`long_nsr_calibration_
      plot.py`) was fitting the population-pooled, global-reference,
      absolute-time regression (same shape as the "global clock", just
      row-filtered to long shedders) instead of a genuine intra-host clock
      -- unlike `plot_sot_sanity_regressions.py`'s already-correct
      `load_intrahost_long_data`, which resolves each host's own founder
      genome and own-infection-relative time. Fixed via a new shared
      `simplicity.tuning.evolutionary_rate.extract_ih_regression_data`
      (own founder genome via `inherited_lineage`, time since own
      infection via `ih_birth - t_infection`), used by both cal_1's
      diagnostic fit and (critically) `impact_long_shedders_cal_2.py`'s
      *separate* re-derivation of the same value, which actually reaches
      `nsr_calibration_table.csv` and had NOT been caught by the first
      pass. Root architecture fixed too: cal_1 now persists its fit result
      (`write_calibrated_long_nsr`) and cal_2 reads it back
      (`read_calibrated_long_nsr`, raising if `--target-osr-long` doesn't
      match) instead of independently recomputing -- one implementation,
      not two kept in sync by hand.
      Also: `edge_case` scenario removed (raw-Hamming saturation at its
      long duration, ~0.34 expected substitutions/site, made its Stage 1
      calibration unreliable regardless of methodology) --
      `NSR_RANGES['cal1_long_nsr']['max']` 0.5 -> 0.1. `detect_sod_outliers`'
      IQR pass now off by default (`use_iqr=False`; hard floor unchanged)
      -- too easily thrown off by the single bad seed it was meant to
      catch, at typical 5-20 seed counts. `CAL1_ISOLATED_FIXED_PARAMS[
      'final_time']` scales per group (`CAL1_FINAL_TIME_MULTIPLIER=3 *
      tau_3_long`: SOT≈145d, HIV≈283d) instead of one uniform 3y.
      **Recalibration required**: run the full pipeline under a fresh
      `--exp-num` (next real production run, exp #4).
- [x] **cal_1 swept the standard NSR in a 100%-long-shedder population** —
      found by run #4, the first full run on the corrected intra-host
      measurement. Stage 1's fit came out flat (R²=0.0013/0.0010; measured
      intra-host OSR did not respond to the swept axis at all), so
      inverting at target produced `NSR_long=0` and every production
      scenario ran with long shedders that could not mutate.
      Cause: `mutations.py` applies `rate_long = population.NSR_long` to
      `long_shedder_i` and `rate_standard = NSR` to everyone else, but
      `build_cal1_settings` swept `nucleotide_substitution_rate` while
      `nucleotide_substitution_rate_long` stayed at its
      `standard_values.json` default. The swept axis governed nobody.
      The old population-pooled fit hid this because
      `infect_long_shedder` only fires on new transmissions, so the
      `infected_individuals_at_start` cohort was always standard; those
      seeds did respond to the swept rate, and their divergence reached
      the global-clock regression through inheritance, faking a
      convincing R²=0.98 dose-response.
      Fixed: cal_1 now sweeps `nucleotide_substitution_rate_long` (and
      `long_nsr_calibration_plot.py` reads/fits/plots that parameter);
      new `start_ls` parameter seeds the initial cohort as long shedders
      (set `True` in `CAL1_ISOLATED_FIXED_PARAMS`) so the isolated
      context is genuinely 100% long shedders from t=0;
      `NSR_RANGES['cal1_long_nsr']` retargeted to 1e-5-1e-2 to bracket
      `NSR_long`. Verified locally: seeds come out 10/10 long shedders
      and intra-host OSR now responds (1e-5 -> 0.00124, 1e-3 -> 0.00257,
      bracketing target 0.00205).
      **Recalibration required**: fresh `--exp-num` (exp #5).

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
- [x] **cal_1's isolated calibration: edge_case sequencing shortfall vs.
      SOT depletion, sharing one parameter** — fixed, via `tests/
      test_long_shedder_r_long_grid.py` (new standalone HPC diagnostic
      tool) rather than more local tuning. Root cause of the sequencing
      shortfall turned out to be `CAL1_ISOLATED_FIXED_PARAMS['final_time']`
      itself (365) -- edge_case's `tau_3_long=350.23` is almost as long as
      that window, so individuals barely complete a cycle regardless of
      `R_long` or starting-cohort size. HPC test #1 confirmed the shortfall
      disappears at `final_time=1095` (3 years) even at the default
      `infected_individuals_at_start=10` -- the `infected_individuals_at_start`
      tuning (10->30->40->60) from the same night was chasing an artifact
      of testing at the wrong duration, not a real fix; reverted back to 10.
      Separately, HPC test #2 (all 3 scenarios x R_long in {1.0,1.3}, narrow
      `NSR_RANGES['cal1_long_nsr']`=0.0005-0.002) found the SOT/HIV weak-fit
      question was real, not diagnostic-seed noise: observed OSR never got
      within 3.6-14x of `target_osr_long`=0.00205 anywhere in that range
      (R²=0.002-0.24 everywhere). HPC test #3 (same grid, NSR=0.001-0.1)
      confirmed a much wider range fixes it: R²=0.54-0.92 across all 6
      cells, SOT/HIV fully bracketed, edge_case only slightly short (needs
      NSR≈0.11-0.13 vs sampled max 0.1, small well-constrained
      extrapolation). Also fixed a bug test #2 surfaced: the
      `"HIV_low/HIV_high"` combined-scenario label contained a `/`, read as
      a directory separator by `fit_observed_substitution_rate_regressor`'s
      file-saving, breaking the fit for both HIV cells -- now joined with
      `+`. Final values adopted for the next real pipeline run (exp #1):
      `R_long=1.1` uniformly across all 4 long-shedder scenarios (not the
      1.0/1.3 split from local-only tuning), `NSR_RANGES['cal1_long_nsr']`
      widened to 0.01-0.5, `final_time=1095` (3 years), seeds=20 for every
      stage. **Recalibration required** (as always when R_long/NSR
      ranges/final_time change).

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



