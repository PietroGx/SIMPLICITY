# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

- [ ] **Production configuration set from the calibrated grid: N=3000,
      I0=90, R=1.05, R_long=1.05 — RECALIBRATION REQUIRED** (2026-10-18,
      v2.4.112).

      The pipeline ran at `population_size` 1000, `I0` 50, `R` 1.03 and
      `R_long` 1.1 — numbers no grid ever selected. The calibrated grid
      (#101-#136) measured 36 configurations at their OWN calibrated rates;
      cell #118 is the only candidate with no artefact of either kind. Over
      its 50 runs: phi=1 at a median of day 284, window 817 d, 48/50 reach
      the 1095-day horizon, 0 saturated, 0 that never reach phi=1.

      `R_long` follows `R` rather than staying at 1.1: the grid set them
      equal in every cell, so #118's numbers were measured that way and only
      hold if production does the same.

      `CAL2_FINAL_TIME` 365 -> 730. cal_2 fits the standard rate above the
      burn-in; at 365 days, 78% of its window sits below phi=1 and it fits
      the neutral phase instead.

      **Every calibration and every production run before this is at a
      different configuration.** Rerun both calibration stages.

- [ ] **The individual reservoir is a fixed 100000, not scaled by
      `population_size`** (2026-10-18, v2.4.112). Not fixed, recorded.

      `extrande.check_stop_conditions` ends a run when
      `len(population.reservoir_i) < 1000` (`extrande.py:257`). The reservoir
      starts at `100000 - N` (`population.py:78`, never overridden) and loses
      one individual per diagnosis (`population_model.py:123`) and one per
      recovery (`population.py:460`), so a run stops at
      `diagnosed + recovered = 100000 - N - 1000`, i.e. after
      `(100000 - N - 1000)/N` population turnovers:

          N=1000 -> 98.0      N=3000 -> 32.0      N=5000 -> 18.8

      This is a hard ceiling that tightens as N grows, and it is invisible in
      the decision table: `test_config_grid.measure()` classifies four of the
      five stop conditions and reports the rest as "other". It truncated 9 of
      50 runs in grid cell #134 (N=5000) at days 1005-1092 — the largest
      epidemics, so the truncation is size-correlated, not random. The chosen
      cell uses about 12 of its 32, so it does not bind; it would at N=5000,
      or at any longer `final_time`.

- [ ] **Consensus mode moved to `distribution` everywhere, including both
      calibration stages — RECALIBRATION REQUIRED** (2026-10-15, v2.4.93).

      Until now only the experiment phase set `consensus`. `build_cal1_settings`
      and `build_cal2_settings` never wrote it, so cal_1 and cal_2 ran at
      whatever `standard_values.json` said — `argmax`. The default path was
      calibrate-under-argmax, run-under-distribution.

      Consensus feeds the immune_waning fitness, so it feeds the dynamics and
      therefore the observed substitution rate the calibration fits. Both
      calibration stages now state `consensus` explicitly and default to
      `distribution`, matching the experiment phase.

      **Every calibration produced before this is a calibration of a different
      model configuration.** `nsr_calibration_table.csv` and
      `calibrated_long_nsr.csv` from unbound #8 and earlier do not carry over.
      Recalibrate before the next run.

      Note on why the explicit setting was needed and flipping the code default
      was not enough: `read_standard_parameters_values` returns
      `{**_STANDARD_VALUES_DEFAULTS, **from_file}` — the file wins, and
      `standard_values.json` is machine-local and gitignored. On any machine
      that already had the file, changing the default in `simplicity/` would
      have changed nothing, and the mode a calibration ran under would have
      stayed a property of the host rather than of the pipeline.

      Guards: tests/test_consensus_pipelines.py checks that all three stages
      default to the same mode and that each calibration builder honours an
      override; run_impact_long_shedders_pipeline refuses two `--consensus`
      modes in one run, since one calibration cannot serve both.

- [x] **The analysis path, a calibration fit and a figure all verified in the
      v2.4.71 layout** (2026-10-08).

      Readers: tests/test_analysis_reads_new_layout.py, 24 checks, against a
      real serial run.

      Fit and figure: the real cal_1 pipeline, --exp-num 1 --runner slurm
      --seeds 5. Two groups, 10 NSR_long sweep points each, 50 repeats per
      group. Both fits converged -- SOT R^2 0.960 NSR 0.000720,
      HIV_low+HIV_high R^2 0.981 NSR 0.001143 -- on 50 data points each, so
      every repeat reached the fit. 04_Output holds exactly SOT and
      HIV_low+HIV_high; each group's repeats are numbered 0-49 in its own
      subtree; the calibration plot and calibrated_long_nsr.csv were written.

      The group names reached the fit-results FILENAMES
      (..._exp_fit_results_SOT.csv, ..._exp_fit_results_HIV_low+HIV_high.csv),
      which is exactly where unbound run #1 wrote tau=350.23,R_long=1.1.

      The magnitudes are in the expected range (~1.1e-3 for the corrected
      pipeline), and HIV (tau 94.23) needing a higher NSR_long than SOT
      (tau 48.23) to hit the same target OSR is consistent with longer
      infections saturating more. R^2 is NOT the evidence here -- a good fit
      proves nothing about which quantity was measured.

      NOT a calibration to keep: 5 seeds, run to exercise the machinery.

- [x] **The calibration plot was inspected** (2026-10-08) and it found
      something: the sweep is geometric while the x-axis was linear, so five of
      ten points piled into the leftmost ~8% of the panel -- exactly where the
      target is crossed and the calibrated value is read off. cal_2's plot had
      had plt.xscale('log') all along; the two siblings simply disagreed. Fixed
      in v2.4.75, both now sample the fitted curve with geomspace too.

- [x] **cal_2 and production have run on this code** (2026-10-09). The full
      bound pipeline, --exp-num 2 --cal-seeds 10 --exp-seeds 10: cal_1 (200
      repeats), cal_2 (400), production per scenario, the frozen-table handoffs
      between them, sanity plots, zero failures. The release cap throttled for
      the first time (400 repeats against a cap of 200).

      It also surfaced a SCIENTIFIC finding, not a refactor one, recorded here
      because it is the quantitative case for the design already chosen:
      cal_2's in-context fit degrades monotonically with long-shedder burden
      (R^2 0.884 -> 0.776 -> 0.612 -> 0.290) and for HIV_high it solved at
      NSR 1.4e-7, seven times BELOW the lowest rate swept -- an extrapolation
      off the end of the data, not a calibration. The realised population clock
      then missed its own target anyway (SOT -21%, HIV_low -25%, HIV_high
      +11%), and HIV_high's standard-individual clock was indistinguishable
      from its mixed clock while standard individuals had a mutation rate of
      essentially zero: they were carrying divergence inherited by
      transmission, not generating it.

- [x] **Slurm parsing verified against 26.05.4** (2026-10-08).
      `tests/check_slurm_interface.py` passes on the cluster: squeue still
      emits the header `release_simulations` drops, `sacct --parsable2 -X` is
      one `<job>_<task>|STATE` row per started task, and it retains terminal
      state after the task leaves the queue. It also FOUND a bug --
      `reconcile_launch_failures` could never match a held task, because
      squeue collapses pending tasks sharing a reason into one row with a
      range in the task field (`1-3`). Fixed in v2.4.72.

      A real 4-repeat array then ran end to end: four array positions resolved
      to four distinct repeats, four distinct `main/<index>` map entries, all
      completed, and the polling loop exited on its own. The array mapping --
      the one thing no local test could settle -- is correct.

      NOT exercised: the release cap. The run had
      SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM=200 and only 4 repeats,
      so everything released in one batch and the throttle never engaged.

      No recalibration required. Nothing changes simulation output for the
      same parameters; this is a validation gap, not a correction.

------------------------------------------------------------------------

## Technical Debt

- [x] **Stage 7c (library)** done 2026-10-09, v2.4.77. The nine functions in
      plots_manager and output_manager that list an experiment's simulations
      take `group=None` and thread it, so a multi-group experiment can be read
      one group at a time. The group is in the OSR csv filename too, or two
      groups collide on one path and the second silently reuses the first's
      numbers.

- [x] **Stage 7c (scripts)** done 2026-10-14, v2.4.83. Smaller than the "~35
      call sites" estimate: outside the archive there are 21, and 16 of them
      were already correct. Production experiments are one group each
      (`impact_long_shedders_unbound_exp.py:81` dispatches one experiment per
      scenario, `..._unbound_config.py:328` gives it one named group), so an
      unscoped listing of one IS that group -- now asserted, not assumed, so
      the day that stops being true a test says so instead of a figure.

      The three real ones partitioned a multi-group experiment by reading
      parameters back and bucketing on the float values:
      `long_nsr_calibration_plot.py` on (tau_3_long, R_long);
      `impact_long_shedders_cal_2.py` on (tau_3_long, long_shedders_ratio)
      matched to the frozen table at 3 decimals; `nsr_sensitivity_analysis.py`
      on (tau_3_long, k_v), with groups that had no names at all. All three now
      list one group at a time and take the cell from the record.

      **No recalibration.** tests/test_group_partition.py computes both
      partitions on the real cal_1, cal_2 and sensitivity grids and they are
      the same partition, so every number produced so far stands. What changes
      is what happens when they would differ: two groups sharing a cohort used
      to merge into one bucket and report the union's fit under one of their
      names, and now raise naming both. The test proves it both ways -- with
      group scoping disabled, all four cal_2 scenarios collapse onto one
      blended value and 9 checks fail.

- [x] **unbound cal_2 could never have run** -- found 2026-10-14 by
      tests/test_stage_dispatch.py, fixed in v2.4.83.
      `impact_long_shedders_unbound_cal_2.py` died with
      `KeyError: 'nucleotide_substitution_rate'` four lines before its
      dispatch: 7a moved the sweep axis inside the named group and the
      progress banner still subscripted the old top-level key. Exactly the
      v2.4.26 failure CLAUDE.md records, surviving only because the unbound
      pipeline has never been run. No recalibration -- nothing downstream of
      the banner was reached, so no number was ever wrong; the stage simply
      stopped.

      The guard is the new test: for every stage of both pipelines it stubs
      the runner, calls `main()`, and requires the dispatch to be reached with
      named, non-empty groups. A printout cannot be caught by reading a config.

- [ ] **Stage 8: fold a pipeline back into one experiment**, its stages becoming
      groups. The structure now supports it -- groups are named, numbered from
      0, and each group's subtree is self-contained -- and the frozen
      `nsr_calibration_table.csv` already externalises the stage-to-stage
      handoff. Gated on 7c: until listing is scoped, one experiment holding
      three stages would hand a mixed listing to anything that reads it.

- [ ] Env variables are set in dir_manager (should be set elsewhere, or rename
      dir_manager).

- [ ] Tree builder (lines 43, 78): a tree node stores only the first of the IH
      lineages (for correct tree lineage coloring). If multiple lineages are
      stored, decide on a coloring strategy.
