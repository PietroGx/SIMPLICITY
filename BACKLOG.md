# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

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
