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

- [ ] **Nobody has looked at the calibration plot.** A png was rendered; that
      is not the same as it being right, and this repo's figure history is bad
      (figure 1's panel G fitting 19 sequences, figure 3's C/D comparing 1,494
      against 6). Open
      05_Plots/..._long_nsr_calibration_fit.png and check: two series, legend
      reading SOT and HIV_low+HIV_high, 50 scatter points each.

- [ ] **cal_2 and production have never run on this code.** cal_1 is one stage
      of three. The frozen nsr_calibration_table.csv handoff from cal_1 to
      cal_2, and cal_2's own fit, are untested in this layout. Same command
      shape: scripts/experiments/impact_long_shedders_cal_2.py.

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

- [ ] **Stage 7c of the data flow refactor: ~60 `get_simulation_output_dirs`
      call sites still pass no group.** The listing API takes one
      (`group=None` returns everything, which is what they all get today), so
      they can migrate one at a time. Blocked on having real output to verify
      each against. Until they do, a listing over an experiment with several
      groups -- which cal_1 and cal_2 already produce -- mixes them, and
      `plots_manager` sorts by parameter value across the lot. Making `group`
      required is the end state.

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
