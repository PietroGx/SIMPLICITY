# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

- [ ] **No FIGURE has been drawn and no calibration FIT run against real
      output in the v2.4.71 layout.** Narrowed 2026-10-08: the readers are now
      covered. `tests/test_analysis_reads_new_layout.py` runs a real
      experiment with the serial runner and drives the actual analysis
      functions over its output -- every output_manager `read_*`, the SSOD
      path helpers against the group level they gained, the resolver,
      `evolutionary_rate.extract_ih_regression_data`, and the
      write/read_OSR_vs_parameter_csv walk. 23 checks, all passing.

      What is still unproven: `plots_manager` actually rendering a panel, and
      a calibration regression converging on data in this layout. Both need a
      run large enough to support a fit, which the probe experiment
      deliberately is not. The cluster's `slurm_iface_check_#1` tree is in the
      right shape but is far too small for either.

      NOTE for whoever does this: `write_OSR_vs_parameter_csv`
      (output_manager.py:627) wraps its per-repeat body in a bare
      `except Exception: continue`. An empty OSR table therefore looks
      identical whether the cause is thin data or a tree it cannot read. If a
      fit comes back empty, do not infer anything from the table alone --
      check the per-repeat reads directly, as that test does.

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
