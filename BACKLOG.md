# SIMPLICITY — Backlog

Working list of blockers, active refactors, features, and technical debt.

------------------------------------------------------------------------

## Blockers (before next production rerun)

- [ ] **The analysis and figure path has never run against real output in the
      v2.4.71 layout.** The data flow refactor is verified end to end for
      setup, dispatch and the resolver, and a real serial run produces correct
      output -- but no figure has been drawn and no calibration fit has been
      run from a `04_Output/<group>/sim_NNN__<label>/seed_NNNN/` tree. That is
      40+ `get_parameter_value_from_simulation_output_dir` call sites, all of
      `plots_manager`, both cal_2 fitters and the figure preprocessors. Their
      signatures are unchanged and they compile, which is not the same thing.
      Watch `long_nsr_calibration_plot` hardest: how it gets group names
      changed, and that is verified only at the settings level, never through
      an actual fit.

- [ ] **Slurm parsing is unverified against 26.05.4.** `tests/fakeslurm/`
      covers the control-plane logic and `tests/test_slurm_lifecycle.py`
      passes, but the shims could be lying about output format. Run
      `python tests/check_slurm_interface.py --raw` on the cluster -- three
      no-op tasks, about a minute -- before trusting a real submission. Then
      one small real experiment through the slurm runner: repeat order is now
      settings order where it used to be readdir order, and only a real array
      exercises that mapping.

      No recalibration required for either. Nothing changes simulation output
      for the same parameters; both are validation gaps, not corrections.

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
