# SIMPLICITY — Backlog

Work not yet done, nothing else. Findings, measurements and decisions already
taken do not belong here: `CHANGELOG` records what shipped, and the numbers
behind a choice live in the code that implements it.

Every item names the file and the fix. **Recalibration required** means doing
it invalidates output produced before it.

------------------------------------------------------------------------

## Bugs

- **`scripts/archive_paper_runs.sh` ships the patient data.** Line 79 copies
  all of `Data/RealWorldData` into its archive, which on a machine that holds
  `patient_global_clock.csv` (non-public, gitignored) would put it in a file
  meant to be shared. Line 78 collects all of `Data/figures`, which includes
  each experiment's `fig1/` (values derived from that file). Exclude the
  gitignored `patient_global_clock.csv` and every `fig1/`, or collect only
  tracked inputs. The script also still uses the pre-stamp bound/unbound
  numbering.

------------------------------------------------------------------------

## To do

- **Rebuild the studio grid at `eccccec8da18`.** v2.5.16 (censored episodes),
  v2.5.20 (the move to `figures_gui/`), v2.5.24 (`mt_seeded_fix`: seeded
  infections are not transmissions), v2.5.25 (`prom_peak_fix`: the
  prominence bar is a fraction of each peak's own height) and v2.6.1 (the
  unused `keep_clades_by_transmission` deleted from `clustering.py`) all
  change the build hash, so no build on disk is readable. On the HPC:
  `python -I figures_gui/prepare.py --stamp 2 --slurm --pack`, which also
  re-exports `runs/` (v2.5.27: its host counts no longer include the 90
  seeded infections, so the browser's clades match the build's). Then run
  `figures_gui/check_explore.py`, `figures_gui/test_censored_episodes.py --real`
  and `check_app.js` against the studio. **Recalibration required** for every
  figure-2 number taken from a build before v2.5.25. Figure 1 rebuilds itself
  locally when its tab is opened (v2.5.27).

- **`check_completed_simulations` counts a seed valid whose end time is
  unknown.** `check_seeded_simulation_output` tests `os.path.isfile` only, so a
  seed with a zero-byte `final_time.csv` is reported under `Valid simulations:
  N/N` while the same run separately lists it as a skip. The count and the skip
  line disagree about the same seed.

- **Packing an experiment destroys its Slurm logs.**
  `scripts/package_simplicity_data.sh:115` runs `rm -rf "$folder/slurm"` on the
  caller's real directories unless `SIMPLICITY_KEEP_SLURM_LOGS=1`. Documented in
  the file's header and presumably intended, but it means the export step is
  destructive to its input, and it is why the packer must never be exercised
  against `Data/` directly.


- **Recalibrate figures 2, 3 and 4 under rule B.** Every per-clade metric is
  computed on the merged frequency frame, so every panel number moves. Rule B
  folds a refused lineage into its surviving ancestor where the old path
  dropped it, which is why the clade counts barely change and the coverage
  does.

- **Label the composition shares as shares of NON-founding clades when the
  figures are regenerated.** `tests/export_fig2_builder.py:208` (deleted in
  v2.5.20; find where the studio applies the same exclusion) keeps only
  clades with a parent and `n_defining > 0`, which excludes exactly the
  founding clade. The exclusion is intended, the labelling was not: the
  figures reported 2026-10-31 (edge_case 33.7 / 17.7 / 26.9 / 21.7) are
  non-founding shares presented as totals. Re-measured under rule B
  2026-11-18, the excluded burden is 0.0-23.0% of a run depending on the seed,
  not the 21-42% measured under the old partition.

------------------------------------------------------------------------

## Technical debt

- **Env variables are set in `dir_manager`.** Set them elsewhere, or rename
  `dir_manager`.

- **A tree node stores only the first of its intra-host lineages.**
  `simplicity/tree/tree_builder.py:47` takes `row.IH_lineages[0]`, which is what
  makes tree lineage colouring well defined. If several
  are stored, a colouring strategy has to be decided.
