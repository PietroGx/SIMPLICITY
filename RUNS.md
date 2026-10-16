# Run register

One row per `--exp-num`, newest last. Numbering restarts at 1 with the
settings/data-flow refactor: a run is now ONE experiment,
`impact_long_shedders_#N`, holding every phase as a named group, so an old
number does not describe the same thing as a new one.

Fill the row when the run is LAUNCHED, not when it finishes — a row that was
never written is a run nobody can trace back to a commit. Update `status` and
`notes` afterwards.

`commit` is what was checked out at launch. Results are only comparable
between runs whose commit and consensus mode match.

| # | launched | commit | consensus | cal/exp seeds | cal_long | cal_std | status | notes |
|---|----------|--------|-----------|---------------|----------|---------|--------|-------|
| 1 | 2026-10-16 | v2.4.98 | distribution | 30 / 30 | 1e-4 .. 1e-2 | 1e-5 .. 1e-3 | **abandoned** | Died in cal_2's fit after all 300 of its simulations completed — `compute_standard_nsr` named the fit from `dm.get_groups`, which since 8e returns all four calibration groups. Fixed in v2.4.103. Left on disk; rerun as #2 rather than cleared. |
| 2 | 2026-10-16 | v2.4.104 | distribution | 30 / 30 | 1e-4 .. 1e-2 | 1e-5 .. 1e-3 | launched | Rerun after the cal_2 fix, under a fresh number rather than clearing #1. Ranges from probes #1–#4 (`tests/run_nsr_range_finder.sh`, v2.4.97). First calibration under `distribution` — the pre-refactor long-side rates were fitted under `argmax`. First run at the profiled resource defaults (1G / 5h, v2.4.100–102); `--slurm-mem-cal2` still overrides cal_2 to 5G. |

## Pre-refactor runs

Kept for provenance only. They used one experiment directory per stage and per
scenario, and `argmax` consensus; their calibration tables do not transfer.

| run | layout | notes |
|-----|--------|-------|
| bound #1–#7 | `impact_long_shedders_<scenario>_#N` | Standard NSR fitted per scenario inside the mixed population. #7 drove HIV_high to 9.3e-7, 150× below its own control. Pipeline deleted in v2.4.90. |
| unbound #8 | `impact_long_shedders_unbound_*_#8` | First end-to-end unbound run. Calibration exact on re-fit; `argmax`. |
