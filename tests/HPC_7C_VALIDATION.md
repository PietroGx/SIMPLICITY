# Validating stage 7c

Cheapest first. Each step says what it proves and what a failure means.

## Which tier runs where

| | local tier | HPC tier |
|---|---|---|
| script | `tests/run_tests_local.sh` | `tests/run_tests_on_compute.sh` |
| where | your own machine | a cluster **compute** node, via `sbatch` |
| needs | a python env | a real Slurm, or a real run's `Data/` |
| runners exercised | serial, multiprocessing, faked Slurm | real Slurm |

The ten local tests build their own experiments and simulate for minutes. A
compute node makes none of them truer, so they do not belong on the cluster —
and a login node is the one place they must never run. `run_tests_local.sh`
refuses to start if it finds `sbatch` on `PATH` and no `SLURM_JOB_ID`.

Run everything from the **repo root**. A script run from inside `scripts/` or
`tests/` leaves a stray `Data/` the root `.gitignore` does not cover.

---

## 1. The local tier — on your machine, ~5 min

```bash
bash tests/run_tests_local.sh
```

`--quick` drops the two that simulate end to end (`analysis reads layout`,
`multiprocessing runner`) and returns in about 40 s; `--list` prints what would
run.

Expect `10 passed, 0 failed`. The two added for 7c are:

- **group partition (7c)** — 50 checks. Builds the real cal_1, cal_2, production
  and sensitivity grids through the configs' own builders, and asserts the old
  value-derived partition and the new group partition are the **same partition**.
- **stage dispatch** — calls `main()` on every stage of both pipelines with the
  runner stubbed. This is the one that found the unbound cal_2 bug.

A failure here is a logic problem and nothing on the cluster should be run until
it is understood.

---

## 2. Get the branch on the cluster

```bash
cd ~/SIMPLICITY
git fetch origin
git checkout settings_data_flow
git pull --ff-only origin settings_data_flow
git log --oneline -1
```

`--ff-only` refuses to invent a merge commit. If it errors, the cluster copy has
local commits — say so rather than forcing it.

---

## 3. Does this cluster's Slurm still parse — ~2 min

```bash
bash tests/run_tests_on_compute.sh
```

With no flags the HPC tier is just this probe. It submits itself with `sbatch`
and tails the log; Ctrl-C stops the tailing, not the job.

`tests/check_slurm_interface.py` submits a small held array and checks the three
replies `simplicity/runners/slurm.py` parses, against the controller that
actually runs the work. The cluster is on Slurm 26.05.4 and the newest Slurm
packaged for Ubuntu 24.04 is 23.11.4, so `test_slurm_lifecycle`'s fake can
confirm the *logic* and say nothing about the *format*. This is the gap it fills.

---

## 4. The proof that 7c changed no number — ~5 min

This is the step that matters. Step 1 runs on grids the test built itself; this
one re-fits **run #2's real output** with the refactored code and diffs it
against the numbers that run froze.

```bash
bash tests/run_tests_on_compute.sh --no-interface --calib 2
```

(It takes the target and model from the run's own saved CSV, so there is
nothing to remember — except
`--target-osr-std`, which the frozen table does not record; it defaults to the
pipeline's 0.0013 and is printed. Drop `--no-interface` to re-run step 3 too.)

It copies the experiment's `Fit_results/` and `05_Plots/` aside, recomputes, and
**puts them back** before exiting. `--keep` leaves the recomputed ones instead.

Expect, per cohort and per scenario:

```
  [ok ] tau=48.23 R=1.1        was 0.000123   now 0.000123   (+0.0000%)
  ...
the refactored code reproduces this run's calibration exactly.
No recalibration is owed.
```

The tolerance is `1e-9` on purpose: a repartition that still reads the same rows
should reproduce the fit bit for bit.

**A difference here is the finding**, not a tolerance to widen. It would mean the
old bucketing and the group assignment disagreed on real data — in which case
run #2's frozen table was fitted on the wrong rows and the recalibration
question is open. Send the output; do not adjust `--rtol`.

`SKIPPED` on both stages means it compared nothing (no finished run under that
`--exp-num` on this machine) and proves nothing. Point it at a run that
completed.

---

## 5. The paper figures, against real output — ~15 min

```bash
bash tests/run_tests_on_compute.sh --no-interface --figures 2
```

Never yet run on real data in its current form: it was rewritten after four
defects were found in it, including one that reported a returned path as
"139 item(s)". It now cross-checks every preprocessor against `settings.json`,
`repeats.json` and the directories on disk, and names the difference.

The two `long_paper` entries read a **different experiment family** — the grid,
`master_grid_log_#N.csv` and `calibration_run_#N`. Against an impact-pipeline
run they are gated off and report `n/a`, each on the file it actually opens.
`n/a` there is correct, not a skip hiding a failure.

---

## 6. Resume and reconciliation, against real Slurm — ~20 min

```bash
bash tests/run_tests_on_compute.sh --no-interface --slurm 95    # any FREE exp-num
```

Pick an `--exp-num` nothing else uses; it writes `slurm_resume_check_#N` and
`slurm_kill_check_#N` and leaves them for inspection.

Phase 3 deliberately lets Slurm kill tasks over a 2-minute walltime, so the
reconciler has something to reconcile. It previously passed while testing
nothing — a 3000-day horizon against 200 individuals finished in twelve seconds,
because `final_time` is *simulated* days. It is now sized by work (4000
individuals) and reports **INCONCLUSIVE with a non-zero exit** if nothing is
killed, rather than going green. INCONCLUSIVE is not a pass: raise the
population or lower `SIMPLICITY_SLURM_TIME` and run it again.

---

## 7. The unbound pipeline, end to end — hours

The pipeline you selected, and the only thing that exercises `edge_case` naming
all the way through. It has **never been run** — step 1's dispatch test is what
found that its stage 2 died on a stale printout before submitting anything
(fixed in v2.4.83), so expect to be the first person through this path.

In tmux, never `nohup`:

```bash
bash scripts/start_simplicity_session.sh

# inside the session, from the repo root:
python scripts/experiments/run_impact_long_shedders_pipeline.py \
    --exp-num <FREE> --runner slurm --cal-seeds 10 --exp-seeds 10
```

Detach with `Ctrl-b d`, reattach with `tmux attach -t simplicity`.

Watch for, in order:

1. **stage 1** dispatches `impact_long_shedders_cal_long_#N` with three
   groups: `SOT`, `HIV_low+HIV_high`, `edge_case`. If `edge_case` appears as
   `tau_350.23_Rl_1.1`, the group is not naming itself and 7a regressed.
2. **stage 2** prints `standard_only  NSR grid: 10 pts [...]` and then
   dispatches. That banner is the line that was broken; it reaching the
   dispatch at all is the fix working.
3. **stage 3** dispatches TWO production families, five scenarios each:
   `impact_long_shedders_<scenario>_#N` (argmax) and
   `impact_long_shedders_dist_<scenario>_#N` (distribution), with a
   sanity-plot job per family.
4. `Fit_results/` under cal_long holds one file per group, named by the group —
   `..._exp_fit_results_edge_case.csv`, not a float tuple — and
   `..._calibrated_long_nsr.csv` now carries a `group` column naming which
   cohort each row came from.

When it finishes, run step 4 against it:

```bash
bash tests/run_tests_on_compute.sh --no-interface --calib <N>
```

Read its stage 2 for what it is. The pipeline fits ONE standard NSR
from the `standard_only` group and applies it to every scenario
(`compute_standard_nsr`, impact_long_shedders_cal_2.py:55), so the five
`[ok]` lines are one number printed five times, not five independent
comparisons. The evidence is the single fit above them. Stage 1, with its three
cohorts, is where that stage carries five numbers' worth of information.

---

## If you want the local tier on the cluster anyway

```bash
bash tests/run_tests_on_compute.sh --with-local
```

This is not part of the validation above. Its only use is confirming that the
cluster's own python and filesystem have not broken something the development
machine cannot see — a different numpy, a `/localscratch` quirk. It asks for 2
CPUs because the multiprocessing test uses a pool of 2.

---

## What is still open after all seven

- **Stage 8** — folding a pipeline back into one experiment, its stages becoming
  groups. The 16 call sites that list a production experiment unscoped are
  correct only while each production scenario is its own experiment; step 1
  asserts that, and stage 8 is what ends it. Those 16 are stage 8's first task,
  not 7c debt.
- **The archive** (`scripts/00_SIMPLICITY_paper_archive/`, 9 sites) is
  deliberately not migrated. Three of its files carry committed merge-conflict
  markers from 98d43d8 and do not compile.
