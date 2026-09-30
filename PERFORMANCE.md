# Performance — what landed, and what is still open

Measurements are wall clock, means of repeated runs, on the machine this was
developed on. Two workloads are used throughout:

  benchmark   population 500,  NSR 0.006,   R 1.3,  150 days
              -> 1127 lineages, 1012 infections
  production  population 1000, NSR 1.4e-4,  R 1.03, 1095 days
              -> 1072 lineages, 14040 infections

The benchmark runs at ~43x the production substitution rate, so it weights
mutation-related work heavily. Judge changes by the production figure.

## Where it stands

                        benchmark   production   ms/infection
    e673e3e                19.23 s     334.2 s         26.96
    v2.4.41 (9c14c68)       9.30 s     118.3 s          8.43
    v2.4.42 (this)          5.12 s      50.6 s          3.60

v2.4.41 -> v2.4.42 is **2.34x on identical work** (same lineage and infection
counts). e673e3e -> v2.4.42 is not like-for-like: the C1-C9 model changes
altered the dynamics, so that run did ~12% less work. Per infection it is 7.5x.

300 production simulations (30 seeds x 5 scenarios x 2 consensus arms):
27.9 CPU-hours at e673e3e, 9.9 at v2.4.41, **4.2 now**.

## Landed in v2.4.41

1. `consensus.build_weighted_consensus_matrix` -- seed every column with the
   total weight on the reference base, correct only where a genome differs.
   Inner steps 9,866,192 -> 113,168; isolated 2.222 s -> 0.046 s (48.5x).
   Verified numerically identical (max diff 7.7e-12).
2. Distance cache per lineage per consensus epoch. Measured redundancy 153.6x
   (434,714 calls, 2,831 distinct answers). Caches the DISTANCE only; the
   fitness arithmetic still recomputes, since phi and n_act change constantly.
3. `distance.hamming_iw` without the set union: 1.4-1.65x. Identical on 500
   random pairs. Keeps the C7 reference guard.

## Landed in v2.4.42

Verified in two stages. A + D + E are bit-identical: a fixed-seed run produced
**0 of 8 output files differing**. B + C change behaviour, deliberately.

A. **`np.average` length-one guard** -- `phenotype/update.py:78,113`.
   100% of the 434,714 calls average a single value, and `np.average([x])` is
   `x` exactly (0 differences over 300,000 values), so the guard is exact for
   every input. 3.05 us -> 0.04 us. Alone: 8.79 -> 7.38 s.
   NOT `sum(v)/len(v)`: that differs in the last bit for 12% of multi-element
   cases (24,476 of 200,000), because numpy sums pairwise.

B. **`round(..., 4)` removed from `fitness_score`** -- `population_model.py:237,299`.
   The same field was written unrounded from `update.py`, the rounding carried
   no comment and dated from v0.17.00/v2.0.0, and nothing downstream needs 4 dp.
   BEHAVIOUR CHANGE: `individuals_data.csv` differs in `fitness_score` only,
   39 of 1012 rows, max 4.90e-05 -- bounded by the 5e-5 half-step of 4 dp
   rounding. No other column moves, so the trajectory itself is unchanged: the
   rounded value was overwritten by `update_fitness_step` before it could
   select. Removing it also made the length-one guard safe at these two sites,
   since numpy and Python round differently at exact ties (87,099 of 200,000
   tie values differ; random values essentially never expose it).

C. **Fitness trajectory behind `population.track_fitness_traj`** (default False).
   An internal constant, not a simulation parameter. When off,
   `update_fitness_trajectory` returns immediately and `simulation.py` skips the
   save. When on it uses one array conversion for mean and std and
   `np.sum(special.entr(a / a.sum()))` for the entropy -- 68x, and the same
   computation `scipy.stats.entropy` performs internally, without its dispatch
   wrapper or the redundant second normalisation the old code paid for by
   pre-dividing.
   Nothing reads the `Entropy` column; `plots_manager` and the archived figures
   read only Time, Mean and Std. `fitness_trajectory.csv` is therefore now
   optional output, and was removed from `required_files` in
   `check_completed_simulations.py` and from both lists in
   `tests/test_local_runme.py` -- the same treatment that file already gives the
   optional FASTAs.

D. **The output path** -- `population.py:506-531`. It built a frame from all
   100,000 reservoir individuals, `iterrows`-ed over them, wrote each result back
   with `.at[]`, then dropped every susceptible -- 99% of the work. Now it filters
   first, normalises the trajectory dicts in `self.individuals` directly (which
   is what the frame was carrying by reference anyway) and builds the frame once.
   Prototyped at the real proportions: 3.11 s -> 0.02 s. Alone: 8.79 -> 5.88 s.

E. **`expm` identity** -- `tuning/diagnosis_rate.py:71-73`.
   `(e^A)^n = e^(nA)`, so `fractional_matrix_power(expm(B), 1000)` is
   `expm(1000 * B)`. 40.81 ms -> 2.92 ms, agreeing to 4.6e-15, and `k_d` comes
   out identical (0.0055 after 54 steps) for both production diagnosis rates.

## Still open

### `_init_individuals` and the reservoir

`Population.__init__(..., reservoir=100000)` against `population_size=1000`.
Every reservoir individual gets a full record at startup: 0.79 s on the
benchmark, and it scales with `reservoir` rather than with the population. Item
D removed the output-side cost of this; the startup side would need lazy
allocation. The value is a modelling choice, so nothing is proposed.

### `get_k_d_from_diagnosis_rate` is a linear scan

`tuning/diagnosis_rate.py:80` steps from 0.0001 by 0.0001, `max_iter` 10000.
For the production parameters it converges in 54 steps (k_d = 0.0055), so it is
~0.12 s per simulation and not worth changing. Cost is proportional to k_d: a
scenario landing on k_d = 0.5 would take 5,000 steps. Bisection would be ~15
regardless, since the function is monotone in k_d.

## Noted, not a performance issue

`plots_manager.py:101-103` reads the in-memory fitness trajectory as a nested
list (`coord[1][0]`) while `update_fitness_trajectory` appends a dict. That
function would raise on its first call. Stale, and now doubly so: with
`track_fitness_traj` off the list is empty.
