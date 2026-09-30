# Performance — measured findings, not yet acted on

Measured 2026-09-30 on the post-C1..C9 code, experiment:
population 500, NSR 0.006, R 1.3, R_long 3, long_shedders_ratio 0.10,
final_time 150 (1127 lineages, 1012 infections). Profiles via cProfile;
timings are wall clock, means of two runs.

## Already done (v2.4.41 work, in the tree)

Three optimisations landed with the C1-C9 change set. Same experiment:

    base e673e3e     20.91 / 19.35 s
    new  argmax      11.89 /  9.93 s
    new  distribution 9.03 / 11.29 s

1. `consensus.build_weighted_consensus_matrix` -- seed every column with the
   total weight on the reference base, correct only where a genome differs.
   Inner steps 9,866,192 -> 113,168; isolated 2.222 s -> 0.046 s (48.5x).
   Verified numerically identical (max diff 7.7e-12).
2. Distance cache per lineage per consensus epoch. Measured redundancy 153.6x
   (434,714 calls, 2,831 distinct answers). Caches the DISTANCE only; the
   fitness arithmetic still recomputes, since phi and n_act change constantly.
3. `distance.hamming_iw` without the set union: 1.4-1.65x. Identical on 500
   random pairs. Keeps the C7 reference guard.

## Open, in order of expected value

### 1. `population.py:506-525` -- the output loop processes 100,000 rows to keep ~1,000

    individuals_data = pd.DataFrame(self.individuals).transpose()   # 100,000 rows
    for idx, row in individuals_data.iterrows():
        ...
        individuals_data.at[idx, 'IH_lineages_trajectory'] = lineage_traj_dic

`individuals_data_to_df` (`:530`) then drops every susceptible -- 99% of what was
just processed. The `.at[]` inside a 100,000-iteration loop is the profile's
100,002 `Series.__init__` and 200,031 `sanitize_array` calls, about 3.4 s.

Two independent fixes: filter to non-susceptible BEFORE the loop, and assign the
column once instead of 100,000 element writes.

Output path only -- cannot affect model behaviour. Biggest single win.

### 2. `np.average` on 1-10 element lists -- 435,706 calls, about 3.0 s

    size  1: np.average 3.12 us   sum/len 0.09 us   34x
    size  5: np.average 3.13 us   sum/len 0.11 us   29x
    size 10: np.average 3.30 us   sum/len 0.13 us   26x

Call sites: `phenotype/update.py:78,113`, `population_model.py:237,299`.
Numpy dispatch overhead dwarfs the arithmetic at these sizes.

### 3. `tuning/diagnosis_rate.py:73` -- a matrix power that is an identity away

    B_ex = scipy.linalg.expm(B_aug)
    Bt   = scipy.linalg.fractional_matrix_power(B_ex, 1000)

Since (e^A)^n = e^(nA) for integer n, this is `expm(1000 * B_aug)`:

    current              40.81 ms  -> 0.1001446159
    expm(1000 * B_aug)    2.92 ms  -> 0.1001446159
    max abs difference 4.6e-15, 14x faster

This is what generates the 10,442 `pade_UV_calc` calls -- `fractional_matrix_power`
calls `expm` internally many times per invocation.

Related but separate: `get_k_d_from_diagnosis_rate` (`:80`) is a LINEAR scan from
0.0001 in steps of 0.0001, max_iter 10000. For the parameters tested it converged
in 54 and 5 steps (0.15 s total), but cost is proportional to k_d: k_d = 0.05
would take 500 steps, k_d = 0.5 would take 5,000. Bisection would be ~15
regardless, since the function is monotone in k_d. CHECK WHAT k_d THE PRODUCTION
SCENARIOS LAND ON before deciding whether this matters.

## Checked and ruled out

The intra-host `expm` table works: 193 misses against 300 cached entries per host
model (`intra_host_model.py:96-106`). Not a problem.

## Not investigated

`_init_individuals` (`population.py:187`) is 0.78 s for 100,000 reservoir
individuals -- one-time per simulation, scales with `reservoir` (default 100,000)
rather than `population_size`.
