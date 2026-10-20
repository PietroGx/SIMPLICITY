#!/usr/bin/env python3
"""Pooling vs per-seed fits, and what a 10% sequencing effort costs.

Two questions, one per clock:

  POOLING  figure 1 fits ONE zero-intercept line through every seed
           concatenated (fig1_preprocess_data.get_model_global_clock /
           get_model_intrahost_clock). Points inside a seed share a genealogy,
           so pooling treats correlated points as independent and the line
           looks far more certain than it is. The two-stage estimator -- one
           slope per seed, then the 30 slopes as the sample -- is what
           output_manager.py:691 already does for the calibration.

  EFFORT   saving keeps the complete diagnosed record; the ratio is a post-hoc
           host-level draw. Whether fitting at 10% costs anything depends on
           which variance dominates:

               Var(per-seed slope at 10%) = Var_between + Var_within
                                            epidemic     sampling

           Var_between is estimated from the full-data per-seed slopes (where
           sampling noise is small); Var_within from repeated draws inside one
           seed. If between dominates, the draw is nearly free and we can fit
           at 10% everywhere for consistency. If within dominates, calibration
           should keep the complete record and only the surveillance
           comparison should draw.

The two clocks are NOT the same measurement and are not treated alike:

  D  global root-to-tip, standard individuals of the control scenario, read
     from sequencing_data_regression.csv. The draw is over diagnosed hosts,
     via sq.subsample -- the surveillance model. Run #3 predates the
     individual_index column, so the host index is joined positionally from
     sequencing_data.csv; that join is test_global_clock_window.seed_points,
     reused here rather than rewritten, and it asserts itself per seed.

  E  intra-host, long shedders, from individuals_data + phylogenetic_data via
     evolutionary_rate.extract_ih_regression_data. It never reads the
     sequencing record, so a surveillance ratio has nothing to bite on. The
     analogous draw is over long-shedder HOSTS, each then followed completely
     -- test_ih_clock_seq_effort.per_host_points, reused here, which validates
     itself against the shared function row for row.

Run from the repo root:
    python tests/test_slope_pooling_and_effort.py --exp-num 3
"""
import os
import sys
import argparse

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))
sys.path.insert(0, os.path.join(HERE, os.pardir, 'scripts', 'experiments'))
sys.path.insert(0, HERE)

import simplicity.dir_manager as dm
import simplicity.sequencing as sq
import simplicity.tuning.evolutionary_rate as er
import impact_long_shedders_config as cfg

from test_global_clock_window import seed_points          # D: the asserted join
from test_ih_clock_seq_effort import per_host_points      # E: validated split

DAYS = 365.25


def slope(x, y):
    """er.tempest_regression -- the zero-intercept fit both panels use."""
    df = pd.DataFrame({'Sequencing_time': np.asarray(x, float),
                       'Distance_from_root': np.asarray(y, float)})
    return float(er.tempest_regression(df).coef_[0])


def decompose(per_seed_full, per_seed_draws, label, target=None):
    """Print the variance decomposition for one clock.

    per_seed_full  : [slope]            one per seed, complete record
    per_seed_draws : [[slope, ...]]     per seed, one per repeated draw
    """
    full = np.asarray(per_seed_full, float)
    draws = np.asarray(per_seed_draws, float)          # (seeds, reps)

    sd_between = full.std(ddof=1)
    sd_within = float(np.mean(draws.std(axis=1, ddof=1)))
    sd_total = draws[:, 0].std(ddof=1)
    frac = sd_between ** 2 / (sd_between ** 2 + sd_within ** 2)

    print(f"\n  --- {label} ---")
    print(f"  seeds                                  {len(full)}")
    print(f"  mean of per-seed slopes   (complete)   {full.mean():.6e}")
    print(f"  mean of per-seed slopes   (10% draw)   {draws.mean():.6e}")
    print(f"  SD BETWEEN seeds          (complete)   {sd_between:.6e}")
    print(f"  SD WITHIN a seed, 10% draw             {sd_within:.6e}")
    print(f"  SD between seeds at 10%   (observed)   {sd_total:.6e}")
    print(f"  between-seed share of variance at 10%  {100 * frac:5.1f}%")
    print(f"  SE of the mean (SD_between/sqrt(n))    "
          f"{sd_between / np.sqrt(len(full)):.6e}")
    if target is not None:
        print(f"  mean vs target {target:.6e}             "
              f"{100 * (full.mean() / target - 1):+.1f}%")
    verdict = ("between-seed dominates -- the 10% draw is nearly free"
               if frac > 0.8 else
               "sampling noise is material -- keep the complete record for "
               "calibration" if frac < 0.5 else
               "mixed -- neither term is negligible")
    print(f"  => {verdict}")
    return sd_between, sd_within


def panel_d(experiment, ratio, reps, rng0):
    """Global clock, control standards. Draw over diagnosed hosts."""
    print("\n=== PANEL D : global clock, control standards ===")
    sods = dm.get_simulation_output_dirs(experiment, group='exp__control')
    if len(sods) != 1:
        print(f"  [warn] {len(sods)} sods in exp__control; using all of them")
    full, draws, pooled_rows = [], [], []
    for sod in sods:
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            reg = seed_points(ssod)                   # asserts its own join
            if reg is None:
                continue
            reg = reg[reg['individual_type'] == 'standard']
            if len(reg) < 30:
                continue
            rows = [{'individual_index': int(i), 't': float(t), 'd': float(d)}
                    for i, t, d in zip(reg['individual_index'],
                                       reg['Sequencing_time'],
                                       reg['Distance_from_root'])]
            full.append(slope([r['t'] for r in rows], [r['d'] for r in rows]))
            pooled_rows.extend(rows)
            per_rep = []
            for rep in range(reps):
                sub = sq.subsample(rows, ratio, seed=rng0 + rep)
                if len(sub) < 10:
                    continue
                per_rep.append(slope([r['t'] for r in sub],
                                     [r['d'] for r in sub]))
            if len(per_rep) == reps:
                draws.append(per_rep)
            else:
                full.pop()
    if not full:
        print("  no usable seeds")
        return
    print(f"  points per seed (complete): "
          f"{len(pooled_rows) // len(full)} median-ish, "
          f"{len(pooled_rows)} pooled")
    pooled = slope([r['t'] for r in pooled_rows], [r['d'] for r in pooled_rows])
    sd_b, _ = decompose(full, draws, "D / control standards")
    print(f"  POOLED slope (what figure 1 draws) {pooled:.6e}")
    print(f"  mean-of-slopes minus pooled        "
          f"{100 * (np.mean(full) / pooled - 1):+.2f}%")
    print(f"  pooled SE if points were independent ~"
          f"{sd_b / np.sqrt(len(pooled_rows)):.6e}  "
          f"(vs the honest {sd_b / np.sqrt(len(full)):.6e})")


def panel_e(experiment, ratio, reps, rng0, scenarios):
    """Intra-host clock, long shedders. Draw over long-shedder hosts."""
    print("\n=== PANEL E : intra-host clock, long shedders ===")
    all_full, all_draws, pooled = [], [], []
    for scenario in scenarios:
        sods = dm.get_simulation_output_dirs(experiment,
                                             group=f'exp__{scenario}')
        if not sods:
            continue
        full, draws = [], []
        for sod in sods:
            for ssod in dm.get_seeded_simulation_output_dirs(sod):
                hosts = per_host_points(ssod)          # validated vs shared fn
                if not hosts:
                    continue
                pts = np.concatenate(hosts)
                if len(pts) < 30:
                    continue
                full.append(slope(pts[:, 0], pts[:, 1]))
                pooled.append(pts)
                per_rep = []
                for rep in range(reps):
                    rng = np.random.default_rng(rng0 + rep)
                    keep = rng.uniform(size=len(hosts)) < ratio
                    kept = [h for h, k in zip(hosts, keep) if k]
                    sub = np.concatenate(kept) if kept else np.empty((0, 2))
                    if len(sub) < 10:
                        continue
                    per_rep.append(slope(sub[:, 0], sub[:, 1]))
                if len(per_rep) == reps:
                    draws.append(per_rep)
                else:
                    full.pop()
                    pooled.pop()
        if full:
            decompose(full, draws, f"E / {scenario}")
            all_full.extend(full)
            all_draws.extend(draws)
    if not all_full:
        print("  no usable seeds")
        return
    decompose(all_full, all_draws, "E / all long-shedder scenarios pooled")
    pts = np.concatenate(pooled)
    print(f"  POOLED slope across cohorts+seeds (figure 1) "
          f"{slope(pts[:, 0], pts[:, 1]):.6e}   n={len(pts)}")
    print(f"  mean of per-seed slopes                      "
          f"{np.mean(all_full):.6e}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    p.add_argument('--ratio', type=float, default=0.10)
    p.add_argument('--reps', type=int, default=20)
    p.add_argument('--seed0', type=int, default=9000)
    p.add_argument('--panels', nargs='+', default=['D', 'E'])
    args = p.parse_args()

    experiment = f"{cfg.exp_name(cfg.DEFAULT_CONSENSUS)}_#{args.exp_num}"
    print(f"experiment : {experiment}")
    print(f"ratio      : {args.ratio}   reps per seed: {args.reps}")

    if 'D' in args.panels:
        panel_d(experiment, args.ratio, args.reps, args.seed0)
    if 'E' in args.panels:
        panel_e(experiment, args.ratio, args.reps, args.seed0,
                [s for s in ('SOT', 'HIV_low', 'HIV_high', 'edge_case')])


if __name__ == '__main__':
    main()
