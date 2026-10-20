#!/usr/bin/env python3
"""What a 10% sequencing effort would mean for cal_1.

cal_1 fits the INTRA-HOST clock, and its points do not come from the
sequencing record at all: extract_ih_regression_data reads individuals_data
+ phylogenetic_data, i.e. every intra-host lineage of every long shedder
(simplicity/tuning/evolutionary_rate.py:62-101). A surveillance ratio applied
to sequencing_data has nothing to bite on there.

It reaches cal_1 in two other ways, and this measures both:

  GATE   long_nsr_calibration_plot.py accepts a seed only when
         len(read_sequencing_data_regression(ssod)) >= min_seq (30). That IS
         the sequencing record, so a 10% effort drops seeds from the fit while
         leaving every surviving seed's points untouched. Reported in
         expectation from the row counts -- cal_1 skips sequencing_data.csv
         (CALIBRATION_SKIPPED_OUTPUTS), so there is no host index to draw on.

  HOSTS  if we decide the ratio should apply to the intra-host record too --
         10% of long-shedder HOSTS, each then followed completely, which is
         what sequencing a long shedder means -- the fit loses hosts. That is
         a drawn subsample, per host, and it is exact.

The per-host split is re-derived here because extract_ih_regression_data
returns a flat frame with no host column. It is checked against the shared
function on every seed, row for row: if the two ever disagree the seed is
reported as DRIFT and excluded, so this file cannot quietly fit something
else than cal_1 does.

om.write_fit_results_csv is stubbed out: the shared fitter writes into the
experiment's Fit_results, and a diagnostic must not overwrite the run's own
calibration output.

Run from the repo root:
    python tests/test_ih_clock_seq_effort.py --exp-num 3
"""
import io
import os
import sys
import argparse
import contextlib

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'scripts', 'experiments'))

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.phenotype.distance as dis
import simplicity.tuning.evolutionary_rate as er
import simplicity.settings_manager as sm
import impact_long_shedders_config as cfg

MIN_SEQ = 30
MIN_LEN = 100


def per_host_points(ssod):
    """[(host_index, DataFrame)] -- the same loop extract_ih_regression_data
    runs, kept host-resolved. Returns None if it does not reproduce the shared
    function exactly."""
    ind_path = os.path.join(ssod, "individuals_data.csv")
    phylo_path = os.path.join(ssod, "phylogenetic_data.csv")
    if not (os.path.isfile(ind_path) and os.path.isfile(phylo_path)):
        return []

    ind_df = om.read_individuals_data(ssod)
    if "type" not in ind_df.columns:
        return []
    ind_df = ind_df[ind_df["type"] == "long_shedder"]
    if ind_df.empty:
        return []

    genome_map = dict(zip(om.read_phylogenetic_data(ssod)["Lineage_name"].astype(str),
                          om.read_phylogenetic_data(ssod)["Genome"]))
    ref_len = len(dis.reference)

    hosts = []
    for host, row in ind_df.iterrows():
        inherited = str(row.get("inherited_lineage"))
        t_inf = row.get("t_infection")
        if inherited not in genome_map or pd.isna(t_inf):
            continue
        traj = row.get("IH_lineages_trajectory")
        if not isinstance(traj, dict):
            continue
        inh_genome = genome_map[inherited]
        pts = []
        for lineage, ev in traj.items():
            lin = str(lineage)
            if lin not in genome_map:
                continue
            ih_birth = ev.get("ih_birth") if isinstance(ev, dict) else None
            if ih_birth is None:
                continue
            pts.append(((float(ih_birth) - float(t_inf)) / 365.25,
                        dis.hamming_iw(genome_map[lin], inh_genome) / ref_len))
        if pts:
            hosts.append(np.asarray(pts, dtype=float))

    shared = er.extract_ih_regression_data(ssod)
    mine = (np.concatenate(hosts) if hosts else np.empty((0, 2)))
    if len(mine) != len(shared) or not np.allclose(
            mine, shared.to_numpy(dtype=float)):
        return None
    return hosts


def osr(points):
    """er.tempest_regression on an (N, 2) array of [time, distance]."""
    return float(er.tempest_regression(pd.DataFrame(
        points, columns=['Sequencing_time', 'Distance_from_root'])).coef_[0])


def draw(hosts, ratio, rng):
    """Keep each HOST with probability `ratio` -- sq.subsample's rule, applied
    to the intra-host record instead of the diagnosis record.
    Returns (points, hosts_kept)."""
    if ratio >= 1.0:
        return np.concatenate(hosts), len(hosts)
    keep = rng.uniform(size=len(hosts)) < ratio
    kept = [h for h, k in zip(hosts, keep) if k]
    return (np.concatenate(kept) if kept else np.empty((0, 2))), len(kept)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exp-num', type=int, required=True)
    p.add_argument('--ratios', type=float, nargs='+', default=[1.0, 0.5, 0.25, 0.10])
    p.add_argument('--reps', type=int, default=10)
    p.add_argument('--target', type=float, default=0.00205)
    p.add_argument('--model-type', default='exp')
    p.add_argument('--seeds-per-point', type=int, default=0,
                   help="0 = every seed.")
    args = p.parse_args()

    om.write_fit_results_csv = lambda *a, **k: None     # never touch the run

    experiment = f"{cfg.exp_name(cfg.DEFAULT_CONSENSUS)}_#{args.exp_num}"
    groups = [g for g in dm.get_groups(experiment) if g.startswith('calib_long__')]
    print(f"experiment : {experiment}")
    print(f"groups     : {groups}")
    print(f"ratios     : {args.ratios}   reps: {args.reps}\n")

    for group in groups:
        print(f"=== {group} ===")
        grid, gate_rows, drift = {}, [], 0
        for sod in dm.get_simulation_output_dirs(experiment, group=group):
            nsr = float(sm.get_parameter_value_from_simulation_output_dir(
                sod, 'nucleotide_substitution_rate_long'))
            ssods = dm.get_seeded_simulation_output_dirs(sod)
            if args.seeds_per_point:
                ssods = ssods[:args.seeds_per_point]
            seeds = []
            for ssod in ssods:
                try:
                    if om.read_final_time(ssod) < MIN_LEN:
                        continue
                    seq_n = len(om.read_sequencing_data_regression(ssod))
                except Exception:
                    continue
                if seq_n < MIN_SEQ:
                    continue
                hosts = per_host_points(ssod)
                if hosts is None:
                    drift += 1
                    continue
                if hosts:
                    seeds.append((hosts, seq_n))
            if seeds:
                grid[nsr] = seeds
                gate_rows.extend(s[1] for s in seeds)

        if drift:
            print(f"  [DRIFT] {drift} seeds did not reproduce "
                  f"extract_ih_regression_data -- excluded")
        if not grid:
            print("  no usable data\n")
            continue

        # --- the gate -------------------------------------------------------
        counts = np.array(gate_rows, dtype=float)
        print(f"  sequencing rows per seed: median {np.median(counts):.0f}, "
              f"min {counts.min():.0f}, max {counts.max():.0f}  (n={len(counts)})")
        for ratio in args.ratios:
            lost = int((counts * ratio < MIN_SEQ).sum())
            print(f"    at {ratio * 100:5.1f}% effort, seeds under the "
                  f"min_seq={MIN_SEQ} gate (in expectation): "
                  f"{lost}/{len(counts)}")

        # --- the hosts ------------------------------------------------------
        print(f"  {'ratio':>6s} {'hosts':>7s} {'points':>8s} {'OSR@top':>10s} "
              f"{'sd':>9s} {'NSR_long':>11s} {'vs 100%':>8s}")
        base_nsr = None
        for ratio in args.ratios:
            reps = 1 if ratio >= 1.0 else args.reps
            nsrs, tops, hostn, ptn = [], [], [], []
            for rep in range(reps):
                rng = np.random.default_rng(1000 + rep)
                rows = []
                kept_hosts = 0
                for nsr, seeds in grid.items():
                    for hosts, _ in seeds:
                        sub, n_kept = draw(hosts, ratio, rng)
                        kept_hosts += n_kept
                        if len(sub) < 2:
                            continue
                        rows.append({'nucleotide_substitution_rate_long': nsr,
                                     'observed_substitution_rate': osr(sub)})
                        if rep == 0 and nsr == max(grid):
                            tops.append(rows[-1]['observed_substitution_rate'])
                hostn.append(kept_hosts)
                df = pd.DataFrame(rows)
                clean = om.detect_sod_outliers(df)
                clean = clean[clean['is_outlier'] == 0]
                try:
                    # the shared fitter prints lmfit's full report per call
                    with contextlib.redirect_stdout(io.StringIO()):
                        fit = er.fit_observed_substitution_rate_regressor(
                            experiment, clean, model_type=args.model_type,
                            parameter_name='nucleotide_substitution_rate_long',
                            experiment_group=group)
                    nsrs.append(float(er.compute_calibrated_parameter(
                        args.model_type, fit, args.target)))
                except Exception as exc:
                    print(f"    {ratio:6.2f}  fit failed: {exc}")
                    continue
                ptn.append(len(clean))
            if not nsrs:
                continue
            mean_nsr = float(np.mean(nsrs))
            if base_nsr is None:
                base_nsr = mean_nsr
            print(f"  {ratio:6.2f} {np.mean(hostn):7.0f} {np.mean(ptn):8.0f} "
                  f"{np.mean(tops) if tops else np.nan:10.6f} "
                  f"{np.std(nsrs):9.6f} {mean_nsr:11.6f} "
                  f"{100 * (mean_nsr / base_nsr - 1):+7.1f}%")
        print()


if __name__ == '__main__':
    main()
