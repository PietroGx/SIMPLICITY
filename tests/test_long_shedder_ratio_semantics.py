# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ============================================================================
# HPC test: what does `long_shedders_ratio` actually control?
#
# population_model.infect_long_shedder uses it for two different things at
# once: a per-new-infection recruitment PROBABILITY, and a concurrent
# PREVALENCE cap (max_long_shedders = population_size * ratio). Those only
# coincide when long shedders clear at the same rate as standard
# individuals, which is exactly what they do not do -- SOT lasts 63 days
# and HIV 109, against ~22 for a standard infection. A recruitment
# probability p should therefore settle at a prevalence well above p, so
# `long_shedders_ratio=0.12` may not mean "12% of infected individuals are
# long shedders". HIV_low vs HIV_high is a pure ratio contrast, so which
# quantity the number denotes changes how that comparison reads.
#
# Measures, per production scenario, over the real production population
# context (USER_FIXED_PARAMS + derive_scenario_params):
#   incidence  -- share of non-founder infections recruited as long
#                 shedders; this is what the probability sets, so it
#                 should track the nominal ratio;
#   prevalence -- time-weighted mean long_shedders/infected from the
#                 trajectory; this is what "12% long shedders" reads as;
#   cap        -- share of simulated time spent at the max_long_shedders
#                 ceiling, i.e. how often the second role binds and
#                 truncates the first.
#
# A large prevalence/nominal ratio, or a non-trivial cap share, means the
# parameter is not doing what its name says. Both are reported per seed and
# aggregated, with the nominal ratio drawn on the figure for reference.
#
# Not part of the impact_long_shedders pipeline -- a standalone diagnostic,
# and it writes nothing the pipeline reads.
#
# NSR: long-shedder recruitment is a coin flip in infect_long_shedder and
# does not depend on the substitution rate, so this defaults to
# standard_values.json's rate rather than a calibrated production one --
# far cheaper (memory tracks cumulative mutation events) and it keeps the
# test decoupled from any particular run's frozen table. The remaining
# NSR-dependence is second-order, through the fitness model feeding back
# into epidemic size; --nsr/--nsr-long are there to confirm the answer
# doesn't move at production rates.
#
# Writes a plain-text recap (Data/tests/test_ls_ratio_out_#{test_n}.txt)
# plus a figure, zipped into Data/{experiment_name}_artifacts.zip.
# ============================================================================

import os
import sys
import zipfile
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt

from simplicity.runme import run_experiment
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm
import simplicity.settings_manager as sm
import simplicity.dir_manager as dm
import simplicity.output_manager as om

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'scripts', 'experiments'))
from impact_long_shedders_config import (
    SCENARIOS, TAU_ROUND, USER_FIXED_PARAMS, derive_scenario_params,
)

EXP_NAME = "long_shedder_ratio_semantics_test"

FINAL_TIME = 365          # long enough for prevalence to equilibrate at tau_3_long=94
DEFAULT_NSR = 0.00008759  # standard_values.json default; see header

RUNNERS = {
    'serial': simplicity.runners.serial,
    'multiprocessing': simplicity.runners.multiprocessing,
    'slurm': simplicity.runners.slurm,
}


def long_scenarios(sp):
    """Frozen params for the long-shedder scenarios, control excluded (its
    ratio is 0, so there is nothing to measure)."""
    return [derive_scenario_params(s, sp) for s in SCENARIOS
            if s["long_shedders_ratio"] > 0.0]


def build_settings(frozen_scenarios, n_seeds, nsr, nsr_long):
    """One group per long-shedder scenario, in the production population
    context. Nothing is swept -- each group is a single parameter point."""
    scenario_groups = [
        {
            'long_shedders_ratio': f["long_shedders_ratio"],
            'tau_3_long': f["tau_3_long"],
            'R_long': f["R_long"],
            'sequence_long_shedders': f["sequence_long_shedders"],
            'nucleotide_substitution_rate_long': nsr_long,
        }
        for f in frozen_scenarios
    ]

    fixed_params = USER_FIXED_PARAMS.copy()
    fixed_params.update({
        'final_time': FINAL_TIME,
        'nucleotide_substitution_rate': nsr,
    })

    def make_settings():
        return ({'_scenario_groups': scenario_groups}, fixed_params, n_seeds)

    return make_settings


def measure_seed(ssod, population_size, ratio):
    """incidence / prevalence / cap-bound share for one seeded simulation,
    or None if its output is unusable."""
    try:
        traj = om.read_simulation_trajectory(ssod)
        ind = om.read_individuals_data(ssod)
    except Exception:
        return None

    if traj.empty or ind.empty:
        return None

    # Incidence: the founder cohort never passes through
    # infect_long_shedder (it only fires on new transmissions), so counting
    # it would dilute the recruitment probability we are trying to read.
    infected_ind = ind[ind["t_infection"].notna()]
    recruited = infected_ind[infected_ind["parent"].astype(str) != "root"]
    if recruited.empty:
        return None
    incidence = float((recruited["type"] == "long_shedder").mean())

    # Prevalence: time-weighted, since Extrande steps are not uniform.
    # Steps with nobody infected carry no ratio and are dropped rather than
    # counted as zero.
    t = traj["time"].to_numpy(dtype=float)
    infected = traj["infected"].to_numpy(dtype=float)
    long_shed = traj["long_shedders"].to_numpy(dtype=float)
    dt = np.diff(t, append=t[-1])
    live = (infected > 0) & (dt > 0)
    if not live.any():
        return None
    weights = dt[live]
    prevalence = float(np.average(long_shed[live] / infected[live], weights=weights))
    peak_prevalence = float((long_shed[live] / infected[live]).max())

    # Cap: how much of that time sat at the max_long_shedders ceiling.
    cap = population_size * ratio
    cap_share = float(weights[long_shed[live] >= cap].sum() / weights.sum())

    return {
        "incidence": incidence,
        "prevalence": prevalence,
        "peak_prevalence": peak_prevalence,
        "cap_share": cap_share,
        "n_recruited": len(recruited),
        "final_time": float(t[-1]),
    }


def analyze_and_report(numbered, frozen_scenarios, test_n, population_size):
    """Aggregate per-scenario incidence/prevalence/cap, write a recap and a
    figure, and zip both into Data/."""
    print(f"--- Analyzing long_shedders_ratio semantics: {numbered} ---")

    simulation_output_dirs = dm.get_simulation_output_dirs(numbered)

    rows = []
    for sod in simulation_output_dirs:
        ratio_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'long_shedders_ratio')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(
            sod, 'tau_3_long')
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            m = measure_seed(ssod, population_size, ratio_val)
            if m is not None:
                m.update({'long_shedders_ratio': ratio_val, 'tau_3_long': tau_val})
                rows.append(m)

    master_df = pd.DataFrame(rows)

    recap = [
        f"=== long_shedders_ratio semantics: {numbered} ===",
        "incidence  = share of non-founder infections recruited as long shedders",
        "prevalence = time-weighted mean long_shedders/infected",
        "cap_share  = share of time at max_long_shedders = population_size * ratio",
        f"population_size={population_size}, final_time={FINAL_TIME}",
        "",
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    x_pos, x_labels = [], []

    for k, frozen in enumerate(frozen_scenarios):
        name = frozen["scenario_name"]
        ratio = frozen["long_shedders_ratio"]
        recap.append(f"--- {name} (nominal ratio={ratio:.3f}, "
                     f"tau_3_long={frozen['tau_3_long']:.2f}) ---")

        df = master_df[
            (master_df['long_shedders_ratio'].round(TAU_ROUND) == round(ratio, TAU_ROUND)) &
            (master_df['tau_3_long'].round(TAU_ROUND) == round(frozen['tau_3_long'], TAU_ROUND))
        ] if not master_df.empty else master_df

        if df.empty:
            recap.append("  no valid data")
            print(f"{name:10s}: no valid data")
            recap.append("")
            continue

        inc, prev = df['incidence'].mean(), df['prevalence'].mean()
        recap.append(f"  n_seeds={len(df)}  mean_final_time={df['final_time'].mean():.1f}")
        recap.append(f"  incidence      = {inc:.4f}  ({inc / ratio:.2f}x nominal)")
        recap.append(f"  prevalence     = {prev:.4f}  ({prev / ratio:.2f}x nominal)"
                     f"  [per-seed sd {df['prevalence'].std():.4f}]")
        recap.append(f"  peak prevalence= {df['peak_prevalence'].mean():.4f}")
        recap.append(f"  cap_share      = {df['cap_share'].mean():.4f}")
        recap.append("")
        print(f"{name:10s}: nominal={ratio:.3f}  incidence={inc:.4f}  "
              f"prevalence={prev:.4f} ({prev / ratio:.2f}x)  "
              f"cap_share={df['cap_share'].mean():.3f}")

        for offset, col, color in ((-0.15, 'incidence', 'tab:blue'),
                                   (0.15, 'prevalence', 'tab:red')):
            ax.scatter(np.full(len(df), k + offset), df[col],
                       color=color, alpha=0.5, s=25,
                       label=col if k == 0 else None)
        ax.hlines(ratio, k - 0.35, k + 0.35, color='black', linestyle='--',
                  linewidth=1.5, label='nominal ratio' if k == 0 else None)
        x_pos.append(k)
        x_labels.append(f"{name}\nratio={ratio:g}")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels)
    ax.set_ylabel("share of infected individuals")
    ax.set_title("long_shedders_ratio: nominal vs realized incidence and prevalence")
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()

    plot_path = os.path.join(dm.get_experiment_plots_dir(numbered),
                             f"{numbered}_ratio_semantics.png")
    plt.savefig(plot_path, dpi=200)
    plt.close()
    print(f"\nPlot saved to: {plot_path}")

    recap_dir = os.path.join(dm.get_data_dir(), "tests")
    os.makedirs(recap_dir, exist_ok=True)
    recap_path = os.path.join(recap_dir, f"test_ls_ratio_out_#{test_n}.txt")
    with open(recap_path, "w") as f:
        f.write("\n".join(recap))
    print(f"Recap written to: {recap_path}")

    zip_path = os.path.join(dm.get_data_dir(), f"{numbered}_artifacts.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(recap_path, arcname=os.path.basename(recap_path))
        zf.write(plot_path, arcname=os.path.basename(plot_path))
    print(f"Artifacts zipped to: {zip_path}")

    return recap_path, plot_path, zip_path


def main():
    parser = argparse.ArgumentParser(
        description="HPC test: measure realized long-shedder incidence and "
                    "prevalence against the nominal long_shedders_ratio, per "
                    "production scenario.")
    parser.add_argument('test_n', type=int,
                        help="Test number (numbers this test's experiment).")
    parser.add_argument('--runner', type=str, default='slurm',
                        choices=['serial', 'multiprocessing', 'slurm'])
    parser.add_argument('--n-seeds', type=int, default=20)
    parser.add_argument('--nsr', type=float, default=DEFAULT_NSR)
    parser.add_argument('--nsr-long', type=float, default=DEFAULT_NSR)
    args = parser.parse_args()

    sp = sm.read_standard_parameters_values()
    frozen_scenarios = long_scenarios(sp)
    population_size = USER_FIXED_PARAMS['population_size']

    experiment_name = f"{EXP_NAME}_#{args.test_n}"
    print("\n=========================================================")
    print(f" long_shedders_ratio semantics test #{args.test_n}")
    for f in frozen_scenarios:
        print(f"  {f['scenario_name']:10s} ratio={f['long_shedders_ratio']:.3f}  "
              f"tau_3_long={f['tau_3_long']:7.2f}  R_long={f['R_long']:.4f}")
    print(f" population_size={population_size}  final_time={FINAL_TIME}  "
          f"NSR={args.nsr}  NSR_long={args.nsr_long}")
    print("=========================================================\n")

    settings_func = build_settings(frozen_scenarios, args.n_seeds,
                                   args.nsr, args.nsr_long)
    run_experiment(experiment_name, settings_func,
                   simplicity_runner=RUNNERS[args.runner], archive_experiment=False)

    analyze_and_report(experiment_name, frozen_scenarios, args.test_n, population_size)


if __name__ == "__main__":
    main()
