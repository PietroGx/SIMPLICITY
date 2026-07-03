#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ============================================================================
# Per-scenario STANDARD-only NSR -> OSR calibration plots.
#
# Reads the standard-only OSR CSVs already produced by the Stage 2 sweeps of
# calibrate_impact_long_shedders.py (imp_ls_cal_{scenario}_#{exp_num}),
# re-fits in memory, and plots per scenario:
#   - standard-individual OSR scatter vs NSR
#   - fitted curve
#   - target OSR line
#   - calibrated NSR intersection (star + annotation)
#   - shaded sweep-range guides (to reveal out-of-range calibrations)
#
# No simulations are run. No core files are modified. Output: PNG.
# ============================================================================

import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

SCENARIOS = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"]
PARAMETER = "nucleotide_substitution_rate"


def plot_scenario(scenario_name, exp_num, target_osr_std, model_type,
                  min_seq, min_len, nsr_min, nsr_max):
    """Produce one standard-only NSR->OSR calibration PNG for a scenario."""
    numbered = f"imp_ls_cal_{scenario_name}_#{exp_num}"
    print(f"\n[Plot] {numbered}")

    # Ensure the standard-only CSV exists (idempotent; cached from File 1).
    om.write_OSR_vs_parameter_csv(
        numbered, PARAMETER, min_seq, min_len, individual_type='standard')

    data = om.read_OSR_vs_parameter_csv(
        numbered, PARAMETER, min_seq, min_len, individual_type='standard')

    if data is None or data.empty:
        print(f"   [skip] No standard-only OSR data for {numbered}.")
        return

    # Re-fit in memory. experiment_group tag prevents clobbering File 1's fit CSV.
    fit = er.fit_observed_substitution_rate_regressor(
        numbered, data, model_type,
        parameter_name=PARAMETER,
        experiment_group='standard_plot',
    )

    try:
        calib_nsr = er.compute_calibrated_parameter(model_type, fit, target_osr_std)
    except Exception as e:
        print(f"   [warn] Could not invert calibration: {e}")
        calib_nsr = None

    x = data[PARAMETER].values
    y = data['observed_substitution_rate'].values

    # Dense curve across the observed x-range (extend to calib_nsr if outside).
    x_lo = min(x.min(), nsr_min)
    x_hi = max(x.max(), nsr_max)
    if calib_nsr is not None:
        x_lo = min(x_lo, calib_nsr)
        x_hi = max(x_hi, calib_nsr)
    x_dense = np.linspace(x_lo, x_hi, 300)
    y_dense = fit.eval(x=x_dense)

    fig, ax = plt.subplots(figsize=(8, 5))

    # Sweep-range guides: shade the region actually simulated.
    ax.axvspan(nsr_min, nsr_max, color='green', alpha=0.06,
               label=f'Swept range [{nsr_min:.1e}, {nsr_max:.1e}]')

    ax.scatter(x, y, s=14, alpha=0.4, color='#0173B2',
               label='Standard OSR (single sims)', zorder=1)
    ax.plot(x_dense, y_dense, color='black', linewidth=1.5,
            label=f'{model_type} fit', zorder=2)
    ax.axhline(target_osr_std, color='red', linestyle='--', linewidth=1.2,
               label=f'Target OSR = {target_osr_std}', zorder=3)

    if calib_nsr is not None:
        ax.plot(calib_nsr, target_osr_std, marker='*', markersize=16,
                color='red', markeredgecolor='black', zorder=4)
        in_range = nsr_min <= calib_nsr <= nsr_max
        flag = "" if in_range else "  (OUT OF SWEPT RANGE)"
        ax.annotate(f"calibrated NSR = {calib_nsr:.3e}{flag}",
                    xy=(calib_nsr, target_osr_std),
                    xytext=(0.02, 0.95), textcoords='axes fraction',
                    fontsize=9, color='red' if not in_range else 'black',
                    va='top')

    ax.set_xscale('log')
    ax.set_xlabel("Nucleotide substitution rate (site/year)")
    ax.set_ylabel("Observed Substitution Rate (standard individuals)")
    ax.set_title(f"Standard OSR calibration — {scenario_name}")
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    out_path = os.path.join(
        dm.get_experiment_plots_dir(numbered),
        f"{numbered}_standard_OSR_calibration.png")
    plt.savefig(out_path, format='png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"   saved: {out_path}")
    if calib_nsr is not None and not (nsr_min <= calib_nsr <= nsr_max):
        print(f"   [note] calibrated NSR {calib_nsr:.3e} is OUTSIDE the swept "
              f"range -- widen --min-nsr/--max-nsr for a reliable fit.")


def main():
    parser = argparse.ArgumentParser(
        description="Per-scenario standard-only NSR->OSR calibration plots.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    parser.add_argument('--min-nsr', type=float, default=3e-5,
                        help="Lower bound of the swept range (for guide shading).")
    parser.add_argument('--max-nsr', type=float, default=5e-4,
                        help="Upper bound of the swept range (for guide shading).")
    parser.add_argument('--only', type=str, default=None,
                        help="Plot only this scenario.")
    args = parser.parse_args()

    scenarios = [args.only] if args.only else SCENARIOS
    if args.only and args.only not in SCENARIOS:
        raise ValueError(f"Unknown scenario '{args.only}'. Choices: {SCENARIOS}")

    for name in scenarios:
        try:
            plot_scenario(name, args.exp_num, args.target_osr_std, args.model,
                          args.min_seq, args.min_len, args.min_nsr, args.max_nsr)
        except Exception as e:
            print(f"[error] {name}: {e}")

    print("\n[Done] Plotting complete.")


if __name__ == "__main__":
    main()
