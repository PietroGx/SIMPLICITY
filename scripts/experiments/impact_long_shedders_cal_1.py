#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np

import simplicity.settings_manager as sm
from experiment_script_runner import run_experiment_script
from impact_long_shedders_config import unique_long_taus, read_nsr_ranges
from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def user_set_experiment_settings(seeds, ranges):
    sp = sm.read_standard_parameters_values()

    def make_settings():
        varying_params = {
            'nucleotide_substitution_rate': np.geomspace(
                ranges['min'], ranges['max'], ranges['steps']).tolist(),
            'tau_3_long': unique_long_taus(sp),  # corrected: D - (tau_1+tau_2+tau_4)
        }

        fixed_params = {
            'long_shedders_ratio': 1.0,
            'R': 1.0,
            'R_long': 1.5,
            'infected_individuals_at_start': 100,
            'final_time': 1200,
            'sequence_long_shedders': True
        }

        return (varying_params, fixed_params, seeds)

    return make_settings


def run_isolated_calibration(exp_num, runner, seeds, target_osr_long):
    print("\n=========================================================")
    print(" PHASE 1: ISOLATED LONG-SHEDDER CALIBRATION (R_long=1.5, 100% LS)")
    print("=========================================================\n")

    exp_name = "calibrate_long_nsr"
    ranges = read_nsr_ranges()['cal1_long_nsr']
    settings_func = user_set_experiment_settings(seeds, ranges)

    print(f"Dispatching grid experiment: {exp_name}_#{exp_num} to {runner}...")
    run_experiment_script(runner, exp_num, settings_func, exp_name)
    print("\n[Success] Calibration sweep submitted via native runner.")

    numbered = f"{exp_name}_#{exp_num}"
    print(f"\n[Plot] Fitting + plotting long-NSR calibration for {numbered}...")
    plot_and_fit_long_nsr_calibration(numbered, target_osr_long)


def main():
    parser = argparse.ArgumentParser(
        description="Isolated long-shedder NSR calibration (Phase 1 of the "
                    "impact_long_shedders pipeline).")
    parser.add_argument('--exp-num', type=int, required=True,
                        help="Experiment number.")
    parser.add_argument('--runner', type=str,
                        choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    parser.add_argument('--seeds', type=int, default=25,
                        help="Seeds per grid point.")
    parser.add_argument('--target-osr-long', type=float, default=0.00205,
                        help="Target long OSR, used only for the diagnostic "
                            "calibration-fit plot (not for the grid itself).")
    args = parser.parse_args()

    run_isolated_calibration(args.exp_num, args.runner, args.seeds, args.target_osr_long)


if __name__ == "__main__":
    main()
