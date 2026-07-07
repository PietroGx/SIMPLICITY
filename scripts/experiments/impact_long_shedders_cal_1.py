#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np

import simplicity.settings_manager as sm
from experiment_script_runner import run_experiment_script
from impact_long_shedders_config import (
    LONG_NSR_EXP_NAME, CAL1_ISOLATED_FIXED_PARAMS, unique_long_taus, derive_r_long,
    read_nsr_ranges, add_slurm_resource_args, set_slurm_resource_env,
)
from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def user_set_experiment_settings(seeds, ranges):
    sp = sm.read_standard_parameters_values()

    def make_settings():
        nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()
        # One group per tau_3_long value (corrected: D - (tau_1+tau_2+tau_4)),
        # each sweeping the same NSR grid. The actual R_long fed to each
        # point is derived from this isolated context's own weekly rate
        # (CAL1_ISOLATED_FIXED_PARAMS['R_long']) scaled by that point's own
        # tau_3_long via derive_r_long -- the same formula every production
        # scenario uses -- rather than a single fixed value shared by every
        # group regardless of duration.
        isolated_r_long_per_week = CAL1_ISOLATED_FIXED_PARAMS["R_long"]
        scenario_groups = [
            {
                'nucleotide_substitution_rate': nsr_values,
                'tau_3_long': tau,
                'R_long': derive_r_long(isolated_r_long_per_week, tau),
            }
            for tau in unique_long_taus(sp)
        ]
        varying_params = {'_scenario_groups': scenario_groups}

        return (varying_params, CAL1_ISOLATED_FIXED_PARAMS.copy(), seeds)

    return make_settings


def run_isolated_calibration(exp_num, runner, seeds, target_osr_long):
    r_long_per_week = CAL1_ISOLATED_FIXED_PARAMS["R_long"]
    print("\n=========================================================")
    print(f" PHASE 1: ISOLATED LONG-SHEDDER CALIBRATION "
         f"(R_long = {r_long_per_week} * tau_3_long/7, 100% LS)")
    print("=========================================================\n")

    exp_name = LONG_NSR_EXP_NAME
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
    add_slurm_resource_args(parser)
    args = parser.parse_args()

    set_slurm_resource_env(args.slurm_mem, args.slurm_time)
    run_isolated_calibration(args.exp_num, args.runner, args.seeds, args.target_osr_long)


if __name__ == "__main__":
    main()
