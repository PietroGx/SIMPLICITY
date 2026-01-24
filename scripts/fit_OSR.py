#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np

from experiments.experiment_script_runner import run_experiment_script
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

def generate_nsr_settings(min_nsr, max_nsr, steps, n_seeds):
    # Log-spaced values for NSR
    nsr_values = np.logspace(np.log10(min_nsr), np.log10(max_nsr), num=steps)
    nsr_list = [float(x) for x in nsr_values]
    
    print(f"\n[Setup] Generated {steps} NSR values from {min_nsr:.2e} to {max_nsr:.2e}")

    def user_set_experiment_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        fixed_params = {
            'infected_individuals_at_start': 100,
            'final_time': 365 * 3,
        }
        return (varying_params, fixed_params, n_seeds)

    return user_set_experiment_settings

def run_fitting(experiment_name, experiment_number, parameter='nucleotide_substitution_rate', 
                      min_seq=30, min_len=100, model_type='exp'):
    print(f"\n[Fitting] Starting OSR fit for: {experiment_name}")
    
    experiment_numbered_name = f'{experiment_name}_#{experiment_number}' 
    
    try:
        # Generate CSVs
        om.write_combined_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len)
        om.write_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len)
        
        # Read Data
        OSR_single = om.read_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len)
        
        if OSR_single.empty:
            print(f"[Warning] No data found for {experiment_numbered_name}.")
            return

        # Fit Model
        print(f'[Fitting] Fitting {model_type} model...')
        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_numbered_name, OSR_single, model_type, weights=None
        )

        # Get Stats
        OSR_combined = om.read_combined_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len)
        OSR_mean = om.get_mean_std_OSR(experiment_numbered_name, parameter, min_seq, min_len)

        # Plot
        print('[Plotting] Saving plot...')
        pm.plot_OSR_fit_figure(
            experiment_numbered_name, fit_result, OSR_single, OSR_combined, OSR_mean,
            model_type, min_seq, min_len
        )
        print("[Fitting] SUCCESS. Plot saved.")

    except Exception as e:
        print(f"[Error] Fitting failed: {e}")

def main():
    parser = argparse.ArgumentParser(description="Simplicity Experiment Manager Pipeline")
    parser.add_argument('--name', type=str, required=True, help="Experiment Name")
    parser.add_argument('--runner', type=str, required=True, help="slurm, serial, or multiprocessing")
    parser.add_argument('--exp-num', type=int, default=1)
    parser.add_argument('--min', type=float, required=True, help="Min NSR")
    parser.add_argument('--max', type=float, required=True, help="Max NSR")
    parser.add_argument('--steps', type=int, default=15)
    parser.add_argument('--seeds', type=int, default=100)

    args = parser.parse_args()

    # 1. Settings
    settings_func = generate_nsr_settings(args.min, args.max, args.steps, args.seeds)

    # 2. Run Simulations
    print(f"\n[Manager] Launching Simulations for {args.name}...")
    run_experiment_script(args.runner, args.exp_num, settings_func, args.name)

    # 3. Run Fitting (Immediately follows because step 2 waited)
    print("\n[Manager] Simulations complete. Starting Fitting...")
    run_fitting(args.name, args.exp_num)

if __name__ == "__main__":
    main()