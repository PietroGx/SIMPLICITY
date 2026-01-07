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

def run_fitting(experiment_name, parameter='nucleotide_substitution_rate', 
                      min_seq=30, min_len=100, model_type='exp'):
    
    # Define the target host type 
    HOST_TYPE = 'normal'
    
    print(f"\n[Fitting] Starting OSR fit for: {experiment_name} (Type: {HOST_TYPE})")

    try:
        # Generate CSVs (Pass target_host_type)
        om.write_combined_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len, 
                                               target_host_type=HOST_TYPE)
        om.write_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len, 
                                      target_host_type=HOST_TYPE)
        
        # Read Data (Pass target_host_type)
        OSR_single = om.read_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len, 
                                                  target_host_type=HOST_TYPE)
        
        if OSR_single.empty:
            print(f"[Warning] No data found for {experiment_name} (Type: {HOST_TYPE}).")
            return

        # Fit Model
        print(f'[Fitting] Fitting {model_type} model...')
        # Note: We must also update write_fit_results_csv inside er if you updated that module. 
        # Assuming er.fit_... returns a result object we can save.
        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_name, OSR_single, model_type, weights=None
        )
        
        # Save fit results with the specific type tag (Updated logic from previous prompt)
        om.write_fit_results_csv(experiment_name, model_type, fit_result, target_host_type=HOST_TYPE)

        # Get Stats (Pass target_host_type)
        OSR_combined = om.read_combined_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len,
                                                             target_host_type=HOST_TYPE)
        OSR_mean = om.get_mean_std_OSR(experiment_name, parameter, min_seq, min_len,
                                       target_host_type=HOST_TYPE)

        # Plot
        print('[Plotting] Saving plot...')
        # You might want to pass the type to the plot function if you updated 
        # plot_OSR_fit_figure to handle titles/filenames differently.
        # Otherwise, it will save using the standard naming convention.
        pm.plot_OSR_fit_figure(
            experiment_name, fit_result, OSR_single, OSR_combined, OSR_mean,
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
    parser.add_argument('--seeds', type=int, default=10)

    args = parser.parse_args()

    # 1. Settings
    settings_func = generate_nsr_settings(args.min, args.max, args.steps, args.seeds)

    # 2. Run Simulations
    print(f"\n[Manager] Launching Simulations for {args.name}...")
    run_experiment_script(args.runner, args.exp_num, settings_func, args.name)

    # 3. Run Fitting (Immediately follows because step 2 waited)
    print("\n[Manager] Simulations complete. Starting Fitting...")
    run_fitting(args.name)

if __name__ == "__main__":
    main()