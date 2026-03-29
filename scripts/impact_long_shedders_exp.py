#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import pandas as pd

# Native SIMPLICITY imports
import simplicity.settings_manager as sm
import simplicity.tuning.diagnosis_rate as dr
import simplicity.tuning.evolutionary_rate as er
import simplicity.plots_manager as pm  # <-- Added correct import for plots
from experiments.experiment_script_runner import run_experiment_script

from fit_OSR_mixed_pop import (
    calculate_expected_lineage_days, 
    fit_calibration_curve, 
    calculate_time_fractions
)

# =============================================================================
# 1. CONFIGURATION & CONSTANTS
# =============================================================================
EXP_NAME = "impact_long_shedders"

# Experiment Overrides (Defaults - will be updated by argparse)
USER_FIXED_PARAMS = {
    "population_size": 1000,
    "infected_individuals_at_start": 100,
    "R": 1.05,
    "final_time": 1095,
    "IH_virus_emergence_rate": 0.1,
}

# =============================================================================
# 2. HELPER FUNCTIONS
# =============================================================================
def generate_homogeneous_settings(min_nsr, max_nsr, steps, n_seeds):
    """Generates the parameter grid for the baseline calibration sweep."""
    nsr_values = np.logspace(np.log10(min_nsr), np.log10(max_nsr), steps)
    nsr_list = [float(x) for x in nsr_values]
    
    def make_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        fixed_params = USER_FIXED_PARAMS.copy()
        fixed_params.update({
            'long_shedders_ratio': 0.0,
            'sequence_long_shedders': False,
            'M_nsr_long': 1.0,
        })
        return (varying_params, fixed_params, n_seeds)
    return make_settings


def dispatch_scenario(scenario_name, duration_D, ratio, r_freq, nsr_eff, sp, args):
    """
    Handles all math derivations, normalizations, and Slurm dispatch for a single scenario.
    Because the Slurm runner loops internally, this function call will BLOCK until the job finishes.
    """
    print(f"\n-> Preparing Scenario: {scenario_name}...")
    
    tau_1, tau_2, tau_3, tau_4 = sp['tau_1'], sp['tau_2'], sp['tau_3'], sp['tau_4']
    k_v = USER_FIXED_PARAMS['IH_virus_emergence_rate']
    
    # Calculate evolutionary multiplier based on user arguments
    m_nsr_long = args.target_osr_long / args.target_osr_std
    
    # 1. Math Derivations based on cohort (Control vs. Long-shedders)
    if ratio == 0.0:
        # Control Logic
        tau_3_long = sp['tau_3_long']
        R_long = sp['R_long']
        adjusted_nsr = nsr_eff
        seq_long = False
    else:
        # Long-shedder Logic
        tau_3_long = duration_D - (tau_1 + tau_2 + tau_4)
        R_long = r_freq * (tau_3_long / 7.0)
        
        # Calculate kinetic diagnosis rates (used only for lineage days math)
        k_ds = dr.get_k_d_from_diagnosis_rate(sp['diagnosis_rate_standard'], tau_1, tau_2, tau_3, tau_4)
        k_dl = dr.get_k_d_from_diagnosis_rate(sp['diagnosis_rate_long'], tau_1, tau_2, tau_3_long, tau_4)
        
        # Calculate lineage days using the rigorous equation
        W_norm = calculate_expected_lineage_days(tau_1, tau_2, tau_3, tau_4, k_ds, k_v, 1, 4)
        W_long = calculate_expected_lineage_days(tau_1, tau_2, tau_3_long, tau_4, k_dl, k_v, 5, 15)
        
        # Calculate time fractions and normalize NSR
        F_norm, F_long = calculate_time_fractions(ratio, W_norm, W_long)
        adjusted_nsr = nsr_eff / (F_norm + (m_nsr_long * F_long))
        seq_long = True

    # 2. Build Settings Dictionary
    scenario_params = USER_FIXED_PARAMS.copy()
    scenario_params.update({
        "long_shedders_ratio": ratio,
        "tau_3_long": tau_3_long,
        "R_long": R_long,
        "nucleotide_substitution_rate": adjusted_nsr,
        "M_nsr_long": m_nsr_long,
        "sequence_long_shedders": seq_long
    })

    # 3. Create Settings Factory & Dispatch
    def make_settings():
        return ({}, scenario_params.copy(), args.seeds)

    print(f"   [Dispatching] {EXP_NAME}_{scenario_name} -> Waiting for {args.runner} completion...")
    run_experiment_script(args.runner, args.exp_num, make_settings, f"{EXP_NAME}_{scenario_name}")
    
    # Return the dictionary for logging
    scenario_params['scenario_name'] = scenario_name
    return scenario_params


# =============================================================================
# 3. MAIN EXECUTION PIPELINE
# =============================================================================
def main():
    parser = argparse.ArgumentParser(description="Grid Manager: Sweeps parameters with calibrated OSR.")
    
    # Required User Inputs
    parser.add_argument('--target-osr-std', type=float, required=True, help="Empirical standard OSR to target (e.g., 0.001306)")
    parser.add_argument('--target-osr-long', type=float, required=True, help="Empirical long-shedders OSR to target (e.g., 0.002053)")
    parser.add_argument('--seeds', type=int, required=True, help="Number of seeds per scenario")
    parser.add_argument('--exp-num', type=int, required=True, help="Experiment number")
    parser.add_argument('--runner', type=str, choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    
    # Population Inputs
    parser.add_argument('--pop-size', type=int, default=1000, help="Total population size")
    parser.add_argument('--start-infections', type=int, default=100, help="Initial number of infected individuals")

    # Sweep Setup
    parser.add_argument('--min-nsr', type=float, default=3e-5, help="Min NSR for baseline sweep")
    parser.add_argument('--max-nsr', type=float, default=3e-4, help="Max NSR for baseline sweep")

    args = parser.parse_args()

    # Dynamically update the fixed parameters with user arguments
    USER_FIXED_PARAMS["population_size"] = args.pop_size
    USER_FIXED_PARAMS["infected_individuals_at_start"] = args.start_infections

    # Create setup directory tagged with the experiment number
    setup_dir = f"Data/{EXP_NAME}_setup_data_#{args.exp_num}"
    os.makedirs(setup_dir, exist_ok=True)

    m_nsr_long = args.target_osr_long / args.target_osr_std

    # Load dynamic background parameters from SIMPLICITY core
    sp = sm.read_standard_parameters_values()
    
    # --- Phase 1: Baseline Sweep & Calibration ---
    experiment_name = f"{EXP_NAME}_calibration_run"
    print(f"\n[Phase 1] Launching homogeneous baseline sweep to establish curve...")
    
    settings_func = generate_homogeneous_settings(
        min_nsr=args.min_nsr, 
        max_nsr=args.max_nsr, 
        steps=10, 
        n_seeds=args.seeds
    ) 
    run_experiment_script(args.runner, args.exp_num, settings_func, experiment_name) 
    
    exp_numbered_name = f"{experiment_name}_#{args.exp_num}"
    print("\n[Phase 1] Baseline sweep finished. Fitting curve...")
    
    fit_result = fit_calibration_curve(exp_numbered_name, model_type='exp')
    if not fit_result:
        print("Error: Calibration fit failed. Exiting script.")
        return
        
    nsr_eff = er.compute_calibrated_parameter('exp', fit_result, args.target_osr_std)
    
    with open(f"{setup_dir}/calibration_results.txt", "w") as f:
        f.write(f"Target OSR Standard: {args.target_osr_std}\n")
        f.write(f"Target OSR Long: {args.target_osr_long}\n")
        f.write(f"Calibrated Base NSR (nsr_eff): {nsr_eff}\n")
        f.write(f"M_nsr_long used: {m_nsr_long}\n")
        f.write(f"Population Size: {args.pop_size}\n")
        f.write(f"Starting Infections: {args.start_infections}\n")

    print("\n[Phase 1] Generating Calibration Plots (Tempest grids and Curve Fit)...")
    try:
        # 1. Plot the grid of root-to-tip tempest regressions
        pm.plot_combined_tempest_regressions(
            experiment_name=exp_numbered_name, 
            parameter='nucleotide_substitution_rate',
            individual_type='standard'
        )
        
        # 2. Plot the 3-panel OSR vs NSR calibration curve fit
        pm.plot_OSR_fit(
            experiment_name=exp_numbered_name,
            fit_result=fit_result,
            model_type='exp',
            min_seq_number=30,    
            min_sim_lenght=300    
        )
        print(f"Plots successfully generated and saved to the experiment plots directory.")
    except Exception as e:
        print(f"Warning: Could not generate calibration plots. Error: {e}")

    # --- Phase 2 & 3: Definition and Sequential Submission ---
    print("\n[Phase 2 & 3] Dispatching Scenarios...")
    
    scenarios = [
        # (Scenario Name, Total Duration D, Prevalence Ratio, Freq for R_long)
        ("control",  19.41, 0.00, None),
        ("SOT",      63.00, 0.01, 0.5),
        ("HIV_low",  109.0, 0.01, 0.5),
        ("HIV_high", 109.0, 0.12, 1.0),
        ("edge_case",365.0, 0.01, 0.5)
    ]
    
    logged_data = []
    
    for name, D, ratio, freq in scenarios:
        result_dict = dispatch_scenario(name, D, ratio, freq, nsr_eff, sp, args)
        logged_data.append(result_dict)

    # Save Parameter Table
    df = pd.DataFrame(logged_data)
    cols = ['scenario_name'] + [c for c in df if c != 'scenario_name']
    df[cols].to_csv(f"{setup_dir}/parameter_table.csv", index=False)
    print(f"\n[Success] All scenarios sequentially completed. Data saved to {setup_dir}/parameter_table.csv")

if __name__ == "__main__":
    main()