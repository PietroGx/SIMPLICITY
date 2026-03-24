#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import argparse
import itertools
import numpy as np
import time
import concurrent.futures

import simplicity.settings_manager as sm
import simplicity.dir_manager as dm
import simplicity.tuning.diagnosis_rate as dr
import simplicity.tuning.evolutionary_rate as er
from experiments.experiment_script_runner import run_experiment_script

# =============================================================================
# 1. IMPORT CORE PIPELINE MATH
# =============================================================================
from fit_OSR_mixed_pop import (
    calculate_expected_lineage_days,
    fit_calibration_curve,
    calculate_time_fractions
)

# =============================================================================
# 2. USER CONFIGURATION DICTIONARY
# =============================================================================
USER_FIXED_PARAMS = {
    'IH_virus_emergence_rate': 0.01,  
    'infected_individuals_at_start': 100,
    'R': 1.03, 
    'final_time': 365 * 3,
    'sequence_long_shedders': True
}

# =============================================================================
# 2b. CALIBRATION SETTINGS GENERATOR
# =============================================================================
def generate_homogeneous_settings(min_nsr, max_nsr, steps, n_seeds):
    nsr_values = np.logspace(np.log10(min_nsr), np.log10(max_nsr), num=steps)
    nsr_list = [float(x) for x in nsr_values]
    
    def user_set_experiment_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        
        # 1. Start with the exact same base universe as the grid!
        fixed_params = USER_FIXED_PARAMS.copy()
        
        # 2. Force it to be a homogeneous standard population
        fixed_params.update({
            'long_shedders_ratio': 0.0,
            'sequence_long_shedders': False # No long-shedders to sequence
        })
        return (varying_params, fixed_params, n_seeds)
        
    return user_set_experiment_settings

def main():
    parser = argparse.ArgumentParser(description="Grid Manager: Sweeps parameters with calibrated OSR.")
    
    # Required User Inputs
    parser.add_argument('--target-osr', type=float, required=True, help="Empirical mixed OSR to target (e.g., 0.0008)")
    parser.add_argument('--seeds', type=int, required=True, help="Number of seeds per grid point")
    parser.add_argument('--exp-num', type=int, required=True, help="Experiment number")
    parser.add_argument('--runner', type=str, choices=['serial', 'multiprocessing', 'slurm'], default='slurm')
    
    # Sweep Setup
    parser.add_argument('--min-nsr', type=float, default=1e-4, help="Min NSR for baseline sweep")
    parser.add_argument('--max-nsr', type=float, default=1e-2, help="Max NSR for baseline sweep")
    parser.add_argument('--diag-std', type=float, default=0.1, help="Diagnosis rate for standard hosts")
    parser.add_argument('--diag-long', type=float, default=0.1, help="Diagnosis rate for long shedders")

    args = parser.parse_args()
    
    # Grid Definition
    M_values = [1.0, 2.0]
    long_shedders_ratio_values = [0.01, 0.12] 
    tau_3_long_values = [30.0, 60.0, 90.0]
    R_long_values = [1.0, 4.0]
    
    grid = list(itertools.product(M_values, long_shedders_ratio_values, tau_3_long_values, R_long_values))
    print(f"\n[Grid Manager] Generated {len(grid)} valid experimental combinations.")

    # --- PHASE 1: Baseline Sweep ---
    experiment_name = "calibration_run"
    settings_func = generate_homogeneous_settings(args.min_nsr, args.max_nsr, steps=10, n_seeds=100) 
    print(f"\n[Phase 1] Launching homogeneous baseline sweep to establish curve...")
    run_experiment_script(args.runner, args.exp_num, settings_func, experiment_name)
    
    exp_numbered_name = f"{experiment_name}_#{args.exp_num}"
    
    fit_result = fit_calibration_curve(exp_numbered_name, model_type='exp')
    if not fit_result: return
    
    NSR_eff = er.compute_calibrated_parameter('exp', fit_result, args.target_osr)
    
    # --- PHASE 2: Dispatch Grid ---
    print("\n[Phase 2] Calculating rates and dispatching grid experiments...")
    
    sp = sm.read_standard_parameters_values()
    k_v = USER_FIXED_PARAMS.get('IH_virus_emergence_rate', 0.01)
    k_ds = dr.get_k_d_from_diagnosis_rate(args.diag_std, sp['tau_1'], sp['tau_2'], sp['tau_3'], sp['tau_4'])
    
    W_norm = calculate_expected_lineage_days(
        sp['tau_1'], sp['tau_2'], sp['tau_3'], sp['tau_4'], k_ds, k_v, 1, 4
    )
    
    import simplicity.dir_manager as dm
    os.makedirs(dm.get_data_dir(), exist_ok=True)
    log_file = os.path.join(dm.get_data_dir(), f"master_grid_log_#{args.exp_num}.csv")
    
    grid_tasks = [] 
    
    with open(log_file, mode='w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Exp_ID", "M", "long_shedders_ratio", "tau_3_long", "R_long", "W_long", "F_norm", "F_long", "NSR_base", "NSR_long"])

        for idx, (M, long_shedders_ratio, tau_3_long, R_long) in enumerate(grid, start=2):
            
            print(f"Calculating grid point M={M}, ratio={long_shedders_ratio}, tau={tau_3_long}...") # <-- ADD THIS LINE

            grid_exp_name = f"long_shedders_exp_M{M}_lsr{long_shedders_ratio}_tau{tau_3_long}_R{R_long}"
            
            k_dl = dr.get_k_d_from_diagnosis_rate(args.diag_long, sp['tau_1'], sp['tau_2'], tau_3_long, sp['tau_4'])
            W_long = calculate_expected_lineage_days(
                sp['tau_1'], sp['tau_2'], tau_3_long, sp['tau_4'], k_dl, k_v, 5, 15
            )
            
            F_norm, F_long = calculate_time_fractions(long_shedders_ratio, W_norm, W_long)
            NSR_base = NSR_eff / (F_norm + (M * F_long))
            
            writer.writerow([idx, M, long_shedders_ratio, tau_3_long, R_long, round(W_long, 2), round(F_norm, 3), round(F_long, 3), NSR_base, NSR_base * M])
            f.flush() 
            
            def make_settings(nsr=NSR_base, m=M, ratio=long_shedders_ratio, t3l=tau_3_long, rl=R_long, s=args.seeds):
                def grid_settings():
                    varying_params = {}
                    fixed_params = USER_FIXED_PARAMS.copy()
                    fixed_params.update({
                        'nucleotide_substitution_rate': nsr,
                        'M_nsr_long': m,
                        'long_shedders_ratio': ratio,
                        'tau_3_long': t3l,
                        'R_long': rl
                    })
                    return (varying_params, fixed_params, s)
                return grid_settings
            
            grid_tasks.append((args.runner, args.exp_num, make_settings(), grid_exp_name))
            
    print(f"\n[Phase 2] Grid generated! Total experiments: {len(grid_tasks)}. Total jobs: {len(grid_tasks) * args.seeds}.")
    
    # --- DIRECT SUBMISSION LOGIC ---
    print(f"\n[Phase 2] Dispatching {len(grid_tasks)} experiments directly to {args.runner}...")
    
    for i, task in enumerate(grid_tasks, 1):
        runner, exp_num, settings_func, exp_name = task
        print(f"Submitting job {i}/{len(grid_tasks)}: {exp_name}")
        
        # Adding a tiny fraction of a second delay prevents overwhelming the head node's socket/scheduler API
        time.sleep(0.1) 
        run_experiment_script(runner, exp_num, settings_func, exp_name)

    print(f"\n[Success] All {len(grid_tasks)} experiments successfully submitted to the queue!")
    print(f"Mathematical log saved to: {log_file}")

if __name__ == "__main__":
    main()