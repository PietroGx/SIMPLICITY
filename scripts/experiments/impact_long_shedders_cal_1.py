#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from experiment_script_runner import run_experiment_script

def user_set_experiment_settings():
    varying_params = {
        'nucleotide_substitution_rate': np.geomspace(0.0001, 0.002, 6).tolist(),
        'tau_3_long': [63.0, 109.0, 365.0]  # Durations for SOT, HIV, EdgeCase
    }
    
    fixed_params = {
        'M_nsr_long': 1.0,           
        'long_shedders_ratio': 1.0,  
        'R': 1.0,                    
        'R_long': 1.5,
        'infected_individuals_at_start': 100,
        'final_time': 1200,
        'sequence_long_shedders': True
    }
    
    n_seeds = 25
    
    return (varying_params, fixed_params, n_seeds)


def run_isolated_calibration():
    print("\n=========================================================")
    print(" PHASE 1: ISOLATED LONG-SHEDDER CALIBRATION (R_long=1.5, 100% LS)")
    print("=========================================================\n")
    
    runner = "slurm" 
    exp_num = 2
    exp_name = "calibrate_long_nsr"
    
    print(f"Dispatching grid experiment: {exp_name}_#{exp_num} to {runner}...")
    run_experiment_script(runner, exp_num, user_set_experiment_settings, exp_name)
    print("\n[Success] Calibration sweep submitted via native runner.")


if __name__ == "__main__":
    run_isolated_calibration()
