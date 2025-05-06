#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro

If you want to change any parameter, you can specify them in the parameters dictionary below. 

Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""
from experiment_script_runner import run_experiment_script
import argparse
import numpy as np 

experiment_name = 'generate_data_R_diagnosis_rate_parameter_space'

def fixture_experiment_settings():
   
    varying_params = {
        'R': [1, 1.25, 1.5, 2],
        'diagnosis_rate': [round(num, 2) for num in np.arange(0.01, 0.11, 0.01)]
    }

    fixed_params = {}

    n_seeds = 100 

    return varying_params, fixed_params, n_seeds

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to generate IH lineages data")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          fixture_experiment_settings,
                          experiment_name)

if __name__ == "__main__":
    main()