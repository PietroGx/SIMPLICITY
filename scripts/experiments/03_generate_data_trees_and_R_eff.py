 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro
   
Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""
from experiment_script_runner import run_experiment_script
import argparse

experiment_name =  'generate_data_trees_and_R_eff'

def fixture_experiment_settings():
    
    varying_params = {
        'phenotype_model': ['distance from wt', 'immune waning']
    }

    fixed_params = {
        'infected_individuals_at_start': 10,
        'final_time': 365 * 3
    }

    n_seeds = 300

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
    
    
