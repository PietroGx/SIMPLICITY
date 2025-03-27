 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro
   
Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""

from experiment_script_runner import run_experiment_script
import simplicity.settings_manager as sm
import argparse

experiment_name =  'Experiment_name'

def user_set_experiment_settings():
    
    # --------- Specify parameter values manually -----------------------------
    
    # parameters value to get combinations from
    varying_params = {
        'phenotype_model': ['distance from wt', 'immune waning']
    }
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {
        'infected_individuals_at_start': 100,
        'final_time': 365*3
    }
    
    # ---------- OR import them from file -------------------------------------
   
    # # leave empty if you only want to import parameters values from file
    # varying_params = {}
    # # import fixed parameters from user geenerated file. You can either create 
    # # it manually or use the provided script: generate_user_set_parameters_file.py
    # filename = 'standard_values.json'
    # fixed_params = sm.read_user_set_parameters_file(filename)
    
    # -------------------------------------------------------------------------
    n_seeds = 100
    
    return (varying_params,fixed_params,n_seeds)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to generate IH lineages data")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          user_set_experiment_settings,
                          experiment_name)

if __name__ == "__main__":
    main()
    
    
