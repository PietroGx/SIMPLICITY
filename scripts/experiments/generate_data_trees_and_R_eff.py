#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro
   
For each parameter, specify a list of values that you would like to use for the 
simulation. If you want to change more than one parameter at the time, consider 
that you need to enter the same number of values for each parameter, e.g. :
    par 1 = [value1, value2]
    par 2 = [value3, value4]
This will run a simulation with par 1 = value1 and par 2 = value 3, and a simulation
with par 1 = value2 and par 2 = value4. 

Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""
from experiment_script_runner import run_experiment_script
import argparse

experiment_name =  'generate_data_trees_and_R_eff'

def fixture_experiment_settings():
    
    parameters      = {'infected_individuals_at_start': [1],
                       'final_time':[365*3]
                       }
    n_seeds = 100

    return (parameters, n_seeds)

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
    
    
