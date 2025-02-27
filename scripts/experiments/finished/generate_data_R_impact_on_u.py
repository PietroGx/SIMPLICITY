#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 

    "population_size": 1000,
    "tau_3": 7.5,
    "infected_individuals_at_start": 10,
    "R": 1.5,
    "diagnosis_rate": 0.1,             # in percentage, will be converted to in model 
    "IH_virus_emergence_rate": 0       # k_v in theoretical model equations
    "evolutionary_rate": 0.0017,       # e in theoretical model equations
    "final_time": 365,
    "max_runtime": 259200, 
    "phenotype_model": 'immune waning',  # or 'distance from wt'
    "sequencing_rate": 0.05,
    "seed": None,
    "F": 1.25
    
If you want to change any, you can specify them in the parameters dictionary below. 
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

R = 1
experiment_name = f'generate_data_R{R}_impact_on_u'

def fixture_experiment_settings():

    R_values       = [R]*6
    evolutionary_rates = [0.00001, 
                        0.0001, 
                        0.001, 
                        0.01, 
                        0.1, 
                        1]
    
    parameters      = {'R': R_values,
                       'evolutionary_rate': evolutionary_rates
                       }
    n_seeds = 100

    return (parameters, n_seeds)
 
def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to explore R impact on u")
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
    
