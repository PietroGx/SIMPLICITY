#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 
    
    STANDARD_VALUES = {
        "population_size": 1000,
        "tau_3": 7.5,
        "infected_individuals_at_start": 10,
        "R": 1.3,
        "diagnosis_rate": 0.1,             # in percentage, will be converted to kd in model 
        "IH_virus_emergence_rate": 0,      # k_v in theoretical model equations
        "evolutionary_rate": 0.0017,       # e in theoretical model equations
        "final_time": 365,
        "max_runtime": 259200, 
        "phenotype_model": 'immune waning',  # or 'distance from wt'
        "sequencing_rate": 0.05,
        "seed": None,
        "F": 1.25
    }
    
y, you can specify them in the parameters dictionary below. 
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
import numpy as np
experiment_name =  'generate_data_OER_fit'

def fixture_experiment_settings():

    # number of values for evolutionary rate
    evolutionary_rate_num_values = 15
    
    # Generate values spaced logarithmically between 10^-5 and 10^-2
    values = np.logspace(np.log10(1e-7), np.log10(0.2), 
                         num=evolutionary_rate_num_values)
    evolutionary_rate_values = values.tolist()
    mapped_sim_lenght = np.logspace(np.log10(3650), np.log10(365), num=evolutionary_rate_num_values).tolist()

    
    parameters      = {'evolutionary_rate': evolutionary_rate_values,
                       'infected_individuals_at_start': [50]*evolutionary_rate_num_values,
                       'final_time': mapped_sim_lenght
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
    
    
