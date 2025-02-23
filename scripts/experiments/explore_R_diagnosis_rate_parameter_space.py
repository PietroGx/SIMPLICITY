#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 15:23:13 2025

@author: pietro
"""

"""

@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 

    "population_size": 1000,
    "infected_individuals_at_start": 100,
    "R": 1.5,
    "diagnosis_rate": 0.0055,
    "IH_virus_emergence_rate": 0.0085,
    "evolutionary_rate": 0.0017,
    "final_time": 365*3 ,
    "max_runtime": 100000000, 
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
from simplicity.runme import run_experiment
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm
import argparse
import itertools
import numpy as np 

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    # Given lists
    # R_values       = [0.9,1.1,1.3,1.5,2]
    diagnosis_rates = [round(num, 2) for num in np.arange(0.01,0.11,0.01)]
    
    # # Repeat the elements in `a` according to the length of `b`
    # R_for_params = list(itertools.chain.from_iterable(
    #                            [[x] * len(diagnosis_rates) for x in R_values]))
    
    # # Repeat the elements in `b` for each element in `a`
    # k_d_for_params = list(itertools.chain.from_iterable([diagnosis_rates] * len(R_values)))
    
    # parameters      = {'R': R_for_params,
    #                    'diagnosis_rate':k_d_for_params
    #                    }
    
    parameters = {'diagnosis_rate': diagnosis_rates}
    n_seeds = 10

    return (parameters, n_seeds)

def explore_R_k_d_space(runner:str, experiment_number:int):
    if runner == 'serial':
        runner_module = simplicity.runners.serial
    elif runner == 'multiprocessing':
        runner_module = simplicity.runners.multiprocessing
    elif runner == 'slurm':
        runner_module = simplicity.runners.slurm
    else:
        raise ValueError('Runner must be either "serial" or "multiprocessing" or "slurm"')
    print('')
    print('##########################################')
    print('testing parameter space for diagnosis rate tuning')
    print('##########################################')
    print('')
    experiment_name = f'test_d_rates_#{experiment_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = False,
                       archive_experiment = False)
    # except RuntimeError:
    #         warnings.warn(f'Experiment {experiment_name} already ran, proceeding to plotting')
    except Exception as e:
        print(f'The simulation failed to run: {e}')
        
    print('')
    print(f'EXPLORATION OF R/k_d PARAM SPACE #{experiment_number} -- COMPLETED.')
    print('##########################################')
 
def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to explore R/k_d parameter space")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script with the provided parameter
    explore_R_k_d_space(args.runner,args.experiment_number)

if __name__ == "__main__":
    main()