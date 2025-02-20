#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
@author: pietro
"""
from simplicity.runme import run_experiment
# import simplicity.runners.slurm
import simplicity.runners.serial
import argparse
from scripts.slurm_diagnostics.slurm_error_summary import print_slurm_error_summary
from scripts.check_completed_simulations import count_completed_simulations
from memory_profiler import profile

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {"evolutionary_rate": [1,
                                             10,
                                             100]}
    n_seeds         = 3
    return (parameters, n_seeds)

##### <actual test>
def test_HPC_OOM_run(test_number:int):

    print('')
    print('##########################################')
    print('testing HPC for OOM errors')
    print('##########################################')
    experiment_name = f'test_HPC_OOM_run_#{test_number}'
    try:
        profile(run_experiment)(experiment_name, 
                           fixture_experiment_settings,             
                           simplicity_runner  = simplicity.runners.serial,
                           plot_trajectory = True,
                           archive_experiment = False)
    except: pass
    return experiment_name

##### </actual test>

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run test experiment on HPC")
    parser.add_argument('test_number', type=int, help="Test number")
    args = parser.parse_args()
    # Run the test 
    experiment_name = test_HPC_OOM_run(args.test_number)
    # count completed simulations
    count_completed_simulations(experiment_name)
    # print summary of slurm errors 
    print_slurm_error_summary(experiment_name)
    

if __name__ == "__main__":
    main()