#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STANDARD_VALUES for SIMPLICITY simulation: 

If you want to change any parameter, you can specify them in the parameters dictionary below. 
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
import simplicity.runners.serial
import argparse
from scripts.check_completed_simulations import count_completed_simulations


## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {"nucleotide_substitution_rate": [0.0001,
                                             0.001,
                                             0.01,
                                             0.1,
                                             1]}
    n_seeds         = 1
    return (parameters, n_seeds)

def diagnosis_OOM_errors(test_number:int):

    print('')
    print('##########################################')
    print('testing HPC for OOM errors')
    print('##########################################')
    experiment_name = f'test_HPC_OOM_run_#{test_number}'
    try:
        run_experiment(experiment_name, 
                           fixture_experiment_settings,             
                           simplicity_runner  = simplicity.runners.serial,
                           plot_trajectory = True,
                           archive_experiment = False)
    except: pass
    return experiment_name

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run test to check OOM relation to NSR parameter")
    parser.add_argument('test_number', type=int, help="Test number")
    args = parser.parse_args()
    # Run 
    experiment_name = diagnosis_OOM_errors(args.test_number) 
    # count completed simulations
    count_completed_simulations(experiment_name)
    
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    