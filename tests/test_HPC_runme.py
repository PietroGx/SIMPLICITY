#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 11:07:15 2025

@author: pietro
"""
from simplicity.runme import run_experiment
from tests.test_local_runme import test_experiment_output
import simplicity.runners.slurm
import argparse

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {}
    n_seeds         = 10
    return (parameters, n_seeds)

##### <actual test>
def test_run_experiment_HPC(test_number:int):

    print('')
    print('##########################################')
    print('testing HPC runner')
    print('##########################################')
    experiment_name = f'test_HPC_experiment_#{test_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings,             
                       simplicity_runner  = simplicity.runners.slurm,
                       plot_trajectory = True,
                       archive_experiment = False)
    except:
        raise RuntimeError('The code did not pass the running test')
    print('')
    print(f'TEST HPC RUNNER #{test_number} -- SUCCESS.')
    print('##########################################')
    return experiment_name
    
def test_experiment_HPC(test_number:int):
    test_experiment_output(test_run_experiment_HPC(test_number))
    
##### </actual test>

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run test experiment on HPC")
    parser.add_argument('param', type=int, help="Test number")
    args = parser.parse_args()
    # Run the test with the provided parameter
    test_experiment_HPC(args.param)

if __name__ == "__main__":
    main()