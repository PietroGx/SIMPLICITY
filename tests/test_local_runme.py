#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 11:07:15 2025

@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 

 
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
from simplicity.dir_manager import get_experiment_output_dir
import simplicity.settings_manager as sm
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import os

def fixture_experiment_settings():
    # --------- Specify parameter values manually -----------------------------
    
    # parameters value to get combinations from
    # varying_params = {
    #     'phenotype_model': ['distance from wt', 'immune waning']
        
    # }
    # # parameters to keep fixed (but different from standard_value) across combinations
    # fixed_params = {
    #     'infected_individuals_at_start': 10,
    #     'final_time': 365
    # }
    varying_params = {
        'R': [0.8,1,1.3,1.5,2],
        'diagnosis_rate': [0,0.05,0.1],
        'phenotype_model': ['linear', 'immune_waning']
    }

    fixed_params = {
        'infected_individuals_at_start': 10,
        'final_time': 365
    }
    
    # ---------- OR import them from file -------------------------------------
    
    # # leave empty if you only want to import parameters values from file
    # varying_params = {}
    # # import fixed parameters from user geenerated file. You can either create 
    # # it manually or use the provided script: generate_user_set_parameters_file.py
    # filename = 'standard_values.json'
    # fixed_params = sm.read_user_set_parameters_file(filename)
   
    # -------------------------------------------------------------------------
    n_seeds = 10
    
    return (varying_params,fixed_params,n_seeds)

def check_output_file(directory, filename): 
    file_path = os.path.join(directory, filename)
    if os.path.isfile(file_path):
        pass
    else:
        raise ValueError(f"'{filename}' does not exist in the directory.")
        
##### <actual test>
def test_run_experiment_local(runner:str, test_number:int):
    if runner == 'serial':
        runner_module = simplicity.runners.serial
    elif runner == 'multiprocessing':
        runner_module = simplicity.runners.multiprocessing
    else:
        raise ValueError('Runner must be either "serial" or "multiprocessing" ')
    print('')
    print('##########################################')
    print(f'testing local {runner} runner')
    print('##########################################')
    experiment_name = f'test_local_experiment_{runner}_#{test_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = True,
                       archive_experiment = False)
    except Exception as e:
        raise RuntimeError(f'The code did not pass the running test. Error: {e}')
    print('')
    print(f'TEST LOCAL {runner} RUNNER {test_number} -- SUCCESS.')
    print('##########################################')
    return experiment_name

def test_experiment_output(experiment_name):
    print('')
    print('##########################################')
    print(f'testing {experiment_name} output')
    print('##########################################')
    # get experiment output directory
    directory = get_experiment_output_dir(experiment_name)
    # loop over each simulation output folder
    for folder in os.listdir(directory):
        simulation_directory = os.path.join(directory,folder)
        # loop over each seeded simulation output folder
        for subfolder in os.listdir(simulation_directory):
            seed_directory = os.path.join(simulation_directory,subfolder)
            check_output_file(seed_directory, 'final_time.csv')
            check_output_file(seed_directory, 'fitness_trajectory.csv')
            check_output_file(seed_directory, 'individuals_data.csv')
            check_output_file(seed_directory, 'lineage_frequency.csv')
            check_output_file(seed_directory, 'phylogenetic_data.csv')
            check_output_file(seed_directory, 'sequencing_data_regression.csv')
            check_output_file(seed_directory, 'sequencing_data.fasta')
            check_output_file(seed_directory, 'simulation_trajectory.csv')
    print('')
    print(f'TEST {experiment_name} OUTPUT -- SUCCESS.')
    print('##########################################')
    
def test_experiment_local(runner:str, test_number:int):
    test_experiment_output(test_run_experiment_local(runner, test_number))
    
##### </actual test>

if __name__ == "__main__":
    # test_experiment_local('serial',4)
    test_experiment_local('multiprocessing',3)