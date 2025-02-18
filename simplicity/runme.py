#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: pietro

STANDARD_VALUES for SIMPLICITY simulation: 
    
    "population_size": 1000
    "infected_individuals_at_start": 100
    "R": 1.5
    "k_d": 0.0055
    "k_v": 0.0085
    "e": 0.0017 (evolutionary rate)
    "final_time": 365*3 
    "max_runtime": 300
    "phenotype_model": 'immune waning' or 'distance from wt'
    "sequencing_rate": 0.05
    "seed": None
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
import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om
import simplicity.runners.serial 
from   simplicity.runners.unit_run import run_seeded_simulation
import types
experiment_name = 'Test_runme'

## experiment settings (sm.write_settings arguments)

def user_set_experiment_settings():
    parameters      = {}
    n_seeds         = 2
    return (parameters, n_seeds)

def run_experiment(experiment_name: str, 
                   experiment_settings: types.FunctionType,
                   simplicity_runner: types.ModuleType,
                   plot_trajectory:bool,
                   archive_experiment=False):
    print('')
    print('##########################################')
    print('')
    # setup experiment files directories
    dm.create_directories(experiment_name)
    # set parameters 
    parameters, n_seeds = experiment_settings()
    # Write experiment settings file
    sm.write_experiment_settings(experiment_name, parameters, n_seeds)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(experiment_name)
    # write seeded simulation parameters files
    sm.write_seeded_simulation_parameters(experiment_name)
    
    # let one of simplicity.runners run each seeded simulation
    simplicity_runner.run_seeded_simulations(experiment_name, 
                                             run_seeded_simulation,
                                             plot_trajectory)
    
    if archive_experiment: 
        om.archive_experiment(experiment_name)
    
    print('')
    print('##########################################')
    print(f'EXPERIMENT {experiment_name} EXECUTED SUCCESSFULLY.')
    print('##########################################')

if __name__ == "__main__":
    
    run_experiment(experiment_name,       
                   user_set_experiment_settings,
                   simplicity_runner  = simplicity.runners.serial,
                   plot_trajectory = True,
                   archive_experiment = False)