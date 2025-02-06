#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 11:07:15 2025

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
import simplicity.plots_manager as pm
import warnings

## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings():
    parameters      = {'evolutionary_rate': [0.00001, 
                                             0.0001,
                                             0.001,
                                             0.01,
                                             0.1],
                       'final_time':[365,
                                     365,
                                     365, 
                                     365,
                                     365],
                       'max_runtime': [1200,
                                       1200,
                                       1200,
                                       1200,
                                       1200
                                       ]
                           }
    n_seeds = 3
    return (parameters, n_seeds)

def plot_regressions_and_export(experiment_name):
    pm.plot_combined_regressions(experiment_name)
    pm.plot_u_vs_parameter(experiment_name,'evolutionary_rate')
    pm.export_u_regression_plots(experiment_name)

##### <actual test>
def test_param_space_local(runner:str, test_number:int):
    if runner == 'serial':
        runner_module = simplicity.runners.serial
    elif runner == 'multiprocessing':
        runner_module = simplicity.runners.multiprocessing
    else:
        raise ValueError('Runner must be either "serial" or "multiprocessing" ')
    print('')
    print('##########################################')
    print('testing parameter space for e/u regression')
    print('##########################################')
    print('Running first batch, t_final=100')
    print('')
    experiment_name = f'test_param_space_{test_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = True,
                       archive_experiment = False)
    except RuntimeError:
            warnings.warn(f'Experiment {experiment_name} already ran, proceeding to plotting')
    except Exception as e:
        print(f'The simulation failed to run: {e}')
    plot_regressions_and_export(experiment_name)
        
    print('')
    print(f'TEST PARAM SPACE #{test_number} -- COMPLETED.')
    print('##########################################')
 

 
##### </actual test>

if __name__ == "__main__":
    plot_regressions_and_export('test_param_space_3')
    # test_param_space_local('multiprocessing',3)