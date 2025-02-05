#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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
import simplicity.plots_manager as pm
import warnings
import argparse

n_seeds_script = 500
## fixture  experiment settings (sm.write_settings arguments)
def fixture_experiment_settings_t100():
    parameters      = {'evolutionary_rate': [0.00001, 
                                             0.0001,
                                             0.001,
                                             0.01,
                                             0.1],
                       'final_time':[100,
                                     100,
                                     100, 
                                     100,
                                     100
                           ]}
    n_seeds = n_seeds_script
    return (parameters, n_seeds)

def fixture_experiment_settings_t500():
    parameters      = {'evolutionary_rate': [0.00001, 
                                             0.0001,
                                             0.001,
                                             0.01,
                                             0.1],
                       'final_time':[500,
                                     500,
                                     500,
                                     500,
                                     500
                           ]}
    n_seeds = n_seeds_script
    return (parameters, n_seeds)

def fixture_experiment_settings_t1000():
    parameters      = {'evolutionary_rate': [0.00001, 
                                             0.0001,
                                             0.001,
                                             0.01,
                                             0.1
                           ]}
    n_seeds = n_seeds_script
    return (parameters, n_seeds)

def plot_regressions_and_export(experiment_name):
    pm.plot_combined_regressions(experiment_name)
    pm.plot_u_vs_parameter(experiment_name,'evolutionary_rate')
    pm.export_u_regression_plots(experiment_name)

def explore_u_e_space(runner:str, experiment_number:int):
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
    print('testing parameter space for e/u regression')
    print('##########################################')
    print('Running first batch, t_final=100')
    print('')
    experiment_name = f'explore_param_space_tfinal=100_#{experiment_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings_t100,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = False,
                       archive_experiment = False)
    except RuntimeError:
            warnings.warn(f'Experiment {experiment_name} already ran, proceeding to plotting')
    except Exception as e:
        print(f'The simulation failed to run: {e}')
    plot_regressions_and_export(experiment_name)
        
    print('##########################################')
    print('Running second batch, t_final=500')
    print('')
    experiment_name = f'explore_param_space_tfinal=500_#{experiment_number}'

    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings_t500,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = False,
                       archive_experiment = False)
    except RuntimeError:
            warnings.warn(f'Experiment {experiment_name} already ran, proceeding to plotting')
    except Exception as e:
        print(f'The simulation failed to run: {e}')
    plot_regressions_and_export(experiment_name)
        
    print('##########################################')
    print('Running third batch, t_final=1000')
    print('')
    experiment_name = f'explore_param_space_tfinal=1000_#{experiment_number}'
    try:
        run_experiment(experiment_name, 
                       fixture_experiment_settings_t1000,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = False,
                       archive_experiment = False)
    except RuntimeError:
            warnings.warn(f'Experiment {experiment_name} already ran, proceeding to plotting')
    except Exception as e:
        print(f'The simulation failed to run: {e}')
    plot_regressions_and_export(experiment_name)
        
    print('')
    print(f'EXPLORATION OF E/U PARAM SPACE #{experiment_number} -- COMPLETED.')
    print('##########################################')
 
def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to explore u/e parameter space")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the test with the provided parameter
    explore_u_e_space(args.runner,args.experiment_number)

if __name__ == "__main__":
    main()
    # explore_u_e_space('serial',1)
