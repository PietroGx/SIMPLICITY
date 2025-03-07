#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 17:54:52 2025

@author: pietro
"""
from simplicity.runme import run_experiment
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm

def run_experiment_script(runner:str, 
                          experiment_number:int, 
                          fixture_experiment_settings,
                          experiment_name):
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
    print(f'Running {experiment_name}')
    print('##########################################')
    print('')
    try:
        run_experiment(f'{experiment_name}_#{experiment_number}', 
                       fixture_experiment_settings,             
                       simplicity_runner  = runner_module,
                       plot_trajectory = True,
                       archive_experiment = False)
    except Exception as e:
        print(f'The simulation failed to run: {e}')
        
    print('')
    print(f'{experiment_name} #{experiment_number} -- COMPLETED.')
    print('##########################################')