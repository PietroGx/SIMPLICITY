#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 14:46:04 2024

@author: pietro
@author: jbescudie
"""
import simplicity.config as config
import simplicity.settings_manager as sm
import simplicity.tuning.evolutionary_rate as er
import os 
import json
import numpy as np
import math

def evaluate_u(experiment_name, target_u: float)->float: 
    """ 
    """
    
    # Do regression and compute error
    df, _ = er.create_joint_df(experiment_name)
    u, _  = er.perform_regression(df)

    error = target_u - u
    return error


def read_next_step(experiment_name):
    # look for maximal step in files with name pattern 'settings_step_{step}.json'
    # return max(-1, step) + 1
    step = -1
    for filename in os.listdir(config.get_data_dir()):
        if filename.startswith(f'settings_step_') and filename.endswith('.json'):
            _step = int(filename[len('settings_step_'):len('.json')]) 
            assert _step == read_step(experiment_name, _step), f"filename '{filename}' content does not match step {_step}"
            step = max(step, 0)
    return step + 1

def write_paremeter_search_settings(experiment_name, threshold: float):
    settings = {
        "threshold": threshold,
    }
    with open(os.path.join(config.get_data_dir(),'parameter_search_settings.json'), 'w') as json_file:
        json.dump(settings, json_file)

def read_parameter_search_settings(experiment_name):
    with open(os.path.join(config.get_data_dir(),'parameter_search_settings.json'), 'r') as json_file:
        return json.load(json_file)

def get_experiment_name_step(experiment_name, step: int):
    return f"{experiment_name}_{step}"
    
def read_settings_step(experiment_name, step: int):
    exp_name = get_experiment_name_step(experiment_name, step)
    return sm.read_experiment_settings(exp_name)
            
def write_settings_step(experiment_name, step: int, settings_step):
    exp_name = get_experiment_name_step(experiment_name, step)
    # setup experiment files directories
    config.create_directories(exp_name)
    # Write experiment settings file
    sm.write_experiment_settings(exp_name)
    # write simulation parameters files
    sm.read_settings_and_write_simulation_parameters(exp_name)
    # write seeded simulation parameters files
    sm.write_seeded_simulation_parameters(exp_name)

def read_threshold(experiment_name, step: int):
    return read_parameter_search_settings(experiment_name)["threshold"]

def read_step_error(experiment_name, step: int):
    with open(os.path.join(config.get_data_dir(),f'error_step_{step}.json'), 'r') as json_file:
        return json.load(json_file)["error"]
    
def write_step_error(experiment_name, step: int, error: float):
    with open(os.path.join(config.get_data_dir(),f'error_step_{step}.json'), 'w') as json_file:        
        json.dump({"error":error})

            

def compute_error_step(experiment_name, step: int):
    target_value = read_settings_step(experiment_name, step)["target_value"]
    error        = evaluate_u(experiment_name, target_value)
    write_step_error(experiment_name, step, error)
    
def update_settings_step(experiment_name, step: int):
    # read previous settings_step.json file
    # and read previous error
    previous_settings = read_settings_step(experiment_name, step)
    previous_error    = read_step_error   (experiment_name, step) if step > 0 else math.inf
              
    # update step
    step += 1
    
    # Calculate step
    learning_rate    = previous_settings["learning_rate"]
    parameter_update = learning_rate * np.sign(error)
    
    # Update parameter
    parameter_name      = list(parameters.keys())[0]
    previous_parameters = previous_settings['parameters']
    updated_parameters  = {**previous_parameters}              # copy parameters
    updated_parameters[parameter_name][0] += parameter_update  # update parameter
    
    # Adjust learning rate based on the change in error
    if abs(error) > abs(previous_error):
        learning_rate /= 2  # If the error increased, reduce the learning rate
    
    # update settings
    updated_settings = {
        "parameters"   : updated_parameters,
        "n_seeds"      : previous_settings['n_seeds'],
        "target_value" : previous_settings['target_value'],
        #"threshold"    : previous_settings['threshold'],
        "learning_rate": learning_rate,
        "error": error,
        "step": step,
        #"experiment_number": experiment_number
    }
    
    # write next settings_step.json file
    write_settings_step(experiment_name, step, updated_settings)
        
               
def converged(experiment_name, step):
    # initial case
    if step < 0:
        return False

    # read threshold, error
    threshold = read_threshold (experiment_name, step)
    error     = read_step_error(experiment_name, step)

    # decide
    converged = abs(error) <= threshold
    return converged


def fit(experiment_name, simplicity_runner, run_seeded_simulation):
    ## <search>
    step  = None
    while step is None or not converged(experiment_name, step):
        step  = read_next_step(experiment_name)
        ## </search>

        # write settings for this step
        write_settings_step(experiment_name, step, )

        # let one of simplicity.runners run each seeded simulation
        sr = simplicity_runner
        sr.run_seeded_simulations(exp_name, run_seeded_simulation)
    
        ## <search>
        # compute error step (parameter search)
        compute_error_step  (experiment_name, step)
        # update step        (parameter search)
        update_settings_step(experiment_name, step)
        print('')
        print('##########################################')
        error = read_step_error(experiment_name, step)
        print(f"#step={step}, error={error}, updating parameter")
        print('')
        ## </search>

        # archive experiment
        om.archive_experiment(exp_name)
            
