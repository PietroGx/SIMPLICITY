# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 12:51:52 2025

@author: pietro
"""

import sys
import simplicity.tuning.evolutionary_rate as er
import simplicity.settings_manager as sm
import simplicity.output_manager as om

def get_NSR_for_model_from_OSR(experiment_name, OSR):
    model_type = 'exp'
    fit_results_params_df = om.read_fit_results_csv(experiment_name, model_type)
    params = fit_results_params_df.to_dict()
    NSR = er.inverse_exp_regressor(OSR, params)
   
    return NSR

def get_NSR_for_model_from_OSR_using_reference_file(OSR):
    fit_results_params_df = sm.read_OSR_NSR_regressor_parameters()
    params = fit_results_params_df.to_dict()
    NSR = er.inverse_exp_regressor(OSR, params)
    return NSR

def get_NSR():
    ''' Prompts user to give observed evolutionary rate for the simulation and returns the correct NSR parameter
    '''
    while True:
        print('')
        print('-----')
        print('Entering setup of substitution rate for SIMPLICITY')
        print('-----')
        print('')
        check_experiment_run = input("Did you run scripts/experiments/00_generate_data_OSR_fit.py? (y/n/q): ").strip().lower()
        check_fit_run = input("Did you run scripts/plots/plot_OSR_fit.py or scripts/fit_observed_evolutionary_rate_regressor.py? (y/n/q): ").strip().lower()
        
        if check_experiment_run == "y" and check_fit_run == "y":
            experiment_name = input("Enter the name of the experiment used for the observed substitution rate (OSR) fit: ")
            print('')
            print('The observed substitution rate (OSR) value should be within the range obtained from the fitting.')
            print(f'If you are unsure, check the fit plot (Data/{experiment_name}/05_Plots/figure4_OSR_exp_fit.png)')

            
            while True:
                OSR = input("Enter the desired observed substitution rate (OSR) value: ")
                try:
                    OSR_value = float(OSR)
                    return get_NSR_for_model_from_OSR(experiment_name, OSR_value)
                except Exception as e:
                    print(e)
                    sys.exit() 
                    
        elif check_experiment_run == "n" or check_fit_run == "n":
            print('')
            print('If you did not run generate_data_OSR_fit.py or fitted the curve, we will use the repo provided values.')
            print('----------------------')
            
            while True:
                OSR = input("Enter the desired observed substitution rate (OSR) value (OSR must be between 0.0001 and 0.001): ")
                try:
                    OSR_value = float(OSR)
                    if 0.0001 <= OSR_value <= 0.001:
                        return get_NSR_for_model_from_OSR_using_reference_file(OSR_value)  
                    else:
                        print("Invalid OSR value. It must be between 0.0001 and 0.001. Please reenter value.")
                except Exception as e:
                    print(e)
                    sys.exit() 

        elif check_experiment_run == "q":
            print("Exiting...")
            sys.exit()

        else:
            print("Invalid input, please select 'y', 'n', or 'q' to quit.")


def get_valid_input(prompt, expected_type, standard_value, min_value=None, max_value=None):
    while True:
        user_input = input(f"{prompt} (Press Enter to keep standard value: {standard_value}) or type 'q' to quit: ").strip()
        if user_input.lower() == 'q':
            print("Process exited.")
            sys.exit()
        if user_input == "":
            print("Keeping standard value.")
            return standard_value
        try:
            if expected_type == int:
                value = int(user_input)
                if min_value is not None and value < min_value:
                    print(f"Value must be at least {min_value}.")
                    continue
                if max_value is not None and value > max_value:
                    print(f"Value must be at most {max_value}.")
                    continue
            elif expected_type == float:
                value = float(user_input)
                if min_value is not None and value < min_value:
                    print(f"Value must be at least {min_value}.")
                    continue
                if max_value is not None and value > max_value:
                    print(f"Value must be at most {max_value}.")
                    continue
            else:
                value = str(user_input)
            return value
        except ValueError:
            print(f"Please enter a valid {expected_type.__name__} value.")


def get_user_input():
    params = {}
    standard_values = sm.read_standard_parameters_values()
    parameter_specs = sm.read_parameter_specs()
    
    no_prompt_params = {"seed","F"}  # Parameters not prompted to the user
    
    print("""
    Welcome! This program will guide you through creating a settings file for running your SIMPLICITY simulations.
    Input the prompted parameters and we will generate a configuration file for you.
    """)
    
    for key, standard_value in standard_values.items():
        param_spec = parameter_specs.get(key, {})
        expected_type = eval(param_spec.get("type", "str"))  # Convert type from string
        min_value = param_spec.get("min")
        max_value = param_spec.get("max")
        
        if key in no_prompt_params:
            params[key] = standard_value  # Automatically set from standard values
        elif key == 'nucleotide_substitution_rate':
            params[key] = get_NSR()
        else:
            print('')
            prompt = f"Please enter {key} value"
            params[key] = get_valid_input(prompt, expected_type, standard_value, min_value, max_value)
    
    return params

def main():
    user_set_parameters = get_user_input()
    filename = input("Enter the filename to save the settings (default: user_parameters.json): ").strip()
    if not filename:
        filename = "user_parameters"
    filename = filename + '.json'
    sm.write_user_set_parameters_file(user_set_parameters, filename)
    
if __name__ == "__main__":
    main()
