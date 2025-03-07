#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: pietro

"""
from experiment_script_runner import run_experiment_script
import argparse

experiment_name = 'generate_data_tau_3_impact_on_OSR_highR'

def fixture_experiment_settings():
    
    parameters      = {'tau_3':[3.25, 
                                7.5,
                                15,
                                30,
                                60,
                                120
                                ],
                       'final_time':[365*3]*6,
                       'R': [2]*6
                       }
    n_seeds = 300

    return (parameters, n_seeds)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to explore tau_3 impact on u")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          fixture_experiment_settings,
                          experiment_name)
    
if __name__ == "__main__":
    main()
    
