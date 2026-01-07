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

@author: pietro
   
Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
"""

from experiment_script_runner import run_experiment_script
import simplicity.settings_manager as sm
import argparse

experiment_name =  'test_long_shedders'

def user_set_experiment_settings():
    
    # --------- Specify parameter values manually -----------------------------
    
    # parameters value to get combinations from
    varying_params = {
        'long_shedders_ratio' : [0,0.01]
    }
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {
        'population_size': 1000,
        "IH_virus_emergence_rate": 0,
        'infected_individuals_at_start': 10,
        'final_time': 365*2
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

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to generate IH lineages data")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          user_set_experiment_settings,
                          experiment_name)

if __name__ == "__main__":
    main()
    
    
