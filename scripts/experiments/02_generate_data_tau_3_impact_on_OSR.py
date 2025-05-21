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

"""
from experiment_script_runner import run_experiment_script
import argparse

experiment_name = 'generate_data_tau_3_impact_on_OSR'

def fixture_experiment_settings():
    
    # parameters value to get combinations from
    varying_params = {
        'tau_3':[3.25, 7.5, 15, 30, 60, 120]
    }
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {'final_time': 365*3,
                    'IH_virus_emergence_rate': 0.1}
    
    n_seeds = 300
    
    return (varying_params,fixed_params,n_seeds)

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
    
