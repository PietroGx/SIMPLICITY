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
Created on Mon Feb 24 17:54:52 2025

@author: pietro
"""
from simplicity.runme import run_experiment
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm

def run_experiment_script(runner:str, 
                          experiment_number:int, 
                          user_set_experiment_settings,
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
                       user_set_experiment_settings,             
                       simplicity_runner  = runner_module,
                       archive_experiment = False)
    except Exception as e:
        print(f'The simulation failed to run: {e}')
        
    print('')
    print(f'{experiment_name} #{experiment_number} -- COMPLETED.')
    print('##########################################')