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
from simplicity.runme import run_experiment, run_experiment_phase
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm


def _runner_module(runner: str):
    if runner == 'serial':
        return simplicity.runners.serial
    if runner == 'multiprocessing':
        return simplicity.runners.multiprocessing
    if runner == 'slurm':
        return simplicity.runners.slurm
    raise ValueError('Runner must be either "serial" or "multiprocessing" '
                     'or "slurm"')


def run_experiment_phase_script(runner: str,
                                experiment_number: int,
                                user_set_experiment_settings,
                                experiment_name,
                                append: bool = True):
    """One phase of ONE experiment (stage 8). The stages call this instead of
    run_experiment_script: every phase goes into experiment_name_#N and only
    its own groups are run."""
    runner_module = _runner_module(runner)
    numbered = f'{experiment_name}_#{experiment_number}'
    print('')
    print('##########################################')
    print(f'Running a phase of {numbered}')
    print('##########################################')
    print('')
    try:
        groups = run_experiment_phase(numbered,
                                      user_set_experiment_settings,
                                      simplicity_runner=runner_module,
                                      append=append)
    except Exception as e:
        # same reasoning as run_experiment_script below: a later stage would
        # otherwise fit against data that was never produced
        print(f'The simulation failed to run: {e}')
        raise

    print('')
    print(f'{numbered} groups {groups} -- COMPLETED.')
    print('##########################################')
    return groups


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
        # Re-raise instead of swallowing: a caller several stages downstream
        # (e.g. a calibration fit, or the next pipeline stage) would otherwise
        # proceed against data that was never actually produced -- silently
        # printing "COMPLETED" regardless of whether anything ran -- and fail
        # later with a confusing, unrelated-looking error instead of this one.
        print(f'The simulation failed to run: {e}')
        raise

    print('')
    print(f'{experiment_name} #{experiment_number} -- COMPLETED.')
    print('##########################################')