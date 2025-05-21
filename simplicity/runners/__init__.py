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
@author: jbescudie


Abstract functions for running the seeded simulations of an experiment given its name.
Each seeded simulation is run by a run_seeded_simulation function passed as argument.

Current implementations:                Target         Description

- simplicity.runners.serial             Single host    Simple pure python for loop.
- simplicity.runners.multiprocessing    Single host    Uses CPython built-in concurrent.futures.ProcessPoolExecutor.
- simplicity.runners.slurm              Cluster        Submits and monitors seeded simulations as a Slurm job array (Slurm must be installed).
"""
import typing


def run_seeded_simulations(self, experiment_name: str, run_seeded_simulation: typing.Callable[[str,str],None]):
    """Abstract function for running the seeded simulations of an experiment given its name.

    experiment_name: str             reference to the experiment for other simplicity component like simplicity.settings_manager and simplicity.output_manager
    run_seeded_simulation: Callable  function with signature (seeded_simulation_parameters_path: str, experiment_name: str) -> None
    """
    import simplicuty.runners
    help(simplicuty.runners)
    raise NotImplementedError("hint: call an implementation instead.")


def run_seeded_simulation(seeded_simulation_parameters_path: str, experiment_name: str) -> None:
    """Abstract function for running a single seeded simulation of an experiment.

    seeded_simulation_parameters_path: str    path to the file containing the simulation's parameters incl. seed
    experiment_name: str                      reference to the experiment for other simplicity component like simplicity.settings_manager and simplicity.output_manager
    """
    import simplicuty.runners
    help(simplicuty.runners)
    raise NotImplementedError("hint: call an implementation instead.")
