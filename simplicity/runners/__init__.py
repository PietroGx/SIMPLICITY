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


Abstract functions for running the repeats of an experiment given its name.
Each repeat is run by a run_seeded_simulation function passed as argument.

A repeat is identified by (group, index), not by a filesystem path. The order is
recorded once in simplicity.jobs' repeats.json and every reader uses it; it used
to come from an unsorted os.walk performed separately by the submitting process
and by each task.

Current implementations:                Target         Description

- simplicity.runners.serial             Single host    Simple pure python for loop.
- simplicity.runners.multiprocessing    Single host    Uses CPython built-in concurrent.futures.ProcessPoolExecutor.
- simplicity.runners.slurm              Cluster        Submits and monitors seeded simulations as a Slurm job array (Slurm must be installed).
"""
import typing


def run_seeded_simulations(self, experiment_name: str,
                           run_seeded_simulation: typing.Callable[[str, str, int], None],
                           groups: typing.Optional[typing.Sequence[str]] = None):
    """Abstract function for running the repeats of an experiment given its name.

    experiment_name: str             reference to the experiment for other simplicity component like simplicity.settings_manager and simplicity.output_manager
    run_seeded_simulation: Callable  function with signature (experiment_name: str, group: str, index: int) -> None
    groups: sequence of str or None  which groups to run; None means all of them.
                                     Dispatch is deliberately not per-group: one
                                     submission may span several, since cal_2
                                     runs all of its (cell, scenario) groups as
                                     a single Slurm array.
    """
    import simplicuty.runners
    help(simplicuty.runners)
    raise NotImplementedError("hint: call an implementation instead.")


def run_seeded_simulation(experiment_name: str, group: str, index: int) -> None:
    """Abstract function for running a single repeat of an experiment.

    experiment_name: str   reference to the experiment for other simplicity component like simplicity.settings_manager and simplicity.output_manager
    group: str             which group of the experiment this repeat belongs to
    index: int             the repeat's position in that group's repeats.json, from 0
    """
    import simplicuty.runners
    help(simplicuty.runners)
    raise NotImplementedError("hint: call an implementation instead.")
