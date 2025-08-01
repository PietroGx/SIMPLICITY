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

"""
import simplicity.settings_manager as sm
from tqdm import tqdm

def run_seeded_simulations(experiment_name, run_seeded_simulation):
    seeded_simulation_parameters_paths = sm.get_seeded_simulation_parameters_paths(experiment_name)
    # run a Simplicity simulation for each seeded parameters
    # for seeded_simulation_parameters_path in seeded_simulation_parameters_paths:
    #     run_seeded_simulation(seeded_simulation_parameters_path, experiment_name, plot_trajectory)

    for path in tqdm(seeded_simulation_parameters_paths, desc="Running seeded simulations", unit="sim", position=0):
        run_seeded_simulation(path, experiment_name)
