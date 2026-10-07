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
from tqdm import tqdm

def run_seeded_simulations(experiment_name, run_seeded_simulation, groups=None):
    import simplicity.jobs as jobs
    repeats = jobs.all_repeats(experiment_name, groups)
    for group, record in tqdm(repeats, desc="Running repeats", unit="sim",
                              position=0):
        run_seeded_simulation(experiment_name, group, record['index'])
