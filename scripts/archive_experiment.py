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
Created on Thu Aug 29 21:32:21 2024

@author: pietro
"""

import simplicity.output_manager
import argparse

# simplicity.output_manager.archive_experiment('test_refactoring')

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Archive an experiment.")
    parser.add_argument("experiment_name", type=str, help="The name of the experiment to archive")

    # Parse the arguments
    args = parser.parse_args()

    # Archive experiment
    simplicity.output_manager.archive_experiment(args.experiment_name)
    print(f"Experiment '{args.experiment_name}' has been archived.")