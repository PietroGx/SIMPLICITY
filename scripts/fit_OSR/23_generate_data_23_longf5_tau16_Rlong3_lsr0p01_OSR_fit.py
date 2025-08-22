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
"""
  long_evo_rate_f = 5
  tau_3_long = 16
  R_long = 3
  long_shedders_ratio = 0.01
  IH_virus_emergence_rate = 0.01
"""
from scripts.experiments.experiment_script_runner import run_experiment_script
import argparse
import numpy as np

experiment_name = "23_longf5_tau16_Rlong3_lsr0p01_OSR_fit"

def fixture_experiment_settings():
    # number of values for NSR
    nucleotide_substitution_rate_num_values = 15

    # Generate values spaced logarithmically between 10^-6 and 3e-4
    values = np.logspace(np.log10(1e-6), np.log10(3e-4), 
                         num=nucleotide_substitution_rate_num_values)
    nucleotide_substitution_rate_values = values.tolist()

    varying_params = {
        'nucleotide_substitution_rate': nucleotide_substitution_rate_values
    }

    fixed_params = {
        'infected_individuals_at_start': 10,
        'final_time': 365*3,
        'R': 1.1,
        'long_evo_rate_f': 5,
        'tau_3_long': 16,
        'R_long': 3,
        'long_shedders_ratio': 0.01,
        'IH_virus_emergence_rate': 0.01
    }

    n_seeds = 100

    return (varying_params, fixed_params, n_seeds)

def main():
    parser = argparse.ArgumentParser(description="Run script to generate IH lineages data (OSR fit)")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    run_experiment_script(args.runner, args.experiment_number, fixture_experiment_settings, experiment_name)

if __name__ == "__main__":
    main()
