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
Created on Fri Feb 14 15:08:05 2025

@author: pietro
"""
import simplicity.plots_manager as pm
import simplicity.output_manager as om
import argparse

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    parameter = 'tau_3'
    min_seq_number = 5
    min_sim_lenght = 100
    
    # build the dataframe needed for the plot
    om.write_combined_OSR_vs_parameter_csv(args.experiment_name, 
                                            parameter, 
                                            min_seq_number,
                                            min_sim_lenght)
    om.write_OSR_vs_parameter_csv(args.experiment_name, 
                                    parameter, 
                                    min_seq_number,
                                    min_sim_lenght)

    # pm.plot_combined_OSR_vs_parameter(args.experiment_name,
    #                                   parameter, 
    #                                   min_seq_number, 
    #                                   min_sim_lenght)
    # # plot ih variability
    # pm.plot_IH_lineage_distribution_grouped_by_simulation(args.experiment_name)
    
    pm.plot_OSR_and_IH_lineages_by_parameter(args.experiment_name, 
                                               parameter, 
                                               min_seq_number, 
                                               min_sim_lenght)
    
if __name__ == "__main__":
    main()
    