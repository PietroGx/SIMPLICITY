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

def plot_regressions(experiment_name, parameter, min_seq_number, min_sim_lenght):
    
    # write csv needed for plotting 
    om.write_combined_OSR_vs_parameter_csv(experiment_name, 
                                parameter, 
                                min_seq_number,
                                min_sim_lenght)
    
    om.write_OSR_vs_parameter_csv(experiment_name, 
                                parameter, 
                                min_seq_number,
                                min_sim_lenght)
    print('##################################################################')
    print('################## Plot combined regressions #####################')
    print('##################################################################')
    print('')
    combined_OSR_vs_parameter_df = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    y_axis_max = max(combined_OSR_vs_parameter_df['observed_substitution_rate'])*1.2
    pm.plot_combined_tempest_regressions(experiment_name, parameter, min_seq_number, 
                                         min_sim_lenght, y_axis_max)
    

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to perform and plot u regression")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    # parser.add_argument('parameter', type=str, help="parameter to compare with u in regression plots")
    # parser.add_argument('min_seq_number', type=float, help="only use datasets with at least this number of sequences")
    # parser.add_argument('min_sim_lenght', type=float, help="only keep data from simulations that lasted at least this number of days")
    args = parser.parse_args()
    
    parameter = 'nucleotide_substitution_rate'
    min_seq_number = 30
    min_sim_lenght = 100
    
    # Run the script with the provided parameter
    plot_regressions(args.experiment_name, 
                                parameter, 
                                min_seq_number,
                                min_sim_lenght)
    
    pm.plot_figure_tempest_regression(args.experiment_name)
    print('Plotting completed.')
    
if __name__ == "__main__":
    main()