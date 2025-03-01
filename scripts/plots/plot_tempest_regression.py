 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:08:05 2025

@author: pietro
"""
import simplicity.plots_manager as pm
import simplicity.output_manager as om
import argparse
import pandas as pd 

def plot_regressions_and_export(experiment_name, parameter, min_sim_lenght):
    print('')
    print('')
    print('##################################################################')
    print('################## Plot e vs u relationship ######################')
    print('##################################################################')
    print('')
   
    pm.plot_combined_observed_evolutionary_rate_vs_parameter(experiment_name,parameter, min_sim_lenght)
    print('')
    print('##################################################################')
    print('################## Plot combined regressions #####################')
    print('##################################################################')
    print('')
    csv_with_values = om.get_combined_observed_evolutionary_rate_vs_parameter_df_file_path(experiment_name,parameter)
    
    y_axis_max = max(pd.read_csv(csv_with_values)['observed_evolutionary_rate'])*1.2
    pm.plot_combined_tempest_regressions(experiment_name, parameter, min_sim_lenght, y_axis_max)
    print('')
    print('##################################################################')
    print('######## Exporting plots to Tempest Regression Folder ############')
    print('##################################################################')
    print('')
    print('')
    pm.export_tempest_regression_plots(experiment_name)
    

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to perform and plot u regression")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    parser.add_argument('parameter', type=str, help="parameter to compare with u in regression plots")
    parser.add_argument('min_sim_lenght', type=float, help="only keep data from simulations that lasted at least this number of days")
    args = parser.parse_args()
    # Run the script with the provided parameter
    plot_regressions_and_export(args.experiment_name, args.parameter, args.min_sim_lenght)
    
if __name__ == "__main__":
    main()