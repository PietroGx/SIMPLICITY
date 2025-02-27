 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:08:05 2025

@author: pietro
"""
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import argparse
import pandas as pd 
import os 

def plot_regressions_and_export(experiment_name, parameter, min_sim_lenght):
    print('')
    print('')
    print('##################################################################')
    print('################## Plot e vs u relationship ######################')
    print('##################################################################')
    print('')
   
    pm.plot_u_vs_parameter(experiment_name,parameter, min_sim_lenght)
    print('')
    print('##################################################################')
    print('################## Plot combined regressions #####################')
    print('##################################################################')
    print('')
    csv_with_u_values = os.path.join(dm.get_experiment_output_dir(experiment_name),
                                     f'{experiment_name}_u_vs_{parameter}_values.csv')
    
    y_axis_max = max(pd.read_csv(csv_with_u_values)['u'])*1.2
    pm.plot_combined_regressions(experiment_name, parameter, min_sim_lenght, y_axis_max)
    print('')
    print('##################################################################')
    print('######## Exporting plots to Tempest Regression Folder ############')
    print('##################################################################')
    print('')
    print('')
    pm.export_u_regression_plots(experiment_name)

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