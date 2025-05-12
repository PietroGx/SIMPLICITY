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
    # # plot combined OSR vs parameter 
    # pm.plot_combined_OSR_vs_parameter(experiment_name, 
    #                                     parameter,  
    #                                     min_seq_number,
    #                                     min_sim_lenght)
    print('')
    print('##################################################################')
    print('################## Plot combined regressions #####################')
    print('##################################################################')
    print('')
    combined_OSR_vs_parameter_df = om.read_combined_OSR_vs_parameter_csv(experiment_name,
                                        parameter,
                                        min_seq_number,
                                        min_sim_lenght)
    
    y_axis_max = max(combined_OSR_vs_parameter_df['observed_substitution_rate'])*1.2
    pm.plot_combined_tempest_regressions(experiment_name, parameter, min_sim_lenght, y_axis_max)
    

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
    
    pm.plot_figure_tempest_regression(experiment_name)
    print('Plotting completed.')
    
if __name__ == "__main__":
    main()