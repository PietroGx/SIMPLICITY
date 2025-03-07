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
    min_seq_number = 10
    min_sim_lenght = 0
    
    # build the dataframe needed for the plot
    om.write_combined_OSR_vs_parameter_csv(args.experiment_name, 
                                            parameter, 
                                            min_seq_number,
                                            min_sim_lenght)
    om.write_OSR_vs_parameter_csv(args.experiment_name, 
                                    parameter, 
                                    min_seq_number,
                                    min_sim_lenght)
    # Run the script with the provided parameter
    pm.plot_combined_OSR_vs_parameter(args.experiment_name,
                                      parameter, 
                                      min_seq_number, 
                                      min_sim_lenght)
    
if __name__ == "__main__":
    main()
    