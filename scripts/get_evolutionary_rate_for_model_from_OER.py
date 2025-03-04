#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 15:34:12 2025

@author: pietro
"""
import simplicity.output_manager as om
import argparse

def get_evolutionary_rate_for_model_from_OER(experiment_name):
    model_type = 'log'
    fit_results_params_df = om.read_fit_results_csv(experiment_name, model_type)
    print(fit_results_params_df)
    # inverse_log_regressor(OER, params)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    get_evolutionary_rate_for_model_from_OER(args.experiment_name)
    
    
if __name__ == "__main__":
    main()