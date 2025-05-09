#!/uparser.add_argument('experiment_name', type=str, help="experiment name")sr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 15:34:12 2025

@author: pietro
"""
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import argparse

def get_NSR_for_model_from_OSR(experiment_name, OSR):
    model_type = 'exp'
    fit_results_params_df = om.read_fit_results_csv(experiment_name, model_type)
    params = fit_results_params_df.to_dict()

    NSR = er.inverse_log_regressor(OSR, params)
    print(f'nucleotide substitution rate for simulations: {NSR:.5f}')
    return NSR

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    parser.add_argument('OSR', type=float, help="desired observed substitution rate for simulation")
    args = parser.parse_args()
    
    get_NSR_for_model_from_OSR(args.experiment_name, args.OSR)
    
    
if __name__ == "__main__":
    main()