#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:08:05 2025

@author: pietro
"""
import simplicity.plots_manager as pm
import argparse

def plot_regressions_and_export(experiment_name):
    pm.plot_combined_regressions(experiment_name)
    pm.plot_u_vs_parameter(experiment_name,'evolutionary_rate')
    pm.export_u_regression_plots(experiment_name)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to perform and plot u regression")
    args = parser.parse_args()
    # Run the script with the provided parameter
    plot_regressions_and_export(args.experiment_name)
    
if __name__ == "__main__":
    main()