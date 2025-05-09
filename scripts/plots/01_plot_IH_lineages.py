#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:08:05 2025

@author: pietro
"""
import simplicity.plots_manager as pm
import argparse

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Plot experiment_name IH_lineage_distribution")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    # Run the script with the provided parameter
    pm.plot_IH_lineage_distribution(args.experiment_name)
    
if __name__ == "__main__":
    main()
    