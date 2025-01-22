#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 21:32:21 2024

@author: pietro
"""

import simplicity.output_manager
import argparse

# simplicity.output_manager.archive_experiment('test_refactoring')

if __name__ == "__main__":
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description="Archive an experiment.")
    parser.add_argument("experiment_name", type=str, help="The name of the experiment to archive")

    # Parse the arguments
    args = parser.parse_args()

    # Archive experiment
    simplicity.output_manager.archive_experiment(args.experiment_name)
    print(f"Experiment '{args.experiment_name}' has been archived.")