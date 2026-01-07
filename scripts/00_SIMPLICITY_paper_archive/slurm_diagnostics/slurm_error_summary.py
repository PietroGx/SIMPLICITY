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
Created on Tue Feb 18 20:16:31 2025

@author: pietro
"""

import os
import re
from collections import Counter
import simplicity.dir_manager as dm
import argparse

def parse_slurm_error_files(experiment_name):
    
    folder_path = dm.get_slurm_logs_dir(experiment_name)
    
    # Dictionary to store job_id and corresponding error type
    job_errors = {}
    error_types = Counter()

    # Iterate over all .err files in the specified folder
    for filename in os.listdir(folder_path):
        if filename.endswith(".err"):
            file_path = os.path.join(folder_path, filename)
            # Only process non-empty files
            if os.path.getsize(file_path) > 0:
                # Extract the slurm job ID from the filename
                match = re.match(r"(.+)-(\d+)_(\d+)\.err", filename)
                job_id = match.group(2)  # The job ID (numeric value after the '-')
                task_id = match.group(3)  # The task ID (numeric value after the '_')
                slurm_job_id = f"{job_id}_{task_id}"  # Combined job identifier
                
                # Read the content of the error file
                with open(file_path, 'r') as f:
                    content = f.read()

                # Check if the job failed 
                if "error" in content.lower():
                    # Determine the error type 
                    if "OOM Killed" in content:
                        error_type = "Out of Memory"
                    elif "DUE TO TIME LIMIT" in content:
                        error_type = "Timeout"
                    else:
                        error_type = "Other"
                    
                    # Update the job_errors dictionary
                    job_errors[slurm_job_id] = error_type
                    error_types[error_type] += 1

    # Return a summary of failed jobs and error types
    return job_errors, error_types

def print_slurm_error_summary(experiment_name):
    job_errors, error_types = parse_slurm_error_files(experiment_name)
    
    if not job_errors:
        print(f"Experiment {experiment_name}: error log's are emtpy.")
    
    else:
        # Print summary
        print("Failed Job Summary:")
        for job_id, error_type in job_errors.items():
            print(f"Job ID: {job_id}, Error: {error_type}")
        
        print("\nError Type Summary:")
        for error_type, count in error_types.items():
            print(f"{error_type}: {count}")

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to get summary of slurm errors")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    print_slurm_error_summary(args.experiment_name)
    
if __name__ == "__main__":
    main()
    