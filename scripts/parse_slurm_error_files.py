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
                match = re.match(r"(.+)-(\d+)_(\d+)\.out", filename)
                if match:
                    job_id = match.group(2)  # The job ID (numeric value after the '-')
                    task_id = match.group(3)  # The task ID (numeric value after the '_')
                    slurm_job_id = f"{job_id}_{task_id}"  # Combined job identifier
                    
                    # Read the content of the error file
                    with open(file_path, 'r') as f:
                        content = f.read()

                    # Check if the job failed (this can be adjusted depending on the content format)
                    if "error" in content.lower() or "failed" in content.lower():
                        # Determine the error type (simplified, adjust based on actual content format)
                        if "out of memory" in content.lower():
                            error_type = "Out of Memory"
                        elif "segmentation fault" in content.lower():
                            error_type = "Segmentation Fault"
                        elif "timeout" in content.lower():
                            error_type = "Timeout"
                        else:
                            error_type = "Other"
                        
                        # Update the job_errors dictionary
                        job_errors[slurm_job_id] = error_type
                        error_types[error_type] += 1

    # Return a summary of failed jobs and error types
    return job_errors, error_types

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to get summary of slurm errors")
    parser.add_argument('experiment_name', type=str, help="experiment name")
    args = parser.parse_args()
    
    job_errors, error_types = parse_slurm_error_files(args.experiment_name)
    
    if not job_errors:
        print(f"The experiment: {args.experiment_name} ran SUCCESSFULLY without errors.")
    
    else:
        # Print summary
        print("Failed Job Summary:")
        print(job_errors.items())
        for job_id, error_type in job_errors.items():
            print(f"Job ID: {job_id}, Error: {error_type}")
        
        print("\nError Type Summary:")
        for error_type, count in error_types.items():
            print(f"{error_type}: {count}")

if __name__ == "__main__":
    main()
    