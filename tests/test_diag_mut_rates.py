#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 13:34:47 2026

@author: pietro
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from simplicity.runme import run_experiment
import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import simplicity.runners.serial
import simplicity.runners.multiprocessing
import simplicity.runners.slurm
from test_local_runme import test_experiment_output 

def generate_settings(args):
    """Generates the settings function needed by run_experiment based on CLI args."""
    def user_set_experiment_settings():
        varying_params = {}
        fixed_params = {
            'population_size': 1000,
            'long_shedders_ratio': args.ls_ratio,
            'infected_individuals_at_start': 100,
            
            'nucleotide_substitution_rate': args.nsr,
            'M_nsr_long': args.m_nsr_long,
        
            'R': 1.03,
            'R_long': 3.0,
            
            'diagnosis_rate_standard': args.diag_std,
            'diagnosis_rate_long': args.diag_long,
            
            'final_time': 600,
            'sequence_long_shedders': True
        }
        return (varying_params, fixed_params, args.seeds)
    return user_set_experiment_settings

def plot_host_type_regressions(sim_dir, seq_df_std, fit_std, seq_df_long, fit_long):
    """Generates a side-by-side plot of Tempest regressions for both host types."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True, sharex=True)
    
    # Standard Hosts Plot
    if seq_df_std is not None and not seq_df_std.empty:
        sns.scatterplot(data=seq_df_std, x='Sequencing_time', y='Distance_from_root', 
                        color='blue', alpha=0.5, ax=axes[0], label='Sequences')
        x_std = seq_df_std['Sequencing_time'].values.reshape(-1, 1)
        axes[0].plot(x_std, fit_std.predict(x_std), color='orange', 
                     linewidth=2, label=f'OSR: {fit_std.coef_[0]:.6f}')
        axes[0].set_title('Standard Hosts Regression')
        axes[0].set_xlabel('Simulation time (y)')
        axes[0].set_ylabel('Genetic distance from root')
        axes[0].legend()

    # Long Shedders Plot
    if seq_df_long is not None and not seq_df_long.empty:
        sns.scatterplot(data=seq_df_long, x='Sequencing_time', y='Distance_from_root', 
                        color='red', alpha=0.5, ax=axes[1], label='Sequences')
        x_long = seq_df_long['Sequencing_time'].values.reshape(-1, 1)
        axes[1].plot(x_long, fit_long.predict(x_long), color='black', 
                     linewidth=2, label=f'OSR: {fit_long.coef_[0]:.6f}')
        axes[1].set_title('Long Shedders Regression')
        axes[1].set_xlabel('Simulation time (y)')
        axes[1].legend()

    plt.tight_layout()
    
    # Save Plot
    sim_name = os.path.basename(sim_dir)
    exp_name = dm.get_experiment_foldername_from_SSOD(sim_dir)
    plots_dir = dm.get_experiment_plots_dir(exp_name)
    out_path = os.path.join(plots_dir, f"{sim_name}_host_type_regressions.png")
    
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"    [Plot] Saved regression comparison to: {out_path}")


def check_rates_by_host_type(experiment_name, do_plot=False):
    """Calculates and verifies the OSR and Effective Diagnosis rates by host type."""
    print('\n====================================================')
    print(' VERIFYING MUTATION AND DIAGNOSIS RATES BY HOST TYPE')
    print('====================================================')
    
    sim_dirs = dm.get_simulation_output_dirs(experiment_name)

    for sim_dir in sim_dirs:
        print(f"\n--- Analyzing Simulation: {os.path.basename(sim_dir)} ---")
        
        # --- 1. CHECK OBSERVED SUBSTITUTION RATES (OSR) ---
        seq_df_std, fit_std = None, None
        seq_df_long, fit_long = None, None
        
        try:
            seq_df_std = om.create_combined_sequencing_df(sim_dir, individual_type='standard')
            if seq_df_std is not None and not seq_df_std.empty:
                fit_std = er.tempest_regression(seq_df_std)
                print(f"  [Mutations] Standard OSR:     {fit_std.coef_[0]:.6f} subs/site/year")
            else:
                print("  [Mutations] Standard OSR:     Not enough data")

            seq_df_long = om.create_combined_sequencing_df(sim_dir, individual_type='long_shedder')
            if seq_df_long is not None and not seq_df_long.empty:
                fit_long = er.tempest_regression(seq_df_long)
                print(f"  [Mutations] Long Shedder OSR: {fit_long.coef_[0]:.6f} subs/site/year")
            else:
                print("  [Mutations] Long Shedder OSR: Not enough data")
                
            if do_plot and fit_std and fit_long:
                plot_host_type_regressions(sim_dir, seq_df_std, fit_std, seq_df_long, fit_long)
                
        except Exception as e:
            print(f"  [Mutations] Error calculating OSR: {e}")

        # --- 2. CHECK EFFECTIVE DIAGNOSIS RATES ---
        try:
            seed_dirs = dm.get_seeded_simulation_output_dirs(sim_dir)
            for seed_dir in seed_dirs:
                ind_file = os.path.join(seed_dir, 'individuals_data.csv')
                if os.path.exists(ind_file):
                    df = pd.read_csv(ind_file)
                    
                    # Standard Hosts
                    std_df = df[df['type'] == 'standard']
                    std_diag = len(std_df[std_df['state'] == 'diagnosed'])
                    std_rec = len(std_df[std_df['state'] == 'recovered'])
                    std_dec = len(std_df[std_df['state'] == 'deceased']) if 'deceased' in std_df['state'].values else 0
                    
                    std_resolved = std_diag + std_rec + std_dec
                    std_eff_rate = std_diag / std_resolved if std_resolved > 0 else 0
                    
                    # Long Shedders
                    long_df = df[df['type'] == 'long_shedder']
                    long_diag = len(long_df[long_df['state'] == 'diagnosed'])
                    long_rec = len(long_df[long_df['state'] == 'recovered'])
                    long_dec = len(long_df[long_df['state'] == 'deceased']) if 'deceased' in long_df['state'].values else 0
                    
                    long_resolved = long_diag + long_rec + long_dec
                    long_eff_rate = long_diag / long_resolved if long_resolved > 0 else 0
                    
                    print(f"  [Diagnosis] Seed: {os.path.basename(seed_dir)}")
                    print(f"      Standard Eff. Rate:     {std_eff_rate:.4f} (from {std_resolved} resolved)")
                    print(f"      Long Shedder Eff. Rate: {long_eff_rate:.4f} (from {long_resolved} resolved)")
        except Exception as e:
            print(f"  [Diagnosis] Error calculating effective rates: {e}")
            
    print('====================================================\n')


def main():
    parser = argparse.ArgumentParser(description="Unified Pipeline Test for SIMPLICITY")
    
    # General execution arguments
    parser.add_argument('--runner', type=str, choices=['serial', 'multiprocessing', 'slurm'], default='serial', help="Runner to use")
    parser.add_argument('--test-num', type=int, default=1, help="Identifier for the experiment name")
    parser.add_argument('--seeds', type=int, default=2, help="Number of simulated seeds")
    parser.add_argument('--plot', action='store_true', help="Generate Tempest regression plots by host type")
    
    # Biological parameters
    parser.add_argument('--nsr', type=float, default=0.0001, help="Base Nucleotide Substitution Rate")
    parser.add_argument('--m-nsr-long', type=float, default=1, help="Multiplier for long shedder mutation rate")
    parser.add_argument('--ls-ratio', type=float, default=0.1, help="Ratio of long shedders in population")
    parser.add_argument('--diag-std', type=float, default=0.1, help="Theoretical diagnosis rate for standard hosts")
    parser.add_argument('--diag-long', type=float, default=0.1, help="Theoretical diagnosis rate for long shedders")

    args = parser.parse_args()

    # Determine Runner Module
    if args.runner == 'serial':
        runner_module = simplicity.runners.serial
    elif args.runner == 'multiprocessing':
        runner_module = simplicity.runners.multiprocessing
    elif args.runner == 'slurm':
        runner_module = simplicity.runners.slurm

    experiment_name = f'test_pipeline_{args.runner}_#{args.test_num}'
    settings_func = generate_settings(args)

    print(f"\n[INIT] Launching {experiment_name} using {args.runner.upper()} runner...")
    start = time.time()
    
    # 1. Run Experiment
    run_experiment(experiment_name, 
                   settings_func,             
                   simplicity_runner = runner_module,
                   archive_experiment = False)
    
    # Note: For SLURM, the script submits jobs and might exit immediately. 
    # The following checks assume the data has been generated (i.e. local runs).
    if args.runner != 'slurm':
        # 2. Check output file integrity
        print("\n[VERIFY] Checking output files integrity...")
        test_experiment_output(experiment_name)
        
        # 3. Check specific biology constraints (Mutations & Diagnosis)
        check_rates_by_host_type(experiment_name, do_plot=args.plot)

        elapsed = time.time() - start
        mins, secs = divmod(elapsed, 60)
        print(f"Pipeline test completed in {int(mins)} min {secs:.2f} sec")
    else:
        print("\n[INFO] Jobs submitted to SLURM. Once they finish, run the rate checks manually.")

if __name__ == "__main__":
    main()