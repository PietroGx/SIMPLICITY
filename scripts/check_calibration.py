#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Ensure the simplicity module can be imported when running from scripts/
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

def main():
    # --- Settings ---
    experiment_name = "calibrate_long_nsr_#2" # Update if your exp number is different
    target_osr = 0.00205
    min_seq = 30
    min_len = 100
    model_type = 'exp'
    
    print(f"--- Analyzing 2D Calibration: {experiment_name} ---")
    
    # ---------------------------------------------------------
    # PHASE 1 & 2: Custom 2D Data Extraction
    # ---------------------------------------------------------
    try:
        simulation_output_dirs = dm.get_simulation_output_dirs(experiment_name)
    except FileNotFoundError:
        print(f"Error: Experiment '{experiment_name}' not found.")
        sys.exit(1)

    all_sod_results = []
    print("Extracting parameters and computing OSR for all seeds (this may take a moment)...")
    
    for sod in simulation_output_dirs:
        nsr_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'nucleotide_substitution_rate')
        tau_val = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')
        
        seeded_dirs = dm.get_seeded_simulation_output_dirs(sod)
        sod_df_rows = []
        
        for ssod in seeded_dirs:
            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
                
                if final_time >= min_len and len(seq_data) >= min_seq:
                    # Calculate OSR for this seed
                    osr_val = er.tempest_regression(seq_data).coef_[0]
                    sod_df_rows.append({
                        'tau_3_long': tau_val,
                        'nucleotide_substitution_rate': nsr_val,
                        'observed_substitution_rate': osr_val,
                        'ssod_path': ssod
                    })
            except Exception:
                continue # Skip failed or incomplete seeds
        
        # Apply outlier detection per specific (NSR, tau_3_long) grid point
        if sod_df_rows:
            sod_df = pd.DataFrame(sod_df_rows)
            sod_df = om.detect_sod_outliers(sod_df)
            all_sod_results.append(sod_df)

    if not all_sod_results:
        print("No valid data found. Check if simulations finished and lengths/seqs meet minimums.")
        sys.exit(1)

    master_df = pd.concat(all_sod_results, ignore_index=True)
    clean_df = master_df[master_df['is_outlier'] == 0]
    
    # ---------------------------------------------------------
    # PHASE 3 & 4: Stratified Regression & Multi-Curve Plotting
    # ---------------------------------------------------------
    plt.figure(figsize=(10, 6))
    colors = {63.0: 'blue', 109.0: 'orange', 365.0: 'red'}
    labels = {63.0: 'SOT (63d)', 109.0: 'HIV (109d)', 365.0: 'Edge Case (365d)'}
    
    clinical_groups = sorted(clean_df['tau_3_long'].unique())
    results = {}
    
    print(f"\n--- Calibration Results (Target OSR = {target_osr}) ---")
    
    for tau in clinical_groups:
        group_df = clean_df[clean_df['tau_3_long'] == tau]
        color = colors.get(tau, 'black')
        label = labels.get(tau, f'tau={tau}')
        
        # Scatter the raw, clean seed data (faded)
        plt.scatter(group_df['nucleotide_substitution_rate'], 
                    group_df['observed_substitution_rate'], 
                    color=color, alpha=0.15, s=10)
        
        # Fit the regressor
        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_name, 
            group_df, 
            model_type=model_type, 
            parameter_name='nucleotide_substitution_rate'
        )
        
        # Interpolate the required input NSR
        try:
            calibrated_nsr = er.compute_calibrated_parameter(model_type, fit_result, target_osr)
            results[tau] = calibrated_nsr
            print(f"{label:<17}: NSR = {calibrated_nsr:.6f}")
        except Exception as e:
            print(f"{label:<17}: Could not compute calibration ({e})")
            continue
            
        # Plot the fit line
        x_vals = np.linspace(group_df['nucleotide_substitution_rate'].min(), 
                             group_df['nucleotide_substitution_rate'].max(), 100)
        y_vals = fit_result.eval(x=x_vals)
        plt.plot(x_vals, y_vals, color=color, linewidth=2, label=f"{label} fit")
        
        # Mark the intersection
        plt.plot(calibrated_nsr, target_osr, marker='*', markersize=12, color=color, markeredgecolor='black')

    # Formatting the Plot
    plt.axhline(target_osr, color='black', linestyle='--', linewidth=1.5, label='Target OSR')
    plt.title('Long-Shedder OSR Calibration by Clinical Group')
    plt.xlabel('Input Nucleotide Substitution Rate (NSR)')
    plt.ylabel('Observed Substitution Rate (OSR)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save output
    plot_path = os.path.join(dm.get_experiment_dir(experiment_name), 'multi_tau_calibration_fit.png')
    plt.savefig(plot_path, dpi=300)
    print(f"\nCalibration plot saved to: {plot_path}")

if __name__ == "__main__":
    main()
