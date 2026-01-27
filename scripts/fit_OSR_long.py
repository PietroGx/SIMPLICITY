#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib.lines import Line2D

from experiments.experiment_script_runner import run_experiment_script
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import simplicity.dir_manager as dm


# ==============================================================================
# 1. SETTINGS GENERATORS
# ==============================================================================

def generate_nsr_settings(min_val, max_val, steps, n_seeds):
    """
    Stage 1: Vary NSR, fixed standard parameters.
    """
    # Log-spaced values for NSR
    nsr_values = np.logspace(np.log10(min_val), np.log10(max_val), num=steps)
    nsr_list = [float(x) for x in nsr_values]
    
    print(f"\n[Setup] Generated {steps} NSR values from {min_val:.2e} to {max_val:.2e}")

    def user_set_experiment_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        fixed_params = {
            'infected_individuals_at_start': 100,
            'final_time': 365 * 3,
            'R': 1.05,
            'sequencing_rate': 0.1
        }
        return (varying_params, fixed_params, n_seeds)

    return user_set_experiment_settings


def generate_ls_factor_settings(calibrated_nsr, min_val, max_val, steps, n_seeds):
    """
    Stage 2: Vary 'long_evo_rate_f', fixed NSR (calibrated).
    """
    # Linear or Log space? Evolution factors usually behave linearly or log-linearly.
    # We'll use linspace for the factor (e.g. 1x to 10x) unless orders of magnitude are needed.
    # If the range is large (1 to 100), logspace is better. 
    # Assuming 1 to 10 -> Linspace is usually fine, but let's stick to logspace 
    # to be consistent with the previous step structure or if the user passes large ranges.
    
    if max_val / min_val > 20:
         vals = np.logspace(np.log10(min_val), np.log10(max_val), num=steps)
    else:
         vals = np.linspace(min_val, max_val, num=steps)
         
    factor_list = [float(x) for x in vals]

    print(f"\n[Setup] Generated {steps} LS Factors from {min_val:.2f} to {max_val:.2f}")
    print(f"[Setup] Fixed NSR: {calibrated_nsr:.4e}")

    def user_set_experiment_settings():
        varying_params = {'long_evo_rate_f': factor_list}
        fixed_params = {
            'nucleotide_substitution_rate': calibrated_nsr,
            'infected_individuals_at_start': 100, 
            'final_time': 365 * 3,
            'R': 1.05,
            'sequencing_rate': 1,
            'long_shedders_ratio': 0.01,
            'sequence_long_shedders':True
        }
        return (varying_params, fixed_params, n_seeds)

    return user_set_experiment_settings


# ==============================================================================
# 2. PLOTTING UTILS
# ==============================================================================

def plot_osr_df_scatter(df, parameter_name, experiment_name, individual_type=None):
    if df.empty:
        print("No data available to plot.")
        return

    plt.figure(figsize=(12, 7))
    
    # Sort for plotting
    df_sorted = df.sort_values(by=parameter_name)
    df_sorted['param_str'] = df_sorted[parameter_name].apply(lambda x: f"{x:.2g}")
    
    unique_params = df_sorted['param_str'].unique()
    n_colors = len(unique_params)
    palette = sns.color_palette("tab20", n_colors=n_colors)
    param_to_color = dict(zip(unique_params, palette))

    # Plot Normal
    sns.stripplot(
        data=df_sorted[df_sorted['is_outlier'] == 0], 
        x='param_str',           
        y='observed_substitution_rate',
        hue='param_str',         
        palette=param_to_color,
        marker='o', 
        jitter=0.25,             
        alpha=0.7, 
        edgecolor='gray', 
        linewidth=0.5,
        legend=False             
    )

    # Plot Outliers
    sns.stripplot(
        data=df_sorted[df_sorted['is_outlier'] == 1], 
        x='param_str', 
        y='observed_substitution_rate',
        hue='param_str',      
        palette=param_to_color,  
        marker='X', 
        s=8, 
        jitter=0.25, 
        alpha=1.0, 
        linewidth=1.5,
        legend=False             
    )
    
    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Normal',
               markerfacecolor='gray', markersize=8),
        Line2D([0], [0], marker='X', color='w', label='Outlier',
               markeredgecolor='gray', markersize=8, linestyle='None')
    ]
    plt.legend(handles=legend_elements, loc='best', title='Status')
    
    # Titles & Labels
    title = f"Distribution of OSR across {parameter_name}"
    if individual_type:
        title += f" (Host: {individual_type})"
        
    plt.title(title, fontsize=14)
    plt.xlabel(f"Parameter: {parameter_name}", fontsize=12)
    plt.ylabel("Observed Substitution Rate", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    
    plt.yscale("log") 
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save
    type_suffix = f"_{individual_type}" if individual_type else ""
    base_filename = f"OSR_scatter_{experiment_name}_{parameter_name}{type_suffix}.png"
    
    plot_dir = dm.get_experiment_plots_dir(experiment_name)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
        
    save_path = os.path.join(plot_dir, base_filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close() 
    print(f"DF plot saved as: {save_path}")


# ==============================================================================
# 3. FITTING LOGIC (Reusable Step)
# ==============================================================================

def perform_model_selection(experiment_numbered_name, OSR_filtered, individual_type, parameter_name):
    """
    Fits available models and returns the best result.
    """
    models_to_test = ['lin', 'log', 'exp', 'tan'] # Spline removed for now
    
    best_result = None
    best_model_name = ""
    best_aic = float('inf')

    print(f"\n[Model Selection] Comparing {len(models_to_test)} models for {individual_type}...")
    print(f"{'Model':<10} | {'AIC':<15} | {'Status':<15}")
    print("-" * 45)

    for model_type in models_to_test:
        try:
            # Note: We pass parameter_name so the fitter knows which column is X
            fit_result = er.fit_observed_substitution_rate_regressor(
                experiment_numbered_name, 
                OSR_filtered, 
                model_type, 
                weights=None,
                parameter_name=parameter_name 
            )
            
            # Save specific CSV
            om.write_fit_results_csv(experiment_numbered_name, model_type, fit_result, individual_type=individual_type)

            aic = fit_result.aic
            status = "Success"
            
            if aic < best_aic:
                best_aic = aic
                best_result = fit_result
                best_model_name = model_type
                status += " (*)"

            print(f"{model_type:<10} | {aic:<15.2f} | {status:<15}")

        except Exception as e:
            # print(e) # Debug
            print(f"{model_type:<10} | {'N/A':<15} | Failed")

    print("-" * 45)
    
    if best_result is None:
        raise RuntimeError("All models failed to fit.")

    print(f"[Result] Best Model: '{best_model_name}' (AIC: {best_aic:.2f})")
    return best_result, best_model_name


def run_fitting_step(experiment_name, experiment_number, parameter_name, individual_type, 
                     min_seq=30, min_len=365):
    """
    Executes the fitting pipeline for a single stage.
    Returns (best_fit_result, best_model_type)
    """
    experiment_numbered_name = f'{experiment_name}_#{experiment_number}'
    
    print(f"\n[Fitting] Starting analysis for: {experiment_name}")
    print(f"          Target Host: {individual_type}")
    print(f"          Parameter:   {parameter_name}")

    try:
        # A. Detect Outliers 
        om.write_OSR_vs_parameter_csv(experiment_numbered_name, parameter_name, min_seq, min_len, 
                                      individual_type=individual_type)
        
        # B. Create Combined Data (All)
        om.write_combined_OSR_vs_parameter_csv(experiment_numbered_name, parameter_name, min_seq, min_len, 
                                               individual_type=individual_type, include_outliers=True)
        
        # C. Load Data for Plotting (All)
        OSR_all = om.read_OSR_vs_parameter_csv(experiment_numbered_name, parameter_name, min_seq, min_len, 
                                                  individual_type=individual_type, include_outliers=True)
                                                  
        plot_osr_df_scatter(df=OSR_all, 
                            parameter_name=parameter_name, 
                            experiment_name=experiment_numbered_name, 
                            individual_type=individual_type)
        
        if OSR_all.empty:
            print(f"[Warning] No data found for {experiment_name}.")
            return None, None

        # D. Model Selection (Fitting)
        best_fit_result, best_model_type = perform_model_selection(
            experiment_numbered_name, OSR_all, individual_type, parameter_name
        )

        # E. Final Plotting
        OSR_combined = om.read_combined_OSR_vs_parameter_csv(
            experiment_numbered_name, parameter_name, min_seq, min_len, individual_type=individual_type
        )
        OSR_mean = om.get_mean_std_OSR(
            experiment_numbered_name, parameter_name, min_seq, min_len, individual_type=individual_type, include_outliers=True
        )

        print('[Plotting] Generating Fit Figures...')
        pm.plot_OSR_fit_figure(
            experiment_numbered_name, parameter_name, best_fit_result, OSR_all, OSR_combined, OSR_mean,
            best_model_type, min_seq, min_len
        )
        pm.plot_combined_tempest_regressions(experiment_numbered_name, parameter_name, 
                                      min_seq, min_len, individual_type)
        
        return best_fit_result, best_model_type

    except Exception as e:
        print(f"[Error] Fitting Step failed: {e}")
        import traceback
        traceback.print_exc()
        return None, None


# ==============================================================================
# 4. MAIN PIPELINE
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Simplicity Two-Stage Calibration Pipeline")
    
    # Experiment Config
    parser.add_argument('--name', type=str, required=True, help="Base Experiment Name")
    parser.add_argument('--runner', type=str, required=True, help="slurm, serial, or multiprocessing")
    parser.add_argument('--exp-num', type=int, default=1)
    parser.add_argument('--seeds', type=int, default=100)
    
    # Stage 1: NSR Sweep
    parser.add_argument('--nsr-min', type=float, required=True, help="Min NSR")
    parser.add_argument('--nsr-max', type=float, required=True, help="Max NSR")
    parser.add_argument('--nsr-steps', type=int, default=15)
    
    # Stage 2: Factor Sweep
    parser.add_argument('--factor-min', type=float, default=1.0, help="Min Long Shedder Factor")
    parser.add_argument('--factor-max', type=float, default=1000.0, help="Max Long Shedder Factor")
    parser.add_argument('--factor-steps', type=int, default=15)

    # Targets
    parser.add_argument('--target-osr-normal', type=float, required=True, help="Target OSR for Normal Hosts")
    parser.add_argument('--target-osr-ls', type=float, required=True, help="Target OSR for Long Shedders")

    args = parser.parse_args()

    # --------------------------------------------------------------------------
    # STAGE 1: CALIBRATE NSR (Normal Hosts)
    # --------------------------------------------------------------------------
    print("\n" + "="*60)
    print("STAGE 1: Calibrating Baseline NSR (Normal Hosts)")
    print("="*60)
    
    exp_name_1 = f"{args.name}_Step1"
    settings_1 = generate_nsr_settings(args.nsr_min, args.nsr_max, args.nsr_steps, args.seeds)
    
    # 1.1 Run Simulation
    print(f"\n[Manager] Launching Simulations for {exp_name_1}...")
    run_experiment_script(args.runner, args.exp_num, settings_1, exp_name_1)
    
    # 1.2 Run Fitting
    fit_res_1, model_1 = run_fitting_step(
        exp_name_1, args.exp_num, 'nucleotide_substitution_rate', 'normal'
    )
    
    if not fit_res_1:
        print("[Error] Stage 1 fitting failed. Aborting.")
        return

    # 1.3 Calculate Calibrated NSR
    calibrated_nsr = er.compute_calibrated_parameter(model_1, fit_res_1, args.target_osr_normal)
    print(f"\n>> STAGE 1 COMPLETE")
    print(f">> Target OSR (Normal): {args.target_osr_normal}")
    print(f">> Calibrated NSR:      {calibrated_nsr:.6e}")
    
    
    # --------------------------------------------------------------------------
    # STAGE 2: CALIBRATE LS FACTOR (Long Shedders)
    # --------------------------------------------------------------------------
    print("\n" + "="*60)
    print("STAGE 2: Calibrating Long Shedder Factor")
    print("="*60)
    
    exp_name_2 = f"{args.name}_Step2"
    settings_2 = generate_ls_factor_settings(calibrated_nsr, args.factor_min, args.factor_max, 
                                             args.factor_steps, args.seeds)
    
    # 2.1 Run Simulation
    print(f"\n[Manager] Launching Simulations for {exp_name_2}...")
    run_experiment_script(args.runner, args.exp_num, settings_2, exp_name_2)
    
    # 2.2 Run Fitting
    fit_res_2, model_2 = run_fitting_step(
        exp_name_2, args.exp_num, 'long_evo_rate_f', 'long_shedder', min_seq=0
    )
    
    if not fit_res_2:
        print("[Error] Stage 2 fitting failed. Aborting.")
        return

    # 2.3 Calculate Calibrated Factor
    calibrated_factor = er.compute_calibrated_parameter(model_2, fit_res_2, args.target_osr_ls)
    
    print(f"\n>> STAGE 2 COMPLETE")
    print(f">> Target OSR (LS):     {args.target_osr_ls}")
    print(f">> Calibrated Factor:   {calibrated_factor:.4f}")

    # --------------------------------------------------------------------------
    # FINAL REPORT
    # --------------------------------------------------------------------------
    print("\n" + "#"*60)
    print("CALIBRATION COMPLETE")
    print("#"*60)
    print(f"1. Baseline NSR:        {calibrated_nsr:.6e}")
    print(f"2. Long Shedder Factor: {calibrated_factor:.4f}")
    print("#"*60 + "\n")

if __name__ == "__main__":
    main()