#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np

from experiments.experiment_script_runner import run_experiment_script
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

def generate_nsr_settings(min_nsr, max_nsr, steps, n_seeds):
    # Log-spaced values for NSR
    nsr_values = np.logspace(np.log10(min_nsr), np.log10(max_nsr), num=steps)
    nsr_list = [float(x) for x in nsr_values]
    
    print(f"\n[Setup] Generated {steps} NSR values from {min_nsr:.2e} to {max_nsr:.2e}")

    def user_set_experiment_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        fixed_params = {
            'infected_individuals_at_start': 10,
            'final_time': 365 * 3,
            'R': 1.05,
            'sequencing_rate': 1
        }
        return (varying_params, fixed_params, n_seeds)

    return user_set_experiment_settings
    
import matplotlib.pyplot as plt
import seaborn as sns
import simplicity.dir_manager as dm
import os 

#def plot_osr_df_scatter(df, parameter_name, experiment_name, target_host_type=None):
#    if df.empty:
#        print("No data available to plot.")
#        return
#
#    plt.figure(figsize=(10, 6))
#    
#    # We use 'hue' to color-code based on the outlier flag
#    # 'palette' allows you to choose specific colors (e.g., blue for 0, red for 1)
#    sns.stripplot(
#        data=df, 
#        x=parameter_name, 
#        y='observed_substitution_rate',
#        hue='is_outlier', 
#        palette={0: 'dodgerblue', 1: 'crimson'}, 
#        jitter=0, # Added back some jitter to see the density of the clusters
#        alpha=0.6, 
#        edgecolor='orange', 
#        linewidth=0.5
#    )
#    
#    # Optional: Adjust legend to be more descriptive
#    handles, labels = plt.gca().get_legend_handles_labels()
#    plt.legend(handles, ['Normal', 'Outlier'], title='Status', loc='best')
#    
#    # Labels and Titles
#    title = f"Distribution of OSR across {parameter_name}"
#    if target_host_type:
#        title += f" (Host: {target_host_type})"
#        
#    plt.title(title, fontsize=14)
#    plt.xlabel(f"Parameter: {parameter_name}", fontsize=12)
#    plt.ylabel("Observed Substitution Rate", fontsize=12)
#    plt.grid(axis='y', linestyle='--', alpha=0.4)
#    
#    plt.xticks(rotation=45)
#    plt.xscale("log")
#    plt.yscale("log")
#    
#    plt.tight_layout()
#
#    # File Path Logic
#    type_suffix = f"_{target_host_type}" if target_host_type else ""
#    base_filename = f"OSR_scatter_{experiment_name}_{parameter_name}{type_suffix}.png"
#    
#    plot_dir = dm.get_experiment_plots_dir(experiment_name)
#    if not os.path.exists(plot_dir):
#        os.makedirs(plot_dir)
#        
#    save_path = os.path.join(plot_dir, base_filename)
#    
#    plt.savefig(save_path, dpi=300, bbox_inches='tight')
#    plt.close() 
#    
#    print(f"DF plot saved as: {save_path}")

from matplotlib.lines import Line2D # Needed for the custom legend

def plot_osr_df_scatter(df, parameter_name, experiment_name, target_host_type=None):
    if df.empty:
        print("No data available to plot.")
        return

    plt.figure(figsize=(12, 7))
    
    # 1. PREPARE CATEGORICAL DATA
    # Create a string version of the parameter to force uniform spacing (Categorical Axis)
    # We sort by the numeric value first to ensure the x-axis order is correct
    df_sorted = df.sort_values(by=parameter_name)
    
    # Format to remove excessive zeros for cleaner labels (e.g., 0.001 instead of 0.001000)
    # This also converts it to string, making the x-axis categorical
    df_sorted['param_str'] = df_sorted[parameter_name].apply(lambda x: f"{x:.2g}")
    
    # 2. DEFINE HIGH-CONTRAST PALETTE
    # 'tab20' has 20 distinct colors. If you have >20 parameters, it repeats.
    # It offers much better distinction for neighbors than a gradient like 'viridis'.
    unique_params = df_sorted['param_str'].unique()
    n_colors = len(unique_params)
    palette = sns.color_palette("tab20", n_colors=n_colors)
    param_to_color = dict(zip(unique_params, palette))

    # 3. Plot NORMAL points (Circles)
    sns.stripplot(
        data=df_sorted[df_sorted['is_outlier'] == 0], 
        x='param_str',           # Plotting against the STRING column
        y='observed_substitution_rate',
        hue='param_str',         # Color by the same STRING column
        palette=param_to_color,
        marker='o', 
        jitter=0.25,             # Uniform jitter in every "lane"
        alpha=0.7, 
        edgecolor='gray', 
        linewidth=0.5,
        legend=False             
    )

    # 4. Plot OUTLIER points (Crosses)
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
    
    # 5. Custom Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Normal',
               markerfacecolor='gray', markersize=8),
        Line2D([0], [0], marker='X', color='w', label='Outlier',
               markeredgecolor='gray', markersize=8, linestyle='None')
    ]
    plt.legend(handles=legend_elements, loc='best', title='Status')
    
    # 6. Formatting
    title = f"Distribution of OSR across {parameter_name}"
    if target_host_type:
        title += f" (Host: {target_host_type})"
        
    plt.title(title, fontsize=14)
    plt.xlabel(f"Parameter: {parameter_name}", fontsize=12)
    plt.ylabel("Observed Substitution Rate", fontsize=12)
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    
    # IMPORTANT: Do NOT set xscale to 'log' because we are now plotting categories (strings).
    # The string labels preserve the order because we sorted df_sorted earlier.
    plt.yscale("log") 
    
    plt.xticks(rotation=45)
    plt.tight_layout()

    # File Path Logic
    type_suffix = f"_{target_host_type}" if target_host_type else ""
    base_filename = f"OSR_scatter_{experiment_name}_{parameter_name}{type_suffix}.png"
    
    plot_dir = dm.get_experiment_plots_dir(experiment_name)
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
        
    save_path = os.path.join(plot_dir, base_filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close() 
    
    print(f"DF plot saved as: {save_path}")

def plot_all_outlier_simulations(df, threshold):
    """
    Iterates through the combined DataFrame and generates a simulation plot
    for every seed marked as an outlier (is_outlier == 1).
    """
    # 1. Filter the DataFrame to get only outliers
    outlier_df = df[df['is_outlier'] == 1]
    
    if outlier_df.empty:
        print("No outliers detected in the provided DataFrame. Nothing to plot.")
        return

    print(f"Detected {len(outlier_df)} outliers. Generating simulation plots...")

    # 2. Iterate through the outlier paths
    for index, row in outlier_df.iterrows():
        ssod = row['ssod_path']
        param_val = row.get(next(iter(df.columns)), "N/A") # Get the first column (parameter)
        
        try:
            print(f"--> Plotting outlier: {os.path.basename(ssod)} (Param: {param_val})")
            
            # Call your existing function from the plot_manager (pm)
            pm.plot_simulation(ssod, threshold)
            
        except Exception as e:
            print(f"Error plotting {ssod}: {e}")

    print("Finished generating outlier plots.")
    
    
def perform_model_selection(experiment_numbered_name, OSR_filtered, HOST_TYPE):
    """
    Iterates through available models, fits them, and returns the best one based on AIC.
    """
    # 1. Define models to test (must match keys in er.factory_model_func)
    models_to_test = ['lin', 'log', 'exp', 'double_log', 'tan', 'spline']
    
    best_result = None
    best_model_name = ""
    best_aic = float('inf')

    print(f"\n[Model Selection] Comparing {len(models_to_test)} models...")
    print("-" * 60)
    print(f"{'Model':<15} | {'AIC':<15} | {'Status':<15}")
    print("-" * 60)

    for model_type in models_to_test:
        try:
            # Fit the model
            # Note: er.fit_... prints a report and saves a generic CSV internally. 
            # We will ignore the internal save and save the specific HOST_TYPE version below.
            fit_result = er.fit_observed_substitution_rate_regressor(
                experiment_numbered_name, OSR_filtered, model_type, weights=None
            )
            
            # Save the specific result for this model variant (with host type)
            om.write_fit_results_csv(experiment_numbered_name, model_type, fit_result, target_host_type=HOST_TYPE)

            aic = fit_result.aic
            status = "Success"
            
            # Check for best
            if aic < best_aic:
                best_aic = aic
                best_result = fit_result
                best_model_name = model_type
                status += " (New Best)"

            print(f"{model_type:<15} | {aic:<15.2f} | {status:<15}")

        except Exception as e:
            print(f"{model_type:<15} | {'N/A':<15} | Failed: {str(e)[:20]}...")

    print("-" * 60)
    
    if best_result is None:
        raise RuntimeError("All models failed to fit.")

    print(f"[Model Selection] WINNER: '{best_model_name}' (AIC: {best_aic:.2f})")
    return best_result, best_model_name


def perform_final_plotting(experiment_numbered_name, fit_result, model_type, 
                           OSR_clean, OSR_combined, OSR_mean, 
                           min_seq, min_len, HOST_TYPE, parameter):
    """
    Generates the final plots using the selected best model.
    """
    print('[Plotting] Generating Fit Figure...')
    pm.plot_OSR_fit_figure(
        experiment_numbered_name, fit_result, OSR_clean, OSR_combined, OSR_mean,
        model_type, min_seq, min_len
    )
    print(f"[Plotting] Fit Figure saved for model '{model_type}'.")
    
    # Combined Tempest Regression Plot
    if not OSR_combined.empty:
        print('[Plotting] Generating Single Regressions Plot...')
        y_axis_max = max(OSR_combined['observed_substitution_rate']) * 1.2
        
        pm.plot_combined_tempest_regressions(
            experiment_numbered_name, parameter, min_seq, min_len, y_axis_max=y_axis_max
        )
        print("[Plotting] Single Regressions Plot saved.")
    else:
        print("[Plotting] Skipped Single Regressions (Combined DF empty).")


def run_fitting(experiment_name, experiment_number, parameter='nucleotide_substitution_rate', 
                      min_seq=30, min_len=365): # Removed model_type (auto-selected now)
    
    experiment_numbered_name = f'{experiment_name}_#{experiment_number}'
    HOST_TYPE = 'normal'
    
    print(f"\n[Fitting] Starting OSR fit for: {experiment_name} (Type: {HOST_TYPE})")

    try:
        # ------------------------------------------------------------------
        # 1. Data Preparation & Outlier Detection
        # ------------------------------------------------------------------
        # A. Detect Outliers (creates Individual CSV)
        om.write_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len, 
                                      target_host_type=HOST_TYPE)
        
        # B. Create Combined Data (Filters out the outliers detected above)
        om.write_combined_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len, 
                                               target_host_type=HOST_TYPE, include_outliers=False)
        
        # C. Load Data for Debug Plotting (All points)
        OSR_all = om.read_OSR_vs_parameter_csv(experiment_numbered_name, parameter, min_seq, min_len, 
                                                  target_host_type=HOST_TYPE, include_outliers=True)
                                                  
        # [DEBUG] Plot with outliers visible
        plot_osr_df_scatter(df=OSR_all, 
                              parameter_name=parameter, 
                              experiment_name=experiment_numbered_name, 
                              target_host_type=HOST_TYPE
                          )
        #plot_all_outlier_simulations(OSR_all,0)
        
        if OSR_all.empty:
            print(f"[Warning] No data found for {experiment_name}.")
            return

        # D. Load Clean Data for Fitting
       

        # ------------------------------------------------------------------
        # 2. Model Selection
        # ------------------------------------------------------------------
        best_fit_result, best_model_type = perform_model_selection(
            experiment_numbered_name, OSR_all, HOST_TYPE
        )

        # ------------------------------------------------------------------
        # 3. Final Plotting
        # ------------------------------------------------------------------
        # Load auxiliary datasets needed for plotting
        OSR_combined = om.read_combined_OSR_vs_parameter_csv(
            experiment_numbered_name, parameter, min_seq, min_len, target_host_type=HOST_TYPE
        )
        OSR_mean = om.get_mean_std_OSR(
            experiment_numbered_name, parameter, min_seq, min_len, target_host_type=HOST_TYPE, include_outliers=False
        )

        perform_final_plotting(
            experiment_numbered_name, best_fit_result, best_model_type,
            OSR_all, OSR_combined, OSR_mean,
            min_seq, min_len, HOST_TYPE, parameter
        )

    except Exception as e:
        print(f"[Error] Fitting Pipeline failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    parser = argparse.ArgumentParser(description="Simplicity Experiment Manager Pipeline")
    parser.add_argument('--name', type=str, required=True, help="Experiment Name")
    parser.add_argument('--runner', type=str, required=True, help="slurm, serial, or multiprocessing")
    parser.add_argument('--exp-num', type=int, default=1)
    parser.add_argument('--min', type=float, required=True, help="Min NSR")
    parser.add_argument('--max', type=float, required=True, help="Max NSR")
    parser.add_argument('--steps', type=int, default=15)
    parser.add_argument('--seeds', type=int, default=100)

    args = parser.parse_args()

    # 1. Settings
    settings_func = generate_nsr_settings(args.min, args.max, args.steps, args.seeds)

    # 2. Run Simulations
    print(f"\n[Manager] Launching Simulations for {args.name}...")
    run_experiment_script(args.runner, args.exp_num, settings_func, args.name)

    # 3. Run Fitting 
    print("\n[Manager] Simulations complete. Starting Fitting...")
    run_fitting(args.name, args.exp_num)
    
if __name__ == "__main__":
    main()