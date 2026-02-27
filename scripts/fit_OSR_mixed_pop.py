#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import scipy.linalg

import simplicity
import simplicity.settings_manager as sm
import simplicity.tuning.diagnosis_rate as dr
import simplicity.plots_manager as pm 
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import simplicity.dir_manager as dm

from simplicity.runme import run_experiment
from experiments.experiment_script_runner import run_experiment_script

def generate_homogeneous_settings(min_nsr, max_nsr, steps, n_seeds):
    """
    Step 2 in LaTeX: Establish Calibration Curve (r=0)
    Generates a parameter sweep of NSR for a completely normal population.
    """
    nsr_values = np.logspace(np.log10(min_nsr), np.log10(max_nsr), num=steps)
    nsr_list = [float(x) for x in nsr_values]
    
    print(f"\n[Calibration] Generating {steps} baseline NSR values from {min_nsr:.2e} to {max_nsr:.2e}")

    def user_set_experiment_settings():
        varying_params = {'nucleotide_substitution_rate': nsr_list}
        fixed_params = {
            'long_shedders_ratio': 0.0,
            'infected_individuals_at_start': 100,
            'final_time': 365 * 3,
            'sequence_long_shedders': False
        }
        return (varying_params, fixed_params, n_seeds)

    return user_set_experiment_settings

def get_expected_lineages(t, k_v, cap_min, cap_max):
    """Calculates E[N_lin(t)] averaged over the uniform distribution of caps."""
    caps = np.arange(cap_min, cap_max + 1)
    # Lineages grow linearly, but cannot exceed the cap
    lineages = np.minimum(1 + k_v * t, caps)
    return np.mean(lineages)

def calculate_expected_lineage_days(tau_1, tau_2, tau_3, tau_4, k_d, k_v, cap_min, cap_max, dt=0.5):
    """
    Numerically integrates W = int N_lin(t) * S(t) dt 
    """
    # Get matrix B
    B = dr.get_B(tau_1, tau_2, tau_3, tau_4, k_d=k_d)
    
    # Setup the numerical integration
    W = 0.0
    p = np.zeros(21) # Note: get_B returns a 21x21 matrix (including the 0 column at the end)
    p[0] = 1.0 
    
    B_dt = scipy.linalg.expm(B * dt)
    t = 0.0
    
    while True:
        S_t = np.sum(p[0:20]) # Only sum the 20 transient infected states
        
        if S_t < 1e-6:
            break
            
        N_t = get_expected_lineages(t, k_v, cap_min, cap_max)
        W += N_t * S_t * dt
        
        p = B_dt @ p
        t += dt
        
    return W

def fit_calibration_curve(experiment_name, parameter='nucleotide_substitution_rate', 
                          min_seq=30, min_len=100, model_type='exp'):
    """Extracts OSR data from the sweep and fits the calibration curve f(NSR)."""
    print(f"\n[Fitting] Establishing calibration curve f(NSR) for: {experiment_name}")
    
    try:
        om.write_combined_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len)
        om.write_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len)
        
        OSR_single = om.read_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len)
        
        if OSR_single.empty:
            raise ValueError(f"No sequence data found for {experiment_name}. Cannot fit curve.")

        # Fit Model
        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_name, OSR_single, model_type, weights=None
        )

        # Plot for visual confirmation
        OSR_combined = om.read_combined_OSR_vs_parameter_csv(experiment_name, parameter, min_seq, min_len)
        OSR_mean = om.get_mean_std_OSR(experiment_name, parameter, min_seq, min_len)
        pm.plot_OSR_fit_figure(
            experiment_name=experiment_name, 
            parameter=parameter,                 
            fit_result=fit_result, 
            OSR_single_sim_data=OSR_single,          
            OSR_combined_sim_data=OSR_combined,      
            OSR_mean_data=OSR_mean,                  
            model_type=model_type, 
            min_seq_number=min_seq, 
            min_sim_lenght=min_len               
        )
        print(f"[Fitting] SUCCESS. Calibration curve '{model_type}' established.")
        return fit_result

    except Exception as e:
        print(f"[Error] Curve fitting failed: {e}")
        return None

def calculate_time_fractions(r, T_norm, T_long):
    """Calculates F_norm and F_long (Equations 10 & 11)."""
    denominator = ((1 - r) * T_norm) + (r * T_long)
    F_norm = ((1 - r) * T_norm) / denominator
    F_long = (r * T_long) / denominator
    return F_norm, F_long

def verify_mixed_simulation(runner_module, NSR_base, M, r, target_osr, n_seeds):
    """Runs a single mixed simulation using the derived parameters to verify the math."""
    print("\n[Verification] Launching verification simulation...")
    exp_name = "Verification_Mixed_OSR_Run"
    
    def verify_settings():
        varying_params = {}
        fixed_params = {
            'nucleotide_substitution_rate': NSR_base,
            'M_nsr_long': M,
            'long_shedders_ratio': r,
            'infected_individuals_at_start': 200,
            'final_time': 365 * 3,
            'sequence_long_shedders': True
        }
        return (varying_params, fixed_params, n_seeds) 

    run_experiment(exp_name, verify_settings, simplicity_runner=runner_module, archive_experiment=False)
    
    # Calculate Global OSR of the mixed run
    sim_dirs = dm.get_simulation_output_dirs(exp_name)
    if not sim_dirs: return
    
    seq_df_mixed = om.create_combined_sequencing_df(sim_dirs[0], individual_type=None) # None = Global Mix
    if seq_df_mixed is not None and not seq_df_mixed.empty:
        fit_mixed = er.tempest_regression(seq_df_mixed)
        actual_osr = fit_mixed.coef_[0]
        error = abs(actual_osr - target_osr) / target_osr * 100
        
        print("\n=========================================================")
        print(" VERIFICATION RESULTS")
        print("=========================================================")
        print(f" Target Empirical OSR : {target_osr:.6f}")
        print(f" Achieved Global OSR  : {actual_osr:.6f}")
        print(f" Margin of Error      : {error:.2f}%")
        print("=========================================================\n")


def main():
    parser = argparse.ArgumentParser(description="Heterogeneous Evolutionary Rate Calibration Pipeline")
    
    # Runner Settings
    parser.add_argument('--name', type=str, default="Calibration_Sweep", help="Experiment Name")
    parser.add_argument('--runner', type=str, choices=['serial', 'multiprocessing', 'slurm'], default='serial')
    parser.add_argument('--exp-num', type=int, default=1)
    
    # Sweep Parameters (Step 2)
    parser.add_argument('--min', type=float, default=1e-5, help="Min NSR for sweep")
    parser.add_argument('--max', type=float, default=1e-3, help="Max NSR for sweep")
    parser.add_argument('--steps', type=int, default=10, help="Number of NSR points")
    parser.add_argument('--seeds', type=int, default=100, help="Seeds per point")
    parser.add_argument('--model', type=str, default='exp', choices=['lin', 'log', 'exp', 'tan'], help="Fit model")
    
    # Target Empirical Parameters (Step 1, 3, & 4)
    parser.add_argument('--target-osr', type=float, required=True, help="Empirical mixed OSR to target (e.g., 0.0008)")
    parser.add_argument('--M', type=float, required=True, help="Empirical rate multiplier (OSR_long / OSR_norm)")
    parser.add_argument('--r', type=float, default=0.01, help="Biological prevalence of long shedders (default 0.01)")
    
    parser.add_argument('--diag-std', type=float, default=0.1, help="Diagnosis rate for standard hosts")
    parser.add_argument('--diag-long', type=float, default=0.1, help="Diagnosis rate for long shedders")
    
    # Actions
    parser.add_argument('--verify', action='store_true', help="Run a verification simulation after calculating NSR_base")

    args = parser.parse_args()
    experiment_numbered_name = f'{args.name}_#{args.exp_num}'

    # --- STEP 2: Homogeneous Sweep ---
    settings_func = generate_homogeneous_settings(args.min, args.max, args.steps, args.seeds)
    print("\n[Pipeline] Phase 1: Launching homogeneous calibration sweep...")
    run_experiment_script(args.runner, args.exp_num, settings_func, args.name)
    
    # --- STEP 2b: Fit Curve ---
    fit_result = fit_calibration_curve(experiment_numbered_name, model_type=args.model)
    if not fit_result: return

    # --- STEP 3: Inverse Math (NSR_eff) ---
    print("\n[Pipeline] Phase 2: Solving for NSR_base...")
    NSR_eff = er.compute_calibrated_parameter(args.model, fit_result, args.target_osr)
        
    # --- STEP 4: Time Weighting (NSR_base) ---
    
    # Fetch ih model params
    standard_params = sm.read_standard_parameters_values()
    tau_1 = standard_params['tau_1']
    tau_2 = standard_params['tau_2']
    tau_3 = standard_params['tau_3']
    tau_3_long = standard_params['tau_3_long']
    tau_4 = standard_params['tau_4']

    # Fetch the k_d values using the user's diagnosis rates
    k_ds = dr.get_k_d_from_diagnosis_rate(args.diag_std, tau_1, tau_2, tau_3, tau_4)
    k_dl = dr.get_k_d_from_diagnosis_rate(args.diag_long, tau_1, tau_2, tau_3_long, tau_4)
    
    # Fetch emergence rate
    k_v = standard_params['IH_virus_emergence_rate']

    # Calculate the expected Lineage evolutionary time (W) 
    W_norm = calculate_expected_lineage_days(
        tau_1=tau_1, tau_2=tau_2, tau_3=tau_3, tau_4=tau_4, k_d=k_ds, 
        k_v=k_v, cap_min=1, cap_max=4
    )
    W_long = calculate_expected_lineage_days(
        tau_1=tau_1, tau_2=tau_2, tau_3=tau_3_long, tau_4=tau_4, k_d=k_dl, 
        k_v=k_v, cap_min=5, cap_max=15
    )

    # Compute base NSR 
    F_norm, F_long = calculate_time_fractions(args.r, W_norm, W_long)
    NSR_base = NSR_eff / (F_norm + (args.M * F_long))
      
    print("\n=========================================================")
    print(" CALIBRATION PIPELINE RESULTS")
    print("=========================================================")
    print(f" Target OSR mixed     : {args.target_osr:.6f}")
    print(f" Target Multiplier (M): {args.M:.2f}")
    print(f" Population ratio (r) : {args.r:.4f}")
    print(f" Time Fractions       : F_norm = {F_norm:.3f}, F_long = {F_long:.3f}")
    print("---------------------------------------------------------")
    print(f" NSR required : {NSR_eff:.6e}")
    print(f" Calculated NSR_base    : {NSR_base:.6e}")
    print(f" Calculated NSR_long    : {NSR_base * args.M:.6e}")
    print("=========================================================\n")

    # --- Verify the Math ---
    if args.verify and args.runner != 'slurm':
        if args.runner == 'serial': runner_module = simplicity.runners.serial
        else: runner_module = simplicity.runners.multiprocessing
        verify_mixed_simulation(runner_module, NSR_base, args.M, args.r, args.target_osr, args.seeds)

if __name__ == "__main__":
    main()