import argparse
import pandas as pd
import numpy as np
import os
import random
import ast

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm
import simplicity.output_manager as om
from simplicity.intra_host_model import Host 

def parse_arguments():
    parser = argparse.ArgumentParser(description="SIMPLICITY Figure 1 Data Preprocessor")
    parser.add_argument('--engine', type=str, choices=['serial', 'multiprocessing', 'slurm'], default='serial')
    parser.add_argument('--force-recompute', action='store_true')
    return parser.parse_args()

def get_clinical_label(scenario):
    if scenario == "control": return "Control"
    if scenario == "SOT": return "SOT"
    if "HIV" in scenario: return "HIV"
    return scenario.replace("_", " ").title()

def get_panel_a_data(exp_num=1, exp_name="impact_long_shedders"):
    scenarios = ["control", "SOT", "HIV_low"]
    df_list = []
    for scenario in scenarios:
        experiment_string = f"{exp_name}_{scenario}_#{exp_num}"
        try:
            sods = dm.get_simulation_output_dirs(experiment_string)
            if not sods: continue
            sod = sods[0]
            tau_1 = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_1')
            tau_2 = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_2')
            tau_4 = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_4')
            tau_3_standard = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3')
            if scenario == "control": tau_3_active = tau_3_standard
            else: tau_3_active = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')
            label = get_clinical_label(scenario)
            host = Host(tau_1=tau_1, tau_2=tau_2, tau_3=tau_3_active, tau_4=tau_4)
            time_max = 800
            step = 1
            p_inf, p_det, p_rec = host.data_plot_ih_solution(state=0, time=time_max, step=step)
            time_array = np.arange(0, time_max, step)
            temp_df = pd.DataFrame({'time': time_array[:len(p_inf)], 'p_infectious': p_inf, 'cohort': label})
            df_list.append(temp_df)
        except Exception: continue
    if not df_list: return pd.DataFrame(columns=['time', 'p_infectious', 'cohort'])
    return pd.concat(df_list, ignore_index=True)

def get_panel_b_data(exp_num=1, exp_name="impact_long_shedders"):
    scenarios = ["control", "SOT", "HIV_low", "HIV_high"]
    df_list = []
    control_exp_string = f"{exp_name}_control_#{exp_num}"
    try:
        control_sods = dm.get_simulation_output_dirs(control_exp_string)
        if not control_sods: return pd.DataFrame(columns=['duration', 'cohort'])
        control_sod = control_sods[0]
        all_control_ssods = dm.get_seeded_simulation_output_dirs(control_sod)
        rng = random.Random(42) 
        sampled_control_ssods = rng.sample(all_control_ssods, min(10, len(all_control_ssods)))
        matched_seed_names = [os.path.basename(ssod) for ssod in sampled_control_ssods]
    except Exception: return pd.DataFrame(columns=['duration', 'cohort'])

    for scenario in scenarios:
        experiment_string = f"{exp_name}_{scenario}_#{exp_num}"
        try:
            sods = dm.get_simulation_output_dirs(experiment_string)
            if not sods: continue
            sod = sods[0]
            all_ssods_for_scenario = dm.get_seeded_simulation_output_dirs(sod)
            for target_seed_name in matched_seed_names:
                target_ssod = next((s for s in all_ssods_for_scenario if os.path.basename(s) == target_seed_name), None)
                if target_ssod is None: continue
                df = om.read_individuals_data(target_ssod)
                if scenario == "control": filtered_df = df[df['type'] == 'standard'].copy()
                else: filtered_df = df[df['type'] == 'long_shedder'].copy()
                label = get_clinical_label(scenario)
                if 't_not_infected' in filtered_df.columns and 't_infection' in filtered_df.columns:
                    filtered_df['duration'] = filtered_df['t_not_infected'] - filtered_df['t_infection']
                    filtered_df = filtered_df[['duration']].dropna()
                    filtered_df['cohort'] = label
                    df_list.append(filtered_df)
        except Exception: continue
    if not df_list: return pd.DataFrame(columns=['duration', 'cohort'])
    return pd.concat(df_list, ignore_index=True)

def get_panel_c_data():
    if not os.path.exists('literature_long_shedders_data.csv'):
        return pd.DataFrame(columns=['duration', 'clinical_category'])
    df = pd.read_csv('literature_long_shedders_data.csv')
    if 'Days (Exact)' in df.columns:
        df['duration'] = df['Days (Exact)']
        df = df.dropna(subset=['clinical_category', 'duration']).copy()
        return df[['duration', 'clinical_category']]
    return pd.DataFrame(columns=['duration', 'clinical_category'])

def get_panel_de_data():
    df_d = pd.DataFrame()
    df_e = pd.DataFrame()
    if os.path.exists('data_fig1_D.csv'):
        df_d = pd.read_csv('data_fig1_D.csv')
        if 'sampling_date' in df_d.columns: df_d['sampling_date'] = pd.to_datetime(df_d['sampling_date'])
    if os.path.exists('data_fig1_E.csv'):
        df_e = pd.read_csv('data_fig1_E.csv')
        if 'sampling_date' in df_e.columns: df_e['sampling_date'] = pd.to_datetime(df_e['sampling_date'])
    return df_d, df_e

def get_panel_fg_data(exp_num=1, exp_name="impact_long_shedders"):
    print(f"\n--- EXTRACTING PANELS F & G ---")
    df_f_list = []
    df_g_list = []
    
    # Extract Panel F (Control Standard)
    control_exp = f"{exp_name}_control_#{exp_num}"
    try:
        sods = dm.get_simulation_output_dirs(control_exp)
        if sods:
            for ssod in dm.get_seeded_simulation_output_dirs(sods[0]):
                seq_path = os.path.join(ssod, 'sequencing_data_regression.csv')
                if not os.path.exists(seq_path): continue
                df = pd.read_csv(seq_path)
                col = 'individual_type' if 'individual_type' in df.columns else ('Individual_type' if 'Individual_type' in df.columns else None)
                if col:
                    filtered_df = df[df[col] == 'standard'].copy()
                    if not filtered_df.empty:
                        filtered_df['cohort'] = 'Control'
                        df_f_list.append(filtered_df)
    except Exception: pass

    # Extract Panel G (Long Shedders)
    scenarios = ["SOT", "HIV_low", "HIV_high"]
    L = 4967 

    for scenario in scenarios:
        exp_str = f"{exp_name}_{scenario}_#{exp_num}"
        try:
            sods = dm.get_simulation_output_dirs(exp_str)
            if not sods: continue
            for ssod in dm.get_seeded_simulation_output_dirs(sods[0]):
                seq_path = os.path.join(ssod, 'sequencing_data.csv')
                ind_path = os.path.join(ssod, 'individuals_data.csv')
                phylo_path = os.path.join(ssod, 'phylogenetic_data.csv')
                
                if not (os.path.exists(seq_path) and os.path.exists(ind_path) and os.path.exists(phylo_path)): 
                    continue
                    
                df_seq = pd.read_csv(seq_path)
                df_ind = pd.read_csv(ind_path)
                df_phylo = pd.read_csv(phylo_path)
                
                df_phylo['genome_length'] = df_phylo['Genome'].apply(lambda x: len(ast.literal_eval(x)) if isinstance(x, str) else 0)
                phylo_dict = dict(zip(df_phylo['Lineage_name'], df_phylo['genome_length']))
                phylo_dict['wt'] = 0 
                phylo_dict['root'] = 0
                
                col_id = 'Unnamed: 0' if 'Unnamed: 0' in df_ind.columns else df_ind.columns[0]
                df_ind_sub = df_ind[[col_id, 'type', 'inherited_lineage']].copy()
                df_ind_sub.rename(columns={col_id: 'individual_index'}, inplace=True)
                
                df_merged = df_seq.merge(df_ind_sub, on='individual_index', how='inner')
                
                if not df_merged.empty:
                    print(f"\n[DIAGNOSTICS - {scenario} - {os.path.basename(ssod)}]")
                    
                    # 1. Map inherited lengths and check for failures
                    df_merged['inherited_length'] = df_merged['inherited_lineage'].map(phylo_dict)
                    
                    unmapped = df_merged[df_merged['inherited_length'].isna()]
                    if not unmapped.empty:
                        print(f"  ⚠️ WARNING: {len(unmapped)} sequences failed to map their inherited_lineage!")
                        print(f"     Example unmapped lineages: {unmapped['inherited_lineage'].unique()[:5]}")
                        
                    # Fillna with 0 temporarily just so we can see the resulting impossible math in the next check
                    df_merged['inherited_length_filled'] = df_merged['inherited_length'].fillna(0)
                    df_merged['intra_host_mutations'] = df_merged['sequence_lenght'] - df_merged['inherited_length_filled']
                    
                    # 2. Check for negative mutations
                    neg_muts = df_merged[df_merged['intra_host_mutations'] < 0]
                    if not neg_muts.empty:
                        print(f"  ⚠️ WARNING: {len(neg_muts)} sequences have NEGATIVE intra-host mutations!")
                        print(f"     Example (Seq Length vs Inherited Length):")
                        print(neg_muts[['sequence_lenght', 'inherited_length_filled', 'intra_host_mutations']].head(3))
                        
                    # 3. Check for suspiciously high mutations very early in infection (e.g. >5 mutations in <10 days)
                    impossible = df_merged[(df_merged['infection_duration'] < 10) & (df_merged['intra_host_mutations'] > 5)]
                    if not impossible.empty:
                        print(f"  ⚠️ WARNING: {len(impossible)} sequences have suspiciously high mutation counts early in infection!")
                        print(f"     Example (Duration vs Mutations):")
                        print(impossible[['infection_duration', 'sequence_lenght', 'inherited_length_filled', 'intra_host_mutations']].head(3))

                    # Proceed with extraction, keeping the raw values so we can plot them as is
                    df_merged['Distance_from_root'] = df_merged['intra_host_mutations'] / L
                    df_merged['Sequencing_time'] = df_merged['infection_duration'] / 365.25
                    df_merged['cohort'] = get_clinical_label(scenario)
                    
                    df_g_list.append(df_merged[['Sequencing_time', 'Distance_from_root', 'cohort']])
        except Exception as e: 
            print(f"❌ Error in {scenario}: {e}")

    df_f = pd.concat(df_f_list, ignore_index=True) if df_f_list else pd.DataFrame(columns=['Sequencing_time', 'Distance_from_root', 'cohort'])
    df_g = pd.concat(df_g_list, ignore_index=True) if df_g_list else pd.DataFrame(columns=['Sequencing_time', 'Distance_from_root', 'cohort'])
    
    return df_f, df_g

def main():
    args = parse_arguments()
    os.makedirs('plot_data_cache', exist_ok=True)
if __name__ == "__main__": main()
