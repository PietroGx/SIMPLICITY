import pandas as pd
import numpy as np
import random
import os

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
import simplicity.settings_manager as sm
from simplicity.phenotype.distance import hamming_iw

def get_clinical_label(scenario):
    if scenario == "control": return "Control"
    if scenario == "SOT": return "SOT"
    if scenario == "HIV_low": return "HIV Low"
    if scenario == "HIV_high": return "HIV High"
    if scenario == "edge_case": return "Edge Case"
    return scenario.replace("_", " ").title()

def get_target_seed_dir(exp_num, scenario, target_seed, exp_name="impact_long_shedders"):
    """Uses dm.get_ssod to cleanly locate the specific Seeded Simulation Output Directory."""
    exp_str = f"{exp_name}_{scenario}_#{exp_num}"
    sods = dm.get_simulation_output_dirs(exp_str)
    
    if not sods:
        print(f"DEBUG [{scenario}]: No output directories found matching '{exp_str}'")
        return None
        
    try:
        seed_int = int(target_seed)
    except ValueError:
        print(f"DEBUG [{scenario}]: target_seed '{target_seed}' cannot be converted to integer.")
        return None

    for sod in sods:
        try:
            ssod = dm.get_ssod(sod, seed_int)
            return ssod
        except ValueError as e:
            continue
            
    return None

def get_shared_valid_seeds(exp_num, scenarios, exp_name="impact_long_shedders", min_days=365, num_seeds=10):
    """Pre-scans the CONTROL group to find valid seeds, establishing a baseline for all groups."""
    print(f"\n--- Establishing baseline: Selecting {num_seeds} seeds from 'control' (run time >= {min_days} days) ---")
    
    control_str = f"{exp_name}_control_#{exp_num}"
    sods = dm.get_simulation_output_dirs(control_str)
    
    if not sods:
        print("DEBUG [control]: No SODs found. Returning empty list.")
        return []
        
    valid_control_seeds = set()
    for sod in sods:
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            try:
                seed_str = dm.get_seed_from_SSOD(str(ssod))
                seed_int = int(seed_str)
                
                t_final = om.read_final_time(ssod)
                if t_final >= min_days:
                    valid_control_seeds.add(seed_int)
            except Exception as e:
                continue
                
    valid_list = sorted(list(valid_control_seeds))
    print(f"DEBUG [control]: Found {len(valid_list)} valid baseline seeds.")
    
    if len(valid_list) == 0:
        print("DEBUG: Zero valid seeds found in control. Violins will be empty.")
        return []
        
    if len(valid_list) < num_seeds:
        print(f"DEBUG: Warning - Only found {len(valid_list)} valid control seeds, proceeding with available.")
        selected = valid_list
    else:
        selected = random.sample(valid_list, num_seeds)
        
    print(f"DEBUG: Selected target seeds for violins (applied universally): {selected}\n")
    return selected

def get_fig2_freq_data(exp_num, scenario, target_seed):
    """Row 1: Extract Lineage Frequencies natively for a specific seed."""
    ssod = get_target_seed_dir(exp_num, scenario, target_seed)
    if not ssod: return pd.DataFrame(), None, 0
    
    try:
        lf = om.read_lineage_frequency(ssod)
        t_final = om.read_final_time(ssod)
        import simplicity.plots_manager as pm
        cmap_df = pm.make_lineages_colormap(ssod)
        return lf, cmap_df, t_final
    except Exception as e:
        print(f"DEBUG [{scenario}]: Error extracting freq data: {e}")
        return pd.DataFrame(), None, 0

def get_fig2_clustered_data(exp_num, scenario, target_seed, cluster_threshold=5):
    """Row 2: Extract Phylogenetic Tree and perform clustering natively for a specific seed."""
    ssod = get_target_seed_dir(exp_num, scenario, target_seed)
    if not ssod: return pd.DataFrame(), {}, {}, None
    
    try:
        phylo_df = om.read_phylogenetic_data(ssod)
        lf = om.read_lineage_frequency(ssod)
        
        clade_to_lineages, _, _, clade_meta_df = cl.cluster_lin_into_clades_with_meta(
            phylo_df, shared_mut_threshold=cluster_threshold
        )
        
        import simplicity.plots_manager as pm
        cmap_df = pm.make_lineages_colormap(ssod)
        
        return lf, clade_to_lineages, clade_meta_df, cmap_df
    except Exception as e:
        print(f"DEBUG [{scenario}]: Error extracting cluster data: {e}")
        return pd.DataFrame(), {}, {}, None

def get_fig2_divergence_data(exp_num, scenario, target_seeds, exp_name="impact_long_shedders"):
    """Row 3: Aggregated Hamming Distance evaluated ONLY on the shared target_seeds."""
    exp_str = f"{exp_name}_{scenario}_#{exp_num}"
    dists = []
    
    sods = dm.get_simulation_output_dirs(exp_str)
    if not sods: return []
    
    target_type = "standard" if scenario == "control" else "long_shedder"
    processed_count = 0
    
    for sod in sods:
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            try:
                seed_str = dm.get_seed_from_SSOD(str(ssod))
                seed_int = int(seed_str)
            except Exception:
                continue
                
            if seed_int not in target_seeds:
                continue
                
            processed_count += 1
            try:
                df_ind = om.read_individuals_data(ssod)
                phylo = om.read_phylogenetic_data(ssod)
            except Exception as e:
                continue
            
            if df_ind.empty or phylo.empty: continue
            
            lin2gen = dict(zip(phylo["Lineage_name"], phylo["Genome"]))
            
            for _, row in df_ind.iterrows():
                if str(row.get("type", "")) != target_type: continue
                inherited = row.get("inherited_lineage", None)
                if inherited not in lin2gen: continue
                
                new_infections = row.get("new_infections", [])
                if not isinstance(new_infections, list): continue
                
                for ev in new_infections:
                    if not isinstance(ev, dict): continue
                    tlin = ev.get("transmitted_lineage", None)
                    if tlin in lin2gen:
                        d = hamming_iw(lin2gen[inherited], lin2gen[tlin])
                        if np.isfinite(d): dists.append(d)
                        
    print(f"[{scenario}] Processed {processed_count}/{len(target_seeds)} baseline seeds: {len(dists)} jumps found.")
    return dists
