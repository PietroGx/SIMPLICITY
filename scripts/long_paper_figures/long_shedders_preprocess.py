#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
from collections import Counter

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
from simplicity.phenotype.distance import hamming_iw

def read_master_log(exp_num):
    """Reads the master grid log to find available parameters."""
    log_path = f"Data/master_grid_log_#{int(exp_num)}.csv"
    
    if not os.path.exists(log_path):
        raise FileNotFoundError(f"Master log not found at {log_path}. Did the experiment run?")
        
    return pd.read_csv(log_path)

def get_baseline_sod(exp_num):
    """Returns the first SOD for the standard calibration baseline."""
    exp_name = f"calibration_run_#{int(exp_num)}"
    sods = dm.get_simulation_output_dirs(exp_name)
    if not sods:
        raise RuntimeError(f"No SOD found for baseline {exp_name}")
    return sods[0]

def get_grid_sod(M, ratio, tau, R, exp_num):
    """Dynamically reconstructs the folder name and returns the SOD."""
    exp_name = f"long_shedders_exp_M{float(M)}_lsr{float(ratio)}_tau{float(tau)}_R{float(R)}_#{int(exp_num)}"
    sods = dm.get_simulation_output_dirs(exp_name)
    if not sods:
        raise RuntimeError(f"No SOD found for {exp_name}")
    return sods[0]

# --- Figure 1 & 2 Extractors ---
def extract_infection_durations(ssod, host_type):
    """Extracts infection durations for a specific host type."""
    df = om.read_individuals_data(ssod)
    df = df.dropna(subset=['t_infectious', 't_not_infectious']).copy()
    sub = df[df["type"] == host_type]
    return (sub['t_not_infectious'] - sub['t_infectious']).to_numpy()

def analyze_long_shedders_trajectories(ssod):
    """Extracts and fixes IH lineage trajectories for individuals."""
    individuals_data = om.read_individuals_data(ssod)
    t_final = om.read_final_time(ssod)
    polished_df = individuals_data.copy()
    polished_df['type'] = polished_df['type'].astype(str)
    
    # Drop infected individuals that did not become infectious yet
    polished_df = polished_df[~((polished_df['state'] == 'infected') & (polished_df['state_t'] < 5))]
    infected_mask = (polished_df['state'] == 'infected') & polished_df['t_not_infected'].isna()
    polished_df.loc[infected_mask, 't_not_infected'] = t_final

    for idx, row in polished_df.iterrows():
        state = row['state']
        t0 = row['t_infection']
        traj = row['IH_lineages_trajectory']
        fixed_traj = {}
        if isinstance(traj, dict):
            for lineage, info in traj.items():
                ih_birth = info.get('ih_birth', t0)
                ih_death = info.get('ih_death')
                if ih_death is None and state == 'infected':
                    ih_death = t_final
                fixed_traj[lineage] = {'ih_birth': ih_birth, 'ih_death': ih_death}
        polished_df.at[idx, 'IH_lineages_trajectory'] = fixed_traj
        
    return polished_df, t_final

# --- Figure 3 Extractors (Hamming Distances) ---
def compute_transmission_hamming_distances(ssod, donor_type_filter=None):
    df = om.read_individuals_data(ssod)
    phylo = om.read_phylogenetic_data(ssod)
    lin2gen = dict(zip(phylo["Lineage_name"], phylo["Genome"]))
    dists = []
    
    for _, row in df.iterrows():
        if donor_type_filter and str(row.get("type", "")) != donor_type_filter:
            continue
        inherited = row.get("inherited_lineage", None)
        if inherited not in lin2gen: continue
        
        for ev in row.get("new_infections", []):
            if not isinstance(ev, dict): continue
            tlin = ev.get("transmitted_lineage", None)
            if tlin in lin2gen:
                d = hamming_iw(lin2gen[inherited], lin2gen[tlin])
                if np.isfinite(d): dists.append(d)
    return dists

def cache_and_aggregate_hamming(exp_num, M, R, ratio, tau_values):
    """Collects genetic distances for baseline and specific tau slices."""
    rows = []
    # 1. Baseline Standard Hosts
    base_sod = get_baseline_sod(exp_num)
    for ssod in dm.get_seeded_simulation_output_dirs(base_sod):
        for d in compute_transmission_hamming_distances(ssod, "standard"):
            rows.append({"group": "Baseline (Standard)", "dist": d})
            
    # 2. Long Shedders across Tau values
    for tau in tau_values:
        try:
            sod = get_grid_sod(M, ratio, tau, R, exp_num)
            for ssod in dm.get_seeded_simulation_output_dirs(sod):
                for d in compute_transmission_hamming_distances(ssod, "long_shedder"):
                    rows.append({"group": f"{float(tau):g}d", "dist": d})
        except Exception as e:
            print(f"[Warning] Could not process tau={tau}: {e}")
            
    return pd.DataFrame(rows)

# --- Figure 4 Extractors (Clade Winners) ---
def _peak_winner(F_clade):
    if F_clade.empty or F_clade.shape[1] == 0: return None
    return F_clade.max(axis=0).idxmax()

def _survival_winner(F_clade, eps=1e-4):
    if F_clade.empty or F_clade.shape[1] == 0: return None
    times = F_clade.index.to_numpy(dtype=float)
    durations = {}
    for clade in F_clade.columns:
        mask = F_clade[clade].to_numpy(dtype=float) > eps
        if not mask.any(): durations[clade] = 0.0
        else: durations[clade] = float(times[len(mask) - 1 - mask[::-1].argmax()] - times[mask.argmax()])
    return max(durations.items(), key=lambda kv: kv[1])[0]

def _fastest_growth_winner(F_clade, lower=0.01, upper=0.50):
    if F_clade.empty or F_clade.shape[1] == 0: return None
    times = F_clade.index.to_numpy(dtype=float)
    candidates = {}
    for clade in F_clade.columns:
        y = F_clade[clade].to_numpy(dtype=float)
        # Simplified crossing logic for speed
        over_lower = times[y >= lower]
        over_upper = times[y >= upper]
        if len(over_lower) > 0 and len(over_upper) > 0 and over_upper[0] > over_lower[0]:
            candidates[clade] = over_upper[0] - over_lower[0]
    return min(candidates.items(), key=lambda kv: kv[1])[0] if candidates else None

def get_clade_winners(ssod, cluster_threshold):
    """Evaluates peak, burden, survival, and growth winners for a single SSOD."""
    phylo_df = om.read_phylogenetic_data(ssod)
    lf = om.read_lineage_frequency(ssod)
    
    clade_to_lineages, _, per_clade_mut_df, clade_meta_df = cl.cluster_lin_into_clades_with_meta(
        phylo_df, shared_mut_threshold=cluster_threshold
    )
    clade_labels_series, _ = cl.label_clades_from_definers(per_clade_mut_df, clade_meta_df)
    clade_labels = clade_labels_series.to_dict() if clade_labels_series is not None else {}
    
    # Build Clade Freq DF
    F_lineage = lf.pivot(index="Time_sampling", columns="Lineage_name", values="Frequency_at_t").sort_index()
    series = []
    for clade, members in clade_to_lineages.items():
        cols = [c for c in members if c in F_lineage.columns]
        if cols: series.append(F_lineage[cols].sum(axis=1).rename(clade))
    F_clade = pd.concat(series, axis=1) if series else pd.DataFrame()

    # Get non-root clades
    non_root = clade_meta_df.loc[(clade_meta_df["parent_clade"].notna()) & (clade_meta_df["n_defining"] > 0), "clade"].tolist()
    F_nr = F_clade.reindex(columns=non_root).dropna(axis=1, how="all") if non_root else pd.DataFrame()

    # Calculate Winners
    peak_winner = _peak_winner(F_nr)
    survival_winner = _survival_winner(F_nr)
    growth_winner = _fastest_growth_winner(F_nr)
    
    # Burden Winner
    lin_to_total = dict(zip(phylo_df["Lineage_name"], phylo_df.get("Total_infections", np.zeros(len(phylo_df)))))
    totals = {clade: sum(int(lin_to_total.get(l, 0)) for l in members) for clade, members in clade_to_lineages.items()}
    burden_winner = max(totals.items(), key=lambda kv: kv[1])[0] if totals else None

    # Resolve labels
    def resolve_label(clade):
        if not clade: return "standard"
        raw = str(clade_labels.get(clade, "standard")).lower()
        if "long" in raw: return "long"
        if "mix" in raw: return "mixed"
        if "founder" in raw: return "founder"
        return "standard"

    return {
        "Peak": resolve_label(peak_winner),
        "Burden": resolve_label(burden_winner),
        "Survival": resolve_label(survival_winner),
        "Growth": resolve_label(growth_winner)
    }

def summarize_sod_pies(sod, cluster_threshold, min_days):
    """Aggregates winners across all seeds in an SOD."""
    results = {"Peak": Counter(), "Burden": Counter(), "Survival": Counter(), "Growth": Counter()}
    valid_seeds = 0
    for ssod in dm.get_seeded_simulation_output_dirs(sod):
        try:
            if om.read_final_time(ssod) >= min_days:
                winners = get_clade_winners(ssod, cluster_threshold)
                for k, v in winners.items(): results[k][v] += 1
                valid_seeds += 1
        except Exception:
            pass
    return results, valid_seeds
    
# --- Figure 5 Extractors (Epidemiological Stats) ---
def extract_epi_stats(sod):
    """Extracts duration, total infections, and total diagnoses across all seeds in an SOD."""
    rows = []
    for ssod in dm.get_seeded_simulation_output_dirs(sod):
        try:
            ft = om.read_final_time(ssod)
            df = om.read_individuals_data(ssod)
            
            std_mask = df['type'] == 'standard'
            lng_mask = df['type'] == 'long_shedder'
            
            # Everyone in the file was infected
            std_inf = std_mask.sum()
            lng_inf = lng_mask.sum()
            
            # Diagnosed individuals end in state == 'diagnosed'
            diag_mask = df['state'] == 'diagnosed'
            std_diag = (std_mask & diag_mask).sum()
            lng_diag = (lng_mask & diag_mask).sum()
                
            rows.append({
                "duration": ft,
                "inf_std": std_inf,
                "inf_long": lng_inf,
                "diag_std": std_diag,
                "diag_long": lng_diag
            })
        except Exception:
            pass
            
    return pd.DataFrame(rows)