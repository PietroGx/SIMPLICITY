#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 12:42:55 2025

@author: pietro
"""

import simplicity.plots_manager as pm
import simplicity.output_manager as om
import simplicity.dir_manager as dm

def analyze_long_shedders(ssod):
    """
    Analyze key metrics for long shedders and update individual data with fixed t_not_infected and ih_death.

    Returns:
        - pd.DataFrame of long shedder metrics:
            - individual_id
            - infection_duration
            - n_IH_lineages_total
            - new_infections_per_lineage (dict)
            - n_IH_lineages_at_end
            - ih_lineage_details (list of dicts)
        - pd.DataFrame of updated individuals_data (with modified t_not_infected and IH_lineages_trajectory)
    """
    import pandas as pd
    from collections import Counter

    # Load original data
    individuals_data = om.read_individuals_data(ssod)
    t_final = om.read_final_time(ssod)

    # Modify a copy 
    polished_df = individuals_data.copy()
    polished_df['type'] = polished_df['type'].astype(str)
    
    # drop infected individuals that did not become infectious yet
    polished_df = polished_df[~((polished_df['state'] == 'infected') & (polished_df['state_t'] < 5))]

    # Apply infection end time fix to infected individuals
    infected_mask = (polished_df['state'] == 'infected') & polished_df['t_not_infected'].isna()
    polished_df.loc[infected_mask, 't_not_infected'] = t_final

    metrics = []

    for idx, row in polished_df.iterrows():
        individual_id = idx
        state = row['state']
        state_t = row['state_t']
        t0 = row['t_infection']
        t1 = row['t_not_infected']
        traj = row['IH_lineages_trajectory']
        new_infs = row['new_infections']

        # Count total IH lineages
        n_lineages = len(traj)

        # Count new infections per lineage
        lineage_counts = Counter()
        for inf in new_infs:
            lineage = inf.get('transmitted_lineage')
            if lineage:
                lineage_counts[lineage] += 1

        # Count lineages at infection end (from column)
        n_at_end = row.get('IH_unique_lineages_number', 0)

        # Adjust lineage trajectory (fix ih_death if missing and still infected)
        ih_lineage_details = []
        fixed_traj = {}
        for lineage, info in traj.items():
            ih_birth = info.get('ih_birth', t0)
            ih_death = info.get('ih_death')

            if ih_death is None and state == 'infected':
                ih_death = t_final

            ih_lineage_details.append({
                'lineage': lineage,
                'ih_birth': ih_birth,
                'ih_death': ih_death
            })

            fixed_traj[lineage] = {'ih_birth': ih_birth, 'ih_death': ih_death}

        # Update the row's trajectory with fixed ih_death
        polished_df.at[individual_id, 'IH_lineages_trajectory'] = fixed_traj
        if row['type'] == 'long_shedder':
            metrics.append({
                'individual_id': individual_id,
                'infection_duration': t1 - t0,
                'n_IH_lineages_total': n_lineages,
                'new_infections_per_lineage': dict(lineage_counts),
                'n_IH_lineages_at_end': n_at_end,
                'ih_lineage_details': ih_lineage_details,
                'is_infected_at_end': (state == 'infected'),
                'is_infectious_at_end': (state_t <19 and state_t >4),
                'time_not_infected': t1
            })


    metrics_df = pd.DataFrame(metrics)
    return metrics_df, polished_df

def prepare_figure_1_data(ssod):
    individuals_data = om.read_individuals_data(ssod)
    
    duration_data = individuals_data.dropna(subset=['t_infectious', 't_not_infectious']).copy()
    duration_data['infection_duration'] = duration_data['t_not_infectious'] - duration_data['t_infectious']

    return duration_data

def prepare_figure_2_data(ssod):

    _, polished_data = analyze_long_shedders(ssod)
    t_final = om.read_final_time(ssod)
    colormap_df = pm.make_lineages_colormap(ssod)
    
    return polished_data, colormap_df, t_final