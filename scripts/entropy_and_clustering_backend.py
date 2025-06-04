# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from collections import defaultdict, deque
from scipy.stats import mannwhitneyu

import simplicity.output_manager as om
import simplicity.dir_manager as dm
import simplicity.plots_manager as pm


pm.apply_plos_rcparams()

# -----------------------------------------------------------------------------
# Utility Functions (no globals)
# -----------------------------------------------------------------------------

def parse_genome(genome):
    """Convert a genome representation to a frozenset of mutation tuples."""
    if genome is None:
        return frozenset()
    try:
        return frozenset(tuple(x) for x in genome)
    except Exception:
        return frozenset()

def filter_seeds(sim_out_dir, min_time, max_time, cut_after_max):
    """
    Return a sorted list of seed indices whose final_time ≥ min_time,
    and (if cut_after_max is False) final_time ≤ max_time.
    If cut_after_max is True, accept final_time > max_time but trim data later.
    """
    filtered = []
    for ssod in dm.get_seeded_simulation_output_dirs(sim_out_dir):
        try:
            final_time = om.read_final_time(ssod)
            if final_time < min_time:
                continue
            if not cut_after_max and final_time > max_time:
                continue
            seed_idx = int(os.path.basename(ssod).split("_")[-1])
            filtered.append(seed_idx)
        except Exception:
            continue
    return sorted(filtered)

def load_seed_data(sim_out_dir, seed, max_time, cut_after_max):
    """
    Load lineage-to-mutations dict and the raw lineage-frequency DataFrame for a given seed.
    If cut_after_max is True, trim freq_df to times ≤ max_time.
    """
    ssod = dm.get_ssod(sim_out_dir, seed)
    phylo = om.read_phylogenetic_data(ssod)
    lineage_to_mutations = dict(zip(
        phylo["Lineage_name"],
        phylo["Genome"].apply(parse_genome)
    ))
    freq_raw = om.read_lineage_frequency(ssod)
    freq_df_lineages = freq_raw.pivot(
        index="Time_sampling",
        columns="Lineage_name",
        values="Frequency_at_t"
    )
    if cut_after_max:
        freq_df_lineages = freq_df_lineages[freq_df_lineages.index <= max_time]
    return lineage_to_mutations, freq_df_lineages

def cluster_lineages_by_shared_mutations(lin2mut: dict, freq_pivot: pd.DataFrame, min_shared: int):
    """
    Cluster lineages based on shared mutations using transitive closure.
    Returns:
      - freq_df_clusters: DataFrame where each column is a cluster's summed frequency
      - clusters: list of sets (each set = lineage names in that cluster)
      - parents: list of representative lineage names (one per cluster)
    """
    adjacency = defaultdict(set)
    items = list(lin2mut.items())
    for i, (l1, m1) in enumerate(items):
        for j, (l2, m2) in enumerate(items):
            if l1 != l2 and len(m1 & m2) >= min_shared:
                adjacency[l1].add(l2)
                adjacency[l2].add(l1)

    visited, clusters, parents = set(), [], []
    for lineage in lin2mut:
        if lineage not in visited:
            cluster = set()
            queue = deque([lineage])
            parents.append(lineage)
            while queue:
                current = queue.popleft()
                if current not in visited:
                    visited.add(current)
                    cluster.add(current)
                    queue.extend(adjacency[current] - visited)
            clusters.append(cluster)

    cluster_freqs = []
    for i, cluster in enumerate(clusters):
        cols = [col for col in cluster if col in freq_pivot.columns]
        if cols:
            summed = freq_pivot[cols].sum(axis=1)
            cluster_freqs.append(summed.rename(f"Cluster_{i}"))

    freq_df_clusters = pd.concat(cluster_freqs, axis=1) if cluster_freqs else pd.DataFrame(index=freq_pivot.index)
    return freq_df_clusters, clusters, parents

def build_clustered_freqs(lineage_to_mutations, freq_df_lineages, threshold):
    """
    Cluster lineages by shared mutations with threshold, returning:
      - freq_df_clusters: DataFrame of cluster-summed frequencies
      - cluster_parents: list of parent lineage names (one per cluster column)
    """
    if threshold == 0:
        return freq_df_lineages.copy(), list(freq_df_lineages.columns)
    freq_df_clusters, clusters, parents = cluster_lineages_by_shared_mutations(
        lineage_to_mutations, freq_df_lineages, threshold
    )
    return freq_df_clusters, parents

def compute_shannon_entropy(freq_df):
    """
    Compute Shannon entropy across columns for each row in freq_df.
    Assumes each row of freq_df is already normalized (sums to 1).
    """
    normed = freq_df.div(freq_df.sum(axis=1), axis=0)
    log_vals = np.log(normed.where(normed > 0))
    entropy = -(normed * log_vals).sum(axis=1, skipna=True)
    entropy[normed.isna().all(axis=1)] = np.nan
    return entropy

def bin_entropy_over_days(entropy_series, max_time):
    """
    Bin an entropy series (indexed by irregular times) into integer days [0, max_time].
    For each day d, take the mean of entropy_series[t] for t in [d, d+1).
    """
    days = np.arange(0, int(max_time) + 1)
    binned = pd.Series(index=days, dtype=float)
    for d in days:
        in_bin = entropy_series[(entropy_series.index >= d) & (entropy_series.index < d + 1)]
        binned[d] = in_bin.mean() if not in_bin.empty else np.nan
    return binned

def run_mwu(x, y, metric_name):
    stat, p = mannwhitneyu(x, y, alternative='two-sided')
    median_diff = np.median(y) - np.median(x)
    direction = "Group2 > Group1" if median_diff > 0 else "Group1 ≥ Group2"
    return {
        "Metric": metric_name,
        "p_value": p,
        "U_statistic": stat,
        "median_diff": median_diff,
        "direction": direction
    }

def build_seed_cache(sim_dir, filtered_seeds, max_time, cut_after_max, cluster_threshold):
    """
    For each seed in filtered_seeds, load and compute once:
      - freq_l: raw lineage-frequency DataFrame
      - freq_c: clustered-frequency DataFrame
      - ent_l: Shannon entropy Series on raw lineages
      - ent_c: Shannon entropy Series on clusters
      - parents: list of cluster‐parent lineage names
    Returns a dict: { seed: { "freq_l":…, "freq_c":…, "ent_l":…, "ent_c":…, "parents":… } }
    """
    cache = {}
    for seed in tqdm(filtered_seeds, desc=f"Building cache for sim_dir={sim_dir}"):
        # 1. Load raw lineage data
        l2m, freq_l = load_seed_data(sim_dir, seed, max_time, cut_after_max)
        # 2. Cluster lineages
        freq_c, parents = build_clustered_freqs(l2m, freq_l, cluster_threshold)
        # 3. Compute entropies
        ent_l = compute_shannon_entropy(freq_l)
        ent_c = compute_shannon_entropy(freq_c)
        cache[seed] = {
            "freq_l": freq_l,
            "freq_c": freq_c,
            "ent_l" : ent_l,
            "ent_c" : ent_c,
            "parents": parents
        }
    return cache

# -----------------------------------------------------------------------------
# Refactored Data Preparation Functions (no reliance on globals)
# -----------------------------------------------------------------------------

def prepare_figure1_data(sim_out_dir_1, sim_out_dir_2,
                         seed,    # example SEED
                         cache_1, cache_2):
    """
    Load only the data needed for Figure 1, pulling from per-seed caches:
      - For 'seed' in each group: raw lineage frequencies, clustered frequencies,
        their Shannon entropies, and the colormap / cluster parents.
    """
    left_dir, right_dir = sim_out_dir_1, sim_out_dir_2

    # Group 1 (linear phenotype model), using 'seed'
    entry1    = cache_1[seed]
    freq_l1   = entry1["freq_l"]
    freq_c1   = entry1["freq_c"]
    ent_l1    = entry1["ent_l"]
    ent_c1    = entry1["ent_c"]
    parents1  = entry1["parents"]
    colormap1 = pm.make_lineages_colormap(dm.get_ssod(left_dir, seed))

    # Group 2 (immune waning phenotype model), using 'seed'
    entry2    = cache_2[seed]
    freq_l2   = entry2["freq_l"]
    freq_c2   = entry2["freq_c"]
    ent_l2    = entry2["ent_l"]
    ent_c2    = entry2["ent_c"]
    parents2  = entry2["parents"]
    colormap2 = pm.make_lineages_colormap(dm.get_ssod(right_dir, seed))

    return {
        # Group 1
        "freq_lineages_ex1":    freq_l1,
        "entropy_lineages_ex1": ent_l1,
        "colormap_ex1":         colormap1,
        "freq_clusters_ex1":    freq_c1,
        "entropy_clusters_ex1": ent_c1,
        "cluster_parents_ex1":  parents1,

        # Group 2
        "freq_lineages_ex2":    freq_l2,
        "entropy_lineages_ex2": ent_l2,
        "colormap_ex2":         colormap2,
        "freq_clusters_ex2":    freq_c2,
        "entropy_clusters_ex2": ent_c2,
        "cluster_parents_ex2":  parents2,
    }

def prepare_figure2_data(sim_out_dir_1, sim_out_dir_2,
                    filtered_seeds_1, filtered_seeds_2,
                    cache_1, cache_2,
                    takeover_threshold):
    """

    """
    min_days_above_threshold = 21
    def count_sweeps(series_bool, min_days_above_threshold):
        arr = series_bool.values.astype(int)
        total = 0
        current_run_length = 0
        for val in arr:
            if val == 1:
                current_run_length += 1
            else:
                if current_run_length >= min_days_above_threshold:
                    total += 1
                current_run_length = 0
        if current_run_length >= min_days_above_threshold:
            total += 1
        return total

    tcl1, tcc1 = [], []
    tcl2, tcc2 = [], []

    # Group 1
    for seed in tqdm(filtered_seeds_1, desc="Fig2 ‐ Group 1 seeds"):
        entry = cache_1[seed]
        freq_l = entry["freq_l"]
        freq_c = entry["freq_c"]

        # Normalize each day to get p_{t,i}
        norm_l = freq_l.div(freq_l.sum(axis=1), axis=0)
        mask_l = (norm_l >= takeover_threshold)
        # Count sustained runs for lineages:
        sweeps_per_lineage = mask_l.apply(count_sweeps, args=(min_days_above_threshold,), axis=0)
        total_sweeps_lineages = int(sweeps_per_lineage.sum())
        tcl1.append(total_sweeps_lineages)

        # Now do the same for clustered‐lineages:
        norm_c = freq_c.div(freq_c.sum(axis=1), axis=0)
        mask_c = (norm_c >= takeover_threshold)
        sweeps_per_cluster = mask_c.apply(count_sweeps, args=(min_days_above_threshold,), axis=0)
        total_sweeps_clusters = int(sweeps_per_cluster.sum())
        tcc1.append(total_sweeps_clusters)

    # Group 2
    for seed in tqdm(filtered_seeds_2, desc="Fig2 ‐ Group 2 seeds"):
        entry = cache_2[seed]
        freq_l = entry["freq_l"]
        freq_c = entry["freq_c"]

        norm_l = freq_l.div(freq_l.sum(axis=1), axis=0)
        mask_l = (norm_l >= takeover_threshold)
        sweeps_per_lineage = mask_l.apply(count_sweeps, args=(min_days_above_threshold,), axis=0)
        total_sweeps_lineages = int(sweeps_per_lineage.sum())
        tcl2.append(total_sweeps_lineages)

        norm_c = freq_c.div(freq_c.sum(axis=1), axis=0)
        mask_c = (norm_c >= takeover_threshold)
        sweeps_per_cluster = mask_c.apply(count_sweeps, args=(min_days_above_threshold,), axis=0)
        total_sweeps_clusters = int(sweeps_per_cluster.sum())
        tcc2.append(total_sweeps_clusters)

    # Mann–Whitney U tests on the new “total_sweeps” lists:
    res_l = run_mwu(tcl1, tcl2, f"Selective sweeps (lineages, ≥{min_days_above_threshold} days)")
    res_c = run_mwu(tcc1, tcc2, f"Selective sweeps (clusters, ≥{min_days_above_threshold} days)")
    summary_df = pd.DataFrame([res_l, res_c])
    print(summary_df)
    return {
        "selective_sweeps_lineages": {1: tcl1, 2: tcl2},
        "selective_sweeps_clusters": {1: tcc1, 2: tcc2},
        "summary_df": summary_df,
    }

def prepare_figure3_data(sim_out_dir_1, sim_out_dir_2,
                         filtered_seeds_1, filtered_seeds_2,
                         cache_1, cache_2,
                         max_time):
    """
    Load the data needed for Figure 3, pulling from seed caches:
      - Per‐seed, binned entropy‐over‐days for raw lineages and clusters
      - Compute group‐mean and group‐std over those binned series
    """
    binned_ent_l1, binned_ent_c1 = {}, {}
    binned_ent_l2, binned_ent_c2 = {}, {}

    # Group 1: bin raw and cluster entropies over days
    for seed in tqdm(filtered_seeds_1, desc="Fig3 ‐ Group 1 seeds"):
        entry = cache_1[seed]
        ent_l = entry["ent_l"]
        ent_c = entry["ent_c"]
        binned_ent_l1[seed] = bin_entropy_over_days(ent_l, max_time)
        binned_ent_c1[seed] = bin_entropy_over_days(ent_c, max_time)

    # Group 2: bin raw and cluster entropies over days
    for seed in tqdm(filtered_seeds_2, desc="Fig3 ‐ Group 2 seeds"):
        entry = cache_2[seed]
        ent_l = entry["ent_l"]
        ent_c = entry["ent_c"]
        binned_ent_l2[seed] = bin_entropy_over_days(ent_l, max_time)
        binned_ent_c2[seed] = bin_entropy_over_days(ent_c, max_time)

    days = np.arange(0, max_time + 1)

    # Group 1 raw
    df_l1  = pd.DataFrame({seed: binned_ent_l1[seed] for seed in filtered_seeds_1}, index=days)
    mean_l1 = df_l1.mean(axis=1, skipna=True)
    std_l1  = df_l1.std(axis=1, skipna=True)

    # Group 2 raw
    df_l2  = pd.DataFrame({seed: binned_ent_l2[seed] for seed in filtered_seeds_2}, index=days)
    mean_l2 = df_l2.mean(axis=1, skipna=True)
    std_l2  = df_l2.std(axis=1, skipna=True)

    # Group 1 clusters
    df_c1  = pd.DataFrame({seed: binned_ent_c1[seed] for seed in filtered_seeds_1}, index=days)
    mean_c1 = df_c1.mean(axis=1, skipna=True)
    std_c1  = df_c1.std(axis=1, skipna=True)

    # Group 2 clusters
    df_c2  = pd.DataFrame({seed: binned_ent_c2[seed] for seed in filtered_seeds_2}, index=days)
    mean_c2 = df_c2.mean(axis=1, skipna=True)
    std_c2  = df_c2.std(axis=1, skipna=True)

    return {
        "days": days,
        "mean_entropy_lineages": {1: mean_l1, 2: mean_l2},
        "std_entropy_lineages":  {1: std_l1,  2: std_l2},
        "mean_entropy_clusters": {1: mean_c1, 2: mean_c2},
        "std_entropy_clusters":  {1: std_c1,  2: std_c2},
    }

def prepare_figure4_data(sim_out_dir_1, sim_out_dir_2,
                         filtered_seeds_1, filtered_seeds_2,
                         cache_1, cache_2):
    """
    Load the data needed for Figure 4, pulling from seed caches:
      - Per‐seed, “max raw entropy” for raw lineages and clusters
      - Mann–Whitney U tests for those two metrics
    """
    max_ent_l1, max_ent_c1 = [], []
    max_ent_l2, max_ent_c2 = [], []

    # Group 1: max raw entropy per seed
    for seed in tqdm(filtered_seeds_1, desc="Fig4 ‐ Group 1 seeds"):
        entry = cache_1[seed]
        ent_l = entry["ent_l"]
        ent_c = entry["ent_c"]
        max_ent_l1.append(ent_l.max())
        max_ent_c1.append(ent_c.max())

    # Group 2: max raw entropy per seed
    for seed in tqdm(filtered_seeds_2, desc="Fig4 ‐ Group 2 seeds"):
        entry = cache_2[seed]
        ent_l = entry["ent_l"]
        ent_c = entry["ent_c"]
        max_ent_l2.append(ent_l.max())
        max_ent_c2.append(ent_c.max())

    # Mann–Whitney U tests
    res_max_l = run_mwu(max_ent_l1, max_ent_l2, "Max entropy (lineages)")
    res_max_c = run_mwu(max_ent_c1, max_ent_c2, "Max entropy (clusters)")
    summary_df = pd.DataFrame([res_max_l, res_max_c])

    return {
        "max_entropy_lineages":  {1: max_ent_l1, 2: max_ent_l2},
        "max_entropy_clusters":  {1: max_ent_c1, 2: max_ent_c2},
        "summary_df":              summary_df,
    }

def prepare_figure5_data(sim_out_dir_1, sim_out_dir_2,
                         filtered_seeds_1, filtered_seeds_2,
                         cache_1, cache_2):
    """
    Load the data needed for Figure 5, pulling from seed caches:
      - Per‐seed, “last observed entropy” = last non‐NaN raw entropy value
        for raw lineages and clusters
      - Mann–Whitney U tests for those two metrics
    """
    last_ent_l1, last_ent_c1 = [], []
    last_ent_l2, last_ent_c2 = [], []

    # Group 1: last observed raw entropy per seed
    for seed in tqdm(filtered_seeds_1, desc="Fig5 ‐ Group 1 seeds"):
        entry = cache_1[seed]
        ent_l = entry["ent_l"]
        ent_c = entry["ent_c"]
        valid_l = ent_l.dropna()
        last_ent_l1.append(valid_l.iloc[-1] if not valid_l.empty else np.nan)
        valid_c = ent_c.dropna()
        last_ent_c1.append(valid_c.iloc[-1] if not valid_c.empty else np.nan)

    # Group 2: last observed raw entropy per seed
    for seed in tqdm(filtered_seeds_2, desc="Fig5 ‐ Group 2 seeds"):
        entry = cache_2[seed]
        ent_l = entry["ent_l"]
        ent_c = entry["ent_c"]
        valid_l = ent_l.dropna()
        last_ent_l2.append(valid_l.iloc[-1] if not valid_l.empty else np.nan)
        valid_c = ent_c.dropna()
        last_ent_c2.append(valid_c.iloc[-1] if not valid_c.empty else np.nan)

    # Mann–Whitney U tests
    res_last_l = run_mwu(last_ent_l1, last_ent_l2, "Last entropy (lineages)")
    res_last_c = run_mwu(last_ent_c1, last_ent_c2, "Last entropy (clusters)")
    summary_df = pd.DataFrame([res_last_l, res_last_c])

    return {
        "last_entropy_lineages":  {1: last_ent_l1, 2: last_ent_l2},
        "last_entropy_clusters":  {1: last_ent_c1, 2: last_ent_c2},
        "summary_df":               summary_df,
    }

# -----------------------------------------------------------------------------
# Plotting Functions 
# -----------------------------------------------------------------------------

def plot_figure1(axs, data):
    """
    Plot 2x2:
    Row 1: raw lineages (example seed) for both groups
    Row 2: clustered lineages (example seed)
    """
    freq_l1 = data["freq_lineages_ex1"]
    ent_l1  = data["entropy_lineages_ex1"]
    cmap1   = data["colormap_ex1"]

    freq_l2 = data["freq_lineages_ex2"]
    ent_l2  = data["entropy_lineages_ex2"]
    cmap2   = data["colormap_ex2"]

    freq_c1 = data["freq_clusters_ex1"]
    ent_c1  = data["entropy_clusters_ex1"]
    parents1 = data["cluster_parents_ex1"]

    freq_c2 = data["freq_clusters_ex2"]
    ent_c2  = data["entropy_clusters_ex2"]
    parents2 = data["cluster_parents_ex2"]

    # Row 1: Lineages
    for col, (freq_df, entropy_series, colormap, title) in enumerate([
        (freq_l1.fillna(0), ent_l1, cmap1, "linear phenotype model (lineages)"),
        (freq_l2.fillna(0), ent_l2, cmap2, "immune waning phenotype model (lineages)")
    ]):
        ax = axs[0, col]
        colors = [pm.get_lineage_color(lin, colormap) for lin in freq_df.columns]
        ax.stackplot(freq_df.index, freq_df.T.values, colors=colors, alpha=0.6,
                     edgecolors='black', linewidth=0.2)
        ax.set_title(title)
        ax.set_ylabel("Frequency")
        ax2 = ax.twinx()
        ax2.plot(entropy_series.index, entropy_series.values, 'k-', linewidth=1)
        ax2.set_ylabel("Entropy")
        ax.set_xlabel("Time")

    # Row 2: Clustered lineages
    for col, (freq_df, entropy_series, colormap, parents, title) in enumerate([
        (freq_c1, ent_c1, cmap1, parents1, "linear phenotype model (clustered lineages)"),
        (freq_c2, ent_c2, cmap2, parents2, "immune waning phenotype model (clustered lineages)")
    ]):
        ax = axs[1, col]
        colors_c = [pm.get_lineage_color(parent, colormap) for parent in parents]
        ax.stackplot(freq_df.index, freq_df.T.values, colors=colors_c, alpha=0.6,
                     edgecolors='black', linewidth=0.2)
        ax.set_title(title)
        ax.set_ylabel("Frequency")
        ax2 = ax.twinx()
        ax2.plot(entropy_series.index, entropy_series.values, 'k-', linewidth=1)
        ax2.set_ylabel("Entropy")
        ax.set_xlabel("Time")
    
    pm.apply_standard_axis_style(ax,True)

def plot_figure2(axs, data):
    """
    Plot 2x1:
    Row 1: violin of takeover_counts_lineages
    Row 2: violin of takeover_counts_clusters
    """
    tcl1 = data["selective_sweeps_lineages"][1]
    tcl2 = data["selective_sweeps_lineages"][2]
    tcc1 = data["selective_sweeps_clusters"][1]
    tcc2 = data["selective_sweeps_clusters"][2]

    # Row 1: Lineages
    df_tl = pd.DataFrame({
        "Group": ["linear phenotype model"] * len(tcl1) +
                 ["immune waning phenotype model"] * len(tcl2),
        "selective_sweeps_count": tcl1 + tcl2
    })
    sns.violinplot(
        ax=axs[0],
        data=df_tl, x="Group", y="selective_sweeps_count",
        inner="box", cut=0,hue="Group"
    )
    axs[0].set_ylabel("Number of takeover events")

    # Row 2: Clusters
    df_tc = pd.DataFrame({
        "Group": ["linear phenotype model"] * len(tcc1) +
                 ["immune waning phenotype model"] * len(tcc2),
        "selective_sweeps_count": tcc1 + tcc2
    })
    sns.violinplot(
        ax=axs[1],
        data=df_tc, x="Group", y="selective_sweeps_count", inner="box", cut=0,hue="Group"
    )
    axs[1].set_ylabel("Number of selective sweeps")
    
    pm.apply_standard_axis_style(axs[0])
    pm.apply_standard_axis_style(axs[1])
    
def plot_figure3(axs, data):
    """
    Plot 2×2 with ±1 std shading around the mean entropy.
    Row 1: mean entropy of lineages (± std) for both groups
    Row 2: mean entropy of clusters (± std) for both groups
    """
    days  = data["days"]
    mel1  = data["mean_entropy_lineages"][1]
    stdl1 = data["std_entropy_lineages"][1]
    mel2  = data["mean_entropy_lineages"][2]
    stdl2 = data["std_entropy_lineages"][2]

    mec1  = data["mean_entropy_clusters"][1]
    stdc1 = data["std_entropy_clusters"][1]
    mec2  = data["mean_entropy_clusters"][2]
    stdc2 = data["std_entropy_clusters"][2]

    # Row 1: Lineages
    ax = axs[0, 0]
    ax.plot(days, mel1, color="black", linewidth=1)
    ax.fill_between(days, mel1 - stdl1, mel1 + stdl1, color="gray", alpha=0.3)
    ax.set_ylabel("Mean entropy")
    ax.set_xlabel("Day")
    pm.apply_standard_axis_style(ax)

    ax = axs[0, 1]
    ax.plot(days, mel2, color="black", linewidth=1)
    ax.fill_between(days, mel2 - stdl2, mel2 + stdl2, color="gray", alpha=0.3)
    ax.set_ylabel("Mean entropy")
    ax.set_xlabel("Day")
    pm.apply_standard_axis_style(ax)

    # Row 2: Clustered lineages
    ax = axs[1, 0]
    ax.plot(days, mec1, color="black", linewidth=1)
    ax.fill_between(days, mec1 - stdc1, mec1 + stdc1, color="gray", alpha=0.3)
    ax.set_ylabel("Mean entropy")
    ax.set_xlabel("Day")
    pm.apply_standard_axis_style(ax)
    
    ax = axs[1, 1]
    ax.plot(days, mec2, color="black", linewidth=1)
    ax.fill_between(days, mec2 - stdc2, mec2 + stdc2, color="gray", alpha=0.3)
    ax.set_ylabel("Mean entropy")
    ax.set_xlabel("Day")
    pm.apply_standard_axis_style(ax)
    
def plot_figure4(axs, data):
    """
    Plot 2x1:
    Row 1: boxplot of max entropy (lineages)
    Row 2: boxplot of max entropy (clusters)
    """
    mel1_list = data["max_entropy_lineages"][1]
    mel2_list = data["max_entropy_lineages"][2]
    mec1_list = data["max_entropy_clusters"][1]
    mec2_list = data["max_entropy_clusters"][2]

    df_me_l = pd.DataFrame({
        "Group": ["linear phenotype model"] * len(mel1_list) +
                 ["immune waning phenotype model"] * len(mel2_list),
        "MaxEntropy": mel1_list + mel2_list
    })
    sns.boxplot(ax=axs[0], data=df_me_l, x="Group", y="MaxEntropy",hue="Group")
    axs[0].set_ylabel("Max entropy")
    pm.apply_standard_axis_style(axs[0])
    
    df_me_c = pd.DataFrame({
        "Group": ["linear phenotype model"] * len(mec1_list) +
                 ["immune waning phenotype model"] * len(mec2_list),
        "MaxEntropy": mec1_list + mec2_list
    })
    sns.boxplot(ax=axs[1], data=df_me_c, x="Group", y="MaxEntropy",hue="Group")
    axs[1].set_ylabel("Max entropy")
    pm.apply_standard_axis_style(axs[1])

def plot_figure5(axs, data):
    """
    Plot 2×1:
    Row 1: boxplot of last-time entropy (lineages)
    Row 2: boxplot of last-time entropy (clusters)
    """
    lel1_list = data["last_entropy_lineages"][1]
    lel2_list = data["last_entropy_lineages"][2]
    lec1_list = data["last_entropy_clusters"][1]
    lec2_list = data["last_entropy_clusters"][2]

    order = ["linear phenotype model", "immune waning phenotype model"]

    df_le_l = pd.DataFrame({
        "Group": ["linear phenotype model"] * len(lel1_list) +
                 ["immune waning phenotype model"] * len(lel2_list),
        "LastEntropy": lel1_list + lel2_list
    }).dropna(subset=["LastEntropy"])

    sns.boxplot(
        ax=axs[0],
        data=df_le_l,
        x="Group",
        y="LastEntropy",
        hue="Group"
    )
    axs[0].set_ylabel("Last entropy")
    pm.apply_standard_axis_style(axs[0])

    df_le_c = pd.DataFrame({
        "Group": ["linear phenotype model"] * len(lec1_list) +
                 ["immune waning phenotype model"] * len(lec2_list),
        "LastEntropy": lec1_list + lec2_list
    }).dropna(subset=["LastEntropy"])

    sns.boxplot(
        ax=axs[1],
        data=df_le_c,
        x="Group",
        y="LastEntropy",
        order=order,
        hue="Group"
    )
    axs[1].set_ylabel("Last entropy")
    pm.apply_standard_axis_style(axs[1])

def summarize_statistical_tests_from_dfs(df2, df4, df5):
    """
    Concatenate three pre‐computed summary_dfs (one per figure) into a single table.
    Each df_{2,4,5} is assumed to have columns:
      ["Metric", "U_statistic", "p_value", "median_diff", "direction"]
    """
    # Tag each with its source figure
    df2l = df2.copy()

    df4l = df4.copy()

    df5l = df5.copy()

    # Concatenate and reorder
    combined = pd.concat([df2l, df4l, df5l], ignore_index=True)
    combined = combined[["Metric", "p_value", "median_diff", "direction"]]
    return combined

# -----------------------------------------------------------------------------
# Main Function (defines all constants and passes them as arguments)
# -----------------------------------------------------------------------------

def main():
    # Define constants locally
    EXPERIMENT_NAME     = "SIMPLICITY_exp_output"
    MIN_SIM_FINAL_TIME  = 300
    MAX_SIM_FINAL_TIME  = 1000
    CUT_AFTER_MAX_TIME  = True
    SEED                = 7
    CLUSTER_THRESHOLD   = 5
    TAKEOVER_THRESHOLD  = 0.50

    # Locate simulation output directories
    sim_dirs = dm.get_simulation_output_dirs(EXPERIMENT_NAME)
    sim_out_dir_1, sim_out_dir_2 = sim_dirs[0], sim_dirs[1]

    # Filter seeds for each group
    filtered_seeds_1 = filter_seeds(sim_out_dir_1, MIN_SIM_FINAL_TIME, MAX_SIM_FINAL_TIME, CUT_AFTER_MAX_TIME)
    filtered_seeds_2 = filter_seeds(sim_out_dir_2, MIN_SIM_FINAL_TIME, MAX_SIM_FINAL_TIME, CUT_AFTER_MAX_TIME)

    # Build per-seed caches (only loads/clusters/entropy once per seed)
    cache_1 = build_seed_cache(sim_out_dir_1, filtered_seeds_1,
                               max_time=MAX_SIM_FINAL_TIME,
                               cut_after_max=CUT_AFTER_MAX_TIME,
                               cluster_threshold=CLUSTER_THRESHOLD)

    cache_2 = build_seed_cache(sim_out_dir_2, filtered_seeds_2,
                               max_time=MAX_SIM_FINAL_TIME,
                               cut_after_max=CUT_AFTER_MAX_TIME,
                               cluster_threshold=CLUSTER_THRESHOLD)

    # Figure 1
    fig1_data = prepare_figure1_data(sim_out_dir_1, sim_out_dir_2,
                                     seed=SEED,
                                     cache_1=cache_1,
                                     cache_2=cache_2)
    fig1, axs1 = plt.subplots(2, 2, figsize=(20, 20), sharex='col',constrained_layout=True)
    plot_figure1(axs1, fig1_data)
    # plt.tight_layout()
    plt.show()

    # Figure 2
    fig2_data = prepare_figure2_data(sim_out_dir_1, sim_out_dir_2,
                                     filtered_seeds_1, filtered_seeds_2,
                                     cache_1, cache_2,
                                     takeover_threshold=TAKEOVER_THRESHOLD)
    fig2, axs2 = plt.subplots(2, 1, figsize=(20, 20), sharex=False,constrained_layout=True)
    plot_figure2(axs2, fig2_data)
    # plt.tight_layout()
    plt.show()

    # # Figure 3
    # fig3_data = prepare_figure3_data(sim_out_dir_1, sim_out_dir_2,
    #                                  filtered_seeds_1, filtered_seeds_2,
    #                                  cache_1, cache_2,
    #                                  max_time=MAX_SIM_FINAL_TIME)
    # fig3, axs3 = plt.subplots(2, 2, figsize=(20, 20), sharex='col',constrained_layout=True)
    # plot_figure3(axs3, fig3_data)
    # # plt.tight_layout()
    # plt.show()

    # # Figure 4
    # fig4_data = prepare_figure4_data(sim_out_dir_1, sim_out_dir_2,
    #                                  filtered_seeds_1, filtered_seeds_2,
    #                                  cache_1, cache_2)
    # fig4, axs4 = plt.subplots(2, 1, figsize=(20, 20), sharex=False,constrained_layout=True)
    # plot_figure4(axs4, fig4_data)
    # # plt.tight_layout()
    # plt.show()

    # # Figure 5
    # fig5_data = prepare_figure5_data(sim_out_dir_1, sim_out_dir_2,
    #                                  filtered_seeds_1, filtered_seeds_2,
    #                                  cache_1, cache_2)
    # fig5, axs5 = plt.subplots(2, 1, figsize=(20, 20), sharex=False,constrained_layout=True)
    # plot_figure5(axs5, fig5_data)
    # # plt.tight_layout()
    # plt.show()
    
    # # test results
    # df2 = fig2_data["summary_df"]
    # df4 = fig4_data["summary_df"]
    # df5 = fig5_data["summary_df"]
    # stats_table = summarize_statistical_tests_from_dfs(df2, df4, df5)
    # print(stats_table)

if __name__ == "__main__":
    main()