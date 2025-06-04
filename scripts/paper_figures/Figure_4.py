# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from tqdm import tqdm
from collections import defaultdict, deque

import simplicity.output_manager as om
import simplicity.dir_manager as dm
import simplicity.plots_manager as pm

pm.apply_plos_rcparams()

# -----------------------------------------------------------------------------
# Utility Functions
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
    Cluster lineages based on shared mutations.
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

    freq_df_clusters = (pd.concat(cluster_freqs, axis=1)
                        if cluster_freqs else pd.DataFrame(index=freq_pivot.index))
    return freq_df_clusters, clusters, parents

def build_clustered_freqs(lineage_to_mutations, freq_df_lineages, threshold):
    """
    Cluster lineages by shared mutations with threshold, returning:
      - freq_df_clusters: DataFrame of cluster‐summed frequencies
      - filtered_parents: list of parent lineage names (one per column in freq_df_clusters)
    """
    if threshold == 0:
        return freq_df_lineages.copy(), list(freq_df_lineages.columns)

    # Get all clusters + their “parent” from the core function
    freq_all, clusters, parents_all = cluster_lineages_by_shared_mutations(
        lineage_to_mutations, freq_df_lineages, threshold
    )

    # re‐build freq_df_clusters and only keep the parents
    cluster_freqs = []
    filtered_parents = []
    for i, cluster in enumerate(clusters):
        cols = [col for col in cluster if col in freq_df_lineages.columns]
        if not cols:
            continue
        summed = freq_df_lineages[cols].sum(axis=1).rename(parents_all[i])
        cluster_freqs.append(summed)
        filtered_parents.append(parents_all[i])

    if cluster_freqs:
        freq_df_clusters = pd.concat(cluster_freqs, axis=1)
    else:
        freq_df_clusters = pd.DataFrame(index=freq_df_lineages.index)

    return freq_df_clusters, filtered_parents

def compute_shannon_entropy(freq_df):
    """
    Compute Shannon entropy across columns for each row in freq_df.
    """
    normed = freq_df.div(freq_df.sum(axis=1), axis=0)
    log_vals = np.log(normed.where(normed > 0))
    entropy = -(normed * log_vals).sum(axis=1, skipna=True)
    entropy[normed.isna().all(axis=1)] = np.nan
    return entropy


def count_sweeps(series_bool, min_days_above_threshold=21):
    """
    Given a boolean Series (True where normalized frequency ≥ threshold),
    count how many sustained runs of length ≥ min_days_above_threshold occur.
    """
    arr = series_bool.values.astype(int)
    total, run_len = 0, 0
    for v in arr:
        if v == 1:
            run_len += 1
        else:
            if run_len >= min_days_above_threshold:
                total += 1
            run_len = 0
    if run_len >= min_days_above_threshold:
        total += 1
    return total


# -----------------------------------------------------------------------------
# Plotting Helpers
# -----------------------------------------------------------------------------

def draw_stack_entropy(ax, freq_df, entropy_series, colormap, data_type):
    """
    ax:             A single Axes object
    freq_df:        DataFrame indexed by time; each column is a lineage (frequency over time)
    entropy_series: Series indexed by time (Shannon entropy at each time point)
    colormap:       result of pm.make_lineages_colormap(...) for coloring each lineage
    data_type:      lineage or cluster, used for axis label
    """

    # Entropy 
    ax.plot(entropy_series.index, entropy_series.values, "k-", linewidth=2)
    ax.set_ylabel("Entropy")
    ax.set_zorder(1)
    ax.patch.set_visible(False) 
    ax.set_xlabel("Time (d)")

    # Stackplot of frequencies
    ax2 = ax.twinx()
    freq_nonan = freq_df.fillna(0)
    colors = [pm.get_lineage_color(lin, colormap) for lin in freq_nonan.columns]
    ax2.stackplot(
        freq_nonan.index,
        freq_nonan.T.values,
        colors=colors,
        alpha=0.5,
        edgecolors="black",
        linewidth=0.1
    )
    ax2.set_ylabel(f"{data_type} Frequency")
    ax2.set_ylim(0,1)
    ax2.set_xlim(0,max(freq_nonan.index))
    ax2.set_zorder(0)
    
    pm.apply_standard_axis_style(ax, True)


def draw_violin(ax, counts_group1, counts_group2, ylabel):
    """
    ax:              A single Axes object
    counts_group1:   list of ints (e.g. tcl1 or tcc1)
    counts_group2:   list of ints (e.g. tcl2 or tcc2)
    ylabel:          string for the y-axis label
    """
    df = pd.DataFrame({
        "Group": ["Linear model  "] * len(counts_group1) + ["  Immune waning model***"] * len(counts_group2),
        "Count": counts_group1 + counts_group2
    })
    sns.violinplot(
        ax=ax,
        data=df,
        x="Group",
        y="Count",
        inner="box",
        cut=0,
        hue="Group",
        legend=False
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    pm.apply_standard_axis_style(ax)


# -----------------------------------------------------------------------------
# Figure 4 
# -----------------------------------------------------------------------------

def plot_figure_4(include_clusters=False):
    """
    Layout when include_clusters=False:
      ┌───────────────────────────────┬───────────────────────────┐
      │ (0,0) Linear raw+entropy      │ (0:2,1) lineage‐level     │
      │                               │         violin            │
      ├───────────────────────────────┼───────────────────────────┤
      │ (1,0) Immune raw+entropy      │                           │
      └───────────────────────────────┴───────────────────────────┘

    Layout when include_clusters=True:
      ┌───────────────────────────────┬───────────────────────────┐
      │ (0,0) Linear raw+entropy      │ (0:2,1) lineage‐level     │
      │                               │         violin            │
      ├───────────────────────────────┤                           │
      │ (1,0) Immune raw+entropy      │                           │
      ├───────────────────────────────┤───────────────────────────┤
      │ (2,0) Linear clustered+entropy│ (2:4,1) cluster‐level     │
      │                               │         violin            │
      ├───────────────────────────────┤                           │
      │ (3,0) Immune clustered+entropy│                           │
      └───────────────────────────────┴───────────────────────────┘
    """
    EXPERIMENT_NAME     = "SIMPLICITY_exp_output"
    MIN_FINAL_TIME      = 300
    MAX_FINAL_TIME      = 450
    CUT_AFTER_MAX_TIME  = True
    SEED                = 7
    CLUSTER_THRESHOLD   = 5
    TAKEOVER_THRESHOLD  = 0.50

    sim_dirs = dm.get_simulation_output_dirs(EXPERIMENT_NAME)
    sim_out_dir_1, sim_out_dir_2 = sim_dirs[0], sim_dirs[1]

    # Filter seeds for sweeps calculations 
    filtered_seeds_1 = filter_seeds(
        sim_out_dir_1,
        MIN_FINAL_TIME,
        MAX_FINAL_TIME,
        CUT_AFTER_MAX_TIME
    )
    filtered_seeds_2 = filter_seeds(
        sim_out_dir_2,
        MIN_FINAL_TIME,
        MAX_FINAL_TIME,
        CUT_AFTER_MAX_TIME
    )

    # Load seed data for the left‐column panels 
    # Group 1: Linear phenotype model
    l2m1, freq_l1 = load_seed_data(
        sim_out_dir_1,
        SEED,
        MAX_FINAL_TIME,
        CUT_AFTER_MAX_TIME
    )
    freq_c1, parents1 = build_clustered_freqs(l2m1, freq_l1, CLUSTER_THRESHOLD)
    freq_c1.columns = parents1
    ent_l1 = compute_shannon_entropy(freq_l1)
    ent_c1 = compute_shannon_entropy(freq_c1)
    cmap1 = pm.make_lineages_colormap(dm.get_ssod(sim_out_dir_1, SEED))

    # Group 2: Immune waning phenotype model
    l2m2, freq_l2 = load_seed_data(
        sim_out_dir_2,
        SEED,
        MAX_FINAL_TIME,
        CUT_AFTER_MAX_TIME
    )
    freq_c2, parents2 = build_clustered_freqs(l2m2, freq_l2, CLUSTER_THRESHOLD)
    freq_c2.columns = parents2
    ent_l2 = compute_shannon_entropy(freq_l2)
    ent_c2 = compute_shannon_entropy(freq_c2)
    cmap2 = pm.make_lineages_colormap(dm.get_ssod(sim_out_dir_2, SEED))

    # Build selective‐sweeps lists for violin plots 
    scl1, scc1 = [], []
    for seed in tqdm(filtered_seeds_1, desc="Computing sweeps for Group 1"):
        l2m, f_l = load_seed_data(
            sim_out_dir_1,
            seed,
            MAX_FINAL_TIME,
            CUT_AFTER_MAX_TIME
        )
        f_c, _ = build_clustered_freqs(l2m, f_l, CLUSTER_THRESHOLD)

        norm_l = f_l.div(f_l.sum(axis=1), axis=0)
        mask_l = (norm_l >= TAKEOVER_THRESHOLD)
        scl1.append(int(mask_l.apply(count_sweeps, axis=0).sum()))

        norm_c = f_c.div(f_c.sum(axis=1), axis=0)
        mask_c = (norm_c >= TAKEOVER_THRESHOLD)
        scc1.append(int(mask_c.apply(count_sweeps, axis=0).sum()))

    scl2, scc2 = [], []
    for seed in tqdm(filtered_seeds_2, desc="Computing sweeps for Group 2"):
        l2m, f_l = load_seed_data(
            sim_out_dir_2,
            seed,
            MAX_FINAL_TIME,
            CUT_AFTER_MAX_TIME
        )
        f_c, _ = build_clustered_freqs(l2m, f_l, CLUSTER_THRESHOLD)

        norm_l = f_l.div(f_l.sum(axis=1), axis=0)
        mask_l = (norm_l >= TAKEOVER_THRESHOLD)
        scl2.append(int(mask_l.apply(count_sweeps, axis=0).sum()))

        norm_c = f_c.div(f_c.sum(axis=1), axis=0)
        mask_c = (norm_c >= TAKEOVER_THRESHOLD)
        scc2.append(int(mask_c.apply(count_sweeps, axis=0).sum()))

    # Build the figure + GridSpec depending on include_clusters ─═══════
    if not include_clusters:
        fig_name = f"Figure_4_{EXPERIMENT_NAME}_seed{SEED}.tiff"
        # Layout: 2 rows × 2 cols
        # Left column: two stack‐+entropy panels
        # Right column: one violin spanning both rows
        fig = plt.figure(figsize=(16, 10))
        gs = gridspec.GridSpec(
            nrows=2, ncols=2,
            height_ratios=[1, 1],
            width_ratios=[3, 2],
            hspace=0.3, wspace=0.3,
            figure=fig
        )

        # Left, top (0,0): Linear raw + entropy
        ax00 = fig.add_subplot(gs[0, 0])
        ax00.set_title("Linear phenotype model")
        draw_stack_entropy(
            ax00,
            freq_l1, ent_l1, cmap1, data_type='Lineage'
        )

        # Left, bottom (1,0): Immune raw + entropy
        ax10 = fig.add_subplot(gs[1, 0])
        ax10.set_title("Immune waning phenotype model")
        draw_stack_entropy(
            ax10,
            freq_l2, ent_l2, cmap2, data_type='Lineage'
        )

        # Right (spans rows 0–1, col 1): lineage‐level violin
        ax_violin = fig.add_subplot(gs[:, 1])
        draw_violin(ax_violin, scl1, scl2, "Number of selective sweeps")

    else:
        fig_name = f"Figure_4S_{EXPERIMENT_NAME}_seed{SEED}.tiff"
        # Layout: 4 rows × 2 cols
        # Left column: four stack‐+entropy panels
        # Right column: two tall violin panels
        fig = plt.figure(figsize=(16, 16))
        gs = gridspec.GridSpec(
        nrows=5, ncols=2,
        height_ratios=[1, 1, 0.05, 1, 1],   # ← row 2 is a thin spacer
        width_ratios=[3, 2],
        hspace=0.3, wspace=0.3,
        figure=fig
        )

        # Left, row 0 (freq_l1 + ent_l1)
        ax00 = fig.add_subplot(gs[0, 0])
        ax00.set_title("Linear phenotype model")
        draw_stack_entropy(
            ax00,
            freq_l1, ent_l1, cmap1, data_type='Lineage'
        )

        # Left, row 1 (freq_c1 + ent_c1)
        ax01 = fig.add_subplot(gs[1, 0])
        ax01.set_title("Immune waning phenotype model")
        draw_stack_entropy(
            ax01,
            freq_l2, ent_l2, cmap2, data_type='Lineage'
        )
        # Right, rows 0–1: lineage‐level violin
        ax_violin1 = fig.add_subplot(gs[0:2, 1])
        draw_violin(ax_violin1, scl1, scl2, "Number of selective sweeps")
        
        
        # Row 2 - SPACER

        # Left, row 2 (freq_l2 + ent_l2)
        ax02 = fig.add_subplot(gs[3, 0])
        ax02.set_title("Linear phenotype model")
        draw_stack_entropy(
            ax02,
            freq_c1, ent_c1, cmap1, data_type='Clustered Lineage'
        )
      
        # Left, row 3 (freq_c2 + ent_c2)
        ax03 = fig.add_subplot(gs[4, 0])
        ax03.set_title("Immune waning phenotype model")
        draw_stack_entropy(
            ax03,
            freq_c2, ent_c2, cmap2, data_type='Clustered Lineage'
        )

        # Right, rows 2–3: cluster‐level violin
        ax_violin2 = fig.add_subplot(gs[3:5, 1])
        draw_violin(ax_violin2, scc1, scc2, "Number of selective sweeps")

   
    output_path = os.path.join("Data", fig_name)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    output_path = os.path.join("Data", f"{fig_name}.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    plt.close(fig)

if __name__ == "__main__":

    plot_figure_4(include_clusters=False)
    plot_figure_4(include_clusters=True)

