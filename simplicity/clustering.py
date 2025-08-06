#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict, deque, Counter
import pandas as pd


def parse_genome(genome):
    """
    Convert a genome (e.g. list of mutations) into a frozenset of mutation tuples.
    Safely handles None or malformed inputs.
    
    Args:
        genome (list or None): Genome data from a row in phylogenetic table.
    
    Returns:
        frozenset: Unique set of mutation tuples.
    """
    if genome is None:
        return frozenset()
    try:
        return frozenset(tuple(x) for x in genome)
    except Exception:
        return frozenset()


def build_lineage_to_mutation_dict(phylo_df):
    """
    Build a dictionary mapping each lineage to its parsed set of mutations.
    
    Args:
        phylo_df (pd.DataFrame): Must include "Lineage_name" and "Genome" columns.
    
    Returns:
        dict: Lineage_name → frozenset of mutations
    """
    return {
        row["Lineage_name"]: parse_genome(row["Genome"])
        for _, row in phylo_df.iterrows()
        if row["Genome"] is not None
    }


def cluster_lineages_by_shared_mutations(lin2mut: dict, freq_pivot: pd.DataFrame, min_shared: int):
    """
    Cluster lineages based on shared mutations using transitive closure.
    
    Args:
        lin2mut (dict): Lineage → frozenset of mutations
        freq_pivot (pd.DataFrame): Pivot table (Time_sampling x Lineage_name)
        min_shared (int): Minimum number of shared mutations for clustering
    
    Returns:
        freq_df_clusters (pd.DataFrame): Clustered frequency table
        clusters (list of sets): Each set = raw lineage names in that cluster
        parents (list): Representative lineage name for each cluster
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
    High-level wrapper for clustering and aggregating frequencies.
    
    Args:
        lineage_to_mutations (dict): Lineage → mutation set
        freq_df_lineages (pd.DataFrame): Pivoted frequency data
        threshold (int): Shared mutation count threshold
    
    Returns:
        freq_df_clusters (pd.DataFrame): Aggregated cluster frequency table
        cluster_parents (list): Parent lineages representing clusters
    """
    if threshold == 0:
        return freq_df_lineages.copy(), list(freq_df_lineages.columns)
    freq_df_clusters, clusters, parents = cluster_lineages_by_shared_mutations(
        lineage_to_mutations, freq_df_lineages, threshold
    )
    return freq_df_clusters, clusters, parents

def assign_cluster_hosttypes_from_parents(parents, cluster_names, host_type_series):
    """
    Assign each cluster the host type of its founder (parent) lineage.
    
    Args:
        parents (list): Representative lineage for each cluster (from clustering)
        cluster_names (list): Corresponding cluster names (e.g., "Cluster_0")
        host_type_series (pd.Series): Lineage_name → host type
    
    Returns:
        pd.Series: Cluster_name → founder's host type
    """
    cluster_hosttypes = {}
    for cname, parent in zip(cluster_names, parents):
        ht = host_type_series.get(parent)
        if ht is not None:
            cluster_hosttypes[cname] = ht
    return pd.Series(cluster_hosttypes)

def assign_cluster_hosttypes_by_max_mutation(clusters, lin2mut, host_type_series):
    """
    Assign cluster to the host type of the lineage with the most mutations.
    
    Args:
        clusters (list of sets): Each cluster as a set of lineage names
        lin2mut (dict): Lineage → mutation set
        host_type_series (pd.Series): Lineage → host type
    
    Returns:
        List of assigned host types per cluster
    """
    hosttypes = []
    for cluster in clusters:
        max_lin = max(cluster, key=lambda l: len(lin2mut.get(l, [])))
        hosttypes.append(host_type_series.get(max_lin, "Unknown"))
    return hosttypes

