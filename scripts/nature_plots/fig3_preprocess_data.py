import os
import sys
import ast
import random
import itertools

import numpy as np
import pandas as pd

import simplicity.dir_manager as dm
import simplicity.output_manager as om
from simplicity.phenotype.distance import hamming_iw

# get_clade_metrics/get_clade_winners/summarize_sod_pies already work against
# plain sod/ssod paths (no dependency on the old M/R/ratio/tau grid), so we
# reuse them here rather than duplicating the clade-clustering logic. That
# module lives in a sibling top-level scripts/ directory, not a package, so
# it needs its own directory on sys.path (same cross-directory sibling-import
# pattern already used by long_nsr_calibration_plot.py).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'long_paper_figures'))
from long_shedders_preprocess import get_clade_metrics, summarize_sod_pies

EXP_NAME = "impact_long_shedders"
SCENARIOS = ["control", "SOT", "HIV_low", "HIV_high"]


def _experiment_sod(exp_num, scenario, exp_name=EXP_NAME):
    """
    First (only) simulation_output_dir for one scenario's experiment, or None
    if that scenario hasn't been run yet at all -- dm.get_simulation_output_dirs
    raises (rather than returning []) when the experiment folder doesn't
    exist, so callers looping over scenarios (e.g. get_panel_b_data) can
    still skip a not-yet-run scenario instead of crashing the whole figure.
    """
    experiment_string = f"{exp_name}_{scenario}_#{exp_num}"
    try:
        sods = dm.get_simulation_output_dirs(experiment_string)
    except ValueError:
        return None
    return sods[0] if sods else None


# =============================================================================
# Panel A -- pies for one scenario
# =============================================================================
def get_panel_a_data(exp_num, group, cluster_threshold=5, min_days=100):
    sod = _experiment_sod(exp_num, group)
    if sod is None:
        return {}, 0
    return summarize_sod_pies(sod, cluster_threshold, min_days)


# =============================================================================
# Panel B -- metrics PCA, one point per seed, all scenarios
# =============================================================================
def get_panel_b_data(exp_num, scenarios=SCENARIOS, cluster_threshold=5, min_days=100):
    rows = []
    for scenario in scenarios:
        sod = _experiment_sod(exp_num, scenario)
        if sod is None:
            continue
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            try:
                if om.read_final_time(ssod) < min_days:
                    continue
                values = get_clade_metrics(ssod, cluster_threshold)["values"]
                rows.append({
                    "scenario": scenario,
                    "seed": dm.get_seed_from_SSOD(ssod),
                    **values,
                })
            except Exception:
                continue
    return pd.DataFrame(rows, columns=["scenario", "seed", "Peak", "Burden", "Survival", "Growth"])


# =============================================================================
# Panel C -- sequence-space PCA + per-seed consistency, single scenario
# =============================================================================
def _standard_sequenced_genomes(ssod):
    """
    [(individual_index, genome), ...] for standard-surveillance sequenced
    reads. sequencing_data.csv's own 'sequence' column already contains the
    genome (simplicity.population_model embeds it at diagnosis-sampling time)
    -- no join to phylogenetic_data.csv needed.
    """
    path = os.path.join(ssod, 'sequencing_data.csv')
    if not os.path.isfile(path):
        return []
    df = pd.read_csv(path)
    df = df[df['individual_type'] == 'standard']
    out = []
    for _, row in df.iterrows():
        try:
            genome = ast.literal_eval(row['sequence']) if isinstance(row['sequence'], str) else row['sequence']
            out.append((row['individual_index'], genome))
        except Exception:
            continue
    return out


def _long_shedder_sequenced_genomes(ssod, max_per_individual=None, rng_seed=42):
    """
    [(individual_index, genome), ...] for every IH lineage of every long
    shedder -- mirrors save_sequencing_dataset's long-shedder export exactly
    (every IH lineage gets "sequenced" at end of infection). Needs the join
    to phylogenetic_data.csv's Genome, since individuals_data.csv only
    stores lineage names.
    """
    ind_df = om.read_individuals_data(ssod)
    phylo_df = om.read_phylogenetic_data(ssod)
    lin2gen = dict(zip(phylo_df["Lineage_name"], phylo_df["Genome"]))

    rng = random.Random(rng_seed)
    out = []
    long_df = ind_df[ind_df['type'] == 'long_shedder']
    for idx, row in long_df.iterrows():
        lineages = [l for l in row['IH_lineages'] if l in lin2gen]
        if max_per_individual is not None and len(lineages) > max_per_individual:
            lineages = rng.sample(lineages, max_per_individual)
        for lin in lineages:
            out.append((idx, lin2gen[lin]))
    return out


def _build_snp_matrix(genomes):
    """
    genomes: list of sparse (position, mutation) lists. Returns a presence/
    absence DataFrame: columns = union of mutated positions observed across
    all genomes, one row per genome (1 = mutated at that position, any
    allele; 0 = matches reference).
    """
    positions = sorted({pos for genome in genomes for pos, _ in genome})
    rows = [[1 if pos in {p for p, _ in genome} else 0 for pos in positions] for genome in genomes]
    return pd.DataFrame(rows, columns=positions)


def get_panel_c_scatter_data(exp_num, group, seed, max_long_per_individual=None):
    """One seed's sequenced genomes as a SNP matrix + individual_type labels."""
    sod = _experiment_sod(exp_num, group)
    if sod is None:
        return pd.DataFrame(), pd.Series(dtype=str)

    ssod = dm.get_ssod(sod, int(seed))
    standard = _standard_sequenced_genomes(ssod)
    long = _long_shedder_sequenced_genomes(ssod, max_long_per_individual)

    genomes = [g for _, g in standard] + [g for _, g in long]
    labels = (["standard"] * len(standard)) + (["long_shedder"] * len(long))

    snp_df = _build_snp_matrix(genomes)
    return snp_df, pd.Series(labels, name="individual_type")


def get_panel_c_consistency_data(exp_num, group, max_long_per_individual=None, max_pairs=200, rng_seed=42):
    """
    Per seed of `group`, independently: mean pairwise hamming_iw(long,
    standard) minus mean pairwise hamming_iw(standard, standard) = "excess
    divergence". Computed per seed (never pooling genomes across seeds) so
    there's no cross-seed batch-effect risk. Pairs are capped/subsampled
    since this is O(n*m) per seed.
    """
    sod = _experiment_sod(exp_num, group)
    if sod is None:
        return pd.DataFrame(columns=["seed", "excess_divergence"])

    rng = random.Random(rng_seed)
    rows = []
    for ssod in dm.get_seeded_simulation_output_dirs(sod):
        try:
            standard = [g for _, g in _standard_sequenced_genomes(ssod)]
            long = [g for _, g in _long_shedder_sequenced_genomes(ssod, max_long_per_individual)]
            if not standard or not long:
                continue

            between_pairs = [(l, s) for l in long for s in standard]
            if len(between_pairs) > max_pairs:
                between_pairs = rng.sample(between_pairs, max_pairs)
            between_dists = [hamming_iw(l, s) for l, s in between_pairs]

            within_pairs = list(itertools.combinations(standard, 2))
            if len(within_pairs) > max_pairs:
                within_pairs = rng.sample(within_pairs, max_pairs)
            within_dists = [hamming_iw(a, b) for a, b in within_pairs]

            if not between_dists or not within_dists:
                continue

            excess = float(np.mean(between_dists) - np.mean(within_dists))
            rows.append({"seed": dm.get_seed_from_SSOD(ssod), "excess_divergence": excess})
        except Exception:
            continue

    return pd.DataFrame(rows, columns=["seed", "excess_divergence"])
