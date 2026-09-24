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
# Scenario lists come from the pipeline config via _scenarios.resolve_scenarios;
# callers pass what they resolved. Nothing here hardcodes a list any more.


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
def get_panel_a_data(exp_num, group, cluster_threshold=5, min_days=100,
                     exp_name=EXP_NAME):
    sod = _experiment_sod(exp_num, group, exp_name=exp_name)
    if sod is None:
        return {}, 0
    return summarize_sod_pies(sod, cluster_threshold, min_days)


# =============================================================================
# Panel B -- metrics PCA, one point per seed, all scenarios
# =============================================================================
_PANEL_B_CACHE = {}


def get_panel_b_data(exp_num, scenarios, cluster_threshold=5, min_days=100,
                     exp_name=EXP_NAME):
    """Clade metrics for every seed of every scenario.

    Cached on disk: this is identical for every --group of a given arm, and
    rendering one figure per scenario would otherwise recompute clade
    clustering for all ~250 seed-scenarios each time.
    """
    key = (exp_name, exp_num, cluster_threshold, min_days, tuple(sorted(scenarios)))
    if key in _PANEL_B_CACHE:
        return _PANEL_B_CACHE[key]
    cache_dir = os.path.join(dm.get_data_dir(), "figures", ".cache")
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(
        cache_dir,
        f"panelB_{exp_name}_{exp_num}_{cluster_threshold}_{min_days}"
        f"_{len(scenarios)}.pkl")
    if os.path.exists(cache_file):
        try:
            df = pd.read_pickle(cache_file)
            _PANEL_B_CACHE[key] = df
            print(f"[fig3] panel B: {len(df)} rows from cache")
            return df
        except Exception:
            pass

    rows = []
    for scenario in scenarios:
        sod = _experiment_sod(exp_num, scenario, exp_name=exp_name)
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
    df = pd.DataFrame(rows, columns=["scenario", "seed", "Peak", "Burden", "Survival", "Growth"])
    try:
        df.to_pickle(cache_file)
    except Exception:
        pass
    _PANEL_B_CACHE[key] = df
    return df


# =============================================================================
# Panel C -- sequence-space PCA + per-seed consistency, single scenario
# =============================================================================
def _ih_lineage_genomes(ssod, individual_type, max_per_individual=None,
                        rng_seed=42):
    """[(individual_index, genome), ...] -- EVERY intra-host lineage of every
    individual of `individual_type`, joined to its genome.

    Both types come through here, which is the point. The previous code took
    all intra-host lineages for long shedders but only the DIAGNOSIS-sequenced
    genomes for standards, and at sequencing_rate=0.05 that is ~6 sequences a
    seed against ~1,494 lineages -- an asymmetry that inflated every
    long-vs-standard comparison built on it.

    Note this measures the diversity each host type CARRIES, not what
    surveillance would DETECT. See BACKLOG: the data-source choice is to be
    re-checked before the final figures.
    """
    ind_df = om.read_individuals_data(ssod)
    phylo_df = om.read_phylogenetic_data(ssod)
    lin2gen = dict(zip(phylo_df["Lineage_name"], phylo_df["Genome"]))

    rng = random.Random(rng_seed)
    out = []
    sub = ind_df[ind_df['type'] == individual_type]
    for idx, row in sub.iterrows():
        lineages = [l for l in row['IH_lineages'] if l in lin2gen]
        if max_per_individual is not None and len(lineages) > max_per_individual:
            lineages = rng.sample(lineages, max_per_individual)
        for lin in lineages:
            out.append((idx, lin2gen[lin]))
    return out


def _standard_sequenced_genomes(ssod, max_per_individual=None):
    return _ih_lineage_genomes(ssod, 'standard', max_per_individual)


def _long_shedder_sequenced_genomes(ssod, max_per_individual=None, rng_seed=42):
    return _ih_lineage_genomes(ssod, 'long_shedder', max_per_individual, rng_seed)


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


def get_panel_c_scatter_data(exp_num, group, seed, max_long_per_individual=None,
                             exp_name=EXP_NAME):
    """One seed's sequenced genomes as a SNP matrix + individual_type labels."""
    sod = _experiment_sod(exp_num, group, exp_name=exp_name)
    if sod is None:
        return pd.DataFrame(), pd.Series(dtype=str)

    ssod = dm.get_ssod(sod, int(seed))
    standard = _standard_sequenced_genomes(ssod)
    long = _long_shedder_sequenced_genomes(ssod, max_long_per_individual)

    genomes = [g for _, g in standard] + [g for _, g in long]
    labels = (["standard"] * len(standard)) + (["long_shedder"] * len(long))

    snp_df = _build_snp_matrix(genomes)
    return snp_df, pd.Series(labels, name="individual_type")


def get_panel_c_consistency_data(exp_num, group, max_long_per_individual=None,
                                 max_pairs=200, rng_seed=42, exp_name=EXP_NAME):
    """Per seed, the THREE mean pairwise Hamming distances:

        standard-standard   how spread the standard population is
        long-long           how spread the long-shedder population is
        long-standard       how far apart the two are

    The panel used to report only (long-standard) minus (standard-standard),
    which invites reading a positive value as "long-shedder virus is
    diverged". It is not: long-long is the LARGEST of the three in most seeds,
    which is the signature of a BROADER cloud around the same centre, not a
    displaced one. Reporting all three makes that visible instead of hiding it
    behind a single difference.

    Computed per seed, never pooling genomes across seeds, so there is no
    cross-seed batch effect. Pair INDICES are sampled rather than materialising
    the product (a seed holds thousands of genomes per type).
    """
    sod = _experiment_sod(exp_num, group, exp_name=exp_name)
    if sod is None:
        return pd.DataFrame(columns=["seed", "comparison", "distance"])

    rng = random.Random(rng_seed)

    def mean_within(genomes):
        n = len(genomes)
        if n < 2:
            return np.nan
        total = n * (n - 1) // 2
        if total > max_pairs:
            vals = []
            while len(vals) < max_pairs:
                i, j = rng.randrange(n), rng.randrange(n)
                if i != j:
                    vals.append(hamming_iw(genomes[i], genomes[j]))
        else:
            vals = [hamming_iw(a, b) for a, b in itertools.combinations(genomes, 2)]
        return float(np.mean(vals)) if vals else np.nan

    def mean_between(a_list, b_list):
        na, nb = len(a_list), len(b_list)
        if na == 0 or nb == 0:
            return np.nan
        total = na * nb
        if total > max_pairs:
            idx = rng.sample(range(total), max_pairs)
            vals = [hamming_iw(a_list[i // nb], b_list[i % nb]) for i in idx]
        else:
            vals = [hamming_iw(a, b) for a in a_list for b in b_list]
        return float(np.mean(vals)) if vals else np.nan

    rows = []
    for ssod in dm.get_seeded_simulation_output_dirs(sod):
        try:
            standard = [g for _, g in _standard_sequenced_genomes(ssod, max_long_per_individual)]
            long = [g for _, g in _long_shedder_sequenced_genomes(ssod, max_long_per_individual)]
            if len(standard) < 2 or len(long) < 2:
                continue
            seed = dm.get_seed_from_SSOD(ssod)
            for label, value in (
                ("standard-standard", mean_within(standard)),
                ("long-long", mean_within(long)),
                ("long-standard", mean_between(long, standard)),
            ):
                if np.isfinite(value):
                    rows.append({"seed": seed, "comparison": label, "distance": value})
        except Exception:
            continue

    return pd.DataFrame(rows, columns=["seed", "comparison", "distance"])
