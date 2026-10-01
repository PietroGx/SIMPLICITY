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
from long_shedders_preprocess import (get_clade_metrics, summarize_sod_pies,
                                     _clade_analysis, _peak_values, _resolve_label,
                                     _survival_durations, _growth_crossing_times)

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

    Cached on disk: this is identical for every --group of a given pipeline, and
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
    positions = sorted({pos for genome in genomes for pos in genome})
    rows = [[1 if pos in genome else 0 for pos in positions] for genome in genomes]
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


# =============================================================================
# Conversion efficiency -- the non-circular measurement
# =============================================================================
# Long shedders generate more substitutions BY CONSTRUCTION (NSR_long is ~13x
# NSR, IH_lineages_max is 5-15 against 1-4, infections last longer). Counting
# their lineages therefore recovers the parameter file, not a result: in
# SIMPLICITY one substitution IS one lineage.
#
# The question that is not circular: does that output convert into ESTABLISHED
# CLADES at a higher rate than the raw output predicts? Clustering is what
# makes a "clade" a stand-in for a real-world lineage -- several co-occurring
# substitutions, not a single SNV -- so it is essential here, not a nuisance.
#
#     efficiency(origin) = successful clades of that origin
#                          --------------------------------
#                          substitutions arising in those hosts
#
# Reported as long / standard. Above 1: long shedders convert BETTER than
# their mutational output predicts. Below 1: worse.
# =============================================================================

# =============================================================================
# Analysis window
# =============================================================================
# phenotype_model is "immune_waning", where a lineage's fitness is
#
#     infected_fraction * distance_from_consensus
#   + (1 - infected_fraction) / active_lineages_n
#
# with infected_fraction = min((diagnosed + recovered)/size, 1). While the
# population is immune-naive that first term vanishes and EVERY lineage scores
# 1/active_lineages_n -- selection on antigenic novelty is off and evolution is
# neutral. Measured over all 250 simulations, immunity reaches 10% at ~35 days,
# 50% at ~110 and 90% at ~200, consistently across scenarios.
#
# So every metric here is measured over the window AFTER the burn-in only.
# Clades are NOT filtered by birth time: a clade present from t=0 (the founder)
# stays in the comparison, but only the infections and frequencies it achieves
# after the cutoff are counted. That is what separates "won because it was
# there first, when nothing could out-compete it" from "won under selection".
#
# Simulations shorter than MIN_FINAL_TIME are early-extinction runs and are
# dropped entirely (14 of 243 in unbound #1).
# =============================================================================

BURNIN_CUTOFF_DAYS = 200
MIN_FINAL_TIME = 365


def _windowed_clade_data(ssod, cluster_threshold, cutoff=BURNIN_CUTOFF_DAYS,
                         ind=None):
    """(F, totals, labels) restricted to the post-burn-in window.

    F keeps only sampling times at or after the cutoff, so peak, survival and
    growth describe the selected phase. `totals` is rebuilt from infection
    times rather than read from the whole-run column, so burden counts only
    infections occurring after the cutoff.

    Note the pools differ by construction, as they do without a window: F holds
    non-root clades only, while totals covers every clade -- which is why
    founder clades compete for burden and for nothing else.
    """
    F, clade_to_lineages, _totals, labels = _clade_analysis(ssod, cluster_threshold)
    if F is not None and not F.empty:
        F = F.loc[F.index.astype(float) >= cutoff]
        F = F.dropna(axis=1, how="all")
        F = F.loc[:, (F.fillna(0) > 0).any(axis=0)]

    if ind is None:
        ind = om.read_individuals_data(ssod)
    post = ind[ind["t_infection"].astype(float) >= cutoff]
    per_lineage = post["inherited_lineage"].value_counts().to_dict()
    totals = {clade: sum(int(per_lineage.get(l, 0)) for l in members)
              for clade, members in clade_to_lineages.items()}
    totals = {c: v for c, v in totals.items() if v > 0}
    return F, totals, labels


def _usable(ssod, min_final_time=MIN_FINAL_TIME):
    """True if this simulation ran long enough to be worth measuring."""
    try:
        return float(om.read_final_time(ssod)) >= min_final_time
    except Exception:
        return False


EFFICIENCY_THRESHOLDS = [3, 5, 8, 10, 15]
EFFICIENCY_PEAKS = [0.10, 0.25, 0.50, 0.75]
REF_THRESHOLD = 5
REF_PEAK = 0.50


def _substitutions_by_type(ind_df, cutoff=None):
    """Substitution events arising in each host type. A lineage in a host's
    IH_lineages that is not the one it was infected with arose there.

    With `cutoff`, only lineages BORN after that time are counted, so the
    denominator covers the same window as the successes it is divided into.
    Birth times come from IH_lineages_trajectory; a lineage with no recorded
    birth (the host's inherited one) is excluded anyway.
    """
    out = {}
    for host_type, key in (('standard', 'standard'), ('long_shedder', 'long')):
        n = 0
        for _, row in ind_df[ind_df['type'] == host_type].iterrows():
            lineages = row['IH_lineages']
            if not isinstance(lineages, (list, tuple)):
                continue
            inherited = row.get('inherited_lineage')
            if cutoff is None:
                n += len([l for l in lineages if l != inherited])
                continue
            traj = row.get('IH_lineages_trajectory')
            t_inf = float(row.get('t_infection', 0.0) or 0.0)
            for l in lineages:
                if l == inherited:
                    continue
                birth = None
                if isinstance(traj, dict) and l in traj:
                    birth = traj[l].get('ih_birth')
                if birth is None:
                    continue
                if t_inf + float(birth) >= cutoff:
                    n += 1
        out[key] = n
    return out


def get_efficiency_data(exp_num, scenarios, exp_name=EXP_NAME, min_days=100):
    """Per seed x scenario x threshold: substitutions and successful clades by
    origin, for every peak criterion. Cached -- this walks every seed's
    individuals_data and re-clusters at five thresholds.
    """
    cache_dir = os.path.join(dm.get_data_dir(), "figures", ".cache")
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f"efficiency_{exp_name}_{exp_num}_w{BURNIN_CUTOFF_DAYS}_m{MIN_FINAL_TIME}.pkl")
    if os.path.exists(cache_file):
        try:
            df = pd.read_pickle(cache_file)
            print(f"[fig3] efficiency: {len(df)} rows from cache")
            return df
        except Exception:
            pass

    rows = []
    for scenario in scenarios:
        sod = _experiment_sod(exp_num, scenario, exp_name=exp_name)
        if sod is None:
            continue
        print(f"[fig3] efficiency: {scenario}...", flush=True)
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            if not _usable(ssod):
                continue
            try:
                ind = om.read_individuals_data(ssod)
                if ind.empty:
                    continue
                subs = _substitutions_by_type(ind, cutoff=BURNIN_CUTOFF_DAYS)
            except Exception:
                continue
            for ct in EFFICIENCY_THRESHOLDS:
                try:
                    F, _tot, labels = _windowed_clade_data(ssod, ct, ind=ind)
                    peaks = _peak_values(F)
                except Exception:
                    continue
                rec = {"scenario": scenario, "seed": dm.get_seed_from_SSOD(ssod),
                       "threshold": ct,
                       "subs_long": subs["long"], "subs_std": subs["standard"]}
                for pk in EFFICIENCY_PEAKS:
                    for origin, key in (("long", "long"), ("standard", "std")):
                        rec[f"succ_{key}_{pk}"] = sum(
                            1 for c in F.columns
                            if _resolve_label(labels, c) == origin
                            and peaks.get(c, 0.0) > pk)
                rows.append(rec)

    df = pd.DataFrame(rows)
    try:
        df.to_pickle(cache_file)
    except Exception:
        pass
    return df


def efficiency_ratio(df, threshold=REF_THRESHOLD, peak=REF_PEAK):
    """Pooled long/standard conversion efficiency per scenario, plus the raw
    per-1000-substitution rates behind it."""
    d = df[df["threshold"] == threshold]
    if d.empty:
        return pd.DataFrame()
    g = d.groupby("scenario").sum(numeric_only=True)
    out = pd.DataFrame(index=g.index)
    # Raw proportion: successes per substitution. No arbitrary scaling -- the
    # axis is labelled "1 in N" instead, which reads directly.
    out["eff_long"] = g[f"succ_long_{peak}"] / g["subs_long"].replace(0, np.nan)
    out["eff_std"] = g[f"succ_std_{peak}"] / g["subs_std"].replace(0, np.nan)
    out["ratio"] = out["eff_long"] / out["eff_std"]
    out["n_succ_long"] = g[f"succ_long_{peak}"]
    out["n_succ_std"] = g[f"succ_std_{peak}"]
    out["subs_long"] = g["subs_long"]
    out["subs_std"] = g["subs_std"]
    return out


def efficiency_robustness(df):
    """Ratio across every (threshold, peak) pair -- the sensitivity grid."""
    rows = []
    for ct in sorted(df["threshold"].unique()):
        for pk in EFFICIENCY_PEAKS:
            r = efficiency_ratio(df, threshold=ct, peak=pk)
            for scenario, rec in r.iterrows():
                rows.append({"scenario": scenario, "threshold": ct, "peak": pk,
                             "ratio": rec["ratio"]})
    return pd.DataFrame(rows)


# =============================================================================
# Win rate -- the within-seed ranking question
# =============================================================================
# The four clade metrics are only comparable BETWEEN CLADES OF ONE SEED: each
# simulation is a different random trajectory, so absolute values do not carry
# across seeds. The question they were built to answer is therefore a ranking
# one -- in this seed, did the top-scoring clade come from a long shedder? --
# which aggregates across seeds as a PROPORTION.
#
# The null is that origin carries no information about score: under
# exchangeability P(top clade is long) equals the long-origin share of the
# clades that were eligible to win.
#
# Eligibility differs by metric. Peak, Burden and Survival rank every clade.
# Growth ranks only clades that crossed BOTH 1% and 50% frequency, in order --
# 3-6 clades per seed against 49-84 for the others -- so its null must use that
# smaller pool or the enrichment is computed against the wrong denominator.
# =============================================================================

WIN_METRICS = ["Peak", "Burden", "Survival", "Growth"]


def _winner_and_pool(F, totals, metric):
    """(winning clade, clades eligible to win) for one metric."""
    if metric == "Peak":
        vals = _peak_values(F)
        return (max(vals, key=vals.get) if vals else None), list(vals)
    if metric == "Survival":
        vals = _survival_durations(F)
        return (max(vals, key=vals.get) if vals else None), list(vals)
    if metric == "Burden":
        return (max(totals, key=totals.get) if totals else None), list(totals)
    if metric == "Growth":
        vals = _growth_crossing_times(F)          # smaller is faster
        return (min(vals, key=vals.get) if vals else None), list(vals)
    raise ValueError(metric)


def get_winrate_data(exp_num, scenarios, exp_name=EXP_NAME, cluster_threshold=5,
                     min_days=100):
    """Per seed and metric: did a long-origin clade win, and what share of the
    eligible clades were long-origin (the null)? Cached."""
    cache_dir = os.path.join(dm.get_data_dir(), "figures", ".cache")
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(
        cache_dir,
        f"winrate_{exp_name}_{exp_num}_{cluster_threshold}"
        f"_w{BURNIN_CUTOFF_DAYS}_m{MIN_FINAL_TIME}.pkl")
    if os.path.exists(cache_file):
        try:
            df = pd.read_pickle(cache_file)
            print(f"[fig3] win rate: {len(df)} rows from cache")
            return df
        except Exception:
            pass

    rows = []
    for scenario in scenarios:
        sod = _experiment_sod(exp_num, scenario, exp_name=exp_name)
        if sod is None:
            continue
        print(f"[fig3] win rate: {scenario}...", flush=True)
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            if not _usable(ssod):
                continue
            try:
                F, totals, labels = _windowed_clade_data(ssod, cluster_threshold)
            except Exception:
                continue
            if F is None or F.empty or not len(F.columns):
                continue
            seed = dm.get_seed_from_SSOD(ssod)
            for metric in WIN_METRICS:
                winner, pool = _winner_and_pool(F, totals, metric)
                if winner is None or not pool:
                    # No eligible clade: the seed cannot vote on this metric.
                    # _resolve_label would silently return "standard" here.
                    continue
                n_long = sum(1 for c in pool
                             if _resolve_label(labels, c) == "long")
                rows.append({
                    "scenario": scenario, "seed": seed, "metric": metric,
                    "won_long": int(_resolve_label(labels, winner) == "long"),
                    "pool": len(pool),
                    "null": n_long / len(pool),
                })

    df = pd.DataFrame(rows)
    try:
        df.to_pickle(cache_file)
    except Exception:
        pass
    return df


def winrate_summary(df):
    """Observed win rate, its Wilson 95% interval, and the mean null, per
    scenario and metric."""
    out = []
    for (scenario, metric), g in df.groupby(["scenario", "metric"]):
        n = len(g)
        if not n:
            continue
        k = int(g["won_long"].sum())
        phat = k / n
        z = 1.959963985
        denom = 1 + z * z / n
        centre = (phat + z * z / (2 * n)) / denom
        half = z * np.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
        out.append({"scenario": scenario, "metric": metric, "n": n,
                    "win": phat, "lo": max(0.0, centre - half),
                    "hi": min(1.0, centre + half),
                    "null": float(g["null"].mean()),
                    "pool": float(g["pool"].mean())})
    return pd.DataFrame(out)
