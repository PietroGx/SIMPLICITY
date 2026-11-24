#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Clade metrics for figures 2, 3 and 4: one clustering pass, every per-clade
number derived from it.

Moved here from scripts/long_paper_figures/long_shedders_preprocess.py in
v2.5.6, when that directory was archived. It was never long-paper-specific --
`get_clade_metrics` and `summarize_sod_pies` are what fig2_preprocess_data
reads, `get_lineage_labels` is what fig4_preprocess_data reads, and the three
underscore helpers are the metric definitions themselves. The figures that use
them live here, so they live here.

`_clade_analysis` is the single clustering pass: everything else derives from
its return value, so the four metrics cannot be computed over four different
partitions.

A clade's origin label comes from `simplicity.clustering.resolve_clade_label`,
which is the only copy of that rule now -- the private `_resolve_label` this
module used to carry was a second statement of it.
"""
import os

import numpy as np
import pandas as pd
from collections import Counter

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.clustering as cl
from scipy.signal import find_peaks

# Skip accounting, from the sibling that already owns it -- every importer of
# this module imports _scenarios too, directly or through
# fig2_preprocess_data, so the import direction costs nothing.
from _scenarios import (EXPECTED_SEED_SKIPS, record_seed_skip,
                        format_seed_skips)



# ---------------------------------------------------------------- episodes
# One definition of an episode and of its rise, shared by the figure-2 builder
# and the studio. They were separate before v2.5.7: the builder had `_rise`
# (the last approach to the global peak) and `_risef` (the first local maximum
# clearing the upper crossing), and the studio measured the episode containing
# the clade's own maximum -- three answers to "how fast did it rise" shipping
# together. `find_peaks` is scipy's, so the prominence test is not reimplemented.

def smooth(y, t, win_days):
    """Centred moving average over a window given in DAYS. Mirrors rise_lab's
    `smooth`: the window is converted to samples with the mean step, forced
    odd, and the edges shrink rather than wrap."""
    n = len(y)
    if n < 2:
        return np.asarray(y, float)
    step = (t[-1] - t[0]) / (n - 1)
    w = max(1, int(round(win_days / max(step, 1e-9))))
    if w % 2 == 0:
        w += 1
    off = w >> 1
    pre_ = np.concatenate([[0.0], np.cumsum(y, dtype=float)])
    a = np.maximum(0, np.arange(n) - off)
    b = np.minimum(n - 1, np.arange(n) + off)
    return (pre_[b + 1] - pre_[a]) / (b - a + 1)


TROUGH_PROM = 1e-4   # an episode start must be this prominent a dip; see episodes_of
CENSOR_TOL = 1e-9    # relative; "at the maximum" for a censored peak, see episodes_of


def episodes_of(ys, prom):
    """(peak height, [(start index, peak index, censored)]). An episode's start
    is the local minimum DIRECTLY before its peak -- not the prominence base,
    which tunnels left to the start of the record on most clades.

    CENSORED: the record ends before the clade comes back down. find_peaks
    never returns a plateau touching the last sample, and a clade that holds
    high to the end has no descent to give it prominence, so a clade that swept
    to fixation had no episode at all. So: the clade's maximum, if it rose to
    it by the prominence bar and never falls back by that bar before the
    record ends, is one more episode, flagged censored. It peaks at the first
    sample within CENSOR_TOL of the maximum -- a smoothed plateau at 1.0
    carries ~1e-16 noise, and a peak finder sees a micro-peak in every wiggle.
    One reaching its maximum on the last sample has no observed peak and
    stays out.

    A trough counts only if it clears TROUGH_PROM of the clade's peak. Without
    that test every shallow wiggle on a long flat approach was a candidate
    start, so which one happened to be LAST could turn on noise far below the
    data's own resolution: rounding a trajectory to five decimals moved one
    clade's rise from 35.02 to 106.06 days, because the rounding erased a
    wiggle at index 600 and promoted the real trough at 496.

    TROUGH_PROM is deliberately tiny, and it was chosen by sweeping it over 150
    clades in 15 runs at k=7, m_T=1% (v2.5.7). Larger bars do remove the
    instability but drag the start leftward, which is the tunnelling the
    prominence base was rejected for:

        bar       median   mean   clades changed   bistable
        none       34.01   32.53        0 of 150      1
        1e-4*h     34.01   33.36        4             0
        1e-3*h     34.99   35.57       14             0
        1e-2*h     36.00   39.09       36             0
        1e-1*h     41.00   52.65       61             0

    1e-4 is the smallest that removes it: the median is unchanged, four clades
    move, and one of those four IS the unstable clade now choosing the deeper
    trough. It sits above the ~1e-5 noise and below the ~1e-4 granularity of a
    frequency, so it discriminates wiggles from dips rather than reshaping the
    metric.
    """
    ys = np.asarray(ys, float)
    n = len(ys)
    h = float(np.max(ys)) if n else 0.0
    if not (h > 0):
        return 0.0, []
    pk, _ = find_peaks(ys, prominence=prom * h)
    peaks = [(int(i), False) for i in pk]
    m = int(np.argmax(ys >= h * (1 - CENSOR_TOL)))
    if (m < n - 1 and h - ys[:m + 1].min() >= prom * h
            and h - ys[m:].min() < prom * h):
        peaks = sorted(peaks + [(m, True)])
    if not peaks:
        return h, []
    mins, _ = find_peaks(-ys, prominence=TROUGH_PROM * h)
    mins = np.concatenate([[0], mins])
    out = []
    for i, c in peaks:
        before = mins[mins < i]
        out.append((int(before[-1]) if len(before) else 0, i, c))
    return h, out


def rise_of(ys, t, a, b, lo, hi):
    """Days from lo to hi of THIS episode's own amplitude. None if not rising."""
    ylo, yhi = float(ys[a]), float(ys[b])
    if not (yhi > ylo):
        return None
    A = ylo + lo * (yhi - ylo)
    B = ylo + hi * (yhi - ylo)
    seg = ys[a:b + 1]
    ia = a + int(np.argmax(seg >= A)) if (seg >= A).any() else a
    ib = a + int(np.argmax(seg >= B)) if (seg >= B).any() else b
    return dict(a=a, b=b, ia=ia, ib=ib,
                days=float(t[ib] - t[ia]), amp=yhi - ylo, top=yhi)


# ---------------------------------------------------------------- per run


def _peak_values(F_clade):
    """{clade: max frequency reached}."""
    if F_clade.empty or F_clade.shape[1] == 0: return {}
    return F_clade.max(axis=0).to_dict()


def _survival_durations(F_clade, eps=1e-4):
    """{clade: duration spent above eps frequency}."""
    if F_clade.empty or F_clade.shape[1] == 0: return {}
    times = F_clade.index.to_numpy(dtype=float)
    durations = {}
    for clade in F_clade.columns:
        mask = F_clade[clade].to_numpy(dtype=float) > eps
        if not mask.any(): durations[clade] = 0.0
        else: durations[clade] = float(times[len(mask) - 1 - mask[::-1].argmax()] - times[mask.argmax()])
    return durations


def _growth_crossing_times(F_clade, lower=0.01, upper=0.50):
    """{clade: time to cross lower->upper frequency} (only clades that cross both)."""
    if F_clade.empty or F_clade.shape[1] == 0: return {}
    times = F_clade.index.to_numpy(dtype=float)
    candidates = {}
    for clade in F_clade.columns:
        y = F_clade[clade].to_numpy(dtype=float)
        # Simplified crossing logic for speed
        over_lower = times[y >= lower]
        over_upper = times[y >= upper]
        if len(over_lower) > 0 and len(over_upper) > 0 and over_upper[0] > over_lower[0]:
            candidates[clade] = over_upper[0] - over_lower[0]
    return candidates


def _clade_analysis(ssod, cluster_threshold):
    """
    Shared clustering + per-clade metric computation for one SSOD. Returns
    (F_nr, clade_to_lineages, totals, clade_labels) so get_clade_winners/
    get_clade_metrics/get_lineage_labels can all derive their result from a
    single clustering pass.
    """
    phylo_df = om.read_phylogenetic_data(ssod)
    lf = om.read_lineage_frequency(ssod)

    # Rule B, exactly as simplicity.clustering.get_clustering_for_ssod applies
    # it: a lineage founds a clade when it has diverged >= cluster_threshold
    # mutations from its clade root AND its subtree carries at least
    # cl.MIN_TX_FRAC of this run's own transmission events. A refused lineage is
    # FOLDED into its surviving ancestor, not dropped.
    # Divergence alone makes a clade of almost every tip -- measured on control
    # at k=5, 402 clades per run of which 86% hold one lineage and 78% infected
    # nobody, against 14 under rule B. Every number below derives from this
    # pass, so this one call is what figures 2 and 4 count.
    individuals_df = om.read_individuals_data(ssod)
    support = cl.subtree_transmission_counts(phylo_df, individuals_df)
    min_support = cl.transmission_floor(individuals_df, cl.MIN_TX_FRAC)

    clade_to_lineages, _, per_clade_mut_df, clade_meta_df = cl.cluster_lin_into_clades_with_meta(
        phylo_df, shared_mut_threshold=cluster_threshold,
        support=support, min_support=min_support
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

    # Burden per clade
    lin_to_total = dict(zip(phylo_df["Lineage_name"], phylo_df.get("Total_infections", np.zeros(len(phylo_df)))))
    totals = {clade: sum(int(lin_to_total.get(l, 0)) for l in members) for clade, members in clade_to_lineages.items()}

    return F_nr, clade_to_lineages, totals, clade_labels


def get_clade_winners(ssod, cluster_threshold):
    """Evaluates peak, burden, survival, and growth winners for a single SSOD."""
    return get_clade_metrics(ssod, cluster_threshold)["winners"]


def get_clade_metrics(ssod, cluster_threshold):
    """
    Like get_clade_winners, but also returns the continuous value achieved by
    long-origin clades for each of the 4 metrics -- 0.0 if no long-origin
    clade exists in this seed (e.g. every 'control' seed: a legitimate "no
    long-shedder impact" data point, not an error).

    Returns {"winners": {Peak, Burden, Survival, Growth -> category label},
             "values":  {Peak, Burden, Survival, Growth -> float}}
    """
    F_nr, _clade_to_lineages, totals, clade_labels = _clade_analysis(ssod, cluster_threshold)

    peak_values = _peak_values(F_nr)
    survival_durations = _survival_durations(F_nr)
    growth_crossings = _growth_crossing_times(F_nr)

    peak_winner = max(peak_values.items(), key=lambda kv: kv[1])[0] if peak_values else None
    survival_winner = max(survival_durations.items(), key=lambda kv: kv[1])[0] if survival_durations else None
    growth_winner = min(growth_crossings.items(), key=lambda kv: kv[1])[0] if growth_crossings else None
    burden_winner = max(totals.items(), key=lambda kv: kv[1])[0] if totals else None

    winners = {
        "Peak": cl.resolve_clade_label(clade_labels, peak_winner),
        "Burden": cl.resolve_clade_label(clade_labels, burden_winner),
        "Survival": cl.resolve_clade_label(clade_labels, survival_winner),
        "Growth": cl.resolve_clade_label(clade_labels, growth_winner),
    }

    def is_long(clade):
        return cl.resolve_clade_label(clade_labels, clade) == "long"

    long_peak = [v for c, v in peak_values.items() if is_long(c)]
    long_survival = [v for c, v in survival_durations.items() if is_long(c)]
    long_growth = [v for c, v in growth_crossings.items() if is_long(c)]
    long_burden = [v for c, v in totals.items() if is_long(c)]

    values = {
        # Peak/Survival: strongest single long-origin clade.
        "Peak": max(long_peak) if long_peak else 0.0,
        "Survival": max(long_survival) if long_survival else 0.0,
        # Burden is additive (infection counts), so sum across long clades.
        "Burden": float(sum(long_burden)) if long_burden else 0.0,
        # Growth: crossing TIME is smaller-is-faster; invert to a rate so
        # "bigger value = more long-shedder impact" holds for all 4 metrics,
        # and "no long clade" consistently means 0.0, not a misleadingly
        # "instant" crossing time.
        "Growth": (1.0 / min(long_growth)) if long_growth else 0.0,
    }

    return {"winners": winners, "values": values}


def get_lineage_labels(ssod, cluster_threshold):
    """
    {lineage_name: normalized_label} for every lineage in this SSOD's tree
    (long/standard/mixed/founder) -- the per-lineage counterpart to
    get_clade_metrics's clade-level aggregates. Used to color Figure 4's
    tree comparison with the same color language as Figure 3's pies/PCA.
    """
    _F_nr, clade_to_lineages, _totals, clade_labels = _clade_analysis(ssod, cluster_threshold)
    return {
        lineage: cl.resolve_clade_label(clade_labels, clade)
        for clade, lineages in clade_to_lineages.items()
        for lineage in lineages
    }


def summarize_sod_pies(sod, cluster_threshold, min_days):
    """Aggregates winners across all seeds in an SOD.

    `valid_seeds` is the pie denominator for figures 2 and 4, so it states
    itself: N of M, with every seed that did not make it and why. Only an
    absent/unreadable or empty/truncated output is a skip; anything structural
    propagates. The accumulation and the increment sit OUTSIDE the try --
    covering them dropped a failing seed from numerator and denominator at
    once, leaving a smaller pie with no indication.
    """
    results = {"Peak": Counter(), "Burden": Counter(), "Survival": Counter(), "Growth": Counter()}
    valid_seeds, skips, short = 0, {}, 0
    ssods = dm.get_seeded_simulation_output_dirs(sod)
    for ssod in ssods:
        try:
            if om.read_final_time(ssod) < min_days:
                short += 1
                continue
            winners = get_clade_winners(ssod, cluster_threshold)
        except EXPECTED_SEED_SKIPS as exc:
            record_seed_skip(skips, exc)
            continue
        for k, v in winners.items(): results[k][v] += 1
        valid_seeds += 1
    print(f"[summarize_sod_pies] {os.path.basename(os.path.dirname(sod))}: "
          f"{valid_seeds} of {len(ssods)} seeds"
          + (f", {short} under {min_days} days" if short else "")
          + f", {format_seed_skips(skips)}", flush=True)
    return results, valid_seeds
    
# --- Figure 5 Extractors (Epidemiological Stats) ---

