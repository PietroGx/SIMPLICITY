#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 3 for long shedders (cached, with progress, PNG output).

Panels:
  A: SOD1 lineage frequency (raw)
  B: SOD1 clustered frequency (shared_mut_threshold = --cluster-threshold)
  C: SOD2 lineage frequency (raw)
  D: SOD2 clustered frequency (shared_mut_threshold = --cluster-threshold)
  E: 5 violins of GENETIC distance per TRANSMISSION EVENT:
     - 7.5d "standard" from BASELINE SOD of --exp-name (SOD with long_shedders_ratio == 0)
     - one violin per tau_3_long from INDEX filtered by evo=5, R=3, ratio=0.01
       (first experiment #1, first SOD; we aggregate ALL seeds in that SOD)

CLI:
  --cluster-threshold  int
  --seed               int
  --exp-name           str (default: generate_data_standard_vs_long_#1)
  --freq-threshold     float (default: 0.01)  # min peak freq to show a lineage in A/C
  --debug              flag to print diagnostics

Notes:
- Caching is automatic: per-SSOD CSVs under Data/cache_hamming/, plus an aggregated CSV under Data/.
- Colors for A-D use your lineage colormap; E uses baseline-first cycle order.
- ASCII-only strings; no grids.
"""

import os
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib import colormaps
from matplotlib.legend_handler import HandlerBase
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from tqdm import tqdm
import seaborn as sns

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.plots_manager as pm
import simplicity.clustering as cl
from simplicity.settings_manager import get_parameter_value_from_simulation_output_dir as get_param
from simplicity.phenotype.distance import hamming_iw

pm.apply_plos_rcparams()

INDEX_CSV = "scripts/long_shedders_experiments/INDEX_exp_scripts.csv"
BASELINE_TAU_LABEL = "7.5d"
CACHE_DIR_NAME = "cache_hamming"  # under Data/
AGG_CSV_NAME = "figure3_long_distances_{exp_name}.csv"

# ---------- Legend ----------

class RainbowLegendBoxHandle:
    """Dummy legend handle for a rainbow gradient box."""
    pass

class HandlerRainbowLegendBox(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, w, h, fs, trans):
        n = 30
        ys = np.linspace(y0, y0 + h, n + 1)
        segs = []
        for i in range(n):
            y = (ys[i] + ys[i + 1]) / 2.0
            segs.append([(x0, y), (x0 + w, y)])
        lc = LineCollection(
            segs,
            cmap=colormaps["gist_rainbow"],
            norm=Normalize(0, 1),
            linewidths=h / n,
            transform=trans,
            zorder=2,
        )
        lc.set_array(np.linspace(0, 1, n))
        border = Rectangle((x0, y0), w, h, transform=trans,
                           edgecolor="black", facecolor="none",
                           linewidth=0.5, zorder=3)
        return [lc, border]

def add_custom_legend_right(fig, violin_colors, violin_labels, anchor=(0.97, 0.97)):

    handles, labels = [], []

    # Rainbow box for the frequency panels A–D
    handles.append(RainbowLegendBoxHandle())
    labels.append("Lineages / Clusters")

    # spacer
    handles.append(Line2D([], [], linestyle="None"))
    labels.append(" ")

    # header + swatches for violin groups
    handles.append(Line2D([], [], linestyle="None"))
    labels.append("Groups (violin)")
    for lab in violin_labels:
        c = violin_colors.get(lab, "0.5")
        handles.append(Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=c))
        labels.append(lab)

    fig.legend(
        handles, labels,
        loc="upper left",
        bbox_to_anchor=anchor,
        frameon=False,
        ncol=1,
        fontsize="small",
        handlelength=1.2,
        handletextpad=0.6,
        labelspacing=0.5,
        borderaxespad=0.0,
        handler_map={RainbowLegendBoxHandle: HandlerRainbowLegendBox()},
    )

    fig.subplots_adjust(right=0.95)

def add_violin_legend_on_ax(ax, labels, color_map, anchor=(0.1, 0.98)):
    """
    Put the violin legend inside ax (axE).
    anchor is (x, y) in axes coords. Lower y to move down; lower x to move left.
    """
    handles = [Rectangle((0, 0), 1, 1, facecolor=color_map[lab], edgecolor=color_map[lab])
               for lab in labels]
    ax.legend(
        handles, labels,
        loc="upper right",
        bbox_to_anchor=anchor,
        frameon=False,
        ncol=1,
        fontsize="small",
        handlelength=1.2,
        handletextpad=0.6,
        labelspacing=0.5,
        borderaxespad=0.0,
    )

# ---------- Utilities ----------

def _is_number(x):
    try:
        float(x)
        return True
    except Exception:
        return False

def format_tau_label(tau):
    x = float(tau)
    if abs(x - int(x)) < 1e-9:
        return "{}d".format(int(x))
    return "{}d".format(x)

def get_sod1_sod2(exp_name, debug=False):
    sods = dm.get_simulation_output_dirs(exp_name)
    if len(sods) < 2:
        raise ValueError("Need at least two simulation_output dirs in experiment '{}'".format(exp_name))
    if debug:
        print("[DEBUG] {} SODs found for exp '{}':".format(len(sods), exp_name))
        for s in sods[:5]:
            print("  [DEBUG] SOD:", s)
    return sods[0], sods[1]

def get_exp_standard_sod(exp_name, tol=1e-12, debug=False):
    """
    Return the SOD from this experiment whose 'long_shedders_ratio' is zero
    (within tolerance). If none equals zero, return the one with the smallest ratio.
    """
    sods = dm.get_simulation_output_dirs(exp_name)
    if not sods:
        raise ValueError("No simulation_output dirs for experiment '{}'".format(exp_name))

    pairs = []
    for sod in sods:
        try:
            r = get_param(sod, "long_shedders_ratio")
            r = float(r) if r is not None and _is_number(r) else float("inf")
        except Exception:
            r = float("inf")
        pairs.append((sod, r))
        if debug:
            print("[DEBUG] SOD ratio long_shedders_ratio =", r, "->", sod)

    for sod, r in pairs:
        if r <= tol:
            if debug:
                print("[DEBUG] Using baseline SOD (ratio==0):", sod)
            return sod

    sod_min, rmin = min(pairs, key=lambda kv: kv[1])
    warnings.warn("No SOD with long_shedders_ratio == 0; using smallest ratio SOD.")
    if debug:
        print("[DEBUG] Using baseline SOD (min ratio {}): {}".format(rmin, sod_min))
    return sod_min

def iter_all_ssods_in_sod(sod_path):
    """Return absolute paths to all SSODs (seed_* folders) under a SOD."""
    ssods = []
    try:
        for name in sorted(os.listdir(sod_path)):
            p = os.path.join(sod_path, name)
            if os.path.isdir(p) and name.startswith("seed_"):
                ssods.append(p)
    except Exception:
        pass
    return ssods

# ---------- GENETIC distance per transmission event ----------

def compute_transmission_hamming_distances(ssod, donor_type_filter=None):
    """
    One Hamming distance per transmission event:
      dist = hamming_iw(G_inherited_by_donor, G_transmitted_by_donor)
    Restrict to donors where row['type'] == donor_type_filter if provided.
    """
    df = om.read_individuals_data(ssod)
    phylo = om.read_phylogenetic_data(ssod)
    lin2gen = dict(zip(phylo["Lineage_name"], phylo["Genome"]))

    dists = []
    for _, row in df.iterrows():
        if donor_type_filter and str(row.get("type", "")) != donor_type_filter:
            continue
        inherited = row.get("inherited_lineage", None)
        if inherited is None:
            continue
        gin = lin2gen.get(inherited, None)
        if gin is None:
            continue

        events = row.get("new_infections", [])
        if not events:
            continue

        for ev in events:
            if not isinstance(ev, dict):
                continue
            tlin = ev.get("transmitted_lineage", None)
            if tlin is None:
                continue
            gout = lin2gen.get(tlin, None)
            if gout is None:
                continue
            try:
                d = hamming_iw(gin, gout)
            except Exception:
                continue
            if np.isfinite(d):
                dists.append(d)
    return dists

# ---------- Caching (automatic) ----------

def cache_dir():
    return os.path.join(dm.get_data_dir(), CACHE_DIR_NAME)

def agg_csv_path(exp_name):
    return os.path.join(dm.get_data_dir(), AGG_CSV_NAME.format(exp_name=exp_name))

def ssod_cache_filename(group_label, donor_type, ssod):
    sodname = os.path.basename(os.path.dirname(ssod))
    seedname = os.path.basename(ssod)
    fn = "{}__{}__{}__{}.csv".format(group_label, donor_type, sodname, seedname)
    return os.path.join(cache_dir(), fn)

def ensure_cached_ssod(group_label, donor_type, ssod, debug=False):
    """
    Ensure per-SSOD CSV exists; compute if missing.
    Returns path to the cache CSV.
    """
    os.makedirs(cache_dir(), exist_ok=True)
    outp = ssod_cache_filename(group_label, donor_type, ssod)
    if os.path.isfile(outp) and os.path.getsize(outp) > 0:
        return outp

    dists = compute_transmission_hamming_distances(ssod, donor_type_filter=donor_type)
    pd.DataFrame({"dist": dists}).to_csv(outp, index=False)
    if debug:
        print("[DEBUG] Wrote {} distances -> {}".format(len(dists), outp))
    return outp

def collect_groups(exp_name, seed, debug=False):
    """
    Build the list of groups for Panel E, each with its SOD and donor type.
    Returns:
      baseline_sod, groups_list
    where groups_list = [(label, donor_type, sod_path), ...] with taus sorted.
    """
    base_sod = get_exp_standard_sod(exp_name, debug=debug)
    groups = [(BASELINE_TAU_LABEL, "standard", base_sod)]

    df = read_index_csv(INDEX_CSV)
    evo = pd.to_numeric(df.get("long_evo_rate_f"), errors="coerce")
    rlg = pd.to_numeric(df.get("R_long"), errors="coerce")
    ratio = pd.to_numeric(df.get("long_shedders_ratio"), errors="coerce")
    sub = df[(evo == 5) & (rlg == 3) & (ratio == 0.01)]
    if sub.empty:
        warnings.warn("INDEX filter returned empty set for evo=5, R=3, ratio=0.01.")
        return base_sod, groups

    taus = sorted(pd.to_numeric(sub["tau_3_long"], errors="coerce").dropna().unique().tolist())
    for t in taus:
        row = sub[pd.to_numeric(sub["tau_3_long"], errors="coerce") == t].iloc[0]
        base = str(row["generated_experiment_name"])
        exp2 = base + "_#1"
        try:
            sods = dm.get_simulation_output_dirs(exp2)
        except Exception as e:
            if debug:
                print("[DEBUG] Missing experiment for tau {}d: {} ({})".format(t, exp2, e))
            continue
        if not sods:
            if debug:
                print("[DEBUG] No SODs for tau {}d in experiment {}".format(t, exp2))
            continue
        lab = format_tau_label(t)
        groups.append((lab, "long_shedder", sods[0]))
        if debug:
            print("[DEBUG] Tau {}d -> exp {}, SOD {}".format(t, exp2, sods[0]))
    return base_sod, groups

def read_index_csv(path):
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="cp1252")

def build_violin_df_with_cache(exp_name, seed, debug=False):
    """
    Ensure per-SSOD caches for all required groups, aggregate to a single DF,
    and save an aggregated CSV under Data/ for future instant loads.
    """
    aggp = agg_csv_path(exp_name)
    if os.path.isfile(aggp) and os.path.getsize(aggp) > 0:
        try:
            df = pd.read_csv(aggp)
            labels = [BASELINE_TAU_LABEL] + sorted(
                [g for g in df["group"].unique() if g != BASELINE_TAU_LABEL],
                key=lambda s: float(s[:-1]) if s.endswith("d") and _is_number(s[:-1]) else 1e9
            )
            return df, labels
        except Exception:
            pass

    baseline_sod, groups = collect_groups(exp_name, seed, debug=debug)

    rows = []
    for lab, dtype, sod in groups:
        ssods = iter_all_ssods_in_sod(sod)
        if debug:
            print("[DEBUG] Group {}, donor_type {}, SOD {}, {} SSODs".format(lab, dtype, sod, len(ssods)))
        for ssod in tqdm(ssods, desc="{} ({})".format(lab, dtype), unit="ssod"):
            cache_csv = ensure_cached_ssod(lab, dtype, ssod, debug=debug)
            try:
                ddf = pd.read_csv(cache_csv)
            except Exception:
                continue
            if "dist" not in ddf.columns:
                continue
            for v in ddf["dist"].dropna().values:
                rows.append({"group": lab, "dist": v, "dtype": dtype, "src_ssod": os.path.basename(ssod)})

    if not rows:
        warnings.warn("No data for Panel E violins.")
        return pd.DataFrame(columns=["group", "dist", "dtype", "src_ssod"]), []

    df = pd.DataFrame(rows)

    try:
        df.to_csv(aggp, index=False)
        if debug:
            print("[DEBUG] Wrote aggregated CSV:", aggp, "rows:", len(df))
    except Exception as e:
        if debug:
            print("[DEBUG] Failed writing aggregated CSV {}: {}".format(aggp, e))

    labels = [BASELINE_TAU_LABEL] + sorted(
        [g for g in df["group"].unique() if g != BASELINE_TAU_LABEL],
        key=lambda s: float(s[:-1]) if s.endswith("d") and _is_number(s[:-1]) else 1e9
    )
    return df, labels

# ---------- Plot helpers for panels A-D ----------

def plot_lineage_frequency_ax(ax, ssod, label_char, freq_threshold):
    lf = om.read_lineage_frequency(ssod)
    t_final = om.read_final_time(ssod)
    cmap_df = pm.make_lineages_colormap(ssod)

    pivot = lf.pivot(index="Time_sampling",
                     columns="Lineage_name",
                     values="Frequency_at_t")

    # keep only columns that ever have nonzero freq
    cols = [c for c in pivot.columns if (pivot[c].fillna(0) > 0).any()]
    pivot = pivot[cols]

    # drop lineages whose peak freq is below the threshold
    if freq_threshold is not None and freq_threshold > 0:
        max_by_col = pivot.fillna(0).max(axis=0)
        keep_cols = [c for c in pivot.columns if max_by_col.get(c, 0.0) >= freq_threshold]
        pivot = pivot[keep_cols]

    # colors aligned to remaining columns
    colors = [pm.get_lineage_color(c, cmap_df) for c in pivot.columns]

    # plot as non-stacked area (overlaid fills)
    pivot.plot(kind="area", stacked=False, color=colors, alpha=0.6, ax=ax, legend=False)

    ax.set_xlim(0, t_final)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Lineage freq.")
    ax.set_xlabel("Time (d)")
    ax.set_title(label_char, loc="left", pad=8, fontsize=16, fontweight="bold")
    pm.apply_standard_axis_style(ax)
    return cmap_df, t_final

def plot_clustered_frequency_ax(ax, ssod, cluster_threshold, colormap_df, t_final, label_char):
    phylo_df = om.read_phylogenetic_data(ssod)
    lf = om.read_lineage_frequency(ssod)
    full_pivot = lf.pivot(index="Time_sampling",
                          columns="Lineage_name",
                          values="Frequency_at_t")

    (clade_to_lineages,
     lineage_to_clade,
     per_clade_mut_df,
     clade_meta_df) = cl.cluster_lin_into_clades_with_meta(
        phylo_df,
        shared_mut_threshold=cluster_threshold
    )

    clade_series = []
    for clade, members in clade_to_lineages.items():
        cols = [c for c in members if c in full_pivot.columns]
        if not cols:
            continue
        s = full_pivot[cols].sum(axis=1)
        s.name = clade
        clade_series.append(s)

    if not clade_series:
        raise ValueError("No clade members present in frequency table.")

    clade_freq = pd.concat(clade_series, axis=1)

    try:
        order = (clade_meta_df
                 .set_index("clade")
                 .loc[clade_freq.columns, "start_time"]
                 .sort_values(kind="mergesort")
                 .index)
        clade_freq = clade_freq[order]
    except Exception:
        pass

    root_by_clade = clade_meta_df.set_index("clade")["root_lineage"].to_dict()
    colors = []
    for clade in clade_freq.columns:
        root_lin = root_by_clade.get(clade, None)
        try:
            colors.append(pm.get_lineage_color(root_lin, colormap_df))
        except Exception:
            cmap = colormaps["tab10"]
            idx = list(clade_freq.columns).index(clade) % 10
            colors.append(cmap(idx))

    clade_freq.plot(kind="area", stacked=True, color=colors, alpha=0.6, ax=ax, legend=False)

    ax.set_xlim(0, t_final)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Clade freq.")
    ax.set_xlabel("Time (d)")
    ax.set_title(label_char, loc="left", pad=8, fontsize=16, fontweight="bold")
    pm.apply_standard_axis_style(ax)

# ---------- Panel E (violins; from cached data) ----------

def plot_violins_panel_E(ax, df, labels):
    if df.empty or not labels:
        ax.axis("off")
        return {}, []

    color_map = build_tau_violin_colors(labels)
    palette = {lab: color_map[lab] for lab in labels}

    sns.violinplot(
        data=df, x="group", y="dist",
        order=labels, hue="group", palette=palette,
        cut=0, inner="quartile", linewidth=1,
        ax=ax, legend=False
    )

    ax.set_xlabel("")
    ax.set_ylabel("Hamming distance")
    ax.set_title("E", loc="left", pad=8, fontsize=16, fontweight="bold")
    pm.apply_standard_axis_style(ax)

    return color_map, labels

def build_tau_violin_colors(labels):
    """
    Baseline first color, then each tau label gets the next colors in the Matplotlib prop cycle.
    """
    cycle = plt.rcParams.get("axes.prop_cycle", None)
    base = cycle.by_key()["color"] if cycle is not None else [colormaps["tab10"](i) for i in range(10)]
    ordered = []
    if BASELINE_TAU_LABEL in labels:
        ordered.append(BASELINE_TAU_LABEL)
    ordered += [lab for lab in labels if lab != BASELINE_TAU_LABEL]
    color_map = {}
    for i, lab in enumerate(ordered):
        color_map[lab] = base[i % len(base)]
    return color_map

# ---------- Build figure ----------

def build_figure(exp_name, cluster_threshold, seed, freq_threshold, debug=False):
    # A-D inputs
    sod1, sod2 = get_sod1_sod2(exp_name, debug=debug)
    ssod1 = dm.get_ssod(sod1, seed)
    ssod2 = dm.get_ssod(sod2, seed)
    if debug:
        print("[DEBUG] SOD1 ->", sod1, "SSOD1 ->", ssod1)
        print("[DEBUG] SOD2 ->", sod2, "SSOD2 ->", ssod2)

    # E data via automatic cache (per-SSOD + aggregated CSV)
    dfE, labelsE = build_violin_df_with_cache(exp_name, seed, debug=debug)

    # figure
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1, 1, 1.1], hspace=0.35, wspace=0.25)
    fig.subplots_adjust(top=0.92)

    axA = fig.add_subplot(gs[0, 0])
    cmap_df1, tfin1 = plot_lineage_frequency_ax(axA, ssod1, label_char="A", freq_threshold=freq_threshold)

    axB = fig.add_subplot(gs[1, 0])
    plot_clustered_frequency_ax(axB, ssod1, cluster_threshold, cmap_df1, tfin1, label_char="B")

    axC = fig.add_subplot(gs[0, 1])
    cmap_df2, tfin2 = plot_lineage_frequency_ax(axC, ssod2, label_char="C", freq_threshold=freq_threshold)

    axD = fig.add_subplot(gs[1, 1])
    plot_clustered_frequency_ax(axD, ssod2, cluster_threshold, cmap_df2, tfin2, label_char="D")

    axE = fig.add_subplot(gs[2, :])
    violin_colors, violin_labels = plot_violins_panel_E(axE, dfE, labelsE)

    # Legend inside violin subplot (upper-left or by your chosen anchor)
    # add_violin_legend_on_ax(axE, violin_labels, violin_colors, anchor=(0.1, 0.98))

    # save to Data/ as PNG
    data_dir = dm.get_data_dir()
    os.makedirs(data_dir, exist_ok=True)
    outfile = os.path.join(
        data_dir, "Figure_3_long_{}_thr{}_seed{}.png".format(exp_name, cluster_threshold, seed)
    )
    print("[SAVE]", outfile)
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)

# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(description="Build Figure 3 (long shedders) with caching and progress")
    ap.add_argument("--cluster-threshold", type=int, default=5, help="Clustering threshold (default 5)")
    ap.add_argument("--seed", type=int, default=1, help="Seed number for A-D (default 1)")
    ap.add_argument("--exp-name", type=str, default="generate_data_standard_vs_long_#1",
                    help="Experiment name (default generate_data_standard_vs_long_#1)")
    ap.add_argument("--freq-threshold", type=float, default=0.05,
                    help="Min peak lineage freq to display in A/C (default 0.01)")
    ap.add_argument("--debug", action="store_true", help="Print debug info (experiments/SODs/SSODs/groups)")
    args = ap.parse_args()

    build_figure(args.exp_name, args.cluster_threshold, args.seed, args.freq_threshold, debug=args.debug)

if __name__ == "__main__":
    main()