#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ============================================================================
# SOT sanity-check: standard vs long-shedder OSR, two clocks
# ----------------------------------------------------------------------------
# LEFT  (global clock):  distance-from-root vs global sequencing time, both
#        cohorts, from sequencing_data_regression.csv.
# RIGHT (intra-host clock): long shedders only. One point per intra-host
#        lineage in each long-shedder's IH_lineages_trajectory:
#          x = (ih_birth - t_infection) / 365.25            [years]
#          y = hamming_iw(genome(lineage), genome(inherited))
#              / len(dis.reference)                          [subs/site]
#
# CRITICAL: intra-host lineage names are unique ONLY within a seed folder.
# Genome lookups are resolved strictly against the SAME ssod's phylogenetic
# data; only the final numeric (x, y) points are pooled across seeds. Names
# and genomes never cross the ssod boundary.
#
# No core files modified. ih_death is never used. Distance is computed only
# via the codebase's dis.hamming_iw (never lineage-name parsing).
# ============================================================================

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless / cluster-safe
import matplotlib.pyplot as plt

import simplicity.dir_manager as dm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er
import simplicity.phenotype.distance as dis

# Cohort styling (matches fig1 D/E colors).
COHORT_STYLE = {
    "standard":     {"color": "#3b528b", "label": "Standard"},
    "long_shedder": {"color": "#e41a1c", "label": "Long shedder"},
}


def set_nature_rcparams():
    """Nature-style rcparams (inlined from fig1_save to avoid its deps)."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7,
        'axes.labelsize': 7,
        'xtick.labelsize': 7,
        'ytick.labelsize': 7,
        'legend.fontsize': 7,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


# ---------------------------------------------------------------------------
# GLOBAL-CLOCK DATA (left subplot) -- unchanged behavior
# ---------------------------------------------------------------------------
def load_regression_data(exp_name, scenario, exp_num):
    """Concatenate sequencing_data_regression across all seeds of the run."""
    experiment_string = f"{exp_name}_{scenario}_#{exp_num}"
    sods = dm.get_simulation_output_dirs(experiment_string)
    if not sods:
        raise RuntimeError(
            f"No simulation output dirs for '{experiment_string}'.")

    frames = []
    for sod in sods:
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            path = os.path.join(ssod, "sequencing_data_regression.csv")
            if not os.path.isfile(path):
                continue
            try:
                df = om.read_sequencing_data_regression(ssod)
            except Exception:
                continue
            if df is not None and not df.empty:
                frames.append(df)

    if not frames:
        raise RuntimeError(
            f"No sequencing_data_regression.csv data for '{experiment_string}'.")

    combined = pd.concat(frames, ignore_index=True)
    required = {"Sequencing_time", "Distance_from_root", "individual_type"}
    missing = required - set(combined.columns)
    if missing:
        raise KeyError(
            f"sequencing_data_regression.csv missing columns: {missing}.")
    return combined


# ---------------------------------------------------------------------------
# INTRA-HOST-CLOCK DATA (right subplot) -- long shedders only
# ---------------------------------------------------------------------------
def parse_ih_trajectory(raw):
    """
    np-aware parse of IH_lineages_trajectory. Pre-fix data serializes values
    as np.float64(...), which defeats ast.literal_eval; a locked-scope eval
    with np available handles both old and fixed data. Returns {} on failure.
    """
    if isinstance(raw, dict):
        return raw
    if not isinstance(raw, str):
        return {}
    try:
        return eval(raw, {"np": np, "__builtins__": {}})
    except Exception:
        return {}


def build_phylo_genome_map(ssod):
    """
    {Lineage_name: Genome(list)} for ONE ssod. Lineage names are unique only
    within a seed folder, so this map must never be shared across ssods.
    """
    ph = om.read_phylogenetic_data(ssod)
    genome_map = {}
    for _, r in ph.iterrows():
        genome_map[str(r["Lineage_name"])] = r["Genome"]
    return genome_map


def load_intrahost_long_data(exp_name, scenario, exp_num):
    """
    Build (Sequencing_time, Distance_from_root, individual_type) points for
    long shedders on the intra-host clock. Names/genomes are resolved strictly
    within each ssod; only numeric points are pooled across seeds.
    """
    experiment_string = f"{exp_name}_{scenario}_#{exp_num}"
    sods = dm.get_simulation_output_dirs(experiment_string)
    if not sods:
        raise RuntimeError(
            f"No simulation output dirs for '{experiment_string}'.")

    ref_len = len(dis.reference)
    points = []            # only (x, y, type) floats cross the ssod boundary
    n_skipped = 0

    for sod in sods:
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            ind_path = os.path.join(ssod, "individuals_data.csv")
            phylo_path = os.path.join(ssod, "phylogenetic_data.csv")
            if not (os.path.isfile(ind_path) and os.path.isfile(phylo_path)):
                continue

            # Per-ssod genome map -- LOCAL SCOPE ONLY, never shared.
            try:
                phylo_map = build_phylo_genome_map(ssod)
            except Exception:
                continue

            # Raw read: IH_lineages_trajectory needs the np-aware parser, so we
            # do not use the native reader (which would ast.literal_eval it).
            try:
                ind_df = pd.read_csv(ind_path, index_col=0)
            except Exception:
                continue
            if "type" not in ind_df.columns:
                continue
            ind_df = ind_df[ind_df["type"] == "long_shedder"]
            if ind_df.empty:
                continue

            for _, row in ind_df.iterrows():
                inherited = str(row.get("inherited_lineage"))
                t_inf = row.get("t_infection")
                if inherited not in phylo_map or pd.isna(t_inf):
                    n_skipped += 1
                    continue
                inh_genome = phylo_map[inherited]

                traj = parse_ih_trajectory(row.get("IH_lineages_trajectory"))
                if not isinstance(traj, dict):
                    continue

                for lineage, ev in traj.items():
                    lin_name = str(lineage)
                    if lin_name not in phylo_map:
                        n_skipped += 1
                        continue
                    ih_birth = ev.get("ih_birth") if isinstance(ev, dict) else None
                    if ih_birth is None:
                        n_skipped += 1
                        continue

                    lin_genome = phylo_map[lin_name]
                    y = dis.hamming_iw(lin_genome, inh_genome) / ref_len
                    x = (float(ih_birth) - float(t_inf)) / 365.25
                    points.append((x, y, "long_shedder"))

    if not points:
        raise RuntimeError(
            f"No intra-host long-shedder points assembled for "
            f"'{experiment_string}'.")

    if n_skipped:
        print(f"[note] intra-host: skipped {n_skipped} record(s) "
              f"(missing lineage/inherited/birth).")

    return pd.DataFrame(
        points, columns=["Sequencing_time", "Distance_from_root", "individual_type"])


# ---------------------------------------------------------------------------
# SHARED FIT / PLOT HELPERS
# ---------------------------------------------------------------------------
def fit_cohort(df_cohort):
    """Fit tempest regression; return (slope=OSR, r2, x, y)."""
    fitted = er.tempest_regression(df_cohort)
    x = df_cohort["Sequencing_time"].values.reshape(-1, 1)
    y = df_cohort["Distance_from_root"].values
    slope = float(fitted.coef_[0])
    r2 = float(fitted.score(x, y))
    return slope, r2, x.ravel(), y


def plot_cohort(ax, df_cohort, cohort_key, x_max):
    """Scatter + fitted zero-intercept line for one cohort. Returns (slope, r2)."""
    style = COHORT_STYLE[cohort_key]
    slope, r2, x, y = fit_cohort(df_cohort)

    if len(x) > 2000:
        idx = np.random.RandomState(42).choice(len(x), 2000, replace=False)
        xs, ys = x[idx], y[idx]
    else:
        xs, ys = x, y

    ax.scatter(xs, ys, s=8, alpha=0.3, color=style["color"],
               edgecolors="none", zorder=2)
    ax.plot([0, x_max], [0, slope * x_max], color=style["color"],
            linewidth=1.5, zorder=3,
            label=f"{style['label']}: {slope:.5f} s/s/y (R²={r2:.2f})")
    return slope, r2


def draw_target(ax, target, x_max, color, label):
    ax.plot([0, x_max], [0, target * x_max], color=color, linestyle="--",
            linewidth=1.0, alpha=0.6, zorder=1, label=label)


def finalize_axis(ax, x_max, title, xlabel):
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Distance (substitutions/site)")
    ax.set_xlim(left=0, right=x_max * 1.02)
    ax.set_ylim(bottom=0)
    ax.set_title(title, fontsize=7, fontweight="bold")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=6, loc="upper left", frameon=False)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="SOT sanity-check: global-clock (both cohorts) vs "
                    "intra-host-clock (long shedders) OSR regressions.")
    parser.add_argument('--exp-num', type=int, default=1)
    parser.add_argument('--scenario', type=str, default='SOT')
    parser.add_argument('--exp-name', type=str, default='impact_long_shedders')
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--outfile', type=str, default=None)
    args = parser.parse_args()

    set_nature_rcparams()

    fig, (ax_global, ax_intra) = plt.subplots(1, 2, figsize=(11, 4.5))

    # ---- LEFT: global clock, both cohorts ----
    combined = load_regression_data(args.exp_name, args.scenario, args.exp_num)
    xg_max = float(combined["Sequencing_time"].max())
    draw_target(ax_global, args.target_osr_std, xg_max, "#3b528b",
                f"Target std = {args.target_osr_std}")
    draw_target(ax_global, args.target_osr_long, xg_max, "#e41a1c",
                f"Target long = {args.target_osr_long}")
    for cohort_key in ("standard", "long_shedder"):
        sub = combined[combined["individual_type"] == cohort_key]
        if sub.empty:
            print(f"[warn] global: no '{cohort_key}' rows; skipping.")
            continue
        slope, r2 = plot_cohort(ax_global, sub, cohort_key, xg_max)
        print(f"[global/{cohort_key}] OSR={slope:.6f} (R²={r2:.3f}, n={len(sub)})")
    finalize_axis(ax_global, xg_max,
                  "Global clock (distance from root)",
                  "Sequencing time (years)")

    # ---- RIGHT: intra-host clock, long shedders only ----
    intra = load_intrahost_long_data(args.exp_name, args.scenario, args.exp_num)
    xi_max = float(intra["Sequencing_time"].max())
    draw_target(ax_intra, args.target_osr_long, xi_max, "#e41a1c",
                f"Target long = {args.target_osr_long}")
    slope, r2 = plot_cohort(ax_intra, intra, "long_shedder", xi_max)
    print(f"[intra/long_shedder] OSR={slope:.6f} (R²={r2:.3f}, n={len(intra)})")
    finalize_axis(ax_intra, xi_max,
                  "Intra-host clock (distance from infecting strain)",
                  "Time since infection (years)")

    plt.tight_layout()

    if args.outfile:
        out_path = args.outfile
    else:
        experiment_string = f"{args.exp_name}_{args.scenario}_#{args.exp_num}"
        out_path = os.path.join(
            dm.get_experiment_plots_dir(experiment_string),
            f"{experiment_string}_global_vs_intrahost.png")

    plt.savefig(out_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()
