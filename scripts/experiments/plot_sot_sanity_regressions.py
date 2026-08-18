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
# One combined figure, one ROW per long-shedder scenario (--scenarios, default
# SOT/HIV_low/HIV_high/edge_case), 2 columns. Every row shares the same x/y
# limits per column, so scenarios are directly visually comparable.
#
# LEFT column (global clock, primary axis): standard individuals'
#        distance-from-root vs global sequencing time, plus a black "mixed"
#        line fit on ALL individuals pooled together (what unstratified
#        surveillance data would show). Secondary (right-hand) y-axis on the
#        same panel overlays the long-shedder intra-host clock fit (line
#        only, no scatter) against its own target, for a direct side-by-side
#        read of both clocks.
# RIGHT column (intra-host clock): long shedders only, full detail
#        (scatter + fit). One point per intra-host lineage in each
#        long-shedder's IH_lineages_trajectory:
#          x = (ih_birth - t_infection) / 365.25            [years]
#          y = hamming_iw(genome(lineage), genome(inherited))
#              / len(dis.reference)                          [subs/site]
#
# Both loaders apply the same final_time/sequence-count gate (--min-len/
# --min-seq, defaults 100/30) cal_2.py's own calibration fits use, so this
# sanity check validates the same population of seeds calibration was
# actually fit on rather than pooling in short/sparse seeds it excluded.
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
# GLOBAL-CLOCK DATA (left subplot)
# ---------------------------------------------------------------------------
def load_regression_data(exp_name, scenario, exp_num, min_seq=30, min_len=100):
    """Concatenate sequencing_data_regression across all seeds of the run.
    Applies the same final_time/sequence-count gate cal_2.py's own
    calibration fits use, so this sanity check validates the same
    population of seeds the calibration was actually fit on."""
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
                final_time = om.read_final_time(ssod)
                df = om.read_sequencing_data_regression(ssod)
            except Exception:
                continue
            if df is None or df.empty:
                continue
            if final_time < min_len or len(df) < min_seq:
                continue
            frames.append(df)

    if not frames:
        raise RuntimeError(
            f"No sequencing_data_regression.csv data for '{experiment_string}' "
            f"passed the min_len={min_len}/min_seq={min_seq} filter.")

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


def load_intrahost_long_data(exp_name, scenario, exp_num, min_seq=30, min_len=100):
    """
    Build (Sequencing_time, Distance_from_root, individual_type) points for
    long shedders on the intra-host clock. Names/genomes are resolved strictly
    within each ssod; only numeric points are pooled across seeds.

    Applies the same final_time/sequence-count gate as load_regression_data
    (and cal_2.py's own calibration fits) per ssod, before extracting any
    intra-host points from it -- same seed population, different columns.
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

            try:
                final_time = om.read_final_time(ssod)
                seq_data = om.read_sequencing_data_regression(ssod)
            except Exception:
                continue
            if seq_data is None or final_time < min_len or len(seq_data) < min_seq:
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


def plot_cohort(ax, df_cohort, cohort_key, x_max, scatter=True, label_suffix=""):
    """Fitted zero-intercept line for one cohort, with optional scatter.
    Returns (slope, r2)."""
    style = COHORT_STYLE[cohort_key]
    slope, r2, x, y = fit_cohort(df_cohort)

    if scatter:
        if len(x) > 2000:
            idx = np.random.RandomState(42).choice(len(x), 2000, replace=False)
            xs, ys = x[idx], y[idx]
        else:
            xs, ys = x, y
        ax.scatter(xs, ys, s=8, alpha=0.3, color=style["color"],
                   edgecolors="none", zorder=2)

    ax.plot([0, x_max], [0, slope * x_max], color=style["color"],
            linewidth=1.5, zorder=3,
            label=f"{style['label']}{label_suffix}: {slope:.5f} s/s/y (R²={r2:.2f})")
    return slope, r2


def draw_target(ax, target, x_max, color, label):
    ax.plot([0, x_max], [0, target * x_max], color=color, linestyle="--",
            linewidth=1.0, alpha=0.6, zorder=1, label=label)


DEFAULT_LONG_SCENARIOS = ["SOT", "HIV_low", "HIV_high", "edge_case"]


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="SOT sanity-check, one row per long-shedder scenario: "
                    "standard global clock + mixed-population global clock "
                    "(left) vs long-shedder intra-host clock (right, and "
                    "overlaid on the left's secondary axis). All rows share "
                    "the same x/y limits per column for direct comparison.")
    parser.add_argument('--exp-num', type=int, default=1)
    parser.add_argument('--scenarios', type=str,
                        default=",".join(DEFAULT_LONG_SCENARIOS),
                        help="Comma-separated long-shedder scenarios, one row each.")
    parser.add_argument('--exp-name', type=str, default='impact_long_shedders')
    parser.add_argument('--target-osr-std', type=float, default=0.0013)
    parser.add_argument('--target-osr-long', type=float, default=0.00205)
    parser.add_argument('--min-seq', type=int, default=30,
                        help="Same gate cal_2.py's own calibration fits use: "
                            "minimum sequencing_data_regression.csv rows for a "
                            "seed to count.")
    parser.add_argument('--min-len', type=int, default=100,
                        help="Same gate cal_2.py's own calibration fits use: "
                            "minimum final_time for a seed to count.")
    parser.add_argument('--outfile', type=str, default=None)
    args = parser.parse_args()

    scenarios = [s.strip() for s in args.scenarios.split(",") if s.strip()]
    if not scenarios:
        raise SystemExit("--scenarios resolved to an empty list.")

    set_nature_rcparams()

    # ---- Pass 1: load every scenario's data up front, to derive shared axis
    # limits across the whole grid before any drawing happens ----
    per_scenario = []
    for scenario in scenarios:
        combined = load_regression_data(args.exp_name, scenario, args.exp_num,
                                        min_seq=args.min_seq, min_len=args.min_len)
        intra = load_intrahost_long_data(args.exp_name, scenario, args.exp_num,
                                         min_seq=args.min_seq, min_len=args.min_len)
        per_scenario.append((scenario, combined, intra))

    x_left_max = max(max(c["Sequencing_time"].max(), i["Sequencing_time"].max())
                     for _, c, i in per_scenario)
    x_right_max = max(i["Sequencing_time"].max() for _, _, i in per_scenario)

    # One shared y-limit across EVERY axis in the grid (left primary, left
    # secondary/twin, and right column all measure the same quantity --
    # substitutions/site -- so a single common scale makes the whole figure
    # directly comparable, not just within a column).
    y_max_std = max(
        combined[combined["individual_type"] == "standard"]["Distance_from_root"].max()
        if not combined[combined["individual_type"] == "standard"].empty else 0.0
        for _, combined, _ in per_scenario)
    y_max_intra = max(i["Distance_from_root"].max() for _, _, i in per_scenario)
    y_max = max(y_max_std, y_max_intra,
               args.target_osr_std * x_left_max, args.target_osr_long * x_left_max,
               args.target_osr_long * x_right_max)

    # ---- Pass 2: draw the grid, one row per scenario, shared limits per column ----
    n = len(per_scenario)
    fig, axes = plt.subplots(n, 2, figsize=(12, 4.2 * n), sharex='col')
    if n == 1:
        axes = axes.reshape(1, 2)  # keep 2D indexing uniform for a single row

    for row, (scenario, combined, intra) in enumerate(per_scenario):
        ax_global, ax_intra = axes[row, 0], axes[row, 1]

        # LEFT, primary axis: standard global clock + mixed (all individuals)
        draw_target(ax_global, args.target_osr_std, x_left_max, "#3b528b",
                    f"Target std = {args.target_osr_std}")
        std_sub = combined[combined["individual_type"] == "standard"]
        if std_sub.empty:
            print(f"[{scenario}][warn] global: no 'standard' rows; skipping.")
        else:
            slope, r2 = plot_cohort(ax_global, std_sub, "standard", x_left_max)
            print(f"[{scenario}][global/standard] OSR={slope:.6f} "
                 f"(R²={r2:.3f}, n={len(std_sub)})")

        slope_mix, r2_mix, _, _ = fit_cohort(combined)
        ax_global.plot([0, x_left_max], [0, slope_mix * x_left_max], color="black",
                       linewidth=1.5, zorder=3,
                       label=f"Mixed (all individuals): {slope_mix:.5f} s/s/y "
                             f"(R²={r2_mix:.2f})")
        print(f"[{scenario}][global/mixed] OSR={slope_mix:.6f} "
             f"(R²={r2_mix:.3f}, n={len(combined)})")

        ax_global.set_ylabel(f"{scenario}\nDistance (subs/site) -- global clock")
        ax_global.set_xlim(left=0, right=x_left_max * 1.02)
        ax_global.set_ylim(0, y_max * 1.05)
        ax_global.spines['top'].set_visible(False)

        # LEFT, secondary axis: long-shedder intra-host clock (line only)
        ax_global2 = ax_global.twinx()
        draw_target(ax_global2, args.target_osr_long, x_left_max, "#e41a1c",
                    f"Target long (intra-host) = {args.target_osr_long}")
        slope_long, r2_long = plot_cohort(ax_global2, intra, "long_shedder", x_left_max,
                                          scatter=False, label_suffix=" (intra-host)")
        print(f"[{scenario}][intra/long_shedder, left panel] OSR={slope_long:.6f} "
             f"(R²={r2_long:.3f}, n={len(intra)})")
        ax_global2.set_ylabel("Distance (subs/site) -- intra-host clock")
        ax_global2.set_ylim(0, y_max * 1.05)
        ax_global2.spines['top'].set_visible(False)

        h1, l1 = ax_global.get_legend_handles_labels()
        h2, l2 = ax_global2.get_legend_handles_labels()
        ax_global.legend(h1 + h2, l1 + l2, fontsize=6, loc="upper left", frameon=False)

        if row == 0:
            ax_global.set_title("Global clock (standard + mixed) & "
                                "long-shedder intra-host clock",
                                fontsize=7, fontweight="bold")
        if row == n - 1:
            ax_global.set_xlabel("Time (years)")

        # RIGHT: intra-host clock, long shedders only, shared column limits
        draw_target(ax_intra, args.target_osr_long, x_right_max, "#e41a1c",
                    f"Target long = {args.target_osr_long}")
        slope, r2 = plot_cohort(ax_intra, intra, "long_shedder", x_right_max)
        print(f"[{scenario}][intra/long_shedder] OSR={slope:.6f} "
             f"(R²={r2:.3f}, n={len(intra)})")
        ax_intra.set_xlim(left=0, right=x_right_max * 1.02)
        ax_intra.set_ylim(0, y_max * 1.05)
        ax_intra.spines['top'].set_visible(False)
        ax_intra.spines['right'].set_visible(False)
        ax_intra.legend(fontsize=6, loc="upper left", frameon=False)
        if row == 0:
            ax_intra.set_title("Intra-host clock (distance from infecting strain)",
                              fontsize=7, fontweight="bold")
        if row == n - 1:
            ax_intra.set_xlabel("Time since infection (years)")

    plt.tight_layout()

    if args.outfile:
        out_path = args.outfile
    else:
        experiment_string = f"{args.exp_name}_sanity_#{args.exp_num}"
        out_path = os.path.join(
            dm.get_experiment_plots_dir(experiment_string),
            f"{experiment_string}_all_scenarios_global_vs_intrahost.png")

    plt.savefig(out_path, format="png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n[saved] {out_path}")


if __name__ == "__main__":
    main()
