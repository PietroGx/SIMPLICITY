
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Figure 1: two subplots
A) Intra-host trajectories: baseline (default tau3) + one curve per unique tau_3_long from INDEX
B) Violin plots: baseline tau3 + one per unique tau_3_long (first experiment #1), same SSOD seed.
   Violin colors match the corresponding line colors in subplot A.
   Violin x-labels are only 'Xd' (e.g., '7.5d').

"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Project-local imports (expected to exist in your environment)
import scripts.long_shedders_figures.preprocess_data as preprocess
import simplicity.dir_manager as dm
import simplicity.intra_host_model as ih

# -----------------------------
# Config
# -----------------------------

DEFAULT_INDEX_CSV = "scripts/long_shedders_experiments/INDEX_exp_scripts.csv"
EXP_NUMBER_SUFFIX = "#1"           # always use first experiment folder
SEED_FOR_SSOD = 1                  # fixed seed for all selections
BASELINE_TAU3 = 7.5                # current default baseline tau3
TIME_GRID_DAYS = 300               # time horizon for IH plot

# -----------------------------
# Helpers
# -----------------------------

def read_index_csv(path):
    """Read CSV; try utf-8, then cp1252."""
    import pandas as pd
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="cp1252")


def get_unique_tau3_long(df):
    """Sorted unique tau_3_long values from INDEX."""
    if "tau_3_long" not in df.columns:
        raise ValueError("INDEX CSV missing column 'tau_3_long'.")
    return sorted(v for v in df["tau_3_long"].dropna().unique())


def first_row_for_tau(df, tau3_value):
    """First INDEX row matching a specific tau_3_long."""
    sub = df[df["tau_3_long"] == tau3_value]
    if sub.empty:
        return None
    return sub.iloc[0]


def resolve_experiment_name(row, exp_number_suffix=EXP_NUMBER_SUFFIX):
    """
    Build experiment folder name:
    prefer 'generated_experiment_name', else 'experiment_name'.
    Append '#1' (or configured suffix).
    """
    if row is None:
        return None
    if "generated_experiment_name" in row.index:
        base = str(row["generated_experiment_name"])
    elif "experiment_name" in row.index:
        base = str(row["experiment_name"])
    else:
        return None
    return f"{base}_{exp_number_suffix}"


def get_first_sod(exp_name):
    """First SOD for the experiment."""
    try:
        sods = dm.get_simulation_output_dirs(exp_name)
    except Exception as e:
        print(f"[WARN] Could not list SODs for experiment '{exp_name}': {e}")
        return None
    if not sods:
        print(f"[WARN] No SODs found for experiment '{exp_name}'.")
        return None
    return sods[0]


def get_ssod_for_seed(sod, seed_number):
    """SSOD for a given seed from a SOD."""
    try:
        return dm.get_ssod(sod, seed_number)
    except Exception as e:
        print(f"[WARN] Could not get SSOD for seed {seed_number} in SOD '{sod}': {e}")
        return None


def durations_for_group_from_ssod(ssod, which_type):
    """
    Use preprocess helper to get durations, then filter by 'type' column.
    Returns np.array of infection_duration.
    """
    try:
        df = preprocess.prepare_figure_1_data(ssod)
    except Exception as e:
        print(f"[WARN] prepare_figure_1_data failed for SSOD '{ssod}': {e}")
        return np.array([])
    if df is None or df.empty:
        return np.array([])
    if "infection_duration" not in df.columns:
        print("[WARN] Data missing 'infection_duration'.")
        return np.array([])
    if "type" not in df.columns:
        print("[WARN] No 'type' column; cannot isolate requested group.")
        return np.array([])

    sub = df[df["type"] == which_type]
    vals = sub["infection_duration"].dropna().to_numpy()
    return vals if vals.size else np.array([])

# -----------------------------
# Plotting
# -----------------------------

def subplot_A_intra_host(ax, tau3_long_values, baseline_tau3=BASELINE_TAU3):
    """
    Subplot A: plot intra-host probability curves using the working API:
    Host(...).data_plot_ih_solution(t0, t1, dt)[0]
    Returns:
        color_map: dict mapping float tau3 -> color (including baseline)
    """
    tvals = np.arange(0, float(TIME_GRID_DAYS), 1.0)
    color_map = {}

    # Baseline (normal) first
    try:
        host = ih.Host(tau_1=2.86, tau_2=3.91, tau_3=float(baseline_tau3), tau_4=8.0)
        y = host.data_plot_ih_solution(0, TIME_GRID_DAYS, 1)[0]
        line = ax.plot(tvals, y, lw=2.0, label=f"tau3={float(baseline_tau3):g}d")[0]
        color_map[float(baseline_tau3)] = line.get_color()
    except Exception as e:
        print(f"[WARN] Failed to plot baseline IH model: {e}")

    # Then one curve per tau_3_long
    for tau3 in tau3_long_values:
        try:
            host = ih.Host(tau_1=2.86, tau_2=3.91, tau_3=float(tau3), tau_4=8.0)
            y = host.data_plot_ih_solution(0, TIME_GRID_DAYS, 1)[0]
            line = ax.plot(tvals, y, lw=1.8, linestyle="--", label=f"tau3={float(tau3):g}d")[0]
            color_map[float(tau3)] = line.get_color()
        except Exception as e:
            print(f"[WARN] Failed plotting IH model for tau3={tau3}: {e}")
            continue

    ax.set_xlabel("time [days]")
    ax.set_ylabel("P(infectious after t)")
    ax.set_xlim(0, TIME_GRID_DAYS)
    ax.set_ylim(0.0, 1.0)
    # remove grid
    ax.grid(False)
    # Only show legend if something got plotted
    if len(ax.lines) > 0:
        ax.legend(title=None, fontsize=9)

    return color_map


def subplot_B_violins(ax, df_index, tau3_long_values, seed_number, baseline_tau3, color_map):
    """
    Subplot B: violins. Order and labels: baseline first, then ascending tau3_long.
    Violin colors match the corresponding line colors from subplot A.
    Labels are only 'Xd' strings.
    """
    labels = []
    data = []
    colors = []

    if df_index.empty:
        print("[WARN] INDEX is empty; cannot build violins.")
        return

    # -------- baseline violin (from the first INDEX row) --------
    first_row = df_index.iloc[0]
    exp_name = resolve_experiment_name(first_row, EXP_NUMBER_SUFFIX)
    if exp_name is None:
        print("[WARN] Could not resolve experiment name for baseline violin.")
    else:
        sod = get_first_sod(exp_name)
        if sod is not None:
            ssod = get_ssod_for_seed(sod, seed_number)
            if ssod is not None:
                dur = durations_for_group_from_ssod(ssod, which_type="normal")
                if dur.size > 0:
                    labels.append(f"{float(baseline_tau3):g}d")
                    data.append(dur)
                    colors.append(color_map.get(float(baseline_tau3), None))
                else:
                    print("[WARN] No durations found for baseline (normal) violin.")
        else:
            print(f"[WARN] Skipping baseline violin; no SOD for '{exp_name}'.")

    # -------- one violin per tau3_long (long_shedder group) --------
    for tau3 in tau3_long_values:
        row = first_row_for_tau(df_index, tau3)
        if row is None:
            print(f"[WARN] No INDEX row found for tau3={tau3}; skipping.")
            continue
        exp_name = resolve_experiment_name(row, EXP_NUMBER_SUFFIX)
        if exp_name is None:
            print(f"[WARN] Could not resolve experiment name for tau3={tau3}; skipping.")
            continue
        sod = get_first_sod(exp_name)
        if sod is None:
            print(f"[WARN] No SOD for experiment '{exp_name}'; skipping tau3={tau3}.")
            continue
        ssod = get_ssod_for_seed(sod, seed_number)
        if ssod is None:
            print(f"[WARN] No SSOD for seed {seed_number} in '{sod}'; skipping tau3={tau3}.")
            continue

        dur = durations_for_group_from_ssod(ssod, which_type="long_shedder")
        if dur.size == 0:
            print(f"[WARN] No durations for tau3={tau3} (type='long_shedder'); skipping.")
            continue

        labels.append(f"{float(tau3):g}d")
        data.append(dur)
        colors.append(color_map.get(float(tau3), None))

    if not data:
        print("[WARN] Nothing to plot in subplot B.")
        return

    # Draw violins
    parts = ax.violinplot(dataset=data, showmeans=False, showmedians=False, showextrema=False)

    # Color each violin to match subplot A
    for i, body in enumerate(parts["bodies"]):
        c = colors[i] if i < len(colors) else None
        if c is not None:
            body.set_facecolor(c)
            body.set_edgecolor(c)
        body.set_alpha(0.8)

    # Add inner quartiles
    for i, vals in enumerate(data, start=1):
        if vals.size < 2:
            continue
        q1, q2, q3 = np.percentile(vals, [25, 50, 75])
        ax.scatter([i], [q2], marker="o", s=12, zorder=3)
        ax.vlines([i], q1, q3, lw=2, zorder=3)

    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, rotation=0)
    ax.set_ylabel("infection duration [days]")
    # remove grid
    ax.grid(False)

# -----------------------------
# Main entry
# -----------------------------

def build_figure(index_csv=DEFAULT_INDEX_CSV, out_png=None):
    # Load INDEX
    df_index = read_index_csv(index_csv)

    # Unique tau_3_long values
    tau3_vals = get_unique_tau3_long(df_index)

    # Limit to at most 4 long tau values for total of up to 5 violins (including baseline)
    if len(tau3_vals) > 4:
        tau3_vals = tau3_vals[:4]

    # Layout
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axA, axB = axes

    # Subplot A (also returns color map)
    color_map = subplot_A_intra_host(axA, tau3_vals, baseline_tau3=BASELINE_TAU3)
    axA.set_title("A. intra-host trajectories", loc="left")

    # Subplot B, colored to match subplot A, labels are only 'Xd'
    subplot_B_violins(
        axB,
        df_index=df_index,
        tau3_long_values=tau3_vals,
        seed_number=SEED_FOR_SSOD,
        baseline_tau3=BASELINE_TAU3,
        color_map=color_map,
    )
    axB.set_title("B. durations by group", loc="left")

    plt.tight_layout()
    if out_png:
        os.makedirs(os.path.dirname(out_png), exist_ok=True)
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
    return fig


def main():
    ap = argparse.ArgumentParser(description="Build Figure 1 (long shedders).")
    ap.add_argument("--index_csv", type=str, default=DEFAULT_INDEX_CSV,
                    help="Path to INDEX_exp_scripts.csv")
    ap.add_argument("--out_png", type=str, default="Data/Figure_1_long.png",
                    help="Optional output PNG path")
    args = ap.parse_args()

    build_figure(index_csv=args.index_csv, out_png=args.out_png)


if __name__ == "__main__":
    main()