
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure 2: segmented IH lineage timelines (all vs only long shedders)
Now parameterized to pick ONE SSOD by:
  long_evo_rate_f, tau_3_long, R_long, long_shedders_ratio, exp_number, seed

Usage example:
  python figure_2_long_shedders.py \
    --evo 5 --tau3 72 --rlong 3 --ratio 0.2 --exp-number 1 --seed 1
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scripts.long_shedders_figures.preprocess_data as preprocess
import simplicity.plots_manager as pm
import simplicity.dir_manager as dm
import simplicity.intra_host_model as ih  # kept for consistency

# ---------- plotting ----------

def plot_segmented_infection_timeline_ax(ax, polished_data, colormap_df, t_final, only_long, label):
    """
    Plot IH lineage segments colored by lineage per individual.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    polished_data : pd.DataFrame (must have IH_lineages_trajectory, type, t_infection)
    colormap_df : pd.DataFrame from pm.make_lineages_colormap(ssod)
    t_final : float
    only_long : bool
    label : str
    """
    jitter_step = 0.15

    # Filter dataset
    data = polished_data[polished_data['type'] == 'long_shedder'] if only_long else polished_data
    data = data.sort_values(by='t_infection').reset_index(drop=True)

    for idx, row in data.iterrows():
        is_long = (row['type'] == 'long_shedder')
        traj = row['IH_lineages_trajectory']
        base_y = idx

        segments = [(lin, times['ih_birth'], times['ih_death']) for lin, times in traj.items()]
        segments.sort(key=lambda x: x[1])
        overlaps = any(segments[i][1] < segments[i - 1][2] for i in range(1, len(segments)))

        for j, (lineage, ih_start, ih_end) in enumerate(segments):
            y_j = base_y + (j * jitter_step - jitter_step / 2) if overlaps else base_y
            color = pm.get_lineage_color(lineage, colormap_df)
            alpha = 1.0 if is_long else 0.3
            zorder = 5 if is_long else 1

            ax.hlines(y_j, ih_start, ih_end, color=color, linewidth=2, alpha=alpha, zorder=zorder)
            if j == 0 and is_long:
                ax.plot(ih_start, y_j, marker='o', color=color, markersize=4, zorder=zorder + 1)

    # ax.axvline(t_final, color='gray', linestyle='--', linewidth=1.5)
    if only_long:
      ax.set_ylabel("Long-shedders inf. timeline")
    else:
      ax.set_ylabel("Individuals inf. timeline")
    ax.set_xlabel("Time")
    ax.set_xlim(0,t_final)
    ax.set_ylim(0)
    ax.grid(False)

    ax.text(-0.05, 1.05, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='left')


def plot_figure_2(ssod):
    """
    Main plot: A (all individuals) full-size.
    Inset: B (only long shedders) as an upper-left box on top of A.
    """
    polished_df, colormap_df, t_final = preprocess.prepare_figure_2_data(ssod)

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # Figure and main axis A
    fig, axA = plt.subplots(1, 1, figsize=(16, 8))

    # Draw A (all individuals)
    plot_segmented_infection_timeline_ax(
        axA, polished_df, colormap_df, t_final,
        only_long=False, label="A"
    )

    # Create inset axis B inside A, upper-left
    axB = axA.inset_axes([0.05, 0.5, 0.5, 0.5])  # [x0, y0, width, height] in Axes A coords
    
    axB.set_zorder(10)         # ensure B draws above A
    axB.set_facecolor("none")

    # Draw B (only long shedders) inside the inset
    plot_segmented_infection_timeline_ax(
        axB, polished_df, colormap_df, t_final,
        only_long=True, label=""  # suppress helper label
    )
    # Small "B" inside the inset
    axB.text(
        0.02, 0.98, "B",
        transform=axB.transAxes, ha="left", va="top",
        fontsize=12, fontweight="bold"
    )

    # Cosmetic tweaks
    for s in ("top", "right"):
        axA.spines[s].set_visible(False)
    axB.tick_params(axis="both", labelsize=8, pad=1)

    # Save
    experiment_name = dm.get_experiment_foldername_from_SSOD(ssod)
    experiment_plots_dir = dm.get_experiment_plots_dir(experiment_name)
    sim_out_name = dm.get_simulation_output_foldername_from_SSOD(ssod)
    seed = os.path.basename(ssod)
    outpath = os.path.join("Data", f"Figure_2_long_{sim_out_name}_{seed}.png")

    print(f"[SAVE] {outpath}")
    plt.savefig(outpath, format="png", dpi=300, bbox_inches="tight")
    plt.close()


# ---------- selection helpers ----------

def read_index_csv(path):
    """Read INDEX robustly, handling possible Windows encodings."""
    try:
        return pd.read_csv(path)
    except UnicodeDecodeError:
        return pd.read_csv(path, encoding="cp1252")

def find_experiment_row(index_df, evo, tau3, r_long, ratio):
    """
    Filter INDEX to the first row matching the parameters.
    Handles numeric vs string ambiguities by coercing to float.
    """
    df = index_df.copy()

    # Column names as used by your grid builder
    # fallbacks are handled just in case
    col_evo = "long_evo_rate_f" if "long_evo_rate_f" in df.columns else "evo_rate_f"
    col_tau = "tau_3_long"
    col_r   = "R_long"
    col_rat = "long_shedders_ratio"

    for c in [col_evo, col_tau, col_r, col_rat]:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in INDEX CSV.")

    df[col_evo] = pd.to_numeric(df[col_evo], errors="coerce")
    df[col_tau] = pd.to_numeric(df[col_tau], errors="coerce")
    df[col_r]   = pd.to_numeric(df[col_r],   errors="coerce")
    df[col_rat] = pd.to_numeric(df[col_rat], errors="coerce")

    sub = df[
        (np.isclose(df[col_evo], float(evo))) &
        (np.isclose(df[col_tau], float(tau3))) &
        (np.isclose(df[col_r],   float(r_long))) &
        (np.isclose(df[col_rat], float(ratio)))
    ]

    if sub.empty:
        raise ValueError(
            "No experiment matched the requested parameters:\n"
            f"  long_evo_rate_f={evo}, tau_3_long={tau3}, R_long={r_long}, long_shedders_ratio={ratio}"
        )
    return sub.iloc[0]

def get_ssod_for_params(index_csv, evo, tau3, r_long, ratio, exp_number, seed):
    """
    Resolve experiment name from INDEX and return the SSOD path for the given seed.
    """
    df = read_index_csv(index_csv)
    row = find_experiment_row(df, evo, tau3, r_long, ratio)

    gen_name = str(row["generated_experiment_name"])
    exp_name = f"{gen_name}_#{int(exp_number)}"

    print(f"[EXP] {exp_name} (from {index_csv})")

    sods = dm.get_simulation_output_dirs(exp_name)
    if not sods:
        raise RuntimeError(f"No Simulation Output Dirs found for experiment '{exp_name}'")

    sim_out_dir = sods[0]
    ssod = dm.get_ssod(sim_out_dir, int(seed))
    print(f"[SSOD] {ssod}")
    return ssod

# ---------- cli ----------

def main():
    ap = argparse.ArgumentParser(description="Plot Figure 2 for ONE SSOD selected by parameters.")
    ap.add_argument("--index-csv",
                    default="scripts/long_shedders_experiments/INDEX_exp_scripts.csv",
                    help="Path to INDEX_exp_scripts.csv")
    ap.add_argument("--evo", type=float, default=5.0, help="long_evo_rate_f value (default: 5)")
    ap.add_argument("--tau3", type=float, default=72.0, help="tau_3_long value (default: 72)")
    ap.add_argument("--rlong", type=float, default=3.0, help="R_long value (default: 3)")
    ap.add_argument("--ratio", type=float, default=0.01, help="long_shedders_ratio (default: 0.2)")
    ap.add_argument("--exp-number", type=int, default=1, help="Experiment number suffix '#N' (default: 1)")
    ap.add_argument("--seed", type=int, default=1, help="Seed number for SSOD (default: 1)")
    args = ap.parse_args()

    ssod = get_ssod_for_params(
        index_csv=args.index_csv,
        evo=args.evo,
        tau3=args.tau3,
        r_long=args.rlong,
        ratio=args.ratio,
        exp_number=args.exp_number,
        seed=args.seed
    )

    plot_figure_2(ssod)

if __name__ == "__main__":
    main()
