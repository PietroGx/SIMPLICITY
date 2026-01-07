# This file is part of SIMPLICITY
# Copyright (C) …
#
# Clade-level metrics across seeds:
#   - Peak frequency winner
#   - Total infections (burden) winner
#   - Longest survival winner
#   - Fastest growth winner (time 1% -> 50%)
#
# Produces 4 pies per SOD summarizing winner labels across seeds.

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import simplicity.output_manager as om
import simplicity.dir_manager as dm
import simplicity.plots_manager as pm
import simplicity.clustering as cl


# ----------------------------- Helpers ---------------------------------

def _pivot_lineage_freqs(lineage_frequency_df: pd.DataFrame):
    f = lineage_frequency_df.pivot(
        index="Time_sampling", columns="Lineage_name", values="Frequency_at_t"
    )
    return f.sort_index()


def _build_clade_freq_df(clade_to_lineages: dict, F_lineage: pd.DataFrame) -> pd.DataFrame:
    series = []
    for clade, members in clade_to_lineages.items():
        cols = [c for c in members if c in F_lineage.columns]
        if not cols:
            continue
        s = F_lineage[cols].sum(axis=1)
        s.name = clade
        series.append(s)
    if not series:
        return pd.DataFrame(index=F_lineage.index)
    return pd.concat(series, axis=1)


def _peak_winner(F_clade: pd.DataFrame):
    if F_clade.empty or F_clade.shape[1] == 0:
        return None, np.nan
    peaks = F_clade.max(axis=0)
    winner = peaks.idxmax()
    return winner, float(peaks.loc[winner])


def _burden_winner_from_phylo(phylo_df: pd.DataFrame, clade_to_lineages: dict):
    if "Total_infections" not in phylo_df.columns:
        raise ValueError("phylogenetic_data is missing 'Total_infections'.")
    lin_to_total = dict(zip(phylo_df["Lineage_name"], phylo_df["Total_infections"]))
    totals = {clade: sum(int(lin_to_total.get(l, 0) or 0) for l in members)
              for clade, members in clade_to_lineages.items()}
    if not totals:
        return None, np.nan
    winner = max(totals.items(), key=lambda kv: kv[1])[0]
    return winner, float(totals[winner])


def _survival_winner(F_clade: pd.DataFrame, eps: float = 1e-4):
    if F_clade.empty or F_clade.shape[1] == 0:
        return None, np.nan
    times = F_clade.index.to_numpy(dtype=float)
    durations = {}
    for clade in F_clade.columns:
        y = F_clade[clade].to_numpy(dtype=float)
        mask = y > eps
        if not mask.any():
            durations[clade] = 0.0
            continue
        t_first = times[mask.argmax()]
        t_last = times[len(mask) - 1 - mask[::-1].argmax()]
        durations[clade] = max(0.0, float(t_last - t_first))
    winner = max(durations.items(), key=lambda kv: kv[1])[0]
    return winner, float(durations[winner])


def _first_sustained_crossing_time(times, values, threshold, sustain_k=1):
    T = np.asarray(times, dtype=float)
    F = np.asarray(values, dtype=float)
    n = len(T)
    if n < 2:
        return None
    if F[0] >= threshold:
        future = F[1:1+sustain_k]
        if len(future) < sustain_k or np.all(future >= threshold):
            return float(T[0])
    for i in range(n - 1):
        f0, f1 = F[i], F[i+1]
        if f0 < threshold <= f1 and (f1 - f0) > 0:
            frac = (threshold - f0) / (f1 - f0)
            t_cross = T[i] + frac * (T[i+1] - T[i])
            future = F[i+1:i+1+sustain_k]
            if len(future) < sustain_k or np.all(future >= threshold):
                return float(t_cross)
    return None


def _fastest_growth_winner(F_clade: pd.DataFrame, lower=0.01, upper=0.50, sustain_k=1):
    if F_clade.empty or F_clade.shape[1] == 0:
        return None, np.nan
    times = F_clade.index.to_numpy(dtype=float)
    candidates = {}
    for clade in F_clade.columns:
        y = F_clade[clade].to_numpy(dtype=float)
        tL = _first_sustained_crossing_time(times, y, lower, sustain_k=sustain_k)
        if tL is None:
            continue
        tail_mask = times > tL
        if not tail_mask.any():
            continue
        tU = _first_sustained_crossing_time(times[tail_mask], y[tail_mask], upper, sustain_k=sustain_k)
        if tU is None:
            continue
        delta = float(tU - tL)
        if delta >= 0:
            candidates[clade] = delta
    if not candidates:
        return None, np.nan
    winner = min(candidates.items(), key=lambda kv: kv[1])[0]
    return winner, float(candidates[winner])


def _label_for_pie(raw_label: str):
    """
    Normalize labels to {'normal','long','mixed','founder'}.
    Safe against None/NaN.
    """
    if raw_label is None:
        return np.nan
    if isinstance(raw_label, float):
        try:
            if np.isnan(raw_label):
                return np.nan
        except Exception:
            pass
    v = str(raw_label).strip().lower()
    if v.startswith("long"):
        return "long"
    if v.startswith("norm"):
        return "normal"
    if v == "mixed":
        return "mixed"
    if v == "founder":
        return "founder"
    return np.nan


# ----------------------- Per-SSOD computation --------------------------

def compute_clade_metric_winners_for_ssod(exp_name: str, sod: str, ssod: str, cluster_threshold: int):
    phylo_df = om.read_phylogenetic_data(ssod)
    lf       = om.read_lineage_frequency(ssod)

    clade_to_lineages, lineage_to_clade, per_clade_mut_df, clade_meta_df = cl.cluster_lin_into_clades_with_meta(
        phylo_df, shared_mut_threshold=cluster_threshold
    )

    # NEW: pass clade_meta_df here
    clade_labels_series, _summary = cl.label_clades_from_definers(per_clade_mut_df, clade_meta_df)
    clade_labels = clade_labels_series.to_dict() if clade_labels_series is not None else {}

    # non-root clades
    non_root = clade_meta_df.loc[
        (clade_meta_df["parent_clade"].notna()) & (clade_meta_df["n_defining"] > 0),
        "clade"
    ].tolist()

    F_lineage = _pivot_lineage_freqs(lf)
    F_clade   = _build_clade_freq_df(clade_to_lineages, F_lineage)

    def label_of(clade):
        if not clade:
            return None
        raw = clade_labels.get(clade)
        if raw in ("normal", "long_shedder", "founder"):
            return _label_for_pie(raw)
        row = clade_meta_df.loc[clade_meta_df["clade"] == clade]
        if not row.empty and (row["parent_clade"].iloc[0] is None or row["n_defining"].iloc[0] == 0):
            return "founder"
        return _label_for_pie(raw if raw is not None else "normal")

    F_nr = F_clade.reindex(columns=non_root).dropna(axis=1, how="all") if non_root else F_clade.iloc[:, 0:0]

    peak_winner, peak_val = _peak_winner(F_nr) if not F_nr.empty else (None, None)
    survival_winner, survival_val = _survival_winner(F_nr) if not F_nr.empty else (None, None)
    growth_winner, growth_val = _fastest_growth_winner(F_nr) if not F_nr.empty else (None, None)
    burden_winner, burden_val = _burden_winner_from_phylo(phylo_df, clade_to_lineages)

    return {
        "ssod": ssod,
        "peak_clade": peak_winner,
        "peak_label": label_of(peak_winner) if peak_winner else None,
        "burden_clade": burden_winner,
        "burden_label": label_of(burden_winner) if burden_winner else None,
        "survival_clade": survival_winner,
        "survival_label": label_of(survival_winner) if survival_winner else None,
        "growth_clade": growth_winner,
        "growth_label": label_of(growth_winner) if growth_winner else None,
        "peak_value": peak_val,
        "burden_value": burden_val,
        "survival_value": survival_val,
        "growth_value": growth_val,
    }

def add_clade_and_panel_legend(fig, fontsize=10, pad=0.012):
    """
    Place a unified legend just OUTSIDE the right edge of the subplot grid,
    top-aligned. Works regardless of how many subplots there are.
    """
    from matplotlib.patches import Patch
    from matplotlib.transforms import Bbox

    colors = {"normal":"#4daf4a","long":"#e41a1c","mixed":"#377eb8","founder":"#8172B2"}
    label_order = ["normal","long","mixed","founder"]
    clade_handles = [Patch(facecolor=colors[k], edgecolor="white", label=k) for k in label_order]

    panel_desc = ["i. Peak frequency",
                  "ii. Total infections (burden)",
                  "iii. Longest survival",
                  "iv. Fastest growth (1%→50%)"]
    panel_handles = [Patch(facecolor="none", edgecolor="none", label=lab) for lab in panel_desc]

    handles = clade_handles + panel_handles
    labels  = [h.get_label() for h in handles]

    # Union of all visible axes in FIGURE coords
    boxes = [ax.get_position() for ax in fig.axes if ax.get_visible()]
    union = Bbox.from_extents(min(b.x0 for b in boxes), min(b.y0 for b in boxes),
                              max(b.x1 for b in boxes), max(b.y1 for b in boxes))

    # Anchor just outside the grid's right edge, top-aligned
    anchor = (union.x1 + pad, union.y1)

    leg = fig.legend(handles=handles, labels=labels,
                     loc="upper left",                # left edge of legend sits on anchor
                     bbox_to_anchor=anchor,
                     bbox_transform=fig.transFigure,  # interpret anchor in figure coords
                     frameon=False, fontsize=fontsize,
                     handlelength=1.2, handletextpad=0.6, labelspacing=0.5)
    return leg



def summarize_sod_to_pies(experiment_name: str,
                          sod_path: str,
                          cluster_threshold: int,
                          min_days: int,
                          axes=None,
                          savefig: bool = True,
                          showlegend: bool = False):
    """
    Draw 4 pie plots for one SOD.

    Behavior:
      - No subplot titles (panel letters i.–iv. only).
      - Optional external legend (set showlegend=True). The legend should be
        provided by a helper function `add_clade_and_panel_legend(...)` if present.
      - Can either save a standalone figure (default) or draw on provided axes and return.

    Returns:
      - If savefig=False: (fig, axes) for further composition.
      - If savefig=True and axes is None: saves PNG and returns None.
    """
    pm.apply_plos_rcparams()

    # Collect SSODs
    ssods = sorted(dm.get_seeded_simulation_output_dirs(sod_path))

    # Filter by minimum duration
    kept_ssods = []
    for ssod in ssods:
        try:
            ft = om.read_final_time(ssod)  # should return a number of days or None
            ft_val = float(ft) if ft is not None else 0.0
        except Exception as e:
            print(f"[SKIP] {ssod} final_time unavailable ({e}); excluded from pies.")
            continue
        if ft_val >= float(min_days):
            kept_ssods.append(ssod)
        else:
            print(f"[SKIP] {ssod} final_time={ft_val} < min_days={min_days}")

    print(f"[INFO] SOD filter by min_days={min_days}: kept {len(kept_ssods)} / {len(ssods)} SSOD(s).")
    if not kept_ssods:
        raise RuntimeError(f"No SSODs meet min_days={min_days} in SOD: {sod_path}")

    # Compute winners per kept seed
    results = []
    for ssod in kept_ssods:
        res = compute_clade_metric_winners_for_ssod(experiment_name, sod_path, ssod, cluster_threshold)
        results.append(res)

    if not results:
        raise RuntimeError("No results produced for any seed.")

    df = pd.DataFrame(results)

    # Normalize labels to {'normal','long','mixed','founder'}; leave NaN as-is
    for col in ["peak_label", "burden_label", "survival_label", "growth_label"]:
        df[col] = df[col].map(_label_for_pie)

    # Counting helper
    def _counts(col_name: str):
        s = df[col_name].dropna()
        order = ["normal", "long", "mixed", "founder"]
        c = s.value_counts().reindex(order, fill_value=0)
        n = int(c.sum())
        return c, n

    # Prepare figure/axes
    created_fig = False
    if axes is None:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        axes = axes.ravel()
        created_fig = True
    else:
        axes = np.ravel(axes)
        fig = axes[0].figure

    # Visual config
    colors = {"normal": "#4daf4a", "long": "#e41a1c", "mixed": "#377eb8", "founder": "#8172B2"}
    label_order = ["normal", "long", "mixed", "founder"]
    metrics = ["peak_label", "burden_label", "survival_label", "growth_label"]
    panel_labels = ["i.", "ii.", "iii.", "iv."]
    panel_titles = {
        "i.": "Peak frequency",
        "ii.": "Total infections (burden)",
        "iii.": "Longest survival",
        "iv.": "Fastest growth (1%->50%)",
    }

    # Draw pies
    for idx, (ax, metric) in enumerate(zip(axes, metrics)):
        c, n = _counts(metric)
        if n == 0:
            ax.text(0.5, 0.5, "No seeds", ha="center", va="center")
            ax.axis("off")
            continue

        ax.pie(
            c.values,
            labels=None,
            autopct=None,
            startangle=90,
            colors=[colors[k] for k in label_order],
            wedgeprops={"linewidth": 0.5, "edgecolor": "white"},
        )
        # Panel letter only (no title)
        ax.text(-0.15, 1.05, panel_labels[idx],
                transform=ax.transAxes, fontsize=12, fontweight="bold",
                va="top", ha="right")

    # Optional external legend (upper right), only if requested.
    if showlegend:
        add_clade_and_panel_legend(
            fig=fig)
        # leave layout room for the legend
        plt.tight_layout()
    else:
        plt.tight_layout()

    # Save or return
    if created_fig and savefig:
        exp_plot_dir = dm.get_experiment_plots_dir(experiment_name)
        os.makedirs(exp_plot_dir, exist_ok=True)
        sod_name = os.path.basename(sod_path.rstrip("/"))
        figfile = os.path.join(exp_plot_dir, f"{sod_name}_clade_winner_pies.png")
        plt.savefig(figfile, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[INFO] Saved pies: {figfile}")
        return None
    else:
        return fig, axes


# ----------------------------- CLI -------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compute clade-level winner stats and plot pies per SOD."
    )
    parser.add_argument("experiment_name", type=str, help="Experiment name containing SODs.")
    parser.add_argument(
        "cluster_threshold",
        type=int,
        help="Shared mutation threshold for clade creation (e.g., 5).",
    )
    parser.add_argument("--min-days", type=int, default=0,
                        help="Only include SSODs with final_time >= this many days (default: 0).")
    args = parser.parse_args()

    experiment_name = args.experiment_name
    cluster_threshold = args.cluster_threshold
    min_days = args.min_days

    sods = dm.get_simulation_output_dirs(experiment_name)
    if not sods:
        raise ValueError(f"No SODs found for experiment '{experiment_name}'.")

    print(f"[INFO] Found {len(sods)} SOD(s) for experiment '{experiment_name}'.")
    for sod in sods:
        print(f"[INFO] Processing SOD: {sod}")
        try:
            summarize_sod_to_pies(experiment_name, sod, cluster_threshold, min_days,
                                        savefig=True, showlegend=True)
        except Exception as e:
            print(f"[WARN] Skipped SOD due to error: {sod}\n  -> {e}")


if __name__ == "__main__":
    main()
