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
import math
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
    """Return time-sorted pivots:
       F_lineage: index=Time_sampling, columns=Lineage_name, values=Frequency_at_t
    """
    f = lineage_frequency_df.pivot(
        index="Time_sampling", columns="Lineage_name", values="Frequency_at_t"
    )
    f = f.sort_index()
    return f


def _build_clade_freq_df(clade_to_lineages: dict, F_lineage: pd.DataFrame) -> pd.DataFrame:
    """Sum lineage frequency columns into clade columns over the same time index."""
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
    """Return (clade, peak_value). None if no columns."""
    if F_clade.empty or F_clade.shape[1] == 0:
        return None, np.nan
    peaks = F_clade.max(axis=0)
    winner = peaks.idxmax()
    return winner, float(peaks.loc[winner])


def _burden_winner_from_phylo(phylo_df: pd.DataFrame, clade_to_lineages: dict):
    """Use phylogenetic_data['Total_infections'] to compute burden per clade.
       Returns (clade, total_infections). Raises if column missing.
    """
    if "Total_infections" not in phylo_df.columns:
        raise ValueError("phylogenetic_data is missing 'Total_infections'.")
    lin_to_total = dict(zip(phylo_df["Lineage_name"], phylo_df["Total_infections"]))
    totals = {}
    for clade, members in clade_to_lineages.items():
        tot = 0
        for lin in members:
            tot += int(lin_to_total.get(lin, 0) or 0)
        totals[clade] = tot
    if not totals:
        return None, np.nan
    winner = max(totals.items(), key=lambda kv: kv[1])[0]
    return winner, float(totals[winner])


def _survival_winner(F_clade: pd.DataFrame, eps: float = 1e-4):
    """Longest time with F_clade(t) > eps. Returns (clade, duration_days)."""
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
        # first & last sample > eps; use sample times (irregular ok)
        t_first = times[mask.argmax()]
        t_last = times[len(mask) - 1 - mask[::-1].argmax()]
        durations[clade] = max(0.0, float(t_last - t_first))
    winner = max(durations.items(), key=lambda kv: kv[1])[0]
    return winner, float(durations[winner])


def _first_sustained_crossing_time(times, values, threshold, sustain_k=1):
    """Linear-interpolated first crossing from below of 'threshold' with a sustain check:
       After the crossing segment (i,i+1), require the next 'sustain_k' points (if exist)
       to remain >= threshold. Returns crossing time or None.
    """
    T = np.asarray(times, dtype=float)
    F = np.asarray(values, dtype=float)

    n = len(T)
    if n < 2:
        return None

    # Handle starting above threshold: treat T0 as crossing if sustained
    if F[0] >= threshold:
        # must remain >= threshold for sustain_k consecutive future samples if they exist
        future = F[1:1+sustain_k]
        if len(future) < sustain_k or np.all(future >= threshold):
            return float(T[0])

    for i in range(n - 1):
        f0, f1 = F[i], F[i+1]
        if f0 < threshold <= f1 and (f1 - f0) > 0:
            # linear interpolation within [Ti, Ti+1]
            frac = (threshold - f0) / (f1 - f0)
            t_cross = T[i] + frac * (T[i+1] - T[i])
            # sustain check on next sustain_k samples
            future = F[i+1:i+1+sustain_k]
            if len(future) < sustain_k or np.all(future >= threshold):
                return float(t_cross)
    return None


def _fastest_growth_winner(F_clade: pd.DataFrame, lower=0.01, upper=0.50, sustain_k=1):
    """Winner is clade with smallest time to go from lower→upper (1%→50%).
       Returns (clade, delta_time). Clades that never reach either threshold are ignored.
    """
    if F_clade.empty or F_clade.shape[1] == 0:
        return None, np.nan

    times = F_clade.index.to_numpy(dtype=float)
    candidates = {}
    for clade in F_clade.columns:
        y = F_clade[clade].to_numpy(dtype=float)
        tL = _first_sustained_crossing_time(times, y, lower, sustain_k=sustain_k)
        if tL is None:
            continue
        # require upper crossing to be after lower crossing: search on the tail
        # Restrict to samples strictly after tL for robustness
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

    # smallest delta wins
    winner = min(candidates.items(), key=lambda kv: kv[1])[0]
    return winner, float(candidates[winner])


def _label_for_pie(raw_label: str) -> str:
    """Normalize labels to exactly {'normal','long','mixed'}."""
    if not isinstance(raw_label, str):
        return "mixed"
    v = raw_label.strip().lower()
    if v.startswith("long"):
        return "long"
    if v.startswith("norm"):
        return "normal"
    if v == "mixed":
        return "mixed"
    return "mixed"


def _tiebreak_by_start_time(candidates, clade_meta_df: pd.DataFrame):
    """Given a list of clade names with equal scores, pick the one with earliest start_time,
       then lexicographic."""
    if not candidates:
        return None
    if clade_meta_df is None or clade_meta_df.empty:
        return sorted(candidates)[0]
    start = clade_meta_df.set_index("clade")["start_time"].to_dict()
    return min(candidates, key=lambda c: (start.get(c, float("inf")), str(c)))


# ----------------------- Per-SSOD computation --------------------------

def compute_clade_metric_winners_for_ssod(ssod: str, cluster_threshold: int):
    """
    For a single seeded run (SSOD):
      - cluster into clades
      - build clade frequency matrix
      - compute winners for: peak, burden, survival, growth
      - label each winner; for burden the root can win and is labeled from its Host_type
    Returns a dict with winners, labels, and values.
    """
    # Read data
    phylo_df = om.read_phylogenetic_data(ssod)
    lf       = om.read_lineage_frequency(ssod)

    # Build clades (tree-based) + per-clade defining-mutation tables + metadata
    clade_to_lineages, lineage_to_clade, per_clade_mut_df, clade_meta_df = cl.cluster_lin_into_clades_with_meta(
        phylo_df, shared_mut_threshold=cluster_threshold
    )

    # Build labels from definers (non-root clades will generally have definers)
    clade_labels_series, _summary = cl.label_clades_from_definers(per_clade_mut_df)
    clade_labels = clade_labels_series.to_dict() if clade_labels_series is not None else {}

    # Identify non-root clades for metrics that exclude the root
    non_root = clade_meta_df.loc[
        (clade_meta_df["parent_clade"].notna()) & (clade_meta_df["n_defining"] > 0),
        "clade"
    ].tolist()
    print(f"[FILTER] non-root clades kept: {len(non_root)} / {clade_meta_df.shape[0]}")

    # Build clade frequency time series
    F_lineage = _pivot_lineage_freqs(lf)                 # index = time, cols = lineage
    F_clade   = _build_clade_freq_df(clade_to_lineages, F_lineage)  # sum lineages → clade

    # Helper: normalize host labels
    def _norm_host(v):
        if pd.isna(v): return None
        s = str(v).strip().lower()
        if s.startswith("long"):   return "long_shedder"
        if s.startswith("norm"):   return "normal"
        if s in ("", "none", "nan", "unknown"): return None
        return s

    # Helper: label for pies with root fallback when definers are empty/unknown
    # - If clade has a clear label from definers, use it (and pass through _label_for_pie).
    # - If it's the ROOT (no definers), label from the root lineage Host_type.
    def label_of(clade):
        if not clade:
            return None
        raw = clade_labels.get(clade)  # may be 'normal' / 'long_shedder' / 'mixed' / 'unknown' / None
        # If we have a clear label already, return a pie-safe label
        if raw in ("normal", "long_shedder"):
            return _label_for_pie(raw)
        # Root fallback: no definers → find root lineage host
        row = clade_meta_df.loc[clade_meta_df["clade"] == clade]
        if not row.empty and (row["parent_clade"].iloc[0] is None or row["n_defining"].iloc[0] == 0):
            root_lin = row["root_lineage"].iloc[0]
            host_raw = phylo_df.loc[phylo_df["Lineage_name"] == root_lin, "Host_type"]
            host = _norm_host(host_raw.iloc[0]) if not host_raw.empty else None
            if host in ("normal", "long_shedder"):
                return _label_for_pie(host)
        # Otherwise fall back to pie-safe mapping of whatever we have
        return _label_for_pie(raw if raw is not None else "unknown")

    # ---------- Metrics ----------
    # Exclude root for PEAK / SURVIVAL / GROWTH
    if len(non_root) == 0:
        F_nr = F_clade.iloc[:, 0:0]  # empty
    else:
        F_nr = F_clade.reindex(columns=non_root).dropna(axis=1, how="all")

    # Peak
    if F_nr.shape[1] == 0:
        peak_winner, peak_val = None, None
    else:
        peak_winner, peak_val = _peak_winner(F_nr)

    # Survival
    if F_nr.shape[1] == 0:
        survival_winner, survival_val = None, None
    else:
        survival_winner, survival_val = _survival_winner(F_nr, eps=1e-4)

    # Growth (1% → 50%)
    if F_nr.shape[1] == 0:
        growth_winner, growth_val = None, None
    else:
        growth_winner, growth_val = _fastest_growth_winner(F_nr, lower=0.01, upper=0.50, sustain_k=1)

    # Burden: INCLUDE ROOT
    # (sum Total_infections across member lineages for all clades, including root)
    burden_winner, burden_val = _burden_winner_from_phylo(phylo_df, clade_to_lineages)

    # ---------- Labels (None if no winner) ----------
    out = {
        "ssod": ssod,
        "peak_clade": peak_winner,
        "peak_label": label_of(peak_winner) if peak_winner else None,
        "burden_clade": burden_winner,
        "burden_label": label_of(burden_winner) if burden_winner else None,  # root gets host-type fallback
        "survival_clade": survival_winner,
        "survival_label": label_of(survival_winner) if survival_winner else None,
        "growth_clade": growth_winner,
        "growth_label": label_of(growth_winner) if growth_winner else None,
        "peak_value": peak_val,
        "burden_value": burden_val,
        "survival_value": survival_val,
        "growth_value": growth_val,
    }

    # Debug prints
    print("[CHK] F_clade shape:", F_clade.shape, "cols:", list(F_clade.columns)[:5], "…")
    print("[CHK] (non-root) F_nr shape:", F_nr.shape, "cols:", list(F_nr.columns)[:5], "…")
    print("[CHK] peak:", peak_winner, peak_val,
          "| burden:", burden_winner, burden_val,
          "| survival:", survival_winner, survival_val,
          "| growth:", growth_winner, growth_val)
    if F_nr.empty:
        print("[WHY] Non-root F_clade is empty — either no non-root clades or all-NaN after filtering.")
    if burden_winner is None:
        missing = "Total_infections" not in phylo_df.columns
        print(f"[WHY] burden winner None — {'missing Total_infections' if missing else 'no clades?'}")

    return out



# ----------------------- SOD-level aggregation -------------------------

def summarize_sod_to_pies(experiment_name: str, sod_path: str, cluster_threshold: int):
    """
    Run metric winners across all seeds (SSODs) in a single SOD and plot four pies
    for the label distribution of winners (normal/long/mixed) for:
      - Peak frequency
      - Burden (total infections)
      - Survival
      - Growth (1% -> 50% time)
    Saves CSVs and a PNG/TIFF figure in out_dir (defaults to SOD's plots dir).
    """
    pm.apply_plos_rcparams()

    ssods = dm.get_seeded_simulation_output_dirs(sod_path)
    ssods = sorted(ssods)

    if not ssods:
        raise ValueError(f"No seeded runs found in SOD: {sod_path}")

    results = []
    for ssod in ssods:
        try:
            res = compute_clade_metric_winners_for_ssod(ssod, cluster_threshold)
            results.append(res)
        except Exception as e:
            print(f"[WARN] Skipping SSOD due to error: {ssod}\n  -> {e}")

    if not results:
        raise RuntimeError("No results produced for any seed; aborting.")

    df = pd.DataFrame(results)

    # Normalize labels (already normalized in compute function, but be safe)
    for col in ["peak_label", "burden_label", "survival_label", "growth_label"]:
        df[col] = df[col].map(_label_for_pie)

    # Aggregates
    def _counts(col):
        c = df[col].value_counts().reindex(["normal", "long", "mixed"], fill_value=0)
        n = int(c.sum())
        pct = (c / max(1, n) * 100.0).round(1)
        return c, pct, n

    agg = {}
    for metric, col in [("Peak", "peak_label"),
                        ("Burden", "burden_label"),
                        ("Survival", "survival_label"),
                        ("Growth", "growth_label")]:
        counts, pct, n = _counts(col)
        agg[metric] = {"counts": counts, "pct": pct, "n": n}
        
    out_dir = dm.get_experiment_output_dir(experiment_name)
    # Save CSVs
    sod_name = os.path.basename(sod_path.rstrip("/"))
    winners_csv = os.path.join(out_dir, f"{sod_name}_clade_winners_by_seed.csv")
    df.to_csv(winners_csv, index=False)
    print(f"[INFO] Saved per-seed winners: {winners_csv}")

    summary_rows = []
    for metric, d in agg.items():
        row = {"SOD": sod_name, "metric": metric, "n_seeds": d["n"]}
        row.update({f"{lab}_count": int(d["counts"].get(lab, 0)) for lab in ["normal", "long", "mixed"]})
        row.update({f"{lab}_pct": float(d["pct"].get(lab, 0.0)) for lab in ["normal", "long", "mixed"]})
        summary_rows.append(row)
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(out_dir, f"{sod_name}_clade_winner_label_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"[INFO] Saved summary: {summary_csv}")

    # Plot pies
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    axes = axes.ravel()
    titles = ["Peak frequency", "Total infections (burden)", "Longest survival", "Fastest growth (1%→50%)"]
    metrics = ["Peak", "Burden", "Survival", "Growth"]
    label_order = ["normal", "long", "mixed"]
    colors = {"normal": "#4daf4a", "long": "#e41a1c", "mixed": "#377eb8"}

    for ax, metric, title in zip(axes, metrics, titles):
        d = agg[metric]
        c = d["counts"].reindex(label_order, fill_value=0)
        if d["n"] == 0:
            ax.text(0.5, 0.5, "No seeds", ha="center", va="center")
            ax.axis("off")
            continue
        labels = [f"{lab} ({int(c[lab])})" if c[lab] > 0 else "" for lab in label_order]
        ax.pie(c.values, labels=labels,
               autopct=lambda p: f"{p:.0f}%" if p > 0.01 else "",
               startangle=90,
               colors=[colors[k] for k in label_order],
               wedgeprops={"linewidth": 0.5, "edgecolor": "white"})
        ax.set_title(f"{title}\n(n={d['n']})")

    plt.tight_layout()
    
    exp_plot_dir = dm.get_experiment_plots_dir(experiment_name)
    figfile = os.path.join(exp_plot_dir, f"{sod_name}_clade_winner_pies.png")
    plt.savefig(figfile, dpi=300, bbox_inches="tight")
    print(f"[INFO] Saved pies: {figfile}")
    plt.close(fig)


# ----------------------------- CLI -------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Compute clade-level winner stats and plot pies per SOD.")
    
    # parser.add_argument("experiment_name", type=str, help="Experiment name containing SODs.")
    # parser.add_argument("cluster_threshold", type=int, help="Shared mutation threshold for clade creation (e.g., 5).")
    
    # args = parser.parse_args()


    experiment_name = 'test_local_experiment_serial_#4'
    cluster_threshold = 5
    
    # Style
    pm.apply_plos_rcparams()

    sods = dm.get_simulation_output_dirs(experiment_name)
    if not sods:
        raise ValueError(f"No SODs found for experiment '{experiment_name}'.")

    print(f"[INFO] Found {len(sods)} SOD(s) for experiment '{experiment_name}'.")
    for sod in sods:
        print(f"[INFO] Processing SOD: {sod}")
        try:
            summarize_sod_to_pies(experiment_name, sod, cluster_threshold)
        except Exception as e:
            print(f"[WARN] Skipped SOD due to error: {sod}\n  -> {e}")


if __name__ == "__main__":
    main()
