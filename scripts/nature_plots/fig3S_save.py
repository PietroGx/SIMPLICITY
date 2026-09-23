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
# Figure 3S -- clade-origin pies for EVERY scenario.
#
# Figure 3 panel A shows the 2x2 pies for one representative group. This is
# the same four metrics for all scenarios side by side, so the trend across
# long-shedder burden is visible rather than asserted.
#
# Built on scripts/long_paper_figures/plot_figures.py::build_figure_4, which
# already gridded these pies -- same summarize_sod_pies / plot_fig4_pies /
# add_global_pie_legend. Only the data access differs: that version indexed an
# old parameter sweep (M, ratio, tau, R) via a master log, this one walks the
# named scenarios of the current pipelines.
# ============================================================================

import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import fig3_preprocess_data as preproc
from _scenarios import BOUND_EXP_NAME, figure_path, resolve_scenarios

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'long_paper_figures'))
from long_shedders_plots import plot_fig4_pies, add_global_pie_legend


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Figure 3S: clade-origin pies for every scenario.")
    parser.add_argument('--exp-num', type=int, required=True)
    parser.add_argument('--exp-name', type=str, default=BOUND_EXP_NAME)
    parser.add_argument('--cluster-threshold', type=int, default=5)
    parser.add_argument('--min-days', type=int, default=100)
    parser.add_argument('--format', type=str, choices=['pdf', 'png'], default='png')
    return parser.parse_args()


def set_nature_rcparams():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 7,
        'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })


def build_figure_3S(exp_num, exp_name, cluster_threshold, min_days, fmt):
    set_nature_rcparams()

    scenarios, _ = resolve_scenarios(exp_name, exp_num)
    if not scenarios:
        raise SystemExit(f"No scenario output found for {exp_name} #{exp_num}.")

    fig = plt.figure(figsize=(180 / 25.4, 62 / 25.4))
    outer = gridspec.GridSpec(1, len(scenarios), figure=fig, wspace=0.28,
                              top=0.70, bottom=0.04, left=0.03, right=0.98)

    for col, scenario in enumerate(scenarios):
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 2, subplot_spec=outer[0, col], wspace=0.08, hspace=0.34)
        axes = [fig.add_subplot(inner[0, 0]), fig.add_subplot(inner[0, 1]),
                fig.add_subplot(inner[1, 0]), fig.add_subplot(inner[1, 1])]
        print(f"  {scenario}...")
        try:
            pie_data, valid_seeds = preproc.get_panel_a_data(
                exp_num, scenario, cluster_threshold, min_days, exp_name=exp_name)
            if valid_seeds:
                plot_fig4_pies(axes, pie_data)
                # plot_fig4_pies titles each cell (Peak/Burden/...), so the
                # scenario name goes above the whole 2x2 block instead of on
                # a cell, where it collided with "Peak".
                box = outer[0, col].get_position(fig)
                fig.text(box.x0 + box.width / 2, box.y1 + 0.055,
                         f"{scenario}  ({valid_seeds} seeds)", ha="center",
                         va="bottom", fontsize=7.5, fontweight="bold")
            else:
                raise ValueError("no valid seeds")
        except Exception as exc:
            for ax in axes:
                ax.axis("off")
            axes[0].text(0.5, 0.5, f"{scenario}\n(no data)", ha="center",
                         va="center", fontsize=6, color="0.45")
            print(f"    [warn] {scenario}: {exc}")

    add_global_pie_legend(fig)

    out = figure_path("3S", exp_name, exp_num, fmt)
    plt.savefig(out, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Generated {out}")


if __name__ == "__main__":
    args = parse_arguments()
    build_figure_3S(args.exp_num, args.exp_name, args.cluster_threshold,
                    args.min_days, args.format)
