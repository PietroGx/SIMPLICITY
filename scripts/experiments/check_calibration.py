#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Thin CLI wrapper around the shared long-shedder NSR calibration fit/plot
# (scripts/experiments/long_nsr_calibration_plot.py), kept as a standalone
# tool to manually re-fit a past impact_long_shedders_#N experiment.
# impact_long_shedders_cal_1.py calls the same shared function automatically
# at the end of its own run -- that call is the calibration stage, and it is
# what persists.
#
# READ-ONLY with respect to the experiment. Nothing is written under its
# directory: no fit-results CSV, no calibrated-long-NSR table (the Stage 1
# result impact_long_shedders_cal_2.py reads back), no per-seed OSR frame,
# and not the figure either. --out says where the figure goes, and is
# required for that reason. Rerunning with a different --target-osr shows
# what that target would give; it does not change what the next cal_2.py run
# under that exp_num will use.

import sys
import argparse

from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def main():
    parser = argparse.ArgumentParser(
        description="Re-fit the long-shedder NSR calibration for a run "
                    "already on disk, without writing into it.")
    parser.add_argument('--experiment-name', type=str, required=True,
                        help="Numbered experiment, e.g. "
                            "impact_long_shedders_#3. Since the pipeline "
                            "merged its phases into one experiment this is "
                            "the whole run; --group picks the cohorts, and "
                            "the default is its calib_long__* ones.")
    parser.add_argument('--group', action='append', default=None,
                        dest='groups',
                        help="Fit only this group; repeatable. Default: the "
                             "experiment's long-calibration groups.")
    parser.add_argument('--out', type=str, required=True,
                        help="Directory the calibration figure is written "
                             "to. Required: this tool writes nothing into "
                             "the experiment, its 05_Plots included.")
    parser.add_argument('--target-osr', type=float, default=0.00205)
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    args = parser.parse_args()

    try:
        plot_and_fit_long_nsr_calibration(
            args.experiment_name, args.target_osr, args.model, args.min_seq,
            args.min_len, groups=args.groups, refit_out_dir=args.out)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
