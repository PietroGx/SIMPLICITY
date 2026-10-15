#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Thin CLI wrapper around the shared long-shedder NSR calibration fit/plot
# (scripts/experiments/long_nsr_calibration_plot.py), kept as a standalone
# tool to manually rerun the fit/plot for a past
# impact_long_shedders_calibration_lng_nsr_#N experiment.
# impact_long_shedders_cal_1.py calls the same shared function automatically
# at the end of its own run.
#
# NOT side-effect-free: this also overwrites that experiment's saved
# calibrated-long-NSR file (write_calibrated_long_nsr), which
# impact_long_shedders_cal_2.py reads back as its Stage 1 result. Rerunning
# this with a different --target-osr changes what the NEXT cal_2.py run
# under that exp_num will use, not just the plot.

import sys
import argparse

from long_nsr_calibration_plot import plot_and_fit_long_nsr_calibration


def main():
    parser = argparse.ArgumentParser(
        description="Rerun the long-shedder NSR calibration fit/plot for a "
                    "run already on disk.")
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
    parser.add_argument('--target-osr', type=float, default=0.00205)
    parser.add_argument('--model', type=str, default='exp',
                        choices=['lin', 'log', 'exp', 'tan'])
    parser.add_argument('--min-seq', type=int, default=30)
    parser.add_argument('--min-len', type=int, default=100)
    args = parser.parse_args()

    try:
        plot_and_fit_long_nsr_calibration(
            args.experiment_name, args.target_osr, args.model, args.min_seq,
            args.min_len, groups=args.groups)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
