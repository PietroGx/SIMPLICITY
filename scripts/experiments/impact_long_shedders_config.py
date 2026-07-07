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
# impact_long_shedders -- SHARED CONFIG
# ----------------------------------------------------------------------------
# Single source of truth for the scenario definitions and fixed parameters
# shared by impact_long_shedders_cal_1.py, impact_long_shedders_cal_2.py and
# impact_long_shedders_exp.py.
#
# NSR sweep ranges are NOT hardcoded here: they live in a user-editable JSON
# file under Data/00_Reference_parameters/, since they get retuned between
# pipeline runs by inspecting the calibration-fit plots. This module only
# knows how to read that file (auto-creating it with sane defaults the first
# time) and write it back.
# ============================================================================

import os
import json

import simplicity.dir_manager as dm

# =============================================================================
# SCENARIOS
# =============================================================================
# One entry per scenario. `inf_duration_long` is the raw clinical long-shedding
# duration (days); tau_3_long is DERIVED from it (see derive_scenario_params)
# as tau_3_long = inf_duration_long - (tau_1 + tau_2 + tau_4). `control` has no
# long side, so inf_duration_long is unused (None).
SCENARIOS = [
    {"name": "control",   "inf_duration_long": None,  "long_shedders_ratio": 0.00},
    {"name": "SOT",       "inf_duration_long": 63.0,  "long_shedders_ratio": 0.01},
    {"name": "HIV_low",   "inf_duration_long": 109.0, "long_shedders_ratio": 0.01},
    {"name": "HIV_high",  "inf_duration_long": 109.0, "long_shedders_ratio": 0.12},
    {"name": "edge_case", "inf_duration_long": 365.0, "long_shedders_ratio": 0.01},
]

TAU_ROUND = 3  # decimals for tau_3_long dict-key / group matching

# Experiment name stems (numbered with _#{exp_num} by run_experiment_script).
# Single source of truth so cal_1, cal_2 and the pipeline orchestrator can
# never disagree on where Stage 1's output lives.
LONG_NSR_EXP_NAME = "impact_long_shedders_calibration_lng_nsr"
STD_NSR_SWEEP_NAME = "impact_long_shedders_calibration_std_nsr"

# Shared color palette for the calibration-fit plots (long-NSR per tau group,
# standard-NSR per scenario), keyed by scenario name so both plots look
# consistent with each other.
DEFAULT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

# Fixed experiment overrides shared by cal_2 and exp (mirrors the population
# context of the impact experiment). Single source of truth -- previously
# duplicated verbatim in both scripts.
USER_FIXED_PARAMS = {
    "population_size": 1000,
    "infected_individuals_at_start": 100,
    "R": 1.05,
    "final_time": 1095,
    "IH_virus_emergence_rate": 0.1,
}

# cal_1's isolated-calibration context: a 100% long-shedder population, used
# ONLY to calibrate the long-shedder NSR in isolation -- deliberately
# different from USER_FIXED_PARAMS (the mixed production population), not a
# duplicate of it. Still belongs in one place rather than inline in cal_1.py.
CAL1_ISOLATED_FIXED_PARAMS = {
    "long_shedders_ratio": 1.0,
    "R": 1.0,
    "R_long": 1.5,
    "infected_individuals_at_start": 100,
    "final_time": 1200,
    "sequence_long_shedders": True,
}


def phase_offset(sp):
    """Sum of the non-long-shedding infection phases (tau_1+tau_2+tau_4)."""
    return sp["tau_1"] + sp["tau_2"] + sp["tau_4"]


def derive_scenario_params(scenario, sp):
    """
    Derive the frozen long-side parameters for one scenario.

    Returns a dict with: scenario_name, long_shedders_ratio, tau_3_long,
    R_long, sequence_long_shedders, is_long.
    """
    name = scenario["name"]
    ratio = scenario["long_shedders_ratio"]

    if ratio == 0.0:
        # Control: no long shedders. Use standard defaults; long NSR unused.
        return {
            "scenario_name": name,
            "long_shedders_ratio": 0.0,
            "tau_3_long": float(sp["tau_3_long"]),
            "R_long": float(sp["R_long"]),
            "sequence_long_shedders": False,
            "is_long": False,
        }

    tau_3_long = float(scenario["inf_duration_long"]) - phase_offset(sp)
    R_long = tau_3_long / 7.0
    return {
        "scenario_name": name,
        "long_shedders_ratio": float(ratio),
        "tau_3_long": tau_3_long,
        "R_long": R_long,
        "sequence_long_shedders": True,
        "is_long": True,
    }


def unique_long_taus(sp):
    """Sorted unique corrected tau_3_long values across long-shedder scenarios."""
    taus = {
        round(derive_scenario_params(scenario, sp)["tau_3_long"], TAU_ROUND)
        for scenario in SCENARIOS
        if scenario["long_shedders_ratio"] > 0.0
    }
    return sorted(taus)


# =============================================================================
# NSR SWEEP RANGES -- user-editable reference file
# =============================================================================
NSR_RANGES_FILENAME = "impact_long_shedders_nsr_ranges.json"

DEFAULT_NSR_RANGES = {
    "cal1_long_nsr": {"min": 0.0001, "max": 0.002, "steps": 6},
    "cal2_standard_nsr": {
        "control":   {"min": 3e-05, "max": 5e-04, "steps": 6},
        "SOT":       {"min": 1e-05, "max": 1e-04, "steps": 6},
        "HIV_low":   {"min": 1e-05, "max": 1e-04, "steps": 6},
        "HIV_high":  {"min": 1e-05, "max": 1e-04, "steps": 6},
        "edge_case": {"min": 1e-05, "max": 1e-04, "steps": 6},
    },
}


def get_nsr_ranges_file_path():
    return os.path.join(dm.get_reference_parameters_dir(), NSR_RANGES_FILENAME)


def write_default_nsr_ranges():
    path = get_nsr_ranges_file_path()
    with open(path, "w") as f:
        json.dump(DEFAULT_NSR_RANGES, f, indent=4)
    print(f"Default NSR ranges written to {path}")


def read_nsr_ranges():
    """
    Read the user-editable NSR-ranges reference file, creating it with the
    defaults above on first use (never silently falls back to a different
    file, unlike simplicity.settings_manager.read_user_set_parameters_file).
    """
    path = get_nsr_ranges_file_path()
    try:
        with open(path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"[note] {path} not found; writing defaults.")
        write_default_nsr_ranges()
        return read_nsr_ranges()


# =============================================================================
# SLURM per-task resource request
# ----------------------------------------------------------------------------
# Read by simplicity.runners.slurm.submit_simulations from the
# SIMPLICITY_SLURM_MEM/SIMPLICITY_SLURM_TIME env vars (same pattern as
# SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM) rather than a function
# argument, so it never has to touch the run_seeded_simulations(...) interface
# shared uniformly by the serial/multiprocessing/slurm runner modules.
#
# The baseline default (2G/1 day) is set exactly once, in
# simplicity.dir_manager, alongside every other simplicity-wide env var --
# not redefined here, so there's a single literal source for it. 2G was set
# from real sacct data on the isolated long-NSR calibration grid
# (population_size=1000, final_time=1200): MaxRSS ranged 0.29-0.45G across
# all 36 grid points -- real headroom over that, well below the old blanket
# 8G every job used to request regardless of actual need. Time limit stays
# at 1 day: observed elapsed time (up to ~88 min) has even more headroom
# already, and unlike memory there's no cost to leaving it generous.
# =============================================================================
def add_slurm_resource_args(parser):
    """Add --slurm-mem/--slurm-time to an argparse parser (shared wiring so
    every stage script exposes the same two flags identically). Defaults come
    from the env vars dir_manager already set on import, not a second
    hardcoded copy of the same values."""
    parser.add_argument('--slurm-mem', type=str, default=os.environ["SIMPLICITY_SLURM_MEM"],
                        help="Per-task SLURM memory request "
                            f"(default: {os.environ['SIMPLICITY_SLURM_MEM']}).")
    parser.add_argument('--slurm-time', type=str, default=os.environ["SIMPLICITY_SLURM_TIME"],
                        help="Per-task SLURM time limit "
                            f"(default: {os.environ['SIMPLICITY_SLURM_TIME']}).")


def set_slurm_resource_env(mem, time):
    """Set the env vars simplicity.runners.slurm reads. Harmless no-op for
    the serial/multiprocessing runners, which never read them."""
    os.environ["SIMPLICITY_SLURM_MEM"] = mem
    os.environ["SIMPLICITY_SLURM_TIME"] = time
