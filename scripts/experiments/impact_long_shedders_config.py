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
