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
# impact_long_shedders -- PIPELINE CONFIG
# ----------------------------------------------------------------------------
# Single place that knows HOW TO BUILD the parameter settings for every
# pipeline stage. impact_long_shedders_cal_1.py / _cal_2.py / _exp.py only
# RUN things: parse CLI args, call a build_*_settings() function below to get
# a ready-to-submit settings callable, dispatch it via run_experiment_script,
# and handle stage-specific I/O (reading calibration fits back, writing the
# frozen table, printing). They never construct simulation parameter dicts
# themselves -- if a parameter dict needs to change, this is the only file
# to touch.
#
# NSR sweep ranges are the one exception: they live in a user-editable JSON
# file under Data/00_Reference_parameters/ (read_nsr_ranges, below), since
# they get retuned between pipeline runs by inspecting the calibration-fit
# plots this pipeline produces.
# ============================================================================

import os
import json

import numpy as np
import pandas as pd

import simplicity.dir_manager as dm
import simplicity.settings_manager as sm

# =============================================================================
# SHARED SCENARIO DATA -- used by every stage
# =============================================================================
# One entry per scenario. `inf_duration_long` is the raw clinical long-shedding
# duration (days); the model parameter tau_3_long is DERIVED from it (see
# derive_tau_3_long) as tau_3_long = inf_duration_long - (tau_1+tau_2+tau_4).
# `control` has no long side, so inf_duration_long is unused (None).
# `beta_multiplier` scales that scenario's long-shedder daily transmission
# rate relative to standard individuals (1.0 = same daily rate, only
# duration differs -- see derive_scenario_params / sm.compute_r_long).
SCENARIOS = [
    {"name": "control",   "inf_duration_long": None,  "long_shedders_ratio": 0.00, "beta_multiplier": 1.0},
    {"name": "SOT",       "inf_duration_long": 63.0,  "long_shedders_ratio": 0.01, "beta_multiplier": 1.0},
    {"name": "HIV_low",   "inf_duration_long": 109.0, "long_shedders_ratio": 0.01, "beta_multiplier": 1.0},
    {"name": "HIV_high",  "inf_duration_long": 109.0, "long_shedders_ratio": 0.12, "beta_multiplier": 1.0},
    {"name": "edge_case", "inf_duration_long": 365.0, "long_shedders_ratio": 0.01, "beta_multiplier": 1.0},
]

TAU_ROUND = 3  # decimals for tau_3_long dict-key / group matching

# Experiment name stems (numbered with _#{exp_num} by run_experiment_script).
# Single source of truth so cal_1, cal_2 and the pipeline orchestrator can
# never disagree on where Stage 1's output lives.
LONG_NSR_EXP_NAME = "impact_long_shedders_calibration_lng_nsr"
STD_NSR_SWEEP_NAME = "impact_long_shedders_calibration_std_nsr"


# =============================================================================
# SHARED MATH -- tau_3_long / R_long derivation
# ----------------------------------------------------------------------------
# R_long comes from simplicity.settings_manager.compute_r_long: long
# shedders transmit at the same daily rate as standard individuals
# (beta_long == beta_standard * beta_multiplier), just for longer. See that
# function's docstring for the formula; single source of truth, shared with
# core.
# =============================================================================
def phase_offset(sp):
    """Sum of the non-long-shedding infection phases (tau_1+tau_2+tau_4)."""
    return sp["tau_1"] + sp["tau_2"] + sp["tau_4"]


def derive_tau_3_long(scenario, sp):
    """Corrected long-shedding duration: inf_duration_long - phase_offset(sp)
    for a long-shedder scenario, or sp's own tau_3_long default for control."""
    if scenario["long_shedders_ratio"] == 0.0:
        return float(sp["tau_3_long"])
    return float(scenario["inf_duration_long"]) - phase_offset(sp)


def derive_scenario_params(scenario, sp, R):
    """
    Derive the frozen long-side parameters for one production scenario.
    R is this experiment's own standard-population R (e.g.
    USER_FIXED_PARAMS['R']) -- NOT sp['R'], simplicity's generic default,
    which this pipeline never actually runs with.

    Returns a dict with: scenario_name, long_shedders_ratio, tau_3_long,
    R_long, sequence_long_shedders, is_long.
    """
    name = scenario["name"]
    ratio = scenario["long_shedders_ratio"]
    tau_3_long = derive_tau_3_long(scenario, sp)
    r_long = sm.compute_r_long(R, sp["tau_2"], sp["tau_3"], tau_3_long,
                                scenario["beta_multiplier"])

    if ratio == 0.0:
        # Control: no long shedders, so R_long is never actually used by the
        # simulation -- still derive it consistently rather than a stray default.
        return {
            "scenario_name": name,
            "long_shedders_ratio": 0.0,
            "tau_3_long": tau_3_long,
            "R_long": r_long,
            "sequence_long_shedders": False,
            "is_long": False,
        }

    return {
        "scenario_name": name,
        "long_shedders_ratio": float(ratio),
        "tau_3_long": tau_3_long,
        "R_long": r_long,
        "sequence_long_shedders": True,
        "is_long": True,
    }


def unique_long_taus(sp):
    """Sorted unique corrected tau_3_long values across long-shedder scenarios."""
    taus = {
        round(derive_tau_3_long(scenario, sp), TAU_ROUND)
        for scenario in SCENARIOS
        if scenario["long_shedders_ratio"] > 0.0
    }
    return sorted(taus)


# =============================================================================
# STAGE 1 (cal_1) -- isolated long-shedder NSR calibration
# ----------------------------------------------------------------------------
# A 100% long-shedder population, used ONLY to calibrate the long-shedder NSR
# in isolation -- deliberately different from USER_FIXED_PARAMS (Stage 2/3's
# mixed production population), not a duplicate of it.
# =============================================================================
CAL1_ISOLATED_FIXED_PARAMS = {
    "long_shedders_ratio": 1.0,
    "R": 1.1,
    "infected_individuals_at_start": 50,
    "final_time": 720,
    "sequence_long_shedders": True,
    # R_long deliberately absent: build_cal1_settings derives it per
    # tau_3_long group via sm.compute_r_long(R, ..., tau, multiplier).
}


def build_cal1_settings(seeds, ranges, R=None, multiplier=1.0):
    """
    Returns a zero-arg make_settings callable ready for run_experiment_script:
    one group per tau_3_long value (see unique_long_taus), each sweeping the
    same NSR grid, with R_long correctly scaled to that group's own duration
    (rather than a single fixed value shared by every group regardless of
    duration). R defaults to CAL1_ISOLATED_FIXED_PARAMS['R'] if not given.
    """
    sp = sm.read_standard_parameters_values()
    fixed_params = CAL1_ISOLATED_FIXED_PARAMS.copy()
    if R is not None:
        fixed_params["R"] = R

    def make_settings():
        nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()
        scenario_groups = [
            {
                'nucleotide_substitution_rate': nsr_values,
                'tau_3_long': tau,
                'R_long': sm.compute_r_long(
                    fixed_params["R"], sp["tau_2"], sp["tau_3"], tau, multiplier),
            }
            for tau in unique_long_taus(sp)
        ]
        varying_params = {'_scenario_groups': scenario_groups}
        return (varying_params, fixed_params, seeds)

    return make_settings


# =============================================================================
# STAGE 2 (cal_2) -- standard NSR calibrated in-context (mixed population)
# =============================================================================
# Fixed experiment overrides shared by cal_2 and exp (mirrors the population
# context of the impact experiment). Single source of truth -- previously
# duplicated verbatim in both scripts.
USER_FIXED_PARAMS = {
    "population_size": 1000,
    "infected_individuals_at_start": 50,
    "R": 1.03,
    "final_time": 1095,
    "IH_virus_emergence_rate": 0.1,
}

# Shared color palette for the calibration-fit plots (long-NSR per tau group,
# standard-NSR per scenario), keyed by scenario name so both plots look
# consistent with each other.
DEFAULT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]


def lookup_long_nsr(long_nsr_by_tau, tau_3_long):
    """Exact (rounded) match lookup. Raises on miss -- never silently guesses."""
    key = round(float(tau_3_long), TAU_ROUND)
    if key not in long_nsr_by_tau:
        raise KeyError(
            f"No calibrated long NSR for tau_3_long={key}. "
            f"Available keys: {sorted(long_nsr_by_tau.keys())}. "
            f"Check that the scenario tau matches the calibration experiment.")
    return long_nsr_by_tau[key]


def build_cal2_scenario_groups(nsr_ranges, long_nsr_by_tau, R, multiplier=None):
    """
    One group dict per scenario: its own nucleotide_substitution_rate sweep
    (list -> that group's local varying axis) plus its fixed tau_3_long/
    long_shedders_ratio/R_long/sequence_long_shedders/long-NSR overrides.

    long_nsr_by_tau is Stage 1's calibration-fit result (computed by cal_2.py
    from already-run simulation output, not parameter setup -- passed in
    rather than recomputed here). multiplier, if given, overrides every
    scenario's own beta_multiplier (SCENARIOS) uniformly.

    Returns (scenario_groups, scenarios_frozen).
    """
    sp = sm.read_standard_parameters_values()
    scenario_groups = []
    scenarios_frozen = []

    for scenario in SCENARIOS:
        if multiplier is not None:
            scenario = {**scenario, "beta_multiplier": multiplier}
        frozen = derive_scenario_params(scenario, sp, R)
        scenarios_frozen.append(frozen)

        name = frozen["scenario_name"]
        r = nsr_ranges[name]
        nsr_list = np.logspace(np.log10(r['min']), np.log10(r['max']), r['steps']).tolist()

        group = {
            "long_shedders_ratio": frozen["long_shedders_ratio"],
            "tau_3_long": frozen["tau_3_long"],
            "R_long": frozen["R_long"],
            "sequence_long_shedders": frozen["sequence_long_shedders"],
            "nucleotide_substitution_rate": [float(x) for x in nsr_list],
        }
        if frozen["is_long"]:
            group["nucleotide_substitution_rate_long"] = lookup_long_nsr(
                long_nsr_by_tau, frozen["tau_3_long"])

        scenario_groups.append(group)

    return scenario_groups, scenarios_frozen


def build_cal2_settings(scenario_groups, seeds, R=None):
    """Zero-arg make_settings callable ready for run_experiment_script,
    submitting every scenario's sweep as ONE combined experiment. R
    defaults to USER_FIXED_PARAMS['R'] if not given -- pass the SAME R used
    to build scenario_groups (build_cal2_scenario_groups) so the frozen
    R_long values stay consistent with what standard individuals actually
    run with."""
    fixed_params = USER_FIXED_PARAMS.copy()
    if R is not None:
        fixed_params["R"] = R

    def make_settings():
        varying_params = {'_scenario_groups': scenario_groups}
        return (varying_params, fixed_params, seeds)
    return make_settings


# =============================================================================
# STAGE 3 (exp) -- production runs
# ----------------------------------------------------------------------------
# Reuses USER_FIXED_PARAMS (Stage 2, above) as-is: the production population
# context must match the one cal_2 calibrated the standard NSR against.
# =============================================================================
def build_exp_scenario_settings(row, n_seeds):
    """
    Zero-arg make_settings callable for one row of the frozen calibration
    table (a single simulation point per scenario -- empty varying_params,
    see generate_experiment_settings: combinations = [()] when there are no
    varying params).
    """
    is_control = float(row["long_shedders_ratio"]) == 0.0

    fixed = USER_FIXED_PARAMS.copy()
    fixed.update({
        "R": float(row["R"]),
        "long_shedders_ratio": float(row["long_shedders_ratio"]),
        "tau_3_long": float(row["tau_3_long"]),
        "R_long": float(row["R_long"]),
        "nucleotide_substitution_rate": float(row["nucleotide_substitution_rate"]),
        "sequence_long_shedders": (not is_control),
    })

    # Long NSR: set only when a long side exists. For control the table value
    # may be empty/NaN; we must not write NaN into the parameters.
    if not is_control:
        long_nsr = row["nucleotide_substitution_rate_long"]
        if pd.isna(long_nsr):
            raise ValueError(
                f"Scenario '{row['scenario_name']}' has long_shedders_ratio>0 "
                f"but no nucleotide_substitution_rate_long in the table.")
        fixed["nucleotide_substitution_rate_long"] = float(long_nsr)

    def make_settings():
        return ({}, fixed.copy(), n_seeds)

    return make_settings


# =============================================================================
# NSR SWEEP RANGES -- user-editable reference file (drives Stage 1 & Stage 2
# grids)
# =============================================================================
NSR_RANGES_FILENAME = "impact_long_shedders_nsr_ranges.json"

# Single source for the grid resolution shared by every sweep below, so
# tuning it doesn't mean touching (and keeping in sync) six separate literals.
_DEFAULT_NSR_STEPS = 10

DEFAULT_NSR_RANGES = {
    "cal1_long_nsr": {"min": 0.0005, "max": 0.002, "steps": _DEFAULT_NSR_STEPS},
    "cal2_standard_nsr": {
        "control":   {"min": 1e-05, "max": 3e-04, "steps": _DEFAULT_NSR_STEPS},
        "SOT":       {"min": 1e-05, "max": 3e-04, "steps": _DEFAULT_NSR_STEPS},
        "HIV_low":   {"min": 1e-05, "max": 3e-04, "steps": _DEFAULT_NSR_STEPS},
        "HIV_high":  {"min": 1e-05, "max": 3e-04, "steps": _DEFAULT_NSR_STEPS},
        "edge_case": {"min": 1e-05, "max": 3e-04, "steps": _DEFAULT_NSR_STEPS},
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
# SLURM per-task resource request -- shared by every stage + the orchestrator
# ----------------------------------------------------------------------------
# Read by simplicity.runners.slurm.submit_simulations from the
# SIMPLICITY_SLURM_MEM/SIMPLICITY_SLURM_TIME env vars (same pattern as
# SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM) rather than a function
# argument, so it never has to touch the run_seeded_simulations(...) interface
# shared uniformly by the serial/multiprocessing/slurm runner modules.
#
# The baseline default (2G/2 days) is set exactly once, in
# simplicity.dir_manager, alongside every other simplicity-wide env var --
# not redefined here, so there's a single literal source for it. 2G was set
# from real sacct data on the isolated long-NSR calibration grid
# (population_size=1000, final_time=1200): MaxRSS ranged 0.29-0.45G across
# all 36 grid points -- real headroom over that, well below the old blanket
# 8G every job used to request regardless of actual need. Time limit was
# bumped from 1 day to 2: some high-R_long tau groups (e.g. edge_case) run
# well past 1 day of wall-clock time to reach final_time and were getting
# killed by the walltime
# mid-simulation -- see simplicity.runners.slurm.reconcile_terminated_tasks,
# which now detects and signals that case regardless, but giving cheap extra
# headroom avoids the kill (and the wasted compute) in the first place.
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
