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
# impact_long_shedders_unbound -- PIPELINE CONFIG
# ----------------------------------------------------------------------------
# Second pipeline, run alongside impact_long_shedders (which is untouched).
# The difference is what calibration is allowed to compensate for:
#
#   impact_long_shedders  calibrates the standard NSR IN CONTEXT, per scenario,
#                         inside the mixed population -- so the standard rate
#                         is pushed down to hold the population clock at
#                         target_osr_std no matter how much divergence long
#                         shedders inject. Run #7 drove HIV_high's standard
#                         NSR to 9.3e-7, 150x below control, to do it.
#
#   impact_long_shedders_ ("unbound") calibrates each cohort ALONE and then
#   unbound              lets them loose together. NSR_long comes from a 100%
#                        long-shedder population, the standard NSR from a 0%
#                        long-shedder population, and production applies both
#                        unchanged. The population-level clock is then an
#                        OUTPUT -- free to rise above target_osr_std -- rather
#                        than a constraint. The two calibrations are
#                        independent: stage 2 does not read stage 1.
#
# Measurement code is deliberately NOT duplicated: stage 1 reuses
# long_nsr_calibration_plot.plot_and_fit_long_nsr_calibration and therefore
# evolutionary_rate.extract_ih_regression_data, exactly as the other pipeline
# does. Parameter construction lives here, in the config, as it does there.
# ============================================================================

import numpy as np
import pandas as pd

import simplicity.settings_manager as sm

# Shared with the bound pipeline on purpose. The tau math must agree, and the
# production population context must be identical or the two pipelines are not
# comparable. Everything below that differs is defined here.
from impact_long_shedders_config import (
    TAU_ROUND, phase_offset, derive_tau_3_long, DEFAULT_COLORS,
    USER_FIXED_PARAMS, CAL1_ISOLATED_FIXED_PARAMS,
    add_slurm_resource_args, set_slurm_resource_env, print_fixed_params,
)

# =============================================================================
# SCENARIOS
# =============================================================================
# edge_case is back. It was dropped in v2.4.24 for raw-Hamming saturation at
# "~0.34 substitutions/site", but that figure dates from before the intra-host
# measurement was corrected, when the calibration was returning NSR_long ~0.116
# -- about 100x too high. At the rate the corrected pipeline actually
# calibrates (~1.1e-3), a 350-day infection accumulates 0.0020 subs/site, some
# 50x below the saturation threshold. Its 0.1% prevalence puts exactly one
# carrier in the initial 1000, which is the point: what does a single
# super-long shedder do to an outbreak?
SCENARIOS = [
    {"name": "control",   "inf_duration_long": None,  "long_shedders_ratio": 0.000, "R_long": None, "susceptibility_long": 1.0},
    {"name": "SOT",       "inf_duration_long": 63.0,  "long_shedders_ratio": 0.010, "R_long": 1.1,  "susceptibility_long": 1.0},
    {"name": "HIV_low",   "inf_duration_long": 109.0, "long_shedders_ratio": 0.010, "R_long": 1.1,  "susceptibility_long": 1.0},
    {"name": "HIV_high",  "inf_duration_long": 109.0, "long_shedders_ratio": 0.120, "R_long": 1.1,  "susceptibility_long": 1.0},
    {"name": "edge_case", "inf_duration_long": 365.0, "long_shedders_ratio": 0.001, "R_long": 1.1,  "susceptibility_long": 1.0},
]

LONG_NSR_EXP_NAME = "impact_long_shedders_unbound_cal_long"
STD_NSR_EXP_NAME = "impact_long_shedders_unbound_cal_std"
PROD_EXP_NAME = "impact_long_shedders_unbound"
SETUP_DIR_TEMPLATE = "Data/impact_long_shedders_unbound_setup_data_#{exp_num}"
TABLE_FILENAME = "nsr_calibration_table.csv"


def derive_scenario_params(scenario, sp):
    """Frozen long-side parameters for one production scenario. Mirrors the
    bound pipeline's version -- same tau derivation, same directly-specified
    R_long, same diagnosis-only sequencing (sequence_long_shedders False)."""
    name = scenario["name"]
    ratio = scenario["long_shedders_ratio"]
    tau_3_long = derive_tau_3_long(scenario, sp)
    is_long = ratio > 0.0

    return {
        "scenario_name": name,
        "long_shedders_ratio": float(ratio),
        "susceptibility_long": float(scenario["susceptibility_long"]),
        "tau_3_long": tau_3_long,
        "R_long": float(scenario["R_long"]) if is_long else 0.0,
        "sequence_long_shedders": False,
        "is_long": is_long,
    }


def unique_long_tau_r_long_pairs(sp):
    """Sorted unique (tau_3_long, R_long) pairs across long-shedder scenarios.
    HIV_low/HIV_high share one; SOT and edge_case each add their own, so
    stage 1 fits three groups here against the bound pipeline's two."""
    pairs = {
        (round(derive_tau_3_long(s, sp), TAU_ROUND),
         round(float(s["R_long"]), TAU_ROUND))
        for s in SCENARIOS if s["long_shedders_ratio"] > 0.0
    }
    return sorted(pairs)


def lookup_long_nsr(long_nsr_by_group, tau_3_long, r_long):
    """Exact (rounded) match on (tau_3_long, R_long). Raises on miss rather
    than guessing."""
    key = (round(float(tau_3_long), TAU_ROUND), round(float(r_long), TAU_ROUND))
    if key not in long_nsr_by_group:
        raise KeyError(
            f"No calibrated long NSR for (tau_3_long, R_long)={key}. "
            f"Available: {sorted(long_nsr_by_group)}.")
    return long_nsr_by_group[key]


# =============================================================================
# STAGE 1 -- long shedders alone
# =============================================================================
# The window is PER GROUP, not global. A group needs room for a few complete
# infection cycles, and infection duration ranges from 63 days (SOT) to 365
# (edge_case). SOT and HIV fit several cycles into the standard 365-day window;
# edge_case does not fit even one, and its census only fires on RECOVERED long
# shedders, so a 365-day window would starve the min_seq gate -- HPC test #1 in
# BACKLOG found exactly that and confirmed it clears at 3 years. The rule below
# gives 365 to SOT and HIV and 3x365 to edge_case, without hardcoding a
# scenario name: any new duration gets a window that fits it.
UNBOUND_CAL1_FINAL_TIME = 365
UNBOUND_CAL1_MIN_CYCLES = 3


def cal1_final_time(tau_3_long, sp):
    """Simulated window for one cal_1 group: the standard one, unless the
    group's own infection duration needs more to complete
    UNBOUND_CAL1_MIN_CYCLES cycles."""
    inf_duration = tau_3_long + phase_offset(sp)
    needed = UNBOUND_CAL1_MIN_CYCLES * inf_duration
    return int(max(UNBOUND_CAL1_FINAL_TIME,
                   np.ceil(needed / UNBOUND_CAL1_FINAL_TIME) * UNBOUND_CAL1_FINAL_TIME))


def build_cal1_settings(seeds, ranges, R=None, ih_virus_emergence_rate=None):
    """One group per (tau_3_long, R_long) pair, each sweeping the same
    NSR_long grid in a 100% long-shedder population, each with its own
    simulated window (see cal1_final_time)."""
    sp = sm.read_standard_parameters_values()
    fixed_params = CAL1_ISOLATED_FIXED_PARAMS.copy()
    fixed_params["final_time"] = UNBOUND_CAL1_FINAL_TIME
    if R is not None:
        fixed_params["R"] = R
    fixed_params["IH_virus_emergence_rate"] = (
        USER_FIXED_PARAMS['IH_virus_emergence_rate']
        if ih_virus_emergence_rate is None else ih_virus_emergence_rate)

    def make_settings():
        nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()
        scenario_groups = [
            {
                'nucleotide_substitution_rate_long': nsr_values,
                'tau_3_long': tau,
                'R_long': r_long,
                # per-group scalar: generate_experiment_settings applies group
                # overrides after fixed_params, so this wins for this group only
                'final_time': cal1_final_time(tau, sp),
            }
            for tau, r_long in unique_long_tau_r_long_pairs(sp)
        ]
        return ({'_scenario_groups': scenario_groups}, fixed_params, seeds)

    return make_settings


# =============================================================================
# STAGE 2 -- standard individuals alone
# =============================================================================
# A 0%-long-shedder population at production settings: the same context every
# scenario runs in, minus the long side. One sweep, one calibrated standard
# NSR, shared by every scenario -- there is nothing scenario-specific left to
# fit once calibration stops compensating for the long side.
#
# final_time 365, matching the standard cal_1 window. Production runs 1095, so
# control will not sit exactly on target_osr_std -- but in this pipeline that
# is not what the standard rate is for: it defines "the rate at which standard
# individuals mutate", measured in a clean population, and the population clock
# is deliberately an output. Worth watching control's production clock all the
# same, since it is the reference every elevation is measured against.
#
# sequencing_rate 1.0, as in the bound pipeline since v2.4.29: pure observation
# (rng6 drives nothing else, sequences are never read back), so it lowers the
# variance of the estimate and nothing else.
UNBOUND_CAL2_FINAL_TIME = 365
UNBOUND_CAL2_SEQUENCING_RATE = 1.0


def build_cal2_settings(seeds, ranges, R=None, ih_virus_emergence_rate=None):
    """Single standard-only sweep: long_shedders_ratio 0, standard NSR swept."""
    fixed_params = USER_FIXED_PARAMS.copy()
    fixed_params.update({
        "long_shedders_ratio": 0.0,
        "sequence_long_shedders": False,
        "final_time": UNBOUND_CAL2_FINAL_TIME,
        "sequencing_rate": UNBOUND_CAL2_SEQUENCING_RATE,
    })
    if R is not None:
        fixed_params["R"] = R
    if ih_virus_emergence_rate is not None:
        fixed_params["IH_virus_emergence_rate"] = ih_virus_emergence_rate

    def make_settings():
        nsr_values = np.geomspace(ranges['min'], ranges['max'], ranges['steps']).tolist()
        varying_params = {'nucleotide_substitution_rate': nsr_values}
        return (varying_params, fixed_params, seeds)

    return make_settings


# =============================================================================
# STAGE 3 -- production
# =============================================================================
def build_exp_scenario_settings(row, n_seeds):
    """One simulation point per frozen-table row. Both rates are applied as
    calibrated; nothing here compensates for their interaction."""
    is_control = float(row["long_shedders_ratio"]) == 0.0

    fixed = USER_FIXED_PARAMS.copy()
    fixed.update({
        "R": float(row["R"]),
        "IH_virus_emergence_rate": float(row["IH_virus_emergence_rate"]),
        "long_shedders_ratio": float(row["long_shedders_ratio"]),
        "susceptibility_long": float(row["susceptibility_long"]),
        "tau_3_long": float(row["tau_3_long"]),
        "R_long": float(row["R_long"]),
        "nucleotide_substitution_rate": float(row["nucleotide_substitution_rate"]),
        "sequence_long_shedders": False,
    })

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
# NSR SWEEP RANGES
# =============================================================================
# cal_long reaches a decade lower than the bound pipeline's 1e-4..1e-2:
# edge_case sheds for 350 days against HIV's 109, so it has far longer to
# accumulate and should calibrate below the others.
# cal_std brackets run #7's control, which landed at 1.38e-4 -- the closest
# thing to this stage the bound pipeline has.
_NSR_STEPS = 10

NSR_RANGES = {
    "cal_long": {"min": 1e-05, "max": 1e-02, "steps": _NSR_STEPS},
    "cal_std":  {"min": 1e-05, "max": 1e-03, "steps": _NSR_STEPS},
}

TABLE_COLUMNS = [
    "scenario_name", "tau_3_long", "long_shedders_ratio", "susceptibility_long",
    "R", "IH_virus_emergence_rate", "R_long",
    "nucleotide_substitution_rate_long", "nucleotide_substitution_rate",
    "target_osr_long", "model_type", "long_calib_source", "std_calib_source",
]
