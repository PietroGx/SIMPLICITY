#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Are long shedders advantaged only in how LONG they shed?

Figure 2 reads every clade-level difference between the scenarios as a
consequence of prolonged shedding. That reading holds only while the two host
types are identical in every other respect that affects transmission, so this
exercises the two manipulation checks on the parameters a production run
actually used:

  equal R                 a long shedder's expected number of secondary
                          infections equals a standard one's. R_long is set to
                          R_STANDARD and the per-day rate is R / infectious
                          duration, so a longer duration means a LOWER daily
                          rate, not more infections
  equal susceptibility    susceptibility_long = 1.0, so a long shedder is no
                          more or less likely to be infected in the first place

Read from `Data/<exp>/setup_data/nsr_calibration_table.csv` -- the frozen table
of a real run -- through `build_exp_scenario_settings`, the same builder the
pipeline calls, rather than from a parameter dict rebuilt here.

The daily rates are computed from the expressions at `simplicity/extrande.py`
lines 145 and 148, which are locals of `get_helpers` and cannot be imported.
Section C asserts those two lines still read as mirrored, so a change there
fails this check instead of silently invalidating it.

Read-only, runs no simulation. Run from the repo root:

    python -I tests/check_manipulation.py [--stamp 1]
"""
import os
import sys
import argparse

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts", "experiments_pipeline"))

import impact_long_shedders_config as C
import simplicity.population_model as pm
import simplicity.settings_manager as sm

PASS = [0]
FAIL = [0]


def ck(cond, msg):
    (PASS if cond else FAIL)[0] += 1
    print(("  ok   " if cond else "  FAIL ") + msg)


# the two expressions this check mirrors, with the file that owns them
MIRRORED = [
    ("simplicity/extrande.py",
     'beta_standard = parameters["R"] / tau_inf_standard'),
    ("simplicity/extrande.py",
     "beta_long = parameters['R_long']/ tau_inf_long"),
    ("simplicity/extrande.py",
     "tau_inf_standard = parameters['tau_2'] + parameters['tau_3']"),
    ("simplicity/extrande.py",
     "tau_inf_long = parameters['tau_2'] + parameters['tau_3_long']"),
    ("simplicity/population_model.py",
     "s = population.susceptibility_long"),
]


class _Stub:
    """Only what SIDR_propensities reads, so the real function can be called."""
    def __init__(self, n_std, n_long):
        self.infectious_standard = n_std
        self.infectious_long = n_long
        self.detectables_standard = 0
        self.detectables_long = 0
        self.infected = n_std + n_long


def section_a(rows):
    print("\nA. the frozen table a run actually used")
    for _, row in rows.iterrows():
        name = row["scenario_name"]
        _, fixed, _ = C.build_exp_scenario_settings(row, 1)()
        ratio = float(fixed["long_shedders_ratio"])
        R, R_long = float(fixed["R"]), float(fixed["R_long"])
        s = float(fixed["susceptibility_long"])
        if ratio == 0.0:
            ck(True, f"{name:10} no long shedders (ratio 0), R={R}, "
                     f"R_long={R_long} unused")
        else:
            ck(R_long == R,
               f"{name:10} R_long == R ({R_long} vs {R})")
        ck(s == 1.0, f"{name:10} susceptibility_long == 1.0 (got {s})")


def section_b(rows):
    print("\nB. equal R means a lower DAILY rate, not more infections")
    print(f"  {'scenario':10} {'tau_inf_std':>12} {'tau_inf_long':>13} "
          f"{'beta_std':>10} {'beta_long':>10} {'beta*tau std':>13} "
          f"{'beta*tau long':>14}")
    for _, row in rows.iterrows():
        name = row["scenario_name"]
        _, f, _ = C.build_exp_scenario_settings(row, 1)()
        if float(f["long_shedders_ratio"]) == 0.0:
            continue
        # the phase durations the run used: settings_manager's own defaults,
        # with anything the frozen table overrode winning
        f = dict(sm._STANDARD_VALUES_DEFAULTS, **f)
        ti_s = float(f["tau_2"]) + float(f["tau_3"])
        ti_l = float(f["tau_2"]) + float(f["tau_3_long"])
        b_s = float(f["R"]) / ti_s
        b_l = float(f["R_long"]) / ti_l
        print(f"  {name:10} {ti_s:12.2f} {ti_l:13.2f} {b_s:10.5f} "
              f"{b_l:10.5f} {b_s * ti_s:13.3f} {b_l * ti_l:14.3f}")
        ck(abs(b_s * ti_s - b_l * ti_l) < 1e-9,
           f"{name:10} expected secondary infections equal "
           f"({b_s * ti_s:.3f} vs {b_l * ti_l:.3f})")
        ck(b_l < b_s,
           f"{name:10} the long daily rate is LOWER "
           f"({b_l:.5f} < {b_s:.5f}), i.e. spread over {ti_l / ti_s:.1f}x longer")

        # the real propensity function, with equal numbers of each type
        props, _ = pm.SIDR_propensities(_Stub(100, 100), b_s, b_l, 0.0, 0.0, 0.0)
        d = dict((p[0], p[1]) for p in props)
        ck(abs(d["infection_long"] / d["infection_standard"] - b_l / b_s) < 1e-12,
           f"{name:10} SIDR_propensities scales the two arms by their own beta "
           f"({d['infection_long'] / d['infection_standard']:.3f})")


def section_c():
    print("\nC. the expressions this check mirrors are still there")
    for rel, frag in MIRRORED:
        src = open(os.path.join(ROOT, rel)).read()
        ck(frag in src, f"{rel}: {frag[:58]}")

    # susceptibility_long == 1.0 takes the uniform branch, which is what
    # "equally likely to be infected" means in code
    src = open(os.path.join(ROOT, "simplicity/population_model.py")).read()
    ck("if s == 1.0:" in src and
       "new_infected_index = population.rng4.choice(susceptibles_list)" in src,
       "population_model: s == 1.0 selects the bare uniform choice")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stamp", default="latest")
    a = ap.parse_args()
    from impact_long_shedders_config import resolve_stamp
    a.stamp = resolve_stamp(a.stamp, 'impact_long_shedders')

    table = os.path.join(ROOT, C.table_path(a.stamp))
    if not os.path.isfile(table):
        print(f"no frozen table at {table} -- run --stamp with a run that is "
              f"present locally")
        sys.exit(2)
    print(f"frozen table: {os.path.relpath(table, ROOT)}")
    rows = pd.read_csv(table)
    print(f"R_STANDARD = {C.R_STANDARD}")

    section_a(rows)
    section_b(rows)
    section_c()
    print("\n  Not a manipulation check, stated so this is not misread: the\n"
          "  substitution rate IS deliberately different. NSR_long is about 13x\n"
          "  the standard rate by construction (scripts/nature_plots/\n"
          "  fig2_preprocess_data.py:263), which is the point of the model --\n"
          "  what this check establishes is that TRANSMISSION is not also\n"
          "  advantaged.")
    print(f"\n{PASS[0]} passed, {FAIL[0]} failed")
    sys.exit(1 if FAIL[0] else 0)
