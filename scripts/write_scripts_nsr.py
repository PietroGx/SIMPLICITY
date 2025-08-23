"""
Generate NSR-run driver scripts from batch-fit results.

Reads:
  - data/OSR_long_batch_fit_summary.csv  (columns: experiment_name,best_model,NSR,OSR)
  - scripts/fit_OSR/INDEX.csv            (columns include: order,experiment_name,long_evo_rate_f,tau_3_long,R_long,long_shedders_ratio)

Writes:
  - scripts/nsr_runs/<NN>_generate_data_<NN_..._NSR_run>.py   (one per scenario)
  - scripts/nsr_runs/INDEX_NSR_runs.csv
"""

import os
import re
import csv
import textwrap
import pandas as pd

SUMMARY_CSV = "Data/OSR_long_batch_fit_summary.csv"
INDEX_CSV   = "scripts/fit_OSR/INDEX.csv"
OUT_DIR     = "scripts/nsr_runs"

# -------- license header exactly as you use --------
LICENSE_HEADER = """# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

# -------- template writer --------
TEMPLATE = """#!/usr/bin/env python3
# -*- coding: utf-8 -*-
\"\"\"

@author: pietro
   
Each simulation will be repeated n_seeds time with a different random seed.

The set of all simulations is what we call an experiment.
\"\"\"

from experiment_script_runner import run_experiment_script
import simplicity.settings_manager as sm
import argparse

experiment_name =  '{experiment_name}'

def user_set_experiment_settings():
    
    # --------- Specify parameter values manually -----------------------------
    
    # parameters value to get combinations from
    varying_params = {{
        # no sweeps for this run; NSR is fixed from the fit
    }}
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {{
        # global constants
        'population_size': 1000,
        'IH_virus_emergence_rate': 0.01,
        'infected_individuals_at_start': 100,
        'R': 1.03,
        'final_time': 365*3,

        # scenario parameters
        'long_evo_rate_f': {long_evo_rate_f},
        'tau_3_long': {tau_3_long},
        'R_long': {R_long},
        'long_shedders_ratio': {long_shedders_ratio},

        # fitted value from batch OSR fit
        'nucleotide_substitution_rate': {NSR:.8g}
    }}
    
    # ---------- OR import them from file -------------------------------------
    # varying_params = {{}}
    # filename = 'standard_values.json'
    # fixed_params = sm.read_user_set_parameters_file(filename)
    # -------------------------------------------------------------------------

    n_seeds = 100
    return (varying_params, fixed_params, n_seeds)

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Run script to generate IH lineages data")
    parser.add_argument('runner', type=str, help="runner")
    parser.add_argument('experiment_number', type=int, help="experiment number")
    args = parser.parse_args()
    # Run the script 
    run_experiment_script(args.runner, 
                          args.experiment_number, 
                          user_set_experiment_settings,
                          experiment_name)

if __name__ == "__main__":
    main()
"""

SCENARIO_COMMENT = '''"""
  long_evo_rate_f = {long_evo_rate_f}
  tau_3_long = {tau_3_long}
  R_long = {R_long}
  long_shedders_ratio = {long_shedders_ratio}
  NSR (fitted) = {NSR}
"""
'''

def ensure_outdir():
    os.makedirs(OUT_DIR, exist_ok=True)

def derive_prefix(exp_name: str) -> str:
    """Return 'NN' two-digit prefix from experiment_name like '01_longf1_...'."""
    m = re.match(r'^(\d{2})_', exp_name)
    if not m:
        # fallback: use 00 if missing
        return "00"
    return m.group(1)

def write_script(row: dict) -> str:
    exp_name = row['experiment_name']
    nsr = float(row['NSR'])
    longf = row['long_evo_rate_f']
    tau = row['tau_3_long']
    rlong = row['R_long']
    lsr = row['long_shedders_ratio']

    # new experiment name to avoid colliding with OSR_fit runs
    new_exp_name = f"{exp_name}_NSR_run"
    prefix = derive_prefix(exp_name)
    filename = f"{prefix}_generate_data_{new_exp_name}.py"
    out_path = os.path.join(OUT_DIR, filename)

    content = (
        LICENSE_HEADER
        + SCENARIO_COMMENT.format(
            long_evo_rate_f=longf, tau_3_long=tau, R_long=rlong,
            long_shedders_ratio=lsr, NSR=f"{nsr:.8g}"
          )
        + textwrap.dedent(TEMPLATE).format(
            experiment_name=new_exp_name,
            long_evo_rate_f=longf,
            tau_3_long=tau,
            R_long=rlong,
            long_shedders_ratio=lsr,
            NSR=nsr
          )
    )

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    return filename, new_exp_name

def main():
    ensure_outdir()

    # Load summary (NSR per experiment) and scenario params from INDEX
    if not os.path.exists(SUMMARY_CSV):
        raise FileNotFoundError(f"Missing summary CSV: {SUMMARY_CSV}")
    if not os.path.exists(INDEX_CSV):
        raise FileNotFoundError(f"Missing INDEX CSV: {INDEX_CSV}")

    df_sum = pd.read_csv(SUMMARY_CSV)               # experiment_name, NSR, ...
    df_idx = pd.read_csv(INDEX_CSV)                 # experiment_name, longf, tau, R_long, lsr, order,...

    # Join on experiment_name to get scenario params + NSR
    df = pd.merge(df_idx, df_sum[['experiment_name','NSR','best_model','OSR']], on='experiment_name', how='inner')

    rows = []
    for _, r in df.iterrows():
        filename, new_exp = write_script(r)
        rows.append({
            "order": int(r.get("order", 0)),
            "source_experiment_name": r['experiment_name'],
            "generated_experiment_name": new_exp,
            "filename": filename,
            "NSR": float(r['NSR']),
            "best_model": r.get('best_model', None),
            "OSR_target": r.get('OSR', None),
            "long_evo_rate_f": r['long_evo_rate_f'],
            "tau_3_long": r['tau_3_long'],
            "R_long": r['R_long'],
            "long_shedders_ratio": r['long_shedders_ratio'],
        })

    # Write an index for the generated NSR-run scripts
    index_out = os.path.join(OUT_DIR, "INDEX_NSR_runs.csv")
    with open(index_out, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[
            "order","source_experiment_name","generated_experiment_name","filename",
            "NSR","best_model","OSR_target",
            "long_evo_rate_f","tau_3_long","R_long","long_shedders_ratio"
        ])
        writer.writeheader()
        # keep order if present
        rows_sorted = sorted(rows, key=lambda x: x["order"])
        writer.writerows(rows_sorted)

    print(f"Generated {len(rows)} scripts in {OUT_DIR}")
    print(f"Wrote index: {index_out}")

if __name__ == "__main__":
    main()
