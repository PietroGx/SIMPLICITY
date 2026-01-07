#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate per-scenario NSR-run driver scripts from batch-fit results.

Reads:
  - Data/OSR_long_batch_fit_summary.csv   (columns: experiment_name,best_model,NSR,OSR)
  - scripts/fit_OSR/INDEX.csv             (columns: order,filename,experiment_name,long_evo_rate_f,tau_3_long,R_long,long_shedders_ratio,IH_virus_emergence_rate)

Writes:
  - scripts/long_shedders_experiments/<NN>_generate_data_<core>.py
    where <core> = experiment_name without numeric prefix and without "_OSR_fit"
  - scripts/long_shedders_experiments/INDEX_exp_scripts.csv
"""

import os
import re
import csv
import textwrap
import pandas as pd

SUMMARY_CSV = "Data/OSR_long_batch_fit_summary.csv"   # capital D
INDEX_CSV   = "scripts/fit_OSR/INDEX.csv"
OUT_DIR     = "scripts/long_shedders_experiments"

# -------- license header --------
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

# -------- experiment script comment (added above the script body) --------
SCENARIO_COMMENT = '''"""
  long_evo_rate_f = {long_evo_rate_f}
  tau_3_long = {tau_3_long}
  R_long = {R_long}
  long_shedders_ratio = {long_shedders_ratio}
  NSR (fitted) = {NSR}
"""
'''

# -------- template body --------
TEMPLATE = """#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from experiment_script_runner import run_experiment_script
import simplicity.settings_manager as sm
import argparse

# By request: experiment_name equals the script filename
experiment_name =  '{experiment_name}'

def user_set_experiment_settings():
    
    # --------- Specify parameter values manually -----------------------------
    
    # parameters value to get combinations from
    varying_params = {{}}
    # parameters to keep fixed (but different from standard_value) across combinations
    fixed_params = {{
        # constants across scenarios
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

def ensure_outdir():
    os.makedirs(OUT_DIR, exist_ok=True)

def canon(name: str) -> str:
    """
    Canonicalize for joining:
      - drop leading NN_ numeric prefix
      - drop trailing _#N suffix if present
    """
    s = str(name).strip()
    s = re.sub(r'^\d+_', '', s)
    s = re.sub(r'_#\d+$', '', s)
    return s

def core_from_index_name(index_exp_name: str) -> str:
    """
    Return the 'core' experiment name:
      - remove leading NN_ numeric prefix
      - remove trailing '_OSR_fit' if present
    """
    s = re.sub(r'^\d+_', '', str(index_exp_name).strip())
    s = re.sub(r'_OSR_fit$', '', s)
    return s

def write_script(row: dict) -> tuple[str, str]:
    """
    Create one script file and return (filename, generated_experiment_name).
    Filename:  NN_generate_data_<core>.py
    experiment_name inside file: NN_generate_data_<core>.py
    """
    order = int(row['order'])
    idx_name = str(row['experiment_name'])      # from INDEX (includes NN_ prefix)
    core = core_from_index_name(idx_name)       # no prefix, no _OSR_fit

    filename = f"{order:02d}_generate_data_{core}.py"
    experiment_name_value = f"{order:02d}_generate_data_{core}"            
    out_path = os.path.join(OUT_DIR, filename)

    nsr = float(row['NSR'])
    content = (
        LICENSE_HEADER
        + SCENARIO_COMMENT.format(
            long_evo_rate_f=row['long_evo_rate_f'],
            tau_3_long=row['tau_3_long'],
            R_long=row['R_long'],
            long_shedders_ratio=row['long_shedders_ratio'],
            NSR=f"{nsr:.8g}",
        )
        + textwrap.dedent(TEMPLATE).format(
            experiment_name=experiment_name_value,
            long_evo_rate_f=row['long_evo_rate_f'],
            tau_3_long=row['tau_3_long'],
            R_long=row['R_long'],
            long_shedders_ratio=row['long_shedders_ratio'],
            NSR=nsr,
        )
    )

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(content)
    return filename, experiment_name_value

def main():
    ensure_outdir()

    if not os.path.exists(SUMMARY_CSV):
        raise FileNotFoundError(f"Missing summary: {SUMMARY_CSV}")
    if not os.path.exists(INDEX_CSV):
        raise FileNotFoundError(f"Missing index: {INDEX_CSV}")

    df_sum = pd.read_csv(SUMMARY_CSV)
    df_idx = pd.read_csv(INDEX_CSV)

    # Canonical keys for a robust join (no prefix, no run-tag suffix)
    df_sum['key'] = df_sum['experiment_name'].map(canon)
    df_idx['key'] = df_idx['experiment_name'].map(canon)

    # Deduplicate summary per key (if multiple runs like _#1/_#2 exist)
    df_sum = df_sum.sort_values('experiment_name').drop_duplicates('key', keep='first')

    df = df_idx.merge(
        df_sum[['key', 'NSR', 'best_model', 'OSR']],
        on='key',
        how='inner',
        validate='one_to_one'
    )

    print(f"INDEX rows:   {len(df_idx)}")
    print(f"SUMMARY rows: {len(df_sum)}")
    print(f"JOIN rows:    {len(df)}")

    # Keep original ordering if present
    df = df.sort_values('order', kind='mergesort')

    rows = []
    for _, r in df.iterrows():
        fn, expname = write_script(r)
        rows.append({
            "order": int(r["order"]),
            "source_experiment_name": r['experiment_name'],  
            "generated_experiment_name": expname,            
            "filename": fn,
            "NSR": float(r['NSR']),
            "best_model": r.get('best_model', None),
            "OSR_target": r.get('OSR', None),
            "long_evo_rate_f": r['long_evo_rate_f'],
            "tau_3_long": r['tau_3_long'],
            "R_long": r['R_long'],
            "long_shedders_ratio": r['long_shedders_ratio'],
        })

    # Write the new index
    index_out = os.path.join(OUT_DIR, "INDEX_exp_scripts.csv")
    with open(index_out, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=[
            "order","source_experiment_name","generated_experiment_name","filename",
            "NSR","best_model","OSR_target",
            "long_evo_rate_f","tau_3_long","R_long","long_shedders_ratio"
        ])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Generated {len(rows)} scripts in {OUT_DIR}")
    print(f"Wrote index: {index_out}")

if __name__ == "__main__":
    main()

