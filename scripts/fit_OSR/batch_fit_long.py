# This file is part of SIMPLICITY
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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch OSR-fitting across generated experiments, then invert the AIC-best model to get NSR.

Hard-coded settings:
  - Reads experiment names from: scripts/fit_OSR/INDEX.csv
  - Uses target OSR = 0.001
  - Uses data_type = 'single_rates'
  - Writes CSV summary to: data/OSR_long_batch_fit_summary.csv

For each experiment_name:
  1) Build the OSR-vs-NSR dataframe (same logic as your per-experiment fitter)
  2) Fit 'log' and 'exp' regressors (AIC selection)
  3) Save fit plots (same plotting functions/kwargs as your original script)
  4) Read saved fit params via om.read_fit_results_csv(experiment_name, best_model)
  5) Invert that model with er.inverse_{best_model}_regressor(OSR, params) to get NSR*
  6) Append a row to CSV: experiment_name, best_model, NSR, OSR


"""

import os
from typing import Dict, List

import pandas as pd

# project imports
import simplicity.plots_manager as pm
import simplicity.output_manager as om
import simplicity.tuning.evolutionary_rate as er

# --------------------------- constants ---------------------------
INDEX_CSV = os.path.join('scripts', 'fit_OSR', 'INDEX.csv')
OUT_CSV = os.path.join('data', 'OSR_long_batch_fit_summary.csv')
TARGET_OSR = 1e-3
DATA_TYPE = 'single_rates'  # or 'combined_rate'
PARAMETER = 'nucleotide_substitution_rate'
MIN_SEQ_NUMBER = 30
MIN_SIM_LENGTH = 100  # note: passed as 'min_sim_lenght' kwarg to match existing API

# --------------------------- helpers ---------------------------
def read_experiment_names_from_index(index_csv_path: str) -> List[str]:
    if not os.path.exists(index_csv_path):
        raise FileNotFoundError(f"INDEX CSV not found: {index_csv_path}")
    df = pd.read_csv(index_csv_path)
    # Prefer experiment_name column; fall back to deriving from filename if needed
    if 'experiment_name' in df.columns:
        names = df['experiment_name'].astype(str).tolist()
    elif 'filename' in df.columns:
        names = [
            os.path.basename(fn).replace('.py', '').split('generate_data_')[-1]
            for fn in df['filename'].astype(str).tolist()
        ]
    else:
        raise KeyError("INDEX CSV must have 'experiment_name' or 'filename' column")
    # Keep stable order: by 'order' if present, else as-is
    if 'order' in df.columns:
        df['_expname_'] = names
        df_sorted = df.sort_values('order', kind='mergesort')
        names = df_sorted['_expname_'].tolist()
    return names


def build_df(experiment_name: str) -> pd.DataFrame:
    if DATA_TYPE == 'combined_rate':
        om.write_combined_OSR_vs_parameter_csv(
            experiment_name, PARAMETER, MIN_SEQ_NUMBER, MIN_SIM_LENGTH
        )
        df = om.read_combined_OSR_vs_parameter_csv(
            experiment_name, PARAMETER, MIN_SEQ_NUMBER, MIN_SIM_LENGTH
        )
    elif DATA_TYPE == 'single_rates':
        om.write_combined_OSR_vs_parameter_csv(
            experiment_name, PARAMETER, MIN_SEQ_NUMBER, MIN_SIM_LENGTH
        )
        om.write_OSR_vs_parameter_csv(
            experiment_name, PARAMETER, MIN_SEQ_NUMBER, MIN_SIM_LENGTH
        )
        df = om.read_OSR_vs_parameter_csv(
            experiment_name, PARAMETER, MIN_SEQ_NUMBER, MIN_SIM_LENGTH
        )
    else:
        raise ValueError("invalid DATA_TYPE; use 'single_rates' or 'combined_rate'")
    return df


def fit_models_and_plot(experiment_name: str, df: pd.DataFrame) -> Dict[str, float]:
    """Fit log/exp models, save plots, and return AICs."""
    aics: Dict[str, float] = {}

    # match your original kwargs (note the 'min_sim_lenght' spelling)
    kwargs = {"min_seq_number": MIN_SEQ_NUMBER, "min_sim_lenght": MIN_SIM_LENGTH}
    plot_fit = pm.plot_combined_OSR_fit if DATA_TYPE == 'combined_rate' else pm.plot_OSR_fit

    for m in ['log', 'exp']:
        print(f"\n##### Fitting {m} for {experiment_name} #####\n")
        fit_result = er.fit_observed_substitution_rate_regressor(
            experiment_name, df, m, None
        )
        aics[m] = fit_result.aic
        # Save plots (your plotters handle saving into the experiment directory)
        try:
            plot_fit(experiment_name, fit_result, m, **kwargs)
            print(f"Saved {m} fit plot for {experiment_name}.")
        except Exception as e:
            print(f"Plotting failed for {m} on {experiment_name}: {e}")
    return aics


def invert_best_model(experiment_name: str, best_model: str, target_osr: float) -> float:
    """Read fit params and invert using er.inverse_{best_model}_regressor(OSR, params)."""
    fit_params_df = om.read_fit_results_csv(experiment_name, best_model)
    params = fit_params_df.to_dict()
    inv_fn_name = f"inverse_{best_model}_regressor"
    inv_fn = getattr(er, inv_fn_name)  # both log & exp inverses exist
    nsr = float(inv_fn(target_osr, params))
    return nsr

# --------------------------- main ---------------------------
def main():
    exp_names = read_experiment_names_from_index(INDEX_CSV)

    rows = []
    for exp in exp_names:
        print(f"\n================== {exp} ==================")
        df = build_df(exp)
        aics = fit_models_and_plot(exp, df)
        best_model = min(aics.items(), key=lambda kv: kv[1])[0]
        nsr_star = invert_best_model(exp, best_model, TARGET_OSR)
        rows.append({
            'experiment_name': exp,
            'best_model': best_model,
            'NSR': nsr_star,
            'OSR': TARGET_OSR,
        })

    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    pd.DataFrame(rows, columns=['experiment_name', 'best_model', 'NSR', 'OSR']).to_csv(
        OUT_CSV, index=False
    )
    print(f"\nWrote summary to {OUT_CSV}")

if __name__ == "__main__":
    main()
