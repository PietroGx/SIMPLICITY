import argparse
import pandas as pd
import numpy as np
import os
import random
import ast

import simplicity.dir_manager as dm
import simplicity.sequencing as sq
import simplicity.settings_manager as sm
import simplicity.tuning.evolutionary_rate as er
from _scenarios import (scenario_names, scenario_output_dirs, scenario_label,
                        one_sod)
import simplicity.output_manager as om
from simplicity.intra_host_model import Host 

def parse_arguments():
    parser = argparse.ArgumentParser(description="SIMPLICITY Figure 1 Data Preprocessor")
    parser.add_argument('--engine', type=str, choices=['serial', 'multiprocessing', 'slurm'], default='serial')
    parser.add_argument('--force-recompute', action='store_true')
    return parser.parse_args()

def get_clinical_label(scenario):
    if scenario == "control": return "Control"
    if scenario == "SOT": return "SOT"
    if "HIV" in scenario: return "HIV"
    if scenario == "edge_case": return "Edge Case"
    return scenario.replace("_", " ").title()

def duration_groups(exp_name):
    """One scenario per DISTINCT clinical label, in config order.

    Panels A and B show infection-duration distributions, so HIV_low and
    HIV_high (which share tau_3_long) would draw the same curve twice. Keying
    on the label collapses them and lets edge_case in wherever the config
    defines it: Control, SOT, HIV, Edge Case.
    """
    seen, out = set(), []
    for name in scenario_names(exp_name):
        label = get_clinical_label(name)
        if label in seen:
            continue
        seen.add(label)
        out.append(name)
    return out


# Sequencing effort differs between the two clocks because the two draws are
# not the same operation, measured in tests/test_slope_pooling_and_effort.py on
# run #3:
#
#   D  global clock, ~2790 diagnosed sequences per seed. A 10% host draw still
#      leaves ~279, and between-seed variance is 99.7% of the total -- the draw
#      moves the mean slope by 0.05%. Free.
#   E  intra-host clock, read from individuals_data + phylogenetic_data, NOT
#      from the sequencing record. The analogous draw is over long-shedder
#      HOSTS, of which there are tens per seed; at 10% sampling noise becomes
#      91-94% of the variance and edge_case loses 28 of its 30 seeds. It also
#      describes a study that sequences 10% of long shedders and then recovers
#      each one's complete within-host history, which no study does.
#
# So the standard/global clock is read at 10% and the long-shedder/intra-host
# clock at 100%. These are not free parameters to tune.
GLOBAL_CLOCK_RATIO = 0.10
INTRAHOST_CLOCK_RATIO = 1.0
_RATIO_SEED = 20261020


def get_panel_a_data(exp_num=1, exp_name="impact_long_shedders"):
    scenarios = duration_groups(exp_name)
    df_list = []
    for scenario in scenarios:
        experiment_string = scenario_label(exp_name, scenario, exp_num)
        # outside the try on purpose: a guard caught by `except Exception`
        # below would turn a loud stop back into a silent skip.
        sod = one_sod(exp_name, scenario, exp_num)
        if sod is None: continue
        try:
            tau_1 = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_1')
            tau_2 = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_2')
            tau_4 = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_4')
            tau_3_standard = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3')
            if scenario == "control": tau_3_active = tau_3_standard
            else: tau_3_active = sm.get_parameter_value_from_simulation_output_dir(sod, 'tau_3_long')
            label = get_clinical_label(scenario)
            host = Host(tau_1=tau_1, tau_2=tau_2, tau_3=tau_3_active, tau_4=tau_4)
            time_max = 800
            step = 1
            p_inf, p_det, p_rec = host.data_plot_ih_solution(state=0, time=time_max, step=step)
            time_array = np.arange(0, time_max, step)
            temp_df = pd.DataFrame({'time': time_array[:len(p_inf)], 'p_infectious': p_inf, 'cohort': label})
            df_list.append(temp_df)
        except Exception: continue
    if not df_list: return pd.DataFrame(columns=['time', 'p_infectious', 'cohort'])
    return pd.concat(df_list, ignore_index=True)

def get_panel_b_data(exp_num=1, exp_name="impact_long_shedders"):
    scenarios = duration_groups(exp_name)
    df_list = []
    control_exp_string = scenario_label(exp_name, "control", exp_num)
    control_sod = one_sod(exp_name, "control", exp_num)
    if control_sod is None: return pd.DataFrame(columns=['duration', 'cohort'])
    try:
        all_control_ssods = dm.get_seeded_simulation_output_dirs(control_sod)
        rng = random.Random(42) 
        sampled_control_ssods = rng.sample(all_control_ssods, min(10, len(all_control_ssods)))
        matched_seed_names = [os.path.basename(ssod) for ssod in sampled_control_ssods]
    except Exception: return pd.DataFrame(columns=['duration', 'cohort'])

    for scenario in scenarios:
        experiment_string = scenario_label(exp_name, scenario, exp_num)
        sod = one_sod(exp_name, scenario, exp_num)      # outside the try
        if sod is None: continue
        try:
            all_ssods_for_scenario = dm.get_seeded_simulation_output_dirs(sod)
            for target_seed_name in matched_seed_names:
                target_ssod = next((s for s in all_ssods_for_scenario if os.path.basename(s) == target_seed_name), None)
                if target_ssod is None: continue
                df = om.read_individuals_data(target_ssod)
                if scenario == "control": filtered_df = df[df['type'] == 'standard'].copy()
                else: filtered_df = df[df['type'] == 'long_shedder'].copy()
                label = get_clinical_label(scenario)
                if 't_not_infected' in filtered_df.columns and 't_infection' in filtered_df.columns:
                    filtered_df['duration'] = filtered_df['t_not_infected'] - filtered_df['t_infection']
                    filtered_df = filtered_df[['duration']].dropna()
                    filtered_df['cohort'] = label
                    df_list.append(filtered_df)
        except Exception: continue
    if not df_list: return pd.DataFrame(columns=['duration', 'cohort'])
    return pd.concat(df_list, ignore_index=True)

# Real-world inputs live in Data/RealWorldData/, resolved from THIS file's
# location rather than the working directory. They used to be bare relative
# paths, so every one of them vanished silently if a figure was run from
# anywhere but the repo root -- and only panel C said so on the plot.
REAL_WORLD_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              os.pardir, os.pardir, 'Data', 'RealWorldData')


def _real_world_csv(filename, needed_by):
    """Read one curated input, or warn by name and return None."""
    path = os.path.normpath(os.path.join(REAL_WORLD_DIR, filename))
    if not os.path.exists(path):
        print(f"[fig1][warn] {needed_by}: missing {path}")
        return None
    return pd.read_csv(path)


def get_panel_c_data():
    """Published long-shedder durations: one row per patient, grouped by
    immunocompromise category. Source table also carries DOI/quotes, which
    the panel ignores."""
    df = _real_world_csv('literature_long_shedders_data.csv', 'panel C')
    empty = pd.DataFrame(columns=['duration', 'clinical_category'])
    if df is None:
        return empty
    if 'Days (Exact)' not in df.columns:
        print("[fig1][warn] panel C: no 'Days (Exact)' column; got "
              f"{list(df.columns)[:8]}")
        return empty
    df['duration'] = df['Days (Exact)']
    df = df.dropna(subset=['clinical_category', 'duration']).copy()
    print(f"[fig1] panel C: {len(df)} patients across "
          f"{df['clinical_category'].nunique()} categories")
    return df[['duration', 'clinical_category']]

def get_panel_de_data():
    df_d = pd.DataFrame()
    df_e = pd.DataFrame()
    loaded_d = _real_world_csv('data_fig1_D.csv', 'panel D')
    if loaded_d is not None:
        df_d = loaded_d
        if 'sampling_date' in df_d.columns: df_d['sampling_date'] = pd.to_datetime(df_d['sampling_date'])
    loaded_e = _real_world_csv('data_fig1_E.csv', 'panel E')
    if loaded_e is not None:
        df_e = loaded_e
        if 'sampling_date' in df_e.columns: df_e['sampling_date'] = pd.to_datetime(df_e['sampling_date'])
    return df_d, df_e

def get_model_global_clock(exp_num=1, exp_name="impact_long_shedders",
                           ratio=GLOBAL_CLOCK_RATIO):
    """Model GLOBAL clock: standard individuals in control, divergence from the
    outbreak root against absolute sequencing time. The model counterpart of
    the real standard-cohort root-to-tip regression."""
    control_exp = scenario_label(exp_name, "control", exp_num)
    rows = []
    sod = one_sod(exp_name, "control", exp_num)
    if sod is None:
        print(f"[fig1][warn] global clock: no output for {control_exp}")
        return pd.DataFrame(columns=['Sequencing_time', 'Distance_from_root', 'cohort'])
    warned = False
    for s, ssod in enumerate(dm.get_seeded_simulation_output_dirs(sod)):
        seq_path = os.path.join(ssod, 'sequencing_data_regression.csv')
        if not os.path.exists(seq_path):
            continue
        df = pd.read_csv(seq_path)
        col = next((c for c in ('individual_type', 'Individual_type') if c in df.columns), None)
        if col is None:
            continue
        sub = df[df[col] == 'standard']
        if sub.empty:
            continue
        if ratio < 1.0:
            # Host-level, via the surveillance model's own rule: a person is
            # sequenced or not, and then all of their genomes are read.
            if 'individual_index' not in sub.columns:
                if not warned:
                    print(f"[fig1][warn] global clock: no individual_index in "
                          f"{os.path.basename(seq_path)} -- this run predates "
                          f"it, so ratio {ratio} CANNOT be applied and the "
                          f"complete record is used instead. Rerun to enable.")
                    warned = True
            else:
                kept = sq.subsample(sub.to_dict('records'), ratio,
                                    seed=_RATIO_SEED + s)
                sub = pd.DataFrame(kept, columns=sub.columns)
                if sub.empty:
                    continue
        keep = sub[['Sequencing_time', 'Distance_from_root']].copy()
        keep['seed'] = os.path.basename(ssod)
        rows.append(keep)
    if not rows:
        print(f"[fig1][warn] global clock: no standard sequences in {control_exp}")
        return pd.DataFrame(columns=['Sequencing_time', 'Distance_from_root', 'cohort'])
    out = pd.concat(rows, ignore_index=True)
    out['cohort'] = 'Control'
    print(f"[fig1] global clock: {len(out)} points from {control_exp} "
          f"({out['seed'].nunique()} seeds, ratio {ratio})")
    return out


def get_model_intrahost_clock(exp_num=1, exp_name="impact_long_shedders",
                              max_seeds=None):
    """Model INTRA-HOST clock for long shedders, via
    evolutionary_rate.extract_ih_regression_data -- the same function cal_1 and
    the sanity plots use.

    This replaces a bespoke extraction that read sequencing_data.csv and
    differenced genome lengths. That file only holds DIAGNOSIS-path sequences,
    and since v2.4.28 production no longer censuses long shedders, so it
    carried 11-19 long-shedder rows per run: the panel reported SOT at
    0.00022 s/s/y where the real intra-host clock reads 0.00202.
    """
    rows = []
    for scenario in duration_groups(exp_name):
        if scenario == "control":
            continue
        exp_str = scenario_label(exp_name, scenario, exp_num)
        sod = one_sod(exp_name, scenario, exp_num)
        if sod is None:
            print(f"[fig1][warn] intra-host clock: no output for {exp_str}")
            continue
        label = get_clinical_label(scenario)
        ssods = dm.get_seeded_simulation_output_dirs(sod)
        if max_seeds:
            ssods = ssods[:max_seeds]
        n = 0
        for ssod in ssods:
            try:
                ih = er.extract_ih_regression_data(ssod)
            except Exception:
                continue
            if ih is None or ih.empty:
                continue
            ih = ih.copy()
            ih['cohort'] = label
            ih['seed'] = os.path.basename(ssod)     # INTRAHOST_CLOCK_RATIO = 1.0
            rows.append(ih)
            n += len(ih)
        print(f"[fig1] intra-host clock: {n} points for {label} ({scenario})")
    if not rows:
        return pd.DataFrame(columns=['Sequencing_time', 'Distance_from_root',
                                     'cohort', 'seed'])
    return pd.concat(rows, ignore_index=True)


def main():
    args = parse_arguments()
    os.makedirs('plot_data_cache', exist_ok=True)
if __name__ == "__main__": main()
