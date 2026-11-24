#!/usr/bin/env python3
'''Wall time and peak RSS of single simulations at a larger scale, per stage.

    python -I tests/benchmark_scale.py --list
    python -I tests/benchmark_scale.py --one exp_HIV_high --root <dir> --out <jsonl>
    bash tests/run_benchmark_scale.sh            # all of them, two at a time, in tmux

Each job is ONE simulation, built by the stage's own builder in
impact_long_shedders_config (build_exp_scenario_settings, build_cal1_settings,
build_cal2_settings) and run through simplicity.runme.run_experiment_phase with
the serial runner, in its own process. Only these are overridden:

    population_size, infected_individuals_at_start, final_time   on the dicts the
                                     builder returns (cal_1 keeps its own window)
    reservoir                        Population.__init__ is wrapped: it is not a
                                     parameter (simplicity/population.py:100)
    max_runtime                      lifted, so the guard cannot end a run early
    the NSR sweep                    cut to the one point being timed

The experiment phase uses the rates frozen by impact_long_shedders_#2. The `ref_`
jobs are the production configuration, the local half of a local-to-HPC speed
ratio. Outputs are deleted after each job; the row keeps what the run reached.
'''
import argparse
import json
import os
import resource
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path[:0] = [REPO, os.path.join(REPO, 'scripts', 'experiments_pipeline')]

NEW = dict(population_size=10000, infected_individuals_at_start=300,
           final_time=1825, reservoir=1_000_000)
REF = dict(population_size=3000, infected_individuals_at_start=90,
           final_time=1095, reservoir=100_000)
TABLE = os.path.join(REPO, 'Data', 'impact_long_shedders_#2', 'setup_data',
                     'nsr_calibration_table.csv')
SCENARIOS = ('control', 'SOT', 'HIV_low', 'HIV_high', 'edge_case')
COHORTS = ('SOT', 'HIV_low+HIV_high', 'edge_case')


def jobs():
    out = {f'ref_{s}': ('exp', s, None, REF) for s in ('control', 'HIV_high')}
    out.update({f'exp_{s}': ('exp', s, None, NEW) for s in SCENARIOS})
    for c in COHORTS:
        for end in ('min', 'max'):
            out[f'cal1_{c}_{end}'] = ('cal1', c, end, NEW)
    for end in ('min', 'max'):
        out[f'cal2_{end}'] = ('cal2', None, end, NEW)
    return out


def settings_for(stage, which, end, scale):
    import pandas as pd
    import impact_long_shedders_config as cfg
    size = {k: scale[k] for k in ('population_size', 'infected_individuals_at_start')}
    if stage == 'exp':
        cfg.set_experiment_output_env()
        row = pd.read_csv(TABLE).set_index('scenario_name').loc[which]
        row = dict(row, scenario_name=which)
        varying, fixed, _ = cfg.build_exp_scenario_settings(row, 1)()
        fixed.update(size, final_time=scale['final_time'])
    else:
        cfg.set_calibration_output_env()
        key = 'cal_long' if stage == 'cal1' else 'cal_std'
        r = cfg.NSR_RANGES[key]
        nsr = r[end]
        build = cfg.build_cal1_settings if stage == 'cal1' else cfg.build_cal2_settings
        varying, fixed, _ = build(1, dict(r, steps=2))()
        fixed.update(size)
        groups = varying['_scenario_groups']
        if stage == 'cal1':
            groups = [g for g in groups if g['name'] == cfg.calib_long_group(which)]
            assert len(groups) == 1, [g['name'] for g in varying['_scenario_groups']]
            groups[0]['nucleotide_substitution_rate_long'] = [nsr]
        else:
            fixed['final_time'] = scale['final_time']    # CAL2_FINAL_TIME = production's
            groups[0]['nucleotide_substitution_rate'] = [nsr]
        varying = {'_scenario_groups': groups}
    fixed['max_runtime'] = 10 ** 7
    return lambda: (varying, fixed, 1)


def run_one(label, root, out):
    stage, which, end, scale = jobs()[label]
    os.environ.setdefault('SIMPLICITY_SLURM_TIME', '2-00:00:00')
    import simplicity.dir_manager as dm
    import simplicity.settings_manager as sm
    import simplicity.population as popmod
    import simplicity.runners.serial as serial
    from simplicity.runme import run_experiment_phase

    seen = {}
    original = popmod.Population.__init__

    def patched(self, *a, **k):
        k['reservoir'] = scale['reservoir']
        original(self, *a, **k)
        seen.update(reservoir=self.reservoir, individuals=len(self.individuals),
                    size=self.size)
    popmod.Population.__init__ = patched

    work = os.path.join(root, label)
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    dm.set_data_dir(work)
    sm._data_dir = work
    name = f'test_benchmark_scale_{label}_#1'
    make = settings_for(stage, which, end, scale)
    t0 = time.time()
    rec = dict(label=label, stage=stage, which=which, nsr_end=end, **scale)
    try:
        run_experiment_phase(name, make, simplicity_runner=serial, append=False)
        rec['outcome'] = 'ok'
    except BaseException as exc:                                # noqa: BLE001
        rec['outcome'] = f'{type(exc).__name__}: {str(exc)[:120]}'
    rec['seconds'] = round(time.time() - t0, 1)
    rec['peak_rss_mb'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    rec.update(seen)
    try:
        import simplicity.output_manager as om
        sod = dm.get_simulation_output_dirs(name)[0]
        ssod = dm.get_seeded_simulation_output_dirs(sod)[0]
        rec.update(reached=float(om.read_final_time(ssod)),
                   infections=len(om.read_individuals_data(ssod)),
                   lineages=len(om.read_phylogenetic_data(ssod)),
                   stem=os.path.basename(sod))
        rec['output_mb'] = round(sum(os.path.getsize(os.path.join(d, f))
                                     for d, _, fs in os.walk(ssod) for f in fs) / 2 ** 20, 1)
    except Exception as exc:                                    # noqa: BLE001
        rec['read_error'] = str(exc)[:120]
    with open(out, 'a') as fh:
        fh.write(json.dumps(rec) + '\n')
    shutil.rmtree(work, ignore_errors=True)
    print(json.dumps(rec), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--list', action='store_true')
    p.add_argument('--one')
    p.add_argument('--root')
    p.add_argument('--out')
    a = p.parse_args()
    if a.list:
        print('\n'.join(jobs()))
        return 0
    run_one(a.one, a.root, a.out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
