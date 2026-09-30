#!/usr/bin/env python3
"""C3 + C2 verification.

In-run, once per simulated day: the frequencies sum to 1.
After the run: they still sum to 1 in the CSV, and long shedders' lineages
appear in it (today they never could).

Run from the repo root:
    python tests/test_frequency.py
"""
import os
import shutil

import simplicity.population as pop
import simplicity.dir_manager as dm
import simplicity.output_manager as om
from simplicity.runme import run_experiment
import simplicity.runners.serial as serial

EXPERIMENT = 'zz_frequency'
_fails = []
_inrun = []

_original = pop.Population.update_lineage_frequency_t


def checked(self, t):
    _original(self, t)
    today = [r for r in self.lineage_frequency if r['Time_sampling'] == t]
    if today:
        s = sum(r['Frequency_at_t'] for r in today)
        if abs(s - 1.0) > 1e-9:
            _inrun.append(f"t={t:.1f} sum={s!r}")


pop.Population.update_lineage_frequency_t = checked


def check(label, got, want):
    ok = got == want
    if not ok:
        _fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:44s} got={got!r} want={want!r}")


def fixture():
    return ({'long_shedders_ratio': [0.10]},
            {'population_size': 300,
             'nucleotide_substitution_rate': 0.0005,
             'infected_individuals_at_start': 20,
             'R': 1.2, 'R_long': 3,
             'diagnosis_rate_standard': 0.1, 'diagnosis_rate_long': 0.01,
             'final_time': 150},
            1)


def after():
    for sod in dm.get_simulation_output_dirs(EXPERIMENT):
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            lf = om.read_lineage_frequency(ssod)
            worst = (lf.groupby('Time_sampling')['Frequency_at_t'].sum() - 1.0).abs().max()
            print(f"  worst deviation from 1 across {lf['Time_sampling'].nunique()} days: {worst:.2e}")
            check("frequencies sum to 1 in the CSV", bool(worst < 1e-9), True)

            ind = om.read_individuals_data(ssod)
            carried = set()
            for _, row in ind[ind['type'] == 'long_shedder'].iterrows():
                if isinstance(row['IH_lineages'], (list, tuple)):
                    carried |= set(row['IH_lineages'])
            present = carried & set(lf['Lineage_name'])
            print(f"  long-shedder lineages carried {len(carried)}, present in CSV {len(present)}")
            check("long shedders reach the frequency table", bool(carried and present), True)


if __name__ == '__main__':
    try:
        run_experiment(EXPERIMENT, fixture, simplicity_runner=serial,
                       archive_experiment=False)
        print(f"\n  in-run days where the sum was not 1: {len(_inrun)}")
        for f in _inrun[:5]:
            print("    ", f)
        check("frequencies sum to 1 every day", len(_inrun), 0)
        after()
    finally:
        d = os.path.join(dm.get_data_dir(), EXPERIMENT)
        if os.path.isdir(d):
            shutil.rmtree(d)
            print(f"\ncleaned up {d}")
    print("\nC3+C2 PASSED" if not _fails else f"\nC3+C2 FAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
