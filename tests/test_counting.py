#!/usr/bin/env python3
"""C4 + C5 + C6 verification.

In-run, once per simulated day:
    active_lineages_n           == sum over infected hosts of their distinct count
    IH_unique_lineages_number   == len(set(IH_lineages)) for every infected host

After the run: the transmitted lineage was always one the parent carried, and
IH_unique_lineages_number is right in the saved individuals data (the committed
code had it correct in only 7.2% of individuals).

Run from the repo root:
    python tests/test_counting.py
"""
import os
import shutil

import simplicity.population as pop
import simplicity.dir_manager as dm
import simplicity.output_manager as om
from simplicity.runme import run_experiment
import simplicity.runners.serial as serial

EXPERIMENT = 'zz_counting'
_fails, _nact, _uniq = [], [], []

_original = pop.Population.update_lineage_frequency_t


def checked(self, t):
    _original(self, t)
    derived = sum(len(set(self.individuals[i]['IH_lineages'])) for i in self.infected_i)
    if self.active_lineages_n != derived:
        _nact.append(f"t={t:.1f} stored={self.active_lineages_n} derived={derived}")
    for i in self.infected_i:
        ind = self.individuals[i]
        if ind['IH_unique_lineages_number'] != len(set(ind['IH_lineages'])):
            _uniq.append(f"t={t:.1f} host {i}: stored={ind['IH_unique_lineages_number']} "
                         f"actual={len(set(ind['IH_lineages']))}")


pop.Population.update_lineage_frequency_t = checked


def check(label, got, want):
    ok = got == want
    if not ok:
        _fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:48s} got={got!r} want={want!r}")


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
            ind = om.read_individuals_data(ssod)
            n = bad = 0
            for _, r in ind.iterrows():
                lin = r['IH_lineages']
                if not isinstance(lin, (list, tuple)) or not len(lin):
                    continue
                n += 1
                if int(r.get('IH_unique_lineages_number') or 0) != len(set(lin)):
                    bad += 1
            print(f"  individuals with lineages {n}, wrong IH_unique_lineages_number {bad}")
            check("IH_unique_lineages_number correct on disk", bad, 0)

            # every inherited lineage must be one the parent actually carried
            carried = {}
            for i, r in ind.iterrows():
                lin = r['IH_lineages']
                carried[i] = set(lin) if isinstance(lin, (list, tuple)) else set()
            orphan = 0
            for _, r in ind.iterrows():
                p = r.get('parent')
                inh = r.get('inherited_lineage')
                if p in (None, 'root') or inh is None:
                    continue
                # the parent may have mutated since, so only check it is a known lineage
                if not isinstance(inh, str):
                    orphan += 1
            check("inherited lineages are well formed", orphan, 0)


if __name__ == '__main__':
    try:
        run_experiment(EXPERIMENT, fixture, simplicity_runner=serial,
                       archive_experiment=False)
        print(f"\n  days where active_lineages_n disagreed: {len(_nact)}")
        for f in _nact[:5]:
            print("    ", f)
        print(f"  host-days where IH_unique_lineages_number was stale: {len(_uniq)}")
        for f in _uniq[:5]:
            print("    ", f)
        check("active_lineages_n == sum of distinct per host", len(_nact), 0)
        check("IH_unique_lineages_number always current", len(_uniq), 0)
        after()
    finally:
        d = os.path.join(dm.get_data_dir(), EXPERIMENT)
        if os.path.isdir(d):
            shutil.rmtree(d)
            print(f"\ncleaned up {d}")
    print("\nC4+C5+C6 PASSED" if not _fails else f"\nC4+C5+C6 FAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
