#!/usr/bin/env python3
'''
Does the incremental consensus give the same answer as rescanning the snapshot?

c.get_consensus rebuilds from the whole snapshot every 5 simulated days, so its
cost grows with elapsed days x circulating lineages and the list itself grows
without bound. ConsensusAccumulator folds each entry in once and evaluates in
O(positions). The maths is exact; this checks the implementation is too, and
measures what it saves.

Four unit checks, then a LIVE check: a real simulation runs with both paths
computed at every rebuild on identical data -- the old result drives the
simulation, so nothing about the run changes -- comparing the matrix, the
column distributions, delta, the consensus sequence, and the time each took.

    python tests/test_consensus_accumulator.py [--final-time 1095]

Nothing in the model is wired to the accumulator yet. That is the next step,
after this passes.
'''
import argparse
import os
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import numpy as np

import simplicity.phenotype.consensus as c
import simplicity.phenotype.distance as dis
import simplicity.population as pop
import simplicity.dir_manager as dm
from simplicity.runme import run_experiment
import simplicity.runners.serial as serial

REF = dis.reference
failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def alt(position, k=0):
    '''A base that differs from the reference at this position.'''
    return [b for b in 'ACGT' if b != REF[position]][k]


def both(snapshot, t):
    '''(old, new) consensus triples from the same entries.'''
    acc = c.ConsensusAccumulator()
    for genome, freq, t_entry in snapshot:
        acc.add(genome, freq, t_entry)
    return c.get_consensus(snapshot, t), acc.consensus(t)


def compare(old, new, label, tol=1e-9):
    chat_o, p_o, delta_o = old
    chat_n, p_n, delta_n = new
    check(f'{label}: consensus sequence identical', chat_o == chat_n, True)
    check(f'{label}: same positions in P', sorted(p_o) == sorted(p_n), True)
    worst = 0.0
    for position, column in p_o.items():
        for base, value in column.items():
            worst = max(worst, abs(value - p_n.get(position, {}).get(base, 0.0)))
    print(f'    max |dP| = {worst:.3e}   |d(delta)| = {abs(delta_o - delta_n):.3e}')
    check(f'{label}: column distributions agree', worst < tol, True)
    check(f'{label}: delta agrees', abs(delta_o - delta_n) < tol, True)


# ---------------------------------------------------------------- unit

def test_matches_on_a_simple_snapshot():
    snapshot = [[{}, 1.0, 0.0],
                [{100: alt(100)}, 0.6, 3.0],
                [{100: alt(100), 250: alt(250)}, 0.4, 5.0],
                [{250: alt(250, 1)}, 0.9, 11.0]]
    compare(*both(snapshot, 12.0), label='simple')


def test_matches_across_a_long_span():
    '''Rebasing has to hold over a production horizon, not just a few days.'''
    rng = np.random.default_rng(7)
    positions = [137, 900, 1801, 3000, 4500]
    snapshot = [[{}, 1.0, 0.0]]
    for day in range(1, 1096):
        genome = {p: alt(p, int(rng.integers(0, 3)))
                  for p in positions if rng.uniform() < 0.3}
        snapshot.append([genome, float(rng.uniform(0.01, 1.0)), float(day)])
    compare(*both(snapshot, 1095.0), label='1095 days')
    print(f'    snapshot entries: {len(snapshot)}')


def test_incremental_evaluation_matches_single_shot():
    '''Evaluating repeatedly (which rescales) must not drift from evaluating
    once at the end.'''
    rng = np.random.default_rng(11)
    positions = [500, 1900, 4100]
    entries = []
    for day in range(1, 1096):
        entries.append([{p: alt(p) for p in positions if rng.uniform() < 0.4},
                        float(rng.uniform(0.1, 1.0)), float(day)])

    stepped = c.ConsensusAccumulator()
    fed = 0
    for t_eval in range(5, 1096, 5):          # the model's rebuild cadence
        while fed < len(entries) and entries[fed][2] <= t_eval:
            stepped.add(*entries[fed])
            fed += 1
        stepped.consensus(float(t_eval))
    while fed < len(entries):
        stepped.add(*entries[fed])
        fed += 1
    stepped_result = stepped.consensus(1095.0)

    oneshot = c.ConsensusAccumulator()
    for e in entries:
        oneshot.add(*e)
    compare(oneshot.consensus(1095.0), stepped_result,
            label='219 rescalings vs one', tol=1e-9)


def test_genome_keying_is_free():
    '''Merging entries that share a genome, summing their frequencies, must
    give the identical answer -- weights enter every sum linearly.'''
    g1 = {300: alt(300)}
    g2 = {300: alt(300), 4700: alt(4700)}
    split = [[{}, 1.0, 0.0],
             [dict(g1), 0.3, 4.0], [dict(g1), 0.25, 4.0], [dict(g1), 0.1, 4.0],
             [dict(g2), 0.2, 4.0], [dict(g2), 0.15, 4.0]]
    keyed = [[{}, 1.0, 0.0],
             [dict(g1), 0.65, 4.0],
             [dict(g2), 0.35, 4.0]]
    old_split, new_split = both(split, 9.0)
    old_keyed, new_keyed = both(keyed, 9.0)
    compare(old_split, old_keyed, label='keying, old path', tol=1e-12)
    compare(new_split, new_keyed, label='keying, accumulator', tol=1e-12)


# ---------------------------------------------------------------- live

class TeeAccumulator(c.ConsensusAccumulator):
    '''The accumulator the model now uses, which additionally keeps the entry
    list the old path needs so both can be compared on identical data.

    The simulation is driven by the accumulator (what ships); the rescan is
    computed alongside purely to check it agrees.
    '''

    def __init__(self, stats):
        super().__init__()
        self.entries = []
        self.stats = stats

    def add(self, genome, frequency, t):
        super().add(genome, frequency, t)
        self.entries.append([genome, frequency, t])

    def consensus(self, t):
        t0 = time.perf_counter()
        new = super().consensus(t)
        t1 = time.perf_counter()
        old = c.get_consensus(self.entries, t)
        t2 = time.perf_counter()

        s = self.stats
        s['calls'] += 1
        s['new_s'] += t1 - t0
        s['old_s'] += t2 - t1
        s['entries'].append(len(self.entries))
        if old[0] != new[0]:
            s['seq_mismatch'] += 1
        for position, column in old[1].items():
            for base, value in column.items():
                s['worst_dP'] = max(
                    s['worst_dP'],
                    abs(value - new[1].get(position, {}).get(base, 0.0)))
        s['worst_ddelta'] = max(s['worst_ddelta'], abs(old[2] - new[2]))
        return new


def live(final_time, population_size, seeds):
    '''Run a real simulation with both paths computed at every rebuild.

    The OLD result is what gets returned, so the simulation is driven exactly
    as it is today and nothing about the trajectory changes.
    '''
    stats = {'calls': 0, 'old_s': 0.0, 'new_s': 0.0, 'entries': [],
             'seq_mismatch': 0, 'worst_dP': 0.0, 'worst_ddelta': 0.0}

    real_create = pop.create_population

    def create_population(parameters):
        population = real_create(parameters)
        population.consensus = TeeAccumulator(stats)
        return population

    name = 'zz_consensus_acc'
    fixed = {'population_size': population_size, 'final_time': final_time,
             'infected_individuals_at_start': 50, 'R': 1.03,
             'nucleotide_substitution_rate': 0.00014278943497365743,
             'long_shedders_ratio': 0.0, 'consensus': 'distribution'}

    pop.create_population = create_population
    try:
        wall0 = time.perf_counter()
        run_experiment(name, lambda: ({}, fixed, seeds),
                       simplicity_runner=serial, archive_experiment=False)
        wall = time.perf_counter() - wall0
    finally:
        pop.create_population = real_create
        try:
            shutil.rmtree(os.path.join(dm.get_data_dir(), name))
        except OSError:
            pass
    return stats, wall


def test_live(final_time, population_size, seeds):
    stats, wall = live(final_time, population_size, seeds)
    if not stats['calls']:
        check('live: rebuilds were intercepted', stats['calls'] > 0, True)
        return
    biggest = max(stats['entries'])
    print(f"    {stats['calls']} rebuilds, snapshot grew to {biggest:,} entries")
    print(f"    max |dP| = {stats['worst_dP']:.3e}   "
          f"|d(delta)| = {stats['worst_ddelta']:.3e}")
    check('live: consensus sequence identical at every rebuild',
          stats['seq_mismatch'], 0)
    check('live: column distributions agree', stats['worst_dP'] < 1e-9, True)
    check('live: delta agrees', stats['worst_ddelta'] < 1e-9, True)

    old_s, new_s = stats['old_s'], stats['new_s']
    speedup = old_s / new_s if new_s else float('inf')
    print(f"\n    rescan  : {old_s:8.3f} s total, "
          f"{1000 * old_s / stats['calls']:7.2f} ms per rebuild")
    print(f"    accumul.: {new_s:8.3f} s total, "
          f"{1000 * new_s / stats['calls']:7.2f} ms per rebuild")
    print(f"    speedup : {speedup:.1f}x   "
          f"(rescan was {100 * old_s / wall:.1f}% of the {wall:.0f}s run)")
    check('live: the accumulator is faster', new_s < old_s, True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--final-time', type=int, default=1095)
    p.add_argument('--population-size', type=int, default=1000)
    p.add_argument('--seeds', type=int, default=1)
    p.add_argument('--skip-live', action='store_true')
    args = p.parse_args()

    for t in (test_matches_on_a_simple_snapshot,
              test_matches_across_a_long_span,
              test_incremental_evaluation_matches_single_shot,
              test_genome_keying_is_free):
        print(f'\n{t.__name__}')
        t()

    if not args.skip_live:
        print(f'\ntest_live  (N={args.population_size}, '
              f'T={args.final_time}, {args.seeds} seed)')
        test_live(args.final_time, args.population_size, args.seeds)

    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
