#!/usr/bin/env python3
"""Verification for the three optimisations: they must change nothing.

1. build_weighted_consensus_matrix   -- reference-seeded form vs the direct form
2. the distance cache                -- same answers, and refreshed on rebuild
3. hamming_iw                        -- set-free form vs the set-union form

Run from the repo root:
    python tests/test_speedups.py
"""
import random

import numpy as np

import simplicity.phenotype.consensus as c
import simplicity.phenotype.distance as dis
import simplicity.phenotype.update as upd

R = dis.reference
_fails = []


def check(label, ok):
    if not ok:
        _fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label}")


def rand_genome(n, rng):
    return {p: rng.choice([b for b in 'ACGT' if b != R[p]])
            for p in rng.sample(range(len(R)), n)}


def direct_matrix(data, positions, bases, letter_index, num_index):
    """The straightforward form: every genome votes in every column."""
    m = np.zeros((len(bases), len(positions)))
    for seq, n_inf, w in data:
        for num in positions:
            ch = seq.get(num, R[num])
            m[letter_index[ch], num_index[num]] += n_inf * w
    return m


def test_matrix():
    rng = random.Random(11)
    for trial in range(5):
        data = [[rand_genome(rng.randint(1, 15), rng), rng.random(), rng.random()]
                for _ in range(40)]
        got, bases, pos = c.build_weighted_consensus_matrix(data)
        li = {b: i for i, b in enumerate(bases)}
        ni = {p: i for i, p in enumerate(pos)}
        want = direct_matrix(data, pos, bases, li, ni)
        check(f"matrix trial {trial} identical to the direct form",
              bool(np.allclose(got, want, atol=1e-9)))
    check("empty data returns an empty matrix",
          c.build_weighted_consensus_matrix([])[0].shape[1] == 0)


def set_union_hamming(a, b):
    """The previous implementation, kept here as the reference."""
    d = 0
    for p in set(a) | set(b):
        rb = R[p]
        if a.get(p, rb) != b.get(p, rb):
            d += 1
    return d


def test_hamming():
    rng = random.Random(12)
    bad = 0
    for _ in range(500):
        a, b = rand_genome(rng.randint(0, 20), rng), rand_genome(rng.randint(0, 20), rng)
        if dis.hamming_iw(a, b) != set_union_hamming(a, b):
            bad += 1
    check(f"hamming_iw matches the set-union form on 500 pairs (mismatches {bad})", bad == 0)
    # the C7 guard must survive
    p = 100
    check("reference-valued entry still scores 0", dis.hamming_iw({p: R[p]}, {}) == 0)


class FakePop:
    """Minimal stand-in exposing what update_fitness touches."""
    def __init__(self, genomes, lineages):
        self._g = genomes
        self.individuals = {0: {'IH_lineages': list(lineages),
                                'IH_lineages_fitness_score': [0.0] * len(lineages),
                                'fitness_score': 0.0}}
        self.diagnosed, self.recovered, self.size = 30, 20, 100
        self.active_lineages_n = 7

    def get_lineage_genome(self, name):
        return self._g[name]


def test_cache():
    rng = random.Random(13)
    names = [f"L{i}" for i in range(6)]
    genomes = {n: rand_genome(rng.randint(1, 8), rng) for n in names}
    consensus_a = ({100: 'A'}, {100: {'A': .6, 'C': .4, 'T': 0., 'G': 0.}}, 0.2)
    # two mutated positions, so the argmax distance differs from consensus_a's
    consensus_b = ({250: 'G', 400: 'C'},
                   {250: {'G': .7, 'T': .3, 'A': 0., 'C': 0.},
                    400: {'C': .8, 'A': .2, 'T': 0., 'G': 0.}}, 0.1)

    for mode in ('argmax', 'distribution'):
        f = upd.update_fitness_factory('immune_waning', mode)
        pop = FakePop(genomes, names)
        f(pop, [0], consensus_a)
        cached = list(pop.individuals[0]['IH_lineages_fitness_score'])
        # uncached reference, same consensus
        want = [upd.immune_waning_fitness_score(pop, genomes[n], consensus_a,
                                                mode == 'distribution') for n in names]
        check(f"{mode}: cached scores equal the uncached computation",
              all(abs(x - y) < 1e-12 for x, y in zip(cached, want)))
        # calling again with the SAME consensus must not change anything
        f(pop, [0], consensus_a)
        check(f"{mode}: second call with the same consensus is stable",
              all(abs(x - y) < 1e-12
                  for x, y in zip(pop.individuals[0]['IH_lineages_fitness_score'], want)))
        # a NEW consensus must refresh the cache
        f(pop, [0], consensus_b)
        got_b = list(pop.individuals[0]['IH_lineages_fitness_score'])
        want_b = [upd.immune_waning_fitness_score(pop, genomes[n], consensus_b,
                                                  mode == 'distribution') for n in names]
        check(f"{mode}: cache refreshed when the consensus is rebuilt",
              all(abs(x - y) < 1e-12 for x, y in zip(got_b, want_b)))
        check(f"{mode}: the two consensuses actually give different scores",
              any(abs(x - y) > 1e-12 for x, y in zip(cached, got_b)))


if __name__ == '__main__':
    print("-- 1. consensus matrix --");   test_matrix()
    print("-- 3. hamming_iw --");         test_hamming()
    print("-- 2. distance cache --");     test_cache()
    print("\nSPEEDUPS PASSED" if not _fails else f"\nSPEEDUPS FAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
