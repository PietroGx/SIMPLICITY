#!/usr/bin/env python3
"""C8 verification.

Unit:  d_P(consensus) == 0, d_P == d_H when every column is pure, and a
       position no snapshot carries costs exactly 1.
Live:  both pipelines run, land in different output paths, and diverge.

Run from the repo root:
    python tests/test_consensus_distance.py
"""
import os
import shutil
import filecmp

import simplicity.phenotype.consensus as c
import simplicity.phenotype.distance as dis
import simplicity.dir_manager as dm
from simplicity.runme import run_experiment
import simplicity.runners.serial as serial

R = dis.reference
P1, P2 = 100, 250
A1 = next(b for b in 'ACGT' if b != R[P1])
A2 = next(b for b in 'ACGT' if b != R[P2])
_fails = []


def check(label, got, want, tol=1e-12):
    ok = abs(got - want) < tol if isinstance(want, float) else got == want
    if not ok:
        _fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:46s} got={got!r} want={want!r}")


def profile(day, t=21.0):
    m, bases, pos = c.build_weighted_consensus_matrix(c.get_seq_weights(day, t))
    P, delta = c.consensus_distribution(m, bases, pos)
    return c.weighted_consensus(m, pos), P, delta


def unit():
    print(f"reference[{P1}]={R[P1]}  reference[{P2}]={R[P2]}")
    # contested: 55/45 at P1, 5/95 at P2
    chat, P, delta = profile([[{P1: A1}, 0.55, 0.0], [{P2: A2}, 0.05, 0.0], [{}, 0.40, 0.0]])
    print(f"\n  contested columns  consensus={chat}  delta={delta:.4f}")
    check("d_P(consensus) == 0", dis.distributional(chat, P, delta), 0.0)
    d_wt = dis.distributional({}, P, delta)
    d_rare = dis.distributional({P2: A2}, P, delta)
    print(f"       d_P(wild type)={d_wt:.4f} vs d_H={dis.hamming_iw({}, chat)}")
    print(f"       d_P(rare base)={d_rare:.4f} vs d_H={dis.hamming_iw({P2: A2}, chat)}")
    check("a rare base costs more than the wild type", d_rare > d_wt, True)
    check("a position outside C costs exactly 1",
          dis.distributional({7777: 'A'}, P, delta) - d_wt, 1.0)

    # pure columns: the two distances must agree exactly
    chat, P, delta = profile([[{P1: A1, P2: A2}, 1.0, 0.0]])
    print(f"\n  pure columns  consensus={chat}  delta={delta:.4f}")
    for g in ({}, {P1: A1}, {P1: A1, P2: A2}, {P2: A2}):
        check(f"d_P == d_H for {g}", dis.distributional(g, P, delta),
              float(dis.hamming_iw(g, chat)))


BASE = {'population_size': 300, 'nucleotide_substitution_rate': 0.0005,
        'infected_individuals_at_start': 20, 'R': 1.2, 'R_long': 3,
        'diagnosis_rate_standard': 0.1, 'diagnosis_rate_long': 0.01,
        'final_time': 120}


def fixture(mode):
    def f():
        fixed = dict(BASE)
        if mode:
            fixed['consensus'] = mode
        return ({'long_shedders_ratio': [0.10]}, fixed, 1)
    return f


def seed_dir(name):
    sod = dm.get_simulation_output_dirs(name)[0]
    return dm.get_seeded_simulation_output_dirs(sod)[0]


def live():
    a, b = 'zz_cdist_argmax', 'zz_cdist_dist'
    run_experiment(a, fixture(None), simplicity_runner=serial, archive_experiment=False)
    run_experiment(b, fixture('distribution'), simplicity_runner=serial, archive_experiment=False)
    da, db = seed_dir(a), seed_dir(b)
    na, nb = os.path.basename(os.path.dirname(da)), os.path.basename(os.path.dirname(db))
    print(f"\n  argmax path: {na}\n  dist   path: {nb}")
    check("default pipeline's path carries no cdist segment", 'cdist' in na, False)
    check("distribution pipeline's path is distinguished", 'cdist_distribution' in nb, True)
    same = filecmp.cmp(os.path.join(da, 'phylogenetic_data.csv'),
                       os.path.join(db, 'phylogenetic_data.csv'), shallow=False)
    check("the two pipelines diverge", same, False)
    return [a, b]


if __name__ == '__main__':
    made = []
    unit()
    try:
        made = live()
    finally:
        for name in made:
            d = os.path.join(dm.get_data_dir(), name)
            if os.path.isdir(d):
                shutil.rmtree(d)
        if made:
            print(f"\ncleaned up {len(made)} scratch experiments")
    print("\nC8 PASSED" if not _fails else f"\nC8 FAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
