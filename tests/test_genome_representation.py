#!/usr/bin/env python3
"""C1 + C7 verification.

Unit level: the distance functions and the consensus vote matrix on dict genomes.
Live run:   no duplicate positions, no reference-valued entries, clades still form.

Run from the repo root:
    python tests/test_genome_representation.py
"""
import os
import shutil

import simplicity.phenotype.distance as dis
import simplicity.phenotype.consensus as c
import simplicity.clustering as cl
import simplicity.dir_manager as dm
import simplicity.output_manager as om
from simplicity.runme import run_experiment
import simplicity.runners.serial as serial

EXPERIMENT = 'zz_genome_representation'
R = dis.reference
P1, P2 = 100, 250
ALT1 = next(b for b in 'ACGT' if b != R[P1])
ALT1b = next(b for b in 'ACGT' if b not in (R[P1], ALT1))   # a second allele at P1
ALT2 = next(b for b in 'ACGT' if b != R[P2])

_fails = []


def check(label, got, want):
    ok = got == want
    if not ok:
        _fails.append(label)
    print(f"  {'ok  ' if ok else 'FAIL'} {label:46s} got={got!r} want={want!r}")


def unit():
    print(f"reference[{P1}]={R[P1]}  reference[{P2}]={R[P2]}\n-- hamming --")
    check("empty genome", dis.hamming({}), 0)
    check("one mutation", dis.hamming({P1: ALT1}), 1)
    check("entry equal to the reference", dis.hamming({P1: R[P1]}), 0)

    print("-- hamming_iw (C7 guard) --")
    check("both empty", dis.hamming_iw({}, {}), 0)
    check("mutation vs empty", dis.hamming_iw({P1: ALT1}, {}), 1)
    check("reference-valued vs empty", dis.hamming_iw({P1: R[P1]}, {}), 0)
    check("identical", dis.hamming_iw({P1: ALT1}, {P1: ALT1}), 0)
    check("different base, same position", dis.hamming_iw({P1: ALT1}, {P1: ALT1b}), 1)

    print("-- parse_genome --")
    check("empty", cl.parse_genome({}), frozenset())
    check("one", cl.parse_genome({P1: ALT1}), frozenset({(P1, ALT1)}))
    try:
        cl.parse_genome(12345)
        check("malformed raises", False, True)
    except Exception:
        check("malformed raises instead of empty set", True, True)

    print("-- consensus matrix --")
    day = [[{P1: ALT1}, 0.6, 0.0], [{P1: ALT2, P2: ALT2}, 0.2, 0.0], [{}, 0.2, 0.0]]
    m, bases, pos = c.build_weighted_consensus_matrix(c.get_seq_weights(day, 21.0))
    check("every column totals 1.0", [round(x, 9) for x in m.sum(axis=0)], [1.0] * len(pos))
    check("consensus is a dict", isinstance(c.weighted_consensus(m, pos), dict), True)


def fixture():
    return ({'long_shedders_ratio': [0.05]},
            {'population_size': 300,
             'nucleotide_substitution_rate': 0.0005,
             'infected_individuals_at_start': 20,
             'R': 1.2, 'R_long': 3,
             'diagnosis_rate_standard': 0.1, 'diagnosis_rate_long': 0.01,
             'final_time': 120},
            1)


def live():
    n_lin = n_ref = n_bad = 0
    clades = []
    for sod in dm.get_simulation_output_dirs(EXPERIMENT):
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            phylo = om.read_phylogenetic_data(ssod)
            for g in phylo['Genome']:
                n_lin += 1
                if not isinstance(g, dict):
                    n_bad += 1
                    continue
                n_ref += sum(1 for p, b in g.items() if R[p] == b)
            clades.append(len(cl.cluster_lin_into_clades_with_meta(phylo, 5)[0]))
    print(f"\n  lineages {n_lin} | not a dict {n_bad} | reference-valued entries {n_ref}"
          f" | clades per simulation {clades}")
    check("every genome is a dict", n_bad, 0)
    check("no reference-valued entries", n_ref, 0)
    check("clades still form", bool(clades) and all(x > 1 for x in clades), True)


if __name__ == '__main__':
    unit()
    try:
        run_experiment(EXPERIMENT, fixture, simplicity_runner=serial,
                       archive_experiment=False)
        live()
    finally:
        d = os.path.join(dm.get_data_dir(), EXPERIMENT)
        if os.path.isdir(d):
            shutil.rmtree(d)
            print(f"\ncleaned up {d}")
    print("\nC1+C7 PASSED" if not _fails else f"\nC1+C7 FAILED: {_fails}")
    raise SystemExit(1 if _fails else 0)
