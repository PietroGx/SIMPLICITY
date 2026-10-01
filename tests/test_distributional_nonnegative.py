'''
The distributional consensus distance must never be negative.

Production run unbound #1 lost most of its distribution-consensus seeds to
"Probabilities are not non-negative" from population_model.py:176, where
rng4.choice is handed per-host fitness_score as p=. The guard there tests only
`fitness_sum > 0`, which a mix of one tiny negative and many positives passes.

The cause was round-off, not the model: d_P is >= 0 by construction -- the
consensus base is each column's argmax, so P(ref) - P(l) >= P(ref) - P(chat) at
every position and the sum cannot fall below -delta. A lineage that IS the
consensus reaches exactly 0 by near-total cancellation, which in doubles landed
at -5.551115123125783e-17 (one ulp near delta=1.70). Once phi = 1 the
(1-phi)/n_act term is exactly 0, so that sign became the whole fitness.

    python tests/test_distributional_nonnegative.py [n_seeds]

Unit checks run always; the live production-parameter run needs the frozen
table and is skipped without it.
'''
import os, sys, shutil, traceback

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import simplicity.phenotype.consensus as c
import simplicity.phenotype.distance as dis

REF = dis.reference
P1, P2, P3 = 100, 250, 400
ALT1 = next(b for b in 'ACGT' if b != REF[P1])
ALT2 = next(b for b in 'ACGT' if b != REF[P2])
ALT3 = next(b for b in 'ACGT' if b != REF[P3])

TABLE = os.path.join(REPO, 'Data',
                     'impact_long_shedders_unbound_setup_data_#1',
                     'nsr_calibration_table.csv')
RERUN_STD_NSR = 0.00014279
REQUIRED = ['final_time.csv', 'individuals_data.csv', 'lineage_frequency.csv',
            'phylogenetic_data.csv', 'simulation_trajectory.csv']

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def profile(snapshot, t=21.0):
    """(chat, column_dist, delta) from a weighted consensus snapshot, through
    the same path the model uses (get_consensus takes the current time, which
    the Bateman kernel weights each entry's age against)."""
    return c.get_consensus(snapshot, t)


def test_consensus_scores_zero():
    '''The case that round-off pushed below zero.'''
    chat, P, delta = profile([[{P1: ALT1, P2: ALT2}, 0.7, 0.0],
                              [{P1: ALT1}, 0.3, 0.0]])
    d = dis.distributional(chat, P, delta)
    print(f'    delta={delta!r}  d(consensus)={d!r}')
    check('consensus distance is not negative', d >= 0.0, True)
    check('consensus distance is zero', abs(d) < 1e-9, True)


def test_never_negative_over_many_profiles():
    '''Sweep weights and genomes; d_P must stay >= 0 throughout.'''
    worst = float('inf')
    n = 0
    for w in (0.05, 0.2, 0.5, 0.7, 0.9, 0.95, 0.99):
        snapshot = [[{P1: ALT1, P2: ALT2, P3: ALT3}, w, 0.0],
                    [{P1: ALT1, P2: ALT2}, 1.0 - w, 0.0]]
        chat, P, delta = profile(snapshot)
        for g in ({}, {P1: ALT1}, {P2: ALT2}, {P3: ALT3},
                  {P1: ALT1, P2: ALT2}, {P1: ALT1, P3: ALT3},
                  {P1: ALT1, P2: ALT2, P3: ALT3}, chat,
                  {P1: ALT1, P2: ALT2, P3: ALT3, 900: ALT1}):
            d = dis.distributional(g, P, delta)
            worst = min(worst, d)
            n += 1
    print(f'    {n} genome/profile pairs, minimum d_P = {worst!r}')
    check('no negative distance anywhere', worst >= 0.0, True)


def test_unmapped_position_still_costs_one():
    '''The clamp must not flatten real distances toward zero.'''
    chat, P, delta = profile([[{P1: ALT1}, 1.0, 0.0]])
    far = dis.distributional({P1: ALT1, 900: ALT2}, P, delta)
    near = dis.distributional({P1: ALT1}, P, delta)
    print(f'    d(consensus)={near!r}  d(consensus + unmapped)={far!r}')
    check('an unmapped position still costs 1', abs(far - near - 1.0) < 1e-9,
          True)
    check('distances are still ordered', far > near, True)


def test_live_production_seeds(n_seeds):
    '''Before the clamp this was 1/6 at these parameters; argmax was 6/6.'''
    if not os.path.isfile(TABLE):
        print('  [skip] frozen calibration table not present')
        return
    import pandas as pd
    import simplicity.dir_manager as dm
    from simplicity.runme import run_experiment
    import simplicity.runners.serial as serial
    import impact_long_shedders_unbound_config as ucfg

    row = pd.read_csv(TABLE)
    row = row[row['scenario_name'] == 'control'].iloc[0].to_dict()
    row['nucleotide_substitution_rate'] = RERUN_STD_NSR

    name = 'zz_distnonneg'
    try:
        run_experiment(name,
                       ucfg.build_exp_scenario_settings(row, n_seeds,
                                                        consensus='distribution'),
                       simplicity_runner=serial, archive_experiment=False)
    except Exception:
        print('  run_experiment raised:')
        traceback.print_exc()
    total = ok = 0
    for sod in dm.get_simulation_output_dirs(name):
        for ssod in dm.get_seeded_simulation_output_dirs(sod):
            total += 1
            ok += all(os.path.isfile(os.path.join(ssod, f)) for f in REQUIRED)
    try:
        shutil.rmtree(os.path.join(dm.get_data_dir(), name))
    except OSError:
        pass
    check(f'all {total} production seeds completed', ok, total)


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    for t in (test_consensus_scores_zero, test_never_negative_over_many_profiles,
              test_unmapped_position_still_costs_one):
        print(f'\n{t.__name__}')
        t()
    print('\ntest_live_production_seeds')
    test_live_production_seeds(n)
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
