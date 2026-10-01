'''
Removing the in-simulation sequencing draw must not perturb the dynamics.

The draw was `population.rng6.uniform(0,1) < seq_rate` inside
population_model.diagnosis, and rng6 was used nowhere else, so dropping it (and
rng6 with it, created last from seeds_generator) should leave every other
output byte-identical. This compares two runs of the same seeds across the
change.

sequencing_data*.csv are EXPECTED to differ: the run now writes the complete
record (every diagnosed host) instead of a `sequencing_rate` sample. The old
file must be a subset of the new one, with identical sequences for the hosts it
sampled -- that is the real check, and it is made here.

    python tests/test_sequencing_removal.py <before_dir> <after_dir>
'''
import ast, csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

# individuals_data.csv is handled separately: it gains the t_diagnosis column,
# so it cannot be byte-identical -- every OTHER column in it must be.
DYNAMICS_FILES = ['final_time.csv', 'lineage_frequency.csv',
                  'phylogenetic_data.csv', 'simulation_trajectory.csv']
NEW_COLUMNS = {'t_diagnosis'}

failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def seed_dirs(root):
    out = {}
    for dirpath, _, files in os.walk(root):
        if 'final_time.csv' in files:
            out[os.path.basename(dirpath)] = dirpath
    return out


def read_rows(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def test_dynamics_bit_identical(before, after):
    b, a = seed_dirs(before), seed_dirs(after)
    check('same seeds present', sorted(b), sorted(a))
    for seed in sorted(set(b) & set(a)):
        for name in DYNAMICS_FILES:
            pb, pa = os.path.join(b[seed], name), os.path.join(a[seed], name)
            if not (os.path.isfile(pb) and os.path.isfile(pa)):
                check(f'{seed}/{name} present on both sides', False, True)
                continue
            with open(pb, 'rb') as f1, open(pa, 'rb') as f2:
                same = f1.read() == f2.read()
            check(f'{seed}/{name} byte-identical', same, True)


def test_individuals_data_unchanged_except_new_column(before, after):
    '''Every pre-existing column must carry the same values; t_diagnosis is
    the only addition.'''
    b, a = seed_dirs(before), seed_dirs(after)
    for seed in sorted(set(b) & set(a)):
        old = read_rows(os.path.join(b[seed], 'individuals_data.csv'))
        new = read_rows(os.path.join(a[seed], 'individuals_data.csv'))
        old_cols, new_cols = set(old[0]), set(new[0])
        check(f'{seed}: only t_diagnosis is added', new_cols - old_cols,
              NEW_COLUMNS)
        check(f'{seed}: no column is removed', old_cols - new_cols, set())
        check(f'{seed}: same number of individuals', len(new), len(old))
        differing = set()
        for orow, nrow in zip(old, new):
            for col in old_cols:
                if orow[col] != nrow[col]:
                    differing.add(col)
        check(f'{seed}: every pre-existing column is unchanged', differing,
              set())


def test_t_diagnosis_persisted(after):
    a = seed_dirs(after)
    for seed in sorted(a):
        rows = read_rows(os.path.join(a[seed], 'individuals_data.csv'))
        check(f'{seed}: t_diagnosis column exists',
              't_diagnosis' in rows[0], True)
        stamped = [r for r in rows if r.get('t_diagnosis') not in (None, '')]
        diagnosed = [r for r in rows if r.get('state') == 'diagnosed']
        check(f'{seed}: every diagnosed individual has a t_diagnosis',
              len(stamped), len(diagnosed))
        check(f'{seed}: some individuals were diagnosed', len(stamped) > 0, True)


def test_complete_record_superset(before, after):
    '''The old sampled file must be contained in the new complete one, with the
    same sequence for every (individual, lineage) it recorded.'''
    b, a = seed_dirs(before), seed_dirs(after)
    for seed in sorted(set(b) & set(a)):
        pb = os.path.join(b[seed], 'sequencing_data.csv')
        pa = os.path.join(a[seed], 'sequencing_data.csv')
        if not os.path.isfile(pb):
            print(f'  [skip] {seed}: no sampled file (rate drew nobody)')
            continue
        old, new = read_rows(pb), read_rows(pa)
        check(f'{seed}: complete record is larger', len(new) >= len(old), True)

        new_by_key = {}
        for r in new:
            new_by_key[(r['individual_index'], r['lineage_name'])] = r

        missing, mismatched = [], []
        for r in old:
            key = (r['individual_index'], r['lineage_name'])
            match = new_by_key.get(key)
            if match is None:
                missing.append(key)
            elif (ast.literal_eval(match['sequence'])
                  != ast.literal_eval(r['sequence'])):
                mismatched.append(key)
        check(f'{seed}: every sampled host is in the complete record',
              len(missing), 0)
        check(f'{seed}: sequences agree for every sampled host',
              len(mismatched), 0)
        old_hosts = {r['individual_index'] for r in old}
        new_hosts = {r['individual_index'] for r in new}
        print(f'    {seed}: {len(old_hosts)} sampled hosts -> '
              f'{len(new_hosts)} in the complete record '
              f'({len(old)} -> {len(new)} rows)')


def test_reconstruction_matches(after):
    '''reconstruct_sequencing_data off the saved CSVs must rebuild the same
    complete record the run wrote, and subsample it per individual.'''
    import simplicity.output_manager as om
    import simplicity.sequencing as sq

    a = seed_dirs(after)
    for seed in sorted(a):
        ssod = a[seed]
        written = read_rows(os.path.join(ssod, 'sequencing_data.csv'))
        rebuilt = sq.reconstruct_sequencing_data(
            om.read_individuals_data(ssod), om.read_phylogenetic_data(ssod),
            sequencing_ratio=1.0)
        check(f'{seed}: reconstruction has the same row count',
              len(rebuilt), len(written))
        same_seq = all(
            ast.literal_eval(w['sequence']) == r['sequence']
            for w, r in zip(written, rebuilt))
        check(f'{seed}: reconstruction reproduces every sequence', same_seq,
              True)

        half = sq.reconstruct_sequencing_data(
            om.read_individuals_data(ssod), om.read_phylogenetic_data(ssod),
            sequencing_ratio=0.5, seed=1)
        again = sq.reconstruct_sequencing_data(
            om.read_individuals_data(ssod), om.read_phylogenetic_data(ssod),
            sequencing_ratio=0.5, seed=1)
        check(f'{seed}: a subsample is reproducible from its seed',
              [r['lineage_name'] for r in half],
              [r['lineage_name'] for r in again])
        check(f'{seed}: a 0.5 subsample is no larger than the whole',
              len(half) <= len(rebuilt), True)
        # the draw is per individual: a host appears whole or not at all
        whole_by_host = {}
        for r in rebuilt:
            whole_by_host.setdefault(r['individual_index'], 0)
            whole_by_host[r['individual_index']] += 1
        sub_by_host = {}
        for r in half:
            sub_by_host.setdefault(r['individual_index'], 0)
            sub_by_host[r['individual_index']] += 1
        partial = [i for i, n in sub_by_host.items() if n != whole_by_host[i]]
        check(f'{seed}: no host is partially sequenced', len(partial), 0)
        break  # one seed is enough for the reconstruction contract


if __name__ == '__main__':
    before, after = sys.argv[1], sys.argv[2]
    print('\ntest_dynamics_bit_identical')
    test_dynamics_bit_identical(before, after)
    print('\ntest_individuals_data_unchanged_except_new_column')
    test_individuals_data_unchanged_except_new_column(before, after)
    print('\ntest_t_diagnosis_persisted')
    test_t_diagnosis_persisted(after)
    print('\ntest_complete_record_superset')
    test_complete_record_superset(before, after)
    print('\ntest_reconstruction_matches')
    test_reconstruction_matches(after)
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
