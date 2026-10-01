'''
Exercises Host._load_or_precompute_exponentials, which nine tasks of the
production array died in (EOFError: Ran out of input at pickle.load).
Covers: cold write, warm read, recovery from a zero-byte file left by a
concurrent writer, and 16 processes racing on one shared path.
'''
import os, sys, time, pickle, tempfile, multiprocessing as mp
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import simplicity.dir_manager as dm

TAUS = dict(tau_1=2.86, tau_2=3.91, tau_3=7.5, tau_4=8)


def _build(data_dir):
    dm.set_data_dir(data_dir)
    import simplicity.intra_host_model as h
    return h.Host(update_mode='matrix', **TAUS)


def _path(data_dir):
    dm.set_data_dir(data_dir)
    import simplicity.output_manager as om
    return om.get_procomputed_matrix_table_filepath(
        TAUS['tau_1'], TAUS['tau_2'], TAUS['tau_3'], TAUS['tau_4'])


def _tables_equal(a, b):
    if set(a) != set(b):
        return False
    return all(np.array_equal(a[k], b[k]) for k in a)


def test_cold_write_leaves_no_temp():
    with tempfile.TemporaryDirectory() as d:
        host = _build(d)
        p = _path(d)
        assert os.path.exists(p), 'table not written'
        assert len(host.exp_table) == 300, len(host.exp_table)
        strays = [f for f in os.listdir(d) if f.endswith('.tmp')]
        assert not strays, f'temp files left behind: {strays}'
        print('cold write ok')


def test_warm_read_matches():
    with tempfile.TemporaryDirectory() as d:
        first = _build(d).exp_table
        second = _build(d).exp_table
        assert _tables_equal(first, second), 'reloaded table differs'
        print('warm read ok')


def test_recovers_from_zero_byte_file():
    with tempfile.TemporaryDirectory() as d:
        good = _build(d).exp_table
        p = _path(d)
        open(p, 'wb').close()               # what a concurrent writer leaves
        assert os.path.getsize(p) == 0
        recovered = _build(d).exp_table     # pre-fix: EOFError
        assert _tables_equal(good, recovered)
        assert os.path.getsize(p) > 0, 'zero-byte file not replaced'
        print('zero-byte recovery ok')


def test_recovers_from_truncated_file():
    with tempfile.TemporaryDirectory() as d:
        good = _build(d).exp_table
        p = _path(d)
        with open(p, 'rb') as f:
            head = f.read(4096)
        with open(p, 'wb') as f:
            f.write(head)                   # mid-stream truncation
        recovered = _build(d).exp_table
        assert _tables_equal(good, recovered)
        print('truncated recovery ok')


def _worker(args):
    # fork, not spawn: a spawned child re-imports simplicity and re-fetches the
    # reference genome from NCBI, which 429s. All that matters here is that the
    # processes are distinct and hit the shared path at the same instant.
    data_dir, start_at = args
    while time.time() < start_at:
        pass
    try:
        t = _build(data_dir).exp_table
        return len(t), sorted(t)[:3]
    except Exception as e:
        return f'{type(e).__name__}: {e}'


def test_concurrent_builders():
    n = 16
    with tempfile.TemporaryDirectory() as d:
        start_at = time.time() + 2.0
        with mp.get_context('fork').Pool(n) as pool:
            results = pool.map(_worker, [(d, start_at) for _ in range(n)])
        bad = [r for r in results if isinstance(r, str)]
        assert not bad, f'{len(bad)}/{n} workers raised: {bad[:3]}'
        assert all(r[0] == 300 for r in results), [r[0] for r in results]
        assert len({tuple(r[1]) for r in results}) == 1, 'tables disagree'
        strays = [f for f in os.listdir(d) if f.endswith('.tmp')]
        assert not strays, f'temp files left behind: {strays}'
        print(f'{n} concurrent builders ok')


if __name__ == '__main__':
    test_cold_write_leaves_no_temp()
    test_warm_read_matches()
    test_recovers_from_zero_byte_file()
    test_recovers_from_truncated_file()
    test_concurrent_builders()
    print('\nall passed')
