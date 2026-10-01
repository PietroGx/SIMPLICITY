'''
Stage 3 submits every scenario at once instead of waiting for each, and the
pipeline packs the run for download at the end.

Before v2.4.49 the five scenarios ran strictly one after another: each
run_seeded_simulations blocks until its own 30 seeds finish, so only 30 tasks
were ever in flight against a cap of 200. dispatch_all submits them together and
splits the cap between them.

Everything Slurm-facing is stubbed; dispatch_all and compress_results are the
real ones.

    python tests/test_parallel_dispatch.py
'''
import os, sys, threading, time, types

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import impact_long_shedders_unbound_exp as exp
import run_impact_long_shedders_unbound_pipeline as pipe

SCENARIOS = ['control', 'SOT', 'HIV_low', 'HIV_high', 'edge_case']
failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def rows():
    return [{'scenario_name': n, 'R': 1.03, 'IH_virus_emergence_rate': 0.1,
             'nucleotide_substitution_rate': 1.4e-4, 'long_shedders_ratio': 0.0}
            for n in SCENARIOS]


def _instrumented(hold=0.25):
    '''Replace dispatch_scenario with one that records overlap.'''
    live = {'now': 0, 'peak': 0}
    order, lock = [], threading.Lock()

    def fake(row, exp_num, runner, n_seeds, consensus='argmax'):
        with lock:
            live['now'] += 1
            live['peak'] = max(live['peak'], live['now'])
            order.append(row['scenario_name'])
        time.sleep(hold)            # stands in for the blocking poll loop
        with lock:
            live['now'] -= 1
    return fake, live, order


def test_parallel_overlaps():
    saved = exp.dispatch_scenario
    fake, live, order = _instrumented()
    try:
        exp.dispatch_scenario = fake
        exp.dispatch_all(rows(), 1, 'slurm', 30, 'argmax', parallel=True)
    finally:
        exp.dispatch_scenario = saved
    check('every scenario dispatched', sorted(order), sorted(SCENARIOS))
    check('all five were in flight at once', live['peak'], len(SCENARIOS))


def test_sequential_does_not_overlap():
    saved = exp.dispatch_scenario
    fake, live, order = _instrumented(hold=0.02)
    try:
        exp.dispatch_scenario = fake
        exp.dispatch_all(rows(), 1, 'slurm', 30, 'argmax', parallel=False)
    finally:
        exp.dispatch_scenario = saved
    check('--sequential keeps one at a time', live['peak'], 1)
    check('order follows the table', order, SCENARIOS)


def test_non_slurm_runners_stay_sequential():
    '''serial/multiprocessing already manage their own parallelism.'''
    saved = exp.dispatch_scenario
    fake, live, _ = _instrumented(hold=0.02)
    try:
        exp.dispatch_scenario = fake
        exp.dispatch_all(rows(), 1, 'multiprocessing', 30, 'argmax', parallel=True)
    finally:
        exp.dispatch_scenario = saved
    check('multiprocessing runner is not threaded', live['peak'], 1)


def test_cap_is_split_not_multiplied():
    saved = exp.dispatch_scenario
    before = os.environ.get('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM')
    seen = []

    def fake(row, exp_num, runner, n_seeds, consensus='argmax'):
        seen.append(int(os.environ[
            'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM']))
    try:
        exp.dispatch_scenario = fake
        os.environ['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'] = '200'
        exp.dispatch_all(rows(), 1, 'slurm', 30, 'argmax', parallel=True)
    finally:
        exp.dispatch_scenario = saved
        if before is None:
            os.environ.pop('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM',
                           None)
        else:
            os.environ['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'] = before
    check('every scenario sees the same budget', len(set(seen)), 1)
    check('200 split five ways', seen[0], 40)
    check('total stays within the cap', seen[0] * len(SCENARIOS) <= 200, True)


def test_cap_is_restored_between_calls():
    '''Two dispatches in one process must not divide the budget twice.'''
    saved = exp.dispatch_scenario
    before = os.environ.get('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM')
    seen = []

    def fake(row, exp_num, runner, n_seeds, consensus='argmax'):
        seen.append(int(os.environ[
            'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM']))
    try:
        exp.dispatch_scenario = fake
        os.environ['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'] = '200'
        exp.dispatch_all(rows(), 1, 'slurm', 30, 'argmax', parallel=True)
        after_first = os.environ[
            'SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM']
        exp.dispatch_all(rows(), 1, 'slurm', 30, 'distribution', parallel=True)
    finally:
        exp.dispatch_scenario = saved
        if before is None:
            os.environ.pop('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM',
                           None)
        else:
            os.environ['SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM'] = before
    check('the cap is put back after dispatching', after_first, '200')
    check('the second dispatch gets the same budget, not half of it',
          set(seen), {40})


def test_failures_are_collected_and_raised():
    saved = exp.dispatch_scenario

    def fake(row, exp_num, runner, n_seeds, consensus='argmax'):
        if row['scenario_name'] in ('SOT', 'edge_case'):
            raise RuntimeError(f"boom in {row['scenario_name']}")
    try:
        exp.dispatch_scenario = fake
        raised = None
        try:
            exp.dispatch_all(rows(), 1, 'slurm', 30, 'argmax', parallel=True)
        except RuntimeError as e:
            raised = e
    finally:
        exp.dispatch_scenario = saved
    check('a failing scenario still raises', raised is not None, True)
    check('the message names the scenario', 'boom in' in str(raised), True)


def test_compress_results_invocation():
    calls = []

    class Result:
        returncode = 0
        stdout = 'Submitting packaging job to SLURM...\nSubmitted batch job 42\n'

    saved = pipe.subprocess.run

    class Sink:
        def write(self, _): pass
        def flush(self): pass

    try:
        pipe.subprocess.run = lambda cmd, **kw: (calls.append((cmd, kw)), Result())[1]
        out = pipe.compress_results(7, Sink(), keep_slurm_logs=True)
    finally:
        pipe.subprocess.run = saved

    cmd, kw = calls[0]
    check('packs the Data directory', cmd[1], 'Data')
    check('suffix targets this run only', cmd[2], '_#7')
    check('packs on a compute node', '--slurm' in cmd, True)
    check('calls package_simplicity_data.sh',
          os.path.basename(cmd[0]), 'package_simplicity_data.sh')
    check('slurm logs are kept by default',
          kw['env']['SIMPLICITY_KEEP_SLURM_LOGS'], '1')
    check('returns the sbatch output', 'Submitted batch job 42' in out, True)

    calls.clear()
    try:
        pipe.subprocess.run = lambda cmd, **kw: (calls.append((cmd, kw)), Result())[1]
        pipe.compress_results(7, Sink(), keep_slurm_logs=False)
    finally:
        pipe.subprocess.run = saved
    check('--drop-slurm-logs is honoured',
          calls[0][1]['env']['SIMPLICITY_KEEP_SLURM_LOGS'], '0')


def test_packaging_failure_does_not_fail_the_run():
    class Result:
        returncode = 1
        stdout = 'pixz not found'

    class Sink:
        def write(self, _): pass
        def flush(self): pass

    saved = pipe.subprocess.run
    try:
        pipe.subprocess.run = lambda cmd, **kw: Result()
        out = pipe.compress_results(7, Sink())
    finally:
        pipe.subprocess.run = saved
    check('a packing failure is reported, not raised', out, None)


if __name__ == '__main__':
    for t in [test_parallel_overlaps, test_sequential_does_not_overlap,
              test_non_slurm_runners_stay_sequential,
              test_cap_is_split_not_multiplied, test_cap_is_restored_between_calls,
              test_failures_are_collected_and_raised,
              test_compress_results_invocation,
              test_packaging_failure_does_not_fail_the_run]:
        print(f'\n{t.__name__}')
        t()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
