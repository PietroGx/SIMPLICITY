'''
Exercises run_seeded_simulations' reporting: the SimulationsStatus line must
print once per actual status change (never on a timer) and carry a timestamp,
and the long-running report's threshold must be reachable by real tasks.

Everything Slurm-facing is stubbed; the loop itself is the real one.
'''
import io, os, re, sys, time, contextlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import simplicity.runners.slurm as slurm

TS = r'\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\](?: \S+)? '
failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def _status(**kw):
    base = dict(total=10, submitted=10, released=0, left=10, pending=0,
                started=0, running=0, completed=0, failed=0)
    base.update(kw)
    return slurm.SimulationsStatus(**base)


def _drive(sequence):
    '''Run the real run_seeded_simulations over a scripted status sequence.'''
    seq = list(sequence)
    saved = {n: getattr(slurm, n) for n in
             ('submit_simulations', 'release_simulations',
              'reconcile_terminated_tasks', 'reconcile_launch_failures',
              'report_long_running_simulations', 'poll_simulations_status')}
    real_sleep = time.sleep
    try:
        slurm.submit_simulations = lambda *a, **k: None
        slurm.release_simulations = lambda *a, **k: None
        slurm.reconcile_terminated_tasks = lambda *a, **k: None
        slurm.reconcile_launch_failures = lambda *a, **k: None
        slurm.report_long_running_simulations = lambda *a, **k: None
        slurm.poll_simulations_status = lambda exp: seq.pop(0)
        time.sleep = lambda s: None
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            slurm.run_seeded_simulations('probe', lambda: None)
        return buf.getvalue()
    finally:
        time.sleep = real_sleep
        for n, fn in saved.items():
            setattr(slurm, n, fn)


def test_thresholds_are_reachable():
    '''unbound #1's slowest of 1,400 tasks took 2,423 s; HIV_high's p95 was
    1,217 s. Post-speedup those are ~370 s and ~184 s.'''
    check('threshold below the observed slowest task',
          slurm.LONG_RUNNING_THRESHOLD_S < 2423, True)
    check('threshold below the post-speedup slowest task (~370s)',
          slurm.LONG_RUNNING_THRESHOLD_S < 370, True)
    check('report fires at least as often as the threshold',
          slurm.LONG_RUNNING_REPORT_INTERVAL_S <= 2 * slurm.LONG_RUNNING_THRESHOLD_S,
          True)


def test_prints_once_per_change():
    first = _status()
    second = _status(released=10, started=3, pending=7, running=3)
    done = _status(released=10, started=10, left=0, completed=10)
    # polled once before submit, then once per while-condition evaluation
    out = _drive([first, first, first, first, second, second, done])
    check('every status line names its experiment',
          all('probe' in l for l in out.splitlines()
              if 'SimulationsStatus' in l), True)
    lines = [l for l in out.splitlines() if 'SimulationsStatus' in l]
    for l in lines:
        print(f'    {l}')
    check('one line per distinct consecutive status', len(lines), 3)
    check('no two consecutive lines identical after the timestamp',
          len({re.sub(TS, '', l) for l in lines}), 3)
    check('every status line is timestamped',
          all(re.match(TS, l) for l in lines), True)
    check('final line reports completion',
          'completed=10' in lines[-1] and 'left=0' in lines[-1], True)


def test_no_timer_reprint():
    '''The old loop reprinted an unchanged status every 17 s. Hold one status
    across many iterations and it must still print exactly once.'''
    stable = _status()
    done = _status(left=0, completed=10, started=10, released=10)
    out = _drive([stable] * 40 + [done])
    lines = [l for l in out.splitlines() if 'SimulationsStatus' in l]
    check('40 unchanged polls produce one line, plus the final', len(lines), 2)


def test_timestamp_format():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        slurm.print_simulations_status(_status())
    line = buf.getvalue().strip()
    print(f'    {line}')
    check('timestamped', bool(re.match(TS, line)), True)
    check('status payload intact', 'SimulationsStatus(total=10' in line, True)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        slurm.print_simulations_status(_status(), 'my_experiment_#3')
    tagged = buf.getvalue().strip()
    print(f'    {tagged}')
    check('experiment name is on the line when given',
          'my_experiment_#3' in tagged, True)
    check('tagged line still matches the format',
          bool(re.match(TS, tagged)), True)
    stamp = line[1:20]
    drift = abs(time.mktime(time.strptime(stamp, '%Y-%m-%d %H:%M:%S')) - time.time())
    check('timestamp is now (<5s drift)', drift < 5, True)


def test_queued_message_is_unambiguous():
    """'submitted N' read as 'N are running'. The line now says what actually
    happens: the whole array goes to Slurm at once, held, and is released a
    capped number at a time."""
    import os
    os.environ.setdefault('SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM', '200')
    stable = _status(total=900, submitted=900, left=900)
    done = _status(total=900, submitted=900, released=900, started=900,
                   completed=900, left=0)
    out = _drive([stable, done])
    line = next(l for l in out.splitlines() if 'queued' in l)
    print(f'    {line}')
    check('names the experiment', 'probe' in line, True)
    check('says they are queued, not running', 'queued 900' in line, True)
    check('says the array is held', 'held Slurm array' in line, True)
    check('states the release cap', 'at a time' in line, True)
    check('no longer says "submitted"', 'submitted' in line, False)


if __name__ == '__main__':
    for t in [test_thresholds_are_reachable, test_prints_once_per_change,
              test_no_timer_reprint, test_timestamp_format,
              test_queued_message_is_unambiguous]:
        print(f'\n{t.__name__}')
        t()
    print('\n' + ('FAILURES:\n  ' + '\n  '.join(failures) if failures
                  else 'all passed'))
    sys.exit(1 if failures else 0)
