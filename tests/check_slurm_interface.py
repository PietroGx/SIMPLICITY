#!/usr/bin/env python3
'''Does this Slurm's output still parse the way the pipeline assumes?

Run this ON THE CLUSTER, after changing simplicity/runners/slurm.py or after a
Slurm upgrade.

    python tests/check_slurm_interface.py              # probe only, ~1 minute
    python tests/check_slurm_interface.py --raw        # and dump every reply
    python tests/check_slurm_interface.py --keep       # leave the array queued

WHY THIS EXISTS

tests/test_slurm_lifecycle.py already drives the whole control plane against
tests/fakeslurm/, and it reaches the cases a real controller will not produce on
demand -- an OOM kill, a requeued-held launch failure, a lost terminal state
write. What a fake cannot tell you is whether it lies about the INTERFACE.

That gap is not hypothetical here. The cluster runs slurm 26.05.4; the newest
Slurm packaged for Ubuntu 24.04 is 23.11.4, five releases behind. So a local
controller would confirm 23.11's output formats and say nothing about the ones
the pipeline actually meets. This check runs where the version is the one that
matters.

WHAT IT CHECKS

The three replies simplicity/runners/slurm.py parses, and the exact shape each
parser needs. This asserts SHAPE, not behaviour -- the logic is the fake's job:

  squeue --Format=ArrayJobID --name=X                     slurm.py:~600
      release_simulations drops line 0 as a header and requires exactly one
      distinct ArrayJobID across the rest. If --Format ever stops emitting that
      header, line 0 becomes a real job id and one task goes unreleased; if it
      emits two headers, the set has two members and the assert fires.

  squeue --Format=ArrayJobID,ArrayTaskID,Reason --noheader   slurm.py:~510
      reconcile_launch_failures does line.split(None, 2) and needs three
      whitespace-separated fields, the third being free text. A requeued-held
      task is only found by substring-matching that text.

  sacct -j ID --format=JobID,State --noheader --parsable2 -X   slurm.py:~450
      reconcile_terminated_tasks splits on the first "|" and expects one line
      per array task, JobID as "<job>_<task>". This is the branch that unblocked
      profile_grid_#910; if the shape changes it silently matches nothing and
      the polling loop waits forever.

It also confirms sacct RETAINS terminal state after the job leaves the queue.
With AccountingStorageType=accounting_storage/none it does not, and the
reconciler has nothing to read.

This submits a held array of three no-op tasks under its own job name. It is not
a pipeline and runs no simulation.
'''
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

JOB_NAME = 'simplicity_iface_probe'

findings = []
raw_replies = []


def check(label, ok, detail=''):
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label}" + (f": {detail}" if detail else ''))
    if not ok:
        findings.append(f'{label}{f": {detail}" if detail else ""}')
    return ok


def run(args, stdin=None, timeout=120):
    try:
        # input is always supplied, never None: a tool that reads stdin would
        # otherwise inherit the terminal's and block forever
        process = subprocess.run(args, input=stdin or b'',
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, timeout=timeout)
    except FileNotFoundError:
        return None, f'{args[0]} not on PATH'
    except subprocess.TimeoutExpired:
        return None, f'{args[0]} timed out after {timeout}s'
    reply = process.stdout.decode()
    raw_replies.append((' '.join(args), process.returncode, reply,
                        process.stderr.decode()))
    if process.returncode != 0:
        return None, f'exit {process.returncode}: {process.stderr.decode().strip()[:200]}'
    return reply, None


def version():
    reply, error = run(['sbatch', '--version'])
    if error:
        check('sbatch is available', False, error)
        return None
    text = reply.strip()
    print(f'  slurm   : {text}')
    return text


def submit():
    """A held array of three no-op tasks, under our own job name."""
    script = '#!/bin/sh\nsleep 2\n'
    reply, error = run(['sbatch', '--parsable', '--array=1-3', '--hold',
                        f'--job-name={JOB_NAME}', '--time=00:05:00',
                        '--mem=100M', '--output=/dev/null',
                        '--error=/dev/null'],
                       stdin=script.encode())
    if error:
        check('sbatch accepts a held array', False, error)
        return None
    job_id = reply.strip().split(';')[0]
    check('sbatch accepts a held array', bool(job_id), f'job {job_id}')
    return job_id


def check_release_listing(job_id):
    """squeue --Format=ArrayJobID --name=X -- the header assumption."""
    print('\nsqueue, release path  (slurm.py release_simulations)')
    reply, error = run(['squeue', '--Format=ArrayJobID', f'--name={JOB_NAME}'])
    if error:
        check('squeue --Format replies', False, error)
        return
    lines = reply.splitlines()
    check('squeue --Format replies', True, f'{len(lines)} line(s)')
    if not lines:
        check('there is a header line to drop', False, 'empty reply')
        return
    header = lines[0].strip()
    check('line 0 is a header, not a job id',
          not header.isdigit(), f'line 0 = {header!r}')
    # exactly what the real code does with the rest
    ids = {line.strip() for line in lines[1:]}
    check('exactly one distinct ArrayJobID in the rest', len(ids) == 1,
          f'{sorted(ids)}')
    check('and it is the job we submitted', ids == {job_id},
          f'{sorted(ids)} vs {job_id}')


def check_reason_listing():
    """squeue --Format=...,Reason --noheader -- three fields, free-text third."""
    print('\nsqueue, launch-failure path  (slurm.py reconcile_launch_failures)')
    reply, error = run(['squeue', '--name', JOB_NAME,
                        '--Format=ArrayJobID,ArrayTaskID,Reason',
                        '--noheader'])
    if error:
        check('squeue --noheader replies', False, error)
        return
    lines = [line for line in reply.splitlines() if line.strip()]
    check('squeue --noheader replies without a header', bool(lines),
          f'{len(lines)} line(s)')
    if not lines:
        return
    parsed = [line.split(None, 2) for line in lines]
    check('every line splits into three whitespace fields',
          all(len(parts) >= 3 for parts in parsed),
          f'first = {parsed[0]!r}')
    check('field 2 parses as an array task id',
          all(parts[1].isdigit() for parts in parsed if len(parts) >= 2),
          f'{[p[1] for p in parsed if len(p) >= 2]}')
    check('a held task reports a reason we could substring-match',
          any(parts[2].strip() for parts in parsed if len(parts) >= 3),
          f'reasons = {sorted({p[2].strip() for p in parsed if len(p) >= 3})}')


def check_sacct(job_id, released_task):
    """sacct --parsable2 -X -- the #910 branch's only source of truth."""
    print('\nsacct, reconciler path  (slurm.py reconcile_terminated_tasks)')
    reply, error = run(['sacct', '-j', job_id, '--format=JobID,State',
                        '--noheader', '--parsable2', '-X'])
    if error:
        check('sacct replies', False, error)
        return
    lines = [line for line in reply.splitlines() if line.strip()]
    if not check('sacct returns rows for the job', bool(lines),
                 f'{len(lines)} row(s)'):
        print('      ^ accounting is probably off. reconcile_terminated_tasks '
              'has nothing to read, which is the branch that unblocked #910.')
        return
    check('every row contains a "|" separator',
          all('|' in line for line in lines), f'first = {lines[0]!r}')
    rows = dict(line.split('|', 1) for line in lines if '|' in line)
    check('JobID is "<job>_<task>" for array tasks',
          any(key.startswith(f'{job_id}_') for key in rows),
          f'{sorted(rows)[:4]}')
    wanted = f'{job_id}_{released_task}'
    state = (rows.get(wanted) or '').strip().split()[:1]
    check('the task we released has a terminal state recorded',
          state and state[0] not in ('PENDING', 'RUNNING'),
          f'{wanted} -> {state}')
    # a trailing suffix is normal, e.g. "CANCELLED by 12345"
    import simplicity.runners.slurm as slurm
    known = (slurm.SLURM_TERMINAL_FAILURE_STATES
             | slurm.SLURM_TERMINAL_SUCCESS_STATES)
    check('that state is one the reconciler recognises',
          bool(state) and state[0] in known,
          f'{state} vs {sorted(known)}')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--raw', action='store_true',
                        help='print every command and its full reply')
    parser.add_argument('--keep', action='store_true',
                        help='leave the probe array queued instead of cancelling')
    args = parser.parse_args()

    print('slurm interface probe')
    if version() is None:
        sys.exit(1)

    job_id = submit()
    if job_id is None:
        sys.exit(1)

    try:
        time.sleep(3)           # let the controller register the array
        check_release_listing(job_id)
        check_reason_listing()

        print('\nreleasing one task, so sacct has a terminal state to keep')
        _, error = run(['scontrol', 'release', f'{job_id}_1'])
        check('scontrol release accepts "<job>_<task>"', error is None,
              error or 'accepted')
        for _ in range(20):
            time.sleep(3)
            reply, _ = run(['sacct', '-j', job_id, '--format=JobID,State',
                            '--noheader', '--parsable2', '-X'])
            if reply and any(s in reply for s in
                             ('COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT')):
                break
        check_sacct(job_id, released_task=1)
    finally:
        if not args.keep:
            run(['scancel', job_id])
            print(f'\ncancelled {job_id}')

    if args.raw:
        print('\n' + '=' * 72 + '\nraw replies\n' + '=' * 72)
        for command, code, out, err in raw_replies:
            print(f'\n$ {command}   (exit {code})')
            for line in (out or '').splitlines():
                print(f'  | {line}')
            for line in (err or '').splitlines():
                print(f'  ! {line}')

    print('\n' + '=' * 72)
    if findings:
        print(f'{len(findings)} assumption(s) do NOT hold on this Slurm:')
        for line in findings:
            print(f'  {line}')
        print('\nRe-run with --raw and compare against the parsers named in '
              'this file\'s docstring before trusting a submission.')
    else:
        print('every parsing assumption in simplicity/runners/slurm.py holds '
              'on this Slurm.')
    sys.exit(1 if findings else 0)


if __name__ == '__main__':
    main()
