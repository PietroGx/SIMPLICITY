#!/usr/bin/env python3
'''Turn the per-process samples from tests/profiling/sitecustomize.py into a
report: what the memory is spent on, what grows, and where the time goes.

    export PYTHONPATH=$PWD/tests/profiling:$PYTHONPATH
    export SIMPLICITY_PROFILE_DIR=$PWD/profiles
    <launch the pipeline as usual>
    python tests/profile_report.py profiles

The profiler samples while the run is alive and fsyncs every row, so a task
Slurm kills for exceeding --mem still leaves its growth curve behind: that
partial curve is the point. Nothing in simplicity/ or scripts/ is involved --
Python imports sitecustomize at startup from PYTHONPATH, and
simplicity/runners/slurm.py passes the whole environment through sbatch.

How to read the three memory numbers, which answer different questions:

  measured MB   sampled directly off the container (section 2). This is the
                one to act on. It OVER-counts anything shared between two
                containers, so the column sums to more than RSS; compare each
                row against RSS, not the total.

  fixed/grows   whether the container's length rose over the run. A fixed
                container is a constant cost paid by every task however short;
                a growing one is what decides whether a long run survives.
                Different fixes, so never pool them.

  bytes/row     slope of RSS against the container's length, fitted within
                each process and then taken as a median (section 3). High R^2
                means that container's growth tracks the process's memory -- it
                does not prove cause, since the containers grow together. Use
                it as a cross-check on the measured figure, not as a second
                opinion with equal weight.
'''
import csv
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

import numpy as np

CONTAINERS = ('individuals', 'lineage_frequency', 'trajectory',
              'phylogenetic_data', 'phylodots', 'consensus_sequences_t',
              'fitness_trajectory', 'R_effective_trajectory',
              'infected_i', 'infectious_i', 'recovered_i', 'reservoir_i')

W = 96


def rule(char='-'):
    print(char * W)


def heading(title):
    rule()
    print(title)
    rule()


def fmt(value, spec='.1f', dash='-'):
    return dash if value is None else format(value, spec)


def read_samples(path):
    rows = []
    with open(path) as handle:
        for row in csv.DictReader(handle):
            parsed = {}
            for key, value in row.items():
                if value == '' or value is None:
                    parsed[key] = None
                else:
                    try:
                        parsed[key] = float(value)
                    except ValueError:
                        parsed[key] = None
            rows.append(parsed)
    return rows


def read_allocations(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) == 3:
                try:
                    rows.append((float(parts[0]), int(parts[1]), parts[2]))
                except ValueError:
                    pass
    return rows


def load(profile_dir):
    """One entry per profiled process: samples, meta, folded stacks, allocs."""
    tasks = []
    for experiment in sorted(os.listdir(profile_dir)):
        exp_dir = os.path.join(profile_dir, experiment)
        if not os.path.isdir(exp_dir):
            continue
        for name in sorted(os.listdir(exp_dir)):
            if not name.endswith('.csv'):
                continue
            tag = name[:-4]
            samples = read_samples(os.path.join(exp_dir, name))
            if not samples:
                continue
            meta = {}
            meta_path = os.path.join(exp_dir, tag + '.meta.json')
            if os.path.isfile(meta_path):
                try:
                    with open(meta_path) as handle:
                        meta = json.load(handle)
                except (OSError, ValueError):
                    pass
            stacks = {}
            folded = os.path.join(exp_dir, tag + '.folded')
            if os.path.isfile(folded):
                with open(folded) as handle:
                    for line in handle:
                        stack, _, count = line.rpartition(' ')
                        if stack:
                            try:
                                stacks[stack] = int(count)
                            except ValueError:
                                pass
            allocs = []
            alloc_path = os.path.join(exp_dir, tag + '.alloc')
            if os.path.isfile(alloc_path):
                allocs = read_allocations(alloc_path)
            tasks.append({'experiment': experiment, 'tag': tag,
                          'samples': samples, 'meta': meta, 'stacks': stacks,
                          'allocs': allocs, 'scenario': None, 'seed': None,
                          'final_time': None})
    return tasks


def join_scenarios(tasks, data_dir):
    """Recover scenario and seed from the experiment's slurm id map, which
    slurm.job() writes under the same job/task ids the profile is named by."""
    for task in tasks:
        mapping = os.path.join(data_dir, task['experiment'], 'slurm',
                               'job_id_mapping',
                               f"{task['experiment']}_{task['tag']}.csv")
        if not os.path.isfile(mapping):
            continue
        try:
            with open(mapping) as handle:
                content = handle.read().strip()
        except OSError:
            continue
        resolved = resolve_repeat(data_dir, task['experiment'], content)
        if resolved is None:
            continue
        group, record = resolved
        task['group'] = group
        task['scenario'] = record['stem']
        task['seed'] = f"seed_{record['seed']:04d}"
        task['repeat'] = f"{group}/{record['stem']}/{task['seed']}"
        params_path = os.path.join(data_dir, task['experiment'],
                                   '02_Simulations', f"{record['stem']}.json")
        if os.path.isfile(params_path):
            try:
                with open(params_path) as handle:
                    task['final_time'] = json.load(handle).get('final_time')
            except (OSError, ValueError):
                pass


def resolve_repeat(data_dir, experiment, content):
    """(group, repeat record) from a slurm id map file's "<group>/<index>".

    The map used to hold an ABSOLUTE seeded-params path written on the cluster
    (/scratch/.../Data/...), which is why relocate() below had to exist: keying
    on it only matched when the report ran on the same machine under the same
    root, so downloading a run and reporting on it locally silently lost every
    cell. A group and an index travel.
    """
    group, _, index = content.rpartition('/')
    if not group:
        return None
    path = os.path.join(data_dir, experiment, '03_Repeats', group,
                        'repeats.json')
    try:
        with open(path) as handle:
            records = json.load(handle)
        return group, records[int(index)]
    except (OSError, ValueError, IndexError, KeyError):
        return None


def relocate(path, data_dir):
    """Re-root an absolute path recorded on the cluster onto this machine.

    The slurm id map stores the path as the task saw it (/scratch/.../Data/...),
    so a downloaded run would otherwise find none of its parameter files.
    """
    if os.path.isfile(path):
        return path
    marker = os.sep + 'Data' + os.sep
    if marker in path:
        return os.path.join(data_dir, path.split(marker, 1)[1])
    return path


def live(task):
    """Samples taken while a Population existed -- the others are startup."""
    return [s for s in task['samples'] if s.get('sim_time') is not None]


def baseline_mb(task):
    rows = live(task)
    return rows[0]['rss_mb'] if rows else None


def peak_mb(task):
    peaks = [s['hwm_mb'] for s in task['samples'] if s.get('hwm_mb')]
    if peaks:
        return max(peaks)
    meta_peak = task['meta'].get('peak_rss_mb')
    if isinstance(meta_peak, (int, float)):
        return meta_peak
    seen = [s['rss_mb'] for s in task['samples'] if s.get('rss_mb')]
    return max(seen) if seen else None


def growth_per_day(task):
    """MB of RSS per simulated day, least squares over the live samples."""
    rows = live(task)
    if len(rows) < 3:
        return None, None
    days = np.array([r['sim_time'] for r in rows], dtype=float)
    rss = np.array([r['rss_mb'] for r in rows], dtype=float)
    span = float(days.max() - days.min())
    if span <= 0:
        return None, None
    return float(np.polyfit(days, rss, 1)[0]), span


def mem_limit_mb(tasks):
    for task in tasks:
        raw = task['meta'].get('slurm_mem')
        if not raw:
            continue
        text = str(raw).strip().upper()
        try:
            if text.endswith('G'):
                return float(text[:-1]) * 1024
            if text.endswith('M'):
                return float(text[:-1])
            return float(text)
        except ValueError:
            return None
    return None


def section_completeness(tasks, limit):
    heading('0. WHAT WAS PROFILED')
    experiments = defaultdict(int)
    for task in tasks:
        experiments[task['experiment']] += 1
    for experiment, count in sorted(experiments.items()):
        print(f'  {experiment}: {count} processes')
    clean = sum(1 for t in tasks if t['meta'].get('exited_cleanly'))
    print(f'  exited cleanly          : {clean} of {len(tasks)}')
    print(f'  killed or still running : {len(tasks) - clean}')
    print('  --mem per task          : '
          + (f'{fmt(limit)} MB' if limit else 'unknown'))
    intervals = sorted({t['meta'].get('sample_interval_s') for t in tasks} - {None})
    print(f'  sample interval(s)      : {intervals}')
    traced = sum(1 for t in tasks if t['allocs'])
    if traced:
        print(f'  tracemalloc captures    : {traced}')
    print()
    print('  A process that did not exit cleanly was killed (walltime, --mem,')
    print('  scancel) or was still alive when the report ran. Its last sample')
    print('  is the last thing it managed to write, which is what makes those')
    print('  rows the interesting ones.')


def section_memory(tasks, limit):
    heading('1. MEMORY PER PROCESS       the baseline is fixed, the growth is not')
    print(f"{'task':<18}{'scenario':<24}{'peak MB':>9}{'base MB':>9}"
          f"{'MB/day':>9}{'days':>7}{'proj MB':>9}  exit")
    rows = []
    for task in tasks:
        slope, span = growth_per_day(task)
        base = baseline_mb(task)
        projected = None
        if slope is not None and base is not None and task['final_time']:
            projected = base + slope * float(task['final_time'])
        rows.append((peak_mb(task) or 0, task, base, slope, span, projected))
    for peak, task, base, slope, span, projected in sorted(rows, key=lambda r: -r[0]):
        flag = 'ok' if task['meta'].get('exited_cleanly') else 'KILLED/LIVE'
        print(f"{task['tag'][:16]:<18}{(task['scenario'] or '')[:22]:<24}"
              f"{fmt(peak):>9}{fmt(base):>9}{fmt(slope, '.3f'):>9}"
              f"{fmt(span, '.0f'):>7}{fmt(projected, '.0f'):>9}  {flag}")
    if limit:
        over = [r for r in rows if r[5] and r[5] > limit]
        near = [r for r in rows if r[0] > 0.8 * limit]
        print()
        print(f'  peaked above 80% of --mem  : {len(near)} of {len(rows)}')
        print(f'  projected to exceed --mem  : {len(over)} of {len(rows)}')
        if over:
            print('  *** a projection extrapolates each task\'s own slope to its')
            print('  *** final_time: trust it only where "days" covers a decent')
            print('  *** share of the run.')


def section_breakdown(tasks):
    heading('2. WHAT THE MEMORY IS IN       measured off the containers themselves')
    worst = {}
    for task in tasks:
        rows = live(task)
        if not rows:
            continue
        first, last = rows[0], rows[-1]
        for name in CONTAINERS:
            megabytes = last.get('mb_' + name)
            length = last.get('n_' + name)
            if megabytes is None or length is None:
                continue
            start = first.get('n_' + name)
            grows = start is not None and length > start * 1.5 and length > start + 10
            current = worst.get(name)
            if current is None or megabytes > current[0]:
                worst[name] = (megabytes, length, grows, last.get('rss_mb'))
    extra = [(t, s) for t in tasks for s in live(t)[-1:]
             if s.get('mb_consensus') is not None]
    if extra:
        best = max(extra, key=lambda ts: ts[1]['mb_consensus'])[1]
        worst['consensus (accumulator)'] = (best['mb_consensus'],
                                            best.get('n_consensus_positions'),
                                            False, best.get('rss_mb'))
    if not worst:
        print('  no size samples (SIMPLICITY_PROFILE_SIZES=0?)')
        return
    print(f"{'container':<28}{'entries':>12}{'MB':>10}{'bytes/entry':>13}"
          f"{'% of RSS':>10}  kind")
    for name, (megabytes, length, grows, rss) in sorted(worst.items(),
                                                        key=lambda kv: -kv[1][0]):
        per = megabytes * 1048576.0 / length if length else None
        share = 100.0 * megabytes / rss if rss else None
        kind = 'GROWS' if grows else 'fixed'
        print(f'{name:<28}{fmt(length, ".0f"):>12}{megabytes:>10.1f}'
              f'{fmt(per, ".0f"):>13}{fmt(share, ".0f"):>10}  {kind}')
    print()
    print('  Worst case over the profiled processes, at each one\'s last live')
    print('  sample. Shared objects are charged to every container that holds')
    print('  them, so these sum to more than RSS -- read each row against RSS.')
    print()
    print('  GROWS = its length rose over the run, so it scales with run length')
    print('  and is what decides whether a long run survives. fixed = a')
    print('  constant every task pays however short it is. The two need')
    print('  different fixes and should never be pooled.')


def section_fit(tasks):
    heading('3. GROWTH PER ROW       RSS against each length, fitted per process')
    print(f"{'container':<28}{'runs':>6}{'max len':>12}{'bytes/row':>12}"
          f"{'R^2':>8}")
    results = []
    for name in CONTAINERS:
        slopes, scores, longest = [], [], 0.0
        for task in tasks:
            rows = live(task)
            pairs = [(s['n_' + name], s['rss_mb']) for s in rows
                     if s.get('n_' + name) is not None and s.get('rss_mb') is not None]
            if len(pairs) < 4:
                continue
            x = np.array([p[0] for p in pairs], dtype=float)
            y = np.array([p[1] for p in pairs], dtype=float)
            longest = max(longest, float(x.max()))
            if x.max() - x.min() <= 0:
                continue
            slope, intercept = np.polyfit(x, y, 1)
            residual = float(((y - (slope * x + intercept)) ** 2).sum())
            total = float(((y - y.mean()) ** 2).sum())
            if total <= 0:
                continue
            slopes.append(float(slope) * 1048576.0)
            scores.append(1 - residual / total)
        if not slopes:
            if longest:
                results.append((name, 0, longest, None, None))
            continue
        results.append((name, len(slopes), longest,
                        float(np.median(slopes)), float(np.median(scores))))
    for name, runs, longest, bytes_row, r2 in sorted(results,
                                                     key=lambda r: -(r[4] or 0)):
        print(f'{name:<28}{runs:>6}{longest:>12.0f}'
              f'{fmt(bytes_row, ".0f"):>12}{fmt(r2, ".3f"):>8}')
    print()
    print('  Median across processes of a fit done WITHIN each process. Pooling')
    print('  the samples instead would fit the gap between processes -- which')
    print('  differ in baseline -- rather than growth, and report near-zero R^2')
    print('  for everything.')
    print()
    print('  A cross-check on section 2, not a second opinion: the containers')
    print('  grow together, so this ranks correlation. Where it disagrees with')
    print('  the measured bytes/entry, the measured figure is the one to trust.')
    print('  A constant-length container cannot be fitted and shows "-".')


def section_time(tasks):
    heading('4. WHERE THE TIME GOES       sampled stacks, pooled across processes')
    leaves = defaultdict(int)
    frames = defaultdict(int)
    total = 0
    for task in tasks:
        for stack, count in task['stacks'].items():
            parts = stack.split(';')
            leaves[parts[-1]] += count
            for part in set(parts):
                frames[part] += count
            total += count
    if not total:
        print('  no stack samples (SIMPLICITY_PROFILE_STACKS=0?)')
        return
    print(f'  {total} samples pooled\n')
    print(f"{'leaf function (where the interpreter actually was)':<62}"
          f"{'samples':>9}{'%':>8}")
    for name, count in sorted(leaves.items(), key=lambda kv: -kv[1])[:20]:
        print(f'{name[:60]:<62}{count:>9}{100.0 * count / total:>7.1f}%')
    print()
    print(f"{'function on stack (inclusive: it, or something it called)':<62}"
          f"{'samples':>9}{'%':>8}")
    shown = 0
    for name, count in sorted(frames.items(), key=lambda kv: -kv[1]):
        if name == 'import':
            continue
        print(f'{name[:60]:<62}{count:>9}{100.0 * count / total:>7.1f}%')
        shown += 1
        if shown == 20:
            break
    print()
    print('  Leaf = self time: these are what to make faster. On stack =')
    print('  inclusive time: use it to see which stage owns the cost.')
    print('  "import" collapses interpreter startup, which is not the model.')


def section_allocations(tasks):
    traced = [t for t in tasks if t['allocs']]
    if not traced:
        return
    heading('5. ALLOCATION SITES       tracemalloc, live bytes by source line')
    pooled = defaultdict(lambda: [0.0, 0])
    for task in traced:
        for megabytes, count, source in task['allocs']:
            entry = pooled[source]
            entry[0] = max(entry[0], megabytes)
            entry[1] = max(entry[1], count)
    print(f"{'source':<56}{'MB':>10}{'objects':>12}")
    for source, (megabytes, count) in sorted(pooled.items(),
                                             key=lambda kv: -kv[1][0])[:25]:
        print(f'{source[:54]:<56}{megabytes:>10.2f}{count:>12}')
    print()
    print(f'  From {len(traced)} process(es) run with')
    print('  SIMPLICITY_PROFILE_TRACEMALLOC=1. This attributes LIVE bytes to')
    print('  the line that allocated them, which is the direct answer to "what')
    print('  should I change". It is a diagnosis-run tool: tracemalloc roughly')
    print('  doubles memory and slows every allocation, so these processes are')
    print('  not comparable on time with the untraced ones.')


def section_slowest(tasks):
    heading('6. SLOWEST PROCESSES       wall, CPU, and cost per simulated day')
    print(f"{'task':<18}{'scenario':<22}{'wall s':>8}{'cpu s':>8}{'cpu/wall':>9}"
          f"{'days':>7}{'s/day':>8}{'ms/react':>10}  exit")
    rows = []
    for task in tasks:
        samples = task['samples']
        wall = task['meta'].get('wall_s') or samples[-1].get('wall_s')
        cpu = task['meta'].get('cpu_s') or samples[-1].get('cpu_s')
        if not isinstance(cpu, (int, float)):
            cpu = None
        reached = max((s['sim_time'] for s in live(task)), default=None)
        reactions = max((s['reactions'] for s in live(task)
                         if s.get('reactions')), default=None)
        rows.append((wall or 0, task, wall, cpu, reached, reactions))
    for _, task, wall, cpu, reached, reactions in sorted(rows,
                                                         key=lambda r: -r[0])[:25]:
        ratio = cpu / wall if cpu and wall else None
        per_day = wall / reached if wall and reached else None
        per_reaction = 1000.0 * wall / reactions if wall and reactions else None
        flag = 'ok' if task['meta'].get('exited_cleanly') else 'KILLED/LIVE'
        print(f"{task['tag'][:16]:<18}{(task['scenario'] or '')[:20]:<22}"
              f"{fmt(wall, '.0f'):>8}{fmt(cpu, '.0f'):>8}{fmt(ratio, '.2f'):>9}"
              f"{fmt(reached, '.0f'):>7}{fmt(per_day, '.2f'):>8}"
              f"{fmt(per_reaction, '.2f'):>10}  {flag}")
    print()
    print('  s/day separates "ran long because it had many days to simulate"')
    print('  from "ran long because each day got dearer"; ms/react does the')
    print('  same against the event count, which is the quantity extrande')
    print('  actually controls.')
    print()
    print('  cpu/wall above 1 means threads: with --cpus-per-task=1 that is')
    print('  BLAS oversubscribing the allocation, which Slurm will account for')
    print('  and the node will contend over. Set OMP_NUM_THREADS=1 to test it.')


def main():
    profile_dir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
        'SIMPLICITY_PROFILE_DIR', 'profiles')
    data_dir = sys.argv[2] if len(sys.argv) > 2 else 'Data'
    if not os.path.isdir(profile_dir):
        sys.exit(f'no profile directory: {profile_dir}')
    tasks = load(profile_dir)
    if not tasks:
        sys.exit(f'no profile CSVs under {profile_dir}')
    join_scenarios(tasks, data_dir)
    limit = mem_limit_mb(tasks)

    rule('=')
    print('SIMULATION PROFILE - memory and time per process')
    rule('=')
    print(f'profile dir : {profile_dir}')
    print(f'data dir    : {data_dir}   (for scenario, seed and final_time)')
    section_completeness(tasks, limit)
    section_memory(tasks, limit)
    section_breakdown(tasks)
    section_fit(tasks)
    section_time(tasks)
    section_allocations(tasks)
    section_slowest(tasks)
    rule('=')


if __name__ == '__main__':
    main()
