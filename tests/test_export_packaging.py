#!/usr/bin/env python3
'''What a run exports must be what that run produced -- nothing else.

    python tests/test_export_packaging.py

WHY

Two bugs of the same shape, both visible on run #2.

1. The artifacts zip walked CONSENSUS_MODES and looked for a sanity grid per
   mode, while the orchestrator has refused to run more than one mode since
   v2.4.93. Every distribution run therefore printed

       [skip] artifact not found: Data/impact_long_shedders_amx_sanity_#2/...

   for a file nothing had been asked to write. Cosmetic, but it trains the eye
   to read skip lines as normal -- and a skip line is how a genuinely missing
   figure would announce itself.

2. The export handed package_simplicity_data.sh `Data` + `_#<n>` and let it
   glob. `Data/*_#2` also matches `nsr_range_finder_test_s1_#2`, so a probe
   directory rode along inside the run's tarball. A run's directories are
   knowable; run_export_dirs names them.

Both are the same mistake: deriving a set by matching when it could be stated.
So the checks here are about the SET -- what goes in and what stays out -- not
about whether the packer compresses.

The last section runs the real shell script over a tree that contains a decoy,
then runs the old glob over the same tree to show the decoy really would have
been swept in. It needs pixz; without it that section reports SKIP and the
rest still runs.
'''
import io
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'scripts'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'experiments'))

import impact_long_shedders_config as cfg
import run_impact_long_shedders_pipeline as pipeline

PKG_SCRIPT = os.path.join(REPO, 'scripts', 'package_simplicity_data.sh')
failures = []


def check(label, got, want):
    ok = got == want
    print(f"  [{'ok ' if ok else 'FAIL'}] {label}: {got!r}")
    if not ok:
        failures.append(f'{label}: got {got!r}, want {want!r}')


def plant(root, name, *files):
    d = os.path.join(root, 'Data', name)
    os.makedirs(d, exist_ok=True)
    for f in files or ('marker.txt',):
        p = os.path.join(d, f)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, 'w') as fh:
            fh.write('x')
    return d


# ===========================================================================
def test_export_dirs():
    print('\n=== run_export_dirs names the one directory a run writes')
    for mode, root in (('distribution', 'Data/impact_long_shedders_#7'),
                       ('argmax', 'Data/impact_long_shedders_amx_#7')):
        dirs = cfg.run_export_dirs(7, mode)
        check(f'{mode}: one directory', len(dirs), 1)
        # it used to name the experiment with a bare EXP_NAME, so an argmax
        # run exported a distribution directory that was never written
        check(f'{mode}: the run directory', dirs[0], root)
        check(f'{mode}: setup data is inside it',
              cfg.setup_dir(7, mode).startswith(root + os.sep), True)
        check(f'{mode}: so is the sanity grid',
              cfg.sanity_plot_path(7, mode).startswith(root + os.sep), True)
    # the whole point: nothing in the list is a pattern
    globby = [d for d in cfg.run_export_dirs(7) if set('*?[') & set(d)]
    check('no wildcards anywhere in the list', globby, [])
    check('default mode is the default consensus',
          cfg.run_export_dirs(7), cfg.run_export_dirs(7, cfg.DEFAULT_CONSENSUS))


# ===========================================================================
def test_artifacts_archive():
    print('\n=== the artifacts zip collects the mode that ran, and only that')
    root = tempfile.mkdtemp(prefix='simp_export_')
    cwd = os.getcwd()
    try:
        os.chdir(root)
        num = 3
        exp = cfg.run_exp_name(num)
        # BOTH modes fully on disk, so the mode asked for is the only thing
        # that can decide which artifacts are picked up. Each mode's sanity
        # grid now nests inside that mode's own run directory.
        for mode in cfg.CONSENSUS_MODES:
            run = cfg.run_exp_name(num, mode)
            plant(root, run,
                  f'05_Plots/{run}_long_nsr_calibration_fit.png',
                  f'05_Plots/{run}_std_nsr_calibration_fit.png')
            rel = os.path.relpath(cfg.sanity_plot_path(num, mode), 'Data')
            plant(root, os.path.dirname(rel), os.path.basename(rel))
        log = os.path.join(root, 'pipeline.log')
        with open(log, 'w') as fh:
            fh.write('log')

        for mode in cfg.CONSENSUS_MODES:
            other = [m for m in cfg.CONSENSUS_MODES if m != mode][0]
            zpath = os.path.join(root, f'artifacts_{mode}.zip')
            pipeline.write_artifacts_archive(zpath, num, log, [mode])
            names = sorted(zipfile.ZipFile(zpath).namelist())
            check(f'{mode}: four artifacts', len(names), 4)
            check(f'{mode}: its own grid is in',
                  os.path.basename(cfg.sanity_plot_path(num, mode)) in names,
                  True)
            check(f'{mode}: the other grid is NOT in',
                  os.path.basename(cfg.sanity_plot_path(num, other)) in names,
                  False)
            check(f'{mode}: both calibration fits are in',
                  sum(1 for n in names if n.endswith('_calibration_fit.png')), 2)

        print('\n=== a genuinely missing figure is still reported')
        os.remove(os.path.join(root, 'Data', exp, '05_Plots',
                               f'{exp}_std_nsr_calibration_fit.png'))
        zpath = os.path.join(root, 'artifacts_missing.zip')
        code = ('import sys, os\n'
                'sys.path[:0] = %r\n'
                'import run_impact_long_shedders_pipeline as p\n'
                'os.chdir(%r)\n'
                'p.write_artifacts_archive(%r, %d, %r, ["distribution"])\n'
                % (sys.path[:3], root, zpath, num, log))
        r = subprocess.run([sys.executable, '-c', code],
                           capture_output=True, text=True)
        skips = [l for l in r.stdout.splitlines() if 'artifact not found' in l]
        check('one skip line, for the file really missing', len(skips), 1)
        check('and it names the std fit',
              'std_nsr_calibration_fit' in (skips[0] if skips else ''), True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(root, ignore_errors=True)


# ===========================================================================
def test_compress_command():
    print('\n=== the export command names directories; no SOURCE_DIR/SUFFIX')
    captured = {}

    class FakeResult:
        returncode = 0
        stdout = 'Submitted batch job 1'

    def fake_run(cmd, **kw):
        captured['cmd'] = cmd
        captured['env'] = kw.get('env')
        return FakeResult()

    log = io.StringIO()
    real_run = pipeline.subprocess.run
    pipeline.subprocess.run = fake_run
    try:
        for mode in cfg.CONSENSUS_MODES:
            pipeline.compress_results(4, log, mode, keep_slurm_logs=True)
            cmd = captured['cmd']
            check(f'{mode}: calls the packer',
                  os.path.basename(cmd[0]), 'package_simplicity_data.sh')
            dirs = [cmd[i + 1] for i, a in enumerate(cmd) if a == '--dir']
            check(f'{mode}: --dir list is run_export_dirs',
                  dirs, cfg.run_export_dirs(4, mode))
            check(f'{mode}: tagged with the experiment',
                  cmd[cmd.index('--tag') + 1], f'{cfg.EXP_NAME}_#4')
            check(f'{mode}: submitted to slurm', '--slurm' in cmd, True)
            # the old call was [script, "Data", "_#4", "--slurm"]
            check(f'{mode}: no bare suffix argument', '_#4' in cmd, False)
            check(f'{mode}: no bare source directory', 'Data' in cmd, False)
        check('slurm logs kept when asked',
              captured['env']['SIMPLICITY_KEEP_SLURM_LOGS'], '1')
        pipeline.compress_results(4, log, cfg.DEFAULT_CONSENSUS,
                                  keep_slurm_logs=False)
        check('and dropped when not',
              captured['env']['SIMPLICITY_KEEP_SLURM_LOGS'], '0')

        # carried over from tests/test_parallel_dispatch.py, retired with the
        # per-scenario threading it was written for: the data is on disk and
        # intact whatever the packer does, so a packing failure is reported
        # rather than raised -- it must not read as a failed run.
        class Failed:
            returncode = 1
            stdout = 'pixz not found'

        pipeline.subprocess.run = lambda cmd, **kw: Failed()
        out = pipeline.compress_results(4, log, cfg.DEFAULT_CONSENSUS)
        check('a packing failure is reported, not raised', out, None)
        check('and the log says how to pack it by hand',
              '--dir' in log.getvalue() and 'packaging exited 1' in log.getvalue(),
              True)
    finally:
        pipeline.subprocess.run = real_run


# ===========================================================================
def _top_level(archive):
    out = subprocess.run(['tar', '-I', 'pixz', '-tf', archive],
                         capture_output=True, text=True)
    return sorted({p.split('/')[1] for p in out.stdout.split()
                   if p.startswith('Data/') and len(p.split('/')) > 1})


def test_packer_end_to_end():
    print('\n=== the real packer over a tree holding a decoy')
    if shutil.which('pixz') is None:
        print('  [SKIP] pixz not installed; the shell packer cannot run here')
        return
    root = tempfile.mkdtemp(prefix='simp_pack_')
    cwd = os.getcwd()
    try:
        os.chdir(root)
        os.makedirs('scripts', exist_ok=True)
        shutil.copy(PKG_SCRIPT, 'scripts/')
        wanted = sorted(os.path.basename(d) for d in cfg.run_export_dirs(2))
        # the decoy is real: tests/run_nsr_range_finder.sh writes this name
        decoy = 'nsr_range_finder_test_s1_#2'
        for name in wanted + [decoy, f'{cfg.EXP_NAME}_#20']:
            plant(root, name, '05_Plots/a.png')

        cmd = ['bash', 'scripts/package_simplicity_data.sh']
        for d in cfg.run_export_dirs(2):
            cmd += ['--dir', d]
        cmd += ['--tag', f'{cfg.EXP_NAME}_#2']
        env = dict(os.environ, SIMPLICITY_KEEP_SLURM_LOGS='1')
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        check('packer succeeded', r.returncode, 0)
        archives = sorted(os.listdir('Data_Export'))
        check('one archive', len(archives), 1)
        check('named for the experiment, not just the date',
              f'{cfg.EXP_NAME}_#2' in archives[0], True)
        packed = _top_level(os.path.join('Data_Export', archives[0]))
        check('packed exactly the run directories', packed, wanted)
        check('the probe directory stayed out', decoy in packed, False)
        check('and so did _#20', f'{cfg.EXP_NAME}_#20' in packed, False)

        print('\n  -- negative control: the glob this replaced --')
        shutil.rmtree('Data_Export')
        subprocess.run(['bash', 'scripts/package_simplicity_data.sh',
                        'Data', '_#2'], capture_output=True, text=True, env=env)
        globbed = _top_level(os.path.join(
            'Data_Export', sorted(os.listdir('Data_Export'))[0]))
        check('the old glob DID sweep the probe in', decoy in globbed, True)
        check('which is the bug this test exists for', globbed != packed, True)
    finally:
        os.chdir(cwd)
        shutil.rmtree(root, ignore_errors=True)


# ===========================================================================
if __name__ == '__main__':
    test_export_dirs()
    test_artifacts_archive()
    test_compress_command()
    test_packer_end_to_end()
    print('\n' + '=' * 70)
    if failures:
        print(f'{len(failures)} FAILED')
        for f in failures:
            print(f'  - {f}')
        sys.exit(1)
    print('all checks passed')
