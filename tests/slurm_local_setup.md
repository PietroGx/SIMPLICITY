# A real Slurm on this machine

Commands for **you** to run — they need `sudo` and start system services, so
they are written out rather than executed.

Target: Ubuntu 24.04, `slurm-wlm 23.11.4`. One node, acting as controller,
database and compute.

## Why bother, given the fake exists

`tests/test_slurm_lifecycle.py` already drives the whole control plane against
`tests/fakeslurm/`, and it reaches cases a real controller will not produce on
demand (`OUT_OF_MEMORY`, `launch failed requeued held`, a lost terminal state
write). It found a real bug doing so — see the note at the end.

What the fake cannot tell you is whether it **lies about the interface**: whether
real `squeue` output parses the way the shims pretend, whether `sacct --parsable2
-X` is really one line per array task, whether `scontrol release` really accepts
a comma-joined `jobid_taskid` list. Those are the assumptions the production code
rests on, and only a real controller settles them.

So: the fake is the test, this is the interface check. Run the pipeline through
it once after a change to `simplicity/runners/slurm.py`.

## 1. Install

```bash
sudo apt update
sudo apt install -y slurm-wlm slurmdbd mariadb-server munge
```

`slurmdbd` and `mariadb-server` are **not optional** for this to be useful.
`reconcile_terminated_tasks` reads terminal-state history from `sacct`, and with
the default `AccountingStorageType=accounting_storage/none` that history does not
exist — `sacct` returns nothing and the single most important branch in the
control plane is untestable. That branch is what unblocked profile_grid_#910.

## 2. Accounting database

```bash
sudo mysql -e "CREATE DATABASE IF NOT EXISTS slurm_acct_db;"
sudo mysql -e "CREATE USER IF NOT EXISTS 'slurm'@'localhost' IDENTIFIED BY 'slurmdbpass';"
sudo mysql -e "GRANT ALL ON slurm_acct_db.* TO 'slurm'@'localhost'; FLUSH PRIVILEGES;"
```

`/etc/slurm/slurmdbd.conf` — **must** be mode 0600 and owned by `slurm`, or
`slurmdbd` refuses to start:

```ini
AuthType=auth/munge
DbdHost=localhost
DbdPort=6819
SlurmUser=slurm
StorageType=accounting_storage/mysql
StorageHost=localhost
StorageUser=slurm
StoragePass=slurmdbpass
StorageLoc=slurm_acct_db
LogFile=/var/log/slurm/slurmdbd.log
PidFile=/run/slurmdbd.pid
```

```bash
sudo chown slurm:slurm /etc/slurm/slurmdbd.conf
sudo chmod 600 /etc/slurm/slurmdbd.conf
sudo mkdir -p /var/log/slurm /var/spool/slurmctld /var/spool/slurmd
sudo chown slurm:slurm /var/log/slurm /var/spool/slurmctld
```

## 3. Controller

`/etc/slurm/slurm.conf`. Replace `CPUS` with `nproc` and `MEM` with a little
under your RAM in MB — Slurm refuses to start a node that claims more than it
has:

```ini
ClusterName=simplicity
SlurmctldHost=localhost
SlurmUser=slurm
AuthType=auth/munge
StateSaveLocation=/var/spool/slurmctld
SlurmdSpoolDir=/var/spool/slurmd
SlurmctldPidFile=/run/slurmctld.pid
SlurmdPidFile=/run/slurmd.pid
ProctrackType=proctrack/linuxproc
TaskPlugin=task/none

# cons_tres so several array tasks share the node, which is the whole point:
# a one-task-at-a-time node never exercises the release cap.
SelectType=select/cons_tres
SelectTypeParameters=CR_CPU_Memory

# this is what makes sacct retain terminal states
AccountingStorageType=accounting_storage/slurmdbd
AccountingStorageHost=localhost
JobAcctGatherType=jobacct_gather/linux

# short, so a walltime kill is reachable in a test rather than tomorrow
MinJobAge=300
KillWait=5

NodeName=localhost CPUs=CPUS RealMemory=MEM State=UNKNOWN
PartitionName=debug Nodes=ALL Default=YES MaxTime=INFINITE State=UP
```

```bash
sudo systemctl enable --now munge slurmdbd
sleep 3                       # slurmdbd must be up before slurmctld registers
sudo systemctl enable --now slurmctld slurmd
sinfo                         # the node should read idle, not down/drain
sacctmgr -i add cluster simplicity
```

If the node comes up `drain`, `scontrol update nodename=localhost state=resume`
and check `RealMemory` against `free -m`.

## 4. Check the interface, not the logic

```bash
cd ~/Documents/SIMPLICITY
scripts/start_simplicity_session.sh         # tmux + conda + SBATCH_QOS
```

Then a small real experiment, slurm runner:

```bash
python - <<'PY'
import simplicity.runme as runme
import simplicity.runners.slurm as slurm
runme.run_experiment(
    'slurm_iface_check_#1',
    lambda: ({'R': [1.05, 1.2]},
             {'population_size': 40, 'final_time': 5,
              'infected_individuals_at_start': 4}, 2),
    simplicity_runner=slurm, archive_experiment=False)
PY
```

Four repeats, two simulations, a held array released under the cap. What to
confirm — these are the shim assumptions, in order of how much they would cost
if wrong:

| Check | How |
|---|---|
| every array position ran its own repeat, exactly once | `04_Output/main/sim_*/seed_*/` — 4 directories, each with its CSVs |
| `sacct` really retains terminal state | `sacct -j <jobid> --format=JobID,State --noheader --parsable2 -X` returns one line per task |
| the id map matches | `Data/slurm_iface_check_#1/slurm/job_id_mapping/*.csv` each contain `main/<index>` |
| the loop terminated on its own | it returned without the "no held task" exception |
| release honoured the cap | `SIMPLICITY_MAX_PARALLEL_SEEDED_SIMULATIONS_SLURM=2` and never more than 2 running in `squeue` |

To reach a walltime kill for real: `SIMPLICITY_SLURM_TIME=00:00:30` with
`final_time=1095`. `sacct` should report `TIMEOUT`, and
`reconcile_terminated_tasks` should mark the repeat failed within
`RECONCILE_INTERVAL_S` (900s — set it lower in the session if you do not want to
wait).

## The bug the fake found

Worth knowing, because it is the reason this file exists in this shape.

`release_simulations` records `RELEASED` for a batch **after** `scontrol release`
returns. Under the old scheme that was harmless: `.released` and `.completed`
were separate files, so touching one next to the other changed nothing. With one
state field per repeat, a task that managed to start *and finish* inside that
window would have its `COMPLETED` overwritten by `RELEASED` — and the polling
loop would then wait on it forever, with its output already on disk.

In production that window is milliseconds and needs a task to start faster than
a function returns, so it would have surfaced rarely and looked like the monitor
hangs already fixed twice. Under a synchronous fake controller it happens every
single time, which is how it was found within a minute.

`jobs.set_state` now refuses to move a repeat backwards (`jobs.STATE_RANK`). A
real Slurm would most likely never have shown this.
