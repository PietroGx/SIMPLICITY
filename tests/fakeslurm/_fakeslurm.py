#!/usr/bin/env python3
'''A Slurm controller small enough to lie about on purpose.

sbatch / squeue / scontrol / sacct shims that implement just enough of a held
job array for simplicity.runners.slurm to drive end to end. Put this directory
first on PATH and the real control plane runs unmodified: it really does submit,
poll, release, resolve and reconcile.

WHY A FAKE RATHER THAN A REAL LOCAL SLURM

A single-node slurmctld would exercise submit/hold/release honestly, but the
branches that have actually cost production runs are the ones you cannot ask a
real controller to produce on demand:

    TIMEOUT / OUT_OF_MEMORY          a walltime or memory kill, which terminates
                                     the task's process so job() never reaches
                                     its except block
    launch failed requeued held      Slurm auto-requeueing a task it could not
                                     start
    COMPLETED with no state written  profile_grid_#910: 120 tasks exited 0,
                                     sacct said COMPLETED for all of them, and
                                     only 116 had written their own terminal
                                     state. The other four matched no branch at
                                     all, and the run polled for hours.

Here each of those is one word in a plan file. A real controller is still worth
having -- see tests/slurm_local_setup.md -- but as a check that these shims do
not lie about the INTERFACE, not as a way to reach these cases.

HOW IT BEHAVES

Deterministic and synchronous: `scontrol release` runs the released tasks to
completion before returning, so the caller's next poll already sees the result.
No daemon, no clock to race.

State lives in $FAKESLURM_STATE (a json file). Per-task behaviour comes from
$FAKESLURM_PLAN, a json object keyed by array task id (1-based):

    ok            run the task; it records its own completion
    raise         the work raises; job() records failed
    hang          the work blocks and the task is killed -- sacct then reports
                  $FAKESLURM_KILL_STATE (default TIMEOUT) and the repeat is left
                  in `started`, exactly as an external kill leaves it
    silent        run the task, then roll its terminal state back to `started`
                  while sacct still reports COMPLETED -- the #910 case
    launch_fail   refuse to start it, with the reason Slurm uses, until it has
                  been released $FAKESLURM_LAUNCH_FAILS times (default 2)
'''
import json
import os
import signal
import subprocess
import sys
import time

JOB_ID = "4242"
KILL_AFTER_S = 1.5


def _state_path():
    path = os.environ.get("FAKESLURM_STATE")
    if not path:
        sys.stderr.write("FAKESLURM_STATE is not set\n")
        raise SystemExit(2)
    return path


def _load():
    try:
        with open(_state_path()) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return {"job_id": JOB_ID, "name": None, "script": None, "env": {},
                "tasks": {}}


def _save(state):
    path = _state_path()
    tmp = f"{path}.tmp"
    with open(tmp, "w") as handle:
        json.dump(state, handle, indent=1)
    os.replace(tmp, path)


def _plan(task_id):
    try:
        plan = json.loads(os.environ.get("FAKESLURM_PLAN", "{}"))
    except ValueError:
        plan = {}
    return plan.get(str(task_id), "ok")


def _arg(args, name, default=None):
    """--name=value or --name value."""
    for index, arg in enumerate(args):
        if arg.startswith(name + "="):
            return arg.split("=", 1)[1]
        if arg == name and index + 1 < len(args):
            return args[index + 1]
    return default


# ------------------------------------------------------------------- sbatch

def sbatch(args):
    if "--version" in args:
        print("slurm 00.00.0-fake")      # answered before reading stdin
        return 0
    script = sys.stdin.read()
    array = _arg(args, "--array", "1-1")
    low, _, high = array.partition("-")
    state = _load()
    state["job_id"] = JOB_ID
    state["name"] = _arg(args, "--job-name")
    state["script"] = script
    # the environment sbatch was handed IS how the task learns its experiment,
    # its submission and the function to call
    state["env"] = {key: value for key, value in os.environ.items()
                    if key.startswith("SIMPLICITY") or key in
                    ("USER_RUN_SEEDED_SIMULATION", "PYTHONPATH",
                     "PYTHONHASHSEED", "FAKESLURM_STATE", "FAKESLURM_PLAN")}
    state["tasks"] = {
        str(task): {"state": "PENDING", "held": True, "reason": "(JobHeldUser)",
                    "releases": 0, "sacct": None}
        for task in range(int(low), int(high or low) + 1)}
    _save(state)
    # --parsable: the job id alone, which is what check_slurm_interface reads
    print(JOB_ID if "--parsable" in args else f"Submitted batch job {JOB_ID}")
    return 0


# -------------------------------------------------------------------- squeue

def _collapse(tasks):
    """[1,2,3,7] -> "1-3,7", the way squeue writes a pending array."""
    parts, start, previous = [], None, None
    for task in tasks + [None]:
        if previous is not None and task == previous + 1:
            previous = task
            continue
        if start is not None:
            parts.append(str(start) if start == previous else f"{start}-{previous}")
        start = previous = task
    return ",".join(parts)


def squeue(args):
    state = _load()
    fmt = _arg(args, "--Format", "") or ""
    wanted = _arg(args, "--name")
    if wanted and state.get("name") not in (None, wanted):
        return 0
    queued = [task for task, info in sorted(state["tasks"].items(),
                                            key=lambda kv: int(kv[0]))
              if info["state"] in ("PENDING", "RUNNING")]
    if "ArrayTaskID" in fmt:
        # launch-failure detection: job id, task id(s), free-text reason.
        #
        # Slurm collapses pending tasks that share a state AND a reason into
        # ONE row, with the task field as a range expression ("1-3", "1,3,5").
        # These shims used to emit one row per task, which is precisely the
        # kind of lie a fake can tell: test_slurm_lifecycle passed while the
        # real parser's int() rejected every row on slurm 26.05.4, and only
        # check_slurm_interface.py against a real controller found it.
        by_reason = {}
        for task in queued:
            info = state["tasks"][task]
            by_reason.setdefault((info["state"], info["reason"]), []).append(int(task))
        for (_state, reason), tasks in by_reason.items():
            print(f'{state["job_id"]} {_collapse(sorted(tasks))} {reason}')
        return 0
    # release path: the real code drops line 0 as a header and expects exactly
    # one distinct ArrayJobID across the rest
    if "--noheader" not in args:
        print("ARRAY_JOB_ID")
    for _ in queued:
        print(state["job_id"])
    return 0


# ------------------------------------------------------------------ scontrol

def _run_task(state, task):
    """Start one task the way Slurm would, then apply its planned fate."""
    info = state["tasks"][task]
    behaviour = _plan(task)

    if behaviour == "launch_fail":
        limit = int(os.environ.get("FAKESLURM_LAUNCH_FAILS", "2"))
        if info["releases"] <= limit:
            info["state"] = "PENDING"
            info["held"] = True
            info["reason"] = "(launch failed requeued held)"
            return
        behaviour = "ok"        # it eventually starts

    script = state.get("script") or ""
    interpreter = sys.executable
    if script.startswith("#!"):
        interpreter = script.splitlines()[0][2:].strip() or interpreter
    body = "\n".join(script.splitlines()[1:])

    env = {**os.environ, **state.get("env", {}),
           "SLURM_ARRAY_JOB_ID": state["job_id"],
           "SLURM_ARRAY_TASK_ID": str(task),
           "SLURM_ARRAY_TASK_MIN": "1",
           "FAKEWORK_BEHAVIOUR": behaviour}

    info["state"] = "RUNNING"
    info["held"] = False
    info["reason"] = "None"

    process = subprocess.Popen([interpreter, "-c", body], env=env,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
    if behaviour == "hang":
        # a walltime or memory kill: SIGKILL, so the task's own process never
        # runs job()'s except block and the repeat stays `started`
        time.sleep(KILL_AFTER_S)
        process.send_signal(signal.SIGKILL)
        process.wait()
        info["state"] = os.environ.get("FAKESLURM_KILL_STATE", "TIMEOUT")
        info["sacct"] = info["state"]
        return

    returncode = process.wait()
    info["state"] = "COMPLETED" if returncode == 0 else "FAILED"
    info["sacct"] = info["state"]

    if behaviour == "silent":
        # the task ran and exited 0, but its terminal state write was lost
        _rollback_to_started(state, task)
        info["sacct"] = "COMPLETED"


def _rollback_to_started(state, task):
    """Leave one repeat's state file saying `started` -- the #910 case.

    Written straight to the file, NOT through jobs.set_state: a lost write means
    the terminal content never reached the file, and set_state is monotonic now,
    so asking it to go backwards is correctly refused. Simulating the loss
    through the API would simulate nothing.
    """
    env = state.get("env", {})
    experiment = env.get("SIMPLICITY_EXPERIMENT_NAME")
    repo = env.get("FAKESLURM_REPO") or os.environ.get("FAKESLURM_REPO")
    if not experiment or not repo:
        return
    sys.path.insert(0, repo)
    try:
        import simplicity.dir_manager as dm
        import simplicity.jobs as jobs
        data_dir = env.get("FAKESLURM_DATA_DIR") or os.environ.get(
            "FAKESLURM_DATA_DIR")
        if data_dir:
            dm.set_data_dir(data_dir)
        submission = env.get("SIMPLICITY_SUBMISSION")
        group, record = jobs.resolve_task(experiment, submission, int(task) - 1)
        path = jobs.state_path(experiment, group, record["index"])
        with open(path) as handle:
            record_on_disk = json.load(handle)
        record_on_disk["state"] = jobs.STARTED
        with open(path, "w") as handle:
            json.dump(record_on_disk, handle)
    except Exception:
        pass


def scontrol(args):
    if not args or args[0] != "release":
        return 0
    state = _load()
    for item in args[1].split(","):
        _, _, task = item.partition("_")
        if task in state["tasks"]:
            state["tasks"][task]["releases"] += 1
    _save(state)
    # run them one at a time, re-reading nothing: synchronous by design, so the
    # caller's next poll already sees the outcome
    for item in args[1].split(","):
        _, _, task = item.partition("_")
        if task in state["tasks"]:
            _run_task(state, task)
            _save(state)
    return 0


# --------------------------------------------------------------------- sacct

def sacct(args):
    state = _load()
    job_id = _arg(args, "-j", state.get("job_id"))
    if job_id != state.get("job_id"):
        return 0
    for task, info in sorted(state["tasks"].items(), key=lambda kv: int(kv[0])):
        reported = info.get("sacct")
        if reported:
            print(f'{state["job_id"]}_{task}|{reported}')
    return 0


def scancel(args):
    state = _load()
    for info in state["tasks"].values():
        if info["state"] in ("PENDING", "RUNNING"):
            info["state"] = "CANCELLED"
            info["sacct"] = info.get("sacct") or "CANCELLED"
    _save(state)
    return 0


def main(tool):
    return {"sbatch": sbatch, "squeue": squeue, "scontrol": scontrol,
            "sacct": sacct, "scancel": scancel}[tool](sys.argv[1:])
