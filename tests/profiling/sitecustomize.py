"""Per-process memory and time profiler for SIMPLICITY simulations.

Python imports `sitecustomize` automatically at interpreter startup if it is
importable, so putting this directory on PYTHONPATH profiles every Slurm task
without a single change to simplicity/ or scripts/. The chain works because
simplicity/runners/slurm.py:243 passes env={**os.environ, ...} to sbatch, and
sbatch exports its environment to the tasks by default.

    export PYTHONPATH=$PWD/tests/profiling:$PYTHONPATH
    export SIMPLICITY_PROFILE_DIR=$PWD/profiles

Unset SIMPLICITY_PROFILE_DIR and this module returns immediately: nothing is
sampled, no thread starts, no file is written.

Three samplers share one daemon thread:

  memory   every SIMPLICITY_PROFILE_INTERVAL_S (default 15): RSS, the kernel's
           own high-water mark, CPU time, and for every container on Population
           that can grow, both its length and an estimate of its bytes. RSS
           says how much, the lengths say which one is growing, and the byte
           estimate says whether that growth is big enough to be the cause.

  speed    the same row carries extrande's own event counters (leaps,
           reactions, thinning) off the live ProgressReporter, so cost can be
           read per event rather than per wall-clock second: a run that is slow
           because it took more steps is a different problem from one that is
           slow because each step got dearer.

  stacks   every SIMPLICITY_PROFILE_STACK_INTERVAL_S (default 0.25) the main
           thread's Python stack, counted in a folded histogram. This is a
           statistical profiler: cheap enough to leave on for a production run,
           and it attributes time to the Python frame that called into numpy
           rather than losing it inside the C call.

Every row is flushed and fsync'd as it is written, deliberately: a task Slurm
kills for exceeding --mem dies on SIGKILL and never runs its atexit handler, so
the partial curve on disk is the whole diagnosis.

Optional, off by default because it roughly doubles memory and slows every
allocation -- for a dedicated diagnosis run, never a production one:

    SIMPLICITY_PROFILE_TRACEMALLOC=1   top allocation sites by source line

Nothing here may break a run: every sampler body is wrapped, and a failure
disables profiling rather than propagating.
"""
import os
import sys

_DIR = os.environ.get("SIMPLICITY_PROFILE_DIR")

# Containers on Population whose length can grow with the run. Anything absent
# on a given version is reported blank rather than failing -- this file must
# never be the reason a production run dies.
_TRACKED = (
    "lineage_frequency",
    "individuals",
    "phylogenetic_data",
    "phylodots",
    "consensus_sequences_t",
    "trajectory",
    "fitness_trajectory",
    "R_effective_trajectory",
    "infected_i",
    "infectious_i",
    "recovered_i",
    "reservoir_i",
)

_COLUMNS = (
    ("wall_s", "cpu_s", "rss_mb", "hwm_mb", "sim_time", "infected",
     "active_lineages", "leaps", "reactions", "thinning", "gc_collected",
     "traced_mb")
    + tuple("n_" + name for name in _TRACKED)
    + tuple("mb_" + name for name in _TRACKED)
    + ("n_consensus_positions", "mb_consensus")
)

# How many elements of a container to measure when estimating its size. The
# rows of these containers are uniform, so a handful is enough; sampling both
# ends guards against a container whose early rows differ from its late ones.
_SIZE_SAMPLE = 24


def _load_shadowed_sitecustomize():
    """If the interpreter already had a sitecustomize, this file shadows it.
    Execute that one too, so this only adds behaviour."""
    here = os.path.dirname(os.path.abspath(__file__))
    for entry in sys.path:
        try:
            candidate = os.path.join(os.path.abspath(entry or "."), "sitecustomize.py")
            if os.path.dirname(candidate) == here or not os.path.isfile(candidate):
                continue
            import importlib.util
            spec = importlib.util.spec_from_file_location("_shadowed_sitecustomize",
                                                          candidate)
            spec.loader.exec_module(importlib.util.module_from_spec(spec))
        except Exception:
            pass
        return


def _approx_size(obj, depth, seen):
    """Shallow-ish deep size, counting each distinct object once.

    `seen` carries object ids across a whole measurement, so the interned key
    strings shared by 100,000 identical dicts are charged once rather than
    100,000 times -- which is the difference between a useful number and a
    wild over-count.
    """
    key = id(obj)
    if key in seen:
        return 0
    seen.add(key)
    size = sys.getsizeof(obj)
    if depth <= 0:
        return size
    try:
        if isinstance(obj, dict):
            for item_key, value in obj.items():
                size += _approx_size(item_key, depth - 1, seen)
                size += _approx_size(value, depth - 1, seen)
        elif isinstance(obj, (list, tuple, set, frozenset)):
            for item in obj:
                size += _approx_size(item, depth - 1, seen)
    except RuntimeError:
        pass   # mutated under us by the simulation thread; the estimate stands
    return size


def _elements(container, limit):
    """Up to `limit` elements, taken from both ends so a container whose early
    rows differ from its late ones is not misread from the front alone."""
    half = max(1, limit // 2)
    if isinstance(container, dict):
        keys = list(container)
        picks = keys[:half] + keys[-half:]
        return [(k, container[k]) for k in picks if k in container]
    if isinstance(container, (list, tuple)):
        return list(container[:half]) + list(container[-half:])
    out = []
    for index, item in enumerate(container):
        if index >= limit:
            break
        out.append(item)
    return out


def _measure(items, seen):
    total = 0
    for item in items:
        if isinstance(item, tuple) and len(item) == 2:
            total += _approx_size(item[0], 1, seen) + _approx_size(item[1], 2, seen)
        else:
            total += _approx_size(item, 2, seen)
    return total


def _container_mb(container):
    """Estimated megabytes held by a container.

    Measures two nested samples and extrapolates on the MARGINAL cost of an
    element -- what one more row adds -- rather than on the average. Structure
    shared across rows is paid for once by the first sample and then correctly
    excluded, so a container of near-identical dicts is not charged for its
    keys over and over.
    """
    try:
        count = len(container)
    except TypeError:
        return ""
    base = sys.getsizeof(container)
    if count == 0:
        return round(base / 1048576.0, 3)
    try:
        picks = _elements(container, _SIZE_SAMPLE)
        if not picks:
            return round(base / 1048576.0, 3)
        split = max(1, len(picks) // 2)
        seen = set()
        first = _measure(picks[:split], seen)
        second = _measure(picks[split:], seen)
    except (RuntimeError, KeyError, IndexError):
        return ""
    sampled = len(picks)
    remaining = sampled - split
    if remaining > 0 and second > 0:
        marginal = second / float(remaining)
    else:
        marginal = (first + second) / float(sampled)
    total = base + first + second + marginal * max(0, count - sampled)
    return round(total / 1048576.0, 3)


def _start():
    import csv
    import gc
    import json
    import threading
    import time
    import weakref

    interval = float(os.environ.get("SIMPLICITY_PROFILE_INTERVAL_S", "15"))
    stack_interval = float(os.environ.get("SIMPLICITY_PROFILE_STACK_INTERVAL_S", "0.25"))
    stacks_on = os.environ.get("SIMPLICITY_PROFILE_STACKS", "1") != "0"
    sizes_on = os.environ.get("SIMPLICITY_PROFILE_SIZES", "1") != "0"
    tracemalloc_on = os.environ.get("SIMPLICITY_PROFILE_TRACEMALLOC", "0") == "1"

    experiment = os.environ.get("SIMPLICITY_EXPERIMENT_NAME", "unknown")
    job = os.environ.get("SLURM_ARRAY_JOB_ID") or os.environ.get("SLURM_JOB_ID")
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    # the name carries job and task so profile_report.py can join a row back to
    # its scenario and seed through the experiment's slurm_id_map, which
    # slurm.job() writes under exactly these two ids.
    tag = f"{job}_{task}" if job and task else f"local_{os.getpid()}"

    out_dir = os.path.join(_DIR, experiment)
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, tag + ".csv")
    meta_path = os.path.join(out_dir, tag + ".meta.json")
    stack_path = os.path.join(out_dir, tag + ".folded")
    alloc_path = os.path.join(out_dir, tag + ".alloc")

    page = os.sysconf("SC_PAGE_SIZE")
    ticks = os.sysconf("SC_CLK_TCK")
    # "seen" latches once a Population has been built: every Python process
    # on PYTHONPATH loads this module, including the pipeline's own driver and
    # any helper script, and only the ones that actually run a simulation
    # should leave a file behind.
    state = {"pop": None, "reporter": None, "hooked": set(), "seen": False}
    stacks = {}
    started = time.time()

    if tracemalloc_on:
        import tracemalloc
        tracemalloc.start(1)

    def rss_mb():
        with open("/proc/self/statm") as handle:
            return int(handle.read().split()[1]) * page / 1048576.0

    def hwm_mb():
        # VmHWM is the kernel's own high-water mark, so a spike between two
        # samples is still recorded. This is the number to compare against
        # --mem and against sacct's MaxRSS.
        try:
            with open("/proc/self/status") as handle:
                for line in handle:
                    if line.startswith("VmHWM:"):
                        return int(line.split()[1]) / 1024.0
        except OSError:
            pass
        return ""

    def cpu_s():
        # utime+stime of the whole process. Compared against wall this says
        # whether a slow task is computing or waiting on the filesystem.
        try:
            with open("/proc/self/stat") as handle:
                fields = handle.read().rsplit(") ", 1)[1].split()
            return round((int(fields[11]) + int(fields[12])) / float(ticks), 2)
        except (OSError, IndexError, ValueError):
            return ""

    def register_population(cls):
        original = cls.__init__

        def patched(self, *args, **kwargs):
            original(self, *args, **kwargs)
            state["pop"] = weakref.ref(self)
            state["seen"] = True

        cls.__init__ = patched
        state["hooked"].add("population")

    def register_reporter(cls):
        original = cls.__init__

        def patched(self, *args, **kwargs):
            original(self, *args, **kwargs)
            state["reporter"] = weakref.ref(self)

        cls.__init__ = patched
        state["hooked"].add("extrande")

    # Which class to grab out of which module once that module finishes loading.
    _HOOKS = {"simplicity.population": ("Population", register_population),
              "simplicity.extrande": ("ProgressReporter", register_reporter)}

    def hook_now():
        """Fallback for a module that was already imported before this ran."""
        for module_name, (class_name, register) in _HOOKS.items():
            key = module_name.rsplit(".", 1)[1]
            if key in state["hooked"]:
                continue
            module = sys.modules.get(module_name)
            cls = getattr(module, class_name, None) if module else None
            if cls is not None:
                register(cls)

    def install_import_hook():
        """Patch each class the instant its module finishes executing.

        sitecustomize runs before simplicity is imported, so the classes cannot
        be patched up front. Polling for them in the sampler thread would race
        the first instance; wrapping the loader cannot.
        """
        import importlib.abc

        class Finder(importlib.abc.MetaPathFinder):
            def find_spec(self, fullname, path=None, target=None):
                entry = _HOOKS.get(fullname)
                if entry is None:
                    return None
                class_name, register = entry
                for finder in list(sys.meta_path):
                    if finder is self:
                        continue
                    found = finder.find_spec(fullname, path, target)
                    if found is not None:
                        break
                else:
                    return None
                loader = getattr(found, "loader", None)
                if loader is None:
                    return found
                original_exec = loader.exec_module

                def exec_module(module):
                    original_exec(module)
                    cls = getattr(module, class_name, None)
                    if cls is not None:
                        register(cls)

                try:
                    loader.exec_module = exec_module
                except (AttributeError, TypeError):
                    pass   # a loader that will not take the wrapper: hook_now covers it
                return found

        sys.meta_path.insert(0, Finder())

    def sample_row():
        row = {"wall_s": round(time.time() - started, 2),
               "cpu_s": cpu_s(),
               "rss_mb": round(rss_mb(), 1)}
        peak = hwm_mb()
        row["hwm_mb"] = round(peak, 1) if peak != "" else ""
        row["gc_collected"] = sum(s.get("collected", 0) for s in gc.get_stats())
        if tracemalloc_on:
            import tracemalloc
            # the interpreter's own count of live traced bytes: the one
            # independent check on the per-container estimates below
            row["traced_mb"] = round(tracemalloc.get_traced_memory()[0] / 1048576.0, 1)

        reporter = state["reporter"]() if state["reporter"] else None
        if reporter is not None:
            row["leaps"] = getattr(reporter, "leap_counter", "")
            row["reactions"] = getattr(reporter, "reactions_counter", "")
            row["thinning"] = getattr(reporter, "thinning_counter", "")

        population = state["pop"]() if state["pop"] else None
        if population is not None:
            row["sim_time"] = round(getattr(population, "time", 0) or 0, 3)
            row["infected"] = getattr(population, "infected", "")
            row["active_lineages"] = getattr(population, "active_lineages_n", "")
            for name in _TRACKED:
                container = getattr(population, name, None)
                if container is None:
                    continue
                try:
                    row["n_" + name] = len(container)
                except TypeError:
                    continue
                if sizes_on:
                    row["mb_" + name] = _container_mb(container)
            accumulator = getattr(population, "consensus", None)
            if accumulator is not None:
                positions = getattr(accumulator, "_s_e", None)
                if positions is not None:
                    row["n_consensus_positions"] = len(positions)
                    if sizes_on:
                        row["mb_consensus"] = _container_mb(positions)
        return row

    def sample_stack():
        frame = sys._current_frames().get(main_thread_id)
        if frame is None:
            return
        parts = []
        while frame is not None and len(parts) < 120:
            code = frame.f_code
            name = os.path.basename(code.co_filename)
            if name.startswith("<frozen importlib"):
                part = "import"          # a chain of these says "importing", nothing more
            elif name == "sitecustomize.py":
                frame = frame.f_back     # the profiler's own hook is not the program
                continue
            else:
                part = f"{name}:{code.co_name}"
            if not parts or parts[-1] != part:
                parts.append(part)       # collapses recursion and import chains
            frame = frame.f_back
        key = ";".join(reversed(parts))
        stacks[key] = stacks.get(key, 0) + 1

    def dump_stacks():
        tmp = stack_path + ".tmp"
        with open(tmp, "w") as handle:
            for key, count in sorted(stacks.items(), key=lambda kv: -kv[1]):
                handle.write(f"{key} {count}\n")
        os.replace(tmp, stack_path)

    def dump_allocations():
        import tracemalloc
        snapshot = tracemalloc.take_snapshot()
        tmp = alloc_path + ".tmp"
        with open(tmp, "w") as handle:
            handle.write(f"# wall_s={round(time.time() - started, 2)}\n")
            handle.write("# size_mb\tcount\tsource\n")
            for stat in snapshot.statistics("lineno")[:40]:
                frame = stat.traceback[0]
                handle.write(f"{stat.size / 1048576.0:.3f}\t{stat.count}\t"
                             f"{os.path.basename(frame.filename)}:{frame.lineno}\n")
        os.replace(tmp, alloc_path)

    def write_meta(final=False):
        meta = {"experiment": experiment, "job": job, "task": task,
                "pid": os.getpid(), "host": os.uname().nodename,
                "python": sys.version.split()[0],
                "slurm_mem": os.environ.get("SIMPLICITY_SLURM_MEM"),
                "slurm_time": os.environ.get("SIMPLICITY_SLURM_TIME"),
                "cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK"),
                "sample_interval_s": interval,
                "stack_interval_s": stack_interval if stacks_on else None,
                "sizes": sizes_on, "tracemalloc": tracemalloc_on,
                "wall_s": round(time.time() - started, 2),
                "cpu_s": cpu_s(),
                "peak_rss_mb": hwm_mb(),
                "hooked": sorted(state["hooked"]),
                "exited_cleanly": final}
        tmp = meta_path + ".tmp"
        with open(tmp, "w") as handle:
            json.dump(meta, handle, indent=1)
        os.replace(tmp, meta_path)

    main_thread_id = threading.get_ident()
    install_import_hook()
    hook_now()   # in case simplicity was already imported

    def loop():
        handle = writer = None
        pending = []
        last_memory = 0.0
        while True:
            try:
                time.sleep(stack_interval if stacks_on else interval)
                hook_now()
                if stacks_on:
                    sample_stack()
                now = time.time()
                if now - last_memory < interval:
                    continue
                last_memory = now
                row = sample_row()
                if writer is None:
                    if not state["seen"]:
                        # no simulation here yet: hold the startup curve in
                        # memory rather than writing a file this process may
                        # never earn. Bounded, so a long-lived driver process
                        # cannot grow it without limit.
                        pending.append(row)
                        del pending[:-240]
                        continue
                    handle = open(csv_path, "w", newline="")
                    writer = csv.DictWriter(handle, fieldnames=_COLUMNS,
                                            restval="")
                    writer.writeheader()
                    for buffered in pending:
                        writer.writerow(buffered)
                    pending = []
                writer.writerow(row)
                handle.flush()
                os.fsync(handle.fileno())   # survives a SIGKILL from Slurm
                if stacks_on:
                    dump_stacks()
                if tracemalloc_on:
                    dump_allocations()
                write_meta()
            except Exception:
                # a profiler must never take a run down with it
                try:
                    if handle is not None:
                        handle.close()
                except Exception:
                    pass
                return

    def at_exit():
        """Each step is isolated: a failure in one must not cost us the
        others. write_meta in particular decides whether the report calls this
        process cleanly exited or killed, so it can never ride on a dump
        succeeding -- tracemalloc can already be torn down by this point."""
        if not state["seen"]:
            return   # this process never ran a simulation
        for action in (_append_final_row, dump_stacks_if_on,
                       dump_allocations_if_on, lambda: write_meta(final=True)):
            try:
                action()
            except Exception:
                pass

    def _append_final_row():
        new = not os.path.exists(csv_path)
        with open(csv_path, "a", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=_COLUMNS, restval="")
            if new:
                writer.writeheader()   # a run too short for the sampler to open it
            writer.writerow(sample_row())

    def dump_stacks_if_on():
        if stacks_on:
            dump_stacks()

    def dump_allocations_if_on():
        if tracemalloc_on:
            dump_allocations()

    import atexit
    atexit.register(at_exit)
    threading.Thread(target=loop, name="simplicity-profiler", daemon=True).start()


if _DIR:
    _load_shadowed_sitecustomize()
    try:
        _start()
    except Exception:
        pass
