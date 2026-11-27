"""Figure-2 studio server. Local only.

Binds 127.0.0.1 by default and is reached from a laptop over an SSH tunnel:

    ssh -N -L 8600:127.0.0.1:8600 <this box>
    then open http://127.0.0.1:8600/

Two classes of endpoint, because one class can occupy the machine for minutes
and this box also serves production vhosts:

  cheap   read caches, aggregate, answer in well under a second. Safe to expose
          behind a password if the app is ever published.
  heavy   run a compute pass. Hold a single-job lock, refuse concurrency, and
          are rejected unless the request arrives on the loopback interface.

Every path the browser uses is origin-relative, so putting nginx in front needs
no code change.
"""
import argparse, json, os, shutil, subprocess, sys, threading, time
import numpy as np
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WEB = os.path.join(HERE, "web")
sys.path.insert(0, HERE)

import compute
import dist
import explore
import fig1
import layout
import panels
import sweep

# Which experiment this server serves. The rebuild endpoints pass these to the
# export scripts, which otherwise default to experiment 1 and would write a
# different experiment's data into this store under identical directory names
# -- compute.params_key carries no experiment identity.
# The store is Data/figures, untracked like the rest of Data/. It was
# ~/fig2_studio, <repo>/fig2_studio, Data/fig2_studio, then Data/figures_studio
# until v2.5.21.
BASE_STORE = os.path.join(ROOT, "Data", "figures")
if os.path.isdir(os.path.join(ROOT, "Data", "figures_studio")):
    sys.exit("Data/figures_studio is the pre-v2.5.21 store: move its contents into Data/figures\n"
             "  mv Data/figures_studio/* Data/figures/ && rmdir Data/figures_studio")
if layout.legacy_stores(BASE_STORE):
    sys.exit("Data/figures holds the pre-v2.5.21 layout: python -I figures_gui/layout.py migrate --apply")
# compute.py, panels.py and export_transmissions.py keep the name --exp-num/EXP_NUM
# (their AST is their cache key) but take a minute stamp as well as a number.
# STAMP / EXP_NAME from the environment, as run_precompute.sh takes them;
# --stamp / --exp-name override
STAMP = os.environ.get("STAMP", "1")
EXP_NAME = os.environ.get("EXP_NAME", compute.EXP_NAME)
STORE = sweep.store_dir(BASE_STORE, STAMP, EXP_NAME)
PY = sys.executable


def source_dir():
    """The run output this store was built from. /api/index carries it so the
    method page can name the experiment instead of a literal, which went stale
    the moment the studio moved to #2. Resolved per request: main() rebinds
    EXP_NAME/STAMP after this module is imported."""
    from impact_long_shedders_config import scenario_experiment
    import simplicity.dir_manager as dm
    name = scenario_experiment(EXP_NAME, STAMP)
    try:
        return dm.get_experiment_dir(name)
    except Exception:
        return name


def sa_levels():
    """The sensitivity tabs' levels (sa_levels.json), the file warm_explore.py
    warms from, so the page asks for exactly what was precomputed."""
    with open(os.path.join(HERE, "sa_levels.json")) as fh:
        return {k: v for k, v in json.load(fh).items() if not k.startswith("_")}


def grid_axes():
    """The values the build grid was swept over, straight from sweep.py.

    A build parameter names a directory on disk, so a slider that can stop
    between grid values is a trap: prom=0.5 has no build, /api/run_daily
    answers 404, and the Sweep tab falls back to the 9.2-day export -- where
    m_F is never applied, because the fallback measures episodes in the
    browser. The rail asks for these and stops only on them.

    Served here rather than from panels.index so that adding it does not move
    panels.code_hash() and discard every memoised panel.
    """
    return dict(k=list(sweep.K), mt=[float(v) for v in sweep.MT],
                win=list(sweep.WIN), prom=list(sweep.PROM),
                lo=[sweep.LO], hi=[sweep.HI], window=[sweep.WINDOW])


# ------------------------------------------------------------------ jobs

class Jobs:
    """One compute pass at a time. A second request is refused, never queued:
    two passes at --procs 2 would put this box into swap, and the OOM killer
    does not necessarily choose the compute job."""

    def __init__(self):
        self.lock = threading.Lock()
        self.cur = None
        self.log = []
        self.history = []

    def status(self):
        with self.lock:
            return dict(running=self.cur, tail=self.log[-40:],
                        history=self.history[-10:])

    def start(self, name, argv):
        with self.lock:
            if self.cur:
                return False, f"{self.cur['name']} is already running"
            self.cur = dict(name=name, argv=argv, started=time.time())
            self.log = []
        threading.Thread(target=self._run, args=(name, argv), daemon=True).start()
        return True, "started"

    def _run(self, name, argv):
        rc = -1
        try:
            p = subprocess.Popen(argv, cwd=ROOT, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, bufsize=1)
            for line in p.stdout:
                with self.lock:
                    self.log.append(line.rstrip())
                    if len(self.log) > 400:
                        del self.log[:200]
            rc = p.wait()
        except Exception as exc:
            with self.lock:
                self.log.append(f"launch failed: {exc}")
        finally:
            with self.lock:
                started = self.cur["started"] if self.cur else time.time()
                self.history.append(dict(name=name, rc=rc,
                                         secs=round(time.time() - started, 1),
                                         when=time.strftime("%H:%M:%S")))
                self.cur = None


JOBS = Jobs()


def cache_cap(limit_gb):
    """Keep the cache under a cap, oldest keyed directory first. The box runs at
    84% full and a full disk would break the production vhosts, not just this."""
    root = os.path.join(STORE, "cache")
    if not os.path.isdir(root):
        return
    def size(p):
        return sum(os.path.getsize(os.path.join(d, f))
                   for d, _, fs in os.walk(p) for f in fs)
    ds = [(os.path.getmtime(os.path.join(root, d)), os.path.join(root, d))
          for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))]
    ds.sort()
    tot = sum(size(p) for _, p in ds)
    cap = limit_gb * 1024 ** 3
    while tot > cap and len(ds) > 1:
        _, p = ds.pop(0)
        tot -= size(p)
        shutil.rmtree(p, ignore_errors=True)


# ------------------------------------------------------------------ presets
# Named sets of the data rules, saved from the exploration tabs (presets.js).
# One file for every experiment: a rule set is not tied to a run.
PRESETS = os.path.join(BASE_STORE, "detection_presets.json")
PRESET_KEYS = {"k", "mt", "originRule", "sdKernel", "sdWin", "sdProm", "sdPromRef", "sdWidth", "sdCens", "sdStart",
               "sdTrough", "sdFrac", "sdPick", "sdGate", "burnin",
               "brkF", "brkH", "brkLo"}


def presets_read():
    try:
        with open(PRESETS) as fh:
            d = json.load(fh)
        return d if isinstance(d, dict) else {}
    except (OSError, ValueError):
        return {}


def presets_write(name, values):
    name = str(name or "").strip()
    if not name or len(name) > 60:
        raise ValueError("a preset needs a name of 1 to 60 characters")
    d = presets_read()
    if values is None:
        d.pop(name, None)
    else:
        if not isinstance(values, dict):
            raise ValueError("a preset's settings must be an object")
        bad = sorted(set(values) - PRESET_KEYS)
        if bad or not all(isinstance(v, (int, float, str, bool)) for v in values.values()):
            raise ValueError(f"not a detection setting: {bad or 'a non-scalar value'}")
        d[name] = values
    os.makedirs(os.path.dirname(PRESETS), exist_ok=True)
    tmp = f"{PRESETS}.{os.getpid()}.tmp"
    with open(tmp, "w") as fh:
        json.dump(d, fh, indent=1, sort_keys=True)
    os.replace(tmp, PRESETS)
    return d


def count_tests(groups):
    """Per-simulation counts by scenario -> Kruskal-Wallis across all groups
    and two-sided Mann-Whitney U for every pair, Holm-adjusted over the pairs.
    The simulation is the unit: one count per run."""
    from scipy.stats import kruskal, mannwhitneyu
    names = [n for n, v in groups.items() if isinstance(v, list) and len(v) >= 2]
    vals = {n: [float(x) for x in groups[n]] for n in names}
    out = dict(kruskal=None, pairs=[])
    if len(names) >= 2:
        try:
            h, p = kruskal(*(vals[n] for n in names))
        except ValueError:
            h, p = float("nan"), float("nan")
        if not (np.isfinite(h) and np.isfinite(p)):   # every value identical: nothing to rank
            h, p = 0.0, 1.0
        out["kruskal"] = dict(H=float(h), p=float(p), groups=len(names))
    pairs = []
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            if len(set(vals[a] + vals[b])) < 2:
                u, p = len(vals[a]) * len(vals[b]) / 2, 1.0
            else:
                u, p = mannwhitneyu(vals[a], vals[b], alternative="two-sided")
            pairs.append(dict(a=a, b=b, U=float(u), p=float(p)))
    order = sorted(range(len(pairs)), key=lambda k: pairs[k]["p"])
    m, run = len(pairs), 0.0
    for rank, k in enumerate(order):
        run = max(run, min(1.0, (m - rank) * pairs[k]["p"]))
        pairs[k]["p_holm"] = run
    out["pairs"] = pairs
    return out


# ------------------------------------------------------------------ handler

class H(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=WEB, **kw)

    def log_message(self, fmt, *args):
        pass

    def end_headers(self):
        # page files revalidate on every load (a 304 when unchanged), so a plain reload shows an edit
        if not self.path.startswith("/api/"):
            self.send_header("Cache-Control", "no-cache")
        super().end_headers()

    def _send(self, obj, code=200):
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _local(self):
        return self.client_address[0] in ("127.0.0.1", "::1")

    # ---- GET
    def do_GET(self):
        u = urlparse(self.path)
        q = {k: v[0] for k, v in parse_qs(u.query).items()}
        try:
            if u.path == "/api/index":
                ix = panels.index(STORE)
                ix["experiment"] = dict(name=EXP_NAME, num=STAMP,
                                        source=source_dir())
                ix["axes"] = grid_axes()
                ix["sa_levels"] = sa_levels()
                return self._send(ix)
            if u.path == "/api/jobs":
                return self._send(JOBS.status())
            if u.path == "/api/panel":
                return self._send(panels.panel(STORE, q))
            if u.path == "/api/dist":
                return self._send(dist.dist(STORE, q))
            if u.path == "/api/episodes":
                return self._send(panels.episode_view(STORE, q))
            if u.path == "/api/across":
                return self._send(panels.across(STORE, q))
            if u.path == "/api/runs":
                return self._send(panels.run_index(STORE))
            if u.path == "/api/run":
                return self._send(panels.run_file(STORE, q))
            if u.path == "/api/build_run":
                return self._send(panels.build_run(STORE, q))
            if u.path == "/api/run_daily":
                return self._send(panels.run_daily(STORE, q))
            if u.path == "/api/explore":
                return self._send(explore.explore(STORE, q))
            if u.path == "/api/explore_panel":
                return self._send(explore.explore_panel(STORE, q))
            if u.path == "/api/fig1":
                return self._send(fig1.load(layout.root_of(STORE)))
            if u.path == "/api/presets":
                return self._send(presets_read())
        except FileNotFoundError as exc:
            return self._send({"error": f"not computed yet: {exc}"}, 404)
        except Exception as exc:
            return self._send({"error": f"{type(exc).__name__}: {exc}"}, 500)
        if u.path == "/":
            self.path = "/index.html"
        return super().do_GET()

    # ---- POST
    def do_POST(self):
        u = urlparse(self.path)
        n = int(self.headers.get("Content-Length") or 0)
        try:
            body = json.loads(self.rfile.read(n) or b"{}")
        except Exception:
            return self._send({"error": "bad json"}, 400)

        if u.path == "/api/count_tests":
            try:
                return self._send({k: count_tests(v) for k, v in (body.get("metrics") or {}).items()})
            except Exception as exc:                                 # noqa: BLE001
                return self._send({"error": f"{type(exc).__name__}: {exc}"}, 400)
        if u.path == "/api/presets":
            if not self._local():
                return self._send({"error": "saving presets is loopback-only"}, 403)
            try:
                return self._send(presets_write(body.get("name"), body.get("settings")))
            except Exception as exc:                                 # noqa: BLE001
                return self._send({"error": f"{type(exc).__name__}: {exc}"}, 400)

        heavy = {"/api/rebuild/clades": "clades",
                 "/api/rebuild/transmissions": "tx",
                 "/api/rebuild/fig1": "fig1"}
        if u.path not in heavy:
            return self._send({"error": "no such endpoint"}, 404)
        if not self._local():
            return self._send({"error": "compute is loopback-only"}, 403)

        which = heavy[u.path]
        try:
            argv, name = self._job(which, body)
        except Exception as exc:
            # do_GET guards; this did not, so a bad parameter killed the handler
            # thread and the Export page saw a dropped connection, not an error
            return self._send({"error": f"{type(exc).__name__}: {exc}"}, 500)

        cache_cap(float(body.get("cap_gb", 4)))
        ok, msg = JOBS.start(name, argv)
        return self._send({"started": ok, "msg": msg, "name": name},
                          200 if ok else 409)

    def _job(self, which, body):
        if which == "clades":
            p = panels.clean_params(body)
            argv = [PY, "-I", os.path.join(HERE, "compute.py"),
                    "--exp-num", str(STAMP), "--exp-name", EXP_NAME,
                    "--out", os.path.join(STORE, "cache"),
                    "--procs", str(int(body.get("procs", 2))),
                    "--seeds", str(int(body.get("seeds", 50))),
                    "--window", p["window"],
                    "--k", str(p["k"]), "--mt", str(p["mt"]),
                    "--win", str(p["win"]),
                    "--prom", str(p["prom"]), "--lo", str(p["lo"]), "--hi", str(p["hi"])]
            name = f"clade build {compute.params_key(p)} [{EXP_NAME} #{STAMP}]"
        elif which == "fig1":
            # here, never on the cluster: it reads patient data (fig1.INPUTS)
            argv = [PY, "-I", os.path.join(HERE, "fig1.py"), "build", "--stamp", str(STAMP),
                    "--exp-name", EXP_NAME, "--store", layout.root_of(STORE)]
            name = f"Figure 1 [{EXP_NAME} #{STAMP}]"
        else:
            argv = [PY, "-I", os.path.join(ROOT, "figures_gui", "export_transmissions.py"),
                    "--exp-num", str(STAMP), "--exp-name", EXP_NAME,
                    "--out", os.path.join(STORE, "tx"),
                    "--procs", str(int(body.get("procs", 2)))]
            name = f"transmissions [{EXP_NAME} #{STAMP}]"
        return argv, name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8600)
    ap.add_argument("--stamp", default=STAMP, help="a stamp, an older number, or latest")
    ap.add_argument("--exp-name", default=EXP_NAME)
    ap.add_argument("--store", default=None,
                    help="override; by default the experiment's own store")
    a = ap.parse_args()
    from _scenarios import resolve_stamp
    a.stamp = resolve_stamp(a.stamp, a.exp_name)
    globals()["STAMP"] = a.stamp
    globals()["EXP_NAME"] = a.exp_name
    # panels.run_daily reads the RAW run output, so it needs the experiment too;
    # compute.ssod_for otherwise defaults to 1 and the Sweep tab falls back.
    panels.EXP_NUM = a.stamp
    panels.EXP_NAME = a.exp_name
    globals()["STORE"] = (os.path.expanduser(a.store) if a.store
                          else sweep.store_dir(BASE_STORE, a.stamp, a.exp_name))
    if a.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"REFUSING to bind {a.host}: this serves compute endpoints and the "
              f"box runs production vhosts. Use an SSH tunnel.", file=sys.stderr)
        return 2
    srv = ThreadingHTTPServer((a.host, a.port), H)
    print(f"studio on http://{a.host}:{a.port}/   {EXP_NAME} #{STAMP}   "
          f"store={STORE}", flush=True)
    for sub in ("cache", "runs", "tx"):
        if not os.path.isdir(os.path.join(STORE, sub)):
            print(f"  WARNING: {sub}/ is missing from this store", file=sys.stderr)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
