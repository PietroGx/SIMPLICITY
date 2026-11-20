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
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from urllib.parse import urlparse, parse_qs

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
WEB = os.path.join(HERE, "web")
sys.path.insert(0, HERE)

import compute
import dist
import fig1
import panels
import sweep

# Which experiment this server serves. The rebuild endpoints pass these to the
# export scripts, which otherwise default to experiment 1 and would write a
# different experiment's data into this store under identical directory names
# -- compute.params_key carries no experiment identity.
# The store is Data/figures_studio, untracked like the rest of Data/. It was
# ~/fig2_studio, then <repo>/fig2_studio, then Data/fig2_studio until v2.5.19.
BASE_STORE = os.path.join(ROOT, "Data", "figures_studio")
EXP_NUM = 1
EXP_NAME = compute.EXP_NAME
STORE = sweep.store_dir(BASE_STORE, EXP_NUM, EXP_NAME)
PY = sys.executable


def source_dir():
    """The run output this store was built from. /api/index carries it so the
    method page can name the experiment instead of a literal, which went stale
    the moment the studio moved to #2. Resolved per request: main() rebinds
    EXP_NAME/EXP_NUM after this module is imported."""
    from impact_long_shedders_config import scenario_experiment
    import simplicity.dir_manager as dm
    name = scenario_experiment(EXP_NAME, EXP_NUM)
    try:
        return dm.get_experiment_dir(name)
    except Exception:
        return name


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


# ------------------------------------------------------------------ handler

class H(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=WEB, **kw)

    def log_message(self, fmt, *args):
        pass

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
                ix["experiment"] = dict(name=EXP_NAME, num=EXP_NUM,
                                        source=source_dir())
                ix["axes"] = grid_axes()
                return self._send(ix)
            if u.path == "/api/jobs":
                return self._send(JOBS.status())
            if u.path == "/api/panel":
                return self._send(panels.panel(STORE, q))
            if u.path == "/api/dist":
                return self._send(dist.dist(STORE, q))
            if u.path == "/api/sens":
                return self._send(dist.sens(STORE, q))
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
            if u.path == "/api/fig1":
                return self._send(fig1.load(STORE))
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

        heavy = {"/api/rebuild/clades": "clades",
                 "/api/rebuild/transmissions": "tx"}
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
                    "--exp-num", str(EXP_NUM), "--exp-name", EXP_NAME,
                    "--out", os.path.join(STORE, "cache"),
                    "--procs", str(int(body.get("procs", 2))),
                    "--seeds", str(int(body.get("seeds", 50))),
                    "--window", p["window"],
                    "--k", str(p["k"]), "--mt", str(p["mt"]),
                    "--win", str(p["win"]),
                    "--prom", str(p["prom"]), "--lo", str(p["lo"]), "--hi", str(p["hi"])]
            name = f"clade build {compute.params_key(p)} [{EXP_NAME} #{EXP_NUM}]"
        else:
            argv = [PY, "-I", os.path.join(ROOT, "tests", "export_transmissions.py"),
                    "--exp-num", str(EXP_NUM), "--exp-name", EXP_NAME,
                    "--out", os.path.join(STORE, "tx"),
                    "--procs", str(int(body.get("procs", 2)))]
            name = f"transmissions [{EXP_NAME} #{EXP_NUM}]"
        return argv, name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8600)
    ap.add_argument("--exp-num", type=int, default=EXP_NUM)
    ap.add_argument("--exp-name", default=EXP_NAME)
    ap.add_argument("--store", default=None,
                    help="override; by default the experiment's own store")
    a = ap.parse_args()
    globals()["EXP_NUM"] = a.exp_num
    globals()["EXP_NAME"] = a.exp_name
    # panels.run_daily reads the RAW run output, so it needs the experiment too;
    # compute.ssod_for otherwise defaults to 1 and the Sweep tab falls back.
    panels.EXP_NUM = a.exp_num
    panels.EXP_NAME = a.exp_name
    globals()["STORE"] = (os.path.expanduser(a.store) if a.store
                          else sweep.store_dir(BASE_STORE, a.exp_num, a.exp_name))
    if a.host not in ("127.0.0.1", "localhost", "::1"):
        print(f"REFUSING to bind {a.host}: this serves compute endpoints and the "
              f"box runs production vhosts. Use an SSH tunnel.", file=sys.stderr)
        return 2
    srv = ThreadingHTTPServer((a.host, a.port), H)
    print(f"studio on http://{a.host}:{a.port}/   {EXP_NAME} #{EXP_NUM}   "
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
