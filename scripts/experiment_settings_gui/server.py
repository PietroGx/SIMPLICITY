#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''Experiment settings panel: pick every parameter of an impact_long_shedders
pipeline run and write it as a settings file into scripts/experiments_settings/.

    python -I scripts/experiment_settings_gui/server.py          # 127.0.0.1:8601
    ssh -N -L 8601:127.0.0.1:8601 <this box>                     # from a laptop

Validation, derived values and the file format are the pipeline's own
(experiments_pipeline/experiment_settings.py and impact_long_shedders_config),
so the panel cannot accept a file the pipeline would refuse. Loopback only.
'''
import argparse
import json
import os
import re
import sys
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, os.pardir, os.pardir))
sys.path[:0] = [REPO, os.path.join(REPO, "scripts", "experiments_pipeline")]

import numpy as np                                               # noqa: E402

import experiment_settings as es                                 # noqa: E402
import impact_long_shedders_config as cfg                        # noqa: E402
import simplicity.settings_manager as sm                         # noqa: E402

WEB = os.path.join(HERE, "web")
FILENAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*\.json$")
LAUNCHER = "scripts/launchers/run_pipeline_and_grid.sh"


def rel(path):
    return os.path.relpath(path, REPO)


def list_files():
    out = []
    for f in sorted(os.listdir(es.SETTINGS_DIR)):
        if not f.endswith(".json"):
            continue
        path = os.path.join(es.SETTINGS_DIR, f)
        row = {"name": f, "mtime": os.path.getmtime(path)}
        try:
            s = es.load(path)
            row.update(description=s["description"], valid=True)
        except Exception as exc:                                 # noqa: BLE001
            row.update(valid=False, error=str(exc).splitlines()[0])
        out.append(row)
    return out


def derive(s):
    """What the pipeline would build from `s`, through the config's own
    functions. Assumes validate(s) returned no errors."""
    sp = sm.read_standard_parameters_values()
    pop, scen = s["population"], s["scenarios"]
    n = pop["population_size"]
    warn = []
    rows = []
    for x in scen:
        tau = cfg.derive_tau_3_long(x, sp)
        rows.append(dict(name=x["name"], tau_3_long=round(tau, 3),
                         initial_long=int(round(n * x["long_shedders_ratio"]))))
        if x["long_shedders_ratio"] > 0 and tau <= 0:
            warn.append(f"{x['name']}: inf_duration_long {x['inf_duration_long']} d leaves "
                        f"tau_3_long = {tau:.2f} d after the other phases "
                        f"({cfg.phase_offset(sp):.2f} d); it must be longer")
        if x["long_shedders_ratio"] > 0 and rows[-1]["initial_long"] == 0:
            warn.append(f"{x['name']}: {n} x {x['long_shedders_ratio']} rounds to 0 long "
                        f"shedders in the initial population")
    c1 = s["cal_1"]
    cohorts = []
    if all(r["tau_3_long"] > 0 for r, x in zip(rows, scen) if x["long_shedders_ratio"] > 0):
        for pair in cfg.unique_long_tau_r_long_pairs(sp, scen):
            cohorts.append(dict(
                group=cfg.calib_long_group(cfg.group_name_for_pair(pair, sp, scen)),
                tau_3_long=pair[0], R_long=pair[1],
                final_time=cfg.cal1_final_time(pair[0], sp, c1["final_time"],
                                               c1["min_cycles"])))
    grid = lambda r: [float(f"{v:.4g}") for v in np.geomspace(r["min"], r["max"], r["steps"])]
    n1, n2 = c1["nsr_long"]["steps"], s["cal_2"]["nsr_std"]["steps"]
    sims = dict(cal_1=len(cohorts) * n1 * s["seeds"]["cal"],
                cal_2=n2 * s["seeds"]["cal"],
                exp=len(scen) * s["seeds"]["exp"])
    sims["total"] = sum(sims.values())

    ref = es.load(es.DEFAULT_SETTINGS)
    scale = ("population_size", "infected_individuals_at_start", "R", "final_time")
    if any(pop[k] != ref["population"][k] for k in scale) and (
            c1["nsr_long"] == ref["cal_1"]["nsr_long"]
            or s["cal_2"]["nsr_std"] == ref["cal_2"]["nsr_std"]):
        warn.append("The NSR ranges are the ones hand-calibrated for "
                    + ", ".join(f"{k}={ref['population'][k]}" for k in scale)
                    + "; with a different population they may not bracket the targets. "
                      "Re-choose them (tests/run_nsr_range_finder.sh) before trusting a run.")
    for x in scen:
        if x["R_long"] is not None and x["R_long"] != pop["R"]:
            warn.append(f"{x['name']}: R_long {x['R_long']} differs from population R {pop['R']}")
    return dict(reservoir=es.reservoir(n),
                infected_at_start_pct=100 * pop["infected_individuals_at_start"] / n,
                phase_offset=cfg.phase_offset(sp), scenarios=rows, cohorts=cohorts,
                cal2_final_time=pop["final_time"],
                nsr_long_grid=grid(c1["nsr_long"]), nsr_std_grid=grid(s["cal_2"]["nsr_std"]),
                simulations=sims, warnings=warn)


def check(s):
    errs = es.validate(s)
    if errs:
        return dict(errors=errs)
    try:
        return dict(errors=[], derived=derive(s))
    except Exception as exc:                                     # noqa: BLE001
        return dict(errors=[f"derived values failed: {type(exc).__name__}: {exc}"])


def launch_command(name):
    return f"bash {LAUNCHER} --settings {rel(os.path.join(es.SETTINGS_DIR, name))}"


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

    def end_headers(self):
        if not self.path.startswith("/api/"):
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_GET(self):
        path, _, query = self.path.partition("?")
        if path == "/api/files":
            return self._send(dict(dir=rel(es.SETTINGS_DIR), files=list_files(),
                                   default=os.path.basename(es.DEFAULT_SETTINGS)))
        if path == "/api/file":
            name = dict(p.split("=", 1) for p in query.split("&") if "=" in p).get("name", "")
            if not FILENAME.match(name):
                return self._send(dict(error="bad file name"), 400)
            target = os.path.join(es.SETTINGS_DIR, name)
            if not os.path.isfile(target):
                return self._send(dict(error=f"no such file: {name}"), 404)
            with open(target) as fh:
                try:
                    return self._send(dict(name=name, settings=json.load(fh)))
                except json.JSONDecodeError as exc:
                    return self._send(dict(error=f"{name}: not valid JSON ({exc})"), 422)
        return super().do_GET()

    def do_POST(self):
        try:
            body = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))) or b"{}")
        except json.JSONDecodeError:
            return self._send(dict(error="request body is not JSON"), 400)
        if self.path == "/api/check":
            return self._send(check(body.get("settings")))
        if self.path == "/api/save":
            name, s = body.get("name", ""), body.get("settings")
            if not FILENAME.match(name or ""):
                return self._send(dict(error="file name: letters, digits, _ . + - and "
                                             "ending in .json"), 400)
            res = check(s)
            if res["errors"]:
                return self._send(dict(error="not saved: the settings have errors",
                                       errors=res["errors"]), 422)
            target = os.path.join(es.SETTINGS_DIR, name)
            if os.path.exists(target) and not body.get("overwrite"):
                return self._send(dict(error=f"{name} exists; tick overwrite to replace it",
                                       exists=True), 409)
            es.write(s, target)
            return self._send(dict(saved=rel(target), command=launch_command(name)))
        return self._send(dict(error="unknown endpoint"), 404)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8601)
    a = ap.parse_args()
    if a.host not in ("127.0.0.1", "localhost", "::1"):
        raise SystemExit("loopback only: it writes files into the repo. "
                         "Reach it over an SSH tunnel.")
    srv = ThreadingHTTPServer((a.host, a.port), H)
    print(f"experiment settings panel: http://{a.host}:{a.port}/  "
          f"(writes to {rel(es.SETTINGS_DIR)}/)", flush=True)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
