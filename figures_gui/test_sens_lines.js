/* Panel 4's sensitivity facets, against a stubbed /api/sens.

     python -I figures_gui/server.py --port 8601 &
     node figures_gui/test_sens_lines.js [http://127.0.0.1:8601]

   /api/sens does not exist yet, so the chart path cannot be exercised by
   check_app.js: sensFor returns null and the tab shows the note instead. This
   loads the studio the same way check_app.js does, replaces sensFor with a
   payload shaped exactly like the contract in web/dist.js, renders the Panel 4
   tab once per axis choice, and counts the draw calls on a per-canvas context.

   Zero draw calls is the failure mode this is here for: linesOverLevels hands
   back a context either way, so a series mapped to the wrong shape draws an
   empty frame and looks like a chart that is merely boring.
*/
const fs = require("fs"), path = require("path"), vm = require("vm"), http = require("http");

const realTimeout = setTimeout;
const wait = ms => new Promise(r => realTimeout(r, ms));
const BASE = (process.argv[2] || "http://127.0.0.1:8601").replace(/\/$/, "");
const WEB = path.join(__dirname, "web");
const F = require(path.join(__dirname, "fakedom.js"));

F.El.prototype.append = function (...cs) { cs.forEach(c => this.children.push(c)); };
F.El.prototype.prepend = function (...cs) { cs.forEach(c => this.children.unshift(c)); };
if (!F.El.prototype.closest) F.El.prototype.closest = function () { return null; };

/* One context per canvas, so a facet that drew nothing is visible as such --
   fakedom hands every canvas the same shared context. The element itself is
   reused across renders (getElementById answers from one registry), so the
   counters are zeroed per render rather than the map rebuilt. */
const CTXS = new Map();
F.El.prototype.getContext = function () {
  if (!this._ctx) this._ctx = F.mkCtx();
  CTXS.set(this.id || "?", this._ctx);
  return this._ctx;
};
function zeroCtxs() {
  CTXS.forEach(c => { c.ops = 0; c.texts.length = 0; });
}

function httpGet(url) {
  return new Promise((res, rej) => {
    http.get(url, r => {
      let b = "";
      r.on("data", d => b += d);
      r.on("end", () => res({ok: r.statusCode < 400, status: r.statusCode, body: b}));
    }).on("error", rej);
  });
}
global.fetch = async u => {
  const url = String(u);
  const r = await httpGet(url.startsWith("http") ? url : BASE + url);
  return {ok: r.ok, status: r.status,
          async json() { return JSON.parse(r.body); },
          async text() { return r.body; }};
};
global.AbortController = class { constructor() { this.signal = {aborted: false}; }
                                 abort() { this.signal.aborted = true; } };
global.URLSearchParams = URLSearchParams;
global.window = global;
global.setInterval = () => 0;
global.Path2D = class { moveTo() {} lineTo() {} arc() {} closePath() {} };

let pass = 0, fail = 0;
const ck = (c, m) => { c ? pass++ : fail++; console.log((c ? "  ok   " : "  FAIL ") + m); };

const SCRIPTS = fs.readFileSync(path.join(WEB, "index.html"), "utf8")
  .match(/<script src="([^"]+)"><\/script>/g)
  .map(s => s.replace(/.*src="|"><\/script>/g, ""));

const METRICS = ["pct_burden_majority_long", "rise_days", "latent_days"];
const SCS = ["control", "SOT", "HIV_low", "HIV_high", "edge_case"];

/* The contract in web/dist.js:

     {kind: "sens", axis, pinned: {...}, levels: [...],
      scenarios: {scenario: {metric: [{n, mean, lo, hi} per level, ...]}}}

   Three deliberate holes, because the real grid has them: a level with no
   build is null, a level under three runs has lo/hi null (panels.py:146-156),
   and one scenario is missing a metric entirely. */
function fakeSens(axis, levels, pinned) {
  const base = {pct_burden_majority_long: 18, rise_days: 22, latent_days: 40};
  const sc = {};
  SCS.forEach((s, si) => {
    const per = {};
    METRICS.forEach(m => {
      if (s === "edge_case" && m === "latent_days") return;     // no metric
      per[m] = levels.map((lv, i) => {
        if (s === "HIV_low" && i === 1) return null;            // no build
        const mean = base[m] * (1 + 0.08 * si) + i * (m === "rise_days" ? 1.4 : 0.9);
        if (s === "SOT" && i === 0) return {n: 2, mean, lo: null, hi: null};
        return {n: 50, mean, lo: mean - 1.3, hi: mean + 1.3};
      });
    });
    sc[s] = per;
  });
  return {kind: "sens", axis, pinned, levels, scenarios: sc};
}

async function main() {
  const up = await httpGet(BASE + "/api/index").catch(() => null);
  if (!up || !up.ok) {
    console.log(`server not answering at ${BASE}`);
    process.exit(2);
  }
  F.REG.clear();
  global.document.body = F.$el("body");
  const ctx = vm.createContext(global);
  for (const s of SCRIPTS)
    new vm.Script(fs.readFileSync(path.join(WEB, s), "utf8"), {filename: s}).runInContext(ctx);
  const peek = expr => new vm.Script(`(${expr})`).runInContext(ctx);
  await peek("loadAxes()");
  const AXES = peek("AXES");
  peek("startStudio()");

  const PINS = {
    k:    {mt: 0.01, win: 15, prom: 0.1},
    mt:   {k: 7, win: 15, prom: 0.1},
    win:  {k: 7, mt: 0.01, prom: 0.1},
    prom: {k: 7, mt: 0.01, win: 15},
  };

  for (const axis of ["k", "mt", "win", "prom"]) {
    const levels = AXES[axis];
    global.sensFor = async ax => fakeSens(ax, levels, PINS[ax]);
    zeroCtxs();
    const moved = peek(`setSetting("axis4", ${JSON.stringify(axis)})`);
    /* The axis selector is a view parameter: the shell has to mark the tab
       stale on its own, or the picker would change nothing until a build
       parameter moved. The render below is then forced, so the debounce in
       shell.js cannot race it. */
    ck(!moved || peek(`staleTabs.has("panel4")`),
       `axis=${axis}: changing settings.axis4 marked the tab stale`);
    peek(`staleTabs.add("panel4")`);
    await peek(`renderTab(tabById("panel4"))`);
    for (let i = 0; i < 40 && !CTXS.has("p4s3"); i++) await wait(25);
    const main4 = F.REG.get("panel4_main");
    const html = main4 ? main4.innerHTML : "";
    const ops = [1, 2, 3].map(i => (CTXS.get("p4s" + i) || {}).ops || 0);
    const texts = [1, 2, 3].map(i => ((CTXS.get("p4s" + i) || {}).texts || []));
    console.log(`\naxis=${axis}  levels=${JSON.stringify(levels)}`);
    console.log(`  draw calls per facet: ${ops.join(" / ")}`);
    ck(ops.every(n => n > 0),
       `axis=${axis}: all three facets drew (${ops.join("/")} ops)`);
    ck(/Sensitivity over/.test(html) && !/not answering yet/.test(html),
       `axis=${axis}: the chart block rendered, not the missing-endpoint note`);
    /* The levels have to reach the x axis as the grid's own stops, and the
       operating point has to be one of them -- an unmatched `mark` is dropped
       silently by linesOverLevels (charts.js:205-212). */
    const want = levels.map(v => peek(`p4Level(${JSON.stringify(axis)}, ${v})`));
    const seen = texts[2];
    /* xLabels thins the row when the stops crowd (charts.js:57-68), so not
       every stop is required -- but the first one is, and a thinned row is not
       an empty one. */
    const got = want.filter(t => seen.includes(t));
    ck(seen.includes(want[0]) && got.length >= Math.min(2, want.length),
       `axis=${axis}: x labels are the grid's own stops, ${got.length}/${want.length} shown `
       + `(${got.join(", ")})`);
    const op = peek(`p4Level(${JSON.stringify(axis)}, SETTING_DEFAULTS.${axis})`);
    ck(levels.map(v => peek(`p4Level(${JSON.stringify(axis)}, ${v})`)).includes(op),
       `axis=${axis}: the operating point ${op} is one of the levels, so the rule is drawn`);
    ck(!/>\s*(undefined|NaN)\s*</i.test(html) && !/(undefined|NaN)\s*(%|d\b)/i.test(html),
       `axis=${axis}: no undefined/NaN in the markup`);
  }

  /* And the other half of the contract: with no endpoint, the note -- never
     by_axis drawn against a duplicated axis. */
  global.sensFor = async () => null;
  CTXS.clear();
  peek(`setSetting("axis4", "k")`);
  peek(`staleTabs.add("panel4")`);
  await peek(`renderTab(tabById("panel4"))`);
  const html = (F.REG.get("panel4_main") || {}).innerHTML || "";
  ck(/not answering yet/.test(html) && !/id="p4s1"/.test(html),
     "sensFor null: the note is shown and no facet canvas is emitted");
  ck([1, 2, 3].every(i => !CTXS.has("p4s" + i)),
     "sensFor null: nothing was drawn");

  /* The band is the caption's claim, so it has to cost something: the same
     series through linesOverLevels with and without `ribbon`. */
  const one = ribbon => {
    const cv = peek("document.createElement('canvas')");
    cv.id = "ribbon_probe_" + (ribbon ? "on" : "off");
    const s = fakeSens("k", AXES.k, PINS.k);
    global.__probe = {cv, ribbon, levels: AXES.k.map(String), mark: "7",
      series: SCS.map(x => ({colour: "#888", label: x,
        vals: (s.scenarios[x] || {}).rise_days || []}))};
    peek(`linesOverLevels(__probe.cv, {levels: __probe.levels,
      series: __probe.series, mark: __probe.mark, ribbon: __probe.ribbon,
      ylab: "days", height: 190})`);
    return cv._ctx.ops;
  };
  const onOps = one(true), offOps = one(false);
  console.log(`\n  ribbon on/off: ${onOps} / ${offOps} draw calls`);
  ck(offOps > 0 && onOps > offOps,
     `the interval band is drawn: ${onOps} ops with it, ${offOps} without`);

  console.log(`\n${pass} passed, ${fail} failed`);
  process.exit(fail ? 1 : 0);
}

main();
