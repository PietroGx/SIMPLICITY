/* S1 (all trees), S2 (clustering sensitivity) and S3 (outbreak detection
   sensitivity) against a running studio.

     python -I figures_gui/server.py --stamp 2 --port 8619 &
     node figures_gui/test_sa_tabs.js [http://127.0.0.1:8619] [--skip-s2]

   The sensitivity tabs are only as good as their claim to be "the figure's own
   measurement at another setting", so the core checks hold them to that:
   saFig5 at no patch must give the per-simulation means of swScenCounts, the
   counting Cryptic by scenario and Figure 5's panel A share (check_app.js pins
   those two together); under a patch it must equal swScenCounts under the same
   patch; and a patch must never leak into the settings. The renders count draw
   calls per canvas, since an empty frame looks like a boring chart.
*/
const fs = require("fs"), path = require("path"), vm = require("vm"), http = require("http");

const realTimeout = setTimeout;
const wait = ms => new Promise(r => realTimeout(r, ms));
const ARGS = process.argv.slice(2), SKIP_S2 = ARGS.includes("--skip-s2");
const BASE = (ARGS.find(a => a.startsWith("http")) || "http://127.0.0.1:8619").replace(/\/$/, "");
const WEB = path.join(__dirname, "web");
const F = require(path.join(__dirname, "fakedom.js"));

F.El.prototype.append = function (...cs) { cs.forEach(c => this.children.push(c)); };
F.El.prototype.prepend = function (...cs) { cs.forEach(c => this.children.unshift(c)); };
if (!F.El.prototype.closest) F.El.prototype.closest = function () { return null; };
/* one context per canvas element, so a cell or facet that drew nothing shows */
const CTXS = [];
F.El.prototype.getContext = function () {
  if (!this._ctx) { this._ctx = F.mkCtx(); this._ctx.el = this; }
  /* a canvas the last render drew is reused by id: count it afresh */
  if (!CTXS.includes(this._ctx)) { this._ctx.ops = 0; CTXS.push(this._ctx); }
  return this._ctx;
};

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

let pass = 0, fail = 0;
const ck = (c, m) => { c ? pass++ : fail++; console.log((c ? "  ok   " : "  FAIL ") + m); };
const near = (a, b) => (a === null && b === null) || Math.abs(a - b) < 1e-9;

const SCRIPTS = fs.readFileSync(path.join(WEB, "index.html"), "utf8")
  .match(/<script src="([^"]+)"><\/script>/g)
  .map(s => s.replace(/.*src="|"><\/script>/g, ""));

async function main() {
  const up = await httpGet(BASE + "/api/index").catch(() => null);
  if (!up || !up.ok) { console.log(`server not answering at ${BASE}`); process.exit(2); }
  F.REG.clear();
  global.document.body = F.$el("body");
  const ctx = vm.createContext(global);
  for (const s of SCRIPTS)
    new vm.Script(fs.readFileSync(path.join(WEB, s), "utf8"), {filename: s}).runInContext(ctx);
  const peek = expr => new vm.Script(`(${expr})`).runInContext(ctx);
  await peek("loadAxes()");
  peek("startStudio()");

  /* ---- the levels are the file's, one copy */
  const file = JSON.parse(fs.readFileSync(path.join(__dirname, "sa_levels.json"), "utf8")).outbreak;
  ck(JSON.stringify(peek("SA_LEVELS.outbreak")) === JSON.stringify(file),
     "S3's levels in the page are sa_levels.json's, through /api/index");
  ck(peek(`Object.keys(SA_OS).every(k => (SA_LEVELS.outbreak || {})[k])`),
     "every rule S3 offers has levels");

  /* ---- Figure 2 and 4 extraction, on payloads shaped like the endpoints' */
  global.__A = {scenarios: {control: {seeds: [
    {seed: "a", per_value: [{m: 0.01, clades: 4, median_size: 3}, {m: 0.1, clades: 1, median_size: 9}]},
    {seed: "b", per_value: [{m: 0.01, clades: 6, median_size: 5}, {m: 0.1, clades: 3, median_size: 7}]}]}}};
  const f2 = peek("saFig2(__A, 0.1)");
  ck(near(f2.control.clades.mean, 2) && near(f2.control.median_size.mean, 8),
     `saFig2 reads the asked m_T only: clades ${f2.control.clades.mean}, size ${f2.control.median_size.mean} at m_T = 10%`);
  ck(JSON.stringify(f2.SOT) === "{}", "saFig2: a scenario with no seeds is empty, not zero");
  global.__H = {scenarios: {HIV_high: {burden: {any_long: {n: 50, mean: 0.25, lo: 0.2, hi: 0.3},
    only_long: {n: 50, mean: 0.05, lo: null, hi: null}}, infections: {any_long: {n: 50, mean: 0.4, lo: 0.35, hi: 0.45}}}}};
  const f4 = peek("saFig4(__H)");
  ck(near(f4.HIV_high.any_long.mean, 25) && near(f4.HIV_high.any_long.lo, 20) && f4.HIV_high.only_long.lo === null
     && near(f4.HIV_high.inf_long.mean, 40), "saFig4: fractions to percentages, a null bound kept null");

  /* ---- Figure 5, live: the figure's own counting */
  global.__D = await peek(`apiGet("/api/explore_panel?" + exploreQuery({mf: 0}))`);
  ck(!!(global.__D && global.__D.scenarios), "Figure 5's panel answers at the settings in force");
  const cmp = patch => peek(`(() => {
    const P = ${JSON.stringify(patch)};
    const C = swWithSettings(P, () => swScenCounts(__D)), S = saFig5(__D, P), m = a => a.reduce((t, x) => t + x, 0) / a.length;
    return SCENARIOS.map(sc => ({sc, n: C[sc].brk.length, brk: [m(C[sc].brk), (S[sc].brk || {}).mean],
                                 cry: [m(C[sc].cry), (S[sc].cry || {}).mean], nS: (S[sc].brk || {}).n}));
  })()`);
  for (const patch of [{}, {brkF: 0.2}, {brkH: 15}, {brkH: 120}, {brkF: 0.1, brkH: 30}]) {
    const rows = cmp(patch);
    const bad = rows.filter(r => !near(r.brk[0], r.brk[1]) || !near(r.cry[0], r.cry[1]) || r.n !== r.nS);
    ck(!bad.length, `saFig5 ${JSON.stringify(patch)} = swScenCounts per simulation: `
       + rows.map(r => `${r.sc} ${r.brk[1].toFixed(2)}/${r.cry[1].toFixed(2)}`).join(", ")
       + (bad.length ? `  -- off: ${bad.map(r => r.sc).join(", ")}` : ""));
  }
  const before = peek("JSON.stringify([settings.brkF, settings.brkH, settings.p3Source])");
  peek(`saFig5(__D, {brkF: 0.5, brkH: 365})`);
  ck(peek("JSON.stringify([settings.brkF, settings.brkH, settings.p3Source])") === before,
     "a patch leaves the settings as they were");
  const sh = peek(`(() => { const S = saFig5(__D, {});
    return SCENARIOS.map(sc => S[sc].share ? S[sc].share.mean : null); })()`);
  ck(sh.every(v => v === null || (v >= 0 && v <= 100)), `cryptic share is a percentage: ${sh.map(v => v === null ? "—" : v.toFixed(0)).join(", ")}`);

  /* ---- S3 renders, one rule at a time */
  for (const ax of ["brkF", "brkH", "brkLo", "sdWin", "sdProm"]) {
    CTXS.length = 0;
    peek(`setSetting("osAxis", ${JSON.stringify(ax)})`);
    peek(`staleTabs.add("saout")`);
    await peek(`renderTab(tabById("saout"))`);
    const html = (F.REG.get("saout_main") || {}).innerHTML || "";
    const drew = CTXS.filter(c => c.ops > 0 && /^sao5_/.test(c.el.id || "")).length;
    ck(drew === 6 && !/class="err"/.test(html) && !/undefined|NaN/.test(html),
       `S3 over ${ax}: ${drew} of 6 facets drew, ${/class="err"/.test(html) ? "with" : "no"} error box`);
  }

  /* ---- S1: one tree per simulation */
  CTXS.length = 0;
  peek(`setSetting("trScen", "HIV_high")`);
  peek(`showTab("trees")`);              // the walk stops when its tab is not active
  peek(`staleTabs.add("trees")`);
  await peek(`renderTab(tabById("trees"))`);
  for (let i = 0; i < 1200 && CTXS.filter(c => c.ops > 0 && !c.el.id).length < 50; i++) await wait(25);
  const th = (F.REG.get("trees_main") || {}).innerHTML || "";
  const cells = (th.match(/class="trcell"/g) || []).length;
  const drewT = CTXS.filter(c => c.ops > 0 && !c.el.id).length;   // the cells' canvases carry no id
  // fakedom hands out a fresh canvas per lookup, so each of the shell's renders adds 50
  ck(cells === 50 && drewT >= 50 && drewT % 50 === 0, `S1: ${cells} cells, every render drew 50 trees (${drewT / 50} render(s)) for HIV high`);

  /* ---- S2: whatever is on disk at k = 5, 7, 10 */
  if (SKIP_S2) { console.log(`\n${pass} passed, ${fail} failed (S2 skipped)`); process.exit(fail ? 1 : 0); }
  CTXS.length = 0;
  peek(`setSetting("saVary", "k")`);
  /* the Figure 5 walk stops when its tab is not the active one, so show it */
  peek(`showTab("saclust")`);
  peek(`staleTabs.add("saclust")`);
  await peek(`renderTab(tabById("saclust"))`);
  /* the shell's debounced redraw may supersede that render: wait for the walk */
  for (let i = 0; i < 2400 && ((F.REG.get("sac_status") || {}).innerHTML || "") !== ""; i++) await wait(25);
  const st = (F.REG.get("sac_status") || {}).innerHTML || "", ms = (F.REG.get("sac_miss") || {}).innerHTML || "";
  ck(st === "" && !/Figure 5/.test(ms), `S2 over k: Figure 5 measured at every level${ms ? " -- " + ms.replace(/<[^>]+>/g, "").replace(/\s+/g, " ").slice(0, 160) : ""}`);
  const sh3 = (F.REG.get("saclust_main") || {}).innerHTML || "";
  const drew3 = CTXS.filter(c => c.ops > 0 && /^sac[245]_/.test(c.el.id || "")).length;
  ck(drew3 === 11, `S2 over k: ${drew3} of 11 facets drew`);
  ck(!/undefined|NaN/.test(sh3), "S2: no undefined or NaN in the page");

  console.log(`\n${pass} passed, ${fail} failed`);
  process.exit(fail ? 1 : 0);
}
main().catch(e => { console.log("crashed: " + (e.stack || e)); process.exit(1); });
