/* Pins the browser's episode detection to compute.py's.

     python -I figures_gui/server.py &
     node figures_gui/check_episodes.js [http://127.0.0.1:8600]

   web/episodes.js and figures_gui/compute.py both implement one definition --
   smooth, find the peaks clearing the prominence bar, measure each rise across
   the band. The browser's copy exists so the sliders answer without a
   twelve-minute rebuild, which means a drift between them would show the
   panels one number and the Sweep time tab another, with nothing to say which
   is right.

   THE 9.2-DAY GAP IS CLOSED (2026-11-18). The browser used to measure on
   runs/, the clade-explorer export, which samples 120 time points about 9.2
   days apart: control seed 0 at k=5, win=1, prom=0.1 gave python 275 emergence
   episodes against the browser's 13, with surviving rise days differing by up
   to 64 days. It now measures on the daily record served by /api/run_daily,
   the same axis compute.py reads, and the browser's episode pick was changed
   to compute.py's rule (the episode containing the clade's own maximum) from
   the artifact's largest-amplitude one.

   It also compares the two DESIGNATIONS, because the browser does its own.
   `rundata.computeClades` takes the export's divergence-only `starts[k]` and
   filters them by subtree transmission support; `compute.run_context`
   designates under rule B, where refusing a split changes the reference root
   for every descendant. Those are not obviously the same: a lineage could
   clear k from the FOLDED root without ever having been a divergence-only
   start, and the browser would miss it. Measured here, 0 mismatches -- which
   is evidence, not a proof, so this section exists to notice the day it stops
   being true.

   So this compares CLADE BY CLADE, not as a set. The join is each clade's root
   lineage index: the build names clades ("Clade_12") while the browser
   identifies them by root lineage, and /api/run_daily carries both. Two
   implementations of one definition agree per clade; a difference in the
   smoothing, the prominence test, the band or the pick moves one of them.
*/
const fs = require("fs"), path = require("path"), vm = require("vm"), http = require("http");

const BASE = (process.argv[2] || "http://127.0.0.1:8600").replace(/\/$/, "");
const WEB = path.join(__dirname, "web");
const TOL_DAYS = 0.02;            // compute.py rounds its day values to 2 dp

function httpGet(url) {
  return new Promise((res, rej) => {
    http.get(url, r => {
      let b = "";
      r.on("data", d => b += d);
      r.on("end", () => res({ok: r.statusCode < 400, status: r.statusCode, body: b}));
    }).on("error", rej);
  });
}

let pass = 0, fail = 0;
const ck = (c, m) => { c ? pass++ : fail++; console.log((c ? "  ok   " : "  FAIL ") + m); };

/* Only the modules the detection needs, in a bare context: no DOM, because
   nothing here draws. */
const NEED = ["settings.js", "api.js", "draw.js", "charts.js", "rundata.js", "paint.js", "episodes.js",
              "tab_panel3.js"];

function makeContext() {
  const g = {
    console, atob: s => Buffer.from(s, "base64").toString("binary"),
    Int32Array, Uint8Array, Uint16Array, Uint32Array, Float32Array, Float64Array,
    Math, Map, Set, JSON, Number, String, Array, Object, isFinite, isNaN,
    Infinity, NaN, parseInt, parseFloat, Error, Promise, URLSearchParams,
    setTimeout, clearTimeout, AbortController: class {
      constructor() { this.signal = {aborted: false}; } abort() {} },
    localStorage: {getItem: () => null, setItem: () => {}, removeItem: () => {}},
    getComputedStyle: () => ({getPropertyValue: () => "#888"}),
    document: {documentElement: {}, createElement: () => ({
      className: "", innerHTML: "", append: () => {}, querySelector: () => null,
      querySelectorAll: () => [], addEventListener: () => {}})},
    registerTab: () => {},
    fetch: async u => {
      const r = await httpGet(String(u).startsWith("http") ? String(u) : BASE + u);
      return {ok: r.ok, status: r.status, async json() { return JSON.parse(r.body); }};
    },
  };
  g.window = g;
  g.globalThis = g;
  const ctx = vm.createContext(g);
  for (const f of NEED)
    new vm.Script(fs.readFileSync(path.join(WEB, f), "utf8"), {filename: f}).runInContext(ctx);
  return ctx;
}

const sorted = a => a.slice().sort((x, y) => x - y);

async function main() {
  const up = await httpGet(BASE + "/api/index").catch(() => null);
  if (!up || !up.ok) {
    console.log(`server not answering at ${BASE} -- start figures_gui/server.py first`);
    process.exit(2);
  }
  const index = JSON.parse(up.body);
  const builds = (index.builds || []).filter(b => b.code === index.code);
  if (!builds.length) {
    console.log("no clade build matches the current compute.py -- nothing to compare against");
    process.exit(2);
  }
  const build = builds.find(b => b.params.k === 7) || builds[0];
  const p = build.params;
  console.log(`comparing against ${build.key}`);
  console.log(`  k=${p.k} win=${p.win} prom=${p.prom} band ${p.lo}-${p.hi} ${p.window}`);

  const ctx = makeContext();
  const run = expr => new vm.Script(`(${expr})`).runInContext(ctx);

  // the browser must be told the build's own parameters, or this compares two
  // different questions rather than two answers to one
  const S = {k: p.k, mt: p.mt, win: p.win, prom: p.prom, lo: p.lo, hi: p.hi,
             window: p.window, burnin: 1};
  run(`setSettings(${JSON.stringify(S)})`);

  let checked = 0;
  for (const scen of ["control", "HIV_high"]) {
    for (const seedIdx of [0, 1]) {
      const qs = new URLSearchParams(Object.assign({}, S, {
        scenario: scen, seed: seedIdx})).toString();
      const got = await httpGet(`${BASE}/api/build_run?${qs}`);
      if (!got.ok) {
        console.log(`  (the build has no ${scen} seed ${seedIdx}, skipped)`);
        continue;
      }
      const pyRun = JSON.parse(got.body).run;

      const dgot = await httpGet(`${BASE}/api/run_daily?${qs}`);
      if (!dgot.ok) {
        ck(false, `${scen} seed ${seedIdx}: /api/run_daily answered ${dgot.status}`);
        continue;
      }
      const daily = JSON.parse(dgot.body);
      const label = `${scen} seed ${seedIdx}`;

      // the join: clade name -> root lineage index, which is what the browser keys on
      const rootOf = new Map(daily.clades.map(c => [c.name, c.root]));

      /* the designations themselves, before any episode is measured */
      const buildRoots = new Set(daily.clades.filter(c => !c.founding)
                                             .map(c => c.root));
      const browserRoots = new Set(await run(
        `runView("${scen}", ${seedIdx}, ${p.k}, ${S.mt}).then(v =>
           [...v.C.alive].filter(r => !v.C.info.get(r).isRoot))`));
      const onlyBuild = [...buildRoots].filter(r => !browserRoots.has(r));
      const onlyBrowser = [...browserRoots].filter(r => !buildRoots.has(r));
      ck(!onlyBuild.length && !onlyBrowser.length,
         `${label}: one designation -- build ${buildRoots.size} non-root `
         + `clades, browser ${browserRoots.size}`
         + (onlyBuild.length ? `, build-only ${onlyBuild.slice(0, 5)}` : "")
         + (onlyBrowser.length ? `, browser-only ${onlyBrowser.slice(0, 5)}` : ""));
      ck(daily.t.length > 300,
         `${label}: the daily record is daily (${daily.t.length} points, `
         + `spacing ${(daily.t[1] - daily.t[0]).toFixed(2)} d)`);

      /* The browser's pool is non-root clades carrying a defining mutation, so
         match it: the founding clade and definer-less groups are not in it. */
      const pyClades = pyRun.clades.filter(c => !c.founding && c.ndef > 0);
      const pyEmgOf = new Map();        // clade index -> emergence rise days
      const pyEpsOf = new Map();        // clade index -> episode count
      for (const e of pyRun.episodes) {
        pyEpsOf.set(e.c, (pyEpsOf.get(e.c) || 0) + 1);
        if (e.emg) pyEmgOf.set(e.c, e.days);
      }

      const jsRows = await run(
        `runView("${scen}", ${seedIdx}, ${p.k}, ${S.mt}).then(async v => {
           const d = await apiGet("/api/run_daily?${qs}");
           const prep = prepareEpisodes(v, d);
           return {daily: prep.daily, unserved: prep.unserved,
                   rows: prep.rows.map(r => ({id: r.id, neps: r.eps.length,
                                              days: r.main ? r.main.days : null}))};
         })`);

      ck(jsRows.daily === true,
         `${label}: the browser measured on the daily record, not the export`);

      const jsOf = new Map(jsRows.rows.map(r => [r.id, r]));
      let bothHave = 0, worst = 0, worstName = "", nepsBad = [];
      for (const c of pyClades) {
        const root = rootOf.get(c.name);
        if (root === undefined) continue;          // not served at this m_T
        const js = jsOf.get(root);
        const py = pyEmgOf.get(c.c);
        if (!js || py === undefined || js.days === null) continue;
        bothHave++;
        const d = Math.abs(js.days - py);
        if (d > worst) { worst = d; worstName = c.name; }
        const pyn = pyEpsOf.get(c.c) || 0;
        if (js.neps !== pyn) nepsBad.push(`${c.name} (browser ${js.neps}, python ${pyn})`);
      }

      ck(bothHave > 0, `${label}: ${bothHave} clade(s) measured by both`);
      ck(worst <= TOL_DAYS,
         `${label}: every clade's emergence rise agrees `
         + `(worst ${worst.toFixed(3)} d${worstName ? " on " + worstName : ""} `
         + `over ${bothHave} clades)`);
      ck(!nepsBad.length,
         `${label}: episode counts agree per clade`
         + (nepsBad.length ? ` -- differ on ${nepsBad.slice(0, 4).join(", ")}` : ""));

      // and the whole multiset, which catches a clade one side dropped
      const pyRise = sorted([...pyEmgOf.entries()]
        .filter(([c]) => pyClades.some(x => x.c === c && rootOf.has(x.name)))
        .map(([, d]) => d));
      const jsRise = sorted(jsRows.rows.filter(r => r.days !== null).map(r => r.days));
      ck(Math.abs(jsRise.length - pyRise.length) <= 1,
         `${label}: emergence-episode count agrees `
         + `(python ${pyRise.length}, browser ${jsRise.length})`);
      checked++;
    }
  }
  ck(checked > 0, `${checked} run(s) compared`);
  await checkBreakouts(run, S);
  console.log(`\n${pass} passed, ${fail} failed`);
  process.exit(fail ? 1 : 0);
}

/* Figure 5 against Sweep time. Figure 5 times each clade's breakout from the
   build's episodes (/api/panel?n=3b, p3Breakouts); Sweep time from
   /api/run_daily's (prepareEpisodes). Both call breakoutOf, so they agree only
   while both are handed the same episodes: the same window and the same m_F.
   Sweep time must draw no breakout at all when its window is not the build's --
   the burn-in toggle off, or a simulation that never reaches the burn-in day,
   whose post-burn-in window is empty in the build. */
async function checkBreakouts(run, S) {
  console.log("\nbreakouts: Figure 5 against Sweep time");
  const qb = new URLSearchParams(S).toString();
  const g3 = await httpGet(`${BASE}/api/panel?n=3b&${qb}`);
  if (!g3.ok) { ck(false, `/api/panel?n=3b answered ${g3.status}`); return; }
  const d3 = JSON.parse(g3.body);

  const sweepOf = async (scen, i, patch) => {
    run(`setSettings(${JSON.stringify(Object.assign({}, S, patch))})`);
    const qs = new URLSearchParams(Object.assign({}, S, patch, {scenario: scen, seed: i})).toString();
    return run(`runView("${scen}", ${i}, ${S.k}, ${S.mt}).then(async v => {
        const d = await apiGet("/api/run_daily?${qs}");
        const prep = prepareEpisodes(v, d);
        const nameOf = new Map(d.clades.map(c => [c.root, c.name]));
        return {asFig5: prep.asFig5, off: prep.fig5Off, rows: prep.rows.map(r =>
          ({name: nameOf.get(r.id), brk: r.brk && {start: r.brk.start, peak: r.brk.peak}}))};
      })`);
  };
  /* the build's episodes are Figure 5's "Clade build" source, which applies m_F;
     its "Detection tabs" source follows Cryptic by scenario, which does not */
  const fig5Of = patch => {
    run(`setSettings(${JSON.stringify(Object.assign({}, S, patch, {p3Source: "build"}))})`);
    return run(`p3Breakouts(${JSON.stringify(d3)})`);
  };

  for (const [what, patch] of [["burn-in on, m_F 0", {burnin: 1, mf: 0, brkF: 0.05}],
                               ["m_F 10% above F 5%", {burnin: 1, mf: 0.1, brkF: 0.05}]]) {
    const B = fig5Of(patch);
    let n = 0, nBrk = 0;
    const bad = [];
    for (const scen of ["SOT", "HIV_high"]) for (const i of [0, 1, 2]) {
      const sw = await sweepOf(scen, i, patch);
      if (!sw.asFig5) continue;
      const swOf = new Map(sw.rows.map(r => [r.name, r]));
      const f5 = B[scen].cl.filter(c => c.run === i);
      const f5Names = new Set(f5.map(c => c.name));
      for (const c of f5) {
        const r = swOf.get(c.name);
        if (!r) continue;                         // past the daily cap
        n++;
        const same = c.ev ? !!r.brk && Math.abs(r.brk.start - c.a) <= TOL_DAYS
                              && Math.abs(r.brk.peak - c.top) <= 2e-4 : !r.brk;
        if (c.ev) nBrk++;
        if (!same) bad.push(`${scen}/${i} ${c.name} (fig5 ${c.ev ? c.a : "none"}, sweep ${r.brk ? r.brk.start : "none"})`);
      }
      sw.rows.filter(r => r.brk && !f5Names.has(r.name))
        .forEach(r => bad.push(`${scen}/${i} ${r.name} only on Sweep time`));
    }
    ck(n > 0 && nBrk > 0 && !bad.length,
       `${what}: ${n} clades, ${nBrk} breakouts, the same start and peak on both`
       + (bad.length ? ` -- ${bad.length} differ: ${bad.slice(0, 3).join("; ")}` : ""));
  }

  const off = await sweepOf("SOT", 0, {burnin: 0, mf: 0});
  ck(off.off === "burnin" && !off.rows.some(r => r.brk),
     `burn-in off: Sweep time draws no breakout (fig5Off ${off.off}, `
     + `${off.rows.filter(r => r.brk).length} drawn)`);

  let found = null;
  for (const [scen, v] of Object.entries(d3.scenarios))
    v.runs.forEach((r, i) => { if (!found && (r.phi === null || r.phi === undefined)) found = [scen, i, r]; });
  if (!found) { console.log("  (no simulation without a burn-in day in this build)"); return; }
  const [scen, i, r3] = found;
  const np = await sweepOf(scen, i, {burnin: 1, mf: 0});
  ck(r3.clades.every(c => !c[3].length),
     `${scen}/${i} never reaches the burn-in day: the build holds no episode for it`);
  ck(np.off === "nophi" && !np.rows.some(r => r.brk),
     `${scen}/${i}: Sweep time draws no breakout (fig5Off ${np.off}, `
     + `${np.rows.filter(r => r.brk).length} drawn)`);
}

main().catch(e => { console.log("FAIL " + e.message); process.exit(1); });
