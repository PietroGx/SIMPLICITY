/* Renders every tab of the studio headlessly against a RUNNING server.

     python -I figures_gui/server.py &
     node figures_gui/check_app.js [http://127.0.0.1:8600]

   Replaces check_pages.js, which rendered nine separate documents. The studio
   is one document now, so this loads the scripts index.html loads, in that
   order, and then renders each registered tab.

   Three things it checks that the per-page version could not:

     no repeated call   a render must not fetch one URL twice, and switching to
                        a tab already rendered must fetch nothing at all
     superseding        changing a setting mid-flight must not leave the answer
                        to the abandoned value on screen
     a pending label    every endpoint slower than the threshold must report
                        something for the indicator to show
*/
const fs = require("fs"), path = require("path"), vm = require("vm"), http = require("http");

const realTimeout = setTimeout;
const wait = ms => new Promise(r => realTimeout(r, ms));

const BASE = (process.argv[2] || "http://127.0.0.1:8600").replace(/\/$/, "");
const WEB = path.join(__dirname, "web");
const F = require(path.join(__dirname, "fakedom.js"));

F.El.prototype.append = function (...cs) { cs.forEach(c => this.children.push(c)); };
F.El.prototype.prepend = function (...cs) { cs.forEach(c => this.children.unshift(c)); };
if (!F.El.prototype.addEventListener)
  F.El.prototype.addEventListener = function () {};
if (!F.El.prototype.closest)
  F.El.prototype.closest = function () { return null; };

function httpGet(url) {
  return new Promise((res, rej) => {
    http.get(url, r => {
      let b = "";
      r.on("data", d => b += d);
      r.on("end", () => res({ok: r.statusCode < 400, status: r.statusCode, body: b}));
    }).on("error", rej);
  });
}

/* Every fetch the app makes, in order, so a repeat is visible. */
let fetchLog = [], postLog = [];
global.fetch = async (u, opt) => {
  const url = String(u);
  /* one exception: Figure 1 starts its own local build when stale or missing
     (tab_fig1.js f1Rebuild). Answered here, never sent: no check builds it. */
  if (opt && opt.method === "POST" && url === "/api/rebuild/fig1") {
    postLog.push(url);
    return {ok: false, status: 409, async json() { return {started: false, msg: "check_app does not build"}; }};
  }
  if (opt && opt.method === "POST") throw new Error("a render must not POST");
  fetchLog.push(url);
  const r = await httpGet(url.startsWith("http") ? url : BASE + (url.startsWith("/") ? "" : "/") + url);
  return {ok: r.ok, status: r.status,
          async json() { return JSON.parse(r.body); },
          async text() { return r.body; }};
};
global.alert = () => {};
global.AbortController = class { constructor() { this.signal = {aborted: false}; }
                                 abort() { this.signal.aborted = true; } };
global.URLSearchParams = URLSearchParams;
global.history = {replaceState: () => {}};
global.location = {hash: "", pathname: "/", search: ""};
global.window = global;
global.setInterval = () => 0;
global.getComputedStyle = () => ({getPropertyValue: () => "#888"});
global.Path2D = class { moveTo() {} lineTo() {} arc() {} closePath() {} };

let store = {};
global.localStorage = {getItem: k => (k in store ? store[k] : null),
                       setItem: (k, v) => { store[k] = String(v); },
                       removeItem: k => { delete store[k]; }};

let pass = 0, fail = 0;
const ck = (c, m) => { c ? pass++ : fail++; console.log((c ? "  ok   " : "  FAIL ") + m); };

const LISTED = fs.readFileSync(path.join(WEB, "index.html"), "utf8")
  .match(/<script src="([^"]+)"><\/script>/g)
  .map(s => s.replace(/.*src="|"><\/script>/g, ""));
/* A tab listed before its file exists 404s in a browser and boots without it;
   do the same here rather than fail. */
const MISSING = LISTED.filter(s => !fs.existsSync(path.join(WEB, s)));
const SCRIPTS = LISTED.filter(s => !MISSING.includes(s));
const GROUPS = ["explore", "figures", "supp", "tools"];

async function main() {
  const up = await httpGet(BASE + "/api/index").catch(() => null);
  if (!up || !up.ok) {
    console.log(`server not answering at ${BASE} -- start figures_gui/server.py first`);
    process.exit(2);
  }

  F.REG.clear();
  const body = F.$el("body");
  global.document.body = body;
  const ctx = vm.createContext(global);
  console.log(`loading ${SCRIPTS.length} scripts`
              + (MISSING.length ? ` (not on disk yet, skipped: ${MISSING.join(", ")})` : ""));
  for (const s of SCRIPTS) {
    try {
      new vm.Script(fs.readFileSync(path.join(WEB, s), "utf8"), {filename: s}).runInContext(ctx);
    } catch (e) {
      ck(false, `${s}: ${e.message}`);
      process.exit(1);
    }
  }
  ck(true, `all ${SCRIPTS.length} scripts loaded without throwing`);

  /* The modules declare `tabs` and `settings` with const, so they live in the
     context's lexical scope, exactly as separate <script> tags share one in a
     browser. They are not properties of globalThis, so reach them by evaluating
     in the same context rather than off `global`. */
  const peek = expr => new vm.Script(`(${expr})`).runInContext(ctx);
  const tabs = peek("tabs");
  ck(Array.isArray(tabs) && tabs.length >= 10,
     `${tabs ? tabs.length : 0} tabs registered`);

  // one source of truth: no tab may define its own settings copy
  ck(typeof peek("settings") === "object", "a single settings object exists");
  const ids = tabs.map(t => t.id);
  ck(new Set(ids).size === ids.length, `tab ids are unique: ${ids.join(", ")}`);
  const views = tabs.map(t => t.viewId);
  ck(new Set(views).size === views.length, "view ids are unique");

  /* index.html awaits this before the first render, so the checker must too:
     the rail's build-parameter sliders take their stops from the grid's axes,
     and with none they would be dead rather than wrong. */
  /* index.html's own inline script, which the checker does not run: both
     wirings are part of the boot, and a selection that never clears or a rise
     band whose two ends drift apart would otherwise go unseen. */
  peek("wireRiseBand()");
  peek("wireCladeSelection()");
  await peek("loadAxes()");
  const axes = peek("AXES");
  const upAxes = (JSON.parse(up.body) || {}).axes || {};
  ck(axes && Object.keys(axes).length > 0,
     `build grid axes loaded: ${Object.keys(axes || {}).join(", ")}`);
  ck(JSON.stringify(axes) === JSON.stringify(upAxes),
     "the rail's axes are the server's, not a second copy");

  /* S2 over m_T measures Figure 5 at seven (k, m_T) pairs, minutes each where
     the runs are not warmed (the cluster warms all 21); over k it is three */
  peek(`settings.saVary = "k"`);
  peek("startStudio()");

  /* Nav groups: the shell's order is the contract, every tab names one of
     them, and the nav draws them in that order with tabs in registration order. */
  ck(JSON.stringify(peek("TAB_GROUPS.map(g => g.id)")) === JSON.stringify(GROUPS),
     `shell defines the groups ${GROUPS.join(", ")} in that order`);
  const badGroup = tabs.filter(t => !GROUPS.includes(t.group)).map(t => t.id);
  ck(badGroup.length === 0, "every tab names a valid group"
     + (badGroup.length ? ` -- not: ${badGroup.join(", ")}` : ""));
  ck(!ids.includes("figure"), "the retired Figure tab is not registered");
  const nav = (F.REG.get("body").children || []).find(c => c.className === "tabs");
  const navHtml = nav ? nav.innerHTML : "";
  const navGroups = [...navHtml.matchAll(/data-group="([^"]+)"/g)].map(m => m[1]);
  ck(JSON.stringify(navGroups) === JSON.stringify(GROUPS),
     `nav shows the four groups in order: ${navGroups.join(", ")}`);
  const labels = [...navHtml.matchAll(/class="glab">([^<]+)</g)].map(m => m[1]);
  ck(JSON.stringify(labels) === JSON.stringify(
       ["Data exploration", "Figures", "Supplementary", "Tools"]),
     `nav group labels: ${labels.join(" | ")}`);
  const navTabs = [...navHtml.matchAll(/data-tab="([^"]+)"/g)].map(m => m[1]);
  const wantTabs = GROUPS.flatMap(g => tabs.filter(t => t.group === g).map(t => t.id));
  ck(JSON.stringify(navTabs) === JSON.stringify(wantTabs),
     `nav order, registration order inside a group: ${navTabs.join(", ")}`);
  const figs = tabs.filter(t => t.group === "figures").map(t => t.id).filter(i => i !== "fig1");
  ck(JSON.stringify(figs) === JSON.stringify(["clades4", "panel2", "panel1", "panel3"]),
     `figures are clades4, panel2, panel1, panel3 after fig1: ${figs.join(", ")}`);

  for (const tab of tabs) {
    fetchLog = [];
    peek(`showTab(${JSON.stringify(tab.id)})`);
    /* panel 4 takes about 8 s, panel + sens at ~4 s each, past the old 4 s */
    for (let i = 0; i < 600; i++) {
      await wait(25);
      const m = F.REG.get(`${tab.id}_main`);
      if (m && m.innerHTML.length > 300 && !/class="loading"/.test(m.innerHTML)) break;
    }
    const main = F.REG.get(`${tab.id}_main`);
    /* errorBox() is appended as an element, so it lands in children rather than
       innerHTML. Fold both in, or a failed render reads as a silent empty one. */
    const kids = main ? (main.children || []) : [];
    const kidText = kids.map(c => (c.className === "err" ? "ERROR: " : "")
                                  + (c.textContent || c.innerHTML || "")).join(" ");
    const html = (main ? main.innerHTML : "") + kidText;
    if (/ERROR:/.test(kidText)) console.log("         " + kidText.slice(0, 200));
    const label = `${tab.id}`;
    ck(!!main, `${label}: built its main region`);
    ck(html.length > 200, `${label}: rendered ${html.length}B`);
    ck(!/class="loading"/.test(html), `${label}: left the loading state`);
    const errs = /class="err"/.test(html) || /ERROR:/.test(kidText);
    ck(!errs, `${label}: no error panel`);
    const LEAK = />\s*(undefined|NaN)\s*</i, LEAK2 = /(undefined|NaN)\s*(%|d\b|days)/i;
    const bad = LEAK.test(html) || LEAK2.test(html) || /\[object Object\]/.test(html);
    if (bad) {
      const m2 = html.match(/.{0,70}(>\s*(undefined|NaN)\s*<|(undefined|NaN)\s*(%|d\b)|\[object Object\]).{0,40}/i);
      console.log("         context: " + (m2 ? m2[0].replace(/\s+/g, " ") : "?"));
    }
    ck(!bad, `${label}: no undefined/NaN leaked into the markup`);

    const dupes = fetchLog.filter((u, i) => fetchLog.indexOf(u) !== i);
    ck(dupes.length === 0,
       `${label}: fetched ${fetchLog.length} URL(s), none twice`
       + (dupes.length ? ` -- repeated ${dupes[0]}` : ""));

    /* A second visit must be free -- unless the first ERRORED, in which case
       retrying is the wanted behaviour and refetching is correct. */
    if (!errs) {
      fetchLog = [];
      peek(`showTab(${JSON.stringify(tab.id)})`);
      await wait(60);
      ck(fetchLog.length === 0, `${label}: revisiting fetched nothing`);
    }
  }

  const f1Stale = !!(peek("f1Last") && peek("f1Last").index && peek("f1Last").index.stale);
  ck(postLog.length === (f1Stale ? 1 : 0),
     `Figure 1 asks for its local build only when stale: ${postLog.length} POST(s), stale=${f1Stale}`);

  // every endpoint the app uses must have a pending label to show
  const labelled = Object.keys(peek("API_LABELS") || {});
  const used = [...new Set(SCRIPTS.flatMap(s =>
    [...fs.readFileSync(path.join(WEB, s), "utf8")
       .matchAll(/["'`](\/api\/[a-z0-9_]+)/g)].map(m => m[1])))];
  const unlabelled = used.filter(u => !labelled.includes(u));
  ck(unlabelled.length === 0,
     `every endpoint has a pending label` + (unlabelled.length ? ` -- missing ${unlabelled}` : ""));

  /* The Muller's band hit test, against the geometry that was actually drawn:
     the middle of each band, at its thickest sample, must map back to that
     band. An off-by-one in the x or y mapping selects a neighbouring clade,
     which looks like a working feature pointing at the wrong lineage. */
  peek('showTab("sweep")');
  for (let i = 0; i < 160; i++) {
    await wait(25);
    if (peek("sweepMuller")) break;
  }
  const mull = peek("sweepMuller");
  let hitBad = 0, hitN = 0;
  for (const b of (mull && mull.bands) || []) {
    if (b.id === null) continue;
    let j = 0, thick = -1;
    for (let k = 0; k < b.lo.length; k++)
      if (b.hi[k] - b.lo[k] > thick) { thick = b.hi[k] - b.lo[k]; j = k; }
    if (!(thick > 0.01)) continue;        // too thin to click at all
    hitN++;
    if (mull.hit(mull.sx(j), mull.sy((b.lo[j] + b.hi[j]) / 2)) !== b.id) hitBad++;
  }
  ck(hitN > 0 && hitBad === 0,
     `Muller hit test: ${hitN - hitBad}/${hitN} bands select themselves`);
  ck(mull && mull.hit(-5, -5) === null && mull.hit(1e6, 1e6) === null,
     "a click outside the stack selects nothing");

  /* Selecting one, end to end: the painter must report it as selected and the
     cell that matches it must be marked, or the band pops out against a grid
     where nothing says which row it is. */
  /* Not every band has a cell: the founding clade and clades with no defining
     mutation are outside fig2's pool by design, and the biggest band is
     usually the founder. So try bands until one lands in the grid. */
  const bandIds = mull.bands.filter(b => b.id !== null).map(b => b.id);
  let popped = null, marked = null;
  for (const id of bandIds) {
    peek(`setSetting("clade", ${JSON.stringify(id)})`);
    await wait(220);
    if (peek("sweepMuller").selected === id) popped = id;
    /* The cells live in #rl_grid, whose innerHTML the tab sets directly; the
       fake DOM does not fold a child's markup back into its parent's. */
    const gr = F.REG.get("rl_grid");
    const h = gr ? gr.innerHTML : "";
    if (/rcell[^"]*\bsel\b/.test(h)) { marked = id; break; }
  }
  ck(popped !== null,
     `the painter pops out the selected clade: ${JSON.stringify(popped)}`);
  ck(marked !== null,
     `the matching cell in the grid is marked: ${JSON.stringify(marked)}`);

  /* And it cannot outlive the run it names. */
  peek("setSetting('seed', 3)");
  ck(peek("settings.clade") === null, "changing the seed clears the selection");
  peek("setSetting('seed', 0)");

  /* The boot path: index.html snaps every build parameter before the first
     render. Run LAST, because it moves the whole rail off the defaults and the
     tabs above would then be measuring a different build. */
  for (const n of ["k", "mt", "win", "prom", "lo", "hi", "window"])
    peek(`setSetting(${JSON.stringify(n)}, ${JSON.stringify(
      Array.isArray(upAxes[n]) && typeof upAxes[n][0] === "number" ? 999 : "nope")})`);
  peek("snapBuildSettings()");
  const offGrid = ["k", "mt", "win", "prom", "lo", "hi", "window"]
    .filter(n => !(upAxes[n] || []).includes(peek(`settings.${n}`)));
  ck(offGrid.length === 0,
     "every build parameter lands on the grid before the first render"
     + (offGrid.length ? ` -- off: ${offGrid.join(", ")}` : ""));

  /* A build parameter may not come to rest between the grid's stops. A saved
     prom=0.5 against a grid of [0.10, 0.20] left the Sweep tab on the 9.2-day
     fallback, where m_F is never applied and nothing greys out. */
  for (const [name, bad, want] of [["prom", 0.5, 0.2], ["k", 15, 10],
                                   ["win", 61, 21], ["lo", 0.25, 0.1]]) {
    peek(`setSetting(${JSON.stringify(name)}, ${bad})`);
    peek(`sliderControl({setting: ${JSON.stringify(name)}, label: "x",
                         values: AXES[${JSON.stringify(name)}]})`);
    const got = peek(`settings.${name}`);
    ck(got === want, `${name}=${bad} snaps to the grid: ${got} (want ${want})`);
  }
  /* back to the operating point: the renders below must not need a build at
     whichever grid corner the checks above left behind (m_T among them) */
  peek(`setSettings(Object.fromEntries(BUILD_SETTINGS.map(n => [n, SETTING_DEFAULTS[n]])))`);

  /* A control that decides what the rail OFFERS has to rebuild it. The
     distribution box is greyed by the selected plot form, and a rail is built
     once per tab -- so without the rebuild the box kept whatever state it had
     when the tab was first opened, and looked tickable under a stacked bar
     chart that cannot carry per-run points. */
  for (const [tabId, setting, form, want] of [
        ["panel1", "p1Pal", "cvd", null],
        ["panel1", "p1Pal", "texture", null],
        ["panel1", "p1Pal", "normal", null],
        ["panel3", "f5Brackets", 0, null],
        ["panel3", "f5Brackets", 1, null],
        ["panel3", "f5Scale", "linear", null],
        ["panel3", "f5Scale", "loglog", null]]) {
    peek(`setSetting(${JSON.stringify(setting)}, ${JSON.stringify(form)})`);
    peek(`showTab(${JSON.stringify(tabId)})`);
    await wait(400);
    /* innerHTML is stored, not parsed (figures_gui/fakedom.js:1), so a
       parent does not carry its appended children's markup -- walk the tree. */
    const deep = el => !el ? "" :
      (el.innerHTML || "") + (el.children || []).map(deep).join("");
    const html = deep(F.REG.get(`${tabId}_rail`));
    const off = /type="checkbox"[^>]*disabled/.test(html);
    const boxes = (html.match(/Show the individual simulations/g) || []).length;
    const mainHtml = deep(F.REG.get(`${tabId}_main`));
    ck(mainHtml.length > 200 && !/class="err"|Not available/.test(mainHtml),
       `${setting}=${form}: renders without an error (${mainHtml.length}B)`
       + (mainHtml.length > 200 ? "" : ` -- at ${peek("buildQuery()")}; main: ${String(F.REG.get(`${tabId}_main`)?.textContent || "").slice(0, 160)}`));
    if (want === null) ck(boxes === 0, `${setting}=${form}: no runs box (${boxes})`);
    else ck(off === want && boxes === 1,
       `${setting}=${form}: the runs box is ${want ? "greyed" : "tickable"}`
       + ` (got ${off ? "greyed" : "tickable"}, ${boxes} box(es))`);
  }

  /* Figure 5 under the detection tabs' source counts what Cryptic by scenario
     counts: same request, and per scenario the same outbreak and cryptic totals
     as swScenCounts. */
  peek(`setSetting("p3Source", "sweep")`);
  const p3d = await peek(`apiGet(p3Url())`);
  ck(peek(`p3Url()`) === "/api/explore_panel?" + peek(`exploreQuery({mf: 0})`),
     "Figure 5 asks /api/explore_panel exactly what Cryptic by scenario asks");
  {
    ctx.__p3d = p3d;
    const cmp = peek(`(() => { const B = p3Breakouts(__p3d), C = swScenCounts(__p3d);
      return SCENARIOS.map(sc => { const s = a => a.reduce((t, x) => t + x, 0);
        return [sc, B[sc].nRuns, B[sc].cl.filter(c => c.ev).length, B[sc].cl.filter(c => c.cry).length,
                C[sc].brk.length, s(C[sc].brk), s(C[sc].cry)]; }); })()`);
    ck(cmp.every(([, a, b, c, d, e, f]) => a === d && b === e && c === f),
       `Figure 5's simulations/outbreaks/cryptic equal Cryptic by scenario's: `
       + cmp.map(r => `${r[0]} ${r[1]}/${r[2]}/${r[3]} vs ${r[4]}/${r[5]}/${r[6]}`).join("; "));
  }

  /* B-F's between-host hidden time runs from first transmission to the first day the
     smoothed curve is above f. f once stayed off the request, so the server
     measured at its own default (0.01) whatever the rail said: set f away from
     that default, and recompute the crossing from the daily series /api/explore
     sends, for every cryptic outbreak of one simulation per scenario */
  {
    peek(`setSettings({brkLo: 0.02})`);
    ctx.__p3d2 = await peek(`apiGet(p3Url())`);
    const sel = peek(`(() => { const B = p3Breakouts(__p3d2), out = [];
      SCENARIOS.forEach(sc => { const c = B[sc].cl.filter(c => c.cry && isNum(c.bh));
        const run = c.length ? c[0].run : null;
        c.filter(x => x.run === run).forEach(x => out.push([sc, x.run, x.name, x.tem, x.tTx, x.bh, x.t, x.ih])); });
      return out; })()`);
    let n = 0, bad = [], sum = 0;
    for (const [sc, run, name, tem, tTx, bh, t, ih] of sel) {
      const d = await peek(`apiGet("/api/explore?" + exploreQuery({scenario: ${JSON.stringify(sc)}, seed: ${run}, mf: 0}))`);
      const c = d.clades.find(x => x.name === name), T = d.t.slice(d.j0 || 0);
      const k = T.findIndex((v, i) => v >= tem - 1e-9 && c.ys[i] > 0.02);
      const want = k < 0 ? null : T[k] - tTx;
      n++;
      if (want === null || Math.abs(bh - want) > 0.02) bad.push(`${sc} ${name}: ${bh} vs ${want}`);
      if (Math.abs(ih + bh - t) > 0.02) sum++;
    }
    ck(n > 0 && !bad.length, `f=2%: between-host hidden time = first transmission to first day above f, `
       + `recomputed from /api/explore: ${bad.length} of ${n} off ${bad.slice(0, 3).join("; ")}`);
    ck(!sum, `intra-host + between-host hidden time = hidden time: ${sum} of ${n} off`);
  }
  peek(`setSettings({brkLo: SETTING_DEFAULTS.brkLo})`);

  /* Every share the server sends is a fraction now. A tab that forgot to
     convert shows 0.36 where 36% belongs, which no other check would catch:
     the markup is valid and the canvas draws. */
  const pctScaled = peek(`asPct({n: 5, mean: 0.36, lo: null, hi: 0.4}).mean`);
  ck(Math.abs(pctScaled - 36) < 1e-9,
     `asPct scales a fraction to a percentage: ${pctScaled}`);
  const keptNull = peek(`asPct({n: 2, mean: 0.5, lo: null, hi: null}).lo`);
  ck(keptNull === null, `asPct keeps a null bound null: ${keptNull}`);

  console.log(`\n${pass} passed, ${fail} failed`);
  process.exit(fail ? 1 : 0);
}

main();
