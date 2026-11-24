/* Figure 5: breakouts. A clade's breakout is its first sweep peaking at or
   above F; cryptic if hidden at least H days (and, under brkQuiet, never above
   its launch level in that time) -- breakoutOf in episodes.js. Clades are
   grouped by the long-born share of their defining mutations (settings.brkColour).
   Under the "sweep" source the records are the Cryptic by scenario tab's
   (swScenCounts): same request, burn-in rule, pick and skipped simulations. */

const P3_PLOTS = [
  ["breakout_counts", "Cryptic breakouts per simulation",
   "Per scenario and origin group, how many cryptic breakouts a simulation has."],
  ["breakout_scatter", "Hidden time against launch frequency",
   "One point per breakout: how long the clade was around before its sweep "
   + "against its frequency when the sweep began; cryptic ones ringed."],
];

/* a stored or linked plot3 from before v2.5.20 names a form that is gone */
const p3Snap = () => { if (!P3_PLOTS.some(([v]) => v === settings.plot3)) settings.plot3 = P3_PLOTS[0][0];
  settings.p3Source = "sweep"; };     // the clade-build source is no longer offered

function rail_panel3(rail) {
  p3Snap();
  rail.innerHTML = "";
  /* innerHTML alone empties it in a browser. The headless harness stores
     innerHTML without parsing (figures_gui/fakedom.js:1), so appended
     cards have to leave the child list too, or a rebuilt rail shows two of
     every control to check_app.js. */
  if (Array.isArray(rail.children)) rail.children.length = 0;
  /* only what the figure looks like: every data rule is set on the exploration tabs */
  rail.append(card("Plot", [
    pickControl({setting: "plot3", label: "Form", options: P3_PLOTS}),
    segmentControl({setting: "brkColour", label: "Colour by", options: P3_COLOUR,
      hint: "The long-born share of a clade's defining mutations. <em>Classes</em>: one "
          + "colour per j of k. <em>Scale</em>: white to red on the exact share. "
          + "<em>Any long-born</em>: two groups, one long-born definer is enough."}),
  ]));
  const m = document.createElement("div");
  m.className = "card";
  m.innerHTML = `<h2>Data</h2><div class="hint" style="display:block">Every rule is set on the
    exploration tabs: <a href="#sweep">Breakout detection</a> (smoothing, episodes, start, F),
    <a href="#cryptic">Cryptic breakouts</a> (hidden time, H, the launch rule) and the burn-in
    rule. The counts are <a href="#crystats">Cryptic by scenario</a>'s.</div>`;
  rail.append(m);
}

/* "Clade build" source only: the build's window, which /api/panel?n=3b measures
   on whatever the burn-in toggle says; the popup asks /api/run_daily for the same */
const p3Window = () => settings.window === "all" ? "whole record" : "from the burn-in day";
const p3Burnin = () => settings.window === "all" ? 0 : 1;

/* ------------------------------------------------------------ breakouts */

const P3_ORIG = [["std", "standard-born", () => jShade(0)], ["long", "long-born", () => jShade(1)]];
const P3_COLOUR = [["class", "Classes"], ["scale", "Scale"], ["origin", "Any long-born"]];
const P3_CLASSES = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"];
const P3_SCALE = ["#ffffff", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d"];
const p3Ramp = (stops, t) => { const x = Math.max(0, Math.min(1, t)) * (stops.length - 1);
  const i = Math.min(stops.length - 2, Math.floor(x)); return mixHex(stops[i], stops[i + 1], x - i); };

/* The groups bars and table rows are split by, and each clade's colour. Class j
   holds clades with j of k definers long-born; a clade with ndef≠k goes to the
   nearest step, never to 0 or k unless none or all of its definers are long-born. */
function p3Groups() {
  const K = +settings.k, mode = settings.brkColour;
  if (mode === "origin") return {key: c => c.o, list: P3_ORIG, head: "origin",
    colourOf: c => c.o === "long" ? jShade(1) : jShade(0)};
  const pal = mode === "scale" ? P3_SCALE : P3_CLASSES;
  return {key: c => String(c.j), head: "long-born definers",
    list: Array.from({length: K + 1}, (_, j) => [String(j), `${j}/${K}`, () => p3Ramp(pal, j / K)]),
    colourOf: c => p3Ramp(pal, mode === "scale" ? c.share : c.j / K)};
}

/* The request each source measures with. "sweep" asks exactly what Cryptic by
   scenario asks (m_F is not applied on the exploration tabs). */
const p3Url = () => settings.p3Source === "sweep"
  ? "/api/explore_panel?" + exploreQuery({mf: 0})
  : "/api/panel?n=3b&" + buildQuery();

/* Per scenario: one record per clade, each knowing its simulation (`run`, the
   seed index; `slot`, its place among the simulations kept). Under "sweep" this
   is swScenCounts's loop -- simulations without a burn-in day skipped while the
   burn-in is on, sdGateDay's gate, settings.sdPick -- kept per clade so it can
   be split by origin; check_app.js pins the totals to swScenCounts. */
function p3Breakouts(d) {
  const K = +settings.k, sweep = settings.p3Source === "sweep";
  const MF = sweep ? 0 : +settings.mf || 0, pick = sweep ? settings.sdPick : "first";
  const out = {};
  SUB_SCENARIOS.forEach(sc => {
    const runs = ((d.scenarios || {})[sc] || {}).runs || [], cl = [];
    let nRuns = 0, skipped = 0;
    runs.forEach((r, run) => {
      if (sweep && settings.burnin && r.phi == null) { skipped++; return; }
      const slot = nRuns++, gate = sweep ? sdGateDay(settings.burnin, r.phi) : null;
      r.clades.forEach(([nl, nd, tem, eps, name]) => {
        const share = nd > 0 ? nl / nd : 0;
        const j = nl <= 0 ? 0 : nl >= nd ? K : Math.min(K - 1, Math.max(1, Math.round(K * share)));
        /* m_F ("build" only) drops an episode before the pick, as /api/run_daily does */
        const ex = eps.filter(e => e[1] >= MF)
          .map(([aDay, top, amp, days, pre]) => ({aDay, top, amp, days, pre}));
        /* a clade whose breakout began on or before the gate day is left out */
        const pe = gate !== null ? breakoutPick(ex, pick) : null;
        if (pe && isNum(tem) && !(Math.max(pe.aDay, tem) > gate)) return;
        const b = breakoutOf(ex, tem, pick, gate);
        const rec = {o: nl > 0 ? "long" : "std", run, slot, seed: r.seed, nl, nd, share, j, tem, name};
        if (b) Object.assign(rec, {ev: 1, t: b.hidden, clamped: b.clamped, a: b.start,
            launch: b.launch, rise: b.rise, top: b.peak, cry: b.cry});
        else Object.assign(rec, {ev: 0, t: Math.max(0, r.end - tem)});
        cl.push(rec);
      });
    });
    out[sc] = {cl, nRuns, skipped};
  });
  return out;
}

/* Seeded bootstrap over simulations, the unit: resample n simulations, apply
   `stat` to the resampled indices. 95% interval from 1000 draws. */
function p3BootBy(n, stat) {
  if (!n) return {n: 0, mean: null, lo: null, hi: null};
  const all = Array.from({length: n}, (_, i) => i), m = stat(all);
  let seed = 12345;
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const bs = [];
  for (let b = 0; b < 1000; b++) { const v = stat(all.map(() => Math.floor(rnd() * n))); if (isNum(v)) bs.push(v); }
  bs.sort((a, b) => a - b);
  return bs.length < 2 ? {n, mean: m, lo: null, hi: null}
    : {n, mean: m, lo: bs[Math.max(0, Math.round(0.025 * bs.length) - 1)], hi: bs[Math.round(0.975 * bs.length) - 1]};
}
/* mean over simulations */
function p3Boot(vals) {
  const v = vals.filter(isNum);
  return p3BootBy(v.length, ix => ix.reduce((t, i) => t + v[i], 0) / ix.length);
}
async function p3RenderBreakout(main, idp, compact) {
  const sweep = settings.p3Source === "sweep";
  let d; try { d = await apiGetSlot(p3Url()); }
  catch (e) { main.innerHTML = ""; main.append(errorBox(e.message)); return; }
  const plot = settings.plot3, B = p3Breakouts(d), G = p3Groups();
  const F = +settings.brkF, H = +settings.brkH, quiet = !!settings.brkQuiet;
  const fPct = v => v < 0.01 ? fmt(100 * v, 2) + "%" : fmt(100 * v, 0) + "%";
  let H_ = compact ? 200 : 240;
  const inG = (sc, gk) => B[sc].cl.filter(c => G.key(c) === gk);
  const sumOf = (sc, gk, f) => inG(sc, gk).filter(c => c[f]).length;
  const perRun = (sc, gk, f) => { const v = new Array(B[sc].nRuns).fill(0);
    inG(sc, gk).forEach(c => { if (c[f]) v[c.slot]++; }); return v; };
  const shown = G.list.filter(([gk]) => SUB_SCENARIOS.some(sc => inG(sc, gk).length));
  const tot = (sc, f) => B[sc].cl.filter(c => c[f]).length;
  const pctOf = (a, b) => b ? fmt(100 * a / b, 0) + "%" : "—";
  const tmax = Math.ceil(Math.max(100, ...SUB_SCENARIOS.flatMap(sc => B[sc].cl.filter(c => c.ev).map(c => c.t))) / 100) * 100;
  const sw = c => `<span><i style="background:${c};box-shadow:inset 0 0 0 1px var(--line)"></i>`;
  const ring = plot === "breakout_scatter"
    ? `<span><i style="background:transparent;border-radius:50%;box-shadow:inset 0 0 0 1.5px var(--ink)"></i>ringed: cryptic breakout</span>` : "";
  const key = settings.brkColour === "scale"
    ? `<div class="key"><span>long-born share of defining mutations: 0%
        <i style="width:160px;background:linear-gradient(90deg,${P3_SCALE.join(",")});box-shadow:inset 0 0 0 1px var(--line)"></i>100%</span>${ring}</div>`
    : settings.brkColour === "origin"
    ? `<div class="key">${P3_ORIG.map(([, t, c]) => `${sw(c())}${t} clades</span>`).join("")}${ring}</div>`
    : `<div class="key">${shown.map(([, t, c]) => `${sw(c())}${t} long-born</span>`).join("")}${ring}</div>`;
  const colourCap = settings.brkColour === "origin"
    ? ` Colour: clade origin; <strong>long-born</strong> has at least one defining
       mutation that arose in a long shedder.`
    : ` Colour: clade origin, the share of its defining mutations that arose in a long
       shedder${settings.brkColour === "scale" ? "" : `, as j of k = ${settings.k}
       (a clade with more or fewer definers than k goes to the nearest step)`}.`;
  const anchorWord = sweep && settings.sdAnchor === "transmission" ? "first transmission" : "first appearance";
  const off = sweep && typeof sweepOffRules === "function" ? sweepOffRules() : [];
  const cryDef = `hidden ≥ ${H} d from ${anchorWord}${quiet ? ", never above its launch frequency in that time" : ""}`;
  const rule = `${sweep && settings.sdPick === "largest" ? "largest" : "first"} sweep peaking ≥ ${fPct(F)} · `
    + `cryptic = ${cryDef}`;

  const note = sc => `${tot(sc, "cry")}/${tot(sc, "ev")} breakouts cryptic · ${B[sc].nRuns} simulations`;
  const titles = {
    breakout_counts: "Cryptic breakouts per simulation, by clade origin",
    breakout_scatter: "Hidden time against launch frequency, one point per breakout",
  };
  const caps = {
    breakout_counts: `Bars: mean number of <strong>cryptic breakouts</strong> per simulation
      (${cryDef}, then a sweep peaking ≥ ${fPct(F)}); whiskers: 95% bootstrap interval over
      simulations. Groups differ in size; the hover gives each group's cryptic share of its breakouts.${colourCap}`,
    breakout_scatter: `One point per breakout. x: days from the clade's ${anchorWord} to the
      start of its breakout sweep. y: its frequency when the sweep began, log scale (floor
      0.001%). Ringed: cryptic.${quiet ? " Unringed points right of H rose above their launch frequency while hidden." : ""}
      Hover a point for its daily frequency; click to keep it open.${colourCap}`,
  };
  const scHead = sc => scenarioHead(sc, note(sc));
  const skipped = SUB_SCENARIOS.reduce((t, sc) => t + B[sc].skipped, 0);
  const scope = sweep
    ? ` Burn-in: ${!settings.burnin ? "ignored, the whole record of every simulation"
        : settings.sdGate === "start" ? "a clade counts only if its breakout starts after the burn-in day"
        : "measured from the burn-in day; clades already circulating then are left out"}${
        skipped ? `; ${skipped} simulation${skipped > 1 ? "s" : ""} that never reach it left out` : ""},
      as on <a href="#crystats">Cryptic by scenario</a>.`
    : ` Measured ${p3Window()}, the clade build's window, whatever the burn-in toggle says${
        settings.window === "all" ? "" : `: clades already circulating at the burn-in day are left
        out, and so are simulations that never reach it`}.${+settings.mf > 0
        ? ` Episodes peaking under m<sub>F</sub> = ${mtPct(settings.mf)}% do not count.` : ""}`;
  main.classList.toggle("fit", !compact);
  main.innerHTML = `<div class="sw-rule">${rule}.${scope}</div>
    <div class="plot"><div class="ph"><h3>${titles[plot]}</h3><span class="sub">k=${settings.k} · m_T=${mtPct(settings.mt)}%</span></div>
    <div class="scgrid">${SUB_SCENARIOS.map(sc => `<div>${scHead(sc)}<canvas id="${idp}_bo_${sc}"></canvas></div>`).join("")}</div>
    ${key}
    ${compact ? "" : `<div class="cap">${caps[plot]}</div>`}</div>`;

  /* the four panels share what is left of the screen: no page scroll */
  if (!compact) {
    const cvs = SUB_SCENARIOS.map(sc => document.getElementById(`${idp}_bo_${sc}`));
    const rows = (window.innerWidth || 1200) >= 720 ? 2 : 4;
    H_ = clampH((fitFree(main, cvs) - 6) / rows, 130, 340);
  }

  /* one bar per origin group present in any scenario, with its interval */
  const bars = (cv, sc, {ylab, stat, fmtY, top, label}) => {
    const {g, w, h} = fit(cv, H_);
    const pad = {l: 46, r: 10, t: 16, b: 30};
    const bs = shown.map(([gk]) => inG(sc, gk).length ? stat(gk) : {mean: null});
    const ymax = top || Math.max(1, ...bs.map(b => b.hi ?? b.mean ?? 0)) * 1.15;
    axes(g, w, h, pad, "", ylab, 0, ymax, fmtY);
    const slot = (w - pad.l - pad.r) / shown.length, few = shown.length <= 2;
    const Y = v => h - pad.b - (h - pad.b - pad.t) * v / ymax, hits = [];
    shown.forEach(([gk, t, c], i) => {
      const b = bs[i], x = pad.l + slot * (i + 0.5), bw = Math.min(56, slot * (few ? 0.5 : 0.7));
      g.fillStyle = css("--ink-2"); g.font = "600 10px " + css("--f-ui"); g.textAlign = "center";
      g.fillText(t, x, h - pad.b + 13);
      const n = inG(sc, gk).length, br = sumOf(sc, gk, "ev"), cr = sumOf(sc, gk, "cry");
      if (isNum(b.mean)) {
        g.fillStyle = c(); g.fillRect(x - bw / 2, Y(b.mean), bw, Y(0) - Y(b.mean));
        g.strokeStyle = css("--ink-3"); g.lineWidth = 0.6;
        g.strokeRect(x - bw / 2 + 0.3, Y(b.mean) + 0.3, bw - 0.6, Y(0) - Y(b.mean) - 0.6);
        if (isNum(b.lo)) { g.strokeStyle = css("--ink"); g.lineWidth = 1.2;
          g.beginPath(); g.moveTo(x, Y(b.lo)); g.lineTo(x, Y(b.hi)); g.stroke(); }
        g.fillStyle = css("--ink"); g.font = (few ? "600 10px " : "9px ") + css("--f-num");
        const lab = label(b, cr, br, few);
        if (lab) g.fillText(lab, x, Math.max(pad.t + 2, Y(b.hi ?? b.mean) - 5));
      }
      if (n) hits.push({x0: x - slot / 2, x1: x + slot / 2, y0: pad.t, y1: h - pad.b,
        text: {text: `${SCENARIO_LABEL[sc]} · ${G.head} ${t}`,
          math: `cryptic per simulation   ${ci(p3Boot(perRun(sc, gk, "cry")), 2)}\n`
              + `breakouts per simulation ${ci(p3Boot(perRun(sc, gk, "ev")), 2)}\n`
              + `cryptic of breakouts     ${cr}/${br}${br ? ` (${pctOf(cr, br)})` : ""}\n`
              + `clades                   ${n}`}});
    });
    cv._axHit = hits; axisTip(cv);
    g.textAlign = "left";
  };

  SUB_SCENARIOS.forEach(sc => {
    const cv = document.getElementById(`${idp}_bo_${sc}`);
    if (plot === "breakout_counts") bars(cv, sc, {ylab: "cryptic breakouts per simulation",
      stat: gk => p3Boot(perRun(sc, gk, "cry")), fmtY: v => v.toFixed(1),
      label: (b, cr, br, few) => few ? `${fmt(b.mean, 2)} · ${pctOf(cr, br)} of breakouts`
                                     : b.mean > 0 ? fmt(b.mean, 2) : ""});
    else {
      const floor = 1e-5;
      const {g} = scatterCI(cv, {height: H_, zero: true, xlim: [0, tmax], ylim: [-5, 0],
        xlab: "hidden days before the breakout sweep", ylab: "frequency at sweep start",
        fmtX: v => fmt(v, 0) + " d", fmtY: v => fPct(Math.pow(10, v)),
        vlines: [{x: H, colour: css("--ink-3"), label: `H ${H} d`}],
        groups: B[sc].cl.filter(x => x.ev).sort((a, b) => a.cry - b.cry || a.share - b.share).map(x => ({
          colour: G.colourOf(x), fill: G.colourOf(x), stroke: css("--ink-3"), r: 2.4, alpha: x.cry ? 0.95 : 0.6,
          x: {mean: x.t}, y: {mean: Math.log10(Math.max(floor, x.launch))},
          tip: {rec: x, text: `${SCENARIO_LABEL[sc]} · ${x.seed} · ${x.nl}/${x.nd} definers long-born`,
            math: `hidden      ${fmt(x.t, 0)} d${x.clamped ? " (start clamped)" : ""}\n`
              + `launched at ${fPct(x.launch)}\npeak        ${fPct(x.top)}\nrise        ${fmt(x.rise, 0)} d\n`
              + (x.cry ? "cryptic breakout" : "not cryptic")}})),
      });
      g.strokeStyle = css("--ink"); g.lineWidth = 0.9; g.globalAlpha = 0.8;
      (cv._ptHit || []).forEach(p => { if (!p.text.rec.cry) return;
        g.beginPath(); g.arc(p.x, p.y, 4.2, 0, 6.2832); g.stroke(); });
      g.globalAlpha = 1;
      p3Peek(cv, sc);
    }
  });
}

/* Hovering a breakout in the scatter: that clade's daily frequency from
   /api/run_daily (the record the build measures on), with its first appearance,
   hidden time, sweep start, launch and peak marked. */
const p3Daily = new Map();
function p3DailyFor(sc, run) {
  const url = settings.p3Source === "sweep"
    ? "/api/explore?" + exploreQuery({scenario: sc, seed: run, mf: 0})
    : "/api/run_daily?" + buildQuery({scenario: sc, seed: run, burnin: p3Burnin()});
  if (!p3Daily.has(url)) p3Daily.set(url, apiGet(url, {label: "reading the daily record"})
    .catch(e => { p3Daily.delete(url); throw e; }));
  return p3Daily.get(url);
}

/* one popup for every scatter: hover previews a breakout, a click pins it until
   the × , Escape, a click on empty plot, or the next render */
const P3_PEEK = {cur: null, pinned: false};
function p3PeekEl() {
  let e = document.getElementById("p3peek");
  if (!e) { e = document.createElement("div"); e.id = "p3peek"; e.className = "axtip p3peek";
    e.innerHTML = `<button type="button" class="pk-x" title="close">×</button><div class="pk-h"></div>
      <div class="pk-c"><canvas></canvas></div>${breakoutKey()}<pre class="axmath"></pre>
      <div class="pk-n">click the point to keep this open</div>`;
    e.querySelector(".pk-x").addEventListener("click", () => p3PeekHide(true));
    document.addEventListener("keydown", ev => { if (ev.key === "Escape") p3PeekHide(true); });
    window.addEventListener("hashchange", () => p3PeekHide(true));
    document.body.append(e); }
  return e;
}
function p3PeekHide(force) {
  if (P3_PEEK.pinned && !force) return;
  P3_PEEK.cur = null; P3_PEEK.pinned = false;
  const e = document.getElementById("p3peek");
  if (e) { e.style.display = "none"; e.classList.remove("pinned"); }
}
function p3PeekShow(ev, hit, sc, pin) {
  const e = p3PeekEl(), W = 420, rec = hit.text.rec;
  P3_PEEK.pinned = !!pin;
  e.classList.toggle("pinned", !!pin);
  e.style.display = "block"; e.style.width = W + "px";
  e.style.left = (ev.clientX + 16 + W > innerWidth ? ev.clientX - W - 16 : ev.clientX + 16) + "px";
  e.style.top = Math.max(8, Math.min(ev.clientY - 60, innerHeight - (e.offsetHeight || 360) - 8)) + "px";
  if (P3_PEEK.cur === rec) return;
  P3_PEEK.cur = rec;
  const key = e.querySelector(".bkey");
  if (key) key.outerHTML = breakoutKey(undefined,
    settings.p3Source === "sweep" && settings.sdAnchor === "transmission"
      ? "first transmission" : "first appearance (emergence in host)");
  e.querySelector(".pk-h").textContent = hit.text.text;
  e.querySelector("pre").textContent = hit.text.math;
  const pcv = e.querySelector("canvas");
  const {g, h} = fit(pcv, 190);
  g.fillStyle = css("--ink-3"); g.font = "11px " + css("--f-ui");
  g.fillText("reading the daily record…", 8, h / 2);
  p3DailyFor(sc, rec.run).then(d => { if (P3_PEEK.cur === rec) p3DrawPeek(pcv, d, rec); },
    err => { if (P3_PEEK.cur !== rec) return; const {g, h} = fit(pcv, 190);
      g.fillStyle = css("--ink-3"); g.font = "11px " + css("--f-ui");
      g.fillText("daily record unavailable: " + err.message, 8, h / 2); });
}

function p3Peek(cv, sc) {
  if (cv._peekBound || !cv.addEventListener) return;
  cv._peekBound = true;
  const hitAt = ev => {
    const r = cv.getBoundingClientRect(), x = ev.clientX - r.left, y = ev.clientY - r.top;
    let hit = null, best = 8;
    (cv._ptHit || []).forEach(p => { const dd = Math.hypot(p.x - x, p.y - y); if (dd < best) { best = dd; hit = p; } });
    return hit && hit.text.rec ? hit : null;
  };
  cv.addEventListener("mouseleave", () => p3PeekHide(false));
  cv.addEventListener("mousemove", ev => {
    const hit = hitAt(ev);
    cv.style.cursor = hit ? "pointer" : "";
    if (P3_PEEK.pinned) return;
    if (!hit) return p3PeekHide(false);
    p3PeekShow(ev, hit, sc, false);
  });
  cv.addEventListener("click", ev => {
    const hit = hitAt(ev);
    if (!hit || (P3_PEEK.pinned && P3_PEEK.cur === hit.text.rec)) return p3PeekHide(true);
    p3PeekShow(ev, hit, sc, true);
  });
}

function p3DrawPeek(cv, d, rec) {
  const {g, w, h} = fit(cv, 190);
  const msg = m => { g.fillStyle = css("--ink-3"); g.font = "11px " + css("--f-ui"); g.fillText(m, 8, h / 2); };
  const c = (d.clades || []).find(x => x.name === rec.name);
  if (!c) return msg("not among this run's clades in the daily record");
  const j0 = d.j0 || 0, t = d.t.slice(j0);
  const ys = c.ys ? Float64Array.from(c.ys) : smooth(c.y.slice(j0), t, +settings.win);
  if (t.length < 2) return msg("too short a record");
  const LG = -5, box = {l: 44, r: w - 10, t: 22, b: h - 20}, t0 = t[0], t1 = t[t.length - 1];
  const X = v => box.l + (box.r - box.l) * (Math.max(t0, Math.min(t1, v)) - t0) / (t1 - t0 || 1);
  const Y = v => box.b - (box.b - box.t) * (Math.log10(Math.max(1e-5, v)) - LG) / -LG;
  g.font = "9px " + css("--f-num");
  for (let p = LG; p <= 0; p++) {
    g.strokeStyle = css("--line-2"); g.lineWidth = 1;
    g.beginPath(); g.moveTo(box.l, Y(10 ** p)); g.lineTo(box.r, Y(10 ** p)); g.stroke();
    g.fillStyle = css("--ink-3"); g.textAlign = "right"; g.fillText(brkPct(10 ** p), box.l - 4, Y(10 ** p) + 3);
  }
  g.textAlign = "left"; g.fillText(`day ${fmt(t0, 0)}`, box.l, h - 5);
  g.textAlign = "right"; g.fillText(`day ${fmt(t1, 0)}`, box.r, h - 5); g.textAlign = "left";
  /* the point's own numbers (the build's), with the peak day from the daily record */
  const bd = breakoutOf(dailyEps(c.eps, t), rec.tem,
    settings.p3Source === "sweep" ? settings.sdPick : "first");
  const b = {tem: rec.tem, start: rec.a, launch: rec.launch, peak: rec.top, peakDay: bd ? bd.peakDay : null};
  drawBreakout(g, X, Y, box, null, {labels: true});
  g.strokeStyle = css("--ink"); g.lineWidth = 1.5; g.beginPath();
  const jE = emergedFrom(t, isNum(c.t_emerge) ? c.t_emerge : c.tem);
  for (let i = jE; i < t.length; i++) i > jE ? g.lineTo(X(t[i]), Y(ys[i])) : g.moveTo(X(t[i]), Y(ys[i]));
  g.stroke();
  drawBreakout(g, X, Y, box, b, {});
  g.font = "600 10px " + css("--f-ui"); g.fillStyle = css("--ink");
  const xa = X(b.tem), xs = X(b.start);
  const lab = `hidden ${fmt(rec.t, 0)} d${rec.tem < t0 ? " (appeared before the window)" : ""}`;
  const lw = g.measureText(lab).width;
  g.fillText(lab, Math.max(box.l, Math.min(box.r - lw, (xa + xs) / 2 - lw / 2)), 12);
  if (isNum(b.peakDay)) { g.font = "600 10px " + css("--f-num");
    const pt = `peak ${brkPct(b.peak)} · day ${fmt(b.peakDay, 0)}`, pw = g.measureText(pt).width, px = X(b.peakDay);
    g.fillText(pt, px + 6 + pw < box.r ? px + 6 : px - 6 - pw, Math.max(box.t + 10, Y(b.peak) + 4)); }
}

async function render_panel3(main, rail, {idp = "p3", compact = false} = {}) {
  p3Snap();
  p3PeekHide(true);
  return p3RenderBreakout(main, idp, compact);
}

registerTab({
  id: "panel3", group: "figures", num: "5", title: "Breakouts", viewId: "viewP3",
  lede: "Do clades circulate unseen, in long shedders or chains of transmission, and then suddenly break out?",
  rail: rail_panel3,
  render: render_panel3,
});
