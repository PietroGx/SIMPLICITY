/* Figure 5: outbreaks. A clade's outbreak is its first sweep peaking at or
   above F; cryptic if hidden at least H days (and, under brkQuiet, never above
   its launch level in that time) -- breakoutOf in episodes.js. One figure:
   A, cryptic outbreaks per simulation in every scenario, control included,
   stacked by Figure 4's origin groups, starred where it differs from control;
   B-F, each cryptic outbreak's intra-host hidden time (first appearance to first
   transmission) against its between-host hidden time (first transmission to the
   first day above f, t_cross from explore.py). The records
   are Cryptic by scenario's (swScenCounts): same request, burn-in rule, pick
   and skipped simulations. */

const P3_PLOTS = [["breakouts", "Cryptic outbreaks and hidden time", ""]];

/* a stored or linked plot3 names a form that is gone; one form since the merge */
const p3Snap = () => { settings.plot3 = P3_PLOTS[0][0]; settings.p3Source = "sweep"; };

function rail_panel3(rail) {
  p3Snap();
  rail.innerHTML = "";
  /* innerHTML alone empties it in a browser. The headless harness stores
     innerHTML without parsing (figures_gui/fakedom.js:1), so appended
     cards have to leave the child list too, or a rebuilt rail shows two of
     every control to check_app.js. */
  if (Array.isArray(rail.children)) rail.children.length = 0;
  /* only what the figure looks like: every data rule is set on the exploration tabs */
  rail.append(card("Display", [
    checkControl({setting: "f5Titles", label: "Titles", hint: "Off: only the panel letters."}),
    checkControl({setting: "f5TextCol", label: "Scenario names in colour"}),
    checkControl({setting: "f5ShowN", label: "n = under each bar in A"}),
    checkControl({setting: "f5Brackets", label: "Stars vs control in A",
      hint: "Over a bar: its scenario differs from control in cryptic outbreaks per simulation. "
          + "Two-sided Mann–Whitney (Cryptic by scenario's tests), Holm-adjusted over the "
          + "comparisons with control; * p &lt; 0.05, ** p &lt; 0.01, *** p &lt; 0.001."}),
    checkControl({setting: "f5HideEmpty", label: "Hide empty B–F panels",
      hint: "Leave out the panel of a scenario with no cryptic outbreak to draw; the rest are lettered in order."}),
    segmentControl({setting: "f5Scale", label: "Axes, B–F", options: P3_SCALE,
      hint: "<em>Linear</em>: days on both. <em>Semi-log</em>: intra-host hidden time on a log axis, "
          + "between-host hidden time linear. <em>Log-log</em>: both log."}),
  ]));
  const m = document.createElement("div");
  m.className = "card";
  m.innerHTML = `<h2>Data</h2><div class="hint" style="display:block">Every rule is set on the
    exploration tabs: <a href="#sweep">Outbreak detection</a> (smoothing, episodes, start, F),
    <a href="#cryptic">Cryptic outbreaks</a> (hidden time, H, the launch rule, f) and the burn-in
    rule. The counts and tests are <a href="#crystats">Cryptic by scenario</a>'s.</div>`;
  rail.append(m);
}

/* "Clade build" source only: the build's window, which /api/panel?n=3b measures
   on whatever the burn-in toggle says; the popup asks /api/run_daily for the same */
const p3Window = () => settings.window === "all" ? "whole record" : "from the burn-in day";
const p3Burnin = () => settings.window === "all" ? 0 : 1;

/* ------------------------------------------------------------ outbreaks */

/* Figure 4's origin groups (tab_panel1.js) in its colour-blind colours, stacked
   bottom to top; the founding clade is not among Figure 5's candidate clades,
   so it is listed only if one shows up */
const P3_CLS = [["founding", "founding clade", "<i>n</i> = 0"],
  ["standard", "standard only", "<i>j</i> = 0"],
  ["mostly_std", "mostly standard", "0 &lt; <i>j</i> ≤ <i>n</i>/2"],
  ["mostly_long", "mostly long", "<i>n</i>/2 &lt; <i>j</i> &lt; <i>n</i>"],
  ["only_long", "only long", "<i>j</i> = <i>n</i>"]];
const p3Class = (nl, nd) => nd <= 0 ? "founding" : nl <= 0 ? "standard" : nl >= nd ? "only_long"
  : 2 * nl > nd ? "mostly_long" : "mostly_std";
const p3Col = k => k === "founding" ? "#9a9a9a" : P1_CVD[k];
const P3_G = {key: c => c.cls, colourOf: c => p3Col(c.cls),
  list: P3_CLS.map(([k, t]) => [k, t, () => p3Col(k)])};
const P3_SCALE = [["linear", "Linear"], ["semilog", "Semi-log"], ["loglog", "Log-log"]];

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
  const sweep = settings.p3Source === "sweep";
  const MF = sweep ? 0 : +settings.mf || 0, pick = sweep ? settings.sdPick : "first";
  const out = {};
  SCENARIOS.forEach(sc => {
    const runs = ((d.scenarios || {})[sc] || {}).runs || [], cl = [];
    let nRuns = 0, skipped = 0;
    runs.forEach((r, run) => {
      if (sweep && settings.burnin && r.phi == null) { skipped++; return; }
      const slot = nRuns++, gate = sweep ? sdGateDay(settings.burnin, r.phi) : null;
      r.clades.forEach(([nl, nd, tem, eps, name, tEm, tTx, tCross]) => {
        const share = nd > 0 ? nl / nd : 0;
        const cls = p3Class(nl, nd);
        /* m_F ("build" only) drops an episode before the pick, as /api/run_daily does */
        const ex = eps.filter(e => e[1] >= MF)
          .map(([aDay, top, amp, days, pre, bDay]) => ({aDay, top, amp, days, pre, bDay}));
        /* a clade whose outbreak began on or before the gate day is left out */
        const pe = gate !== null ? breakoutPick(ex, pick) : null;
        if (pe && isNum(tem) && !(Math.max(pe.aDay, tem) > gate)) return;
        const b = breakoutOf(ex, tem, pick, gate, sweep ? tCross : null);
        const rec = {run, slot, seed: r.seed, nl, nd, share, cls, tem, name, tTx};
        if (b) Object.assign(rec, {ev: 1, t: b.hidden, clamped: b.clamped, a: b.start, hEnd: b.hEnd,
            ceil: b.ceil, launch: b.launch, rise: b.rise, top: b.peak, cry: b.cry, tCross,
            ih: isNum(tEm) && isNum(tTx) ? tTx - tEm : null, bh: isNum(tTx) && isNum(tCross) ? tCross - tTx : null,
            pk: b.peakDay, tp: isNum(tTx) && isNum(b.peakDay) ? b.peakDay - tTx : null});
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
  let d; try { d = await apiGetSlot(p3Url()); }
  catch (e) { main.innerHTML = ""; main.append(errorBox(e.message)); return; }
  const B = p3Breakouts(d), G = P3_G, SC = SCENARIOS.filter(sc => B[sc]);
  const F = +settings.brkF, H = +settings.brkH;
  const fPct = v => v < 0.01 ? fmt(100 * v, 2) + "%" : fmt(100 * v, 0) + "%";
  const inG = (sc, gk) => B[sc].cl.filter(c => G.key(c) === gk);
  /* per simulation: how many of its clades have flag f, optionally in group gk */
  const perRun = (sc, f, gk) => { const v = new Array(B[sc].nRuns).fill(0);
    B[sc].cl.forEach(c => { if (c[f] && (gk === undefined || G.key(c) === gk)) v[c.slot]++; }); return v; };
  const tot = (sc, f) => B[sc].cl.filter(c => c[f]).length;
  const pctOf = (a, b) => b ? fmt(100 * a / b, 0) + "%" : "—";
  const shown = G.list.filter(([gk]) => SC.some(sc => inG(sc, gk).some(c => c.ev)));
  const anchorWord = "first appearance";
  /* log days, 0.1 to 1000 in every panel; shorter times, and the outbreaks that
     began before their first transmission, sit on the floor */
  const DMIN = 0.1, LG = v => Math.log10(Math.max(DMIN, v));
  /* B-F: cryptic outbreaks only, intra-host against between-host hidden time */
  const yOf = c => c.bh;
  const bfOf = sc => B[sc].cl.filter(c => c.cry && isNum(c.ih) && isNum(yOf(c)));
  const evs = SC.flatMap(bfOf);
  /* f5HideEmpty: B-F only for scenarios with a cryptic outbreak to draw, lettered in order */
  const SCB = settings.f5HideEmpty ? SC.filter(sc => bfOf(sc).length) : SC;
  const bfHidden = SC.filter(sc => !SCB.includes(sc));
  const bfRange = SCB.length > 1 ? `B–${FIG_LETTERS[SCB.length]}` : "B";
  const nNoPk = SC.reduce((t, sc) => t + tot(sc, "cry"), 0) - evs.length;
  /* y: whole orders of magnitude around every panel's time to peak, the same in all five */
  /* f5Scale: linear, semi-log (log x) or log-log; the same axes in every panel */
  const logX = settings.f5Scale !== "linear", logY = settings.f5Scale === "loglog";
  /* what sits on an axis: under DMIN on a log axis, below 0 on a linear y */
  const nFloorX = logX ? evs.filter(c => c.ih < DMIN).length : 0;
  const nFloorY = evs.filter(c => logY ? yOf(c) < DMIN : yOf(c) < 0).length;
  const floorNote = [nFloorX ? `intra-host hidden times under ${DMIN} d (${nFloorX})` : "",
    nFloorY ? `between-host hidden times ${logY ? `under ${DMIN} d` : "below 0"} (${nFloorY})` : ""].filter(Boolean).join(" and ");
  const nice = m => { const e = Math.pow(10, Math.floor(Math.log10(Math.max(m, 1e-9) / 4)));
    return [1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 10].map(k => k * e).find(st => 4 * st >= m); };
  const linAx = vals => { const st = nice(Math.max(1, ...vals)); return {lim: [0, 4 * st], ticks: [st, 2 * st, 3 * st]}; };
  const fx = v => logX ? LG(v) : Math.max(0, v), fy = v => logY ? LG(v) : Math.max(0, v);
  const yv = evs.map(yOf), tps = yv.map(LG);
  /* x per panel, from that scenario's own intra-host hidden times; y shared */
  const xAxOf = sc => { const v = bfOf(sc).map(c => c.ih);
    if (!logX) { const m = Math.max(1, ...v), e = Math.pow(10, Math.floor(Math.log10(m / 5)));
      const st = [1, 2, 5, 10].map(k => k * e).find(q => m / q <= 5), n = Math.max(1, Math.ceil(m / st - 1e-9));
      return {lim: [0, n * st], ticks: Array.from({length: n - 1}, (_, k) => (k + 1) * st), lin: 1}; }
    const lg = v.map(LG), lo = lg.length ? Math.floor(Math.min(...lg)) : -1;
    const hi = Math.max(lo + 1, lg.length ? Math.ceil(Math.max(...lg)) : 3);
    /* a centred tick at every order of magnitude; the chart's own corner labels are blanked (see fmtX) */
    return {lim: [lo, hi], ticks: [lo + 1e-9, ...Array.from({length: hi - lo - 1}, (_, k) => lo + 1 + k), hi - 1e-9]}; };
  const yAx = logY ? {lim: tps.length ? [Math.floor(Math.min(...tps)), Math.max(3, Math.ceil(Math.max(...tps)))] : [-1, 3]}
    : linAx(yv);
  const fmtLog = v => Math.abs(v - Math.round(v)) > 1e-6 ? "" : fmt(Math.pow(10, v), v < 0 ? 1 : 0);
  const cryDef = `hidden ≥ ${H} d from ${anchorWord}${cryHiddenText() ? ", " + cryHiddenText() : ""}`;
  const rule = `outbreak: the ${settings.sdPick === "largest" ? "largest" : "first"} sweep peaking ≥ ${fPct(F)}; `
    + `cryptic: ${cryDef}`;

  /* A's stars: Cryptic by scenario's own Mann–Whitney tests, on the same
     per-simulation totals, the pairs with control Holm-adjusted among themselves */
  let T = null, terr = null;
  if (settings.f5Brackets && !compact) {
    try { const m = await swCryMeasure(null); T = m.T ? m.T.cryptic : null; terr = m.terr; }
    catch (e) { terr = e.message; }
  }
  const vsC = (T ? T.pairs || [] : []).filter(x => (x.a === "control") !== (x.b === "control")
    && SC.includes(x.a) && SC.includes(x.b) && isNum(x.p)).sort((u, v) => u.p - v.p);
  const pC = {};
  vsC.reduce((run, x, i) => { const p = Math.max(run, Math.min(1, (vsC.length - i) * x.p));
    pC[x.a === "control" ? x.b : x.a] = p; return p; }, 0);
  const stars = p => !isNum(p) ? "" : p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : "";
  const nStar = Object.values(pC).filter(p => p < 0.05).length;

  /* -------------------------------------------------------------- markup */
  const head = (sc, i) => `<h4>${figLetter(sc, FIG_LETTERS[i + 1], false)}${settings.f5Titles
    ? figName(sc, cap1(SCENARIO_LABEL[sc]), settings.f5TextCol) : ""}</h4>`;
  const sw = c => `<i class="c4sw" style="background:${c};box-shadow:inset 0 0 0 1px var(--line)"></i>`;
  const row = (s, t) => `<div class="r">${s}<span>${t}</span></div>`;
  const originRows = shown.slice().reverse().map(([k, t, c]) => row(sw(c()),
      `${t} <span class="nw">(${(P3_CLS.find(x => x[0] === k) || [])[2]})</span>`)).join("")
    + `<div class="r" style="margin-top:6px"><span><i>n</i> = clade-defining substitutions</span></div>
       <div class="r"><span><i>j</i> = substitutions that arose in a long shedder</span></div>`;
  const legend = `<aside class="c4leg">${originRows}</aside>`;

  const skipped = SC.reduce((t, sc) => t + B[sc].skipped, 0);
  const scope = ` Burn-in: ${!settings.burnin ? "ignored, the whole record of every simulation"
      : settings.sdGate === "start" ? "a clade counts only if its outbreak starts after the burn-in day"
      : "measured from the burn-in day; clades already circulating then are left out"}${
      skipped ? `; ${skipped} simulation${skipped > 1 ? "s" : ""} that never reach it left out` : ""},
    as on <a href="#crystats">Cryptic by scenario</a>.`;
  const colourCap = ` Colour: Figure 4's origin groups, from a clade's n defining mutations, j of which
      arose in a long shedder: only long, j = n; mostly long, n/2 &lt; j &lt; n; mostly standard,
      0 &lt; j ≤ n/2; standard only, j = 0.`;
  const kw = T && T.kruskal ? ` Kruskal–Wallis across the ${T.kruskal.groups} scenarios: H = ${fmt(T.kruskal.H, 2)}, p ${
    T.kruskal.p < 0.001 ? "&lt; 0.001" : "= " + T.kruskal.p.toFixed(3)}.` : "";
  const testCap = !settings.f5Brackets ? ""
    : T ? ` Stars: the scenario differs from control in cryptic outbreaks per simulation (two-sided
        Mann–Whitney, Holm-adjusted over the ${vsC.length} comparisons with control; * p &lt; 0.05,
        ** p &lt; 0.01, *** p &lt; 0.001; ${nStar} of ${vsC.length} differ).${kw}`
    : ` <strong>Tests not available:</strong> ${terr || "no answer from /api/count_tests"}.`;
  const notes = `<strong>A:</strong> mean number of cryptic outbreaks per simulation in each scenario,
      stacked by origin group: each segment is that group's own mean, so the segments add up to the bar.
      Beside it: each segment's share of the column.${testCap}
      ${bfHidden.length ? `No panel for ${bfHidden.map(sc => SCENARIO_LABEL[sc]).join(", ")}: no cryptic outbreak to draw. ` : ""}
      <strong>${bfRange}:</strong> cryptic outbreaks only, one point each, pooled over the scenario's
      simulations, so a scenario with more simulations has more points; A is the per-simulation rate.
      x, intra-host hidden time: days from the clade's first appearance to its first transmission, inside the
      host where it arose. y, between-host hidden time: days from that first transmission to the first day the
      clade's smoothed frequency rises above f = ${fPct(+settings.brkLo)}; negative where it rose above f
      before the first transmission. ${{linear: "Linear scales", semilog: "Log x, linear y", loglog: "Log scales"}[settings.f5Scale] || "Log scales"},
      y the same in every panel, x fitted to each panel's own intra-host hidden times${floorNote ? `; ${floorNote} sit on the axis` : ""}.${nNoPk ? ` ${nNoPk} cryptic outbreaks with no first transmission on record, or never above f, are not drawn.` : ""}
      Hover a point for its daily frequency; click to keep it open.${colourCap}
      Rules: ${rule}.${scope} Per scenario: ${SC.map(sc =>
        `${SCENARIO_LABEL[sc]} ${tot(sc, "cry")} of ${tot(sc, "ev")} outbreaks cryptic over ${B[sc].nRuns} simulations`).join("; ")}.`;

  main.classList.remove("fit");     // fit mode hides captions; the notes box must show
  main.innerHTML = `<div class="plot"><div class="c4wrap"><div class="f5col">
      <div><h4>${figLetter(null, "A", false)}</h4>
        <canvas id="${idp}_bars"></canvas></div>
      <div class="f5grid">${SCB.map((sc, i) => `<div>${head(sc, i)}<canvas id="${idp}_bo_${sc}"></canvas></div>`).join("")}</div>
    </div>${legend}</div></div>
    ${compact ? "" : `<div class="plot f5notes"><div class="ph"><h3>Panel notes</h3></div><div class="cap">${notes}</div></div>`}`;

  /* columns for B-F from the figure's own width, not the window's: the rail and
     the legend take a share the window width does not know about */
  const col = main.querySelector(".f5col"), grid = main.querySelector(".f5grid");
  const cw = (col && col.clientWidth) || 900, cols = Math.max(1, Math.min(SCB.length, cw >= 900 ? 5 : cw >= 480 ? 3 : 1));
  if (grid && grid.style) grid.style.gridTemplateColumns = `repeat(${cols},minmax(0,1fr))`;

  /* A and the scatter rows share what is left of the screen: no page scroll */
  let barH = compact ? 190 : 230, scH = compact ? 190 : 240;
  if (!compact) {
    const cvs = [document.getElementById(`${idp}_bars`), ...SCB.map(sc => document.getElementById(`${idp}_bo_${sc}`))];
    const rows = Math.max(1, Math.ceil(SCB.length / cols));
    /* the figure fills the window; the notes box sits below it, one scroll away */
    const nb = main.querySelector(".f5notes"), lg = main.querySelector(".c4leg");
    if (lg) lg.style.display = "none";
    const free = fitFree(main, cvs) - 8 + (nb && nb.getBoundingClientRect ? nb.getBoundingClientRect().height + 14 : 0);
    if (lg) lg.style.display = "";
    barH = clampH(free * (rows === 1 ? 0.46 : 0.32), 150, 360);
    scH = clampH((free - barH) / rows, 140, 360);
  }

  /* ---------------------------------------------------------------- A */
  const drawBars = cv => {
    if (!cv) return;
    const {g, w, h} = fit(cv, barH);
    const totals = Object.fromEntries(SC.map(sc => [sc, p3Boot(perRun(sc, "cry"))]));
    const pad = {l: 46, r: 10, t: nStar ? 22 : 14, b: settings.f5ShowN ? 34 : 22};
    const top = t => t.mean ?? 0;
    const yraw = Math.max(0.5, ...SC.map(sc => top(totals[sc]))) * 1.18;
    const ystep = [0.125, 0.25, 0.5, 1, 1.5, 2, 2.5, 5, 10, 20, 25, 50].find(s => 4 * s >= yraw), ymax = ystep ? 4 * ystep : yraw;
    axes(g, w, h, pad, "", "average cryptic outbreaks per simulation", 0, ymax, v => v.toFixed(ystep < 0.25 ? 3 : ystep < 0.5 ? 2 : 1));
    const slot = (w - pad.l - pad.r) / SC.length, bw = Math.min(64, slot * 0.5);
    const X = i => pad.l + slot * (i + 0.5), Y = v => h - pad.b - (h - pad.b - pad.t) * v / ymax, hits = [];
    SC.forEach((sc, i) => {
      const x = X(i), t = totals[sc];
      let y0 = Y(0);
      const parts = G.list.map(([gk, lab, c]) => ({lab, c, m: B[sc].nRuns ? swMean(perRun(sc, "cry", gk)) : 0}));
      parts.forEach(({c, m}) => { if (!(m > 0)) return;
        const y1 = y0 - (Y(0) - Y(m));
        g.fillStyle = c(); g.fillRect(x - bw / 2, y1, bw, y0 - y1);
        g.strokeStyle = css("--ink-3"); g.lineWidth = 0.6; g.strokeRect(x - bw / 2 + 0.3, y1 + 0.3, bw - 0.6, y0 - y1 - 0.6);
        y0 = y1; });
      /* each segment's number, popped out to the right of the column on a leader line */
      const segs = [], sum = parts.reduce((a, p) => a + (p.m > 0 ? p.m : 0), 0);
      let yb = Y(0);
      parts.forEach(p => { if (!(p.m > 0)) return; const y1 = yb - (Y(0) - Y(p.m));
        segs.push({my: (yb + y1) / 2, txt: fmt(100 * p.m / sum, 0) + "%"}); yb = y1; });
      if (segs.length > 1) {
        const gap = 11, lx = x + bw / 2 + 14;
        segs.forEach((q, k) => q.ly = k ? Math.min(q.my, segs[k - 1].ly - gap) : Math.min(q.my, Y(0) - 4));
        const over = pad.t + 6 - segs[segs.length - 1].ly;
        if (over > 0) segs.forEach(q => q.ly += over);
        g.strokeStyle = css("--ink-3"); g.lineWidth = 0.8; g.fillStyle = css("--ink-2");
        g.font = "10px " + css("--f-num"); g.textAlign = "left";
        segs.forEach(q => { g.beginPath(); g.moveTo(x + bw / 2 - 6, q.my); g.lineTo(x + bw / 2 + 5, q.my);
          g.lineTo(lx - 3, q.ly); g.stroke(); g.fillText(q.txt, lx, q.ly + 3.5); });
      }
      g.textAlign = "center"; g.font = "600 11px " + css("--f-ui");
      g.fillStyle = settings.f5TextCol ? scenarioRaw(sc) : css("--ink-2");
      g.fillText(cap1(SCENARIO_LABEL[sc]), x, h - pad.b + 14);
      if (settings.f5ShowN) { g.fillStyle = css("--ink-3"); g.font = "10px " + css("--f-num");
        g.fillText(`n = ${B[sc].nRuns}`, x, h - pad.b + 27); }
      const br = tot(sc, "ev"), cr = tot(sc, "cry");
      hits.push({x0: x - slot / 2, x1: x + slot / 2, y0: pad.t, y1: h - pad.b,
        text: {text: `${SCENARIO_LABEL[sc]} · cryptic outbreaks`,
          math: `cryptic per simulation   ${ci(t, 2)}\n`
              + `outbreaks per simulation ${ci(p3Boot(perRun(sc, "ev")), 2)}\n`
              + `cryptic of outbreaks     ${cr}/${br}${br ? ` (${pctOf(cr, br)})` : ""}\n`
              + `simulations              ${B[sc].nRuns}${B[sc].skipped ? ` (${B[sc].skipped} left out)` : ""}\n`
              + (isNum(pC[sc]) ? `vs control, Holm p       ${pC[sc] < 0.001 ? pC[sc].toExponential(1) : fmt(pC[sc], 3)} ${stars(pC[sc])}\n` : "")
              + parts.filter(p => p.m > 0).map(p => `  ${p.lab}`.padEnd(25) + fmt(p.m, 2)).join("\n")}});
    });
    /* stars over the column of each scenario that differs from control */
    g.fillStyle = css("--ink"); g.textAlign = "center"; g.font = "700 13px " + css("--f-ui");
    SC.forEach((sc, i) => { const s = stars(pC[sc]), t = totals[sc];
      if (s) g.fillText(s, X(i), Math.max(12, Y(top(t)) - 6)); });
    g.textAlign = "left";
    cv._axHit = hits; axisTip(cv);
  };

  /* -------------------------------------------------------------- B-F */
  const drawScatter = (cv, sc, i) => {
    if (!cv) return;
    const xAx = xAxOf(sc);
    scatterCI(cv, {height: scH, xlim: xAx.lim, ylim: yAx.lim,
      xlab: "intra-host hidden time (days)", ylab: i % cols === 0 ? "between-host hidden time (days)" : "",
      fmtX: logX ? v => v === xAx.lim[0] || v === xAx.lim[1] ? "" : fmtLog(v) : v => fmt(v, 0),
      fmtY: logY ? fmtLog : v => Math.abs(v) < 1e-9 ? "" : fmt(v, 0),
      vlines: xAx.ticks.map(x => ({x, colour: css("--line-2"), dash: [1, 0]})), xticks: xAx.ticks,
      groups: bfOf(sc).sort((a, b) => a.share - b.share).map(x => ({
        colour: G.colourOf(x), fill: G.colourOf(x), stroke: css("--ink-3"), r: 2.6, alpha: 0.9,
        x: {mean: fx(x.ih)}, y: {mean: fy(yOf(x))},
        tip: {rec: x, text: `${SCENARIO_LABEL[sc]} · ${x.seed} · ${(P3_CLS.find(k => k[0] === x.cls) || [, x.cls])[1]}, j = ${x.nl} of n = ${x.nd}`,
          math: `intra-host hidden    ${fmt(x.ih, 1)} d\n`
            + `between-host hidden  ${fmt(x.bh, 1)} d  (to above f = ${fPct(+settings.brkLo)})\n`
            + `time to peak         ${fmt(x.tp, 1)} d\nhidden               ${fmt(x.t, 0)} d${x.clamped ? " (start clamped)" : ""}\n`
            + `launched at          ${fPct(x.launch)}\npeak                 ${fPct(x.top)}\nrise                 ${fmt(x.rise, 0)} d`}})),
    });
    p3Peek(cv, sc);
  };

  drawBars(document.getElementById(`${idp}_bars`));
  SCB.forEach((sc, i) => drawScatter(document.getElementById(`${idp}_bo_${sc}`), sc, i));
}

/* Hovering an outbreak in the scatter: that clade's daily frequency from
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

/* one popup for every scatter: hover previews an outbreak, a click pins it until
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
  if (key) key.outerHTML = breakoutKey(undefined, "first appearance (emergence in host)", true);
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
  if (!c) return msg("not among this simulation's clades in the daily record");
  const j0 = d.j0 || 0, t = d.t.slice(j0);
  const ys = c.ys ? Float64Array.from(c.ys) : smooth(c.y.slice(j0), t, +settings.win);
  if (t.length < 2) return msg("too short a record");
  const bd = breakoutOf(dailyEps(c.eps, t), rec.tem,
    settings.p3Source === "sweep" ? settings.sdPick : "first", null, rec.tCross);
  const b = {tem: rec.tem, start: rec.a, hEnd: rec.hEnd, ceil: rec.ceil, launch: rec.launch, peak: rec.top,
             peakDay: bd ? bd.peakDay : null};
  /* only the outbreak's own stretch: first appearance to its peak, a margin either side */
  const end = isNum(b.peakDay) ? b.peakDay : b.start, span = Math.max(1, end - b.tem);
  const t0 = Math.max(t[0], b.tem - Math.max(5, 0.08 * span)), t1 = Math.min(t[t.length - 1], end + Math.max(5, 0.08 * span));
  const inW = i => t[i] >= t0 && t[i] <= t1;
  /* linear frequency, from 0 to just above the larger of the peak and F */
  let top = Math.max(+settings.brkF, b.peak || 0);
  for (let i = 0; i < t.length; i++) if (inW(i)) top = Math.max(top, ys[i]);
  const step = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25].find(s => top * 1.08 / s <= 5) || 0.25;
  const ymax = Math.ceil(top * 1.08 / step) * step;
  const box = {l: 44, r: w - 10, t: 22, b: h - 20};
  const X = v => box.l + (box.r - box.l) * (Math.max(t0, Math.min(t1, v)) - t0) / (t1 - t0 || 1);
  const Y = v => box.b - (box.b - box.t) * Math.max(0, Math.min(ymax, v)) / ymax;
  g.font = "9px " + css("--f-num");
  for (let v = 0; v <= ymax + 1e-9; v += step) {
    g.strokeStyle = css("--line-2"); g.lineWidth = 1;
    g.beginPath(); g.moveTo(box.l, Y(v)); g.lineTo(box.r, Y(v)); g.stroke();
    g.fillStyle = css("--ink-3"); g.textAlign = "right"; g.fillText(fmt(100 * v, 0) + "%", box.l - 4, Y(v) + 3);
  }
  g.textAlign = "left"; g.fillText(`day ${fmt(t0, 0)}`, box.l, h - 5);
  g.textAlign = "right"; g.fillText(`day ${fmt(t1, 0)}`, box.r, h - 5); g.textAlign = "left";
  g.save(); g.beginPath(); g.rect(box.l, box.t - 2, box.r - box.l, box.b - box.t + 4); g.clip();
  drawBreakout(g, X, Y, box, null, {labels: true});
  /* Figure 5's two axes: first appearance to first transmission (dotted), then
     to the first day above f (dash-dot), with f itself */
  if (isNum(rec.tTx) && rec.tTx > b.tem) {
    g.strokeStyle = css("--ink-2"); g.lineWidth = 1.2; g.setLineDash([1, 2]);
    g.beginPath(); g.moveTo(X(rec.tTx), box.t); g.lineTo(X(rec.tTx), box.b); g.stroke(); g.setLineDash([]); }
  if (isNum(rec.tCross)) {
    g.strokeStyle = css("--ink-2"); g.lineWidth = 1; g.setLineDash([6, 2, 1, 2]);
    g.beginPath(); g.moveTo(X(rec.tCross), box.t); g.lineTo(X(rec.tCross), box.b);
    /* rule 2 draws f already, as its ceiling */
    if (!cryUnder()) { g.moveTo(box.l, Y(+settings.brkLo)); g.lineTo(box.r, Y(+settings.brkLo)); }
    g.stroke(); g.setLineDash([]); }
  g.strokeStyle = css("--ink"); g.lineWidth = 1.5; g.beginPath();
  const jE = emergedFrom(t, isNum(c.t_emerge) ? c.t_emerge : c.tem);
  let on = false;
  for (let i = jE; i < t.length; i++) { if (!inW(i)) continue;
    on ? g.lineTo(X(t[i]), Y(ys[i])) : g.moveTo(X(t[i]), Y(ys[i])); on = true; }
  g.stroke();
  drawBreakout(g, X, Y, box, b, {});
  g.restore();
  g.font = "600 10px " + css("--f-ui"); g.fillStyle = css("--ink");
  const xa = X(b.tem), xs = X(isNum(b.hEnd) ? b.hEnd : b.start);
  const lab = `hidden ${fmt(rec.t, 0)} d${rec.tem < t[0] ? " (appeared before the window)" : ""}`;
  const lw = g.measureText(lab).width;
  g.fillText(lab, Math.max(box.l, Math.min(box.r - lw, (xa + xs) / 2 - lw / 2)), 12);
  if (isNum(b.peakDay)) { g.font = "600 10px " + css("--f-num");
    const pt = `peak ${brkPct(b.peak)}`, pw = g.measureText(pt).width, px = X(b.peakDay);
    g.fillText(pt, px + 6 + pw < box.r ? px + 6 : px - 6 - pw, Math.max(box.t + 10, Y(b.peak) + 4)); }
}

async function render_panel3(main, rail, {idp = "p3", compact = false} = {}) {
  p3Snap();
  p3PeekHide(true);
  return p3RenderBreakout(main, idp, compact);
}

registerTab({
  id: "panel3", group: "figures", num: "5", title: "Outbreaks", viewId: "viewP3",
  lede: "Do clades circulate unseen, in long shedders or chains of transmission, and then suddenly break out?",
  rail: rail_panel3,
  render: render_panel3,
});
