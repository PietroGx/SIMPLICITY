/* Panel 3: ported from panel3.html.

   The page became a tab: `boot()` and its own draw() loop are gone, the
   shell calls render_panel3 when a setting it reads changes, and every
   fetch goes through apiGetSlot so a repeat costs nothing.

   The four fixed bar charts became one picked form, settings.plot3. The
   line-across-j, three-counts and dot-and-interval forms are gone (v2.5.19);
   rise against latency is the only one.

   Canvas ids carry a p3_ prefix. The shell builds every tab view in ONE
   document (shell.js, `tabs.forEach(viewElement)`), so a bare "c1" resolved to
   panel 2's canvas of the same name and this panel drew into it: panel 1's
   uniquely named canvases took 50/50/36 fillRect calls and panel 3 took zero.

   railDistribution greys its box from settings.plot3, so the rail has to be
   rebuilt when the form changes -- which settings.js's RAIL_SETTINGS and
   shell.js do. Rebuilding it from render() instead would rebuild it on every
   build-slider drag and take the handle's focus with it.
*/

const P3_PLOTS = [
  ["breakout_curve", "Breakout curves",
   "Per scenario, the share of clades that have had their first big sweep, by "
   + "days since the clade first appeared: standard-born against long-born."],
  ["breakout_counts", "Cryptic breakouts per simulation",
   "Per scenario, how many clades per simulation stayed hidden for H days at "
   + "under f of the population and then swept past F."],
  ["breakout_scatter", "Hidden time against launch frequency",
   "One point per breakout: how long the clade was around before its sweep "
   + "against how rare it was when the sweep began."],
  ["rise_latency", "Rise against latency",
   "One point per scenario: the latent interval on x against the rise time on "
   + "y, each arm that axis's own interval. Puts the waiting time and the "
   + "sweep speed on one pair of axes."],
];

function rail_panel3(rail) {
  rail.innerHTML = "";
  /* innerHTML alone empties it in a browser. The headless harness stores
     innerHTML without parsing (tests/artifacts/fakedom.js:1), so appended
     cards have to leave the child list too, or a rebuilt rail shows two of
     every control to check_app.js. */
  if (Array.isArray(rail.children)) rail.children.length = 0;
  rail.append(card("Plot", [pickControl({setting: "plot3", label: "Form", options: P3_PLOTS}),
    railDistribution("dist3", settings.plot3)]));
  if (String(settings.plot3).startsWith("breakout")) rail.append(card("Breakout", [
    sliderControl({setting: "brkF", label: "Big sweep: peak at least <i>F</i>", values: BRK_F,
      format: v => fmt(100 * v, 0) + "%",
      hint: "A clade's breakout is its FIRST episode peaking at or above F of infected hosts."}),
    sliderControl({setting: "brkH", label: "Cryptic: hidden at least <i>H</i> days", values: BRK_H,
      format: v => v + " d",
      hint: "Days from the clade's first appearance to the start of its breakout sweep."}),
    sliderControl({setting: "brkLo", label: "Cryptic: launched below <i>f</i>", values: BRK_LO,
      format: v => fmt(100 * v, 2) + "%",
      hint: "The clade's frequency when the breakout sweep began."}),
  ]));
  rail.append(railFocus());
  rail.append(railWindow());
  const r = document.createElement("div");
  r.className = "card";
  r.innerHTML = `<h2>Episode rules</h2><div class="hint">Set on the
  <a href="#sweep">Sweep time</a> tab. Currently
  <code>win=${settings.win}d prom=${settings.prom} ${settings.lo}–${settings.hi} ${settings.window}</code>.</div>`;
  rail.append(r);
}

/* distRuns holds the scenario fixed, so a chart whose CATEGORIES are the
   scenarios asks once per scenario and takes the one column back. */
const p3_runs = (dist, sc, key) => (distRuns(dist, sc, [key], key) || [null])[0];

function p3_runsBySc(dist, key) {
  const out = SCENARIOS.map(sc => p3_runs(dist, sc, key));
  return out.some(r => r && r.length) ? out : undefined;
}

const p3_classLab = c => c === "wt" ? "wild type" : c === "other" ? "ndef≠k" : c;

/* ------------------------------------------------------------ breakouts */

const BRK_F = [0.01, 0.02, 0.05, 0.1, 0.2];
const BRK_H = [0, 15, 30, 60, 90, 120, 180];
const BRK_LO = [0.0001, 0.0003, 0.001, 0.003, 0.01];
const P3_ORIG = [["std", "standard-born", () => jShade(0)], ["long", "long-born", () => jShade(1)]];

/* Per scenario and origin: one record per clade, and per-simulation counts.
   The sweep start is clamped at first appearance: an episode's start is the
   dip before its peak, which for a clade at zero until it appears lands at the
   start of the run of zeros, up to ~250 days before the clade exists. */
function p3Breakouts(d) {
  const F = +settings.brkF, H = +settings.brkH, LO = +settings.brkLo;
  const out = {};
  SUB_SCENARIOS.forEach(sc => {
    const S = {std: {cl: [], runs: []}, long: {cl: [], runs: []}};
    ((d.scenarios[sc] || {}).runs || []).forEach(r => {
      const cnt = {std: {n: 0, brk: 0, cry: 0}, long: {n: 0, brk: 0, cry: 0}};
      r.clades.forEach(([nl, nd, tem, eps]) => {
        const o = nl > 0 ? "long" : "std", e = eps.find(x => x[1] >= F);
        const rec = {o, seed: r.seed, nl, nd};
        cnt[o].n++;
        if (e) {
          Object.assign(rec, {ev: 1, t: Math.max(e[0], tem) - tem, clamped: e[0] < tem,
            launch: Math.max(0, e[1] - e[2]), rise: e[3], top: e[1]});
          rec.cry = rec.t >= H && rec.launch < LO;
          cnt[o].brk++; if (rec.cry) cnt[o].cry++;
        } else Object.assign(rec, {ev: 0, t: Math.max(0, r.end - tem)});
        S[o].cl.push(rec);
      });
      S.std.runs.push(cnt.std); S.long.runs.push(cnt.long);
    });
    out[sc] = S;
  });
  return out;
}

/* Kaplan-Meier: share broken out by day t, a clade that never breaks out
   counting until its simulation ends and no further. */
function p3KM(cl) {
  const a = cl.slice().sort((x, y) => x.t - y.t || y.ev - x.ev);
  let n = a.length, S = 1;
  const pts = [[0, 0]];
  for (let i = 0; i < a.length;) {
    const t = a[i].t; let dd = 0, cc = 0;
    while (i < a.length && a[i].t === t) { a[i].ev ? dd++ : cc++; i++; }
    if (dd) { S *= 1 - dd / n; pts.push([t, 100 * (1 - S)]); }
    n -= dd + cc;
  }
  return pts;
}
const p3At = (pts, t) => { let v = 0; for (const [x, y] of pts) { if (x > t) break; v = y; } return v; };

/* mean over simulations with a seeded bootstrap interval: the simulation is the unit */
function p3Boot(vals) {
  const v = vals.filter(isNum);
  if (!v.length) return {n: 0, mean: null, lo: null, hi: null};
  const m = v.reduce((a, b) => a + b, 0) / v.length;
  let seed = 12345;
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const bs = [];
  for (let b = 0; b < 1000; b++) { let t = 0; for (let i = 0; i < v.length; i++) t += v[Math.floor(rnd() * v.length)]; bs.push(t / v.length); }
  bs.sort((a, b) => a - b);
  return {n: v.length, mean: m, lo: bs[24], hi: bs[974]};
}
const p3Med = v => { const a = v.filter(isNum).sort((x, y) => x - y);
  return a.length ? (a.length % 2 ? a[a.length >> 1] : (a[a.length / 2 - 1] + a[a.length / 2]) / 2) : null; };

async function p3RenderBreakout(main, idp, compact) {
  let d; try { d = await apiGetSlot("/api/panel?n=3b&" + buildQuery()); }
  catch (e) { main.innerHTML = ""; main.append(errorBox(e.message)); return; }
  const plot = settings.plot3, B = p3Breakouts(d);
  const F = +settings.brkF, H = +settings.brkH, LO = +settings.brkLo;
  const fPct = v => v < 0.01 ? fmt(100 * v, 2) + "%" : fmt(100 * v, 0) + "%";
  const H_ = compact ? 200 : 240;
  const sumOf = (sc, o, k) => B[sc][o].runs.reduce((a, r) => a + r[k], 0);
  const km = {};
  SUB_SCENARIOS.forEach(sc => P3_ORIG.forEach(([o]) => { km[sc + o] = p3KM(B[sc][o].cl); }));
  const tmax = Math.ceil(Math.max(100, ...SUB_SCENARIOS.flatMap(sc => P3_ORIG.flatMap(([o]) =>
    B[sc][o].cl.map(c => c.t)))) / 100) * 100;
  const rule = `first sweep peaking ≥ ${fPct(F)} · cryptic = hidden ≥ ${H} d and launched < ${fPct(LO)}`;

  const notes = {
    breakout_curve: sc => P3_ORIG.map(([o, t]) => `${t} ${fmt(p3At(km[sc + o], 1e9), 0)}% ever`).join(" · "),
    breakout_counts: sc => P3_ORIG.map(([o, t]) => `${t} ${sumOf(sc, o, "cry")}/${sumOf(sc, o, "brk")} cryptic`).join(" · "),
    breakout_scatter: sc => P3_ORIG.map(([o, t]) => `${t} ${sumOf(sc, o, "cry")}/${sumOf(sc, o, "brk")} cryptic`).join(" · "),
  };
  const titles = {
    breakout_curve: "Share of clades broken out, by days since first appearance",
    breakout_counts: "Cryptic breakouts per simulation",
    breakout_scatter: "Hidden time against launch frequency, one point per breakout",
  };
  const caps = {
    breakout_curve: `A clade's <strong>breakout</strong> is its first sweep whose peak
      reaches ${fPct(F)} of infected hosts. Each curve is the share of clades that have
      broken out by a given number of days after they first appeared — green for
      clades whose defining mutations all arose in standard hosts, red for clades
      with at least one long-born definer. A curve that lags means clades of that
      origin stay hidden longer; its end height is how many ever break out. A clade
      that never breaks out counts only until its simulation ends (Kaplan–Meier),
      so late-appearing clades are not read as failures. Clades pooled over the
      scenario's simulations. The dashed rule is H = ${H} d.`,
    breakout_counts: `Per simulation, the number of <strong>cryptic breakouts</strong>:
      clades that existed for at least ${H} days, were under ${fPct(LO)} of infected
      hosts when their breakout sweep began, and then peaked at ${fPct(F)} or more.
      Bars are the mean over simulations with a 95% bootstrap interval over
      simulations; the label is how many of all breakouts of that origin were
      cryptic. Long-born clades are more numerous in most scenarios, so read the
      share as well as the count.`,
    breakout_scatter: `One point per breakout. x: days from the clade's first
      appearance to the start of its breakout sweep. y: the clade's frequency when
      the sweep began, on a log scale (a clade that rose the moment it appeared sits
      at the floor, 0.001%). The cryptic corner is right of H = ${H} d and below
      f = ${fPct(LO)}. Hover a point for its simulation and numbers.`,
  };
  const scHead = sc => scenarioHead(sc, notes[plot](sc));
  const tbl = compact ? "" : `<div class="plot"><div class="ph"><h3>Breakouts by origin</h3>
    <span class="sub">${rule}</span></div>
    <table class="t"><thead><tr><th>scenario</th><th>origin</th><th>clades</th><th>broke out</th>
      <th>hidden days, median</th><th>launch frequency, median</th><th>rise days, median</th>
      <th>cryptic per simulation</th><th>cryptic share of breakouts</th><th>start clamped</th></tr></thead><tbody>
    ${SUB_SCENARIOS.flatMap(sc => P3_ORIG.map(([o, t]) => {
      const cl = B[sc][o].cl, ev = cl.filter(c => c.ev), br = sumOf(sc, o, "brk");
      if (!cl.length) return "";
      return `<tr><td>${SCENARIO_LABEL[sc]}</td><td>${t}</td><td class="n">${cl.length}</td>
        <td class="n">${fmt(100 * ev.length / cl.length, 0)}%</td>
        <td class="n">${fmt(p3Med(ev.map(c => c.t)), 0)}</td>
        <td class="n">${ev.length ? fPct(p3Med(ev.map(c => c.launch))) : "—"}</td>
        <td class="n">${fmt(p3Med(ev.map(c => c.rise)), 0)}</td>
        <td class="n">${ci(p3Boot(B[sc][o].runs.map(r => r.cry)), 2)}</td>
        <td class="n">${br ? fmt(100 * sumOf(sc, o, "cry") / br, 0) + "%" : "—"}</td>
        <td class="n">${ev.filter(c => c.clamped).length}</td></tr>`; })).join("")}</tbody></table>
    <div class="cap"><strong>Start clamped</strong> counts breakouts whose sweep start
      the build placed before the clade existed — the dip before the peak, for a
      clade at zero until it appears, lands at the start of its run of zeros. Their
      start is set to the first appearance, so their hidden time is 0: they rose as
      soon as they appeared.</div></div>`;

  main.innerHTML = `<div class="note">${rule}. Clades already circulating at the
    burn-in day are left out: their first appearance is outside the window.</div>
    <div class="plot"><div class="ph"><h3>${titles[plot]}</h3><span class="sub">k=${settings.k} · m_T=${mtPct(settings.mt)}%</span></div>
    <div class="scgrid">${SUB_SCENARIOS.map(sc => `<div>${scHead(sc)}<canvas id="${idp}_bo_${sc}"></canvas></div>`).join("")}</div>
    ${keyRow(P3_ORIG.map(([, t, c]) => [c(), t + " clades"]))}
    ${compact ? "" : `<div class="cap">${caps[plot]}</div>`}</div>${tbl}`;

  SUB_SCENARIOS.forEach(sc => {
    const cv = document.getElementById(`${idp}_bo_${sc}`);
    if (plot === "breakout_curve") {
      const {g, w, h} = fit(cv, H_);
      const pad = {l: 46, r: 10, t: 10, b: 34};
      axes(g, w, h, pad, "days since the clade first appeared", "% of clades broken out", 0, 100, v => v.toFixed(0));
      const X = t => pad.l + (w - pad.l - pad.r) * t / tmax, Y = v => h - pad.b - (h - pad.b - pad.t) * v / 100;
      g.fillStyle = css("--ink-3"); g.font = "10px " + css("--f-num"); g.textAlign = "center";
      for (let t = 0; t <= tmax; t += tmax > 600 ? 200 : 100) g.fillText(String(t), X(t), h - pad.b + 12);
      g.textAlign = "left";
      g.strokeStyle = css("--ink-3"); g.setLineDash([4, 3]); g.lineWidth = 1;
      g.beginPath(); g.moveTo(X(H), pad.t); g.lineTo(X(H), h - pad.b); g.stroke(); g.setLineDash([]);
      P3_ORIG.forEach(([o, , c]) => {
        const pts = km[sc + o]; if (pts.length < 2 && !B[sc][o].cl.length) return;
        g.strokeStyle = c(); g.lineWidth = 2; g.beginPath(); g.moveTo(X(0), Y(0));
        let yv = 0;
        pts.forEach(([t, v]) => { g.lineTo(X(t), Y(yv)); g.lineTo(X(t), Y(v)); yv = v; });
        g.lineTo(X(Math.max(...B[sc][o].cl.map(x => x.t), 0)), Y(yv)); g.stroke();
      });
    } else if (plot === "breakout_counts") {
      const {g, w, h} = fit(cv, H_);
      const pad = {l: 46, r: 10, t: 16, b: 30};
      const bs = P3_ORIG.map(([o]) => p3Boot(B[sc][o].runs.map(r => r.cry)));
      const top = Math.max(1, ...bs.map(b => b.hi ?? b.mean ?? 0)) * 1.15;
      axes(g, w, h, pad, "", "cryptic breakouts per simulation", 0, top, v => v.toFixed(1));
      const slot = (w - pad.l - pad.r) / 2, Y = v => h - pad.b - (h - pad.b - pad.t) * v / top;
      P3_ORIG.forEach(([o, t, c], i) => {
        const b = bs[i], x = pad.l + slot * (i + 0.5), bw = Math.min(56, slot * 0.5);
        g.fillStyle = css("--ink-2"); g.font = "600 10px " + css("--f-ui"); g.textAlign = "center";
        g.fillText(t, x, h - pad.b + 13);
        if (!isNum(b.mean)) return;
        g.fillStyle = c(); g.fillRect(x - bw / 2, Y(b.mean), bw, Y(0) - Y(b.mean));
        if (isNum(b.lo)) { g.strokeStyle = css("--ink"); g.lineWidth = 1.2;
          g.beginPath(); g.moveTo(x, Y(b.lo)); g.lineTo(x, Y(b.hi)); g.stroke(); }
        const br = sumOf(sc, o, "brk");
        g.fillStyle = css("--ink"); g.font = "600 10px " + css("--f-num");
        g.fillText(`${fmt(b.mean, 2)} · ${br ? fmt(100 * sumOf(sc, o, "cry") / br, 0) : "—"}% of breakouts`, x, Math.max(pad.t + 2, Y(b.hi ?? b.mean) - 5));
      });
      g.textAlign = "left";
    } else {
      const floor = 1e-5;
      scatterCI(cv, {height: H_, zero: true, xlim: [0, tmax], ylim: [-5, 0],
        xlab: "hidden days before the breakout sweep", ylab: "frequency at sweep start",
        fmtX: v => fmt(v, 0) + " d", fmtY: v => fPct(Math.pow(10, v)),
        vlines: [{x: H, colour: css("--ink-3"), label: `H ${H} d`}],
        hlines: [{y: Math.log10(LO), colour: css("--ink-3"), label: `f ${fPct(LO)}`}],
        groups: P3_ORIG.flatMap(([o, t, c]) => B[sc][o].cl.filter(x => x.ev).map(x => ({
          colour: c(), fill: c(), stroke: css("--ink-3"), r: 2.4, alpha: 0.6,
          x: {mean: x.t}, y: {mean: Math.log10(Math.max(floor, x.launch))},
          tip: {text: `${SCENARIO_LABEL[sc]} · ${x.seed} · ${t} (${x.nl}/${x.nd} long-born)`,
            math: `hidden      ${fmt(x.t, 0)} d${x.clamped ? " (start clamped)" : ""}\n`
              + `launched at ${fPct(x.launch)}\npeak        ${fPct(x.top)}\nrise        ${fmt(x.rise, 0)} d`
              + (x.cry ? "\ncryptic breakout" : "")}}))),
      });
    }
  });
}

async function render_panel3(main, rail, {idp = "p3", compact = false} = {}){
  if (String(settings.plot3).startsWith("breakout")) return p3RenderBreakout(main, idp, compact);

  let d; try { d = await apiGetSlot("/api/panel?n=3&" + buildQuery()); }
  catch(e){ main.innerHTML=""; main.append(errorBox(e.message +
    " — build these rules on the Export page.")); return; }

  const plot = settings.plot3;
  const points = !!settings.dist3 && distAllowed(plot);
  const dist = points ? await distFor(3) : null;
  const S = sc => d.scenarios[sc] || {};
  const scLabs = SCENARIOS.map(sc => SCENARIO_LABEL[sc]);
  const scKey = keyRow(SCENARIOS.map(sc=>[scenarioColour(sc),SCENARIO_LABEL[sc]]));
  const classes = [...new Set(SCENARIOS.flatMap(sc =>
    Object.keys(S(sc).by_class_rise || {})))]
    .sort((a,b)=> (a==="wt")-(b==="wt") || (a==="other")-(b==="other") || (+a)-(+b));

  const tbl = SCENARIOS.map(sc => { const s = S(sc);
    return `<tr><td>${SCENARIO_LABEL[sc]} <small style="color:var(--ink-3)">${fmt(asPctNum(s.prevalence),1)}% · ${s.duration ?? "—"}d</small></td>
      <td class="n">${s.runs ?? 0}</td>
      <td class="n">${ci(s.rise_days,1," d")}</td>
      <td class="n">${s.n_with_rise ?? 0}</td>
      <td class="n">${ci(s.latent_days,1," d")}</td>
      <td class="n">${s.n_with_latency ?? 0}</td>
      <td class="n">${s.n_truncated ?? 0}</td>
      <td class="n">${s.n_censored ?? "—"}</td></tr>`; }).join("");

  const head = `<div class="note">${simWords(d.note)}</div>`;
  const table = `
  <div class="plot"><div class="ph"><h3>Summary</h3></div>
    <table class="t"><thead><tr><th>scenario</th><th>simulations</th><th>rise</th>
      <th>clades with a rise</th><th>latency</th><th>clades with a latency</th>
      <th>truncated at burn-in</th><th>still high at the end</th></tr></thead><tbody>${tbl}</tbody></table>
    <div class="cap">The truncation columns are why the two metrics cannot share
      a denominator: a clade already rising at the burn-in day keeps its rise
      time and loses its latency. <strong>Still high at the end</strong> counts
      emergence episodes the record cuts off before the clade comes back down —
      a sweep to fixation among them. Their rise is fully observed and is
      counted; only their peak is a lower bound. A dash means the build predates
      the flag.</div></div>`;

  {
    const groups = SCENARIOS.map(sc => ({colour: scenarioColour(sc),
      label: SCENARIO_LABEL[sc], x: S(sc).latent_days, y: S(sc).rise_days,
      xruns: p3_runs(dist, sc, "latent_days") || [],
      yruns: p3_runs(dist, sc, "rise_days") || []}));
    const box = `
  <div class="plot"><div class="ph"><h3>Rise time against latent interval</h3>
    <span class="sub">one point per scenario · arms are each axis's own bootstrap CI of the mean</span></div>
    <canvas id="${idp}_c1"></canvas>
    ${scKey}
    <div class="cap"><strong>Measured:</strong> two means per scenario — days
      from first appearance to the start of the sweep on x, days from 10% to 90%
      of that episode's own amplitude on y — each with a bootstrap interval over
      simulations.<br>
      <strong>Read:</strong> down and to the right would be a clade that waited
      longer <em>and then</em> swept faster, which is what chronic shedding
      predicts: more time inside one host accumulating divergence, then a
      lineage that arrives already fit. That is the reading, not the
      measurement. Five scenarios are five design points, not a trend, and the
      two axes do not share a denominator — x is over clades with a latency,
      y over clades with a rise (Summary below).</div></div>`;
    main.innerHTML = compact ? box.replace(/<div class="cap">[\s\S]*?<\/div><\/div>$/, "</div>")
                             : head + box + table;
    scatterCI(document.getElementById(idp + "_c1"), {groups, points, height: compact ? 260 : 300,
      xlab: "latent interval (days) · clades with a latency",
      ylab: "rise time (days) · clades with a rise"});

  }

}


registerTab({
  id: "panel3", group: "figures", num: "5", title: "Rise vs latency", viewId: "viewP3",
  lede: "How long a clade lurked, then how fast it went. The <strong>emergence episode</strong> is the one containing the clade's highest frequency; its rise is measured from the dip immediately before the sweep, which is what makes it the speed of a growth event rather than a disguised lifetime.",
  rail: rail_panel3,
  render: render_panel3,
});
