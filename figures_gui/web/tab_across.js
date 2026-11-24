/* m_T effect on clades: what the clustering does to all 250 runs, not to one.

   The artifact's tab 2. Five metrics, each drawn the same way: one faint line
   per run across the m_T values, the median over runs in bold, and a dashed
   rule at the m_T currently set. All of them arrive in one response, so
   moving the m_T slider only moves the rule -- nothing is refetched.

   The m_T values and the doubling size bins are the artifact's, so these axes are
   directly comparable with the ones it drew; the numbers behind them now come
   from the clade build, so they carry the same clade definition as the panels.
*/

/* the scenarios picked on the rail, or all of them: here a pick hides the rest */
const acrossShown = () => { const f = focusSet(); return f.size ? SCENARIOS.filter(s => f.has(s)) : SCENARIOS; };

const ACROSS_PANELS = [
  {key: "clades", label: "Clades", log: true,
   cap: "How many clades the definition yields. Falls steeply with m_T: most "
      + "groups a bare threshold produces carry almost no transmission."},
  {key: "singleton_frac", label: "Singleton share", pct: true,
   cap: "Clades holding a single lineage. At m_T = 0 most clades are "
      + "singletons; the transmission minimum is what removes them."},
  {key: "median_size", label: "Median clade size", log: true,
   cap: "Member lineages in the middle clade."},
  {key: "max_div", label: "Max divergence",
   cap: "Mutations between the furthest member and its clade root. A bound, "
      + "so the simulation's worst case rather than its typical one."},
  {key: "never_tx_frac", label: "Never transmitted", pct: true,
   cap: "Clades that infected nobody. Non-zero only at m_T = 0: above it the "
      + "definition refuses them, which is the whole point of the minimum."},
];

function acrossRail(rail) {
  railCompact(rail);
  rail.append(presetCard("ac"));
  rail.append(railClade({origin: true}));
  rail.append(railFocus({hide: true}));
}

/* One metric across the m_T values: faint per-run lines, bold medians, and the
   current m_T. Ported from the artifact's drawSweep -- same geometry, same
   log/percent handling, same dashed rule. */
function paintValuesPanel(cv, panel, data, height) {
  const {g, w, h} = fit(cv, height);
  const mtValues = data.values;
  const nm = mtValues.length;
  const L = 38, R = 6, T = 8, B = 16;
  const sx = i => L + (nm < 2 ? 0 : i / (nm - 1)) * (w - L - R);
  /* a share arrives as a fraction and is drawn as a percentage */
  const val = r => panel.pct ? asPctNum(r[panel.key]) : r[panel.key];

  let lo = Infinity, hi = -Infinity;
  acrossShown().forEach(sc => (data.scenarios[sc] || {seeds: []}).seeds.forEach(s =>
    s.per_value.forEach(r => {
      const v = val(r);
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    })));
  if (!isFinite(lo)) { lo = 0; hi = 1; }
  if (panel.pct) { lo = 0; hi = 100; }
  const logY = panel.log && lo > 0;
  if (logY) lo = Math.max(lo, 1); else if (!panel.pct) lo = Math.min(lo, 0);
  if (hi <= lo) hi = lo + 1;
  const tr = v => logY ? Math.log10(Math.max(v, lo)) : v;
  const y0 = tr(lo), y1 = tr(hi);
  const sy = v => T + (1 - (tr(v) - y0) / (y1 - y0)) * (h - T - B);

  g.font = "9px " + css("--f-num");
  g.strokeStyle = css("--line"); g.fillStyle = css("--ink-3"); g.lineWidth = 1;
  const ticks = logY
    ? (() => { const t = [];
               for (let e = Math.ceil(y0); e <= y1; e++) t.push(Math.pow(10, e));
               return t.length > 1 ? t : [lo, hi]; })()
    : [lo, lo + (hi - lo) / 2, hi];
  ticks.forEach(t => {
    g.globalAlpha = 0.45;
    g.beginPath(); g.moveTo(L, sy(t)); g.lineTo(w - R, sy(t)); g.stroke();
    g.globalAlpha = 1;
    g.fillText(panel.pct ? Math.round(t) + "%"
      : t >= 1000 ? (t / 1000).toFixed(t >= 10000 ? 0 : 1) + "k"
      : String(Math.round(t)), 2, sy(t) + 3);
  });
  mtValues.forEach((m, j) => {
    if (j && j !== nm - 1 && j % 4) return;
    g.fillText(m === 0 ? "off" : mtPct(m) + "%", Math.max(0, sx(j) - 8), h - 4);
  });

  // where m_T sits now: snapped to the nearest value, as the artifact did
  let mi = 0;
  mtValues.forEach((m, j) => {
    if (Math.abs(m - settings.mt) < Math.abs(mtValues[mi] - settings.mt)) mi = j;
  });
  g.save(); g.setLineDash([4, 3]); g.strokeStyle = css("--accent"); g.globalAlpha = 0.9;
  g.beginPath(); g.moveTo(sx(mi), T); g.lineTo(sx(mi), h - B); g.stroke(); g.restore();

  acrossShown().forEach(sc => {
    const d = data.scenarios[sc];
    if (!d || !d.seeds.length) return;
    g.strokeStyle = scenarioColour(sc); g.lineWidth = 1; g.globalAlpha = 0.13;
    d.seeds.forEach(s => {
      g.beginPath();
      s.per_value.forEach((r, j) =>
        j ? g.lineTo(sx(j), sy(val(r))) : g.moveTo(sx(j), sy(val(r))));
      g.stroke();
    });
  });
  acrossShown().forEach(sc => {
    const d = data.scenarios[sc];
    if (!d || !d.median.length) return;
    g.strokeStyle = scenarioColour(sc); g.lineWidth = 1.8; g.globalAlpha = 1;
    g.beginPath();
    d.median.forEach((r, j) =>
      j ? g.lineTo(sx(j), sy(val(r))) : g.moveTo(sx(j), sy(val(r))));
    g.stroke();
  });
  g.globalAlpha = 1;
}

/* the built m_T value nearest the one set */
const acrossMi = data => data.values.reduce((b, m, j) =>
  Math.abs(m - settings.mt) < Math.abs(data.values[b] - settings.mt) ? j : b, 0);

/* Clade sizes at the set m_T, as a share of clades at or above each doubling
   bin. The x axis spans every value, so it holds still as the slider moves. */
function paintSizeBins(cv, data, height) {
  const {g, w, h} = fit(cv, height);
  const L = 38, R = 8, T = 10, B = 20;
  const mi = acrossMi(data);
  const bins = acrossShown().flatMap(sc => (data.scenarios[sc] || {size_bins: []}).size_bins.flat());
  const maxBin = bins.reduce((m, b) => Math.max(m, b.bin), 0);
  const sx = b => L + (maxBin < 1 ? 0 : b / maxBin) * (w - L - R);
  const sy = p => T + (1 - p / 100) * (h - T - B);

  g.font = "9px " + css("--f-num");
  [0, 25, 50, 75, 100].forEach(p => {
    g.strokeStyle = css("--line-2"); g.globalAlpha = 0.6;
    g.beginPath(); g.moveTo(L, sy(p)); g.lineTo(w - R, sy(p)); g.stroke();
    g.globalAlpha = 1;
    g.fillStyle = css("--ink-3");
    g.fillText(p + "%", 2, sy(p) + 3);
  });
  for (let b = 0; b <= maxBin; b += Math.max(1, Math.round(maxBin / 6))) {
    g.fillStyle = css("--ink-3");
    g.fillText(String(1 << b), Math.max(0, sx(b) - 6), h - 5);
  }
  acrossShown().forEach(sc => {
    const d = data.scenarios[sc];
    const sb = d && d.size_bins[mi];
    if (!sb || !sb.length) return;
    const total = sb.reduce((a, x) => a + x.n, 0);
    let cum = 0;
    const pts = [];
    for (let b = 0; b <= maxBin; b++) {
      const hit = sb.find(x => x.bin === b);
      pts.push([b, 100 * (1 - cum / Math.max(total, 1))]);
      if (hit) cum += hit.n;
    }
    g.strokeStyle = scenarioColour(sc); g.lineWidth = 1.6;
    g.beginPath();
    pts.forEach(([b, p], i) => i ? g.lineTo(sx(b), sy(p)) : g.moveTo(sx(b), sy(p)));
    g.stroke();
  });
}

function whatIsLeftTable(data) {
  const mi = acrossMi(data);
  const rows = acrossShown().map(sc => {
    const d = data.scenarios[sc];
    if (!d || !d.median.length) return `<tr><td>${SCENARIO_LABEL[sc]}</td>
      <td colspan="4">not built</td></tr>`;
    const off = d.median[0], on = d.median[mi];
    return `<tr><td>${SCENARIO_LABEL[sc]}</td>
      <td class="n">${fmt(off.clades, 0)}</td>
      <td class="n">${fmt(on.clades, 0)}</td>
      <td class="n">${fmt(100 * on.clades / Math.max(off.clades, 1), 1)}%</td>
      <td class="n">${fmt(asPctNum(off.singleton_frac), 1)}% → ${fmt(asPctNum(on.singleton_frac), 1)}%</td></tr>`;
  }).join("");
  return `<div class="plot"><div class="ph"><h3>What is left</h3>
    <span class="sub">median over simulations · m_T = ${mtPct(settings.mt)}% (value ${data.values[mi]})</span></div>
    <table class="t"><thead><tr><th>scenario</th><th>groups at m_T = 0</th>
      <th>clades now</th><th>kept</th><th>singleton share</th></tr></thead>
      <tbody>${rows}</tbody></table>
    <div class="cap">The left column is every group the threshold produced; the
      next is the clade definition at the current m_T. The gap is what the
      transmission minimum removes, and it is most of them.</div></div>`;
}

async function renderAcross(main) {
  main.classList.add("fit");
  const data = await apiGet("/api/across?" + buildQuery(), {slot: "across",
    label: "reading the clade build"});
  /* which values exist: the tab reads one build per value now, so a thin set
     is a gap in the grid rather than a result. */
  const gap = (data.values_missing || []).length
    ? `<div class="note"><strong>${data.values_missing.length} of
       ${(data.values_full || []).length} m<sub>T</sub> values have no clade
       build.</strong> Every value is its own build under rule B, so what is
       drawn below is only the values on disk: ${(data.values || []).join(", ")}%.
       Missing: ${data.values_missing.join(", ")}%. Build them with
       <code>figures_gui/slurm_precompute.sh</code> on the cluster, or one at a
       time from the <strong>Export</strong> tab.</div>`
    : "";
  main.innerHTML = gap + whatIsLeftTable(data) + `
    <div class="setgrid">
    <div class="plot"><div class="ph"><h3>Clade sizes</h3>
      <span class="sub">x: clade size · k = ${settings.k} · m_T = ${data.values.length ? mtPct(data.values[acrossMi(data)]) + "%" : "—"}</span></div>
      <canvas id="cl_ecdf"></canvas>
      <div class="cap">Share of clades with at least n lineages, at the m<sub>T</sub>
        set now: the drop from 1 to 2 is the singleton share. The only panel here
        whose x axis is not m<sub>T</sub>.</div></div>` + ACROSS_PANELS.map(p =>
      `<div class="plot"><div class="ph"><h3>${p.label}</h3>
        <span class="sub">x: m_T · per simulation</span></div>
        <canvas data-panel="${p.key}"></canvas>
        <div class="cap">${p.cap}</div></div>`).join("") + `</div>
    ${keyRow(acrossShown().map(sc => [scenarioColour(sc), SCENARIO_LABEL[sc]]))}
    <div class="cap">Each faint line is one simulation across the m<sub>T</sub>
      values; the bold line is the median <strong>over simulations</strong>. The
      dashed rule is where m<sub>T</sub> is set now. Scenarios are never pooled.
      Every value is its own clade build: under rule B, m<sub>T</sub> is part of
      designation, so one cannot be had by filtering another.</div>`;

  /* one height for every panel, from the window's leftover over the grid's rows */
  const ecdf = document.getElementById("cl_ecdf");
  const cvs = ACROSS_PANELS.map(p => main.querySelector(`canvas[data-panel="${p.key}"]`));
  const free = fitFree(main, [ecdf, ...cvs]);
  const cells = [...(main.querySelectorAll(".setgrid > .plot") || [])];
  const rows = Math.max(1, new Set(cells.map(c => c.offsetTop)).size);
  const ph = clampH(free / rows, 80, 260);
  paintSizeBins(ecdf, data, ph);
  ACROSS_PANELS.forEach((p, i) => paintValuesPanel(cvs[i], p, data, ph));
}

registerTab({
  id: "across", group: "explore", title: "m_T effect on clades", viewId: "viewSet",
  lede: "The same clustering over all 250 simulations, across the m_T values. One "
      + "faint line per simulation, the median over simulations in bold.",
  rail: acrossRail,
  render: renderAcross,
});
