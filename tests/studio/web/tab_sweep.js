/* Sweep time: every clade of one run, each rise episode highlighted.

   The artifact's tab 3. Episodes are computed in the browser from the run's own
   trajectories (see episodes.js), so the smoothing window, the prominence bar
   and the rise band answer immediately instead of needing a rebuild.

   The stacked plot at the top is drawn from the per-seed export; the cells
   below are measured and drawn on the daily record. Both cover the same days,
   so a peak in a cell is the same day as a peak above it -- but they are not
   the same sample grid, and only the cells carry Figure 5's numbers.
*/

const SWEEP_MULLER_H = 184;
const SWEEP_MAX_CELLS = 60;
const DAILY_CAP = 200;          // compute.MAX_DAILY_CLADES

/* The clade the grid was last scrolled to, so a redraw that did not change the
   selection does not scroll again. */
let sweepScrolledTo = null;

/* The last Muller painted here, kept so check_app.js can pin the band hit test
   to the geometry that was actually drawn. */
let sweepMuller = null;

function sweepRail(rail) {
  rail.append(railSimulation(50));
  rail.append(card("Clade definition", [
    (() => { const d = document.createElement("div");
      d.className = "hint";
      d.innerHTML = `Set on the <strong>Clustering</strong> tab and shared with
        every panel: k = ${settings.k}, m<sub>T</sub> = ${mtPct(settings.mt)}%.`;
      return d; })(),
    buttonControl({label: "Change it on the Clustering tab →",
      onClick: () => showTab("clustering")}),
  ]));
  rail.append(railEpisodes());
  rail.append(card("View", [
    segmentControl({setting: "episodeView", label: "Layout",
      options: [["grid", "Per clade"], ["one", "Overlay"]]}),
    segmentControl({setting: "yScale", label: "Height",
      options: [["shared", "Shared"], ["own", "Own peak"]],
      hint: "A shared ceiling lets the cells be read against each other; each "
          + "clade's own peak shows the SHAPE of a small clade's rise."}),
    segmentControl({setting: "colourBy", label: "Colour",
      options: [["origin", "By origin"], ["clade", "By clade"]]}),
    sliderControl({setting: "episodeCols", label: "Columns", min: 1, max: 6, step: 1}),
    checkControl({setting: "showRaw", label: "Show the unsmoothed trajectory",
      hint: "Grey, behind the smoothed line: how much the window is doing."}),
    checkControl({setting: "burnin", label: "Measure from the burn-in day",
      hint: "The stretch before it is washed out rather than cropped — the "
          + "shape to the left is the reason the window is there."}),
  ]));
}

function sweepStats(prep) {
  const rows = prep.rows;
  const withEp = rows.filter(r => r.main);
  const nEps = rows.map(r => r.eps.length);
  const mains = withEp.map(r => r.value);
  const multi = withEp.filter(r => r.eps.length > 1);
  const shares = multi.map(r => r.ampShare);
  const median = a => {
    if (!a.length) return NaN;
    const s = a.slice().sort((x, y) => x - y), m = s.length >> 1;
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
  };
  const cells = [
    ["clades", rows.length.toLocaleString(), ""],
    ["with an episode", withEp.length.toLocaleString(),
     rows.length ? `${Math.round(100 * withEp.length / rows.length)}%` : ""],
    ["episodes / clade", rows.length
      ? (nEps.reduce((a, b) => a + b, 0) / rows.length).toFixed(2) : "—", ""],
    ["more than one", rows.length
      ? `${Math.round(100 * nEps.filter(v => v > 1).length / rows.length)}%` : "—", ""],
    ["median clade rise", isFinite(median(mains))
      ? Math.round(median(mains)) : "—", "days"],
    ["largest ep holds", multi.length && isFinite(median(shares))
      ? `${Math.round(100 * median(shares))}%` : "—", "of the rise"],
  ];
  return `<div class="stats">` + cells.map(([k, v, s]) =>
    `<div class="stat"><div class="k">${k}</div><div class="v">${v}` +
    (s ? ` <small>${s}</small>` : "") + `</div></div>`).join("") + `</div>`;
}

function sweepKey(anyOther) {
  if (settings.colourBy === "origin")
    return keyRow([[css("--long"), "long origin"], [css("--std"), "standard origin"],
                   [css("--mixed"), "mixed (tie)"], [css("--founder"), "founding clade"]]
      .concat(anyOther ? [[css("--ink-3"), "smaller clades"]] : []));
  return `<div class="key"><span>Each colour is one clade; hues are arbitrary
    and reassigned when the partition changes.</span>` +
    (anyOther ? `<span><i style="background:${css("--ink-3")}"></i>smaller clades</span>` : "")
    + `</div>`;
}

async function renderSweep(main) {
  /* The daily record is what panel 3 measures on, so this tab asks for it and
     measures there too. It comes from the clade build, so it is unavailable
     while no build matches the current compute.py -- then the tab falls back to
     the ~9.2-day export and says so, rather than showing nothing. */
  const [runv, daily] = await Promise.all([
    runView(settings.scenario, settings.seed, settings.k, settings.mt),
    apiGetSlot("/api/run_daily?"
               + buildQuery({scenario: settings.scenario, seed: settings.seed}))
      .then(d => d, e => { if (isSuperseded(e)) throw e; return null; }),
  ]);
  const prep = prepareEpisodes(runv, daily);
  const rows = prep.rows;
  const shown = rows.slice(0, SWEEP_MAX_CELLS);
  const gridMode = settings.episodeView === "grid";

  const axisNote = prep.daily
    ? `<div class="note"><strong>These are Figure 5's numbers.</strong> The
       episodes here are measured on the daily record, the same axis
       <code>compute.py</code> uses, so a rise day on this tab is the rise day
       Figure 5 reports. The stacked plot below is still drawn from the per-seed
       export, which samples about every 9.2 days — a shape, not a number — and
       both cover the same days.
       ${prep.unserved ? `<br>${prep.unserved} clade${prep.unserved === 1 ? "" : "s"}
       of this simulation carry no daily series and are not shown: the build packs the
       largest ${DAILY_CAP} that clear m<sub>T</sub>.` : ""}</div>`
    : `<div class="note err"><strong>Falling back to the 9.2-day export.</strong>
       The daily record comes from the clade build, and no build matches the
       current <code>compute.py</code> — so every crossing day below is
       quantised to roughly 9 days and a 30-day rise is 3 samples. Measured gap,
       control seed 0 at win=1: 275 emergence episodes on the daily record
       against 13 here. Rebuild on the <strong>Export</strong> tab, then these
       become Figure 5's numbers.</div>`;
  main.innerHTML = axisNote + sweepStats(prep) + `
    <div class="plot"><div class="ph"><h3>Lineage frequency by clade</h3>
      <span class="sub" id="rl_mullerSub"></span></div>
      <canvas id="rl_muller"></canvas>
      <div class="axis"><span id="rl_mX0"></span><span id="rl_mXmid"></span>
        <span id="rl_mX1"></span></div>
      <div id="rl_mullerKey"></div></div>
    <div class="plot"><div class="ph">
      <h3>${gridMode ? "Clade frequency, each rise episode highlighted"
                     : "Every clade on one axis, each rise episode highlighted"}</h3>
      <span class="sub">${rows.length > SWEEP_MAX_CELLS
        ? `${shown.length} largest of ${rows.length} clades`
        : `${rows.length} clades`} · k = ${settings.k}, m_T = ${mtPct(settings.mt)}%</span></div>
      ${gridMode ? `<div class="rgrid" id="rl_grid"></div>`
                 : `<canvas id="rl_one"></canvas>`}
      <div class="cap">The thick orange stretch is the measured rise: the
        low-crossing to the high-crossing of that episode's own amplitude. Dots
        are the trough it rose from and the peak it reached. A faded stretch is
        an episode that is not the clade's largest.</div></div>`;

  const mcv = document.getElementById("rl_muller");
  const m = paintMuller(mcv,
    {D: runv.D, C: runv.C, ordered: prep.ordered, colMap: prep.colMap,
     height: SWEEP_MULLER_H, washBurnin: !!settings.burnin,
     selected: settings.clade,
     excluded: new Set(rows.filter(r => !r.eps.length).map(r => r.id))});
  /* Click a band to pop it out of the stack; click it again, or click the grey
     overflow band or outside the stack, to clear. The grey band is many clades
     at once, so hit() returns null there rather than picking one. */
  sweepMuller = m;
  mcv.style.cursor = "pointer";
  mcv.onclick = ev => {
    const b = mcv.getBoundingClientRect();
    const id = m.hit(ev.clientX - b.left, ev.clientY - b.top);
    setSetting("clade", id !== null && id === settings.clade ? null : id);
  };
  /* A band can be clicked that has no cell below it, and the reason matters:
     the founding clade and clades with no defining mutation are outside fig2's
     pool by design, and the daily record packs only the largest DAILY_CAP. Say
     which, rather than leaving an empty grid to explain itself. */
  const selRow = rows.find(r => r.id === settings.clade);
  const selInf = settings.clade === null ? null
               : runv.C.info.get(settings.clade);
  const selNote =
    settings.clade === null ? ""
    : selRow ? ` · showing ${labelOf(selRow.inf)} · ${selRow.inf.size}L`
    : !selInf ? " · the selected clade is not in this partition"
    : selInf.isRoot ? " · the founding clade, which carries no metric"
    : !selInf.ndef ? " · no defining mutations, so it carries no metric"
    : " · no daily series for this clade";
  document.getElementById("rl_mullerSub").textContent =
    `${prep.ordered.length.toLocaleString()} clades · ${rows.length} carry the metric`
    + selNote;
  document.getElementById("rl_mullerKey").innerHTML = sweepKey(m.anyOther);
  document.getElementById("rl_mX0").textContent = `day ${m.day0}`;
  document.getElementById("rl_mX1").textContent = `day ${m.day1}`;
  document.getElementById("rl_mXmid").textContent = settings.burnin && prep.phi !== null
    ? "the whole simulation, washed out where the metric does not look — the cells below share this axis"
    : "the whole simulation — the cells below share this axis";

  if (gridMode) {
    const grid = document.getElementById("rl_grid");
    grid.style.gridTemplateColumns = `repeat(${settings.episodeCols}, minmax(0, 1fr))`;
    grid.innerHTML = shown.map((r, i) => {
      const lab = labelOf(r.inf);
      const d = r.value !== null && r.value !== undefined
        ? `${Math.round(r.value)} d` : "no peak";
      const rest = r.eps.filter(e => e !== r.main).map(e => Math.round(e.days));
      const tail = r.eps.length === 0 ? "no episode"
        : r.eps.length === 1 ? "1 episode"
        : `${r.eps.length} eps · also ${rest.join(", ")} d`;
      /* no episode for one of two reasons, and the label says which: no peak
         cleared the prominence bar, or every one that did is under m_F */
      const why = r.eps.length ? (r.main && r.main.cens ? " · still high at the end" : "")
        : r.nFound ? ` · under m_F = ${mtPct(settings.mf)}%` : "";
      return `<div class="rcell${r.eps.length ? "" : " out"}${
        r.id === settings.clade ? " sel" : ""}" data-i="${i}">
        <div class="rh"><b title="${lab}">${lab} · ${r.inf.size}L</b>
          <span class="d">${d}</span></div>
        <canvas></canvas>
        <div class="rf"><span>peak ${(100 * r.peak).toFixed(1)}%</span>
          <span title="every episode's own rise">${tail}${why}</span></div></div>`;
    }).join("");
    const cells = grid.querySelectorAll(".rcell canvas");
    const cellH = settings.episodeCols >= 5 ? 66 : settings.episodeCols >= 3 ? 86 : 132;
    shown.forEach((r, i) => { if (cells[i]) paintOne(cells[i], r, cellH, false, prep.view); });
    /* The same selection from the other end, so a cell and its band can be
       matched in either direction. */
    grid.onclick = ev => {
      const c = ev.target.closest && ev.target.closest(".rcell");
      const r = c ? shown[+c.dataset.i] : null;
      if (r) setSetting("clade", r.id === settings.clade ? null : r.id);
    };
    /* Only when the selection CHANGED: scrolling on every redraw would fight
       with the slider the redraw came from. */
    if (settings.clade === null) sweepScrolledTo = null;
    else if (settings.clade !== sweepScrolledTo) {
      const el = grid.querySelector(".rcell.sel");
      if (el && el.scrollIntoView) el.scrollIntoView({block: "nearest"});
      sweepScrolledTo = settings.clade;
    }
  } else {
    const cv = document.getElementById("rl_one");
    const {g, w, h} = fit(cv, 300);
    const V = prep.view;
    const sx = j => 2 + (V.nt < 2 ? 0 : j / (V.nt - 1)) * (w - 4);
    shown.forEach(r => {
      const out = !r.eps.length;
      const dim = settings.clade !== null && r.id !== settings.clade;
      g.strokeStyle = (out || dim) ? css("--ink-3") : (r.col || css("--ink-3"));
      g.globalAlpha = dim ? 0.12 : (out ? 0.22 : 0.55);
      g.lineWidth = r.id === settings.clade ? 2 : 1;
      g.beginPath();
      for (let j = 0; j < r.ys.length; j++) {
        const x = sx(V.j0 + j), y = h - 5 - (r.ys[j] / V.top) * (h - 10);
        j ? g.lineTo(x, y) : g.moveTo(x, y);
      }
      g.stroke();
    });
    g.globalAlpha = 1;
    shown.forEach(r => r.eps.forEach(e => {
      const dim = settings.clade !== null && r.id !== settings.clade;
      g.strokeStyle = css("--rise"); g.lineWidth = 2.4; g.lineCap = "round";
      g.globalAlpha = dim ? 0.08 : ((r.main && e === r.main) ? 0.95 : 0.35);
      g.beginPath();
      for (let j = e.ia; j <= e.ib; j++) {
        const x = sx(V.j0 + j), y = h - 5 - (r.ys[j] / V.top) * (h - 10);
        j > e.ia ? g.lineTo(x, y) : g.moveTo(x, y);
      }
      g.stroke(); g.lineCap = "butt";
    }));
    g.globalAlpha = 1;
    if (settings.burnin && V.phiIdx !== null) {
      g.save(); g.globalAlpha = 0.62; g.fillStyle = css("--bg");
      g.fillRect(0, 0, sx(V.phiIdx), h); g.restore();
      drawBurnin(g, sx(V.phiIdx), w, h, V.phi, false);
    }
  }
}

registerTab({
  id: "sweep", group: "explore", title: "Sweep time", viewId: "viewRise",
  lede: "Every clade of one simulation, with each rise episode highlighted. The "
      + "episode rules set here are what Figure 5 reads.",
  rail: sweepRail,
  render: renderSweep,
});
