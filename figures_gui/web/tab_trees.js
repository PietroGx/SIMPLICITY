/* S1 · All trees: every simulation of one scenario as a phylogeny, one figure
   per scenario. The clade definition is the Clustering tab's (k, m_T, origin
   rule); the colouring, tree axis and palette are its Display controls. Reads
   only the per-seed export (/api/run, runs/), as the Clustering overview does. */

const TR_ASPECT = 0.62;              // cell height over cell width

function treesRail(rail) {
  railCompact(rail);
  rail.append(card("Scenario", [
    segmentControl({setting: "trScen", label: "Scenario",
      options: SCENARIOS.map(sc => [sc, SCENARIO_LABEL[sc]]),
      hint: "One figure per scenario: export each in turn."}),
  ]));
  rail.append(card("Display", [
    segmentControl({setting: "colourBy", label: "Colour by",
      options: [["origin", "Origin"], ["clade", "Clade"]],
      hint: "Origin compares across simulations: clade colours are per-simulation "
          + "and mean nothing from one cell to the next."}),
    segmentControl({setting: "treeAxis", label: "Tree x axis",
      options: [["time", "Time"], ["div", "Divergence"]]}),
    paletteControl(),
    segmentControl({setting: "trCols", label: "Per row",
      options: [[5, "5"], [7, "7"], [10, "10"]]}),
    checkControl({setting: "trTags", label: "Seed numbers"}),
    checkControl({setting: "trBurnin", label: "Burn-in day",
      hint: "A dashed rule at the simulation's burn-in day, on the time axis only."}),
  ]));
  rail.append(dataCard(DATA_CLADE));
}

/* the legend of the colouring in force, the Figure 2 legend's rows */
function trLegend() {
  if (settings.colourBy !== "origin")
    return `<span>one colour per clade, per simulation</span>`;
  const rows = settings.originRule === "frac"
    ? [[css("--std"), "all definers standard-born"], [css("--long"), "all definers long-born"],
       [css("--founder"), "founding clade"]]
    : [[css("--std"), "standard"], [css("--long"), "long"], [css("--mixed"), "mixed"],
       [css("--founder"), "founding clade"]];
  return rows.map(([c, t]) => `<span><i style="background:${c}"></i>${t}</span>`).join("");
}

async function renderTrees(main) {
  const sc = SCENARIOS.includes(settings.trScen) ? settings.trScen : SCENARIOS[0];
  const runs = ((await apiGet("/api/runs", {slot: "runs", label: "simulation index"}))
    .scenarios || {})[sc] || [];
  const cols = +settings.trCols || 10, k = settings.k, mt = settings.mt;
  main.classList.remove("fit");
  /* every simulation first, under the loading note, then the grid in one go:
     a walk left running after the tab is closed would keep fetching into
     whatever tab is open next. The answers stay cached, so a revisit refetches
     nothing it already has. */
  for (let i = 0; i < runs.length; i++) {
    if (activeTabId !== "trees") return;
    main.innerHTML = `<div class="loading">loading ${SCENARIO_LABEL[sc]}: simulation ${i + 1} of ${runs.length}…</div>`;
    try { await loadRun(sc, i); } catch (e) { if (isSuperseded(e)) return; }
  }
  main.innerHTML = `<div class="plot"><div class="ph">
      <h3>${cap1(SCENARIO_LABEL[sc])}: all ${runs.length} simulations</h3>
      <span class="sub">k = ${k} · m_T = ${mtPct(mt)}% · x: ${settings.treeAxis === "time" ? "time" : "divergence"}</span></div>
    <div class="trgrid" style="grid-template-columns:repeat(${cols},minmax(0,1fr))">${runs.map((r, i) =>
      `<div class="trcell" data-seed="${i}" title="${sdEsc(r.seed)}">${settings.trTags
        ? `<b>${String(i).padStart(2, "0")}</b>` : ""}<canvas></canvas></div>`).join("")}</div>
    <div class="key trkey">${trLegend()}${settings.trBurnin && settings.treeAxis === "time"
      ? `<span><i class="c4bi">${[0, 8, 16].map(x => `<i style="left:${x}px"></i>`).join("")}</i>burn-in day</span>` : ""}</div>
    <div class="cap">Each cell is one simulation of ${SCENARIO_LABEL[sc]}, its phylogeny
      coloured as on <a href="#clustering">Clustering</a>, at k = ${k} and
      m<sub>T</sub> = ${mtPct(mt)}%. Every cell has its own x range, the
      simulation's last day or deepest lineage, so branch lengths compare within
      a cell, not across cells.</div></div>`;
  const grid = main.querySelector(".trgrid");
  const h = Math.max(30, Math.round((grid.clientWidth / cols - 6) * TR_ASPECT));
  for (let i = 0; i < runs.length; i++) {
    const cell = grid.querySelector(`.trcell[data-seed="${i}"]`);
    if (!cell || cell.isConnected === false) return;      // re-rendered under us
    try {
      const v = await runView(sc, i, k, mt);
      const {colMap} = colourMapFor(v.T, v.C);
      const {g, w, h: hh} = fit(cell.querySelector("canvas"), h);
      const {xmax, sx} = paintTree(g, w, hh, v.D, v.T, v.C, colMap, 0.4);
      const day = burninDayOf(v.D);
      if (settings.trBurnin && day !== null && settings.treeAxis === "time" && day <= xmax)
        drawBurnin(g, sx(day), w, hh, null, true);
    } catch (e) {
      if (isSuperseded(e)) return;
      cell.classList.add("empty");
      cell.title = e.message;
    }
  }
}

registerTab({
  id: "trees", group: "supp", num: "S1", title: "All trees", viewId: "viewTrees",
  lede: "Every simulation of one scenario as a phylogeny, coloured by the clade definition set on Clustering.",
  rail: treesRail,
  render: renderTrees,
});
