/* Clustering panel: the Clustering tab's phylogeny and clade frequency, once
   per scenario with long shedders, all at the same seed. Seed index i is
   seed_000i in every scenario (/api/runs), so the four cells share the random
   seed and differ only by scenario. The clade definition, colouring and tree
   axis are the Clustering tab's own settings. */

const C4_TREE_H = 220;
const C4_MULLER_H = 120;

function clades4Rail(rail) {
  railCompact(rail);
  rail.append(card("Display", [
    segmentControl({setting: "colourBy", label: "Colour by",
      options: [["origin", "Origin"], ["clade", "Clade"]],
      hint: "Origin compares across the four: clade colours are per-simulation "
          + "and mean nothing from one cell to the next."}),
    segmentControl({setting: "treeAxis", label: "Tree x axis",
      options: [["time", "Time"], ["div", "Divergence"]]}),
    paletteControl(),
  ]));
  rail.append(dataCard(`The seed (the same in all four scenarios) is the one picked on
    <a href="#clustering">Clustering</a>. ${DATA_CLADE}`));
}

async function renderClades4(main) {
  const seed = settings.seed;
  /* all four first, one at a time (the "run" slot supersedes a concurrent
     load), then the grid in one go: a grid drawn early fills in cell by cell,
     and a second visit mid-load refetched the cell still in flight */
  const views = {};
  for (const sc of SUB_SCENARIOS) {
    try { views[sc] = await runView(sc, seed, settings.k, settings.mt); }
    catch (e) { if (isSuperseded(e)) return; views[sc] = e; }
  }
  main.innerHTML = `<div class="plot"><div class="ph">
      <h3>Phylogeny and clade frequency, by scenario</h3>
      <span class="sub">seed_${String(seed).padStart(4, "0")} · k = ${settings.k} ·
        m_T = ${mtPct(settings.mt)}% · coloured by ${settings.colourBy === "origin" ? "origin" : "clade"}</span></div>
    <div class="scgrid">${SUB_SCENARIOS.map(sc => `<div>${scenarioHead(sc)}
      <div class="c4sub" id="c4_${sc}_sub">${views[sc] instanceof Error ? views[sc].message : ""}</div>
      <canvas id="c4_${sc}_tree"></canvas>
      <div class="axis"><span id="c4_${sc}_tx0"></span><span id="c4_${sc}_tx1"></span></div>
      <canvas id="c4_${sc}_mul"></canvas>
      <div class="axis"><span id="c4_${sc}_mx0"></span><span id="c4_${sc}_mx1"></span></div>
      </div>`).join("")}</div>
    ${settings.colourBy === "origin" ? cladeKeyRow({}) : ""}
    <div class="cap">Each cell is one simulation: the scenario named above it, at
      seed ${seed}. Top, the phylogeny; below, lineage frequency stacked by
      clade, with the burn-in day marked on the tree. The seed is shared, but the
      scenarios differ in who sheds for how long, so the runs diverge from the
      first long shedder on — same seed is not the same epidemic. Each cell has
      its own x range, the simulation's own last day or deepest lineage.</div></div>`;

  SUB_SCENARIOS.forEach(sc => {
    const view = views[sc];
    if (view instanceof Error) return;
    const el = id => document.getElementById(`c4_${sc}_${id}`);
    const {D, T, C} = view;
    const {ordered, colMap} = colourMapFor(T, C);
    const {g, w, h} = fit(el("tree"), C4_TREE_H);
    const {xmax, sx} = paintTree(g, w, h, D, T, C, colMap, 0.6);
    const day = burninDayOf(D);
    if (day !== null && settings.treeAxis === "time" && day <= xmax)
      drawBurnin(g, sx(day), w, h, day, false);
    const m = paintMuller(el("mul"), {D, C, ordered, colMap, height: C4_MULLER_H,
      washBurnin: false, selected: null});
    el("sub").textContent = `${D.n.toLocaleString()} lineages · ${C.alive.size.toLocaleString()} clades`
      + (m.nClades > m.shown ? ` · ${m.shown} largest drawn` : "");
    el("tx0").textContent = settings.treeAxis === "time" ? "day 0" : "0 mutations";
    el("tx1").textContent = settings.treeAxis === "time" ? `day ${xmax}` : `${xmax} mutations`;
    el("mx0").textContent = `day ${m.day0}`;
    el("mx1").textContent = `day ${m.day1}`;
  });
}

registerTab({
  id: "clades4", group: "figures", num: "2", title: "Clustering", viewId: "viewC4",
  lede: "The phylogeny and clade frequency of one seed in each scenario with long shedders, side by side.",
  rail: clades4Rail,
  render: renderClades4,
});
