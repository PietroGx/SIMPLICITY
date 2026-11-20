/* Clustering: one run at a time, as a phylogeny and as clade frequency.

   This is the artifact's tab 1. What it draws is the clustering's own decision
   -- the browser rebuilds clade membership from the `starts` the exporter took
   from `cluster_lin_into_clades_with_meta`, the same function the clade pass
   calls -- so the colours here and the numbers on the panel tabs mean the same
   thing at the same k.

   The element ids are the artifact's (`cl_tree`, `cl_muller`, `cl_ovGrid`),
   because studio.css carries its rules scoped to them.
*/

const CLUSTER_TREE_H = 300;
const CLUSTER_MULLER_H = 190;
const OVERVIEW_CELL_H = 58;

function clusteringRail(rail) {
  rail.append(railSimulation(50));
  rail.append(railClade());
  rail.append(railOrigin());
  rail.append(card("Display", [
    segmentControl({setting: "colourBy", label: "Colour by",
      options: [["clade", "Clade"], ["origin", "Origin"]],
      hint: "Colour by <em>Origin</em> to compare against the panels — clade "
          + "colours are per-simulation and mean nothing across simulations."}),
    segmentControl({setting: "treeAxis", label: "Tree x axis",
      options: [["time", "Time"], ["div", "Divergence"]],
      hint: "Divergence is the lineage's mutation count, so a long flat branch "
          + "in time can be a short one here."}),
    checkControl({setting: "overview", label: "Overview: every simulation as a small tree",
      hint: "All 50 simulations of this scenario, drawn with the settings above."}),
    buttonControl({label: "Reset to defaults", onClick: () =>
      resetSettings(["k", "mt", "originRule", "colourBy", "treeAxis", "burnin"])}),
  ]));
}

function clusteringStats(view) {
  const {D, T, C} = view;
  const sizes = [...C.alive].map(r => C.info.get(r).size);
  const singles = sizes.filter(n => n <= 1).length;
  const cells = [
    ["lineages", D.n.toLocaleString()],
    ["tips", T.nTips.toLocaleString()],
    ["clades", C.alive.size.toLocaleString()],
    ["singletons", `${(100 * singles / Math.max(sizes.length, 1)).toFixed(0)}<small>%</small>`],
    ["burn-in day", burninDayOf(D) === null ? "—" : fmt(burninDayOf(D), 0)],
  ];
  return `<div class="stats">` + cells.map(([k, v]) =>
    `<div class="stat"><div class="k">${k}</div><div class="v">${v}</div></div>`).join("")
    + `</div>`;
}

function cladeKeyRow(view) {
  if (settings.colourBy === "clade")
    return `<div class="cap">One colour per clade, ordered down the tree. Clade
      colours are per-simulation: the same hue in two simulations is not the same clade.</div>`;
  if (settings.originRule === "frac")
    return keyRow([[css("--std"), "all definers standard-born"],
                   [css("--long"), "all definers long-born"],
                   [css("--founder"), "founding clade"]])
      + `<div class="cap">A blend: the hue is the long-born share of the clade's
         defining mutations.</div>`;
  return keyRow([[css("--std"), "standard"], [css("--long"), "long"],
                 [css("--mixed"), "mixed"], [css("--founder"), "founding clade"]]);
}

async function renderClustering(main) {
  const view = await runView(settings.scenario, settings.seed, settings.k, settings.mt);
  const {D, T, C} = view;
  const {ordered, colMap} = colourMapFor(T, C);

  if (settings.overview) {
    main.innerHTML = clusteringStats(view) + `
      <div class="plot"><div class="ph"><h3>Every simulation of ${SCENARIO_LABEL[settings.scenario]}</h3>
        <span class="sub">k = ${settings.k} · m_T = ${mtPct(settings.mt)}%</span></div>
        <div class="ovgrid" id="cl_ovGrid"></div>
        <div class="cap">Each cell is one simulation's phylogeny at these settings.
          Click a cell to open that simulation on its own.</div></div>`;
    await paintOverviewGrid(document.getElementById("cl_ovGrid"));
    return;
  }

  main.innerHTML = clusteringStats(view) + `
    <div class="plot"><div class="ph">
      <h3>Phylogeny, coloured by ${settings.colourBy === "origin" ? "origin" : "clade"}</h3>
      <span class="sub" id="cl_treeSub"></span></div>
      <canvas id="cl_tree"></canvas>
      <div class="axis"><span id="cl_treeX0"></span><span id="cl_treeX1"></span></div></div>
    <div class="plot"><div class="ph"><h3>Lineage frequency by clade</h3>
      <span class="sub" id="cl_mullerSub"></span></div>
      <canvas id="cl_muller"></canvas>
      <div class="axis"><span id="cl_mX0"></span><span id="cl_mX1"></span></div>
      ${cladeKeyRow(view)}</div>`;

  const treeCv = document.getElementById("cl_tree");
  const {g, w, h} = fit(treeCv, CLUSTER_TREE_H);
  const {xmax, sx} = paintTree(g, w, h, D, T, C, colMap);
  const day = burninDayOf(D);
  if (day !== null && settings.treeAxis === "time" && day <= xmax)
    drawBurnin(g, sx(day), w, h, day, false);

  document.getElementById("cl_treeSub").textContent =
    `${D.n.toLocaleString()} lineages · ${T.nTips.toLocaleString()} tips`;
  document.getElementById("cl_treeX0").textContent =
    settings.treeAxis === "time" ? "day 0" : "0 mutations";
  document.getElementById("cl_treeX1").textContent =
    settings.treeAxis === "time" ? `day ${xmax}` : `${xmax} mutations`;

  const mcv = document.getElementById("cl_muller");
  const m = paintMuller(mcv,
    {D, C, ordered, colMap, height: CLUSTER_MULLER_H, washBurnin: false,
     selected: settings.clade});
  /* The same selection as the Sweep tab: one run, one partition, so the id
     means the same thing on both. */
  mcv.style.cursor = "pointer";
  mcv.onclick = ev => {
    const b = mcv.getBoundingClientRect();
    const id = m.hit(ev.clientX - b.left, ev.clientY - b.top);
    setSetting("clade", id !== null && id === settings.clade ? null : id);
  };
  document.getElementById("cl_mullerSub").textContent = (m.nClades > m.shown
    ? `${m.shown} largest of ${m.nClades.toLocaleString()} clades drawn · rest in grey`
    : `${m.nClades.toLocaleString()} clades`)
    + (m.selected === null ? "" : ` · showing ${labelOf(C.info.get(m.selected))}`);
  document.getElementById("cl_mX0").textContent = `day ${m.day0}`;
  document.getElementById("cl_mX1").textContent = `day ${m.day1}`;
}

/* The overview grid. Runs are loaded one at a time and painted as each arrives,
   so the grid fills in rather than blocking on fifty fetches. */
async function paintOverviewGrid(host) {
  const runs = (await apiGet("/api/runs", {slot: "runs", label: "simulation index"}))
    .scenarios[settings.scenario] || [];
  host.innerHTML = runs.map((r, i) =>
    `<div class="ovcell pending" data-seed="${i}" title="${r.seed}">
       <canvas></canvas></div>`).join("");
  host.querySelectorAll(".ovcell").forEach(cell => cell.addEventListener("click", () => {
    setSettings({seed: +cell.dataset.seed, overview: 0});
  }));
  for (let i = 0; i < runs.length; i++) {
    const cell = host.querySelector(`.ovcell[data-seed="${i}"]`);
    if (!cell || !document.body.contains(cell)) return;   // tab changed under us
    try {
      const v = await runView(settings.scenario, i, settings.k, settings.mt);
      const {colMap} = colourMapFor(v.T, v.C);
      const cv = cell.querySelector("canvas");
      const {g, w, h} = fit(cv, OVERVIEW_CELL_H);
      paintTree(g, w, h, v.D, v.T, v.C, colMap, 0.4);
      cell.classList.remove("pending");
    } catch (e) {
      if (isSuperseded(e)) return;
      cell.classList.remove("pending");
      cell.classList.add("empty");
    }
  }
}

registerTab({
  id: "clustering", group: "explore", title: "Clustering", viewId: "viewCl",
  lede: "One simulation as a phylogeny and as clade frequency. The clade definition "
      + "set here is what every panel downstream reads.",
  rail: clusteringRail,
  render: renderClustering,
});
