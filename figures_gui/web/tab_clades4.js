/* Clustering panel: the Clustering tab's phylogeny and clade frequency, once
   per scenario with long shedders, all at the same seed. Seed index i is
   seed_000i in every scenario (/api/runs), so the four cells share the random
   seed and differ only by scenario. The clade definition, colouring and tree
   axis are the Clustering tab's own settings. */

const C4_LETTERS = ["A", "B", "C", "D"];

/* Two sliders on one seed: every seed, and a shortlist the user picks, so a few
   candidates can be stepped through without the ones in between. The shortlist
   is kept in settings without notifying (editing it redraws nothing). */
const C4_SEED_MAX = 49;
const c4Pad = v => String(v).padStart(4, "0");
const c4List = () => [...new Set(String(settings.f2Seeds || "").split(/[^0-9]+/).filter(Boolean)
  .map(Number).filter(v => v >= 0 && v <= C4_SEED_MAX))].sort((a, b) => a - b);
function c4SeedControl() {
  const el = document.createElement("div");
  el.className = "ctl c4seeds";
  el.innerHTML = `
    <div class="lab"><label>Seed</label><span class="val" data-r="all"></span></div>
    <input type="range" data-s="all" min="0" max="${C4_SEED_MAX}" step="1">
    <div class="lab" style="margin-top:8px"><label>Shortlist</label><span class="val" data-r="list"></span></div>
    <input type="range" data-s="list" min="0" step="1">
    <div class="c4srow"><input type="text" data-t="list" placeholder="e.g. 3, 7, 15" spellcheck="false">
      <button type="button" class="sw-reset" data-b="tog"></button></div>
    <div class="hint">Pick candidates on the first slider and add them, or type them in; the second
      slider then steps through those alone. Seed is shared with Clustering.</div>`;
  const q = s => el.querySelector(s), all = q('[data-s="all"]'), lst = q('[data-s="list"]');
  const txt = q('[data-t="list"]'), tog = q('[data-b="tog"]');
  const save = L => { settings.f2Seeds = L.join(", "); mirrorSettings(); };
  const sync = () => {
    const L = c4List(), s = +settings.seed, i = L.indexOf(s);
    all.value = s; q('[data-r="all"]').textContent = c4Pad(s);
    lst.max = Math.max(0, L.length - 1); lst.disabled = L.length < 2;
    if (i >= 0) lst.value = i;
    q('[data-r="list"]').textContent = !L.length ? "empty" : i >= 0 ? `${i + 1} of ${L.length}` : `— of ${L.length}`;
    if (document.activeElement !== txt) txt.value = L.join(", ");
    tog.textContent = i >= 0 ? `Remove ${c4Pad(s)}` : `Add ${c4Pad(s)}`;
  };
  all.addEventListener("input", () => { setSetting("seed", +all.value); sync(); });
  lst.addEventListener("input", () => { const L = c4List(); if (L.length) setSetting("seed", L[+lst.value]); sync(); });
  txt.addEventListener("change", () => { settings.f2Seeds = txt.value; save(c4List()); sync(); });
  tog.addEventListener("click", () => {
    const L = c4List(), s = +settings.seed;
    save(L.includes(s) ? L.filter(v => v !== s) : [...L, s].sort((a, b) => a - b)); sync(); });
  onSettingChange(n => { if (n === "seed" && el.isConnected !== false) sync(); });
  sync();
  return el;
}

function clades4Rail(rail) {
  railCompact(rail);
  rail.append(card("Display", [
    c4SeedControl(),
    checkControl({setting: "f2Titles", label: "Titles", hint: "Off: only the panel letters."}),
    checkControl({setting: "f2Dots", label: "Colour behind letters, panels"}),
    checkControl({setting: "f2DotsLeg", label: "Colour behind letters, legend"}),
    checkControl({setting: "f2TextCol", label: "Scenario names in colour"}),
    segmentControl({setting: "colourBy", label: "Colour by",
      options: [["origin", "Origin"], ["clade", "Clade"]],
      hint: "Origin compares across the four: clade colours are per-simulation "
          + "and mean nothing from one cell to the next."}),
    segmentControl({setting: "treeAxis", label: "Tree x axis",
      options: [["time", "Time"], ["div", "Divergence"]]}),
    paletteControl(),
  ]));
  rail.append(dataCard(DATA_CLADE));
}

/* one tick per simulated year under a frequency plot, placed as paintMuller places days */
function c4YearTicks(el, times) {
  const nt = times.length, last = times[nt - 1];
  const at = d => { const j = dayToIndex(times, d); return j === null || nt < 2 ? null : j / (nt - 1); };
  let html = "";
  for (let d = 0; d <= last + 0.5; d += 365) {
    const f = at(d); if (f === null) continue;
    const pos = `calc(2px + ${(100 * f).toFixed(3)}% - ${(4 * f).toFixed(3)}px)`;
    const al = f < 0.02 ? "translateX(0)" : f > 0.98 ? "translateX(-100%)" : "translateX(-50%)";
    html += `<i style="left:${pos}"></i><span style="left:${pos};transform:${al}">${d.toLocaleString()}</span>`;
  }
  el.innerHTML = html + `<em>days</em>`;
}

async function renderClades4(main) {
  const seed = settings.seed, titles = !!settings.f2Titles;
  /* all four first, one at a time (the "run" slot supersedes a concurrent
     load), then the grid in one go: a grid drawn early fills in cell by cell,
     and a second visit mid-load refetched the cell still in flight */
  const views = {};
  for (const sc of SUB_SCENARIOS) {
    try { views[sc] = await runView(sc, seed, settings.k, settings.mt); }
    catch (e) { if (isSuperseded(e)) return; views[sc] = e; }
  }
  main.classList.add("fit");
  /* the panel letter, on a disc in the scenario's colour unless switched off */
  const let_ = (sc, i, on) => `<b class="c4let${on ? " on" : ""}"${on
    ? ` style="background:${scenarioRaw(sc)}"` : ""}>${C4_LETTERS[i]}</b>`;
  const name = (sc, t) => settings.f2TextCol ? `<span style="color:${scenarioRaw(sc)}">${t}</span>` : t;
  const cap1 = t => t.charAt(0).toUpperCase() + t.slice(1);
  main.innerHTML = `<div class="plot">
    <div class="c4wrap"><div class="scgrid c4grid">${SUB_SCENARIOS.map((sc, i) => `<div>
      <h4>${let_(sc, i, settings.f2Dots)}${titles ? name(sc, cap1(SCENARIO_LABEL[sc])) : ""}</h4>
      ${views[sc] instanceof Error ? `<div class="err">${sdEsc(views[sc].message)}</div>` : ""}
      <canvas id="c4_${sc}_tree"></canvas>
      <canvas id="c4_${sc}_mul"></canvas>
      <div class="c4ticks" id="c4_${sc}_ticks"></div>
      </div>`).join("")}</div>
    <aside class="c4leg" id="c4_leg"></aside></div>
    <div class="cap">Each cell is one simulation: the scenario named above it, at
      seed ${seed}. Top, the phylogeny; below, lineage frequency stacked by
      clade, with the burn-in day marked on both. The seed is shared, but the
      scenarios differ in who sheds for how long, so the runs diverge from the
      first long shedder on — same seed is not the same epidemic. Each cell has
      its own x range, the simulation's own last day or deepest lineage.</div></div>`;

  /* the window's height, shared by the grid's rows: no scrolling */
  const cvs = SUB_SCENARIOS.flatMap(sc => [`c4_${sc}_tree`, `c4_${sc}_mul`].map(id => document.getElementById(id)));
  const grid = main.querySelector(".c4grid");
  const cols = grid && (getComputedStyle(grid).gridTemplateColumns || "").split(" ").filter(Boolean).length || 1;
  const per = fitFree(main, cvs) / Math.ceil(SUB_SCENARIOS.length / cols);
  const treeH = clampH(per * 0.64, 80, 420), mulH = clampH(per * 0.36, 50, 240);

  const counts = [];
  SUB_SCENARIOS.forEach((sc, i) => {
    const view = views[sc];
    if (view instanceof Error) return;
    const el = id => document.getElementById(`c4_${sc}_${id}`);
    const {D, T, C} = view;
    const {ordered, colMap} = colourMapFor(T, C);
    const {g, w, h} = fit(el("tree"), treeH);
    const {xmax, sx} = paintTree(g, w, h, D, T, C, colMap, 0.6);
    const day = burninDayOf(D);
    if (day !== null && settings.treeAxis === "time" && day <= xmax)
      drawBurnin(g, sx(day), w, h, null, false);
    const m = paintMuller(el("mul"), {D, C, ordered, colMap, height: mulH,
      washBurnin: false, selected: null, burninLabel: false});
    c4YearTicks(el("ticks"), D.times);
    counts.push(`<div class="r stack"><span class="c4head">${let_(sc, i, settings.f2DotsLeg)}<b>${name(sc, SCENARIO_LABEL[sc])}</b></span><span>${D.n.toLocaleString()} lineages →
      ${C.alive.size.toLocaleString()} clades${m.nClades > m.shown ? ` (${m.shown} largest drawn)` : ""}</span></div>`);
  });
  const colours = settings.colourBy === "origin"
    ? (settings.originRule === "frac"
        ? [[css("--std"), "all definers standard-born"], [css("--long"), "all definers long-born"], [css("--founder"), "founding clade"]]
        : [[css("--std"), "standard"], [css("--long"), "long"], [css("--mixed"), "mixed"], [css("--founder"), "founding clade"]])
      .map(([c, t]) => `<div class="r"><i class="c4sw" style="background:${c}"></i><span>${t}</span></div>`).join("")
    : `<div class="r"><span>one colour per clade, per simulation</span></div>`;
  document.getElementById("c4_leg").innerHTML = `
    <div class="r"><span><b>Simulation shown:</b> <b>${String(seed).padStart(4, "0")}</b></span></div>
    <div class="r"><i class="c4bi">${[0, 8, 16].map(x => `<i style="left:${x}px"></i>`).join("")}</i><span>burn-in day</span></div>
    <div class="r"><span>Clustering parameters: <span class="nw">k = ${settings.k}, m_T = ${mtPct(settings.mt)}%</span></span></div>
    <h5>Scenarios</h5>${counts.join("")}${settings.colourBy === "origin" ? `<div class="c4legfoot">${colours}</div>` : ""}`;
}

registerTab({
  id: "clades4", group: "figures", num: "2", title: "Clustering", viewId: "viewC4",
  lede: "The phylogeny and clade frequency of one seed in each scenario with long shedders, side by side.",
  rail: clades4Rail,
  render: renderClades4,
});
