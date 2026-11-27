/* S2 · Clustering sensitivity and S3 · Outbreak detection sensitivity.

   Each metric is the one its figure shows, computed by that figure's own code
   at each level, one parameter moved at a time around the settings in force:

     Figure 2  /api/across (the clade builds): clades per simulation, median clade size
     Figure 4  /api/panel?n=1h (the panel Figure 4 draws): frequency AUC and
               infections in clades with long-born mutations
     Figure 5  /api/explore_panel under the detection tabs' rules, then
               p3Breakouts -- Figure 5's own per-clade pass -- with the level
               patched over the settings (swWithSettings)

   Figure 3 reads the transmission table, which no clade definition touches,
   and Figure 1 involves no clustering, so neither has a clustering axis.
   Bands are the bootstrap 95% interval over simulations, the figures' own
   estimator (p3BootBy, panels._boot). */

const SA_F2 = [["clades", "Clades per simulation", "count"],
               ["median_size", "Median clade size (lineages)", "count"]];
const SA_F4 = [["any_long", "Frequency AUC, clades with long-born mutations (%)", "pct"],
               ["only_long", "Frequency AUC, clades with only long-born mutations (%)", "pct"],
               ["inf_long", "Infections, clades with long-born mutations (%)", "pct"]];
const SA_F5 = [["brk", "Outbreaks per simulation", "count"],
               ["cry", "Cryptic outbreaks per simulation", "count"],
               ["cryLong", "Cryptic outbreaks per simulation, clades with long-born mutations", "count"],
               ["share", "Cryptic share of outbreaks (%)", "pct"],
               ["ih", "Median intra-host hidden time, cryptic outbreaks (d)", "days"],
               ["bh", "Median between-host hidden time, cryptic outbreaks (d)", "days"]];

const saMedian = v => {
  if (!v.length) return null;
  const s = [...v].sort((a, b) => a - b), m = s.length >> 1;
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
};
let saToken = 0;                      // a newer render abandons an older one mid-fetch
const saGet = (url, label) => apiGet(url, {label}).catch(e => ({error: e.message}));

/* Figure 5's metrics per scenario from one /api/explore_panel answer, under
   `patch` over the settings in force. p3Source is pinned to the detection
   tabs' rules, the source the answer came from. */
function saFig5(d, patch) {
  const B = swWithSettings(Object.assign({}, patch, {p3Source: "sweep"}), () => p3Breakouts(d));
  const out = {};
  SCENARIOS.forEach(sc => {
    const b = B[sc];
    if (!b || !b.nRuns) { out[sc] = {}; return; }
    const n = b.nRuns, zero = () => new Array(n).fill(0), lists = () => Array.from({length: n}, () => []);
    const brk = zero(), cry = zero(), cryL = zero(), ih = lists(), bh = lists();
    b.cl.forEach(c => {
      if (c.ev) brk[c.slot]++;
      if (!c.cry) return;
      cry[c.slot]++;
      if (c.cls !== "standard" && c.cls !== "founding") cryL[c.slot]++;
      if (isNum(c.ih)) ih[c.slot].push(c.ih);
      if (isNum(c.bh)) bh[c.slot].push(c.bh);
    });
    const sum = (a, ix) => ix.reduce((t, i) => t + a[i], 0);
    out[sc] = {brk: p3Boot(brk), cry: p3Boot(cry), cryLong: p3Boot(cryL),
      share: p3BootBy(n, ix => { const t = sum(brk, ix); return t ? 100 * sum(cry, ix) / t : null; }),
      ih: p3BootBy(n, ix => saMedian(ix.flatMap(i => ih[i]))),
      bh: p3BootBy(n, ix => saMedian(ix.flatMap(i => bh[i])))};
  });
  return out;
}

/* Figure 4's metrics per scenario from /api/panel?n=1h, as percentages */
function saFig4(d) {
  const out = {};
  SCENARIOS.forEach(sc => {
    const s = (d.scenarios || {})[sc];
    out[sc] = !s ? {} : {any_long: asPct((s.burden || {}).any_long),
                         only_long: asPct((s.burden || {}).only_long),
                         inf_long: asPct((s.infections || {}).any_long)};
  });
  return out;
}

/* Figure 2's metrics per scenario at m_T = m from one /api/across answer */
function saFig2(d, m) {
  const out = {};
  SCENARIOS.forEach(sc => {
    const seeds = ((d.scenarios || {})[sc] || {}).seeds || [];
    const at = key => seeds.map(s => ((s.per_value || []).find(x => Math.abs(x.m - m) < 1e-12) || {})[key]);
    out[sc] = seeds.length ? {clades: p3Boot(at("clades")), median_size: p3Boot(at("median_size"))} : {};
  });
  return out;
}

/* One block of facets: a figure's metrics over the levels, scenarios as lines.
   `per[i]` is level i's {scenario: {metric: band}}, or an Error. */
function saBlock(idp, title, metrics, levels) {
  return `<div class="safig">${title}</div><div class="sagrid">${metrics.map(([key, label], j) =>
    `<div><h4>${label}</h4><canvas id="${idp}_${key}"></canvas></div>`).join("")}</div>`;
}
function saDraw(idp, metrics, levels, per, {mark, xlab}) {
  metrics.forEach(([key, label, kind]) => {
    const cv = document.getElementById(`${idp}_${key}`);
    if (!cv) return;
    linesOverLevels(cv, {levels, mark, ribbon: true, height: 170, xlab,
      zero: kind !== "days", pct: kind === "pct" ? "fit" : false, ylab: "",
      fmtY: v => fmt(v, kind === "days" ? 0 : Math.abs(v) < 10 ? 1 : 0),
      series: SCENARIOS.map(sc => ({colour: scenarioColour(sc), label: SCENARIO_LABEL[sc],
        vals: per.map(p => p instanceof Error ? null : ((p || {})[sc] || {})[key] || null)}))});
  });
}
const saMissing = (levels, per) => levels.filter((_, i) => per[i] instanceof Error)
  .map((lv, i) => lv).join(", ");
const saErr = d => d && d.error ? new Error(d.error) : null;

/* ------------------------------------------------- S2: clustering sensitivity */

function saClustRail(rail) {
  rail.append(railFocus());
  rail.append(card("Sensitivity", [
    segmentControl({setting: "saVary", label: "Vary",
      options: [["mt", "m_T"], ["k", "k"]],
      hint: "m<sub>T</sub> walks the grid's m<sub>T</sub> values at the k in force; k walks "
          + "the grid's k values at the m<sub>T</sub> in force. Everything else stays as set."}),
  ]));
  rail.append(dataCard(`${DATA_CLADE} Figure 5's detection rules are set on
    <a href="#sweep">Outbreak detection</a> and <a href="#cryptic">Cryptic outbreaks</a>.`));
}

async function renderSaClust(main) {
  const tok = ++saToken, vary = settings.saVary === "k" ? "k" : "mt";
  const ks = (AXES.k || [5, 7, 10]).map(Number), mts = (AXES.mt || []).map(Number);
  const pairs = vary === "mt" ? mts.map(m => ({k: +settings.k, mt: m}))
                              : ks.map(k => ({k, mt: +settings.mt}));
  const levels = pairs.map(p => vary === "mt" ? mtPct(p.mt) + "%" : String(p.k));
  const mark = vary === "mt" ? mtPct(settings.mt) + "%" : String(settings.k);
  const xlab = vary === "mt" ? "minimum transmissions m_T" : "shared-mutation threshold k";
  main.innerHTML = `<div class="loading">reading ${pairs.length} clade builds…</div>`;

  /* Figures 2 and 4 read the clade builds, fast; Figure 5 is measured level by
     level and drawn as each arrives. A level whose runs are not warmed on disk
     costs minutes (warm_explore.py --pairs all fills them), so the walk stops
     when the tab is left and resumes from the cache on the next visit. */
  const across = {};
  for (const k of new Set(pairs.map(p => p.k)))
    across[k] = await saGet("/api/across?" + buildQuery({k}), "Figure 2 over the clade builds");
  const f4 = await Promise.all(pairs.map(p =>
    saGet("/api/panel?" + buildQuery({n: "1h", k: p.k, mt: p.mt}), "Figure 4 at each clade build")));
  if (tok !== saToken) return;
  const per2 = pairs.map(p => saErr(across[p.k]) || saFig2(across[p.k], p.mt));
  const per4 = f4.map(d => saErr(d) || saFig4(d));
  const per5 = pairs.map(() => null);
  const pinned = vary === "mt" ? `k = ${settings.k}` : `m<sub>T</sub> = ${mtPct(settings.mt)}%`;
  main.innerHTML = `<div class="plot"><div class="ph"><h3>Over ${vary === "mt" ? "m<sub>T</sub>" : "k"}</h3>
      <span class="sub">${pinned} · burn-in ${settings.burnin ? "excluded" : "included"}</span></div>
    ${saBlock("sac2", "Figure 2 · Clustering", SA_F2, levels)}
    ${saBlock("sac4", "Figure 4 · Frequency AUC", SA_F4, levels)}
    ${saBlock("sac5", "Figure 5 · Outbreaks", SA_F5, levels)}
    ${keyRow(SCENARIOS.map(s => [scenarioColour(s), SCENARIO_LABEL[s]]))}
    <div class="note" id="sac_status"></div><div id="sac_miss"></div>
    <div class="cap">Each line is one scenario; the band is the bootstrap 95% interval
      of the mean (of the median, for hidden times) over simulations, the dashed
      rule the setting in force. Every point is its figure's own measurement
      recomputed at that clade definition: the clade builds for Figures 2 and 4,
      Figure 5's pass under the detection rules set on the exploration tabs.
      Figure 3 is left out because no clade definition enters it.</div></div>`;
  const o = {mark, xlab};
  saDraw("sac2", SA_F2, levels, per2, o);
  saDraw("sac4", SA_F4, levels, per4, o);
  saDraw("sac5", SA_F5, levels, per5, o);
  const say = (id, html) => { const el = document.getElementById(id); if (el) el.innerHTML = html; };
  const showMissing = () => {
    const miss = [["Figure 2", per2], ["Figure 4", per4], ["Figure 5", per5]]
      .map(([f, per]) => [f, saMissing(levels, per)]).filter(([, m]) => m);
    say("sac_miss", miss.length ? `<div class="err">Not available at every level — ${miss.map(([f, m]) =>
      `${f}: ${m}`).join("; ")}. The clade builds and the warmed runs come from the HPC
      (<code>figures_gui/prepare.py --slurm</code>).</div>` : "");
  };
  showMissing();
  for (let i = 0; i < pairs.length; i++) {
    if (tok !== saToken || activeTabId !== "saclust") return;
    say("sac_status", `Figure 5: measuring ${levels[i]}, ${i} of ${pairs.length} levels done…`);
    const p = pairs[i];
    const d = await saGet("/api/explore_panel?" + exploreQuery({k: p.k, mt: p.mt, mf: 0}),
                          `Figure 5 at k = ${p.k}, m_T = ${mtPct(p.mt)}%`);
    if (tok !== saToken) return;
    per5[i] = saErr(d) || saFig5(d, {k: p.k, mt: p.mt});
    saDraw("sac5", SA_F5, levels, per5, o);
  }
  say("sac_status", "");
  showMissing();
}

/* -------------------------------------------- S3: outbreak detection sensitivity */

/* axis -> [label, level format, patch for a level, the level in force] */
const SA_OS = {
  brkF: ["outbreak peak F", v => brkPct(v), v => ({brkF: v}), () => +settings.brkF],
  brkH: ["cryptic hidden time H", v => v + " d", v => ({brkH: v}), () => +settings.brkH],
  brkLo: ["ceiling f", v => brkPct(v), v => ({brkLo: v}), () => +settings.brkLo],
  sdWin: ["smoothing window", v => v + " d", v => ({sdWin: v}), () => +settings.sdWin],
  sdProm: ["peak prominence", v => fmt(v, 2), v => ({sdProm: v}), () => +settings.sdProm],
};
const SA_OS_HINT = {
  brkF: "An outbreak is a sweep peaking at or above F.",
  brkH: "An outbreak is cryptic if it stayed hidden at least H days from first appearance.",
  brkLo: "The cryptic ceiling: hidden time, and Figure 5's between-host hidden time, end the first day above f.",
  sdWin: "The moving average the frequency is read through before peaks are found.",
  sdProm: "How far a peak must stand above its dip, as a share of the reference set on Outbreak detection.",
};

function saOsRail(rail) {
  rail.append(railFocus());
  rail.append(card("Sensitivity", [
    pickControl({setting: "osAxis", label: "Rule to vary",
      options: Object.entries(SA_OS).map(([k, [label]]) => [k, label, SA_OS_HINT[k]]),
      hint: "One rule at a time; every other rule stays as set on the exploration tabs."}),
  ]));
  rail.append(dataCard(`The detection and cryptic rules are set on <a href="#sweep">Outbreak
    detection</a> and <a href="#cryptic">Cryptic outbreaks</a>; ${DATA_CLADE}`));
}

async function renderSaOs(main) {
  const tok = ++saToken, ax = SA_OS[settings.osAxis] ? settings.osAxis : "brkF";
  const [label, fmtL, patchOf, now] = SA_OS[ax];
  const base = ((SA_LEVELS.outbreak || {})[ax] || []).slice();
  const cur = now();
  const vals = [...new Set([...base.map(Number), cur])].sort((a, b) => a - b);
  const levels = vals.map(fmtL), mark = fmtL(cur);
  main.innerHTML = `<div class="loading">measuring Figure 5 at ${vals.length} levels…</div>`;
  /* the browser-side rules (F and H) reuse one answer: the url is the same */
  const urls = vals.map(v => "/api/explore_panel?" + swWithSettings(patchOf(v), () => exploreQuery({mf: 0})));
  const got = {};
  for (const u of new Set(urls)) {
    if (tok !== saToken) return;
    got[u] = await saGet(u, `Figure 5, ${label}`);
  }
  if (tok !== saToken) return;
  const per = vals.map((v, i) => saErr(got[urls[i]]) || saFig5(got[urls[i]], patchOf(v)));
  const miss = saMissing(levels, per);
  main.innerHTML = `<div class="plot"><div class="ph"><h3>Over the ${label}</h3>
      <span class="sub">k = ${settings.k} · m<sub>T</sub> = ${mtPct(settings.mt)}% · ${new Set(urls).size} measurement${new Set(urls).size === 1 ? "" : "s"}</span></div>
    ${saBlock("sao5", "Figure 5 · Outbreaks", SA_F5, levels)}
    ${keyRow(SCENARIOS.map(s => [scenarioColour(s), SCENARIO_LABEL[s]]))}
    ${miss ? `<div class="err">Not measured at ${miss}.</div>` : ""}
    <div class="cap">Each line is one scenario; the band is the bootstrap 95% interval
      of the mean (of the median, for hidden times) over simulations, the dashed rule
      the setting in force. Every point is Figure 5's own pass with that one rule
      changed and every other as set on the exploration tabs. ${SA_OS_HINT[ax]}</div></div>`;
  saDraw("sao5", SA_F5, levels, per, {mark, xlab: label});
}

registerTab({
  id: "saclust", group: "supp", num: "S2", title: "Sensitivity to clustering", viewId: "viewSaClust",
  lede: "How Figures 2, 4 and 5 move with the clade definition: k or m<sub>T</sub>, one at a time.",
  rail: saClustRail,
  render: renderSaClust,
});
registerTab({
  id: "saout", group: "supp", num: "S3", title: "Sensitivity to outbreak detection", viewId: "viewSaOut",
  lede: "How Figure 5 moves with the outbreak and cryptic rules, one rule at a time.",
  rail: saOsRail,
  render: renderSaOs,
});
