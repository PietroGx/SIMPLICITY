/* Panel 2, Transmission: per scenario with long shedders, the jump standard
   hosts and long shedders add to what they pass on (fig3 row 3's violin)
   beside where the transmissions go (the arrow graph). The matrix, stacked and intensity forms
   went in v2.5.19; intensity stays in the Summary table.

   Canvas ids are <idp>_* : shell.js builds every tab view into ONE document,
   so bare ids here and in another tab would be the same element.
*/

/* Where the virus goes, as a graph: a standard node and a long node, ss and ll
   as self-loops, sl over the top and ls underneath. Line width is sqrt(share)
   on one scale for every subpanel, so ~99% standard->standard does not reduce
   the other three to hairlines; the label carries the number. */
const P2_ARROW_W = 14;
function p2Arrows(cv, {cell, height, longCol}) {
  const {g, w, h} = fit(cv, height || 220);
  /* the loops reach ~2.8 r past their node and carry a label: keep both on canvas */
  const r = Math.max(12, Math.min(22, h * 0.1, w * 0.075));
  const m = 2.9 * r + 8;
  const cy = h / 2 + 4, xs = Math.max(w * 0.3, m), xl = Math.min(w * 0.7, w - m);
  const col = {s: css("--std"), l: longCol || css("--long")};
  const at = (cx, a, rr = r) => [cx + rr * Math.cos(a), cy + rr * Math.sin(a)];
  const lw = b => b && isNum(b.mean) ? 0.8 + P2_ARROW_W * Math.sqrt(Math.max(0, b.mean) / 100) : 0;
  /* Each share in two lines -- 3 significant digits, then ± half its 95% CI
     to 2 -- set a fixed gap outside its own arrow's extreme, so no two are
     lined up by construction. The full interval is the label's hover. */
  const hits = [], NAME = {ss: "standard → standard", ll: "long → long", sl: "standard → long", ls: "long → standard"};
  const sig = (v, k) => v === 0 ? "0" : v.toFixed(Math.max(0, k - 1 - Math.floor(Math.log10(Math.abs(v)))));
  const GAP = 5, LH = 12;
  const lab = (b, ext, below, which) => {
    if (!b || !isNum(b.mean) || !ext) return;
    const t1 = sig(b.mean, 3) + "%", t2 = isNum(b.lo) && isNum(b.hi) ? `± ${sig((b.hi - b.lo) / 2, 2)}` : "";
    const [x, y] = ext, n = t2 ? 2 : 1;
    const y1 = below ? y + GAP + 9 : y - GAP - (n - 1) * LH;
    g.textAlign = "center"; g.fillStyle = css("--ink"); g.font = "600 10.5px " + css("--f-num");
    g.fillText(t1, x, y1);
    if (t2) { g.font = "10.5px " + css("--f-num"); g.fillText(t2, x, y1 + LH); }
    hits.push({x, y: y1 + (n - 1) * LH / 2 - 4, text: {text: `${NAME[which]}: share of each simulation's transmission events, averaged over simulations.`,
      math: `share    ${sig(b.mean, 3)}%` + (t2 ? `\n95% CI   ${sig(b.lo, 3)} – ${sig(b.hi, 3)}%` : "")}});
  };
  /* the topmost and bottommost points of a drawn stroke, its half-width included */
  const extremes = (p0, c1, c2, e, W) => {
    let top = null, bot = null;
    for (let i = 0; i <= 40; i++) {
      const t = i / 40, u = 1 - t, k = [u * u * u, 3 * u * u * t, 3 * u * t * t, t * t * t];
      const x = k[0] * p0[0] + k[1] * c1[0] + k[2] * c2[0] + k[3] * e[0];
      const y = k[0] * p0[1] + k[1] * c1[1] + k[2] * c2[1] + k[3] * e[1];
      if (!top || y < top[1]) top = [x, y]; if (!bot || y > bot[1]) bot = [x, y];
    }
    return {top: [top[0], top[1] - W / 2], bot: [bot[0], bot[1] + W / 2]};
  };
  /* a cubic from p0 to p3; the stroke stops short of the tip so the head covers it */
  const edge = (c, b, p0, c1, c2, p3) => {
    const W = lw(b);
    const dx = p3[0] - c2[0], dy = p3[1] - c2[1], dl = Math.hypot(dx, dy) || 1;
    const ux = dx / dl, uy = dy / dl, hl = 5 + W * 0.9, hw = 3.5 + W * 0.75;
    const e = [p3[0] - ux * hl * 0.85, p3[1] - uy * hl * 0.85];
    const ext = extremes(p0, c1, c2, e, W);
    if (!W) return ext;
    g.strokeStyle = c; g.fillStyle = c; g.lineWidth = W; g.lineCap = "butt"; g.globalAlpha = 0.85;
    g.beginPath(); g.moveTo(...p0); g.bezierCurveTo(...c1, ...c2, ...e); g.stroke();
    g.beginPath(); g.moveTo(...p3);
    g.lineTo(p3[0] - ux * hl - uy * hw, p3[1] - uy * hl + ux * hw);
    g.lineTo(p3[0] - ux * hl + uy * hw, p3[1] - uy * hl - ux * hw);
    g.closePath(); g.fill(); g.globalAlpha = 1;
    return ext;
  };
  const loop = (cx, side, c, b, which) => {
    const a0 = side < 0 ? Math.PI : 0, s = 0.6;
    const p0 = at(cx, a0 + side * s), p3 = at(cx, a0 - side * s);
    const ext = edge(c, b, p0, [cx + side * r * 3.4, cy + r * 1.9], [cx + side * r * 3.4, cy - r * 1.9], p3);
    lab(b, ext.top, false, which);
  };
  const ll = cell("ll"), ls = cell("ls"), sl = cell("sl"), ss = cell("ss");
  loop(xs, -1, col.s, ss, "ss");
  loop(xl, +1, col.l, ll, "ll");
  const bend = r * 2.2, mx = (xs + xl) / 2;
  lab(sl, edge(col.s, sl, at(xs, -0.55), [xs + r, cy - bend], [xl - r, cy - bend], at(xl, Math.PI + 0.55)).top, false, "sl");
  lab(ls, edge(col.l, ls, at(xl, Math.PI - 0.55), [xl - r, cy + bend], [xs + r, cy + bend], at(xs, 0.55)).bot, true, "ls");
  [[xs, "s"], [xl, "l"]].forEach(([cx, t]) => {
    g.fillStyle = col[t]; g.beginPath(); g.arc(cx, cy, r, 0, 2 * Math.PI); g.fill();
    g.strokeStyle = css("--panel"); g.lineWidth = 2; g.stroke();
    g.fillStyle = "#fff"; g.font = `600 ${Math.min(10, Math.floor(r * 0.62))}px ` + css("--f-ui"); g.textAlign = "center";
    g.fillText(t === "s" ? "std" : "long", cx, cy + 3.5);
  });
  g.textAlign = "left";
  cv._ptHit = hits; axisTip(cv);
}

/* fig3's row-3 violin, one per transmitter type -- standard in green, long in
   the scenario's colour, as the arrows -- from the pooled per-integer counts /api/panel sends:
   Gaussian KDE at Scott's bandwidth, cut at the data's range, width-normalised
   (seaborn's density_norm="width"), quartiles dashed. `ymax` is per subpanel,
   shared by its two violins. */
/* At most P2_SUB events per violin, drawn without replacement from the pooled
   per-integer counts with a fixed seed: the bandwidth shrinks as n^-1/5, so a
   violin on 1.2 million events and one on 877 would be smoothed ~4x apart.
   Quartiles are the sample's, numpy's linear rule. */
const P2_SUB = 800;
function p2Sample(hist, k, seed) {
  const n = hist.reduce((a, c) => a + c, 0);
  if (n <= k) return null;
  let t = seed >>> 0;
  const rnd = () => { t = (t + 0x6D2B79F5) >>> 0; let x = Math.imul(t ^ (t >>> 15), 1 | t);
    x ^= x + Math.imul(x ^ (x >>> 7), 61 | x); return ((x ^ (x >>> 14)) >>> 0) / 4294967296; };
  const pick = new Set();                         /* Floyd: k distinct ranks in [0, n) */
  for (let j = n - k; j < n; j++) { const r = Math.floor(rnd() * (j + 1)); pick.add(pick.has(r) ? j : r); }
  const cum = []; hist.reduce((a, c, v) => (cum[v] = a + c), 0);
  const out = new Array(hist.length).fill(0), vals = [];
  pick.forEach(r => { let lo = 0, hi = cum.length - 1;
    while (lo < hi) { const m = (lo + hi) >> 1; if (cum[m] > r) hi = m; else lo = m + 1; }
    out[lo]++; vals.push(lo); });
  vals.sort((a, b) => a - b);
  const qt = p => { const x = p * (k - 1), i = Math.floor(x); return vals[i] + (vals[Math.min(k - 1, i + 1)] - vals[i]) * (x - i); };
  return {hist: out, n: k, q: [qt(0.25), qt(0.5), qt(0.75)]};
}
const P2_VIOLINS = [["standard", "std", "--std", "a standard host"],
                    ["long", "long", "--long", "a long shedder"]];
function p2Violin(cv, {jump, ymax, height, label, longCol, showN = true, sub = true, seed = 0}) {
  const paint = c => c === "--long" && longCol ? longCol : css(c);
  const {g, w, h} = fit(cv, height || 220);
  const T = 14, B = 32, L = 42, R = 4;
  const Y = v => h - B - (h - B - T) * v / ymax;
  g.font = "10px " + css("--f-num"); g.fillStyle = css("--ink-3"); g.textAlign = "right";
  g.strokeStyle = css("--line"); g.lineWidth = 1;
  const step = ymax > 16 ? 5 : ymax > 8 ? 2 : 1;
  for (let v = 0; v <= ymax; v += step) {
    g.fillText(String(v), L - 4, Y(v) + 3);
    g.globalAlpha = 0.4; g.beginPath(); g.moveTo(L, Y(v)); g.lineTo(w - R, Y(v)); g.stroke();
    g.globalAlpha = 1;
  }
  g.save(); g.translate(11, T + (h - B - T) / 2); g.rotate(-Math.PI / 2);
  g.textAlign = "center"; g.font = "600 11px " + css("--f-ui"); g.fillStyle = css("--ink-2");
  [].concat(label).forEach((t, i, a) => g.fillText(t, 0, (i - (a.length - 1) / 2) * 13 + 2)); g.restore();
  const slot = (w - L - R) / P2_VIOLINS.length, half = slot / 2 - 4, hits = [];
  P2_VIOLINS.forEach(([key, short, colour, who], k) => {
    const jv = (jump || {})[key] || {}, nAll = jv.n || 0;
    const S = sub && nAll ? p2Sample(jv.hist || [], P2_SUB, 42 + 7919 * seed + k) : null;
    const hist = S ? S.hist : jv.hist || [], n = S ? S.n : nAll, qs = S ? S.q : jv.q || [];
    const cx = L + slot * (k + 0.5);
    g.textAlign = "center"; g.font = "600 11px " + css("--f-ui"); g.fillStyle = paint(colour);
    g.fillText(short, cx, h - 19);
    g.font = "10px " + css("--f-num"); g.fillStyle = css("--ink-3");
    if (showN || !n) g.fillText(n ? `n=${n.toLocaleString()}` : "none", cx, h - 7);
    if (!n) return;
    const lo = hist.findIndex(c => c > 0);
    let hi = hist.length - 1; while (hi > lo && !hist[hi]) hi--;
    let m = 0, ss = 0;
    hist.forEach((c, v) => { m += c * v; ss += c * v * v; });
    m /= n;
    const sd = Math.sqrt(Math.max(0, ss / n - m * m));
    const bw = sd > 0 ? sd * Math.pow(n, -0.2) : 0.25;
    const dens = y => hist.reduce((a, c, v) => a + (c ? c * Math.exp(-0.5 * ((y - v) / bw) ** 2) : 0), 0);
    const ys = [], ds = [], N = 80;
    for (let i = 0; i <= N; i++) { const y = lo + (hi - lo) * i / N; ys.push(y); ds.push(dens(y)); }
    const dmax = Math.max(...ds) || 1;
    g.save(); g.beginPath(); g.rect(L, T - 1, w - L - R, h - B - T + 2); g.clip();
    g.beginPath();
    if (hi === lo) g.rect(cx - half, Y(lo) - 1.5, 2 * half, 3);
    else {
      ys.forEach((y, i) => { const x = cx + half * ds[i] / dmax; i ? g.lineTo(x, Y(y)) : g.moveTo(x, Y(y)); });
      for (let i = N; i >= 0; i--) g.lineTo(cx - half * ds[i] / dmax, Y(ys[i]));
      g.closePath();
    }
    g.fillStyle = paint(colour); g.globalAlpha = 0.8; g.fill(); g.globalAlpha = 1;
    g.strokeStyle = css("--ink-2"); g.lineWidth = 0.8; g.stroke();
    qs.forEach((q, i) => {
      const hw = hi === lo ? half : half * dens(q) / dmax;
      g.setLineDash(i === 1 ? [6, 2] : [2, 2]); g.lineWidth = i === 1 ? 2.2 : 0.8;   /* the median bold */
      g.beginPath(); g.moveTo(cx - hw, Y(q)); g.lineTo(cx + hw, Y(q)); g.stroke();
    });
    g.setLineDash([]); g.restore();
    const pct = v => fmt(100 * (hist[v] || 0) / n, 1) + "%";
    hits.push({x0: cx - slot / 2, x1: cx + slot / 2, y0: T, y1: h - B, text: {
      text: `Genome positions at which a lineage ${who} transmits differs from the one it was infected with — one value per transmission event, pooled over ${jv.runs} simulations.`,
      math: `events            ${n.toLocaleString()}${S ? ` drawn at random of ${nAll.toLocaleString()}` : ""}\n`
        + `quartiles         ${qs.map(q => fmt(q, 1)).join(" · ")}\n`
        + `mean per sim.     ${ci(jv.mean, 2)}\n`
        + hist.map((c, v) => c ? `jump = ${v}`.padEnd(18) + pct(v) : null).filter(Boolean).slice(0, 8).join("\n"),
      tex: [String.raw`\Delta_e = \#\{\,x : g_{\text{transmitted}}(x) \ne g_{\text{inherited}}(x)\,\}`]}});
  });
  g.textAlign = "left";
  cv._axHit = hits;
  axisTip(cv);
}

function p2Ratio(s) {
  const L = ((s || {}).intensity || {}).long;
  const S = ((s || {}).intensity || {}).standard;
  const miss = v => !v || v.mean === null || v.mean === undefined;
  if (miss(L) || miss(S) || !S.mean) return null;
  return {n: Math.min(L.n || 0, S.n || 0), mean: L.mean / S.mean,
          lo: null, hi: null};
}

/* settings.js lists burnin here; nothing in this rail rebuilds it. */
function rail_panel2(rail) {
  if (!rail) return;
  rail.append(card("Display", [
    checkControl({setting: "f3Titles", label: "Titles", hint: "Off: only the panel letters."}),
    checkControl({setting: "f3ShowN", label: "n = under the violins",
      hint: "Transmission events behind each violin, pooled over simulations."}),
    checkControl({setting: "f3Sub", label: `Violins on ${P2_SUB} events each`,
      hint: `Each violin drawn from ${P2_SUB} of its transmission events at random (fixed seed), so all eight are smoothed alike. Off: every event.`}),
    checkControl({setting: "f3Dots", label: "Colour behind letters, panels"}),
    checkControl({setting: "f3LegLetters", label: "Letters in the legend", hint: "Off: the coloured dots alone."}),
    checkControl({setting: "f3TextCol", label: "Scenario names in colour"}),
  ]));
  rail.append(dataCard(`Whether events before the burn-in day count is set on
    <a href="#clustering">Clustering</a> (the burn-in box of the clade definition).`));
}

/* `idp`/`compact` are unused since the Figure tab went; `compact` draws the chart alone. */
async function render_panel2(main, rail, {idp = "p2", compact = false} = {}) {

  let d; try { d = await apiGetSlot("/api/panel?n=2&" + buildQuery()); }
  catch (e) { main.innerHTML = ""; main.append(errorBox(e.message)); return; }
  const sOf = sc => d.scenarios[sc] || {};
  const cellOf = (sc, c) => asPct((sOf(sc).where || {})[c]) || null;
  const H = compact ? 170 : 220;
  const noJump = SUB_SCENARIOS.some(sc => !sOf(sc).jump);
  /* each subpanel its own scale, 8% over its largest jump; fig3 shares one
     across scenarios, which left SOT and HIV in the bottom fifth of edge case's */
  const ymax1 = sc => Math.ceil(1.08 * Math.max(1, ...P2_VIOLINS.map(([t]) =>
    (((sOf(sc).jump || {})[t] || {}).hist || []).length - 1)));
  /* A-C on one fixed scale, 0-6, so they compare at a glance; edge case, whose
     long shedders reach ~25, keeps its own. Anything above 6 is clipped. */
  const ymaxOf = sc => sc === "edge_case" ? ymax1(sc) : 6;

  const body = `<div class="scgrid">${SUB_SCENARIOS.map((sc, i) => `<div>
    <h4>${figLetter(sc, FIG_LETTERS[i], settings.f3Dots)}${settings.f3Titles ? figName(sc, cap1(SCENARIO_LABEL[sc]), settings.f3TextCol) : ""}</h4>
    <div class="p2pair"><div><canvas id="${idp}_vi_${sc}"></canvas></div>
      <div><canvas id="${idp}_ar_${sc}"></canvas></div></div></div>`).join("")}</div>`;
  const row = (sw, t) => `<div class="r">${sw}<span>${t}</span></div>`;
  const sq = c => `<i class="c4sw" style="background:${c}"></i>`;
  const key = `<aside class="c4leg">
    <h5>Transmitter</h5>
    <div class="r"><span class="c4head"><b class="c4let on${settings.f3LegLetters ? " lg" : ""}" style="background:${css("--std")}"></b><b${settings.f3TextCol ? ` style="color:${css("--std")}"` : ""}>standard host</b></span></div>
    ${SUB_SCENARIOS.map((sc, i) => `<div class="r"><span class="c4head">${figLetter(sc,
      settings.f3LegLetters ? FIG_LETTERS[i] : "", true)}<b>${figName(sc, SCENARIO_LABEL[sc], settings.f3TextCol)}</b></span></div>`).join("")}
    <h5>Violins</h5>
    ${row(dashSwatch(6, 2, 2.2), "median")}
    ${row(dashSwatch(2, 2, 0.8), "quartiles")}
  </aside>`;
  const cap = `One subpanel per scenario with long shedders; control has none.
    <strong>Left, the violins</strong> are what figure 3's bottom row draws: for
    every transmission, the number of genome positions at which the lineage
    passed on differs from the one the transmitter was infected with — how much
    variation a host adds to what it transmits. Standard transmitters in green,
    long shedders in their scenario's colour; the number under each is its events. Events are pooled
    over simulations, so a simulation with more of them weighs more; the dashed
    lines are the quartiles, and each subpanel has its own y scale,
    shared by its two violins — read the tick labels before comparing across
    scenarios. Hover a
    violin for its counts. <strong>Right, the arrows:</strong> each is one
    transmitter→recipient pair. Its number is that pair's share of each
    simulation's events, averaged over simulations, so every simulation weighs the
    same however many events it had; the four in a subpanel still sum to 100%. The loops are
    long→long and standard→standard; the arc over the top is standard→long, the
    one underneath long→standard. An arrow's colour is its transmitter: green
    from a standard host, the scenario's colour from a long shedder.
    <strong>Width is the square root of the share</strong>, on one scale for all
    four subpanels: on a linear scale standard→standard, near 100%, left the
    other three as hairlines. These shares count events, so they rise with
    prevalence — the per-host rates are in the Summary table.`;
  const stale = noJump ? `<div class="note"><strong>The transmission table on disk
    predates the jump export</strong>, so the violins are empty. Re-export it with
    <code>figures_gui/slurm_tx.sh</code>.</div>` : "";
  const paint = () => SUB_SCENARIOS.forEach(sc => {
    p2Violin(document.getElementById(`${idp}_vi_${sc}`), {jump: sOf(sc).jump,
      ymax: ymaxOf(sc), height: H, label: ["new substitutions", "on transmission"], longCol: scenarioRaw(sc), showN: !!settings.f3ShowN,
      sub: !!settings.f3Sub, seed: SCENARIOS.indexOf(sc)});
    p2Arrows(document.getElementById(`${idp}_ar_${sc}`),
      {cell: c => cellOf(sc, c), height: H, longCol: scenarioRaw(sc)});
  });

  const plot = `<div class="plot"><div class="c4wrap">${body}${key}</div></div>`;
  if (compact) { main.innerHTML = plot; paint(); return; }

  const rows = SCENARIOS.map(sc => { const s = sOf(sc), j = s.jump; return `<tr>
    <td>${SCENARIO_LABEL[sc]}</td>
    <td class="n">${s.runs ?? 0}</td>
    <td class="n">${ci(asPct(s.who_transmits), 2, "%")}</td>
    <td class="n">${ci((s.intensity || {}).long, 2)}</td>
    <td class="n">${ci((s.intensity || {}).standard, 2)}</td>
    <td class="n">${ci(p2Ratio(s), 2)}</td>
    ${P2_VIOLINS.map(([t]) => { const v = (j || {})[t];
      return `<td class="n">${v && v.n ? ci(v.mean, 2) : "—"}</td>`; }).join("")}</tr>`; }).join("");
  main.innerHTML = `${stale}
  ${plot}
  <div class="plot"><div class="ph"><h3>Panel notes</h3></div>
    <div class="cap">${cap} ${cap1(simWords(d.note))}.</div></div>
  <div class="plot"><div class="ph"><h3>Summary</h3>
    <span class="sub">per scenario, all simulations</span></div>
    <table class="t"><thead><tr><th>scenario</th><th>simulations</th>
      <th>caused by long shedders</th><th>events per long shedder</th>
      <th>events per standard host</th><th>long ÷ standard</th>
      <th>mean jump, standard</th><th>mean jump, long</th></tr></thead><tbody>${rows}</tbody></table>
    <div class="cap">Events per host divide by hosts, not events, so unlike the
      arrow shares they do not rise with prevalence. long ÷ standard is one mean
      over another and carries no interval: dividing two CIs is not a CI of a
      ratio. Mean jump is each simulation's mean over that transmitter type's
      transmissions, bootstrapped over simulations. control has no long
      shedders, so its per-long-shedder values are absent, not zero.</div></div>`;
  paint();

}


registerTab({
  id: "panel2", group: "figures", num: "3", title: "Transmission", viewId: "viewP2",
  lede: "A <strong>transmission event</strong> is one infection: a transmitter, a recipient, a day and the lineage passed on. Per scenario, how much variation standard hosts and long shedders add to what they transmit, beside where the transmissions go.",
  rail: rail_panel2,
  render: render_panel2,
});
