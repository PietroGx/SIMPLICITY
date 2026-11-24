/* Figure 4: share of hosts against share of infections and frequency AUC, its one
   chart since v2.5.20 (the clade-points, prevalence and class-scatter forms
   are gone). The table and the note sit under it. */

function rail_panel1(rail) {
  rail.append(dataCard(DATA_CLADE));
  rail.append(card("Colours", [segmentControl({setting: "p1Pal", label: "Colour scheme", options: P1_SCHEMES,
    hint: "<em>Colour-blind colours</em>: lightness steps down from red to yellow and the "
        + "green turns blue-green, so no two classes merge under red–green colour blindness. "
        + "<em>Texture</em>: the same colours, plus a pattern per class, for greyscale print."})]));
}

const P1_SCHEMES = [["normal", "Normal"], ["cvd", "Colour-blind colours"], ["texture", "Colour-blind texture"]];
/* checked with the dataviz validator, all pairs: worst colour-blind ΔE 8.4, worst normal-vision ΔE 19.0 */
const P1_CVD = {only_long: "#a50026", mostly_long: "#f46d43", mostly_std: "#fdbf4f", standard: "#1b7f6a"};
const P1_TEX = {only_long: "cross", mostly_long: "horiz", mostly_std: "dots", founding: "hatch"};
function p1Scheme() {
  const pal = settings.p1Pal === "normal"
    ? {only_long: css("--long"), mostly_long: css("--mostly-long"), mostly_std: css("--mostly-std"), standard: css("--std")}
    : P1_CVD;
  pal.founding = pal.standard;
  const tex = k => settings.p1Pal === "texture" || k === "founding" ? P1_TEX[k] || null : null;
  const light = k => k === "mostly_std";
  const ink = k => light(k) ? "rgba(0,0,0,.45)" : "rgba(255,255,255,.7)";
  /* a textured swatch is a small drawn image: the export's html2canvas does not
     paint CSS repeating gradients */
  const drawn = k => { try {
      const cv = document.createElement("canvas"); cv.width = cv.height = 24;
      const g = cv.getContext("2d"); g.scale(2, 2); g.fillStyle = pal[k]; g.fillRect(0, 0, 12, 12);
      p1Texture(g, tex(k), 0, 0, 12, 12, ink(k));
      return `url(${cv.toDataURL("image/png")}) 0 0/100% 100%`;
    } catch (e) { return null; } };
  const swatch = k => { const c = pal[k], i = ink(k);
    if (tex(k)) { const u = drawn(k); if (u && u.startsWith("url(data:image/png")) return u; }
    switch (tex(k)) {
      case "hatch": return `repeating-linear-gradient(135deg,${c} 0 3px,rgba(255,255,255,.75) 3px 4.5px)`;
      case "cross": return `repeating-linear-gradient(45deg,${i} 0 1px,transparent 1px 5px),repeating-linear-gradient(135deg,${i} 0 1px,transparent 1px 5px),${c}`;
      case "horiz": return `repeating-linear-gradient(0deg,${i} 0 1px,transparent 1px 4px),${c}`;
      case "dots": return `radial-gradient(${i} 1px,transparent 1.3px) 0 0/4px 4px,${c}`;
      default: return c;
    } };
  return {pal, tex, light, ink, swatch};
}
/* a class's texture over a filled segment */
function p1Texture(g, kind, x, y, w, h, ink) {
  if (!kind) return;
  g.save(); g.beginPath(); g.rect(x, y, w, h); g.clip();
  const line = (x0, y0, x1, y1) => { g.beginPath(); g.moveTo(x0, y0); g.lineTo(x1, y1); g.stroke(); };
  if (kind === "hatch") { g.strokeStyle = "rgba(255,255,255,.75)"; g.lineWidth = 1.5;
    for (let u = x - h; u < x + w; u += 4.5) line(u, y + h, u + h, y); }
  else if (kind === "cross") { g.strokeStyle = ink; g.lineWidth = 1;
    for (let u = x - h; u < x + w; u += 6) { line(u, y + h, u + h, y); line(u, y, u + h, y + h); } }
  else if (kind === "horiz") { g.strokeStyle = ink; g.lineWidth = 1;
    for (let v = y + 2.5; v < y + h; v += 4) line(x, v, x + w, v); }
  else if (kind === "dots") { g.fillStyle = ink;
    for (let v = y + 2; v < y + h; v += 4) for (let u = x + 2; u < x + w; u += 4) { g.beginPath(); g.arc(u, v, 0.9, 0, 6.2832); g.fill(); } }
  g.restore();
}

async function panel1Form(d, {idp = "p1", compact = false} = {}) {
  const of = sc => d.scenarios[sc] || {};
  const SUBS = SUB_SCENARIOS;
  const sub = `k=${settings.k} · m_T=${mtPct(settings.mt)}% · ${settings.burnin?"post burn-in":"whole simulation"}`;
  const scNote = sc => `${fmt(asPctNum(of(sc).prevalence), 1)}% · ${of(sc).duration ?? "—"}d`;
  const box = (tag, title, key, cap, note = scNote) => `<div class="plot">
    <div class="ph"><h3>${title}</h3><span class="sub">${sub}</span></div>
    <div class="scgrid">${SUBS.map(sc => `<div>
      ${scenarioHead(sc, note(sc))}
      <canvas id="${idp}${tag}_${sc}"></canvas></div>`).join("")}</div>
    ${key || ""}${compact ? "" : `<div class="cap">${cap}</div>`}</div>`;
  const each = (tag, f) => () => SUBS.forEach(sc => f(document.getElementById(`${idp}${tag}_${sc}`), sc));

  const ph = await apiGetSlot("/api/panel?n=1h&" + buildQuery()).catch(e => ({error: e.message}));
  const phOf = sc => ((ph || {}).scenarios || {})[sc] || {};
  /* red first; the founding clade joins the green segment, hatched */
  const S = p1Scheme();
  const HB = [["only_long", "only long (all definers long-born)"],
              ["mostly_long", "mostly long (> half, not all)"],
              ["mostly_std", "mostly standard (some, ≤ half)"],
              ["standard", "standard only (no long-born definer)"],
              ["founding", "founding clade"]];
  const infWhole = SUBS.some(sc => phOf(sc).infections_window !== "post");
  const infWin = infWhole ? "whole simulation" : "on or after the burn-in day";
  const pctTxt = v => !isNum(v) ? "—" : fmt(v, 1) + "%";

  return {
    html: box("hb", "Long shedders: share of hosts against share of infections and frequency AUC",
      keyRow([[S.swatch("only_long"), "long shedders (hosts bar)"], [S.swatch("standard"), "standard hosts (hosts bar)"]]
        .concat(HB.map(([k, t]) => [S.swatch(k), k === "founding" ? "of which founding clade (hatched)" : t]))),
      `Three bars per scenario, each out of 100%. <strong>Hosts:</strong> the share
       of infected individuals who are long shedders, in each simulation as
       realized. <strong>Infections:</strong> every infection in the simulation,
       split by the clade of the lineage it passed on — the plainest measure of
       spread. ${infWhole
         ? "These clade builds count infections over the whole simulation only, so this bar ignores the burn-in setting until the next rebuild."
         : "Counted on or after the burn-in day, as frequency AUC is."} <strong>Frequency AUC:</strong> each
       simulation's frequency AUC (the area under the clade's frequency curve: its share
       of infected hosts, summed over the days in the window) split by where
       the clade's defining mutations arose — <strong>only long</strong>: every definer
       arose in a long shedder; <strong>mostly long</strong>: more than half, not all;
       <strong>mostly standard</strong>: at least one, at most half; <strong>standard
       only</strong>: none. The green segment also holds the founding clade (it has no
       definers), hatched, and its label is the two together.
       Infections and frequency AUC share the same segments.${settings.p1Pal === "texture"
         ? " Textures: crossed lines, only long; horizontal lines, mostly long; dots, mostly standard; plain, standard only."
         : ""} The headline compares
       the hosts share with the infections carried by clades with <em>at least
       one</em> long-born definer (red, orange and yellow together).
       Means over simulations; hover a segment for its 95% interval.
       A long shedder stays infected far longer than a standard host, so part
       of the over-representation is longer infection itself, not fitter
       mutations.${(ph || {}).error ? ` <strong>Not available:</strong> ${ph.error}` : ""}`,
      sc => scNote(sc)),
    draw: each("hb", (cv, sc) => {
      const v = phOf(sc), b = v.burden || {}, inf = v.infections || {}, hm = asPctNum((v.hosts || {}).mean);
      const {g, w, h} = fit(cv, compact ? 144 : 158);
      const L = 62, R = 6, bh = 20, y1 = 30, y2 = 66, y3 = 102;
      const W = w - L - R, any = asPctNum((inf.any_long || {}).mean);
      g.font = "600 11px " + css("--f-ui"); g.fillStyle = css("--ink"); g.textAlign = "left";
      g.fillText(isNum(hm) && isNum(any)
        ? `${pctTxt(hm)} of hosts → ${pctTxt(any)} of infections${hm > 0 ? ` · ${fmt(any / hm, 0)}×` : ""}`
        : "not available", 0, 13);
      const hits = [];
      const bar = (y, label, segs) => {
        g.font = "600 10px " + css("--f-ui"); g.fillStyle = css("--ink-2"); g.textAlign = "right";
        g.fillText(label, L - 6, y + bh / 2 + 4);
        let x = L;
        const runs = [];      // a hatched segment joins the run before it: one outline, one label
        /* a non-zero segment is never thinner than 5 px; the widest one gives up the difference */
        const ws = segs.map(([val]) => isNum(val) && val > 0 ? Math.max(5, W * val / 100) : 0);
        const big = ws.indexOf(Math.max(...ws)), tot = ws.reduce((t, u) => t + u, 0);
        if (tot > W && big >= 0) ws[big] -= tot - W;
        segs.forEach(([val, k, name, ci], i) => {
          if (!isNum(val) || val <= 0) return;
          const sw = ws[i], col = S.pal[k];
          g.fillStyle = col; g.fillRect(x, y, sw, bh);
          p1Texture(g, S.tex(k), x, y, sw, bh, S.ink(k));
          if (k === "founding" && runs.length && runs[runs.length - 1].end === x && runs[runs.length - 1].col === col) {
            const r = runs[runs.length - 1]; r.end = x + sw; r.val += val; r.plain = false;
          } else runs.push({x, end: x + sw, solid: x + sw, val, col, plain: !S.tex(k),
                            ink: S.light(k) ? css("--ink") : "#fff"});
          hits.push({x0: x, x1: x + sw, y0: y, y1: y + bh, text: {text: `${SCENARIO_LABEL[sc]} · ${name}`,
            math: `mean over ${v.runs} simulations  ${pctTxt(val)}` + (ci ? `\n95% interval  ${ci}` : "")}});
          x += sw;
        });
        const above = [];         // labels too wide for their segment, set just above the bar
        runs.forEach(r => { const sw = r.end - r.x, txt = pctTxt(r.val);
          g.strokeStyle = css("--panel"); g.lineWidth = 1;
          if (sw >= 8) g.strokeRect(r.x + 0.5, y + 0.5, sw - 1, bh - 1);
          else if (r.end < L + W - 1) { g.beginPath(); g.moveTo(r.end - 0.5, y); g.lineTo(r.end - 0.5, y + bh); g.stroke(); }
          g.font = "600 10px " + css("--f-num"); g.textAlign = "center";
          const tw = g.measureText(txt).width + 6, lw = r.solid - r.x > tw ? r.solid - r.x : sw;
          if (sw > tw) {
            if (!r.plain) { g.fillStyle = r.col; g.fillRect(r.x + lw / 2 - tw / 2, y + 4, tw, bh - 8); }   // a clear patch under the label
            g.fillStyle = r.ink; g.fillText(txt, r.x + lw / 2, y + bh / 2 + 4); return; }
          g.font = "600 9px " + css("--f-num");
          above.push({txt, cx: r.x + sw / 2, hw: g.measureText(txt).width / 2 + 2}); });
        /* no overlaps, inside the bar's span: push right, then back left from the end */
        above.forEach((a, i) => { a.cx = Math.max(a.cx, L + a.hw, i ? above[i - 1].cx + above[i - 1].hw + a.hw : -Infinity); });
        for (let i = above.length - 1; i >= 0; i--)
          above[i].cx = Math.min(above[i].cx, i < above.length - 1 ? above[i + 1].cx - above[i + 1].hw - above[i].hw : L + W - above[i].hw);
        g.font = "600 9px " + css("--f-num"); g.fillStyle = css("--ink-2"); g.textAlign = "center";
        above.forEach(a => g.fillText(a.txt, a.cx, y - 3));
        g.strokeStyle = css("--line"); g.strokeRect(L + 0.5, y + 0.5, W - 1, bh - 1);
      };
      const ciOf = o => o && isNum(o.lo) ? `${pctTxt(asPctNum(o.lo))} – ${pctTxt(asPctNum(o.hi))}` : null;
      bar(y1, "hosts", [[hm, "only_long", "long shedders among infected", ciOf(v.hosts)],
                        [isNum(hm) ? 100 - hm : null, "standard", "standard hosts", null]]);
      bar(y2, "infections", HB.map(([k, t]) => [asPctNum((inf[k] || {}).mean), k, `${t} · share of infections, ${infWin}`, ciOf(inf[k])]));
      bar(y3, "freq. AUC", HB.map(([k, t]) => [asPctNum((b[k] || {}).mean), k, `${t} · share of frequency AUC`, ciOf(b[k])]));
      g.font = "9px " + css("--f-num"); g.fillStyle = css("--ink-3");
      g.textAlign = "left"; g.fillText("0%", L, y3 + bh + 13);
      g.textAlign = "right"; g.fillText("100%", L + W, y3 + bh + 13); g.textAlign = "left";
      cv._axHit = hits; axisTip(cv);
    }),
  };
}

async function render_panel1(main){

  let d; try { d = await apiGetSlot("/api/panel?n=1&" + buildQuery()); }
  catch(e){ main.innerHTML=""; main.append(errorBox(e.message + " — build this k on the Export page.")); return; }
  const cls = d.classes;
  const lab = c => c === "wt" ? "wild type" : c === "other" ? "ndef≠k" : c;
  const of = sc => d.scenarios[sc] || {};
  const form = await panel1Form(d);
  const tbl = SCENARIOS.map(sc => {
    const s = of(sc);
    const cc = s.classes || {};
    return `<tr><td>${SCENARIO_LABEL[sc]} <small style="color:var(--ink-3)">${fmt(asPctNum(s.prevalence),1)}% · ${s.duration ?? "—"}d</small></td>`
      + cls.map(c => `<td class="n">${fmt(asPctNum((cc[c]||{}).frac_burden?.mean),1)}</td>`).join("")
      + `</tr>`;
  }).join("");
  main.innerHTML = `
  <div class="note">${simWords(d.note)}.
    Percentages are formed inside each simulation and then averaged over
    simulations; whiskers are a bootstrap over simulations, never over clades.</div>
  ${form.html}
  <div class="plot"><div class="ph"><h3>Percentage of frequency AUC, by class</h3></div>
    <table class="t"><thead><tr><th>scenario</th>
      ${cls.map(c => `<th>${lab(c)}</th>`).join("")}</tr></thead>
      <tbody>${tbl}</tbody></table></div>`;
  form.draw();

}


registerTab({
  id: "panel1", group: "figures", num: "4", title: "Frequency AUC", viewId: "viewP1",
  lede: "How much of the viral population's history traces to mutations that arose in long shedders, against how rare long shedders are.",
  rail: rail_panel1,
  render: render_panel1,
});
