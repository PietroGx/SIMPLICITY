/* Figure 1 · Model: the six panels of figures_gui/prep/fig1_save.py, drawn
   from what that script drew. figures_gui/fig1.py runs build_figure_1 and
   records each painter's frames, limits, outlines and fits, and stops if its
   clock fits differ from the lines the figure drew; nothing here computes a
   clock, a fit or a density. Canvas ids are f1_*. */

API_LABELS["/api/fig1"] = "reading Figure 1";

const F1_H = 240;
let f1Last = null;
/* the Display switches live in settings (settings.js), so a reload keeps them */
const f1On = k => !!settings[k];
/* bottom padding each canvas needs, and the floor its row imposes, so the plot
   areas of one row share their x axis line (f1Align) */
const f1PadNeed = {}, f1PadMin = {};

const f1Data = () => css("--fig1-data");
const f1Col = (d, cohort) => scenarioRaw((d.cohort_scenario || {})[cohort] || "control");
const f1Role = (d, role) => role === "unmodelled" ? css("--fig1-unmod") : f1Col(d, role);
/* the figure's own cohort label for HIV: it is one scenario (HIV_low) standing for the clinical group */
const f1Name = (d, c) => { const sc = (d.cohort_scenario || {})[c];
  return c === "HIV" ? "HIV" : sc ? SCENARIO_LABEL[sc] : c; };
const f1Code = s => `<code>${s}</code>`;

function f1Ticks(a, b, n = 5) {
  const raw = (b - a) / n, p = Math.pow(10, Math.floor(Math.log10(raw || 1)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * p).find(s => s >= raw) || raw || 1;
  const out = [];
  for (let v = Math.ceil(a / step) * step; v <= b + step * 1e-9; v += step) out.push(+v.toPrecision(10));
  return out;
}

/* axes() for y; x ticks, or category labels with a second line, here. Leaves
   the context clipped to the plot area: the caller restores. */
function f1Frame(cv, {x0, x1, y0, y1, xlab, ylab, fmtX, fmtY, cats, sub}) {
  const {g, w, h} = fit(cv, F1_H);
  /* category names tilt, as the original's do, when they do not fit flat */
  const catFont = "600 10px " + css("--f-ui");
  g.font = catFont;
  const catW = cats ? Math.max(...cats.map(c => g.measureText(String(c)).width || 0)) : 0;
  const tilt = cats && catW > 0.95 * (w - 58) / Math.max(1, x1 - x0);
  const need = !cats ? 34 : tilt ? Math.ceil(16 + (sub ? 12 : 0) + catW * 0.66) : (sub ? 40 : 30);
  f1PadNeed[cv.id] = need;
  const pad = {l: 46, r: 12, t: 10, b: Math.max(need, f1PadMin[cv.id] || 0)};
  /* a row's deeper padding would leave the label far below a plain axis: draw it under the ticks */
  const lift = pad.b - need;
  axes(g, w, h, pad, lift ? null : xlab, ylab, y0, y1, fmtY);
  if (lift && xlab) { g.font = "600 11px " + css("--f-ui"); g.fillStyle = css("--ink-2"); g.textAlign = "center";
    g.fillText(xlab, (pad.l + w - pad.r) / 2, h - 3 - lift); g.textAlign = "left"; }
  const X = v => pad.l + (w - pad.l - pad.r) * (v - x0) / ((x1 - x0) || 1);
  const Y = v => h - pad.b - (h - pad.b - pad.t) * (v - y0) / ((y1 - y0) || 1);
  g.strokeStyle = css("--line"); g.textAlign = "center";
  const ticks = cats ? cats.map((c, i) => [i, c]) : f1Ticks(x0, x1).map(v => [v, fmtX ? fmtX(v) : String(v)]);
  ticks.forEach(([v, t], i) => {
    const x = X(v);
    g.beginPath(); g.moveTo(x, h - pad.b); g.lineTo(x, h - pad.b + 3); g.stroke();
    if (!cats && v === 0 && y0 === 0) return;   /* the y axis already prints the origin's 0 */
    g.fillStyle = css(cats ? "--ink-2" : "--ink-3");
    g.font = cats ? catFont : "10px " + css("--f-num");
    if (tilt) {
      if (sub) { g.fillStyle = css("--ink-3"); g.font = "9px " + css("--f-num");
                 g.fillText(sub[i], x, h - pad.b + 12); }
      g.save(); g.translate(x + 3, h - pad.b + (sub ? 22 : 11)); g.rotate(-0.72);
      g.textAlign = "right"; g.fillStyle = css("--ink-2"); g.font = catFont;
      g.fillText(t, 0, 0); g.restore();
      return;
    }
    /* a label wider than the space right of its tick ends at the canvas edge instead of being cut */
    const tw = g.measureText(t).width;
    if (x + tw / 2 > w - 1) { g.textAlign = "right"; g.fillText(t, w - 1, h - pad.b + 13); g.textAlign = "center"; }
    else g.fillText(t, x, h - pad.b + 13);
    if (sub) { g.fillStyle = css("--ink-3"); g.font = "9px " + css("--f-num");
               g.fillText(sub[i], x, h - pad.b + 25); }
  });
  g.textAlign = "left";
  g.save(); g.beginPath(); g.rect(pad.l, pad.t - 1, w - pad.l - pad.r, h - pad.t - pad.b + 1); g.clip();
  return {g, w, h, pad, X, Y};
}

const f1Path = (g, xs, ys, X, Y) => {
  g.beginPath();
  xs.forEach((x, i) => i ? g.lineTo(X(x), Y(ys[i])) : g.moveTo(X(x), Y(ys[i])));
};

/* ------------------------------------------------------------ A */
function f1PaintA(cv, d) {
  const A = d.A;
  const {g, X, Y} = f1Frame(cv, {x0: A.xlim[0], x1: A.xlim[1], y0: A.ylim[0], y1: A.ylim[1],
    xlab: "days since infection", ylab: "probability of being infectious", fmtY: v => v.toFixed(2)});
  A.curves.forEach(c => {
    g.strokeStyle = f1Col(d, c.cohort); g.lineWidth = 1.7; g.globalAlpha = 0.85;
    f1Path(g, c.t, c.p, X, Y); g.stroke();
  });
  g.globalAlpha = 1; g.restore();
}
/* the first day after its peak a curve is below one half, read off the drawn curve */
const f1Half50 = c => {
  const top = c.p.indexOf(Math.max(...c.p));
  const i = c.p.findIndex((p, j) => j > top && p < 0.5);
  return i < 0 ? null : c.t[i];
};

/* ------------------------------------------------------------ B */
/* a violin's half-width at y, from the left edge of the outline seaborn drew */
function f1HalfWidth(v, y) {
  for (let j = 1; j < v.y.length; j++) {
    const a = v.y[j - 1], b = v.y[j];
    if (v.dx[j - 1] <= 0 && v.dx[j] <= 0 && (a - y) * (b - y) <= 0 && a !== b)
      return -(v.dx[j - 1] + (v.dx[j] - v.dx[j - 1]) * (y - a) / (b - a));
  }
  return 0.1;
}
function f1PaintB(cv, d) {
  const B = d.B, V = B.violins;
  const {g, X, Y} = f1Frame(cv, {x0: B.xlim[0], x1: B.xlim[1], y0: B.ylim[0], y1: B.ylim[1],
    ylab: "realized duration (days)", fmtY: v => v.toFixed(0),
    cats: V.map(v => f1Name(d, v.cohort)), sub: f1On("f1ShowN") ? V.map(v => `n=${v.n.toLocaleString()}`) : null});
  const hits = [];
  V.forEach((v, i) => {
    f1Path(g, v.dx.map(dx => i + dx), v.y, X, Y); g.closePath();
    g.fillStyle = f1Col(d, v.cohort); g.globalAlpha = 0.85; g.fill(); g.globalAlpha = 1;
    g.strokeStyle = css("--ink"); g.lineWidth = 0.8; g.stroke();
    v.q.forEach((q, k) => {
      const hw = f1HalfWidth(v, q);
      g.setLineDash(k === 1 ? [5, 2] : [2, 2]); g.lineWidth = k === 1 ? 1.1 : 0.8;
      g.beginPath(); g.moveTo(X(i - hw), Y(q)); g.lineTo(X(i + hw), Y(q)); g.stroke();
    });
    g.setLineDash([]);
    hits.push({x: X(i), y: Y(v.q[1]), text: {
      text: `${f1Name(d, v.cohort)}: realized duration of ${v.cohort === "Control" ? "standard individuals" : "long shedders"}, t_not_infected − t_infection.`,
      math: `individuals  ${v.n.toLocaleString()}${v.n_all > v.n ? ` drawn at random of ${v.n_all.toLocaleString()}` : ""}\nquartiles    ${v.q.map(q => fmt(q, 1)).join(" · ")} d\nmean         ${fmt(v.mean, 1)} d`}});
  });
  g.restore();
  cv._ptHit = hits; axisTip(cv);
}

/* ------------------------------------------------------------ C */
function f1PaintC(cv, d) {
  const C = d.C, K = C.categories;
  const {g, X, Y} = f1Frame(cv, {x0: C.xlim[0], x1: C.xlim[1], y0: C.ylim[0], y1: C.ylim[1],
    ylab: "literature duration (days)", fmtY: v => v.toFixed(0),
    cats: K.map(k => k.cat), sub: f1On("f1ShowN") ? K.map(k => `n=${k.n}`) : null});
  const hits = [];
  let j = 0;
  K.forEach((k, i) => {
    const b = k.box, hw = 0.2, col = f1Role(d, k.role);
    g.fillStyle = col; g.globalAlpha = 0.35;
    g.fillRect(X(i - hw), Y(b.q3), X(i + hw) - X(i - hw), Y(b.q1) - Y(b.q3));
    g.globalAlpha = 1; g.strokeStyle = css("--ink-3"); g.lineWidth = 1;
    g.strokeRect(X(i - hw), Y(b.q3), X(i + hw) - X(i - hw), Y(b.q1) - Y(b.q3));
    g.beginPath();
    g.moveTo(X(i), Y(b.q3)); g.lineTo(X(i), Y(b.whishi)); g.moveTo(X(i - hw / 2), Y(b.whishi)); g.lineTo(X(i + hw / 2), Y(b.whishi));
    g.moveTo(X(i), Y(b.q1)); g.lineTo(X(i), Y(b.whislo)); g.moveTo(X(i - hw / 2), Y(b.whislo)); g.lineTo(X(i + hw / 2), Y(b.whislo));
    g.stroke();
    g.strokeStyle = css("--ink"); g.lineWidth = 1.4;
    g.beginPath(); g.moveTo(X(i - hw), Y(b.med)); g.lineTo(X(i + hw), Y(b.med)); g.stroke();
    k.pts.forEach(([v, role]) => {
      g.fillStyle = f1Role(d, role); g.globalAlpha = 0.8;
      g.beginPath(); g.arc(X(i + jitter(j++, 0.2)), Y(v), 2.6, 0, 6.2832); g.fill();
    });
    g.globalAlpha = 1;
    hits.push({x: X(i), y: Y(b.med), text: {text: `${k.cat}: published shedding durations, one per patient.`,
      math: `patients   ${k.n}\nmedian     ${fmt(b.med, 0)} d\nquartiles  ${fmt(b.q1, 0)} – ${fmt(b.q3, 0)} d\n`
        + `whiskers   ${fmt(b.whislo, 0)} – ${fmt(b.whishi, 0)} d\nmean       ${fmt(k.mean, 0)} d`}});
  });
  g.restore();
  cv._ptHit = hits; axisTip(cv);
}

/* ------------------------------------------------------------ D, F */
function f1PaintClock(cv, d, P, xlab) {
  const {g, X, Y} = f1Frame(cv, {x0: P.xlim[0], x1: P.xlim[1], y0: P.ylim[0], y1: P.ylim[1],
    xlab, ylab: "divergence (10⁻³ subs/site)", fmtY: v => (v * 1e3).toFixed(2)});
  const xm = P.x_max;
  if (f1On("f1Cloud")) P.model.forEach(m => {
    g.fillStyle = f1Col(d, m.cohort); g.globalAlpha = 0.25;
    m.thin.x.forEach((x, i) => { g.beginPath(); g.arc(X(x), Y(m.thin.y[i]), 1.5, 0, 6.2832); g.fill(); });
  });
  if (P.data) {
    g.fillStyle = f1Data(); g.globalAlpha = 0.2;
    P.data.x.forEach((x, i) => { g.beginPath(); g.arc(X(x), Y(P.data.y[i]), 1.8, 0, 6.2832); g.fill(); });
  }
  P.model.forEach(m => {
    const c = f1Col(d, m.cohort);
    if (m.sd > 0) {
      g.fillStyle = c; g.globalAlpha = 0.16; g.beginPath(); g.moveTo(X(0), Y(0));
      g.lineTo(X(xm), Y((m.slope + m.sd) * xm)); g.lineTo(X(xm), Y((m.slope - m.sd) * xm));
      g.closePath(); g.fill();
    }
  });
  g.globalAlpha = 1;
  P.model.forEach(m => {
    g.strokeStyle = f1Col(d, m.cohort); g.lineWidth = 1.6;
    g.beginPath(); g.moveTo(X(0), Y(0)); g.lineTo(X(xm), Y(m.slope * xm)); g.stroke();
  });
  if (P.data) {
    g.strokeStyle = f1Data(); g.lineWidth = 1.6; g.setLineDash([6, 3]);
    g.beginPath(); g.moveTo(X(0), Y(0)); g.lineTo(X(xm), Y(P.data.slope * xm)); g.stroke();
    g.setLineDash([]);
  }
  g.restore();
}

function f1ClockTable(d, P) {
  const row = (name, e, sd, seeds) => `<tr><td>${name}</td><td class="n">${fmt(e.rate, 5)}</td>
    <td class="n">${sd}</td><td class="n">${fmt(e.r2, 3)}</td><td class="n">${e.n.toLocaleString()}</td>
    <td class="n">${seeds}</td></tr>`;
  return `<table class="t" style="margin-top:8px"><thead><tr><th>series</th><th>rate (s/s/y)</th>
    <th>± SD over seeds</th><th>R²</th><th>points</th><th>seeds</th></tr></thead><tbody>
    ${P.data ? row(P.data_label, P.data, "—", "—") : ""}
    ${P.model.map(m => row(`model, ${f1Name(d, m.cohort)}`, m, fmt(m.sd_rate, 5),
      m.fitted_seeds === m.seeds ? m.seeds : `${m.fitted_seeds} of ${m.seeds}`)).join("")}</tbody></table>`;
}

/* ------------------------------------------------------------ E */
function f1PaintE(cv, d) {
  const E = d.E, S = E.series || {};
  const {g, X, Y} = f1Frame(cv, {x0: E.xlim[0], x1: E.xlim[1], y0: E.ylim[0], y1: E.ylim[1],
    /* fig1_plots stores the density divided by 1,000; shown here in plain units */
    xlab: "divergence (subs/site)", ylab: "density (per subs/site)",
    fmtX: v => v.toFixed(3), fmtY: v => (v * 1e3).toFixed(0)});
  const col = k => k === "data" ? f1Data() : f1Col(d, "Control");
  ["model", "data"].forEach(k => {
    const h = (S[k] || {}).hist;
    if (!f1On("f1Hist") || !h || !h.x.length) return;
    g.strokeStyle = col(k); g.globalAlpha = 0.35; g.lineWidth = 0.8;
    g.beginPath(); g.moveTo(X(h.x[0]), Y(h.y[0]));
    for (let i = 1; i < h.x.length; i++) {
      const mx = (h.x[i - 1] + h.x[i]) / 2;
      g.lineTo(X(mx), Y(h.y[i - 1])); g.lineTo(X(mx), Y(h.y[i]));
    }
    g.lineTo(X(h.x[h.x.length - 1]), Y(h.y[h.y.length - 1])); g.stroke();
  });
  g.globalAlpha = 1;
  const km = (S.model || {}).kde, kd = (S.data || {}).kde;
  if (km) {
    g.fillStyle = col("model"); g.globalAlpha = 0.1;
    f1Path(g, km.x, km.y, X, Y); g.lineTo(X(km.x[km.x.length - 1]), Y(0)); g.lineTo(X(km.x[0]), Y(0));
    g.closePath(); g.fill(); g.globalAlpha = 1;
    g.strokeStyle = col("model"); g.lineWidth = 1.6; f1Path(g, km.x, km.y, X, Y); g.stroke();
  }
  if (kd) {
    g.strokeStyle = col("data"); g.lineWidth = 1.8; g.lineCap = "round"; g.setLineDash([0.1, 3.5]);
    f1Path(g, kd.x, kd.y, X, Y); g.stroke(); g.setLineDash([]); g.lineCap = "butt";
  }
  /* each series' own values as a rug of short ticks along the floor, model
     lowest and patients above: a value has an x but no y of its own, the y axis
     being the whole sample's density (fig1.py stores up to 1,500 of each) */
  const y0 = E.ylim[0], span = E.ylim[1] - E.ylim[0];
  [["model", f1On("f1DotsModel"), 0.0], ["data", f1On("f1DotsData"), 0.05]].forEach(([k, on, lo]) => {
    const v = (S[k] || {}).pts;
    if (!on || !v) return;
    g.strokeStyle = col(k); g.globalAlpha = 0.25; g.lineWidth = 0.7;
    const ya = Y(y0 + span * lo), yb = Y(y0 + span * (lo + 0.04));
    g.beginPath(); v.forEach(x => { g.moveTo(X(x), ya); g.lineTo(X(x), yb); }); g.stroke();
  });
  g.globalAlpha = 1;
  g.restore();
}

/* ------------------------------------------------------------ tab */
function f1PaintClocks() {
  const d = f1Last; if (!d) return;
  const D = document.getElementById("f1_D"), F = document.getElementById("f1_F");
  if (D) f1PaintClock(D, d, d.D, "days since first sample");
  if (F) f1PaintClock(F, d, d.F, "days since infection");
}

function f1Rail(rail) {
  const tog = (key, label, hint) => {
    const c = document.createElement("div");
    c.className = "ctl";
    c.innerHTML = `<label class="chk"><input type="checkbox"${f1On(key) ? " checked" : ""}><span>${label}</span></label>`
      + (hint ? `<div class="hint">${hint}</div>` : "");
    c.querySelector("input").addEventListener("change", ev => {
      setSetting(key, ev.target.checked ? 1 : 0); renderTab(tabById("fig1")); });
    return c;
  };
  rail.append(card("Display", [
    tog("f1Titles", "Panel titles", "Off: only the panel letters."),
    tog("f1ShowN", "n = under B and C",
        "The count under each violin and box: simulated individuals in B, patients in C."),
    tog("f1Hist", "Histogram in D", "The faint step outline of each series' histogram behind its density."),
    tog("f1DotsData", "Patient sequences in D",
        "Each patient sequence's divergence as a tick along the floor: the values the dotted density is drawn from."),
    tog("f1DotsModel", "Model sequences in D",
        "The time-matched model sample the solid density is drawn from, 1,500 of its sequences at random, as ticks below the patients'."),
    tog("f1Cloud", "Model points in E and F",
        "Figure 1 draws only the model's fitted lines: its thousands of points buried the patients'. "
        + "Ticked, each scenario's points are added behind, thinned to 600 with <code>fig1_plots._thin</code> "
        + "as <code>plot_clock_overlay(model_scatter=True)</code> would draw them."),
  ]));
}

/* ------------------------------------------------------------ legend
   One key for the whole figure, right of the grid: colours once, then what each
   panel's marks mean. Swatches are CSS, not SVG, so the Export tab's
   html2canvas draws them as the page does. */
const F1_CATS = {HIV: "HIV/AIDS", SOT: "solid organ transplant", "B-CN": "B-cell neoplasm",
  AuIm: "rheumatological / autoimmune", PID: "primary immunodeficiency",
  aHSCT: "haematopoietic stem-cell transplant", OHem: "other haematological disease",
  ImmC: "immunocompromised, unspecified", Other: "other"};
const f1Sw = {
  line: (c, dash) => `<i class="sw" style="border-top:2px ${dash === "dot" ? "dotted" : dash ? "dashed" : "solid"} ${c};height:0;margin-top:1px"></i>`,
  band: c => `<i class="sw"><i style="position:absolute;inset:0;background:${c};opacity:.18"></i><i style="position:absolute;left:0;right:0;top:4px;border-top:2px solid ${c}"></i></i>`,
  area: c => `<i class="sw"><i style="position:absolute;inset:0;background:${c};opacity:.12"></i><i style="position:absolute;left:0;right:0;top:0;border-top:2px solid ${c}"></i></i>`,
  dot: (c, a = 0.8) => `<i class="sw"><i style="position:absolute;left:7px;top:2px;width:6px;height:6px;border-radius:50%;background:${c};opacity:${a}"></i></i>`,
  tick: c => `<i class="sw">${[4, 8, 9, 13, 16].map(x => `<i style="position:absolute;left:${x}px;top:1px;height:8px;border-left:1px solid ${c};opacity:.6"></i>`).join("")}</i>`,
  fill: c => `<i class="sw"><i style="position:absolute;inset:0;background:${c}"></i></i>`,
};
const f1Row = (sw, text) => `<div class="r">${sw}<span>${text}</span></div>`;
/* a dashed segment as spans, so its dash lengths match the canvas's and html2canvas draws it */
const f1Dash = (on, off, wpx, c = "var(--ink)") => {
  let out = "", x = 0;
  while (x < 20) { const l = Math.min(on, 20 - x);
    out += `<i style="position:absolute;left:${x}px;top:4px;width:${l}px;border-top:${wpx}px solid ${c}"></i>`; x += on + off; }
  return `<i class="sw">${out}</i>`;
};
function f1Legend(d) {
  const cs = d.cohort_scenario || {}, dur = d.cohort_duration || {}, C = d.C, pd = f1Data(), mc = f1Col(d, "Control");
  const who = c => c === "Control" ? "standard individuals"
    : `long shedders${isNum(dur[c]) ? ` – ${fmt(dur[c], 0)} days` : ""}`;
  const txt = t => `<div class="r"><span>${t}</span></div>`;
  const clock = (P, dot, data, model) => [
    f1Row(f1Sw.dot(pd, 0.35), dot),
    f1Row(f1Sw.line(pd, true), data),
  ].concat(P.model.map(m => f1Row(f1Sw.band(f1Col(d, m.cohort)), model(m))),
    f1On("f1Cloud") ? [f1Row(f1Sw.dot("var(--ink-2)", 0.3), "model sequences, thinned to 600")] : []).join("");
  const cats = C.categories.map(k => txt(`<b>${k.cat}</b> ${F1_CATS[k.cat] || ""}`)).join("");
  return `<aside class="f1leg">
    <h5>A · intra-host model</h5>
    ${["Control", "SOT", "HIV", "Edge Case"].filter(c => cs[c])
      .map(c => f1Row(f1Sw.fill(f1Col(d, c)), `<b>${f1Name(d, c)}</b> (${who(c)})`)).join("")}
    <h5>B · realized duration</h5>
    ${f1Row(f1Dash(5, 2, 1.5), "median")}
    ${f1Row(f1Dash(2, 2, 1), "quartiles")}
    <h5>C · clinical category</h5>
    ${cats}
    ${f1Row(f1Sw.dot(f1Col(d, "Edge Case")), `any patient beyond ${C.cutoff ?? 300} d`)}
    <h5>D · divergence</h5>
    ${f1Row(f1Sw.line(pd, "dot"), "patients data Gaussian KDE")}
    ${f1Row(f1Sw.area(mc), "model Gaussian KDE")}
    ${f1On("f1DotsData") ? f1Row(f1Sw.tick(pd), "one patient sequence") : ""}
    ${f1On("f1DotsModel") ? f1Row(f1Sw.tick(mc), "one model sequence, time-matched") : ""}
    <h5>E · global clock</h5>
    ${clock(d.D, "standard patient sequences", "standard patients – root-to-tip regression",
            () => "model – root-to-tip regression on standard individuals")}
    <h5>F · intra-host clock</h5>
    ${clock(d.F, "long shedder sequences", "intra-host long shedder longitudinal data – root-to-tip regression",
            m => `model ${f1Name(d, m.cohort)} – root-to-tip regression`)}
    <div class="r f1foot"><span>model line: mean of single sim. fits ±1 SD</span></div>
  </aside>`;
}

/* the original's 2×3 grid (fig1_save.py: GridSpec(2, 6), two columns each):
   a panel is letter, title, canvas and key; captions and tables go below */
/* printed letter per panel: the divergence distribution (data key E) sits
   left of the global clock (data key D), so in reading order they swap letters */
const F1_LETTER = {A: "A", B: "B", C: "C", D: "E", E: "D", F: "F"};
const f1Box = (id, title, sub, key) => `<div class="f1p">
  <div class="ph"><h3><b>${F1_LETTER[id]}</b>${f1On("f1Titles") ? " " + title : ""}</h3></div>
  <canvas id="f1_${id}"></canvas>${key}</div>`;
const f1Note = (id, cap, extra = "") => `<div class="f1n"><h4>${id}</h4><div class="cap">${cap}</div>${extra}</div>`;

/* Figure 1 is built on this machine, never on the cluster: it reads patient
   data that stays here (fig1.INPUTS). Missing or stale data starts the build
   (POST /api/rebuild/fig1, one job at a time) and the tab redraws when it ends.
   Another job holding the lock is waited out; any other refusal is an error. */
let f1Job = null;
const f1Sleep = ms => new Promise(r => setTimeout(r, ms));
async function f1Idle() {
  for (;;) { const j = await apiGet("/api/jobs", {fresh: true});
    if (!j.running) return j; await f1Sleep(4000); }
}
function f1Rebuild() {
  if (f1Job) return f1Job;
  f1Job = (async () => {
    for (;;) {
      const r = await apiPost("/api/rebuild/fig1", {});
      if (r.started || /^Figure 1/.test(r.msg || "")) break;
      if (!/already running/.test(r.msg || "")) throw new Error(r.error || r.msg || "the build was refused");
      await f1Idle();
    }
    const j = await f1Idle(), last = (j.history || []).filter(h => /^Figure 1/.test(h.name)).pop();
    if (!last || last.rc !== 0) throw new Error(`the Figure 1 build failed (exit ${last ? last.rc : "?"}): `
      + (j.tail || []).slice(-3).join(" "));
    apiInvalidate("/api/fig1");
  })().finally(() => { f1Job = null; });
  return f1Job;
}

async function renderFig1(main) {
  let d;
  try { d = await apiGetSlot("/api/fig1"); }
  catch (e) {
    if (isSuperseded(e) || !/not computed yet/.test(e.message)) throw e;
    main.innerHTML = `<div class="note">Figure 1 has not been built for this experiment. Building it
      here now with ${f1Code("figures_gui/fig1.py build")}, about 5 minutes; it is built on this
      machine only, because it reads patient data that is not on the cluster.</div>`;
    await f1Rebuild();
    return renderFig1(main);
  }
  f1Last = d;
  const ix = d.index || {}, ex = d.experiment || {}, src = d.sources || {};
  const A = d.A, B = d.B, C = d.C, D = d.D, E = d.E, F = d.F, S = E.series || {};
  const hiv = (d.cohort_scenario || {}).HIV;

  const why = (ix.stale_why || []).map(w => w === "code" ? `the code changed (${ix.code_hash} → ${ix.current_hash})`
    : "the real-world input files changed").join(" and ");
  const stale = (ix.stale ? `<div class="note" id="f1stale">Rebuilding Figure 1 here: ${why || "it is stale"}.
      These panels are the old build until it finishes, about 5 minutes.</div>` : "")
    + ((ix.inputs_missing || []).length ? `<div class="err">Not on this machine:
      ${ix.inputs_missing.map(f => f1Code("Data/RealWorldData/" + f)).join(", ")}. Without them the
      panels that compare with patients have no patient series.</div>` : "");
  if (ix.stale) f1Rebuild().then(() => { if (activeTabId === "fig1") renderTab(tabById("fig1")); },
    e => { const b = document.getElementById("f1stale"); if (b) { b.className = "err"; b.textContent = e.message; } });
  const note = `<div class="note">Experiment ${f1Code(`${ex.name} #${ex.num}`)}, read back off
    ${f1Code("figures_gui/prep/fig1_save.py")}'s own axes by ${f1Code("figures_gui/fig1.py")}
    on ${ix.built || "—"} (code ${ix.code_hash || "—"}). Every number here is one the figure
    draws or prints; the build stops if its clock fits differ from the figure's lines.
    The HIV series is ${f1Code(hiv || "—")} alone: ${f1Code("duration_groups")} keeps one
    scenario per clinical label, so ${f1Code("HIV_high")} is not in this figure.</div>`;

  const half = A.curves.map(c => `${f1Name(d, c.cohort)} ${f1Half50(c) ?? "—"} d`).join(", ");
  const capA = `The within-host model's probability that a host is still infectious,
    by days since infection. Reads no simulation output beyond each scenario's
    durations: ${src.A || ""}. These are the model's design, not a measurement.
    Read off the curves, the first day past its peak each is below one half: ${half}.`;

  const vq = B.violins.map(v => `${f1Name(d, v.cohort)} ${fmt(v.q[1], 0)} d (${fmt(v.q[0], 0)}–${fmt(v.q[2], 0)}, n=${v.n.toLocaleString()})`).join("; ");
  const capB = `Realized infection durations: standard individuals in control,
    long shedders in every other scenario. Source: ${src.B || ""}. Each violin rests on
    at most 150 individuals drawn at random (seed 42, ${f1Code("fig1_plots.violin_sample")}), so a
    scenario with ~280,000 hosts and one with under 200 are drawn from the same number. Width-normalised
    kernel densities cut at the data's range, outlines exactly as seaborn drew them;
    the dashed lines are the quartiles. Median (quartiles, individuals): ${vq}.`;

  const capC = `Published long-shedder shedding durations, one point per patient
    (${C.n} patients, ${C.categories.length} categories), from ${f1Code(src.C || "")},
    column ${f1Code("Days (Exact)")}; categories ordered by mean duration. Boxes are the
    quartiles and 1.5×IQR whiskers handed to matplotlib's ${f1Code("bxp")}; outliers are
    not marked, but every patient is a point. A category is coloured as the scenario
    that stands in for it, grey if none does; any patient beyond ${C.cutoff ?? 300} d
    is coloured as the edge case whatever the category. Jitter is cosmetic.`;

  const clockLine = P => [P.data ? `${P.data_label} ${fmt(P.data.rate, 5)} s/s/y (R² ${fmt(P.data.r2, 2)})` : null]
    .concat(P.model.map(m => `model ${f1Name(d, m.cohort)} ${fmt(m.rate, 5)} ± ${fmt(m.sd_rate, 5)} s/s/y (R² ${fmt(m.r2, 2)})`))
    .filter(Boolean).join("; ");
  const clockCap = (P, dataSrc, modelSrc, what) => `${what} Patient points
    (${P.data ? P.data.n.toLocaleString() : 0}) read ${f1Code(dataSrc)}; their x is
    ${f1Code(P.real_x)} minus its smallest value${P.data ? ` (${fmt(P.data.x_offset, 0)} d)` : ""}, as
    the figure shifts it. The model reads ${modelSrc}. Every line is a regression
    through the origin (${f1Code("fig1_plots._through_origin")}). The model's line is the
    <strong>mean of per-seed slopes</strong>, its band ±1 SD across seeds — the spread of
    single epidemics, since the patient cohort is itself one draw; its R² is of the whole
    cloud pooled. Rates: ${clockLine(P)}.`;
  const capD = clockCap(D, src.D_data || "", src.D_model || "",
    `The global clock: divergence from the root against sampling time, standard
     individuals only.`);
  const capF = clockCap(F, src.F_data || "", src.F_model || "",
    `The intra-host clock: each long shedder's lineages against that host's own founder
     genome and time since its own infection.`);

  const sd_ = S.data || {}, sm = S.model || {};
  const capE = `Panel E's two series with time collapsed: the distribution of divergence.
    Patients: ${(sd_.n || 0).toLocaleString()} sequences from ${f1Code(src.D_data || "")}. Model:
    control's ${(sm.n || 0).toLocaleString()} standard sequences, resampled with replacement to
    ${(sm.matched || 0).toLocaleString()} so their sampling days follow the patients' in
    ${E.n_time_bins} quantile bins (${f1Code("fig1_plots._time_match")}) — otherwise this would
    compare two sampling designs as much as two clocks. Curves are Gaussian KDEs, each bandwidth floored at
    1.2× its own series' step: 1/genome length for the model, 1/median overlap length
    for the patients (${f1Code("_smoothed")}); the faint
    steps are histograms, each binned at twice its own resolution. Read the shape, not the
    rate: E shows the rates; here the patients are bimodal, the model is not.`;


  main.innerHTML = stale + note + `<div class="plot"><div class="f1wrap"><div class="f1grid">`
    + f1Box("A", "Infectiousness by day", "within-host model", "")
    + f1Box("B", "Realized infection duration", "simulated, 10 seeds", "")
    + f1Box("C", "Published shedding durations", `${C.n} patients`, "")
    + f1Box("E", E.title || "Divergence distribution", "", "")
    + f1Box("D", D.title || "Global clock", "patients vs control", "")
    + f1Box("F", F.title || "Intra-host clock", "patients vs long shedders", "")
    + `</div>${f1Legend(d)}</div></div><div class="plot"><div class="ph"><h3>Panel notes</h3>
      <span class="sub">what each panel reads and the numbers it draws</span></div>`
    + f1Note("A", capA) + f1Note("B", capB) + f1Note("C", capC)
    + f1Note("D", capE) + f1Note("E", capD, f1ClockTable(d, D)) + f1Note("F", capF, f1ClockTable(d, F))
    + `</div>`;

  f1PaintA(document.getElementById("f1_A"), d);
  f1PaintB(document.getElementById("f1_B"), d);
  f1PaintC(document.getElementById("f1_C"), d);
  f1PaintE(document.getElementById("f1_E"), d);
  f1PaintClocks();
  f1Align(d);
}

/* Panels side by side share their x axis line: each row takes the largest
   bottom padding any of its canvases needs, then all six repaint with it */
function f1Align(d) {
  const ids = ["A", "B", "C", "D", "E", "F"], rows = {};
  ids.forEach(id => { const cv = document.getElementById("f1_" + id); if (!cv) return;
    const top = Math.round(cv.getBoundingClientRect().top);
    (rows[top] = rows[top] || []).push(cv.id); });
  let moved = false;
  Object.values(rows).forEach(r => { const b = Math.max(...r.map(i => f1PadNeed[i] || 0));
    r.forEach(i => { if (f1PadMin[i] !== b) { f1PadMin[i] = b; moved = true; } }); });
  if (!moved) return;
  f1PaintA(document.getElementById("f1_A"), d);
  f1PaintB(document.getElementById("f1_B"), d);
  f1PaintC(document.getElementById("f1_C"), d);
  f1PaintE(document.getElementById("f1_E"), d);
  f1PaintClocks();
}

registerTab({
  id: "fig1", group: "figures", num: "1", title: "Model", viewId: "viewFig1",
  lede: "The paper's Figure 1, panel for panel: how long hosts stay infectious in the model and in the literature, and the model's global and intra-host clocks against patient data. Drawn from what <code>figures_gui/prep/fig1_save.py</code> itself drew.",
  rail: f1Rail,
  rerenderOn: [],
  render: renderFig1,
});
