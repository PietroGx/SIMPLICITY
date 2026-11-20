/* Figure 1 · Model: the six panels of scripts/nature_plots/fig1_save.py, drawn
   from what that script drew. tests/studio/fig1.py runs build_figure_1 and
   records each painter's frames, limits, outlines and fits, and stops if its
   clock fits differ from the lines the figure drew; nothing here computes a
   clock, a fit or a density. Canvas ids are f1_*. */

API_LABELS["/api/fig1"] = "reading Figure 1";

const F1_H = 240;
let f1Cloud = false, f1Last = null;

const f1Data = () => css("--fig1-data");
const f1Col = (d, cohort) => scenarioRaw((d.cohort_scenario || {})[cohort] || "control");
const f1Role = (d, role) => role === "unmodelled" ? css("--fig1-unmod") : f1Col(d, role);
const f1Name = (d, c) => { const sc = (d.cohort_scenario || {})[c]; return sc ? SCENARIO_LABEL[sc] : c; };
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
  const pad = {l: 46, r: 12, t: 10,
               b: !cats ? 34 : tilt ? Math.ceil(16 + (sub ? 12 : 0) + catW * 0.66) : (sub ? 40 : 30)};
  axes(g, w, h, pad, xlab, ylab, y0, y1, fmtY);
  const X = v => pad.l + (w - pad.l - pad.r) * (v - x0) / ((x1 - x0) || 1);
  const Y = v => h - pad.b - (h - pad.b - pad.t) * (v - y0) / ((y1 - y0) || 1);
  g.strokeStyle = css("--line"); g.textAlign = "center";
  const ticks = cats ? cats.map((c, i) => [i, c]) : f1Ticks(x0, x1).map(v => [v, fmtX ? fmtX(v) : String(v)]);
  ticks.forEach(([v, t], i) => {
    const x = X(v);
    g.beginPath(); g.moveTo(x, h - pad.b); g.lineTo(x, h - pad.b + 3); g.stroke();
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
    g.fillText(t, x, h - pad.b + 13);
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
    cats: V.map(v => f1Name(d, v.cohort)), sub: V.map(v => `n=${v.n.toLocaleString()}`)});
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
      math: `individuals  ${v.n.toLocaleString()}\nquartiles    ${v.q.map(q => fmt(q, 1)).join(" · ")} d\nmean         ${fmt(v.mean, 1)} d`}});
  });
  g.restore();
  cv._ptHit = hits; axisTip(cv);
}

/* ------------------------------------------------------------ C */
function f1PaintC(cv, d) {
  const C = d.C, K = C.categories;
  const {g, X, Y} = f1Frame(cv, {x0: C.xlim[0], x1: C.xlim[1], y0: C.ylim[0], y1: C.ylim[1],
    ylab: "literature duration (days)", fmtY: v => v.toFixed(0),
    cats: K.map(k => k.cat), sub: K.map(k => `n=${k.n}`)});
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
  if (f1Cloud) P.model.forEach(m => {
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
    xlab: "divergence (10⁻³ subs/site)", ylab: "density (10³ per subs/site)",
    fmtX: v => fmt(v * 1e3, 1), fmtY: v => v.toFixed(2)});
  const col = k => k === "data" ? f1Data() : f1Col(d, "Control");
  ["model", "data"].forEach(k => {
    const h = (S[k] || {}).hist;
    if (!h || !h.x.length) return;
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
    g.strokeStyle = col("data"); g.lineWidth = 1.6; g.setLineDash([6, 3]);
    f1Path(g, kd.x, kd.y, X, Y); g.stroke(); g.setLineDash([]);
  }
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
  const el = document.createElement("div");
  el.className = "ctl";
  el.innerHTML = `<label class="chk"><input type="checkbox"${f1Cloud ? " checked" : ""}>
      <span>Model points in D and F</span></label>
    <div class="hint">Figure 1 draws only the model's fitted lines: its thousands of
      points buried the patients'. Ticked, each scenario's points are added behind,
      thinned to 600 with <code>fig1_plots._thin</code> as
      <code>plot_clock_overlay(model_scatter=True)</code> would draw them.</div>`;
  el.querySelector("input").addEventListener("change", ev => {
    f1Cloud = !!ev.target.checked; f1PaintClocks(); });
  rail.append(card("Display", [el]));
}

/* the original's 2×3 grid (fig1_save.py: GridSpec(2, 6), two columns each):
   a panel is letter, title, canvas and key; captions and tables go below */
const f1Box = (id, title, sub, key) => `<div class="f1p">
  <div class="ph"><h3><b>${id}</b> ${title}</h3><span class="sub">${sub}</span></div>
  <canvas id="f1_${id}"></canvas>${key}</div>`;
const f1Note = (id, cap, extra = "") => `<div class="f1n"><h4>${id}</h4><div class="cap">${cap}</div>${extra}</div>`;

async function renderFig1(main) {
  const d = await apiGetSlot("/api/fig1");
  f1Last = d;
  const ix = d.index || {}, ex = d.experiment || {}, src = d.sources || {};
  const cohortKey = cs => keyRow(cs.map(c => [f1Col(d, c), f1Name(d, c)]));
  const A = d.A, B = d.B, C = d.C, D = d.D, E = d.E, F = d.F, S = E.series || {};
  const hiv = (d.cohort_scenario || {}).HIV;

  const stale = ix.stale ? `<div class="err">Built by code that has since changed
    (hash ${ix.code_hash}, now ${ix.current_hash}): these panels may no longer be the
    figure. Rebuild with ${f1Code(`python tests/studio/fig1.py build --exp-num ${ex.num} --store &lt;store&gt;`)}.</div>` : "";
  const note = `<div class="note">Experiment ${f1Code(`${ex.name} #${ex.num}`)}, read back off
    ${f1Code("scripts/nature_plots/fig1_save.py")}'s own axes by ${f1Code("tests/studio/fig1.py")}
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
    long shedders in every other scenario. Source: ${src.B || ""}. Width-normalised
    kernel densities cut at the data's range, outlines exactly as seaborn drew them;
    the dashed lines are the quartiles. Median (quartiles, individuals): ${vq}.`;

  const capC = `Published long-shedder shedding durations, one point per patient
    (${C.n} patients, ${C.categories.length} categories), from ${f1Code(src.C || "")},
    column ${f1Code("Days (Exact)")}; categories ordered by mean duration. Boxes are the
    quartiles and 1.5×IQR whiskers handed to matplotlib's ${f1Code("bxp")}; outliers are
    not marked, but every patient is a point. A category is coloured as the scenario
    that stands in for it, grey if none does; any patient beyond ${C.cutoff ?? 300} d
    is coloured as the edge case whatever the category. Jitter is cosmetic.`;
  const keyC = keyRow([["SOT", "SOT"], ["HIV", "HIV"]].filter(([c]) => (d.cohort_scenario || {})[c])
    .map(([c, t]) => [f1Col(d, c), t]).concat([[f1Col(d, "Edge Case"), `beyond ${C.cutoff ?? 300} d`],
    [css("--fig1-unmod"), "no modelled scenario"]]));

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
  const capE = `Panel D's two series with time collapsed: the distribution of divergence.
    Patients: ${(sd_.n || 0).toLocaleString()} sequences from ${f1Code(src.D_data || "")}. Model:
    control's ${(sm.n || 0).toLocaleString()} standard sequences, resampled with replacement to
    ${(sm.matched || 0).toLocaleString()} so their sampling days follow the patients' in
    ${E.n_time_bins} quantile bins (${f1Code("fig1_plots._time_match")}) — otherwise this would
    compare two sampling designs as much as two clocks. Curves are Gaussian KDEs with the
    bandwidth floored at 1.2× the model's 1/genome step (${f1Code("_smoothed")}); the faint
    steps are histograms, each binned at twice its own resolution. Read the shape, not the
    rate: D shows the rates; here the patients are bimodal, the model is not.`;

  const keyClock = P => keyRow([[f1Data(), `${P.data_label} (dashed)`]]
    .concat(P.model.map(m => [f1Col(d, m.cohort), `model, ${f1Name(d, m.cohort)}`])));

  main.innerHTML = stale + note + `<div class="plot"><div class="f1grid">`
    + f1Box("A", "Infectiousness by day", "within-host model", cohortKey(A.curves.map(c => c.cohort)))
    + f1Box("B", "Realized infection duration", "simulated, 10 seeds", cohortKey(B.violins.map(v => v.cohort)))
    + f1Box("C", "Published shedding durations", `${C.n} patients`, keyC)
    + f1Box("D", D.title || "Global clock", "patients vs control", keyClock(D))
    + f1Box("E", E.title || "Divergence distribution", "D's points, time collapsed",
      keyRow([[f1Data(), `${sd_.label || "Patient data"} (dashed)`], [f1Col(d, "Control"), "model, control"]]))
    + f1Box("F", F.title || "Intra-host clock", "patients vs long shedders", keyClock(F))
    + `</div></div><div class="plot"><div class="ph"><h3>Panel notes</h3>
      <span class="sub">what each panel reads and the numbers it draws</span></div>`
    + f1Note("A", capA) + f1Note("B", capB) + f1Note("C", capC)
    + f1Note("D", capD, f1ClockTable(d, D)) + f1Note("E", capE) + f1Note("F", capF, f1ClockTable(d, F))
    + `</div>`;

  f1PaintA(document.getElementById("f1_A"), d);
  f1PaintB(document.getElementById("f1_B"), d);
  f1PaintC(document.getElementById("f1_C"), d);
  f1PaintE(document.getElementById("f1_E"), d);
  f1PaintClocks();
}

registerTab({
  id: "fig1", group: "figures", num: "1", title: "Figure 1 · Model", viewId: "viewFig1",
  lede: "The paper's Figure 1, panel for panel: how long hosts stay infectious in the model and in the literature, and the model's global and intra-host clocks against patient data. Drawn from what <code>scripts/nature_plots/fig1_save.py</code> itself drew.",
  rail: f1Rail,
  rerenderOn: [],
  render: renderFig1,
});
