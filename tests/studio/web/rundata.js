/* One run's data, decoded once.

   The server passes a per-seed file through `/api/run`: packed base64 arrays
   holding the tree topology (`parent`, `depth`, `nedge`), the host type of each
   lineage, the clade starts at every one of the k values, and the per-lineage
   frequency trajectories. `decodeRun` turns those into typed arrays.

   `runView` is the reason this file exists. The tree, the Muller plot and the
   episode grid all need the same three derived things -- the topology layout,
   the clade assignment, and the colour map -- and the artifact recomputed them
   in each painter. Here they are computed once per (run, k, m_T) and shared, so
   changing a VIEW setting such as the tree's x axis redraws without touching
   any of it.

   `computeClades` is rule B in the browser, ported unchanged from
   tests/artifacts/fig2_studio.html. Its own comment records the check that
   matters: it agrees with the python designation over 90 scenario x seed x k x
   m combinations with 0 mismatches, which is why the tree's colouring and the
   clade build's numbers mean the same thing.
*/

const RUN_CACHE_MAX = 6;             // ~218 KB each; six is a comfortable budget
const runCache = new Map();
const runViewCache = new Map();

function decodeRun(o) {
  return {
    n: o.n, scenario: o.scenario, seed: o.seed,
    parent: dec(o.parent, Int32Array), host: dec(o.host, Uint8Array),
    nedge: dec(o.nedge, Uint8Array), depth: dec(o.depth, Uint16Array),
    tem: dec(o.tem, Uint16Array), hosts: dec(o.hosts, Uint32Array),
    starts: o.starts, thresholds: o.thresholds,
    times: dec(o.times, Uint16Array), f_start: dec(o.f_start, Int32Array),
    f_len: dec(o.f_len, Uint16Array), f_vals: dec(o.f_vals, Uint16Array),
    phi: o.phi === undefined ? null : o.phi,
  };
}

async function loadRun(scenario, seedIndex) {
  const key = `${scenario}/${seedIndex}`;
  if (runCache.has(key)) return runCache.get(key);
  const raw = await apiGet(`/api/run?scenario=${encodeURIComponent(scenario)}&seed=${seedIndex}`,
                           {slot: "run", label: `loading ${scenario} seed ${seedIndex}`});
  const D = decodeRun(raw);
  if (runCache.size >= RUN_CACHE_MAX) runCache.delete(runCache.keys().next().value);
  runCache.set(key, D);
  return D;
}

/* The decoded run plus everything derived from the clade definition, memoised
   on exactly the inputs that change it. */
async function runView(scenario, seedIndex, k, mt) {
  const key = `${scenario}/${seedIndex}/k${k}/m${mt}`;
  if (runViewCache.has(key)) return runViewCache.get(key);
  const D = await loadRun(scenario, seedIndex);
  if (!D.starts || !(String(k) in D.starts))
    throw new Error(`this simulation carries no clade starts at k=${k} `
                  + `(it has ${Object.keys(D.starts || {}).join(", ")})`);
  const T = buildTopology(D);
  const C = computeClades(D, T, k, mt);
  const view = {D, T, C};
  if (runViewCache.size > RUN_CACHE_MAX * 2)
    runViewCache.delete(runViewCache.keys().next().value);
  runViewCache.set(key, view);
  return view;
}

/* The burn-in day for the run on screen, or null when the record has none.
   Every chart that shades a burn-in asks this one function. */
function burninDayOf(D) {
  return (D && D.phi !== null && D.phi !== undefined && isFinite(D.phi)) ? D.phi : null;
}

function dec(b64, Type){
  const bin = atob(b64), len = bin.length, buf = new ArrayBuffer(len), v = new Uint8Array(buf);
  for (let i = 0; i < len; i++) v[i] = bin.charCodeAt(i);
  return new Type(buf);
}

function buildTopology(D){
  const n = D.n, par = D.parent;
  const kids = new Int32Array(n).fill(-1), next = new Int32Array(n).fill(-1);
  const last = new Int32Array(n).fill(-1);
  let root = 0;
  for (let i = 0; i < n; i++){
    const p = par[i];
    if (p < 0){ root = i; continue; }
    if (kids[p] < 0) kids[p] = i; else next[last[p]] = i;
    last[p] = i;
  }
  // DFS order and tip-based y layout
  const order = new Int32Array(n); let oi = 0;
  const stack = [root];
  while (stack.length){
    const v = stack.pop(); order[oi++] = v;
    const buf = [];
    for (let c = kids[v]; c >= 0; c = next[c]) buf.push(c);
    for (let j = buf.length - 1; j >= 0; j--) stack.push(buf[j]);
  }
  const y = new Float32Array(n); let tip = 0;
  for (let j = n - 1; j >= 0; j--){
    const v = order[j];
    if (kids[v] < 0){ y[v] = tip++; }
    else { let s = 0, c2 = 0;
      for (let c = kids[v]; c >= 0; c = next[c]){ s += y[c]; c2++; }
      y[v] = s / c2; }
  }
  const rank = new Int32Array(n);
  for (let j = 0; j < n; j++) rank[order[j]] = j;
  return {root, kids, next, order, y, nTips: tip, rank};
}

function computeClades(D, T, k, mt){
  const n = D.n, par = D.parent;
  const starts = dec(D.starts[String(k)], Uint32Array);

  // Rule B: diverging by k is not enough. A lineage founds a clade only if it
  // and its descendants also carry >= m_T of the run's transmission events.
  // Cumulative subtree transmissions fall monotonically down the tree, so a
  // refusal is hereditary: the surviving starts are exactly the unfiltered
  // starts that clear the minimum, checked against the python designation over
  // 90 scenario x seed x k x m combinations with 0 mismatches. A refused
  // lineage stays in its ancestor's clade and brings its divergence with it,
  // which is why the k-1 bound does not survive m_T.
  // m_T is a FRACTION of the run's own transmissions, the unit settings.mt
  // holds since v2.5.15 -- an absolute count would be a different filter in
  // each scenario, which run 10k to 29k transmissions. This read `mTxPct / 100`
  // against a parameter every caller fed settings.mt, so the floor was a
  // hundredth of what the rail asked for: 0.01% where 1% was shown, and the
  // browser designated 59 clades where the build designated 16.
  const stx = new Float64Array(n);
  let totTx = 0;
  for (let i = 0; i < n; i++){ stx[i] = D.hosts[i]; totTx += D.hosts[i]; }
  for (let i = n - 1; i >= 0; i--){ const p = par[i]; if (p >= 0) stx[p] += stx[i]; }
  const mTx = mt > 0 ? mt * totTx : 0;

  const isStart = new Uint8Array(n);
  const nBase = starts.length;
  for (let i = 0; i < starts.length; i++){
    const v = starts[i];
    if (mTx <= 0 || par[v] < 0 || stx[v] >= mTx) isStart[v] = 1;
  }

  // clade id == index of its root lineage; parents precede children, so one pass
  const base = new Int32Array(n);
  for (let i = 0; i < n; i++){
    const p = par[i];
    base[i] = (isStart[i] || p < 0) ? i : base[p];
  }
  // clade tree + sizes
  const roots = [];
  for (let i = 0; i < n; i++) if (base[i] === i) roots.push(i);
  const pClade = new Map();
  for (const r of roots) pClade.set(r, par[r] >= 0 ? base[par[r]] : -1);
  // two units of account: member lineages, and the individuals they infected
  const szL = new Map(), szH = new Map();
  for (let i = 0; i < n; i++){
    const c = base[i];
    szL.set(c, (szL.get(c) || 0) + 1);
    szH.set(c, (szH.get(c) || 0) + D.hosts[i]);
  }

  // m_T was applied at designation, so every surviving root IS a clade:
  // nothing to fold, and a lineage's clade is simply its base.
  const alive = new Set(roots);
  const L = szL, H = szH;
  const clade = base;
  const pAlive = pClade;

  // defining mutations: walk from the clade root up to its parent clade's root
  const host = D.host, nedge = D.nedge;
  const info = new Map();
  for (const r of alive){
    const pr = pAlive.get(r);
    let nl = 0, ns = 0, v = r, steps = 0;
    while (v >= 0 && steps++ < 100000){
      // the parent clade's root edge already defines the PARENT; counting it here
      // too would let one mutation vote in two clades' ballots. Excluded always.
      const atParentRoot = (v === pr);
      if (!atParentRoot){
        if (host[v]) nl += nedge[v]; else ns += nedge[v];
      }
      if (atParentRoot) break;
      v = par[v];
    }
    const tot = nl + ns;
    let lab;
    if (pr < 0 && tot === 0) lab = "founder";
    else if (nl > ns) lab = "long";
    else if (ns > nl) lab = "standard";
    else lab = "mixed";
    info.set(r, {nl, ns, p: tot ? nl / tot : 0, vote: lab,
                 root: host[r] ? "long" : "standard",
                 size: L.get(r), hosts: H.get(r), isRoot: pr < 0, ndef: tot});
  }
  return {clade, alive, info, pAlive, nBase};
}

function dayToIndex(times, day){
  const nt = times.length;
  if (nt < 2 || day < times[0] || day > times[nt - 1]) return null;
  let lo = 0, hi = nt - 1;
  while (hi - lo > 1){ const mid = (lo + hi) >> 1; if (times[mid] <= day) lo = mid; else hi = mid; }
  const span = times[hi] - times[lo];
  return span > 0 ? lo + (day - times[lo]) / span : lo;
}

function cladeSeries(D, C, ids){
  const times = D.times, nt = times.length;
  const idxOf = new Map(); ids.forEach((r, i) => idxOf.set(r, i));
  const F = new Float32Array(ids.length * nt);
  const fs = D.f_start, fl = D.f_len, fv = D.f_vals;
  let k = 0;
  for (let i = 0; i < D.n; i++){
    const len = fl[i];
    if (!len) continue;
    const ci = idxOf.get(C.clade[i]);
    const t0 = fs[i];
    if (ci !== undefined) for (let j = 0; j < len; j++) F[ci * nt + t0 + j] += fv[k + j] / 65535;
    k += len;
  }
  return F;
}