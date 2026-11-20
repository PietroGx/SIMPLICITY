/* Method: ported from method.html.

   Static page plus a live state block read from /api/index. Its page-local
   style block moved into studio.css, so this tab carries no CSS of its own.
*/

/* The definition text, from the first of method.html's two script blocks. */


/* The metric definitions, module scope: both the rail's contents list and
   the rendered page read them, so they cannot live inside either. */
const DEFS = `
<div class="sec" id="clade"><h2>Forming a clade</h2>

<div class="def"><h3>Clade membership</h3>
<p>Clade formation is not reimplemented here. The studio calls the pipeline's own
function, the same one <code>tests/export_fig2_builder.py</code> uses, so a clade
in the studio is a clade in the figure.</p>
<pre>c2l, l2c, pcm, meta = cluster_lin_into_clades_with_meta(
                          phylogenetic_data, shared_mut_threshold = k)</pre>
<p class="src"><b>tests/studio/compute.py:136</b> · <b>simplicity/clustering.py</b></p></div>

<div class="def"><h3>The transmission floor</h3>
<p>The second half of the clade definition, and part of the definition rather than
a display filter. A group of lineages is a clade when it has diverged at least
<i>k</i> mutations from its clade root <strong>and</strong> it and its descendants
carry at least <i>m</i><sub>T</sub> of the run's transmission events.</p>
<pre>n_tx_run    = number of individuals carrying an inherited_lineage
ntx(clade)  = sum of that count over the clade's member lineages
floor       = (m_T / 100) * n_tx_run
keep        = ntx >= floor</pre>
<p><i>m</i><sub>T</sub> is a percentage of that run's own transmissions, never an
absolute count: a fixed count would be a different filter in each scenario. Kept
regardless of the floor: the founding clade. Dropped: any group that has a parent
and no defining mutations, which is not a clade.</p>
<p class="src"><b>tests/studio/compute.py:142-143</b> · <b>tests/studio/panels.py:154-168</b></p></div>

<div class="def"><h3>Frequency frame</h3>
<p>Times down, clades across. A clade's frequency on a day is the sum of its
member lineages' frequencies that day.</p>
<pre>F[t, clade] = sum over the clade's lineages of Frequency_at_t</pre>
<p>Read at <strong>daily</strong> resolution from the raw run output. The
clade-explorer export is downsampled to 120 time points about 9 days apart, and
rise times land at a few tens of days, so on that axis a whole rise would be three
or four samples.</p>
<p class="src"><b>tests/studio/compute.py:148-162</b>, mirroring
<b>tests/export_fig2_builder.py:190-209</b></p></div>
</div>

<div class="sec" id="burden"><h2>Burden and composition — panel 1</h2>

<div class="def"><h3>Burden <span class="u">days</span></h3>
<p>A clade's frequency integrated over time. Frequency is a proportion and so
carries no unit; integrating it over the day axis leaves days. It reads as how
many days' worth of the entire viral population the clade accounts for, so 100%
for 10 days and 10% for 100 days both score 10 days.</p>
<pre>burden = integral of F[t, clade] dt
       = np.trapezoid(y, t)        # trapezoidal rule on the day axis</pre>
<p>Two values are stored per clade: <code>burden</code> over the post-burn-in
window and <code>burden_all</code> over the whole record. The burn-in checkbox
picks which one every panel-1 percentage is formed from.</p>
<p class="src"><b>tests/studio/compute.py:186-188</b></p></div>

<div class="def"><h3>Composition class</h3>
<p>Clades grouped by how many of their defining mutations arose in a long-shedder
host. Two classes are not a count at all and are labelled as themselves.</p>
<pre>nlong = defining mutations whose host_type starts with "long"

class = "wt"        if the clade has no parent
        "other"     if n_defining != k
        str(nlong)  otherwise</pre>
<p><strong>wt</strong> is the founding clade. Its root lineage is named
<code>wt</code> and carries no mutations, so its long-born ratio is 0/0, which is
not zero but has no value at all — it is labelled rather than merged into class 0.
<strong>other</strong> holds clades whose defining-mutation count is not <i>k</i>:
<code>n_defining</code> counts mutation events while <i>k</i> thresholds genome
divergence, so the segment can overshoot and a handful carry 8 at <i>k</i>=7. That
is the clustering rule working, not a defect. Together with 0…<i>k</i> the two
account for every clade, which is why the classes sum to 100% of a run.</p>
<p class="src"><b>tests/studio/compute.py:180-182</b> · <b>tests/studio/panels.py:138-152</b></p></div>

<div class="def"><h3>The three measures panel 1 plots</h3>
<p>Each is formed <strong>inside one run</strong>, then averaged over runs.</p>
<pre>frac_burden[class] = burden(class) / burden(all kept clades)   [fraction]
frac_clades[class] = count(class)  / count(all kept clades)    [fraction]
per_clade[class]   = burden(class) / count(class)              [days]</pre>
<p>Every share is stored and served as a fraction, 0 to 1, and multiplied by 100
only where the studio draws or prints it.</p>
<p>The three are one decomposition: burden percentage is a product of how many
clades a class holds and how big each one is, and those can pull in opposite
directions. A class holding far more of the burden than of the count is made of
many small clades, not individually fitter ones — which is why all three are on
the page rather than the first alone.</p>
<p class="src"><b>tests/studio/panels.py:205-207</b></p></div>
</div>

<div class="sec" id="tx"><h2>Transmission events — panel 2</h2>

<div class="def"><h3>The 2x2 table</h3>
<p>One transmission event is one infection: a transmitter, a recipient, a day and
the lineage passed on. Each event falls in exactly one cell, by the host type of
each side.</p>
<pre>ll  long     -> long        ls  long     -> standard
sl  standard -> long        ss  standard -> standard</pre>
<p>The initial seeded infections have <code>parent</code> = <code>root</code>: they
have no transmitter, so they are not transmission events. They are counted
separately per run and excluded from every percentage.</p>
<p class="src"><b>tests/export_transmissions.py:35,66-90</b></p></div>

<div class="def"><h3>The three claims, from that one table</h3>
<pre>who_transmits      = (ll + ls) / (ll + ls + sl + ss)     [fraction]
where              = cell       / (ll + ls + sl + ss)     [fraction]
intensity_long     = (ll + ls) / n_long_hosts          [events per host]
intensity_standard = (sl + ss) / n_standard_hosts      [events per host]</pre>
<p>Read these against the config before reading them as a result.
<code>R_long</code> equals the standard reproduction number and
<code>susceptibility_long</code> is 1.0 in every scenario — long shedders have no
transmission advantage and are no easier to infect by construction; only shedding
duration and mutation rate differ. So <code>who_transmits</code> and the two
intensities are manipulation checks with answers the config already implies. The
<code>where</code> row is the one that can say something the config does not.</p>
<p class="src"><b>tests/studio/panels.py:260-268</b> ·
<b>scripts/experiments/impact_long_shedders_config.py:105-109</b></p></div>

<div class="def"><h3>Why this table is trusted</h3>
<p>It is assembled twice per run — once from each individual's <code>parent</code>
(who infected me) and once from each individual's <code>new_infections</code> (who
I infected) — and the run fails unless the two agree on every cell and every day.
The totals must also close: events plus seeded equals the number infected.</p>
<p class="src"><b>tests/export_transmissions.py:92-99</b></p></div>
</div>

<div class="sec" id="epi"><h2>Episodes, rise and latency — panel 3</h2>

<div class="def"><h3>Smoothing <span class="u">days</span></h3>
<p>A centred moving average whose window is given in days and converted to samples
with the mean step. The edges shrink rather than wrap, so no value is invented
beyond the record.</p>
<pre>step  = (t[-1] - t[0]) / (n - 1)
w     = round(win_days / step), forced odd, at least 1
ys[i] = mean of y over [i - w//2, i + w//2], clipped to the record</pre>
<p class="src"><b>tests/studio/compute.py:59-75</b></p></div>

<div class="def"><h3>Episode</h3>
<p>One rise and fall: a local maximum whose prominence clears a bar set as a
fraction of the clade's own peak height, starting at the local minimum
<strong>directly before it</strong>.</p>
<pre>peaks = find_peaks(ys, prominence = prom * max(ys))
mins  = find_peaks(-ys), with index 0 prepended
start(peak i) = the last element of mins below i
amplitude     = ys[peak] - ys[start]</pre>
<p>The start is deliberately not the prominence base. The base tunnels left to
wherever the series was last that low, which on almost every clade is the start of
the record — it would make every episode begin on day one.</p>
<p class="src"><b>tests/studio/compute.py:77-94</b></p></div>

<div class="def"><h3>Emergence episode</h3>
<p>A clade can rise more than once, and the question the panel asks is about the
sweep that mattered, not the first stirring. The emergence episode is the one
containing the clade's highest frequency.</p>
<pre>imax      = argmax(ys)
emergence = the episode whose [start, peak] contains imax
            (ties broken by the larger amplitude)</pre>
<p class="src"><b>tests/studio/compute.py:210-211</b></p></div>

<div class="def"><h3>Emergence rise time <span class="u">days</span></h3>
<p>How fast that sweep happened: the days spent crossing the middle of its own
climb. The band is a fraction of <strong>this episode's</strong> amplitude, not of
the clade's peak, so a late episode on a clade that was once far higher is still
measured against its own rise.</p>
<pre>ylo, yhi = ys[start], ys[peak]
A  = ylo + lo * (yhi - ylo)        B  = ylo + hi * (yhi - ylo)
ia = first index in [start, peak] with ys >= A
ib = first index in [start, peak] with ys >= B
rise_days = t[ib] - t[ia]</pre>
<p class="src"><b>tests/studio/compute.py:96-114</b></p></div>

<div class="def"><h3>Latent interval <span class="u">days</span></h3>
<p>How long the clade sat around before that sweep began.</p>
<pre>latent_days = start_day(emergence episode) - first_appearance_day</pre>
<p><code>first_appearance_day</code> is the clade's <code>start_time</code> from the
clustering metadata. Reported only for clades that are not truncated, and only
when it comes out at or above zero.</p>
<p class="src"><b>tests/studio/panels.py:318-323</b></p></div>

<div class="def"><h3>Truncated</h3>
<p>A clade already above a trace frequency on the burn-in day. Its first
appearance is outside the window, so its latency cannot be measured — but its
rise, which happens inside the window, still can.</p>
<pre>trunc = phi exists  and  start_time &lt; phi  and  F[phi, clade] &gt; 1e-4</pre>
<p>Such a clade keeps its rise time and loses its latent interval. Rise and
latency therefore have <strong>different denominators</strong>, and panel 3 reports
both counts alongside the number truncated rather than quietly dividing by the
same n twice.</p>
<p class="src"><b>tests/studio/compute.py:190-192</b></p></div>

<div class="def"><h3>Burn-in day <span class="u">φ</span></h3>
<p>The first day on which <code>(diagnosed + recovered) / N</code> reaches 1, read
per run from that run's own trajectory. Some runs never reach it and have no φ.</p>
<p>It matters because selection strength scales with the fraction of the population
already infected and reaches full weight there. Before φ a clade's trajectory is
substantially drift; after it, the fitness term is doing work. The window control
on the Clustering page switches every frequency measure between the post-φ window
and the whole record — and the <code>wt</code> class is the one that moves most,
which is the measurement that would have been invisible had the founding clade
been dropped.</p>
<p class="src"><b>scripts/nature_plots/fig2_preprocess_data.py:344</b> ·
<b>simplicity/phenotype/update.py:39-42</b></p></div>
</div>

<div class="sec" id="agg"><h2>Aggregating over runs</h2>

<div class="def"><h3>Interval</h3>
<p>A percentile bootstrap over runs. The per-run values go in; the interval comes
from resampling those run values, never the clades inside them.</p>
<pre>v      = the per-run values, one number per run
means  = 2000 resamples of v, with replacement, each of size len(v),
         each reduced to its mean
mean   = v.mean()
lo, hi = 2.5th and 97.5th percentile of means</pre>
<p>The resampler is seeded (20261116), so reloading a page gives the same interval
rather than a slightly different one each time. Below three runs no interval is
quoted at all — with one or two there is nothing to quote and a number would
invent precision.</p>
<p class="src"><b>tests/studio/panels.py:22-23,123-134</b></p></div>
</div>

<div class="sec" id="rules"><h2>Two rules that hold everywhere</h2>

<div class="def"><h3>The run is the independent unit</h3>
<p>Every percentage is formed inside a run and then averaged across runs, and every
interval bootstraps over runs. A run that happens to produce 3000 clades does not
outvote one that produces 300. Pooling clades across runs would treat each clade as
an independent observation, which it is not: clades within a run share a history,
a burn-in and a population.</p></div>

<div class="def"><h3>Scenarios are never pooled</h3>
<p>The five scenarios are a designed pair of contrasts. SOT against HIV low holds
long-shedder prevalence at 1.0% and moves shedding duration from 63 to 109 days;
HIV low against HIV high holds duration at 109 days and moves prevalence from 1.0%
to 12%. A statistic averaged across the five destroys both contrasts and answers
nothing the experiment was built to ask.</p>
<p class="src"><b>scripts/experiments/impact_long_shedders_config.py:105-109</b></p></div>
</div>
`;

function rail_method(rail) {
const toc = document.createElement("div");
toc.className = "card";
toc.innerHTML = `<h2>Contents</h2><div class="toc">
  <a href="#state">What is precomputed</a>
  <a href="#flow">The four passes</a>
  <a href="#clade">Forming a clade</a>
  <a href="#burden">Burden and composition</a>
  <a href="#tx">Transmission events</a>
  <a href="#epi">Episodes, rise, latency</a>
  <a href="#agg">Aggregating over runs</a>
  <a href="#rules">Two rules that hold everywhere</a>
</div>`;
rail.append(toc);

}

async function render_method(main){

  main.innerHTML = `<div class="loading">loading…</div>`;
  let ix = null, note = "";
  try { ix = await apiGetSlot("/api/index"); }
  catch(e){ note = e.message; }

  let state = "", flow = "";
  if (ix){
    const b3 = ix.builds || [];
    const tx = ix.transmissions;
    const runs = b3.length ? Object.values(b3[0].counts || {}).reduce((a,c)=>a+c,0) : 0;
    const txruns = tx ? (tx.scenarios || []).length : 0;
    state = `<div class="stats">
      <div class="stat"><div class="k">clade passes cached</div>
        <div class="v">${b3.length}</div></div>
      <div class="stat"><div class="k">runs per pass</div>
        <div class="v">${runs} <small>of 250</small></div></div>
      <div class="stat"><div class="k">transmission export</div>
        <div class="v">${tx ? txruns + " " : "—"}<small>${tx ? "scenarios" : "missing"}</small></div></div>
      <div class="stat"><div class="k">compute code</div>
        <div class="v" style="font-size:.78rem">${ix.code}</div></div>
    </div>`;
  } else {
    state = `<div class="err">could not read the build state: ${note}</div>`;
  }

  /* Paths come from /api/index, never a literal: this diagram named
     impact_long_shedders_#1 for a day after the studio moved to #2. */
  const store = ix ? ix.store : "the store";
  const src = (ix && ix.experiment) ? ix.experiment.source : "the run output";

  flow = `<div class="flow">
    <div class="step cold"><div class="hd"><h3>1 · Raw run output</h3>
      <span class="when">read only during a rebuild</span></div>
      <p>The simulation's own output, one directory per seed. Nothing in the studio
      reads this while you click: it is opened only by a rebuild, which is why a
      rebuild takes minutes and a page takes milliseconds.</p>
      <div class="io"><span>reads</span><span>${src}</span>
      <span>per run</span><span>phylogenetic_data · lineage_frequency · individuals_data</span></div></div>

    <div class="step"><div class="hd"><h3>2 · Transmission export</h3>
      <span class="when">about 10 min · no phylogeny, no clustering</span></div>
      <p>Counts every transmission event into the 2x2 table, per day, per run.
      Independent of <i>k</i> and of the episode rules, so it is built once and
      never rebuilt when you move a slider.</p>
      <div class="io"><span>runs</span><span>tests/export_transmissions.py</span>
      <span>writes</span><span>${store}/tx/</span>
      <span>feeds</span><span>panel 2</span></div></div>

    <div class="step"><div class="hd"><h3>3 · Clade pass</h3>
      <span class="when">about 12 min per parameter set</span></div>
      <p>Clusters each run's lineages at one <i>k</i>, builds the daily frequency
      frame, and records every clade's composition, burden and episodes. One pass
      per parameter set; the cache key carries the parameters <strong>and a hash of
      the computing code</strong>, so a change to the maths cannot be served out of
      a stale build.</p>
      <div class="io"><span>runs</span><span>tests/studio/compute.py</span>
      <span>writes</span><span>${store}/cache/&lt;key&gt;/l3_&lt;scenario&gt;.json</span>
      <span>key</span><span>k · smoothing window · prominence · rise band · window · code hash</span>
      <span>feeds</span><span>the Frequency page and panels 1, 3, 4</span></div></div>

    <div class="step"><div class="hd"><h3>4 · Panel aggregation</h3>
      <span class="when">milliseconds warm · a few seconds cold</span></div>
      <p>Applies the transmission floor, forms each percentage inside a run,
      bootstraps over runs, and answers the page. Cached on disk too, keyed by the
      request, the code hash and <strong>which clade passes exist</strong> — panel 4
      reads every pass, so adding one has to invalidate it.</p>
      <div class="io"><span>runs</span><span>tests/studio/panels.py, in the server process</span>
      <span>writes</span><span>${store}/cache/_panels/</span>
      <span>serves</span><span>/api/panel · /api/episodes · /api/index</span></div></div>
  </div>
  <div class="note">Only pass 1 touches the simulation output. Passes 2 and 3 are
    jobs you start from the Export page — one at a time, refused rather than queued,
    and reachable only over the loopback interface. Everything a panel draws comes
    off disk.</div>`;

  main.innerHTML = `
    <div class="sec" id="state"><h2>What is precomputed</h2>
      <p class="lede" style="margin-bottom:10px">Live, from this server's own store
        — not a note that can go stale.</p>
      ${state}</div>
    <div class="sec" id="flow"><h2>The four passes</h2>${flow}</div>
    ${DEFS}`;

}

registerTab({
  id: "method", num: 0, title: "Method", viewId: "viewMethod",
  lede: "The flow from a run's raw output to a panel, and every metric with the calculation that produces it. Each definition names the file and line that computes it, so a definition here can be checked against the code rather than trusted.",
  rail: rail_method,
  render: render_method,
});
