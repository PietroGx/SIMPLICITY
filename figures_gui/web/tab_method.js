/* Method: the maths behind every plotted quantity, figure by figure.

   Each entry is one quantity: a line saying what is drawn, its formula
   (KaTeX, which index.html loads), its symbols, and the file:line the formula
   mirrors. Written from the code that computes the number, not from captions.
   Above the figures, a live block read from /api/index and the stores a
   number passes through. */

const mthT = String.raw;
const mthEsc = s => String(s).replace(/[&<>]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;"})[c]);

/* One entry: [what, [display formulas], [[symbol tex, meaning html], ...], source].
   Shared definitions, then one section per figure, each a list of panels. */
const MTH_SECTIONS = [
{id: "mth-shared", title: "Used by every figure", panels: [
  {h: "Bootstrap interval, Figures 3, 4 and S1", items: [
    ["Every 95% interval in Figures 3, 4 and S1: a percentile bootstrap of the mean over simulations, never over clades or events.",
     [mthT`\bar v = \frac{1}{S}\sum_{s=1}^{S} v_s`,
      mthT`\bar v^{*}_b = \frac{1}{S}\sum_{i=1}^{S} v_{s_{bi}},\qquad s_{bi}\sim\mathcal U\{1,\dots,S\},\quad b = 1,\dots,2000`,
      mthT`[\,\text{lo},\ \text{hi}\,] = \big[\,P_{2.5}(\bar v^{*}),\ P_{97.5}(\bar v^{*})\,\big]`],
     [[mthT`v_s`, "one simulation's value; non-finite values are dropped"],
      [mthT`S`, "simulations with a value; below 3 the mean is given and no interval"],
      [mthT`P_q`, "numpy's linear-interpolated percentile; generator seeded 20261116, so a reload repeats"]],
     "figures_gui/panels.py:25-26, 234-244"]]},
  {h: "Bootstrap interval, Figure 5", items: [
    ["The interval in the hover over each of Figure 5's A columns: its own bootstrap in the browser, 1,000 draws, read off the sorted draws by rank. Not drawn on the figure.",
     [mthT`\text{lo} = \bar v^{*}_{(r_1)},\quad r_1 = \max\!\big(1,\operatorname{round}(0.025B)\big),\qquad \text{hi} = \bar v^{*}_{(r_2)},\quad r_2 = \operatorname{round}(0.975B)`],
     [[mthT`\bar v^{*}_{(r)}`, "the r-th smallest of the B bootstrap means over simulations"],
      [mthT`B`, "1,000 draws, Park–Miller generator seeded 12345"]],
     "figures_gui/web/tab_panel3.js:115-132, 277, 317"]]},
  {h: "Burn-in day", items: [
    ["The burn-in day φ, read per simulation from its own trajectory. A simulation that never reaches it has none.",
     [mthT`\varphi = \min\Big\{\,t \;:\; \frac{D(t) + R(t)}{N} \ge 1\,\Big\},\qquad N = S(0) + I(0)`],
     [[mthT`D(t),\,R(t)`, "diagnosed and recovered counts in the simulation trajectory"],
      [mthT`N`, "the population, from the trajectory's first row"]],
     "figures_gui/prep/burnin.py:42-67"]]},
  {h: "Two rules that hold everywhere", items: [
    ["The simulation is the independent unit: a percentage is formed inside a simulation, then averaged over simulations. Scenarios are never pooled: they are a designed pair of contrasts (SOT against HIV_low moves shedding duration 63 → 109 d at prevalence 1%; HIV_low against HIV_high moves prevalence 1% → 12% at 109 d).",
     [mthT`\bar y = \frac{1}{S}\sum_{s} y_s,\qquad y_s = \frac{\text{part}_s}{\text{whole}_s}`],
     [[mthT`S`, "simulations of one scenario"]],
     "scripts/experiments_settings/impact_long_shedders_N3000_T1095.json (scenarios)"]]},
]},

{id: "mth-f1", title: "Figure 1 · Model", note: "Built by <code>figures_gui/fig1.py</code>, which runs <code>figures_gui/prep/fig1_save.build_figure_1</code> and records what each painter drew. The studio swaps two letters (<code>F1_LETTER</code> in tab_fig1.js): the divergence distribution, data key <code>E</code>, is shown as <b>D</b>; the global clock, data key <code>D</code>, as <b>E</b>.", panels: [
  {h: "A · Infectiousness by day", items: [
    ["The probability that a host infected on day 0 is infectious t days later: a chain of phases solved exactly by the matrix exponential.",
     [mthT`\mathbf p(t) = e^{At}\,\mathbf e_0,\qquad p_{\text{inf}}(t) = \sum_{i=5}^{18} p_i(t)`,
      mthT`A_{ii} = -\lambda_i,\qquad A_{i,i-1} = \lambda_{i-1},\qquad \lambda_i = \frac{n_m}{\tau_m}\ \ (i \text{ in phase } m)`],
     [[mthT`i`, "states 0–4 pre-detection (n=5), 5 pre-symptomatic (n=1), 6–18 infectious (n=13), 19 post-infectious (n=1), 20 recovered, absorbing"],
      [mthT`\tau_m`, "tau_1, tau_2, tau_4 and tau_3 (control) or tau_3_long (long shedders), read from each scenario's output"],
      [mthT`\mathbf e_0`, "starts in state 0; t = 0…799 d; drawn to 800 d with an edge case, else 400; the studio keeps every second day"]],
     "simplicity/intra_host_model.py:60-85, 157-188 · figures_gui/prep/fig1_preprocess_data.py:78-101 · figures_gui/fig1.py _panel_a"]]},
  {h: "B · Realized infection duration", items: [
    ["Each simulated individual's realized duration: standard individuals in control, long shedders elsewhere.",
     [mthT`d_i = t^{\text{not infected}}_i - t^{\text{infection}}_i`],
     [[mthT`i`, "rows of individuals_data.csv from 10 seeds, drawn with random.Random(42) from control's, the same seed names in every scenario; HIV is HIV_low alone"]],
     "figures_gui/prep/fig1_preprocess_data.py:31-46, 108-149"],
    ["The subsample: at most 150 individuals per scenario, so every violin rests on the same n.",
     [mthT`\mathcal D^{\text{shown}} = \begin{cases}\text{150 drawn without replacement} & |\mathcal D| > 150\\ \mathcal D & \text{otherwise}\end{cases}`],
     [[mthT`\mathcal D`, "one scenario's durations; pandas sample, random_state 42"]],
     "figures_gui/prep/fig1_plots.py:94-102"],
    ["The violin: seaborn's Gaussian KDE at Scott's factor, cut at the data's range, every violin scaled to the same maximum width; dashed lines are the quartiles.",
     [mthT`\hat f(y) = \frac{1}{n h}\sum_i \phi\!\Big(\frac{y - d_i}{h}\Big),\qquad h = \hat\sigma\, n^{-1/5},\qquad y \in [\min d,\ \max d]`,
      mthT`\text{half-width}(y) = w\,\frac{\hat f(y)}{\max \hat f},\qquad Q_{25},\,Q_{50},\,Q_{75}\ \text{of the shown } d_i`],
     [[mthT`\hat\sigma`, "sample SD of the shown durations"],
      [mthT`\phi`, "standard normal density"]],
     "figures_gui/prep/fig1_plots.py:105-114 · figures_gui/fig1.py _panel_b"]]},
  {h: "C · Published shedding durations", items: [
    ["One point per patient, the column Days (Exact); categories ordered by mean duration. Boxes are matplotlib's box statistics; outliers are not marked.",
     [mthT`\text{box} = [Q_1,\,Q_3],\quad \text{line} = Q_2,\quad \text{whiskers} = \big[\min\{x \ge Q_1 - 1.5\,\mathrm{IQR}\},\ \max\{x \le Q_3 + 1.5\,\mathrm{IQR}\}\big]`],
     [[mthT`\mathrm{IQR}`, "Q3 − Q1"],
      [mthT`x`, "one patient; coloured as the edge case beyond 300 d whatever the category; jitter is cosmetic"]],
     "figures_gui/prep/fig1_preprocess_data.py:168-184 · figures_gui/prep/fig1_plots.py:124-171 · figures_gui/fig1.py _panel_c"]]},
  {h: "D · Divergence distribution (data key E)", items: [
    ["Patient sequences of the global-clock cohort: divergence, and day since the first sample.",
     [mthT`y_i = \frac{\text{hamming}_i}{\text{overlap}_i},\qquad x_i = \text{days\_since\_ref}_i - \min_j \text{days\_since\_ref}_j`],
     [[mthT`i`, "rows of patient_global_clock.csv with overlap_coverage ≥ 0.95"]],
     "figures_gui/prep/fig1_preprocess_data.py:197-208 · figures_gui/prep/fig1_plots.py:484-488"],
    ["Model sequences resampled so their sampling days follow the patients', in 10 quantile bins of the patients' days.",
     [mthT`n_b = \#\{\,i : x^{\text{pat}}_i \in [e_b, e_{b+1})\,\},\qquad e_0,\dots,e_{10} = \text{deciles of } x^{\text{pat}}`,
      mthT`\text{draw } n_b \cdot \max\!\big(1,\lfloor N_{\text{mod}}/N_{\text{pat}}\rfloor\big) \text{ model sequences from bin } b,\ \text{with replacement}`],
     [[mthT`x^{\text{mod}}`, "Sequencing_time × 365.25 of control's standard individuals, sequencing_data_regression.csv, 10% of hosts per seed"],
      [mthT`e_0, e_{10}`, "opened to ±∞ for the model; generator seeded 4242"]],
     "figures_gui/prep/fig1_plots.py:435-457, 494-499 · figures_gui/prep/fig1_preprocess_data.py:217-266"],
    ["The density curves: Gaussian KDE at Scott's factor, the bandwidth floored at 1.2 grid steps, each series at its own step.",
     [mthT`\hat f(y) = \frac{1}{n h}\sum_i \phi\!\Big(\frac{y - y_i}{h}\Big),\qquad h = \hat\sigma \cdot \max\!\Big(n^{-1/5},\ \frac{1.2\,\Delta}{\sigma}\Big)`,
      mthT`\Delta_{\text{model}} = \min_{u_{j+1} \ne u_j}(u_{j+1} - u_j)\ \ \text{if} \ge \tfrac{\text{span}}{500},\qquad \Delta_{\text{pat}} = \frac{1}{\operatorname{median}(\text{overlap})}`],
     [[mthT`\sigma,\ \hat\sigma`, "population SD (in the floor) and scipy's sample SD"],
      [mthT`u_j`, "the model's distinct values, sorted"],
      [mthT`y`, "600 points on [0, largest value of either series]"]],
     "figures_gui/prep/fig1_plots.py:389-405, 418-432, 488, 500, 509-511"],
    ["The faint steps: a histogram per series on bins twice its own step, normalised to a density.",
     [mthT`\hat f_b = \frac{c_b}{n\,\omega},\qquad \omega = 2\Delta\ \ \big(\text{or } 2\,y_{\max}/120 \text{ with no step}\big)`],
     [[mthT`c_b`, "values in bin b"]],
     "figures_gui/prep/fig1_plots.py:408-415, 512-517"],
    ["The y axis is the density in plain units: the pipeline draws it ÷ 1,000 and the studio multiplies back. Rug ticks are up to 1,500 of each series' values at random (seed 42).",
     [mthT`\text{shown} = 10^{3}\cdot\big(10^{-3}\hat f\big) = \hat f\quad [\text{per subs/site}]`],
     [],
     "figures_gui/prep/fig1_plots.py:516-525 · figures_gui/web/tab_fig1.js f1PaintE · figures_gui/fig1.py _panel_e"]]},
  {h: "E · Global clock (data key D)", items: [
    ["Every line is a regression through the origin; its rate and R².",
     [mthT`b = \frac{\sum_i x_i y_i}{\sum_i x_i^2},\qquad R^2 = 1 - \frac{\sum_i (y_i - b\,x_i)^2}{\sum_i (y_i - \bar y)^2},\qquad \text{rate} = 365.25\,b\ \ [\text{s/s/y}]`],
     [[mthT`R^2`, "against the centred total sum of squares, although the fit has no intercept"]],
     "figures_gui/prep/fig1_plots.py:275-286"],
    ["Patients: the points and the dashed line, x and y as in D.",
     [mthT`y = b_{\text{pat}}\,x`],
     [],
     "figures_gui/prep/fig1_plots.py:314-326"],
    ["Model: the line is the mean of per-seed slopes, the band ±1 SD across seeds; R² is one fit through every seed pooled.",
     [mthT`\bar b = \frac{1}{S}\sum_{s} b_s,\qquad \mathrm{SD} = \sqrt{\frac{1}{S-1}\sum_s (b_s - \bar b)^2}`,
      mthT`\text{line } y = \bar b\,x,\qquad \text{band } y \in \big[(\bar b - \mathrm{SD})\,x,\ (\bar b + \mathrm{SD})\,x\big]`],
     [[mthT`b_s`, "through-origin slope of seed s"],
      [mthT`x,\ y`, "Sequencing_time × 365.25 and Distance_from_root of control's standard individuals, 10% of hosts per seed (sequencing.subsample, seed 20261020 + s)"],
      [mthT`\text{points}`, "the model's cloud, when ticked: 600 at random, seed 42"]],
     "figures_gui/prep/fig1_plots.py:289-290, 328-365 · figures_gui/prep/fig1_preprocess_data.py:217-266 · figures_gui/fig1.py _clock"]]},
  {h: "F · Intra-host clock", items: [
    ["Patients: patient_intrahost_clock.csv, no coverage cut.",
     [mthT`x_i = \text{time\_point}_i - \min_j \text{time\_point}_j,\qquad y_i = \text{hamming\_distance\_per\_overlap}_i`],
     [],
     "figures_gui/prep/fig1_save.py:100-104 · figures_gui/prep/fig1_plots.py:314-326"],
    ["Model: one point per lineage in every long shedder's intra-host record, every seed; fitted and drawn as in E.",
     [mthT`x = t^{\text{ih birth}}_{\ell} - t^{\text{infection}}_{h},\qquad y = \frac{\operatorname{hamming}\big(g_\ell,\ g_{\text{inherited}(h)}\big)}{L_{\text{ref}}}`],
     [[mthT`h,\ \ell`, "a long shedder and a lineage in its IH_lineages_trajectory"],
      [mthT`L_{\text{ref}}`, "reference genome length"]],
     "simplicity/tuning/evolutionary_rate.py:39-104 · figures_gui/prep/fig1_preprocess_data.py:269-314"]]},
]},

{id: "mth-f2", title: "Figure 2 · Clustering", note: "One seed per scenario with long shedders, drawn from the per-seed export (<code>/api/run</code>); the clade partition is recomputed in the browser. It agrees with the clade build when <code>runs/</code> was written by <code>figures_gui/export_clade_explorer.py</code> since v2.5.24, which counts only transmission events per lineage; an older export also counts the seeded infections, so near the m<sub>T</sub> floor its partition can differ.", panels: [
  {h: "Clade (rule B)", items: [
    ["A lineage founds a clade when it has diverged k mutations from its clade root and its subtree carries m_T of the simulation's transmission events. Seeded infections (parent = root) had no transmitter, so they count in neither. A refused lineage stays in its ancestor's clade.",
     [mthT`\ell \text{ founds a clade} \iff \big|M(\ell)\setminus M(r)\big| \ge k \ \ \wedge\ \ T(\ell) \ge m_T\,T_{\text{sim}}`,
      mthT`T(\ell) = \sum_{\ell' \in \text{subtree}(\ell)} \#\{\,h : \text{inherited\_lineage}(h) = \ell',\ \text{parent}(h) \ne \text{root}\,\},\qquad T_{\text{sim}} = \#\{\,h : \text{parent}(h) \ne \text{root}\,\}`],
     [[mthT`M(\ell)`, "the set of mutations in lineage ℓ's genome"],
      [mthT`r`, "the root of the clade ℓ's parent belongs to"],
      [mthT`k,\ m_T`, "shared-mutation threshold and minimum transmissions, a fraction; m_T = 0 drops the second term"],
      [mthT`\text{parent}(h) = \text{root}`, "a seeded infection: no transmitter, so not a transmission event (<code>transmitted_lineages</code>, since v2.5.24)"]],
     "simplicity/clustering.py:188-230, 269-288 · figures_gui/compute.py:139-145 · figures_gui/web/rundata.js:113-142"]]},
  {h: "Phylogeny", items: [
    ["Each lineage drawn as an elbow from its parent; tips stacked in depth-first order.",
     [mthT`x_\ell = \text{Time\_emergence}(\ell)\ \ \text{or}\ \ |M(\ell)|,\qquad y_{\text{tip}} = 0,1,2,\dots,\qquad y_{\text{node}} = \operatorname{mean}_{c \in \text{children}}\, y_c`],
     [[mthT`x`, "the Tree x axis: Time, or Divergence (mutations in the genome)"]],
     "figures_gui/web/paint.js:77-98 · figures_gui/web/rundata.js:80-110"],
    ["Colour by origin (Ratio rule): a mix of the standard and long colours by the long-born share of the clade's defining mutations.",
     [mthT`p_c = \frac{n^{\text{long}}_c}{n^{\text{long}}_c + n^{\text{std}}_c},\qquad \text{colour} = \operatorname{mix}(\text{std},\ \text{long};\ p_c)`],
     [[mthT`n_c`, "mutation counts on the edges from the clade root up to, not including, the parent clade's root"]],
     "figures_gui/web/rundata.js:170-195 · figures_gui/web/paint.js:58-67"]]},
  {h: "Clade frequency", items: [
    ["A clade's frequency is the sum of its lineages' frequencies; clades are stacked in tree order, capped at the whole population.",
     [mthT`F_c(t_j) = \sum_{\ell \in c} f_\ell(t_j),\qquad \text{band}_c(t_j) = \Big[\min\big(1,\textstyle\sum_{c' < c} F_{c'}\big),\ \min\big(1,\sum_{c' \le c} F_{c'}\big)\Big]`],
     [[mthT`f_\ell(t_j)`, "lineage frequency in the per-seed export: 120 time points about 9.2 d apart, stored in steps of 1/65535"],
      [mthT`c`, "the 70 clades with the highest peak get a band; the rest share one grey band"]],
     "figures_gui/web/rundata.js:217-231 · figures_gui/web/paint.js:21, 302-349"],
    ["The burn-in day: a dashed line at φ (see above), on the frequency plot always and on the tree only on the Time axis.",
     [mthT`x = \varphi`],
     [],
     "figures_gui/web/tab_clades4.js:132-136 · figures_gui/web/paint.js:372-382"]]},
]},

{id: "mth-f3", title: "Figure 3 · Transmission", note: "From the transmission table (<code>/api/panel?n=2</code>). With the burn-in box ticked, only days from ⌊φ⌋ count and simulations with no φ are left out.", panels: [
  {h: "Transmission events", items: [
    ["Each infection with a transmitter, counted by the host type of each side and the day it happened. Seeded infections (parent = root) are not events.",
     [mthT`e_{xy}(d) = \#\{\,h : \text{type}(\text{parent}(h)) = x,\ \text{type}(h) = y,\ \lfloor t^{\text{inf}}_h \rfloor = d\,\},\qquad E_{xy} = \sum_{d \ge d_0} e_{xy}(d)`],
     [[mthT`x,\,y`, "l long shedder, s standard"],
      [mthT`d_0`, "⌊φ⌋ with the burn-in box, else 0"]],
     "figures_gui/export_transmissions.py:79-91 · figures_gui/panels.py:436-450"]]},
  {h: "Violins · new substitutions on transmission", items: [
    ["The jump Δ of one event: genome positions where the lineage passed on differs from the one the transmitter was infected with.",
     [mthT`\Delta_e = \operatorname{hamming}\big(g_{\text{inherited}(x)},\ g_{\text{transmitted}(e)}\big)`],
     [[mthT`x`, "the transmitter, split by its individuals_data type, long_shedder or standard"]],
     "figures_gui/export_transmissions.py:105-128"],
    ["Events pooled over the scenario's simulations into one integer histogram per transmitter type.",
     [mthT`H(j) = \sum_{s}\ \#\{\,e \in s : \Delta_e = j,\ d_e \ge d_0\,\},\qquad n = \sum_j H(j)`],
     [],
     "figures_gui/panels.py:441-446, 488-496"],
    ["The subsample: 800 of the n events, drawn without replacement with a fixed seed; quartiles by the linear rule.",
     [mthT`Q_p = v_{(i)} + (z - i)\big(v_{(i+1)} - v_{(i)}\big),\qquad z = p\,(n' - 1),\quad i = \lfloor z \rfloor`],
     [[mthT`v_{(i)}`, "the drawn values sorted, from 0; n' = 800, or n with the subsample off (the server's quartiles)"],
      [mthT`\text{seed}`, "42 + 7919 × scenario index + violin index"]],
     "figures_gui/web/tab_panel2.js:95-116, 138 · figures_gui/panels.py:499-506"],
    ["The outline: a Gaussian KDE at Scott's bandwidth on the (sub)sample's histogram, cut at its range, scaled to the slot's width.",
     [mthT`\hat f(y) = \sum_j H(j)\,\exp\!\Big(-\tfrac12\Big(\frac{y - j}{h}\Big)^{2}\Big),\qquad h = \sigma\, n^{-1/5}\ \ (0.25 \text{ if } \sigma = 0)`,
      mthT`\sigma^2 = \frac{1}{n}\sum_j H(j)\,j^2 - \bar j^{\,2},\qquad y \in [j_{\min},\ j_{\max}],\qquad \text{half-width}(y) = w\,\frac{\hat f(y)}{\max \hat f}`],
     [[mthT`y`, "81 points; the y axis is 0–6 except the edge case, which runs to 1.08 × its largest Δ"]],
     "figures_gui/web/tab_panel2.js:146-163, 224-228"]]},
  {h: "Arrows · where transmissions go", items: [
    ["Each arrow's number: the mean over simulations of that cell's share of the simulation's events; the ± under it is half the bootstrap 95% interval.",
     [mthT`\text{share}_{xy} = \frac{1}{S}\sum_{s} \frac{E_{xy,s}}{E_{ll,s} + E_{ls,s} + E_{sl,s} + E_{ss,s}},\qquad \pm = \tfrac12\,(\text{hi} - \text{lo})`],
     [[mthT`[\text{lo},\text{hi}]`, "the bootstrap interval above"]],
     "figures_gui/panels.py:453-456, 472 · figures_gui/web/tab_panel2.js:32"],
    ["Arrow width grows with the square root of the share, on one scale for all four subpanels.",
     [mthT`w_{xy} = 0.8 + 14\,\sqrt{\text{share}_{xy}}\ \ \text{px}`],
     [[mthT`\text{share}`, "as a fraction, 0 to 1"]],
     "figures_gui/web/tab_panel2.js:14, 23"]]},
  {h: "Summary table", items: [
    ["Caused by long shedders; events per host of each type (intensity); their ratio; the mean jump.",
     [mthT`\text{caused} = \frac{1}{S}\sum_s \frac{E_{ll,s} + E_{ls,s}}{\sum_{xy} E_{xy,s}},\qquad I_{\text{long}} = \frac{1}{S}\sum_s \frac{E_{ll,s} + E_{ls,s}}{N_{\text{long},s}},\qquad I_{\text{std}} = \frac{1}{S}\sum_s \frac{E_{sl,s} + E_{ss,s}}{N_{\text{std},s}}`,
      mthT`\rho = \frac{I_{\text{long}}}{I_{\text{std}}},\qquad \bar\Delta_x = \frac{1}{S}\sum_s \operatorname{mean}_{e \in s,\ x}\,\Delta_e`],
     [[mthT`N_{\text{long},s}`, "long shedders infected in the simulation, the whole simulation even with the burn-in box"],
      [mthT`\rho`, "one mean over another, no interval"]],
     "figures_gui/panels.py:448-460 · figures_gui/export_transmissions.py:144-145 · figures_gui/web/tab_panel2.js:187-193"]]},
]},

{id: "mth-f4", title: "Figure 4 · Frequency AUC", note: "From the clade build at the rail's k and m<sub>T</sub> (<code>/api/panel?n=1h</code>). A simulation counts if its total frequency AUC is above 0 and it is in the transmission table.", panels: [
  {h: "Frequency AUC of a clade", items: [
    ["The area under a clade's daily frequency curve, from φ with the burn-in box, else over the whole record. Unit: days.",
     [mthT`f_\ell(t) = \frac{1}{H_t}\sum_{h} \frac{\mathbf 1[\ell \in L_h(t)]}{|L_h(t)|},\qquad F_c(t) = \sum_{\ell \in c} f_\ell(t)`,
      mthT`B_c = \int_{t \ge \varphi} F_c(t)\,dt \quad (\text{trapezoid rule on the daily record})`],
     [[mthT`H_t,\ L_h(t)`, "infected undiagnosed hosts on day t, and the distinct lineages host h carries"]],
     "simplicity/population.py:598-643 · figures_gui/compute.py:168-174, 347-350"]]},
  {h: "Origin groups (the bar segments)", items: [
    ["A clade's group is set by the long-born share of its own defining mutations: n is the clade's n_defining, not the clustering k.",
     [mthT`x_c = \frac{j_c}{n_c},\qquad \text{group}(c) = \begin{cases}\text{founding} & \text{no parent, or } n_c = 0\\ \text{standard only} & x_c = 0\\ \text{mostly standard} & 0 < x_c \le \tfrac12\\ \text{mostly long} & \tfrac12 < x_c < 1\\ \text{only long} & x_c = 1\end{cases}`],
     [[mthT`j_c`, "defining mutations whose host_type starts with “long”"],
      [mthT`n_c`, "n_defining, the clade's own count"]],
     "figures_gui/panels.py:377-379 · figures_gui/compute.py:337-341"]]},
  {h: "The three bars", items: [
    ["Hosts: long shedders among the simulation's infected individuals.",
     [mthT`\text{hosts} = \frac{1}{S}\sum_s \frac{N_{\text{long},s}}{N_{\text{infected},s}}`],
     [],
     "figures_gui/panels.py:359-360, 389-390"],
    ["Infections: every transmission event split by the clade of the lineage it passed on; seeded infections are not events and are left out (since v2.5.24); from φ when the build carries ntx_post and the burn-in box is ticked.",
     [mthT`\text{inf}_g = \frac{1}{S}\sum_s \frac{\sum_{c \in g} T_c}{\sum_{c} T_c},\qquad T_c = \#\{\,h : \text{inherited\_lineage}(h) \in c,\ \text{parent}(h) \ne \text{root}\,\}`],
     [[mthT`g`, "an origin group; any long = mostly standard + mostly long + only long, the headline's second number"]],
     "figures_gui/panels.py:366-384 · figures_gui/compute.py:154, 328-344"],
    ["Frequency AUC: each group's share of the simulation's total.",
     [mthT`\text{AUC}_g = \frac{1}{S}\sum_s \frac{\sum_{c \in g} B_c}{\sum_c B_c}`],
     [[mthT`\text{interval}`, "every segment carries the bootstrap above; a segment is drawn at least 5 px wide"]],
     "figures_gui/panels.py:370, 380 · figures_gui/web/tab_panel1.js:141, 183-186"]]},
  {h: "Table · percentage of frequency AUC by class", items: [
    ["The table under the chart uses the composition classes, which do use k.",
     [mthT`\text{share}_j = \frac{1}{S}\sum_s \frac{\sum_{c : \text{class}(c) = j} B_c}{\sum_c B_c},\qquad \text{class}(c) = \begin{cases}\text{wt} & \text{founding clade}\\ \text{other} & n_c \ne k\\ j_c & n_c = k\end{cases}`],
     [],
     "figures_gui/panels.py:257-262, 265-305"]]},
]},

{id: "mth-f5", title: "Figure 5 · Outbreaks", note: "Measured by <code>/api/explore_panel</code> under the exploration tabs' rules (<code>exploreQuery</code>), on every simulation's daily record, with m<sub>F</sub> = 0. The defaults below are those rules as shipped.", panels: [
  {h: "Candidates", items: [
    ["Every clade except the founding one with at least one defining mutation, among the 200 with the highest peak in the simulation's daily record.",
     [mthT`c \in \text{candidates} \iff \text{parent}(c) \ne \varnothing\ \wedge\ n_c > 0`],
     [],
     "figures_gui/explore.py:329-349 · figures_gui/compute.py:223, 280-284"]]},
  {h: "Smoothing", items: [
    ["A centred moving average, its window in days, edges shrinking rather than wrapping (kernel box, the default).",
     [mthT`\tilde y_i = \frac{1}{|J_i|}\sum_{j \in J_i} y_j,\qquad J_i = \big[\max(0,\,i - \lfloor w/2 \rfloor),\ \min(n-1,\,i + \lfloor w/2 \rfloor)\big],\qquad w = \operatorname{odd}\!\big(\operatorname{round}(W/\delta)\big)`],
     [[mthT`W`, "window in days, 15"],
      [mthT`\delta`, "mean step, 1 d"],
      [mthT`\text{other kernels}`, "Gaussian σ = w/2.3548 (full width at half maximum W), median; always on the frequency itself"]],
     "figures_gui/prep/clade_preprocess.py:48-63 · figures_gui/explore.py:123-146"]]},
  {h: "Episodes", items: [
    ["A peak counts when its prominence clears a bar set as a fraction of its own height (Prominence relative to: Each peak, the clade build's), or of the clade's highest frequency (Clade maximum, the rule before v2.5.25).",
     [mthT`\beta_i = \max(p\,\tilde y_i,\ p_{\text{abs}})\ \text{(each peak)},\qquad \beta_i = \max(p\,h,\ p_{\text{abs}}),\ h = \max_j \tilde y_j\ \text{(clade maximum)},\qquad i \text{ is a peak} \iff \operatorname{prom}(i) \ge \beta_i`,
      mthT`\operatorname{prom}(i) = \tilde y_i - \max\!\Big(\min_{[L_i,\,i]} \tilde y,\ \min_{[i,\,R_i]} \tilde y\Big)`],
     [[mthT`p`, "prominence, 0.10; p_abs is held at 0"],
      [mthT`L_i,\,R_i`, "the nearest higher point on each side, or the record's end (scipy find_peaks)"]],
     "figures_gui/explore.py:149-159"],
    ["A censored peak: the clade's maximum, when it rose to it by the bar and never falls back by the bar before the record ends.",
     [mthT`m = \min\{\,i : \tilde y_i \ge h\,(1 - 10^{-9})\,\},\qquad h - \min_{i \le m}\tilde y_i \ge \beta\ \ \wedge\ \ h - \min_{i \ge m}\tilde y_i < \beta`],
     [[mthT`m`, "may be the last day here; the clade build leaves that case out"]],
     "figures_gui/explore.py:182-185 · figures_gui/prep/clade_preprocess.py:121-124"],
    ["An episode starts at the last dip before its peak that lies at least the bar below it.",
     [mthT`a(i) = \max\{\,m \in \{0\} \cup \mathcal T : m < i,\ \tilde y_m \le \tilde y_i - \beta_i\,\}\ \ (0 \text{ if none}),\qquad \mathcal T = \text{peaks of } -\tilde y \text{ with prominence} \ge \theta h`],
     [[mthT`\theta`, "dip depth, 10⁻⁴"],
      [mthT`\tilde y_m \le \tilde y_i - \beta_i`, "the “low” dip rule, which exploreQuery always sends"]],
     "figures_gui/explore.py:188-197 · figures_gui/web/episodes.js:322-324"],
    ["An episode's amplitude and top; its rise is days from 10% to 90% of that amplitude.",
     [mthT`\text{top} = \tilde y_{b},\qquad \text{amp} = \tilde y_{b} - \tilde y_{a},\qquad \text{rise} = t_{i_{90}} - t_{i_{10}},\quad i_q = \min\{\,i \ge a : \tilde y_i \ge \tilde y_a + q\cdot\text{amp}\,\}`],
     [[mthT`a,\ b`, "start and peak index"]],
     "figures_gui/prep/clade_preprocess.py:136-147"]]},
  {h: "Outbreak, hidden time, cryptic outbreak", items: [
    ["The outbreak: the clade's first episode, by start day, peaking at or above F.",
     [mthT`e^{*} = \text{first } e \text{ with } \text{top}_e \ge F`],
     [[mthT`F`, "5%; the outbreak pick “largest” takes the tallest instead"]],
     "figures_gui/web/episodes.js:258-262 · figures_gui/explore.py:485-488"],
    ["Hidden time runs from first appearance to the outbreak's start, clamped at first appearance; the launch frequency is the curve's value at the start.",
     [mthT`t_{\text{start}} = \max\big(t_{a^*},\ t_{\text{em}}\big),\qquad \text{hidden} = t_{\text{start}} - t_{\text{em}},\qquad \text{launch} = \max\big(0,\ \text{top}^{*} - \text{amp}^{*}\big)`],
     [[mthT`t_{\text{em}}`, "first appearance (emergence in host), the root lineage's Time_emergence; hidden time always starts here (exploreQuery sends anchor=emergence)"]],
     "figures_gui/web/episodes.js:274-279 · figures_gui/explore.py:266-295, 478"],
    ["Cryptic, rule 1: hidden at least H days and, under the launch rule, never above the launch frequency while hidden.",
     [mthT`\text{cryptic} \iff \text{hidden} \ge H\ \wedge\ \neg\Big(\text{quiet} \wedge \max_{t_{\text{em}} \le t_i < t_{a^*}} \tilde y_i > \text{launch}\Big)`],
     [[mthT`H`, "60 d"],
      [mthT`\text{quiet}`, "the launch rule, on"]],
     "figures_gui/web/episodes.js:280-291 · figures_gui/explore.py:373-380"],
    ["Cryptic, rule 2: hidden means never above f. Hidden time ends the first day the smoothed curve is above f, and the outbreak is cryptic if that is at least H days.",
     [mthT`t_{\text{cross}} = \min\{\,t_i \ge t_{\text{em}} : \tilde y_i > f\,\},\qquad \text{hidden} = t_{\text{cross}} - t_{\text{em}},\qquad \text{cryptic} \iff \text{hidden} \ge H`],
     [[mthT`f`, "the ceiling, set on Cryptic outbreaks (brkLo)"],
      [mthT`t_{\text{cross}}`, "the outbreak's start when the curve never exceeds f"]],
     "figures_gui/web/episodes.js:284-287 · figures_gui/explore.py:363-370"],
    ["The burn-in rule (“start”, the default): measured on the whole record; a simulation with no φ is left out, and a clade only counts if its outbreak starts after φ.",
     [mthT`\text{counted} \iff \varphi \text{ exists}\ \wedge\ t_{\text{start}} > \varphi`],
     [[mthT`\text{trunc}`, "the alternative: measure from φ and leave out clades with t_em < φ and frequency above 10⁻⁴ at φ"]],
     "figures_gui/web/episodes.js:267-268, 278 · figures_gui/web/tab_panel3.js:89-99 · figures_gui/explore.py:383-392"]]},
  {h: "Colour", items: [
    ["Fixed: Figure 4's origin groups in its colour-blind colours (<code>P1_CVD</code>), from the clade's long-born share of its defining mutations. Every candidate has n_c > 0, so the founding group does not occur.",
     [mthT`x_c = \frac{j_c}{n_c},\qquad \text{class}(c) = \begin{cases}\text{standard only} & x_c = 0\\ \text{mostly standard} & 0 < x_c \le \tfrac12\\ \text{mostly long} & \tfrac12 < x_c < 1\\ \text{only long} & x_c = 1\end{cases}`],
     [[mthT`n_c`, "n_defining"], [mthT`j_c`, "defining mutations that arose in a long shedder"]],
     "figures_gui/web/tab_panel3.js:55-67"]]},
  {h: "Plots", items: [
    ["A: average cryptic outbreaks per simulation in every scenario, control included, stacked by origin group. Each group's segment is its own mean, so the segments add up to the column. No whisker and no total are drawn; the column's bootstrap interval is only in its hover.",
     [mthT`\text{bar} = \frac{1}{S}\sum_{s} \#\{\,c \in s : \text{cryptic}\,\} = \sum_g \text{seg}_g,\qquad \text{seg}_g = \frac{1}{S}\sum_{s} \#\{\,c \in s \cap g : \text{cryptic}\,\}`],
     [[mthT`S`, "simulations kept"], [mthT`g`, "an origin group, as under Colour"]],
     "figures_gui/web/tab_panel3.js:141-142, 274-323"],
    ["A's stars: each scenario against control, the Mann–Whitney tests below, Holm-adjusted over the comparisons with control only.",
     [mthT`p^{\text{Holm}}_{(i)} = \max_{l \le i} \min\big(1,\ (m - l + 1)\, p_{(l)}\big),\qquad {*}\ p < 0.05,\ {**}\ p < 0.01,\ {***}\ p < 0.001`],
     [[mthT`m`, "comparisons with control, 4"]],
     "figures_gui/web/tab_panel3.js:188-201"],
    ["B–F: one point per cryptic outbreak. x, intra-host hidden time: first appearance to first transmission, inside the host where the clade arose; y, between-host hidden time: first transmission to the first day the smoothed frequency rises above f (<code>t_cross</code>, with each clade from <code>/api/explore_panel</code>). y shared by all five panels, x scaled to each.",
     [mthT`x = t_{\text{tx}} - t_{\text{em}},\qquad y = t_{>f} - t_{\text{tx}},\qquad \text{rule 2: } x + y = \text{hidden}`],
     [[mthT`t_{\text{tx}}`, "first transmission, the earliest infection with a transmitter that inherited any of the clade's lineages"],
      [mthT`t_{>f}`, "the first day from first appearance on that the smoothed frequency is above f"]],
     "figures_gui/web/tab_panel3.js:104, 150-152, 333-351 · figures_gui/explore.py:282-292, 484"]]},
  {h: "Tests (Cryptic by scenario, the same counts)", items: [
    ["Do scenarios differ? Kruskal–Wallis across scenarios, then every pair by two-sided Mann–Whitney U, Holm-adjusted; one count per simulation.",
     [mthT`H_{\text{KW}} \text{ on } \{\text{counts}_s\}_{\text{scenarios}},\qquad p^{\text{Holm}}_{(r)} = \max_{r' \le r}\ \min\big(1,\ (m - r' + 1)\,p_{(r')}\big),\qquad \text{significant} \iff p^{\text{Holm}} < 0.05`],
     [[mthT`p_{(r)}`, "the r-th smallest of the m pairwise p"],
      [mthT`\text{ties}`, "all values equal: H = 0, p = 1; U = n_a n_b / 2, p = 1"],
      [mthT`\text{table}`, "mean ± SD (n − 1) per simulation; cryptic share pools Σ cryptic / Σ outbreaks"]],
     "figures_gui/server.py count_tests · figures_gui/web/tab_sweep.js:504-527, 542-557, 602, 606-607"]]},
]},

{id: "mth-s1", title: "Supplementary S1 · Clustering sensitivity", note: "One row per clade build on disk made by the current code. Rise and latency are the clade build's emergence-episode measures (<code>panels.panel3</code>), not Figure 5's outbreaks.", panels: [
  {h: "Frequency AUC in the majority-long classes", items: [
    ["Per simulation, the frequency AUC share of classes with more than half of k defining mutations long-born; clades whose n_defining is not k are not in it.",
     [mthT`y = \frac{1}{S}\sum_s \frac{\sum_{c : n_c = k,\ j_c > k/2} B_c}{\sum_c B_c}`],
     [],
     "figures_gui/dist.py:80, 115-120 · figures_gui/panels.py:660-668"]]},
  {h: "Mean emergence rise time", items: [
    ["The clade build's episodes (box smoothing W, prominence p, dips at any depth over 10⁻⁴h), its emergence episode the one holding the clade's maximum.",
     [mthT`e^{\text{emg}} = \arg\min_e \big(\mathbf 1[\,i_{\max} \notin [a_e, b_e]\,],\ -\text{amp}_e\big),\qquad T_c = t_{i_{90}} - t_{i_{10}}`,
      mthT`y = \frac{1}{S}\sum_s \operatorname{mean}_{c \in s} T_c`],
     [[mthT`i_{\max}`, "index of the clade's smoothed maximum"]],
     "figures_gui/compute.py:367-390 · figures_gui/prep/clade_preprocess.py:70-147 · figures_gui/panels.py:511-528, 572-606"]]},
  {h: "Mean latent interval", items: [
    ["Days from first appearance to the emergence episode's start; none for a truncated clade, or when negative.",
     [mthT`L_c = t_{a^{\text{emg}}} - t_{\text{em}}(c),\qquad y = \frac{1}{S}\sum_s \operatorname{mean}_{c \in s,\ L_c \ge 0} L_c`,
      mthT`\text{truncated} \iff \varphi \text{ exists}\ \wedge\ t_{\text{em}} < \varphi\ \wedge\ F_c(\varphi) > 10^{-4}`],
     [],
     "figures_gui/panels.py:589-603, 617-619 · figures_gui/compute.py:354"]]},
  {h: "Sensitivity lines and summary", items: [
    ["Lines: one axis moved (k, m_T, window or prominence), every other build parameter pinned at the operating point; each level carries the bootstrap above.",
     [mthT`\text{levels} = \{\,b : b_{a'} = \text{op}_{a'}\ \ \forall a' \ne a\,\}`],
     [[mthT`\text{op}`, "k = 7, m_T = 1%, window 15 d, prominence 0.10, rise band 10–90%, post burn-in"]],
     "figures_gui/dist.py:274-335"],
    ["Summary: the operating point's value, beneath it the range over every build on disk.",
     [mthT`\text{span} = \max_b y_b - \min_b y_b`],
     [],
     "figures_gui/panels.py:740-784"]]},
]},
];

/* the stores a number passes through, GLOSSARY.md's names */
const MTH_STORES = `<div class="def"><h3>Where the numbers come from</h3>
  <p><b>transmission table</b> · <code>figures_gui/export_transmissions.py</code> → <code>tx/</code>: Figure 3, the hosts bar of Figure 4.</p>
  <p><b>per-seed export</b> · <code>figures_gui/export_clade_explorer.py</code> → <code>runs/</code>: Figure 2's phylogeny and frequency plot.</p>
  <p><b>clade build</b> · <code>figures_gui/compute.py</code> → <code>cache/&lt;parameters&gt;__&lt;hash&gt;/clades_&lt;scenario&gt;.json</code>: Figures 4 and S1, named for k, m<sub>T</sub>, window, prominence, rise band, burn-in window and the code hash.</p>
  <p><b>daily record</b> · the simulation's own <code>lineage_frequency</code>, read on first request and memoised: Figure 5 and the exploration tabs.</p>
  <p><b>panel cache</b> · <code>figures_gui/panels.py</code> → <code>cache/_panels/</code>: each panel's answer, keyed by the request and the code.</p></div>`;

function mthEntry([what, tex, sym, src]) {
  const syms = sym.length ? `<div class="msym">${sym.map(([s, t]) =>
    `<span class="mtx">${mthEsc(s)}</span><span>${t}</span>`).join("")}</div>` : "";
  return `<div class="mq"><p>${what}</p>${tex.map(t => `<div class="mtex">${mthEsc(t)}</div>`).join("")}
    ${syms}<div class="msrc">${mthEsc(src)}</div></div>`;
}
const mthSection = s => `<div class="sec" id="${s.id}"><h2>${s.title}</h2>
  ${s.note ? `<p class="mnote">${s.note}</p>` : ""}
  ${s.panels.map(p => `<div class="def"><h3>${p.h}</h3>${p.items.map(mthEntry).join("")}</div>`).join("")}</div>`;

/* KaTeX over every formula and symbol; the source stays as text if it is absent */
function mthTypeset(root) {
  const run = () => (root.querySelectorAll(".mtex, .mtx") || []).forEach(el => {
    if (el.dataset && el.dataset.done) return;
    const src = el.textContent;
    try { window.katex.render(src, el, {displayMode: el.classList.contains("mtex"), throwOnError: false});
          if (el.dataset) el.dataset.done = "1"; }
    catch (e) { el.textContent = src; }
  });
  if (typeof window === "undefined") return;
  if (window.katex) run(); else window.addEventListener("load", () => { if (window.katex) run(); }, {once: true});
}

function rail_method(rail) {
  const toc = document.createElement("div");
  toc.className = "card";
  toc.innerHTML = `<h2>Contents</h2><div class="toc">
    <a href="#method" data-to="mth-state">What is precomputed</a>
    ${MTH_SECTIONS.map(s => `<a href="#method" data-to="${s.id}">${s.title}</a>`).join("")}</div>`;
  /* the hash names the tab, so a link scrolls instead of navigating */
  toc.addEventListener("click", ev => {
    const a = ev.target && ev.target.closest ? ev.target.closest("a[data-to]") : null;
    if (!a) return;
    ev.preventDefault();
    const el = document.getElementById(a.dataset.to);
    if (el && el.scrollIntoView) el.scrollIntoView({behavior: "smooth", block: "start"});
  });
  rail.append(toc);
}

async function render_method(main) {
  main.innerHTML = `<div class="loading">loading…</div>`;
  let ix = null, note = "";
  try { ix = await apiGetSlot("/api/index"); }
  catch (e) { note = e.message; }

  let state;
  if (ix) {
    const b3 = ix.builds || [];
    const tx = ix.transmissions;
    const runs = b3.length ? Object.values(b3[0].counts || {}).reduce((a, c) => a + c, 0) : 0;
    const txruns = tx ? (tx.scenarios || []).length : 0;
    state = `<div class="stats">
      <div class="stat"><div class="k">clade builds cached</div>
        <div class="v">${b3.length}</div></div>
      <div class="stat"><div class="k">simulations per build</div>
        <div class="v">${runs} <small>of 250</small></div></div>
      <div class="stat"><div class="k">transmission table</div>
        <div class="v">${tx ? txruns + " " : "—"}<small>${tx ? "scenarios" : "missing"}</small></div></div>
      <div class="stat"><div class="k">compute code</div>
        <div class="v" style="font-size:.78rem">${ix.code}</div></div>
    </div>`;
  } else {
    state = `<div class="err">could not read the build state: ${note}</div>`;
  }

  main.innerHTML = `
    <div class="sec" id="mth-state"><h2>What is precomputed</h2>
      <p class="lede" style="margin-bottom:10px">Live, from this server's own store.</p>
      ${state}${MTH_STORES}</div>
    ${MTH_SECTIONS.map(mthSection).join("")}`;
  mthTypeset(main);
}

registerTab({
  id: "method", group: "tools", title: "Method", viewId: "viewMethod",
  lede: "The maths behind every plotted quantity, figure by figure: what is drawn, its formula, its symbols, and the file and line that computes it, so each can be checked against the code rather than trusted.",
  rail: rail_method,
  render: render_method,
});
