/* Export: one figure tab as an image or PDF at a printed width and resolution.

   The figure is redrawn, not screenshotted: its column is laid out at a width
   that keeps the studio's 10 px labels legible in print, every canvas is
   redrawn at the matching pixel density (EXPORT_DPR in draw.js), and the
   column is then composed into one image with html2canvas. Rule lines,
   captions, notes and tables are the studio's own reading aids and are left
   out. The plots are raster at the chosen dpi, also inside the PDF. */

const EX_FORMATS = [["png", "PNG"], ["jpeg", "JPEG"], ["pdf", "PDF"]];
const EX_WIDTHS = [[89, "89 mm · 1 column"], [120, "120 mm · 1.5 columns"], [183, "183 mm · 2 columns"]];
const EX_DPIS = [[300, "300 dpi"], [600, "600 dpi"], [1200, "1200 dpi"]];
const EX_LIBS = {
  html2canvas: "https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js",
  jspdf: "https://cdnjs.cloudflare.com/ajax/libs/jspdf/2.5.1/jspdf.umd.min.js",
};
/* one CSS pixel prints at this many mm, so a 10 px label is about 5.7 pt */
const EX_MM_PER_PX = 0.2;
const EX_SKIP = ".cap,.note,.sw-rule,table.t,.loading,.err,.hinttog";

const exFigures = () => navOrder().filter(t => t.group === "figures" || t.group === "supp");
const exLayoutPx = mm => Math.round(Math.max(440, Math.min(1400, mm / EX_MM_PER_PX)));
const exPt = px => px * (+settings.exWidth / exLayoutPx(+settings.exWidth)) / 0.3528;
const exSlug = t => `${t.num != null ? "fig" + t.num : t.id}_${t.title}`.toLowerCase().replace(/[^a-z0-9]+/g, "_");

function rail_export(rail) {
  rail.append(card("Figure", [
    pickControl({setting: "exFig", label: "Figure",
      options: exFigures().map(t => [t.id, tabHeading(t), t.lede ? String(t.lede).replace(/<[^>]+>/g, "") : ""])}),
  ]));
  rail.append(card("File", [
    segmentControl({setting: "exFmt", label: "Format", options: EX_FORMATS,
      hint: "PNG and JPEG carry the resolution in the file; the PDF is one page the size of the figure."}),
    segmentControl({setting: "exWidth", label: "Printed width", options: EX_WIDTHS,
      hint: "Nature's column widths. The figure is laid out for this width, so labels print at the same size whichever you pick."}),
    segmentControl({setting: "exDpi", label: "Resolution", options: EX_DPIS,
      hint: "Line art wants 600 dpi or more; 300 dpi is enough for a draft."}),
  ]));
}

function exLoad(name) {
  if (name === "html2canvas" && window.html2canvas) return Promise.resolve();
  if (name === "jspdf" && window.jspdf) return Promise.resolve();
  return new Promise((ok, fail) => {
    const s = document.createElement("script");
    s.src = EX_LIBS[name]; s.onload = ok;
    s.onerror = () => fail(new Error(`could not load ${name} from cdnjs`));
    document.head.append(s);
  });
}

/* PNG's pHYs chunk, so an image viewer or a journal's checker reads the dpi */
const EX_CRC = (() => { const t = new Uint32Array(256);
  for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; t[n] = c >>> 0; }
  return t; })();
function exCrc(bytes) { let c = 0xffffffff; for (const b of bytes) c = EX_CRC[(c ^ b) & 0xff] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; }
async function exPngDpi(blob, dpi) {
  const src = new Uint8Array(await blob.arrayBuffer()), ppm = Math.round(dpi / 0.0254);
  const chunk = new Uint8Array(21), v = new DataView(chunk.buffer);
  v.setUint32(0, 9); chunk.set([0x70, 0x48, 0x59, 0x73], 4);
  v.setUint32(8, ppm); v.setUint32(12, ppm); chunk[16] = 1;
  v.setUint32(17, exCrc(chunk.subarray(4, 17)));
  const at = 33;                                   // 8-byte signature + the 25-byte IHDR chunk
  return new Blob([src.subarray(0, at), chunk, src.subarray(at)], {type: "image/png"});
}

function exSave(blob, name) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 5000);
}

/* Redraw one figure tab at print density, compose it, and hand back a canvas. */
async function exCapture(tab, mm, dpi) {
  const back = activeTabId, layout = exLayoutPx(mm), scale = (mm / 25.4 * dpi) / layout;
  const main = document.getElementById(`${tab.id}_main`), col = main.closest(".maincol");
  const was = col.style.cssText;
  /* width and density first: showTab may start its own render of a stale tab */
  col.style.cssText = `${was};width:${layout}px;max-width:none;flex:none`;
  EXPORT_DPR = scale;
  showTab(tab.id);
  try {
    await renderTab(tab);
    await exLoad("html2canvas");
    return await window.html2canvas(main, {scale, backgroundColor: css("--panel") || "#ffffff", logging: false,
      ignoreElements: el => el.matches && (el.matches(EX_SKIP) || (el.matches(".plot") && !el.querySelector("canvas")))});
  } finally {
    EXPORT_DPR = 0;
    col.style.cssText = was;
    staleTabs.add(tab.id);
    showTab(back);
  }
}

async function exRun(ids, say) {
  const mm = +settings.exWidth, dpi = +settings.exDpi, fmt = settings.exFmt;
  for (const id of ids) {
    const t = tabById(id);
    say(`${tabHeading(t)}: drawing at ${dpi} dpi…`);
    const cv = await exCapture(t, mm, dpi);
    const hmm = mm * cv.height / cv.width, name = `SIMPLICITY_${exSlug(t)}_${mm}mm_${dpi}dpi.${fmt === "jpeg" ? "jpg" : fmt}`;
    if (fmt === "pdf") {
      await exLoad("jspdf");
      const doc = new window.jspdf.jsPDF({orientation: mm > hmm ? "landscape" : "portrait", unit: "mm", format: [mm, hmm]});
      doc.addImage(cv.toDataURL("image/png"), "PNG", 0, 0, mm, hmm, undefined, "FAST");
      exSave(doc.output("blob"), name);
    } else {
      const blob = await new Promise(ok => cv.toBlob(ok, `image/${fmt}`, 0.95));
      exSave(fmt === "png" ? await exPngDpi(blob, dpi) : blob, name);
    }
    say(`${tabHeading(t)}: saved ${name}, ${cv.width} × ${cv.height} px, ${fmt.toUpperCase()}.`);
  }
}

async function render_export(main) {
  const figs = exFigures();
  if (!figs.some(t => t.id === settings.exFig)) settings.exFig = figs[0].id;
  const mm = +settings.exWidth, dpi = +settings.exDpi, px = Math.round(mm / 25.4 * dpi);
  const t = tabById(settings.exFig);
  const row = (k, v) => `<tr><td>${k}</td><td class="n">${v}</td></tr>`;
  main.innerHTML = `<div class="plot"><div class="ph"><h3>${tabHeading(t)}</h3>
      <span class="sub">${settings.exFmt.toUpperCase()} · ${mm} mm · ${dpi} dpi</span></div>
    <table class="t"><tbody>
      ${row("width", `${mm} mm = ${px.toLocaleString()} px`)}
      ${row("layout", `the figure is laid out ${exLayoutPx(mm)} px wide, so 10 px labels print at ${fmt(exPt(10), 1)} pt`)}
      ${row("height", "whatever the layout needs at that width")}
      ${row("file", `SIMPLICITY_${exSlug(t)}_${mm}mm_${dpi}dpi.${settings.exFmt === "jpeg" ? "jpg" : settings.exFmt}`)}
    </tbody></table>
    <div class="grid2">
      <button class="wide prim" id="ex_one">Export ${tabHeading(t)}</button>
      <button class="wide" id="ex_all">Export all ${figs.length} figures</button></div>
    <pre class="log" id="ex_log">Data settings come from the exploration tabs and the look from each
figure's own tab: the export draws exactly what that tab shows, without its rule
line, caption, notes or tables.</pre></div>`;
  const log = document.getElementById("ex_log"), say = s => { log.textContent += "\n" + s; };
  const go = ids => async () => {
    main.querySelectorAll("button").forEach(b => b.disabled = true);
    try { await exRun(ids, say); } catch (e) { say("failed: " + e.message); }
    main.querySelectorAll("button").forEach(b => b.disabled = false);
  };
  document.getElementById("ex_one").onclick = go([t.id]);
  document.getElementById("ex_all").onclick = go(figs.map(f => f.id));
}

registerTab({
  id: "export", group: "tools", title: "Export", viewId: "viewExport",
  lede: "Save a figure as PNG, JPEG or PDF at a journal column width and resolution.",
  rail: rail_export,
  render: render_export,
});
