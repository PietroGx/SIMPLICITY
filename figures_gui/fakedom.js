/* A DOM small enough to run the builder's UI layer headlessly: innerHTML is
   stored, not parsed, so querySelector answers by selector shape. */
const CTX = {cur:null};
function mkCtx(){
  const noop=()=>{};
  const g={ops:0,texts:[],
    beginPath:()=>{g.ops++},moveTo:()=>{g.ops++},lineTo:()=>{g.ops++},arc:()=>{g.ops++},
    stroke:()=>{g.ops++},fill:()=>{g.ops++},fillRect:()=>{g.ops++},strokeRect:()=>{g.ops++},
    clearRect:noop,save:noop,restore:noop,translate:noop,rotate:noop,setTransform:noop,
    setLineDash:noop,drawImage:noop,closePath:noop,clip:noop,
    bezierCurveTo:()=>{g.ops++},quadraticCurveTo:()=>{g.ops++},rect:()=>{g.ops++},ellipse:()=>{g.ops++},
    createLinearGradient:()=>({addColorStop:noop}),
    createImageData:(w,h)=>({data:new Uint8ClampedArray(4*Math.max(1,w|0)*Math.max(1,h|0))}),
    getImageData:(x,y,w,h)=>({data:new Uint8ClampedArray(4*Math.max(1,w|0)*Math.max(1,h|0))}),
    putImageData:noop,scale:noop,transform:noop,
    fillText:(t)=>{g.ops++;g.texts.push(String(t))},
    measureText:t=>({width:String(t).length*5.6})};
  return g;
}
class El{
  constructor(id){ this.id=id||""; this.dataset={}; this.style={}; this._html="";
    this.textContent=""; this.hidden=false; this.value=""; this.placeholder="";
    this.children=[]; this.handlers={}; this.width=0; this.height=0;
    this.clientWidth=620; this.parentElement={clientWidth:620}; this.className="";
    this.offsetHeight=0; this._top=180; }
  getBoundingClientRect(){ return {top:this._top, left:0, width:this.clientWidth,
    height:this.offsetHeight, bottom:this._top+this.offsetHeight, right:this.clientWidth}; }
  /* a real DOM replaces the element's children here, and a rebuilt rail
     stacked a second copy of itself while this did not */
  set innerHTML(v){ this._html=String(v); this.children=[]; }
  get innerHTML(){ return this._html; }
  getContext(){ return CTX.cur || (CTX.cur=mkCtx()); }
  addEventListener(t,f){ (this.handlers[t]=this.handlers[t]||[]).push(f); }
  appendChild(c){ this.children.push(c); return c; }
  setAttribute(k,v){ this[k]=v; }
  getAttribute(k){ return this[k]; }
  querySelector(sel){
    if (sel==="canvas"){ const c=new El(); c.parentElement={clientWidth:620}; return c; }
    const e=new El(); e.addEventListener=()=>{}; return e;
  }
  querySelectorAll(sel){
    /* a grid built by innerHTML: hand back one canvas per <canvas> in it, so
       the paint loop runs instead of silently drawing nothing */
    if (/canvas\s*$/.test(String(sel||""))){
      const n=(this._html.match(/<canvas/g)||[]).length;
      return Array.from({length:n},()=>{ const c=new El(); c.parentElement={clientWidth:240}; return c; });
    }
    return [];
  }
  closest(){ return null; }
  select(){}
  focus(){}
  /* a removed element really leaves the registry, so a later getElementById
     returns null the way a browser's does. Without this the harness hands back
     a fresh stand-in for anything asked of it and a wire-after-remove bug only
     shows up in a real browser -- which is how one shipped. */
  remove(){ this.removed=true; if (this.id){ REG.delete(this.id); GONE.add(this.id); } }
  get classList(){ return {add:()=>{},remove:()=>{},toggle:()=>{}}; }
}
global.Path2D=class Path2D{ constructor(){ this.n=0 } moveTo(){this.n++} lineTo(){this.n++}
  rect(){this.n++} arc(){this.n++} ellipse(){this.n++} closePath(){} bezierCurveTo(){this.n++}
  quadraticCurveTo(){this.n++} };
const REG=new Map(), GONE=new Set();
function $el(id){
  /* the checks start a fresh page by emptying REG, which also retires the
     tombstones: an element removed by the previous boot exists again */
  if (!REG.size) GONE.clear();
  if (GONE.has(id)) return null;
  if(!REG.has(id)) REG.set(id,new El(id));
  return REG.get(id);
}
global.document={
  getElementById:id=>$el(id),
  createElement:()=>new El(),
  querySelector:sel=>$el("sel:"+sel),
  querySelectorAll:()=>[],
  documentElement:{scrollHeight:1000, style:{}}, body:new El("body"),
};
const VARS={"--long":"#e41a1c","--std":"#4daf4a","--mixed":"#377eb8","--founder":"#8172b2",
  "--sc0":"#9aa0a6","--sc1":"#6fa96f","--sc2":"#4c86c6","--sc3":"#d1495b","--sc4":"#9b6bb3",
  "--ink":"#1c1b19","--ink-2":"#55524d","--ink-3":"#8b8780","--line":"#e3e1dd","--line-2":"#efedea",
  "--panel":"#fff","--accent":"#4a6fa5","--f-num":"monospace","--f-ui":"sans-serif"};
global.getComputedStyle=()=>({getPropertyValue:n=>VARS[n]||"#888",gridTemplateColumns:"1fr 1fr"});
const STORE=new Map();
global.localStorage={getItem:k=>STORE.has(k)?STORE.get(k):null,setItem:(k,v)=>STORE.set(k,String(v)),
  removeItem:k=>STORE.delete(k)};
global.window={devicePixelRatio:1,innerHeight:900,innerWidth:1440,scrollY:0,
  matchMedia:()=>({matches:false,addEventListener(){}}),claude:null};
global.innerHeight=900;
global.matchMedia=global.window.matchMedia;
global.location={hash:"", pathname:"/", search:""};
global.history={replaceState:(a,b,url)=>{ const i=String(url).indexOf("#"); global.location.hash=i<0?"":String(url).slice(i); }};
global.addEventListener=()=>{};
global.navigator={clipboard:{writeText:async()=>{}}};
const timers=new Map(); let tid=0;
global.setTimeout=(f,ms)=>{ const id=++tid; timers.set(id,f); setImmediate(()=>{ if(timers.has(id)){ timers.delete(id); try{f()}catch(e){} } }); return id; };
global.clearTimeout=id=>{ timers.delete(id); };
/* serve the two real export directories from disk; anything else is still "no network",
   so a stray outbound fetch fails loudly instead of silently resolving */
global.fetch=async(u)=>{
  const fs=require("fs");
  const p=String(u).replace(/^[a-z]+:\/\/[^/]+\//i,"").replace(/^\.?\//,"").split(/[?#]/)[0];
  if(!/^(fig2b|data)\//.test(p)) throw new Error("no network in the harness: "+p);
  if(!fs.existsSync(p)) return {ok:false,status:404,
    async json(){throw new Error("404 "+p)}, async text(){return ""}};
  const txt=fs.readFileSync(p,"utf8");
  return {ok:true,status:200,async json(){return JSON.parse(txt)},async text(){return txt}};
};
function flushTimers(){ while(timers.length){ const f=timers.shift(); try{f()}catch(e){} } }
function loadApp(file, names0){
  const fs=require("fs");
  const html=fs.readFileSync(file,"utf8");
  let src=[...html.matchAll(/<script>([\s\S]*?)<\/script>/g)][0][1];
  const i=src.indexOf("(async function init(){");
  src=src.slice(0,i);
  const names=["S","PLOTS","PLOT_ORDER","METRICS","METRIC_ORDER","GLYPH","PRESETS","GROUPS","SCEN_ORDER",
    "TX_PCT","TX_PCT_LEGACY","txIdx","decodeTable","vals","computeVals","readCl","applyCl","clChips","encodeState",
    "parseState","applyState","stored","persist","useCtx","withCtx","ctxLabel","copyCl","copyAn",
    "savePanel","deletePanel","editPanel","toggleFig","moveFig","drawFig","drawPreview","selectPlot",
    "syncBuild","renderLib","buildPalette","buildScenList","buildNotes","wire","setView","figRecs",
    "ensureFor","updateCount","applyPreset","loadLib","saveLib","validRec","downloadPng","snapK",
    "defaultTitle","sameCl","amax","quantile","finite","okSeeds","activeScens","seedsOf","VCACHE","concatTables",
    "aucPerSeed","rhoPerSeed","ciOverSeeds","WIN_MIN","RHO_MIN","spearman","STAT","st","perSeed","auc","tieFrac","LADDER_PLOTS","thrLadder",
    "K_SET","ladderKs","NEEDS","metricOk","whyNot","okMetrics","paletteHtml","syncPalette","noteKey","noteLabel",
    "noteSync","noteWrite","notesText","fmlCard","metricHint",
    "scanM","fitM","fitRec","SUBS","subsOf","subsTot","fmtAx",
    "seedLabels","SEEDAGG","AGG_ORDER","agg","axLab"];
  const api=new Function(src+"\nreturn {"+(names0||names).join(",")+"};")();
  return {api, CTX, mkCtx, REG, STORE, flushTimers, El};
}
/* a table may be exported as several row ranges; read it the way the page does */
function readTable(api, dir, f){
  const fs=require("fs");
  const parts=Array.isArray(f.file)?f.file:[f.file];
  return api.concatTables(parts.map(p=>
    api.decodeTable(JSON.parse(fs.readFileSync(`${dir}/${p}`,"utf8")))));
}
module.exports={loadApp, CTX, mkCtx, REG, GONE, STORE, flushTimers, El, $el, readTable};
