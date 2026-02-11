// Particle-only p5.js (mobile-friendly). PNG frames are sampled offscreen only.

// ---- Config ----
// Runtime error overlay (helps when the dev server keeps reloading / canvas stays black).
(() => {
  if (window.__fatalOverlayInstalled) return;
  window.__fatalOverlayInstalled = true;
  const el = document.createElement("pre");
  el.id = "fatal-overlay";
  el.style.cssText = [
    "position:fixed",
    "left:12px",
    "top:12px",
    "right:12px",
    "max-height:45vh",
    "overflow:auto",
    "z-index:9999",
    "padding:10px 12px",
    "border-radius:10px",
    "background:rgba(0,0,0,0.78)",
    "color:#ffb3b3",
    "font:12px/1.35 ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, 'Liberation Mono', 'Courier New', monospace",
    "white-space:pre-wrap",
    "display:none",
  ].join(";");
  const attach = () => {
    try {
      if (el.parentNode) return;
      if (document.body) document.body.appendChild(el);
    } catch (_) {}
  };
  // Scripts in <head> (like CDN p5.js) can delay DOMContentLoaded if they hang.
  // We attach immediately when possible so errors are visible even if DOMContentLoaded is late.
  attach();
  document.addEventListener("DOMContentLoaded", attach);
  const show = (msg) => {
    try {
      attach();
      window.__fatalError = msg;
      el.textContent = msg;
      el.style.display = "block";
    } catch (_) {}
  };
  // If CDN scripts fail (common on phones/offline), p5 never loads and the screen stays black.
  // Detect that early and show a clear message.
  window.setTimeout(() => {
    try {
      if (typeof window.p5 === "undefined") {
        show("[ERROR]\np5.js failed to load.\nCheck internet/CDN access, or host p5 locally.");
      }
    } catch (_) {}
  }, 1500);
  window.addEventListener("error", (e) => {
    // Ignore errors from browser extensions; they can appear on GitHub Pages / Live Server
    // and shouldn't stop the sketch.
    try {
      const file = String(e?.filename || "");
      if (file.startsWith("chrome-extension://")) return;
    } catch (_) {}
    const msg = String(e?.message || e || "Unknown error");
    const src = e?.filename ? `\n@ ${e.filename}:${e.lineno || 0}:${e.colno || 0}` : "";
    show(`[ERROR]\n${msg}${src}`);
  });
  window.addEventListener("unhandledrejection", (e) => {
    const msg = String(e?.reason?.stack || e?.reason || "Unhandled rejection");
    // Same: ignore extension-origin stacks so they don't freeze the sketch.
    if (msg.includes("chrome-extension://")) return;
    show(`[REJECTION]\n${msg}`);
  });
})();
const ACTS = [1, 2, 3, 4];
const SRC_COUNT = { 1: 192, 2: 172, 3: 192, 4: 136 };
// Source frames per "sound-second" of internal time `t`.
// With the dt-based time step, 1.0 sound-second ~= 1 real second at full input speed.
// Chosen so: totalActsDuration + totalTransitionsDuration ~= 60s for a full 1→2→3→4→1 loop.
const FPS_EFFECTIVE = { 1: 14, 2: 14, 3: 14, 4: 14 };
// Frame stepping: set to 1 to use every source frame (full fidelity).
const FRAME_STEP = { 1: 1, 2: 1, 3: 1, 4: 1 };

const TARGET_W = 160, TARGET_H = 284;
const CELL_SIZE = 18; // square grid size (px) (tune: larger -> fewer cells -> denser silhouettes)
// Fixed simulation grid (keeps particle redistribution stable across different viewport sizes/hosts).
// This matches the common "good" resolution from local runs.
const GRID_COLS = 137;
const GRID_ROWS = 74;
const N = 4200, BLACK_PCT = 0.1;
const TRANSITION_DURATION = 2.6, DANCE_PORTION = 0.78;
const MIC_THRESHOLD = 0.03;
// Auto mode speed in sound-seconds per real second (applied with dt).
const AUTO_SPEED = 1.0;
// Deterministic seed so local vs GitHub Pages initialize identically.
const SIM_SEED = 1337;
// Backing-store density cap for performance (still uses devicePixelRatio up to this cap).
// Keeping this >1 avoids the "upscaled/soft" look on high-DPR hosts (e.g., GitHub Pages at 125%/150% zoom).
const PIXEL_DENSITY_CAP = 1.5;
function applyPixelDensity() {
  try {
    const dpr = max(1, (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1));
    pixelDensity(min(PIXEL_DENSITY_CAP, dpr));
  } catch (_) {}
}
function applySimSeed() {
  try { randomSeed(SIM_SEED); } catch (_) {}
  try { noiseSeed(SIM_SEED); } catch (_) {}
}
// Page zoom compensation:
// GitHub Pages can end up with a different per-site zoom than localhost (even on the same device/browser),
// which changes `windowWidth/windowHeight` and therefore the grid resolution + particle look.
// We measure zoom and compensate the grid cell size so the simulation stays consistent.
let pageZoom = 1;
function _measurePageZoom() {
  try {
    const d = document.createElement("div");
    d.style.cssText = "position:absolute;left:-1000px;top:-1000px;width:1in;height:1in;";
    document.body.appendChild(d);
    const px = d.getBoundingClientRect().width || 96;
    d.remove();
    // 1in should be 96 CSS px at 100% zoom.
    const z = px / 96;
    if (isFinite(z) && z > 0.2 && z < 6) return z;
  } catch (_) {}
  // Fallback (less accurate for zoom, but better than nothing)
  try {
    const dpr = (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1);
    return isFinite(dpr) && dpr > 0 ? dpr : 1;
  } catch (_) {}
  return 1;
}
function updatePageZoom() {
  // Clamp to avoid wild values on some mobile browsers.
  pageZoom = constrain(_measurePageZoom(), 0.5, 3);
}

// Scene visibility (fade in/out based on mic level; independent of sound-time `t` so it can fade out on silence)
const VIS_IN = 0.16;
const VIS_OUT = 0.08;
// Use a separate (lower) threshold than MIC_THRESHOLD so visuals still react on quiet mics.
const VIS_LEVEL_ON = 0.006;
const VIS_LEVEL_FULL = 0.030;
const INTERP_SHARPNESS = 1.6; // >1 reduces "double exposure" trails during motion

// Density contrast (tune): higher GAMMA -> stronger emphasis on bright areas (more contrast)
const DENSITY_THR = 0.07;
const DENSITY_GAMMA = 1.25;
const DENSITY_GREY_LIFT = 0.10;
const DENSITY_SMOOTH = 0.28;

// Particle look (tune)
const PARTICLE_SIZE = 3;
const LIGHT_ALPHA = 180;


// "Void travel" cleanup: when particles must cross empty (alpha=0) areas, guide them along 3 clean lanes
// (1 straight + 2 arced) between their current position and their destination.
const VOID_MIN_DIST_CELLS = 4.0;  // only engage lanes if far enough (in cell units)
const VOID_SAMPLES = 3;           // line samples to detect crossing empty space
const VOID_EMPTY_NEED = 2;        // how many samples must be empty to engage lanes
const LANE_LOOKAHEAD = 0.08;      // progress lookahead along lane (0..1)
const LANE_DT_REF = 0.03;         // reference sound-time step (used to normalize lane speed)

// ---- State ----
let mic, amp, audioStarted = false;
let micRunning = false;
let t = 0, _prevT = 0, tAdvanced = false, tDelta = 0;
let debugOn = false, imagesDrawnToCanvas = false;
let showGrid = false;
let autoRun = false;
let sceneVis = 0; // 0..1
let sceneA = 0;   // eased alpha used by renderers
let catchUp = 0;
let debugMeanCellDist = 0;
let debugTransport = 0;

let COLS = 1, ROWS = 1, CELLS = 1;
let gridX0 = 0, gridY0 = 0;
let cellW = 1, cellH = 1, invCellW = 1, invCellH = 1;

function resizeGrid() {
  // Letterbox the simulation to a fixed grid resolution so the computed density field
  // (and therefore particle look) stays consistent across different viewport sizes/aspects.
  COLS = GRID_COLS;
  ROWS = GRID_ROWS;
  CELLS = COLS * ROWS;

  const aspect = COLS / max(1, ROWS);
  let gw = min(width, height * aspect);
  let gh = gw / aspect;
  if (gh > height) { gh = height; gw = gh * aspect; }

  gridX0 = (width - gw) * 0.5;
  gridY0 = (height - gh) * 0.5;
  cellW = gw / COLS;
  cellH = gh / ROWS;
  invCellW = 1 / max(1e-6, cellW);
  invCellH = 1 / max(1e-6, cellH);
  Sampler.realloc();
  Particles.realloc();
  // NOTE: We intentionally do NOT reset the sampler cache on resize.
  // The sampler now computes density in grid-cell space (COLSxROWS), so cached frames remain valid.
}
function posToCell(x, y) {
  let c = ((x - gridX0) * invCellW) | 0, r = ((y - gridY0) * invCellH) | 0;
  if (c < 0) c = 0; else if (c >= COLS) c = COLS - 1;
  if (r < 0) r = 0; else if (r >= ROWS) r = ROWS - 1;
  return r * COLS + c;
}
function cellOriginX(cell) { return gridX0 + (cell % COLS) * cellW; }
function cellOriginY(cell) { return gridY0 + ((cell / COLS) | 0) * cellH; }

function fitRect(sw, sh, dw, dh, out) {
  const s = min(dw / sw, dh / sh), w = sw * s, h = sh * s;
  out.x = (dw - w) * 0.5; out.y = (dh - h) * 0.5; out.w = w; out.h = h; return out;
}
function hash01(n) {
  n = (n | 0) ^ 0x9e3779b9;
  n = Math.imul(n ^ (n >>> 16), 0x85ebca6b);
  n = Math.imul(n ^ (n >>> 13), 0xc2b2ae35);
  n = (n ^ (n >>> 16)) >>> 0;
  return n / 4294967296;
}
function nextActWithFrames(fromAct) {
  let a = fromAct;
  for (let k = 0; k < ACTS.length; k++) {
    a++; if (a > 4) a = 1;
    if ((SRC_COUNT[a] || 0) > 0) return a;
  }
  return fromAct;
}

function drawGridOverlay() {
  push();
  stroke(255, 22);
  strokeWeight(1);
  for (let c = 0; c <= COLS; c++) line(gridX0 + c * cellW, gridY0, gridX0 + c * cellW, gridY0 + ROWS * cellH);
  for (let r = 0; r <= ROWS; r++) line(gridX0, gridY0 + r * cellH, gridX0 + COLS * cellW, gridY0 + r * cellH);
  pop();
}

