// Particle-only p5.js (mobile-friendly). PNG frames are sampled offscreen only.

// ---- Config ----
const ACTS = [1, 2, 3, 4];
const SRC_COUNT = { 1: 192, 2: 172, 3: 192, 4: 0 };
const FPS_EFFECTIVE = { 1: 24, 2: 24, 3: 24, 4: 24 }; // source frames per sound-second

const TARGET_W = 160, TARGET_H = 284;
const CELL_SIZE = 18; // square grid size (px) (tune: larger -> fewer cells -> denser silhouettes)
const N = 4200, BLACK_PCT = 0.1;
const TRANSITION_DURATION = 2.0, DANCE_PORTION = 0.6;
const MIC_THRESHOLD = 0.03;
const AUTO_SPEED = 0.03;
const INTERP_SHARPNESS = 1.6; // >1 reduces "double exposure" trails during motion

// "Void travel" cleanup: when particles must cross empty (alpha=0) areas, guide them along 3 clean lanes
// (1 straight + 2 arced) between their current position and their destination.
const VOID_MIN_DIST_CELLS = 4.0;  // only engage lanes if far enough (in cell units)
const VOID_SAMPLES = 3;           // line samples to detect crossing empty space
const VOID_EMPTY_NEED = 2;        // how many samples must be empty to engage lanes
const LANE_LOOKAHEAD = 0.08;      // progress lookahead along lane (0..1)

// ---- State ----
let mic, amp, audioStarted = false;
let micRunning = false;
let t = 0, _prevT = 0, tAdvanced = false;
let debugOn = true, imagesDrawnToCanvas = false;
let showGrid = false;
let autoRun = false;
let catchUp = 0;
let debugMeanCellDist = 0;
let debugTransport = 0;
let actFrames = { 1: [], 2: [], 3: [], 4: [] };

let COLS = 1, ROWS = 1, CELLS = 1;
let gridX0 = 0, gridY0 = 0;
let cellW = 1, cellH = 1, invCellW = 1, invCellH = 1;

function resizeGrid() {
  cellW = CELL_SIZE;
  cellH = CELL_SIZE;
  COLS = max(1, floor(width / cellW));
  ROWS = max(1, floor(height / cellH));
  CELLS = COLS * ROWS;
  gridX0 = (width - COLS * cellW) * 0.5;
  gridY0 = (height - ROWS * cellH) * 0.5;
  invCellW = 1 / max(1e-6, cellW);
  invCellH = 1 / max(1e-6, cellH);
  Sampler.realloc();
  Particles.realloc();
  Sampler.resetAllCaches();
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
    if ((actFrames[a] || []).length > 0) return a;
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

// ---- Sampler: small buffer -> per-cell desired particle counts ----
const Sampler = {
  g: null, rG: { x: 0, y: 0, w: 0, h: 0 }, rC: { x: 0, y: 0, w: 0, h: 0 },
  w: TARGET_W, h: TARGET_H,
  weights: new Float32Array(1), frac: new Float32Array(1),
  cdf: new Float32Array(1),
  tmp: new Float32Array(1),
  pixCount: new Uint16Array(1),
  mask: new Uint8Array(1),
  visited: new Uint8Array(1),
  stack: new Int32Array(1),
  cache: { 1: null, 2: null, 3: null, 4: null },
  hits: 0, misses: 0, hitsF: 0, missesF: 0,
  init() {
    this.g = createGraphics(this.w, this.h);
    this.g.pixelDensity(1); this.g.noSmooth();
  },
  realloc() {
    this.weights = new Float32Array(CELLS);
    this.frac = new Float32Array(CELLS);
    this.cdf = new Float32Array(CELLS);
    this.tmp = new Float32Array(CELLS);
    this.pixCount = new Uint16Array(CELLS);
    this.mask = new Uint8Array(CELLS);
    this.visited = new Uint8Array(CELLS);
    this.stack = new Int32Array(CELLS);
  },
  resetFrameStats() { this.hitsF = 0; this.missesF = 0; },
  resetAllCaches() { for (let a of ACTS) this.cache[a] = null; },
  ensureActCache(a, cycle) {
    const need = cycle * CELLS, cur = this.cache[a];
    if (cur && cur.cycle === cycle && cur.counts.length === need) return;
    this.cache[a] = { cycle, done: new Uint8Array(cycle), counts: new Uint16Array(need) };
  },
  _weightAt(px) {
    const a = this.g.pixels[px + 3] / 255; if (a <= 0.05) return 0;
    const br = (this.g.pixels[px] + this.g.pixels[px + 1] + this.g.pixels[px + 2]) / 765;
    return br * a;
  },
  _compute(img, outCounts) {
    fitRect(img.width, img.height, this.w, this.h, this.rG);
    fitRect(img.width, img.height, width, height, this.rC);
    const rG = this.rG, rC = this.rC;

    this.g.clear(); this.g.background(0); this.g.imageMode(CORNER);
    this.g.image(img, rG.x, rG.y, rG.w, rG.h);
    this.g.loadPixels();

    this.weights.fill(0);
    this.pixCount.fill(0);
    const gwInv = 1 / max(1, rG.w), ghInv = 1 / max(1, rG.h);
    const x0 = max(0, floor(rG.x)), y0 = max(0, floor(rG.y));
    const x1 = min(this.w, ceil(rG.x + rG.w)), y1 = min(this.h, ceil(rG.y + rG.h));

    const gx1 = gridX0 + COLS * cellW;
    const gy1 = gridY0 + ROWS * cellH;

    for (let y = y0; y < y1; y++) {
      const row = y * this.w;
      const iy = (y + 0.5 - rG.y) * ghInv, cy = rC.y + iy * rC.h;
      if (cy < gridY0 || cy >= gy1) continue;
      for (let x = x0; x < x1; x++) {
        const p = 4 * (row + x);
        const w = this._weightAt(p);
        if (w <= 0) continue;
        const ix = (x + 0.5 - rG.x) * gwInv, cx = rC.x + ix * rC.w;
        if (cx < gridX0 || cx >= gx1) continue;
        const cell = posToCell(cx, cy);
        this.weights[cell] += w;
        if (this.pixCount[cell] < 65535) this.pixCount[cell]++;
      }
    }

    // Convert per-cell brightness to raw weights, normalized per-frame so gray-heavy frames still fill.
    // Then systematic-resample to exactly N particles.
    const BR_THR = 0.06;     // threshold in normalized brightness space (0..1)
    const GAMMA = 0.75;      // <1 boosts midtones so grays get particles
    const GREY_LIFT = 0.22;  // extra linear weight to further support midtones

    this.mask.fill(0);
    let maxBr = 0;
    for (let i = 0; i < CELLS; i++) {
      const pc = this.pixCount[i];
      if (pc === 0) { this.weights[i] = 0; continue; }
      const br = this.weights[i] / pc; // ~0..1
      this.weights[i] = br;
      if (br > maxBr) maxBr = br;
    }

    if (maxBr <= 1e-6) {
      outCounts.fill(0); outCounts[posToCell(width * 0.5, height * 0.5)] = N; return;
    }

    // Normalize and build weights
    let rawSum = 0;
    for (let i = 0; i < CELLS; i++) {
      const brN = this.weights[i] / maxBr;
      const x = (brN - BR_THR) / max(1e-6, 1 - BR_THR);
      if (x <= 0) { this.weights[i] = 0; continue; }
      const w = pow(x, GAMMA) + GREY_LIFT * x;
      this.weights[i] = w;
      this.mask[i] = 1;
      rawSum += w;
    }

    // Light smoothing to fill gray regions and reduce "missing rows" from sampling quantization
    this.tmp.fill(0);
    for (let r = 0; r < ROWS; r++) {
      for (let c = 0; c < COLS; c++) {
        const i = r * COLS + c;
        let sum = 0;
        let cnt = 0;
        for (let dr = -1; dr <= 1; dr++) {
          const rr = r + dr; if (rr < 0 || rr >= ROWS) continue;
          for (let dc = -1; dc <= 1; dc++) {
            const cc = c + dc; if (cc < 0 || cc >= COLS) continue;
            const j = rr * COLS + cc;
            const w = this.weights[j];
            if (w <= 0) continue;
            sum += w;
            cnt++;
          }
        }
        if (cnt > 0) this.tmp[i] = (this.weights[i] * 0.6) + (sum / cnt) * 0.4;
      }
    }
    rawSum = 0;
    for (let i = 0; i < CELLS; i++) { this.weights[i] = this.tmp[i]; rawSum += this.weights[i]; }

    if (rawSum <= 1e-9) {
      outCounts.fill(0); outCounts[posToCell(width * 0.5, height * 0.5)] = N; return;
    }

    // Systematic resampling: assigns exactly N particles with low variance and no scan-order bias
    outCounts.fill(0);
    let total = 0;
    for (let i = 0; i < CELLS; i++) { total += this.weights[i]; this.cdf[i] = total; }
    const step = total / N;
    let u = random(step);
    let i = 0;
    for (let k = 0; k < N; k++) {
      const th = u + k * step;
      while (i < CELLS - 1 && this.cdf[i] < th) i++;
      outCounts[i]++;
    }
  },
  ensure(a, idx, cycle, outOff) {
    this.ensureActCache(a, cycle);
    const c = this.cache[a]; if (!c || idx < 0 || idx >= c.cycle) return false;
    if (c.done[idx]) { this.hits++; this.hitsF++; outOff.off = idx * CELLS; return true; }
    const img = (actFrames[a] || [])[idx] || null; if (!img) return false;
    this.misses++; this.missesF++;
    const off = idx * CELLS;
    this._compute(img, c.counts.subarray(off, off + CELLS));
    c.done[idx] = 1; outOff.off = off; return true;
  },
  counts(a) { return this.cache[a] ? this.cache[a].counts : null; },
};

// ---- Particles: grid-density with local rebalancing + cheap separation ----
const Particles = {
  x: new Float32Array(N), y: new Float32Array(N), vx: new Float32Array(N), vy: new Float32Array(N),
  kind: new Uint8Array(N), cell: new Uint16Array(N), u: new Float32Array(N), v: new Float32Array(N),
  order: new Uint16Array(N),
  cellCounts: new Uint16Array(1), desired: new Uint16Array(1),
  // transport lanes (for clean motion across empty space)
  trA: new Uint8Array(N),
  trSX: new Float32Array(N), trSY: new Float32Array(N),
  trEX: new Float32Array(N), trEY: new Float32Array(N),
  trNX: new Float32Array(N), trNY: new Float32Array(N),
  trAmp: new Float32Array(N),
  // spatial hash
  binSize: 6, binsX: 1, binsY: 1, binHead: new Int32Array(1), binNext: new Int32Array(N),

  realloc() {
    this.cellCounts = new Uint16Array(CELLS);
    this.desired = new Uint16Array(CELLS);
  },
  resizeBins() {
    this.binsX = max(1, ceil(width / this.binSize));
    this.binsY = max(1, ceil(height / this.binSize));
    this.binHead = new Int32Array(this.binsX * this.binsY);
  },
  init() {
    if (this.cellCounts.length !== CELLS) this.realloc();
    this.cellCounts.fill(0);
    const cx = width * 0.5, cy = height * 0.5, R = min(width, height) * 0.22;
    for (let i = 0; i < N; i++) {
      this.order[i] = i;
      const ang = random(TWO_PI), rad = sqrt(random()) * R;
      this.x[i] = cx + cos(ang) * rad; this.y[i] = cy + sin(ang) * rad;
      this.vx[i] = random(-0.6, 0.6); this.vy[i] = random(-0.6, 0.6);
      this.kind[i] = random() < BLACK_PCT ? 1 : 0;
      const c = posToCell(this.x[i], this.y[i]);
      this.cell[i] = c; this.cellCounts[c] = min(65535, this.cellCounts[c] + 1);
      this.u[i] = random(); this.v[i] = random();
      this.trA[i] = 0;
    }
    for (let i = N - 1; i > 0; i--) {
      const j = floor(random(i + 1)), tmp = this.order[i];
      this.order[i] = this.order[j]; this.order[j] = tmp;
    }
  },
  setDesired(counts, off0, off1, alpha) {
    const a = constrain(alpha, 0, 1);
    let w1 = pow(a, INTERP_SHARPNESS);
    let w0 = pow(1 - a, INTERP_SHARPNESS);
    const ws = w0 + w1;
    if (ws > 1e-9) { w0 /= ws; w1 /= ws; } else { w0 = 1; w1 = 0; }
    let sum = 0;
    for (let i = 0; i < CELLS; i++) {
      const v = (counts[off0 + i] * w0 + counts[off1 + i] * w1 + 0.5) | 0;
      this.desired[i] = v;
      sum += v;
    }
    let diff = N - sum;
    if (diff > 0) { for (let i = 0; i < CELLS && diff > 0; i++) if (this.desired[i] > 0) { this.desired[i]++; diff--; } if (diff > 0) this.desired[posToCell(width * 0.5, height * 0.5)] += diff; }
    else if (diff < 0) { diff = -diff; for (let i = CELLS - 1; i >= 0 && diff > 0; i--) if (this.desired[i] > 0) { this.desired[i]--; diff--; } }
  },
  _bestNeighborDeficit(cell) {
    const r = (cell / COLS) | 0, c = cell - r * COLS;
    let best = -1, bestNeed = -1e9;

    // Ring 1 (8-neighborhood)
    for (let dr = -1; dr <= 1; dr++) {
      const rr = r + dr; if (rr < 0 || rr >= ROWS) continue;
      for (let dc = -1; dc <= 1; dc++) {
        if (dr === 0 && dc === 0) continue;
        const cc = c + dc; if (cc < 0 || cc >= COLS) continue;
        const nb = rr * COLS + cc;
        const need = (this.desired[nb] | 0) - (this.cellCounts[nb] | 0);
        if (need > bestNeed) { bestNeed = need; best = nb; }
      }
    }

    // Ring 2 (local hop) helps prevent particles getting stuck in "all-surplus" islands
    if (bestNeed <= 0) {
      for (let dr = -2; dr <= 2; dr++) {
        const rr = r + dr; if (rr < 0 || rr >= ROWS) continue;
        for (let dc = -2; dc <= 2; dc++) {
          if (dr === 0 && dc === 0) continue;
          if (abs(dr) <= 1 && abs(dc) <= 1) continue;
          const cc = c + dc; if (cc < 0 || cc >= COLS) continue;
          const nb = rr * COLS + cc;
          const need = (this.desired[nb] | 0) - (this.cellCounts[nb] | 0);
          if (need > bestNeed) { bestNeed = need; best = nb; }
        }
      }
    }

    // If no deficit nearby, take a guided 1-cell hop toward the strongest deficit in a wider radius.
    // This fixes "left behind" islands without teleporting across the screen.
    if (bestNeed <= 0) {
      if (_hotCount <= 0) return -1;
      let target = -1;
      let bestScore = 0;
      const maxK = min(_hotCount, 24);
      for (let k = 0; k < maxK; k++) {
        const cell2 = _hotCells[k];
        const need2 = _hotNeeds[k];
        if (cell2 < 0 || need2 <= 0) continue;
        const r2 = (cell2 / COLS) | 0;
        const c2 = cell2 - r2 * COLS;
        const dist = abs(r2 - r) + abs(c2 - c);
        const score = need2 / (1 + dist);
        if (score > bestScore) { bestScore = score; target = cell2; }
      }
      if (target < 0) return -1;
      const tr = (target / COLS) | 0;
      const tc = target - tr * COLS;
      const stepR = tr > r ? 1 : tr < r ? -1 : 0;
      const stepC = tc > c ? 1 : tc < c ? -1 : 0;
      const nr = constrain(r + stepR, 0, ROWS - 1);
      const nc = constrain(c + stepC, 0, COLS - 1);
      const nb = nr * COLS + nc;
      if (nb === cell) return -1;
      return nb;
    }
    return best;
  },
  rebalancePasses(passes) {
    let moved = 0;
    for (let p = 0; p < passes; p++) {
      for (let oi = 0; oi < N; oi++) {
        const i = this.order[oi], from = this.cell[i];
        if ((this.cellCounts[from] | 0) <= (this.desired[from] | 0)) continue;
        const to = this._bestNeighborDeficit(from); if (to < 0) continue;
        this.cellCounts[from]--; this.cellCounts[to]++; this.cell[i] = to; moved++;
        const h = (i + 1) * 73856093 ^ (to + 1) * 19349663;
        this.u[i] = hash01(h); this.v[i] = hash01(h ^ 0x68bc21eb);
        this.vx[i] *= 0.6; this.vy[i] *= 0.6;
      }
    }
    return moved;
  },
  hardRetargetToDesired() {
    // Assign cells deterministically to exactly match desired counts (prevents stuck stacks on act switch)
    this.cellCounts.set(this.desired);
    let cell = 0, cum = this.desired[0] | 0;
    for (let oi = 0; oi < N; oi++) {
      const i = this.order[oi], target = oi + 0.5;
      while (cum <= target && cell < CELLS - 1) { cell++; cum += this.desired[cell] | 0; }
      this.cell[i] = cell;
      const h = (i + 1) * 73856093 ^ (cell + 1) * 19349663;
      this.u[i] = hash01(h); this.v[i] = hash01(h ^ 0x68bc21eb);
      this.vx[i] = random(-0.4, 0.4); this.vy[i] = random(-0.4, 0.4);
    }
    for (let i = 0; i < N; i++) {
      const ox = cellOriginX(this.cell[i]), oy = cellOriginY(this.cell[i]);
      this.x[i] = ox + this.u[i] * cellW; this.y[i] = oy + this.v[i] * cellH;
    }
  },
  _buildBins() {
    this.binHead.fill(-1);
    for (let i = 0; i < N; i++) {
      const bx = constrain((this.x[i] / this.binSize) | 0, 0, this.binsX - 1);
      const by = constrain((this.y[i] / this.binSize) | 0, 0, this.binsY - 1);
      const b = by * this.binsX + bx;
      this.binNext[i] = this.binHead[b]; this.binHead[b] = i;
    }
  },
  separate(strength, radius) {
    this._buildBins();
    const r2 = radius * radius;
    for (let by = 0; by < this.binsY; by++) for (let bx = 0; bx < this.binsX; bx++) {
      const b = by * this.binsX + bx;
      for (let i = this.binHead[b]; i !== -1; i = this.binNext[i]) {
        for (let ny = max(0, by - 1); ny <= min(this.binsY - 1, by + 1); ny++) for (let nx = max(0, bx - 1); nx <= min(this.binsX - 1, bx + 1); nx++) {
          const bb = ny * this.binsX + nx;
          for (let j = this.binHead[bb]; j !== -1; j = this.binNext[j]) {
            if (j <= i) continue;
            const dx = this.x[i] - this.x[j], dy = this.y[i] - this.y[j], d2 = dx * dx + dy * dy;
            if (d2 > 1e-6 && d2 < r2) {
              const d = sqrt(d2), k = (1 - d / radius) * strength, ux = dx / d, uy = dy / d;
              this.vx[i] += ux * k; this.vy[i] += uy * k; this.vx[j] -= ux * k; this.vy[j] -= uy * k;
            }
          }
        }
      }
    }
  },
  updateInsideCells(level, snap01, flowAmt, catchUp01) {
    const cu = constrain(catchUp01 || 0, 0, 1);
    const damp = lerp(0.80, 0.72, cu);
    const maxSpd = lerp(3.6, 5.2, cu);
    const att = (0.10 + map(level, 0, 0.2, 0, 0.26, true)) * (1 + cu * 0.9);
    const wobbleAmp = (0.28 + map(level, 0, 0.2, 0, 0.55, true)) * (1 - cu * 0.75);
    const wobbleF = 0.35;
    const flow = 0.006 + map(level, 0, 0.2, 0, 0.03, true);
    const swirlBase = (0.04 + map(level, 0, 0.2, 0, 0.08, true)) * flowAmt * (1 - cu * 0.9);
    const cx = width * 0.5, cy = height * 0.5;

    debugTransport = 0;

    for (let i = 0; i < N; i++) {
      const cell = this.cell[i];
      const ox = cellOriginX(cell), oy = cellOriginY(cell);

      let snapX = 0, snapY = 0;
      let snapBase = 0.10, snapScale = 0.62, snapVelScale = 0.65;
      let hasSnap = false;

      // Base destination inside the assigned cell
      const baseTx = ox + this.u[i] * cellW;
      const baseTy = oy + this.v[i] * cellH;

      // Detect if we're crossing a large empty region (desired==0) between current position and destination.
      const posCell = posToCell(this.x[i], this.y[i]);
      const posEmpty = (this.desired[posCell] | 0) === 0;
      const dxBase = baseTx - this.x[i], dyBase = baseTy - this.y[i];
      const distBase = sqrt(dxBase * dxBase + dyBase * dyBase);

      let useLane = false;
      if (posEmpty && distBase > cellW * VOID_MIN_DIST_CELLS) {
        let emptyHits = 0;
        // Sample along the straight segment; if mostly empty, guide along a clean lane
        for (let s = 1; s <= VOID_SAMPLES; s++) {
          const tt = s / (VOID_SAMPLES + 1);
          const sx = this.x[i] + dxBase * tt;
          const sy = this.y[i] + dyBase * tt;
          const sc = posToCell(sx, sy);
          if ((this.desired[sc] | 0) === 0) emptyHits++;
        }
        useLane = emptyHits >= VOID_EMPTY_NEED;
      }

      // If lane is engaged, route via 1 of 3 deterministic trajectories: straight / arc+ / arc-
      if (useLane) {
        const lane = i % 3;
        const sign = lane === 1 ? 1 : lane === 2 ? -1 : 0;

        // Start lane segment on first entry (or if destination moved too far)
        if (!this.trA[i] || abs(baseTx - this.trEX[i]) + abs(baseTy - this.trEY[i]) > cellW * 6) {
          this.trA[i] = 1;
          this.trSX[i] = this.x[i]; this.trSY[i] = this.y[i];
          this.trEX[i] = baseTx; this.trEY[i] = baseTy;
          const dx = this.trEX[i] - this.trSX[i];
          const dy = this.trEY[i] - this.trSY[i];
          const len = sqrt(dx * dx + dy * dy) + 1e-6;
          this.trNX[i] = -dy / len;
          this.trNY[i] = dx / len;
          this.trAmp[i] = min(cellW * 7.5, len * 0.22);
        } else {
          // Follow destination smoothly so the lane remains stable
          this.trEX[i] = lerp(this.trEX[i], baseTx, 0.25);
          this.trEY[i] = lerp(this.trEY[i], baseTy, 0.25);
          const dx = this.trEX[i] - this.trSX[i];
          const dy = this.trEY[i] - this.trSY[i];
          const len = sqrt(dx * dx + dy * dy) + 1e-6;
          this.trNX[i] = -dy / len;
          this.trNY[i] = dx / len;
          this.trAmp[i] = min(cellW * 7.5, len * 0.22);
        }

        const sx = this.trSX[i], sy = this.trSY[i];
        const ex = this.trEX[i], ey = this.trEY[i];
        const dx = ex - sx, dy = ey - sy;
        const len2 = dx * dx + dy * dy + 1e-6;
        // progress along the segment based on projection
        let prog = ((this.x[i] - sx) * dx + (this.y[i] - sy) * dy) / len2;
        prog = constrain(prog, 0, 1);
        const prog2 = min(1, prog + LANE_LOOKAHEAD);

        const sin1 = sign === 0 ? 0 : sin(prog * PI);
        const sin2 = sign === 0 ? 0 : sin(prog2 * PI);

        const amp = this.trAmp[i] * sign;
        const nx = this.trNX[i], ny = this.trNY[i];

        const px = sx + dx * prog + nx * amp * sin1;
        const py = sy + dy * prog + ny * amp * sin1;
        const ax = sx + dx * prog2 + nx * amp * sin2;
        const ay = sy + dy * prog2 + ny * amp * sin2;

        // Pull toward lane and push forward along it (orderly trajectories)
        const pull = 0.030 + cu * 0.010;
        const push = 0.42 + cu * 0.28;
        const txLane = px;
        const tyLane = py;
        snapX = txLane; snapY = tyLane;
        snapBase = 0.08; snapScale = 0.55; snapVelScale = 0.55;
        hasSnap = true;
        let tdx = ax - px, tdy = ay - py;
        const tl = sqrt(tdx * tdx + tdy * tdy) + 1e-6;
        tdx /= tl; tdy /= tl;

        // Reduce noisy motions while in lane
        const laneSwirl = swirlBase * 0.15;
        const laneFlow = flow * 0.25;
        const ang = noise(this.x[i] * 0.0018, this.y[i] * 0.0018, t * 0.22) * TWO_PI * 2;
        const dxC = this.x[i] - cx, dyC = this.y[i] - cy, dC = sqrt(dxC * dxC + dyC * dyC) + 0.001;
        const swirl = laneSwirl * constrain(1 - dC / (min(width, height) * 0.7), 0, 1);
        this.vx[i] += (cos(ang) * laneFlow + (-dyC / dC) * swirl) * 0.18;
        this.vy[i] += (sin(ang) * laneFlow + (dxC / dC) * swirl) * 0.18;

        const dxT = txLane - this.x[i], dyT = tyLane - this.y[i];
        this.vx[i] += dxT * pull;
        this.vy[i] += dyT * pull;
        this.vx[i] += tdx * push * 0.03;
        this.vy[i] += tdy * push * 0.03;

        // Rejoin normal behavior near the destination
        if (prog > 0.92 || !posEmpty) this.trA[i] = 0;
        else debugTransport++;

      } else {
        this.trA[i] = 0;

        const nx1 = noise((i + 1) * 0.017, (i + 7) * 0.031, t * wobbleF);
        const ny1 = noise((i + 9) * 0.019, (i + 3) * 0.029, t * wobbleF);
        const tx = baseTx + (nx1 - 0.5) * 2 * wobbleAmp;
        const ty = baseTy + (ny1 - 0.5) * 2 * wobbleAmp;
        snapX = tx; snapY = ty;
        snapBase = 0.10; snapScale = 0.62; snapVelScale = 0.65;
        hasSnap = true;

        const ang = noise(this.x[i] * 0.0018, this.y[i] * 0.0018, t * 0.22) * TWO_PI * 2;
        const dxC = this.x[i] - cx, dyC = this.y[i] - cy, dC = sqrt(dxC * dxC + dyC * dyC) + 0.001;
        const swirl = swirlBase * constrain(1 - dC / (min(width, height) * 0.7), 0, 1);
        this.vx[i] += (cos(ang) * flow + (-dyC / dC) * swirl) * 0.18;
        this.vy[i] += (sin(ang) * flow + (dxC / dC) * swirl) * 0.18;

        const dx = tx - this.x[i], dy = ty - this.y[i];
        this.vx[i] += dx * att; this.vy[i] += dy * att;

        // Small continuous micro-motion inside the cell (reduced during catch-up so it doesn't leave trails)
        const jig = (0.03 + map(level, 0, 0.2, 0, 0.06, true)) * (1 - cu * 0.8);
        this.vx[i] += (noise(i * 0.013, 9.1, t * 0.5) - 0.5) * jig;
        this.vy[i] += (noise(i * 0.017, 3.7, t * 0.5) - 0.5) * jig;
      }

      this.vx[i] *= damp; this.vy[i] *= damp;

      const spd = sqrt(this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i]);
      if (spd > maxSpd) { const k = maxSpd / max(0.001, spd); this.vx[i] *= k; this.vy[i] *= k; }
      this.x[i] += this.vx[i]; this.y[i] += this.vy[i];

      // No hard cell clamping (it creates visible empty grid lines). Just keep within canvas bounds.
      if (this.x[i] < 0) { this.x[i] = 0; this.vx[i] *= -0.25; }
      else if (this.x[i] > width) { this.x[i] = width; this.vx[i] *= -0.25; }
      if (this.y[i] < 0) { this.y[i] = 0; this.vy[i] *= -0.25; }
      else if (this.y[i] > height) { this.y[i] = height; this.vy[i] *= -0.25; }

      if (snap01 > 0 && hasSnap) {
        const k = snapBase + snap01 * snapScale;
        this.x[i] = lerp(this.x[i], snapX, k);
        this.y[i] = lerp(this.y[i], snapY, k);
        const vk = 1 - snap01 * snapVelScale;
        this.vx[i] *= vk;
        this.vy[i] *= vk;
      }
    }
  },
  dance(level, s) {
    const flow = 0.20 + map(level, 0, 0.2, 0, 0.65, true);
    const swirlBase = 0.85 + map(level, 0, 0.2, 0, 1.15, true);
    const driftX = lerp(0.25, 0.05, s), driftY = lerp(-0.05, -0.15, s);
    const damp = 0.92, maxSpd = 5.8, cx = width * 0.5, cy = height * 0.5;
    for (let i = 0; i < N; i++) {
      const ang = noise(this.x[i] * 0.0021, this.y[i] * 0.0021, t * 0.35) * TWO_PI * 2.4;
      const dx = this.x[i] - cx, dy = this.y[i] - cy, d = sqrt(dx * dx + dy * dy) + 0.001;
      const swirl = swirlBase * constrain(1 - d / (min(width, height) * 0.55), 0, 1);
      this.vx[i] += (cos(ang) * flow + driftX + (-dy / d) * swirl) * 0.18;
      this.vy[i] += (sin(ang) * flow + driftY + (dx / d) * swirl) * 0.18;
      this.vx[i] *= damp; this.vy[i] *= damp;
      const spd = sqrt(this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i]);
      if (spd > maxSpd) { const k = maxSpd / max(0.001, spd); this.vx[i] *= k; this.vy[i] *= k; }
      this.x[i] += this.vx[i]; this.y[i] += this.vy[i];
      if (this.x[i] < 0) { this.x[i] = 0; this.vx[i] *= -0.25; } else if (this.x[i] > width) { this.x[i] = width; this.vx[i] *= -0.25; }
      if (this.y[i] < 0) { this.y[i] = 0; this.vy[i] *= -0.25; } else if (this.y[i] > height) { this.y[i] = height; this.vy[i] *= -0.25; }
    }
  },
  draw() {
    noStroke();
    for (let i = 0; i < N; i++) { this.kind[i] ? fill(0, 210) : fill(255, 235); circle(this.x[i], this.y[i], 2.5); }
  },
};

function estimateMismatch() {
  // Sum of surplus particles over all cells; higher means more redistribution needed
  let surplus = 0;
  for (let i = 0; i < CELLS; i++) {
    const d = (Particles.cellCounts[i] | 0) - (Particles.desired[i] | 0);
    if (d > 0) surplus += d;
  }
  return surplus;
}

function estimateMeanCellDist() {
  const samples = min(220, N);
  const stride = max(1, (N / samples) | 0);
  let sum = 0;
  for (let s = 0; s < samples; s++) {
    const i = s * stride;
    const cell = Particles.cell[i];
    const ox = cellOriginX(cell), oy = cellOriginY(cell);
    const tx = ox + Particles.u[i] * cellW;
    const ty = oy + Particles.v[i] * cellH;
    const dx = Particles.x[i] - tx, dy = Particles.y[i] - ty;
    sum += sqrt(dx * dx + dy * dy);
  }
  return sum / samples;
}

const _hotCells = new Int32Array(48);
const _hotNeeds = new Int32Array(48);
let _hotCount = 0;

function computeDeficitHotspots() {
  // Top-k deficit cells by (desired - current). O(CELLS * K) with small K.
  const K = _hotCells.length;
  let count = 0;
  for (let i = 0; i < K; i++) { _hotCells[i] = -1; _hotNeeds[i] = 0; }
  for (let cell = 0; cell < CELLS; cell++) {
    const need = (Particles.desired[cell] | 0) - (Particles.cellCounts[cell] | 0);
    if (need <= 0) continue;
    let j = count < K ? count++ : K - 1;
    if (count === K && need <= _hotNeeds[j]) continue;
    while (j > 0 && need > _hotNeeds[j - 1]) {
      _hotNeeds[j] = _hotNeeds[j - 1];
      _hotCells[j] = _hotCells[j - 1];
      j--;
    }
    _hotNeeds[j] = need;
    _hotCells[j] = cell;
  }
  _hotCount = count;
}

function catchUpNudge(strength01) {
  if (_hotCount <= 0) return;
  const cu = constrain(strength01, 0, 1);
  // Gentle "hard catch-up" (no teleport): keep this subtle to avoid feeling too fast.
  const kick = 0.00025 * (0.20 + cu * 0.60);
  const stride = 10; // much fewer nudged particles per update
  const seed = (t * 60) | 0;
  for (let i = 0; i < N; i += stride) {
    const from = Particles.cell[i];
    const surplus = (Particles.cellCounts[from] | 0) - (Particles.desired[from] | 0);
    if (surplus <= 0) continue;
    const idx = (i + seed) % _hotCount;
    const targetCell = _hotCells[idx];
    const tx = cellOriginX(targetCell) + cellW * 0.5;
    const ty = cellOriginY(targetCell) + cellH * 0.5;
    const dx = tx - Particles.x[i];
    const dy = ty - Particles.y[i];
    Particles.vx[i] += dx * kick;
    Particles.vy[i] += dy * kick;
  }
}

// ---- Act controller (single source of truth for timing) ----
const Acts = {
  mode: "ACT", act: 1, next: 2, actStartT: 0, transitionStartT: 0,
  elapsed: 0, dur: 0, fps: 24, cycle: 0, frameFloat: 0, idx: 0, src0: 0, src1: 0, alpha: 0, cycles: 0,
  compute() {
    const loaded = (actFrames[this.act] || []).length;
    const src = SRC_COUNT[this.act] || loaded;
    this.cycle = min(loaded, src);
    this.fps = FPS_EFFECTIVE[this.act] || 24;
    this.elapsed = max(0, t - this.actStartT);
    this.frameFloat = this.elapsed * this.fps;
    this.dur = this.cycle > 0 && this.fps > 0 ? this.cycle / this.fps : 0;
    if (this.cycle <= 0) { this.idx = this.src0 = this.src1 = 0; this.alpha = 0; this.cycles = 0; return; }
    const base = floor(this.frameFloat);
    this.idx = base % this.cycle;
    this.src0 = this.idx;
    this.src1 = (this.idx + 1) % this.cycle;
    this.alpha = constrain(this.frameFloat - base, 0, 1);
    this.cycles = floor(this.frameFloat / this.cycle);
  },
  startTransition() {
    this.mode = "TRANSITION";
    this.transitionStartT = t;
    this.next = nextActWithFrames(this.act);
  },
  endTransition() {
    this.act = this.next;
    this.actStartT = t;
    this.mode = "ACT";
    this.compute();
    if (this.cycle > 0) {
      Sampler.resetFrameStats();
      if (Sampler.ensure(this.act, 0, this.cycle, _off0)) {
        const counts = Sampler.counts(this.act);
        Particles.setDesired(counts, _off0.off, _off0.off, 0);
        Particles.hardRetargetToDesired();
      }
    }
  },
  update() {
    this.compute();
    if (this.mode === "ACT") {
      if (this.cycle <= 0 || this.dur <= 0) return this.startTransition();
      if (this.cycles >= 1 || this.elapsed >= this.dur) return this.startTransition();
    } else {
      if (t - this.transitionStartT >= TRANSITION_DURATION) return this.endTransition();
    }
  },
};

const _off0 = { off: 0 }, _off1 = { off: 0 };

// ---- p5 lifecycle ----
function preload() {
  for (let a of ACTS) {
    const count = SRC_COUNT[a] || 0;
    actFrames[a] = new Array(count);
    for (let i = 0; i < count; i++) actFrames[a][i] = loadImage(`assets/act${a}/act${a}_${nf(i, 5)}.png`);
  }
}

function setup() {
  createCanvas(windowWidth, windowHeight);
  mic = new p5.AudioIn();
  amp = new p5.Amplitude();
  resizeGrid();
  Sampler.init();
  Particles.resizeBins();
  Particles.init();

  const btnGrid = createButton("Grid: OFF");
  btnGrid.position(12, 12);
  btnGrid.style("position", "fixed");
  btnGrid.style("z-index", "10");
  btnGrid.style("padding", "8px 10px");
  btnGrid.style("font-size", "14px");
  btnGrid.style("border", "1px solid rgba(255,255,255,0.25)");
  btnGrid.style("background", "rgba(0,0,0,0.55)");
  btnGrid.style("color", "#fff");
  btnGrid.style("border-radius", "8px");
  btnGrid.mousePressed(() => {
    showGrid = !showGrid;
    btnGrid.html(showGrid ? "Grid: ON" : "Grid: OFF");
  });

  const btnAuto = createButton("Auto: OFF");
  btnAuto.position(110, 12);
  btnAuto.style("position", "fixed");
  btnAuto.style("z-index", "10");
  btnAuto.style("padding", "8px 10px");
  btnAuto.style("font-size", "14px");
  btnAuto.style("border", "1px solid rgba(255,255,255,0.25)");
  btnAuto.style("background", "rgba(0,0,0,0.55)");
  btnAuto.style("color", "#fff");
  btnAuto.style("border-radius", "8px");
  btnAuto.mousePressed(() => {
    autoRun = !autoRun;
    if (autoRun) audioStarted = true;
    btnAuto.html(autoRun ? "Auto: ON" : "Auto: OFF");
  });
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
  resizeGrid();
  Particles.resizeBins();
  Particles.init();
  Sampler.resetAllCaches();
}

function draw() {
  background(0);
  if (!audioStarted && !autoRun) return drawStartScreen();

  const level = autoRun ? 0.09 : (micRunning ? amp.getLevel() : 0);
  if (autoRun) t += AUTO_SPEED;
  else if (micRunning && level > MIC_THRESHOLD) t += map(level, MIC_THRESHOLD, 0.2, 0.01, 0.05, true);
  tAdvanced = t !== _prevT; _prevT = t;

  Sampler.resetFrameStats();
  Acts.update();

  let desiredSum = 0, moved = 0;
  let mismatch = 0;

  if (Acts.mode === "ACT") {
    if (Acts.cycle > 0) {
      const ok0 = Sampler.ensure(Acts.act, Acts.src0, Acts.cycle, _off0);
      const ok1 = Sampler.ensure(Acts.act, Acts.src1, Acts.cycle, _off1);
      if (ok0 || ok1) {
        const counts = Sampler.counts(Acts.act);
        if (ok0 && ok1) Particles.setDesired(counts, _off0.off, _off1.off, Acts.alpha);
        else if (ok0) Particles.setDesired(counts, _off0.off, _off0.off, 0);
        else Particles.setDesired(counts, _off1.off, _off1.off, 0);
        for (let i = 0; i < CELLS; i++) desiredSum += Particles.desired[i];
        if (tAdvanced) {
          mismatch = estimateMismatch();
          debugMeanCellDist = estimateMeanCellDist();
          // Position-based catch-up helps prevent "stuck in space", but keep it mild.
          const posCU = constrain((debugMeanCellDist - cellW * 0.85) / (cellW * 3.2), 0, 1);
          const targetCU = constrain((mismatch - N * 0.10) / (N * 0.45), 0, 1);
          catchUp = max(catchUp * 0.94, targetCU, posCU * 0.35);
          const passes = (mismatch > 900 ? 10 : mismatch > 600 ? 8 : mismatch > 350 ? 6 : mismatch > 180 ? 4 : 3) + floor(catchUp * 2);
          // IMPORTANT: compute hotspots BEFORE rebalancing so surplus "islands" can find long-range deficits
          // in the same update tick (otherwise particles may look "stuck" until the next t advance).
          computeDeficitHotspots();
          moved += Particles.rebalancePasses(passes);
          if (catchUp > 0.32) catchUpNudge(catchUp);
          // Catch-up increases snap/attraction but stays smooth (no teleport)
          Particles.updateInsideCells(level, 0.18 + catchUp * 0.26, 0.15, catchUp * 0.7);
          Particles.separate(0.018, 3.2);
          catchUp *= 0.97;
        }
      }
    }
    Particles.draw();
  } else {
    const p = constrain((t - Acts.transitionStartT) / max(1e-6, TRANSITION_DURATION), 0, 1);
    const stageA = p < DANCE_PORTION;
    const s = stageA ? p / max(1e-6, DANCE_PORTION) : (p - DANCE_PORTION) / max(1e-6, 1 - DANCE_PORTION);

    if (tAdvanced) {
      if (stageA) {
        Particles.dance(level, s);
      } else {
        const a = Acts.next;
        const loaded = (actFrames[a] || []).length;
        const cycle = min(loaded, SRC_COUNT[a] || loaded);
        if (cycle > 0) {
          const ok0 = Sampler.ensure(a, 0, cycle, _off0);
          const ok1 = Sampler.ensure(a, 1 % cycle, cycle, _off1);
          if (ok0 && ok1) {
            const counts = Sampler.counts(a);
            Particles.setDesired(counts, _off0.off, _off1.off, s);
            for (let i = 0; i < CELLS; i++) desiredSum += Particles.desired[i];
            debugMeanCellDist = estimateMeanCellDist();
            computeDeficitHotspots();
            moved += Particles.rebalancePasses(7);
            if (p > 0.74) catchUpNudge(0.25);
            Particles.updateInsideCells(level, map(p, 0.88, 1.0, 0, 1, true), 0.0, 0.65);
            Particles.separate(0.020, 3.2);
          }
        }
      }
    }
    Particles.draw();
  }

  if (debugOn) drawDebug(level, desiredSum, moved, mismatch);
  if (showGrid) drawGridOverlay();
}

function drawStartScreen() {
  fill(255); noStroke(); textAlign(CENTER, CENTER); textSize(18);
  text("Tap to activate microphone\nSound drives the animation", width / 2, height / 2);
}

function mousePressed() {
  if (audioStarted) return;
  userStartAudio();
  mic.start(() => { amp.setInput(mic); audioStarted = true; micRunning = true; Acts.actStartT = t; });
}

function keyPressed() { if (key === "d" || key === "D") debugOn = !debugOn; }

function drawDebug(level, desiredSum, moved, mismatch) {
  const loaded = (actFrames[Acts.act] || []).length;
  push(); noStroke(); fill(0, 160); rect(10, 10, 470, 106, 8); fill(255); textSize(12); textAlign(LEFT, TOP);
  text(`mode:${Acts.mode} act:${Acts.act} next:${Acts.next}  t:${t.toFixed(2)} el:${Acts.elapsed.toFixed(2)} dur:${Acts.dur.toFixed(2)} fps:${Acts.fps}`, 18, 16);
  text(`cycle:${Acts.cycle} SRC:${SRC_COUNT[Acts.act] || 0} loaded:${loaded}  src0/src1:${Acts.src0}/${Acts.src1} a:${Acts.alpha.toFixed(2)} cycles:${Acts.cycles}`, 18, 32);
  text(`grid:${COLS}x${ROWS} desiredSum:${desiredSum} moved:${moved} mismatch:${mismatch} meanDist:${debugMeanCellDist.toFixed(1)} catchUp:${catchUp.toFixed(2)} hot:${_hotCount} lane:${debugTransport}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} level:${level.toFixed(3)}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  pop();
}
