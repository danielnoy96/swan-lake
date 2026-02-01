// Particle-only p5.js (mobile-friendly). PNG frames are sampled offscreen only.

// ---- Config ----
const ACTS = [1, 2, 3, 4];
const SRC_COUNT = { 1: 192, 2: 172, 3: 192, 4: 0 };
const FPS_EFFECTIVE = { 1: 24, 2: 24, 3: 24, 4: 24 }; // source frames per sound-second

const TARGET_W = 160, TARGET_H = 284;
const CELL_SIZE = 18; // square grid size (px) (tune: larger -> fewer cells -> denser silhouettes)
const N = 3200, BLACK_PCT = 0.1;
const TRANSITION_DURATION = 2.0, DANCE_PORTION = 0.6;
const MIC_THRESHOLD = 0.03;
const AUTO_SPEED = 0.03;

// ---- State ----
let mic, amp, audioStarted = false;
let micRunning = false;
let t = 0, _prevT = 0, tAdvanced = false;
let debugOn = true, imagesDrawnToCanvas = false;
let showGrid = false;
let autoRun = false;
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
    const b = br - 0.08;
    if (b <= 0) return 0;
    return b * a;
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

    // Convert per-cell brightness to raw weights (later resampled to exactly N particles)
    // This keeps particles spread across all bright cells instead of collapsing to a few.
    const BR_THR = 0.10;
    const GAMMA = 1.10;

    this.mask.fill(0);
    let rawSum = 0;
    for (let i = 0; i < CELLS; i++) {
      const pc = this.pixCount[i];
      if (pc === 0) { this.weights[i] = 0; continue; }
      let br = this.weights[i] / pc; // ~0..1
      br = (br - BR_THR) / max(1e-6, 1 - BR_THR);
      if (br <= 0) { this.weights[i] = 0; continue; }
      br = pow(br, GAMMA);
      const raw = br;
      this.weights[i] = raw;
      this.mask[i] = 1;
      rawSum += raw;
    }

    // Remove isolated single-cell noise (one pass)
    for (let i = 0; i < CELLS; i++) {
      if (this.mask[i] === 0) continue;
      const r = (i / COLS) | 0;
      const c = i - r * COLS;
      let n = 0;
      for (let dr = -1; dr <= 1; dr++) {
        const rr = r + dr; if (rr < 0 || rr >= ROWS) continue;
        for (let dc = -1; dc <= 1; dc++) {
          if (dr === 0 && dc === 0) continue;
          const cc = c + dc; if (cc < 0 || cc >= COLS) continue;
          n += this.mask[rr * COLS + cc];
        }
      }
      if (n < 2) { rawSum -= this.weights[i]; this.weights[i] = 0; this.mask[i] = 0; }
    }

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
    }
    for (let i = N - 1; i > 0; i--) {
      const j = floor(random(i + 1)), tmp = this.order[i];
      this.order[i] = this.order[j]; this.order[j] = tmp;
    }
  },
  rebuildCellsFromPos() {
    this.cellCounts.fill(0);
    for (let i = 0; i < N; i++) {
      const c = posToCell(this.x[i], this.y[i]);
      this.cell[i] = c; this.cellCounts[c] = min(65535, this.cellCounts[c] + 1);
    }
  },
  setDesired(counts, off0, off1, alpha) {
    const a = constrain(alpha, 0, 1), ia = 1 - a;
    let sum = 0;
    for (let i = 0; i < CELLS; i++) { const v = (counts[off0 + i] * ia + counts[off1 + i] * a + 0.5) | 0; this.desired[i] = v; sum += v; }
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
      let target = -1;
      let bestScore = 0;
      const R = 10;
      for (let dr = -R; dr <= R; dr++) {
        const rr = r + dr; if (rr < 0 || rr >= ROWS) continue;
        for (let dc = -R; dc <= R; dc++) {
          if (dr === 0 && dc === 0) continue;
          const cc = c + dc; if (cc < 0 || cc >= COLS) continue;
          const nb = rr * COLS + cc;
          const score = (this.desired[nb] | 0) - (this.cellCounts[nb] | 0);
          if (score > bestScore) { bestScore = score; target = nb; }
        }
      }
      if (target < 0 || bestScore <= 0) return -1;
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
  updateInsideCells(level, snap01, flowAmt) {
    const damp = 0.78, maxSpd = 3.6;
    const att = 0.10 + map(level, 0, 0.2, 0, 0.26, true);
    const wobbleAmp = 0.28 + map(level, 0, 0.2, 0, 0.55, true), wobbleF = 0.35;
    const flow = 0.006 + map(level, 0, 0.2, 0, 0.03, true);
    const swirlBase = (0.04 + map(level, 0, 0.2, 0, 0.08, true)) * flowAmt;
    const cx = width * 0.5, cy = height * 0.5;

    for (let i = 0; i < N; i++) {
      const cell = this.cell[i];
      const ox = cellOriginX(cell), oy = cellOriginY(cell);
      const nx1 = noise((i + 1) * 0.017, (i + 7) * 0.031, t * wobbleF);
      const ny1 = noise((i + 9) * 0.019, (i + 3) * 0.029, t * wobbleF);
      const tx = ox + this.u[i] * cellW + (nx1 - 0.5) * 2 * wobbleAmp;
      const ty = oy + this.v[i] * cellH + (ny1 - 0.5) * 2 * wobbleAmp;

      const ang = noise(this.x[i] * 0.0018, this.y[i] * 0.0018, t * 0.22) * TWO_PI * 2;
      const dxC = this.x[i] - cx, dyC = this.y[i] - cy, dC = sqrt(dxC * dxC + dyC * dyC) + 0.001;
      const swirl = swirlBase * constrain(1 - dC / (min(width, height) * 0.7), 0, 1);
      this.vx[i] += (cos(ang) * flow + (-dyC / dC) * swirl) * 0.18;
      this.vy[i] += (sin(ang) * flow + (dxC / dC) * swirl) * 0.18;

      const dx = tx - this.x[i], dy = ty - this.y[i];
      this.vx[i] += dx * att; this.vy[i] += dy * att;

      // Small continuous micro-motion inside the cell (prevents "stuck" feeling even when shape is static)
      const jig = 0.03 + map(level, 0, 0.2, 0, 0.06, true);
      this.vx[i] += (noise(i * 0.013, 9.1, t * 0.5) - 0.5) * jig;
      this.vy[i] += (noise(i * 0.017, 3.7, t * 0.5) - 0.5) * jig;

      this.vx[i] *= damp; this.vy[i] *= damp;

      const spd = sqrt(this.vx[i] * this.vx[i] + this.vy[i] * this.vy[i]);
      if (spd > maxSpd) { const k = maxSpd / max(0.001, spd); this.vx[i] *= k; this.vy[i] *= k; }
      this.x[i] += this.vx[i]; this.y[i] += this.vy[i];

      // No hard cell clamping (it creates visible empty grid lines). Just keep within canvas bounds.
      if (this.x[i] < 0) { this.x[i] = 0; this.vx[i] *= -0.25; }
      else if (this.x[i] > width) { this.x[i] = width; this.vx[i] *= -0.25; }
      if (this.y[i] < 0) { this.y[i] = 0; this.vy[i] *= -0.25; }
      else if (this.y[i] > height) { this.y[i] = height; this.vy[i] *= -0.25; }

      if (snap01 > 0) {
        const k = 0.10 + snap01 * 0.62;
        this.x[i] = lerp(this.x[i], tx, k); this.y[i] = lerp(this.y[i], ty, k);
        this.vx[i] *= 1 - snap01 * 0.65; this.vy[i] *= 1 - snap01 * 0.65;
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
    for (let i = 0; i < N; i++) { this.kind[i] ? fill(0, 210) : fill(255, 235); circle(this.x[i], this.y[i], 2.1); }
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
    if (tAdvanced) Particles.rebuildCellsFromPos();
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

  if (Acts.mode === "ACT") {
    if (Acts.cycle > 0) {
      const ok0 = Sampler.ensure(Acts.act, Acts.src0, Acts.cycle, _off0);
      const ok1 = Sampler.ensure(Acts.act, Acts.src1, Acts.cycle, _off1);
      if (ok0 && ok1) {
        const counts = Sampler.counts(Acts.act);
        Particles.setDesired(counts, _off0.off, _off1.off, Acts.alpha);
        for (let i = 0; i < CELLS; i++) desiredSum += Particles.desired[i];
        if (tAdvanced) {
          const mismatch = estimateMismatch();
          const passes = mismatch > 900 ? 10 : mismatch > 600 ? 8 : mismatch > 350 ? 6 : mismatch > 180 ? 4 : 3;
          moved += Particles.rebalancePasses(passes);
          // Mild snap in ACT improves clarity while still allowing motion inside cells
          Particles.updateInsideCells(level, 0.18, 0.15);
          Particles.separate(0.018, 3.2);
          Particles.rebuildCellsFromPos();
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
        Particles.rebuildCellsFromPos();
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
            moved += Particles.rebalancePasses(7);
            Particles.updateInsideCells(level, map(p, 0.88, 1.0, 0, 1, true), 0.0);
            Particles.separate(0.020, 3.2);
            Particles.rebuildCellsFromPos();
          }
        }
      }
    }
    Particles.draw();
  }

  if (debugOn) drawDebug(level, desiredSum, moved);
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

function drawDebug(level, desiredSum, moved) {
  const loaded = (actFrames[Acts.act] || []).length;
  push(); noStroke(); fill(0, 160); rect(10, 10, 470, 106, 8); fill(255); textSize(12); textAlign(LEFT, TOP);
  text(`mode:${Acts.mode} act:${Acts.act} next:${Acts.next}  t:${t.toFixed(2)} el:${Acts.elapsed.toFixed(2)} dur:${Acts.dur.toFixed(2)} fps:${Acts.fps}`, 18, 16);
  text(`cycle:${Acts.cycle} SRC:${SRC_COUNT[Acts.act] || 0} loaded:${loaded}  src0/src1:${Acts.src0}/${Acts.src1} a:${Acts.alpha.toFixed(2)} cycles:${Acts.cycles}`, 18, 32);
  text(`grid:${COLS}x${ROWS} desiredSum:${desiredSum} moved:${moved}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} level:${level.toFixed(3)}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  pop();
}
