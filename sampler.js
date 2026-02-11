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
  // On-demand loader (avoids preloading hundreds of PNGs -> crashes/reload loops on mobile)
  queue: [],
  inFlight: 0,
  maxInFlight: 3,
  hits: 0, misses: 0, hitsF: 0, missesF: 0,
  init() {
    this.g = createGraphics(this.w, this.h);
    this.g.pixelDensity(1); this.g.noSmooth();
    try { if (this.g.drawingContext) this.g.drawingContext.imageSmoothingEnabled = false; } catch (_) {}
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
  resetAllCaches() {
    for (let a of ACTS) this.cache[a] = null;
    this.queue.length = 0;
  },
  ensureActCache(a, cycle) {
    const need = cycle * CELLS, cur = this.cache[a];
    if (cur && cur.cycle === cycle && cur.counts.length === need) return;
    this.cache[a] = {
      cycle,
      done: new Uint8Array(cycle),
      loading: new Uint8Array(cycle),
      failed: new Uint8Array(cycle),
      counts: new Uint16Array(need),
      pixHash: new Uint32Array(cycle),
      countsHash: new Uint32Array(cycle),
      imgW: new Uint16Array(cycle),
      imgH: new Uint16Array(cycle),
    };
  },
  _fnv1aBytes(arr) {
    let h = 2166136261 >>> 0;
    for (let i = 0; i < arr.length; i++) {
      h ^= arr[i] & 0xff;
      h = Math.imul(h, 16777619) >>> 0;
    }
    return h >>> 0;
  },
  _fnv1aU16(arr, stride) {
    let h = 2166136261 >>> 0;
    for (let i = 0; i < arr.length; i += stride) {
      h ^= arr[i] & 0xffff;
      h = Math.imul(h, 16777619) >>> 0;
    }
    return h >>> 0;
  },
  _weightAt(px) {
    const a = this.g.pixels[px + 3] / 255; if (a <= 0.05) return 0;
    const br = (this.g.pixels[px] + this.g.pixels[px + 1] + this.g.pixels[px + 2]) / 765;
    return br * a;
  },
  _compute(img, outCounts, seed01) {
    fitRect(img.width, img.height, this.w, this.h, this.rG);
    // IMPORTANT: map the source into *grid cell space* (not canvas pixels).
    // This makes cell assignment deterministic across hosts/viewports/DPR.
    fitRect(img.width, img.height, COLS, ROWS, this.rC);
    const rG = this.rG, rC = this.rC;

    // Quantize draw rect to integer pixels in the small buffer to avoid subpixel raster differences.
    // (Subpixel scaling can produce slightly different results across hosts/GPUs even in the same browser.)
    rG.x = Math.round(rG.x);
    rG.y = Math.round(rG.y);
    rG.w = Math.max(1, Math.round(rG.w));
    rG.h = Math.max(1, Math.round(rG.h));
    if (rG.x < 0) rG.x = 0;
    if (rG.y < 0) rG.y = 0;
    if (rG.x + rG.w > this.w) rG.w = Math.max(1, this.w - rG.x);
    if (rG.y + rG.h > this.h) rG.h = Math.max(1, this.h - rG.y);

    this.g.clear(); this.g.background(0); this.g.imageMode(CORNER);
    // Ensure nearest-neighbor in this buffer for stable sampling.
    try { this.g.noSmooth(); } catch (_) {}
    try {
      if (this.g.drawingContext) {
        this.g.drawingContext.imageSmoothingEnabled = false;
        this.g.drawingContext.globalCompositeOperation = "source-over";
        this.g.drawingContext.globalAlpha = 1;
        if (this.g.drawingContext.setTransform) this.g.drawingContext.setTransform(1, 0, 0, 1, 0, 0);
      }
    } catch (_) {}
    this.g.image(img, rG.x, rG.y, rG.w, rG.h);
    this.g.loadPixels();

    // Diagnostic: stable hash of the rasterized small buffer pixels.
    // If this differs across hosts for the same (act, src), the divergence is in PNG decode/draw/rasterization.
    const pixHash = this._fnv1aBytes(this.g.pixels);

    this.weights.fill(0);
    this.pixCount.fill(0);
    const gwInv = 1 / max(1, rG.w), ghInv = 1 / max(1, rG.h);
    const x0 = max(0, floor(rG.x)), y0 = max(0, floor(rG.y));
    const x1 = min(this.w, ceil(rG.x + rG.w)), y1 = min(this.h, ceil(rG.y + rG.h));

    for (let y = y0; y < y1; y++) {
      const row = y * this.w;
      const iy = (y + 0.5 - rG.y) * ghInv;
      const rCell = rC.y + iy * rC.h;
      if (rCell < 0 || rCell >= ROWS) continue;
      const rr = rCell | 0;
      for (let x = x0; x < x1; x++) {
        const p = 4 * (row + x);
        const w = this._weightAt(p);
        if (w <= 0) continue;
        const ix = (x + 0.5 - rG.x) * gwInv;
        const cCell = rC.x + ix * rC.w;
        if (cCell < 0 || cCell >= COLS) continue;
        const cell = rr * COLS + (cCell | 0);
        this.weights[cell] += w;
        if (this.pixCount[cell] < 65535) this.pixCount[cell]++;
      }
    }

    // Convert per-cell brightness to raw weights, normalized per-frame so gray-heavy frames still fill.
    // Then systematic-resample to exactly N particles.
    const BR_THR = DENSITY_THR;           // threshold in normalized brightness space (0..1)
    const GAMMA = DENSITY_GAMMA;          // >1 increases contrast (bright areas get more particles)
    const GREY_LIFT = DENSITY_GREY_LIFT;  // keep some midtone support
    const SMOOTH = DENSITY_SMOOTH;        // neighbor smoothing (too high reduces contrast)

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
      outCounts.fill(0); outCounts[((ROWS >> 1) * COLS + (COLS >> 1)) | 0] = N; return;
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
        if (cnt > 0) this.tmp[i] = (this.weights[i] * (1 - SMOOTH)) + (sum / cnt) * SMOOTH;
      }
    }
    rawSum = 0;
    for (let i = 0; i < CELLS; i++) { this.weights[i] = this.tmp[i]; rawSum += this.weights[i]; }

    if (rawSum <= 1e-9) {
      outCounts.fill(0); outCounts[((ROWS >> 1) * COLS + (COLS >> 1)) | 0] = N; return;
    }

    // Systematic resampling: assigns exactly N particles with low variance and no scan-order bias
    outCounts.fill(0);
    let total = 0;
    for (let i = 0; i < CELLS; i++) { total += this.weights[i]; this.cdf[i] = total; }
    const step = total / N;
    // Deterministic phase so results are stable across different async load orders (local vs GitHub Pages).
    // `seed01` should be in [0,1). Fallback to random if not provided.
    let u = (typeof seed01 === "number" && seed01 >= 0) ? (seed01 % 1) * step : random(step);
    let i = 0;
    for (let k = 0; k < N; k++) {
      const th = u + k * step;
      while (i < CELLS - 1 && this.cdf[i] < th) i++;
      outCounts[i]++;
    }

    // Diagnostic: stable hash of the resulting per-cell counts.
    // Hash a stride so it's fast but still comparable.
    const countsHash = this._fnv1aU16(outCounts, max(1, floor(CELLS / 512)));

    return { pixHash, countsHash };
  },
  _path(a, idx) {
    // Keep filenames consistent with existing export: actX_00000.png ...
    return `assets/act${a}/act${a}_${nf(idx, 5)}.png`;
  },
  queueLen() { return this.queue.length; },
  ready(a) {
    const c = this.cache[a];
    if (!c) return 0;
    let n = 0;
    for (let i = 0; i < c.done.length; i++) n += c.done[i] ? 1 : 0;
    return n;
  },
  _pump() {
    while (this.inFlight < this.maxInFlight && this.queue.length > 0) {
      const job = this.queue.shift();
      const a = job.a, idx = job.idx, cycle = job.cycle;
      const c = this.cache[a];
      if (!c || c.cycle !== cycle) continue;
      if (idx < 0 || idx >= c.cycle) continue;
      if (c.done[idx] || c.failed[idx]) { c.loading[idx] = 0; continue; }

      this.inFlight++;
      const path = this._path(a, idx);
      loadImage(
        path,
        (img) => {
          try {
            const cc = this.cache[a];
            if (!cc || cc.cycle !== cycle) return;
            const off = idx * CELLS;
            const seed01 = hash01(((a + 1) * 73856093) ^ ((idx + 1) * 19349663) ^ (cycle * 83492791));
            const out = cc.counts.subarray(off, off + CELLS);
            const diag = this._compute(img, out, seed01);
            // Persist diagnostics for cross-host comparisons.
            cc.pixHash[idx] = (diag && diag.pixHash ? diag.pixHash : 0) >>> 0;
            cc.countsHash[idx] = (diag && diag.countsHash ? diag.countsHash : 0) >>> 0;
            cc.imgW[idx] = (img.width | 0) & 0xffff;
            cc.imgH[idx] = (img.height | 0) & 0xffff;
            cc.done[idx] = 1;
          } catch (e) {
            try {
              const cc = this.cache[a];
              if (cc && cc.cycle === cycle) cc.failed[idx] = 1;
              console.warn("[Sampler] compute failed", a, idx, e);
            } catch (_) {}
          } finally {
            try {
              const cc = this.cache[a];
              if (cc && cc.cycle === cycle) cc.loading[idx] = 0;
            } catch (_) {}
            this.inFlight--;
            this._pump();
          }
        },
        (err) => {
          try {
            const cc = this.cache[a];
            if (cc && cc.cycle === cycle) {
              cc.failed[idx] = 1;
              cc.loading[idx] = 0;
            }
            console.warn("[Sampler] loadImage failed", path, err);
          } catch (_) {}
          this.inFlight--;
          this._pump();
        }
      );
    }
  },
  ensure(a, idx, cycle, outOff) {
    this.ensureActCache(a, cycle);
    const c = this.cache[a]; if (!c || idx < 0 || idx >= c.cycle) return false;
    const off = idx * CELLS;
    outOff.off = off;
    if (c.done[idx]) { this.hits++; this.hitsF++; return true; }
    this.misses++; this.missesF++;
    if (c.failed[idx]) return false;
    if (!c.loading[idx]) {
      c.loading[idx] = 1;
      this.queue.push({ a, idx, cycle });
      this._pump();
    }
    return false;
  },
  counts(a) { return this.cache[a] ? this.cache[a].counts : null; },
  framePixHash(a, idx) {
    const c = this.cache[a];
    if (!c || !c.pixHash || idx < 0 || idx >= c.cycle) return 0;
    return c.pixHash[idx] >>> 0;
  },
  frameCountsHash(a, idx) {
    const c = this.cache[a];
    if (!c || !c.countsHash || idx < 0 || idx >= c.cycle) return 0;
    return c.countsHash[idx] >>> 0;
  },
  frameImgWH(a, idx) {
    const c = this.cache[a];
    if (!c || !c.imgW || !c.imgH || idx < 0 || idx >= c.cycle) return { w: 0, h: 0 };
    return { w: c.imgW[idx] | 0, h: c.imgH[idx] | 0 };
  },
};

// ---- Particles: grid-density with local rebalancing + cheap separation ----
