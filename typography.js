// Typography overlay rendered as particles (template-driven).
// Uses two separate templates:
// - assets/typography/bigtypography.png
// - assets/typography/smalltypography.png

const Typography = {
  bigImg: null,
  smallImg: null,

  bigRect: { x: 0, y: 0, w: 0, h: 0 },
  smallRect: { x: 0, y: 0, w: 0, h: 0 },

  // Tunables
  TOTAL: 9000,
  BIG_N: 3200,
  SMALL_N: 5800,
  BIG_SIZE: 2.0,
  SMALL_SIZE: 0.95,
  BASE_ALPHA: 215,

  // Particles (in normalized template space, 0..1)
  n: 0,
  group: new Uint8Array(0), // 0=big, 1=small
  x: new Float32Array(0),
  y: new Float32Array(0),
  vx: new Float32Array(0),
  vy: new Float32Array(0),
  seed: new Float32Array(0),
  kind: new Uint8Array(0), // 0=dim, 1=hi
  anchor: new Uint32Array(0), // packed: (group<<24) | idx

  // Cached small-text twinkle + jitter (updated incrementally for performance).
  _tw: new Float32Array(0),
  _sx: new Float32Array(0),
  _sy: new Float32Array(0),
  _stampOx: new Float32Array(0),
  _stampOy: new Float32Array(0),
  _cursor: 0,
  _chunk: 1600,

  // Template data
  _big: null,
  _small: null,
  ready: false,
  _stats: { bigCand: 0, smallCand: 0 },

  preload() {
    this.bigImg = loadImage("assets/typography/bigtypography.png");
    this.smallImg = loadImage("assets/typography/smalltypography.png");
  },

  init() {
    this.ready = false;
    this._alloc();
    this.resize();
  },

  resize() {
    this.ready = false;
    // Huge source PNGs can blow up memory (especially on mobile) and make everything stutter.
    // We only need enough resolution to build low-res masks, so downscale aggressively.
    this._downscaleSources();
    this._updateRects();
    this._buildTemplates();
    this._seedParticles();
    this.ready = true;
  },

  debugString() {
    const bw = this.bigImg ? this.bigImg.width : 0;
    const bh = this.bigImg ? this.bigImg.height : 0;
    const sw = this.smallImg ? this.smallImg.width : 0;
    const sh = this.smallImg ? this.smallImg.height : 0;
    return `typ big:${bw}x${bh} small:${sw}x${sh} bigCand:${this._stats.bigCand} smallCand:${this._stats.smallCand}`;
  },

  // ---- internals ----
  _downscaleSources() {
    const maxDim = 2400;
    const down = (img) => {
      if (!img || !img.width) return;
      const w = img.width, h = img.height;
      const m = max(w, h);
      if (m <= maxDim) return;
      if (w >= h) img.resize(maxDim, 0);
      else img.resize(0, maxDim);
    };
    down(this.bigImg);
    down(this.smallImg);
  },
  _alloc() {
    const total = max(1, this.TOTAL | 0);
    const bigN = constrain(this.BIG_N | 0, 0, total);
    const smallN = constrain(this.SMALL_N | 0, 0, total - bigN);
    this.n = bigN + smallN;
    this.group = new Uint8Array(this.n);
    this.x = new Float32Array(this.n);
    this.y = new Float32Array(this.n);
    this.vx = new Float32Array(this.n);
    this.vy = new Float32Array(this.n);
    this.seed = new Float32Array(this.n);
    this.kind = new Uint8Array(this.n);
    this.anchor = new Uint32Array(this.n);

    this._tw = new Float32Array(this.n);
    this._sx = new Float32Array(this.n);
    this._sy = new Float32Array(this.n);
    this._stampOx = new Float32Array(this.n);
    this._stampOy = new Float32Array(this.n);
    this._cursor = 0;

    for (let i = 0; i < this.n; i++) this.group[i] = i < bigN ? 0 : 1;
  },

  _updateRects() {
    if (this.bigImg && this.bigImg.width) {
      fitRect(this.bigImg.width, this.bigImg.height, width, height, this.bigRect);
    } else {
      this.bigRect.x = 0; this.bigRect.y = 0; this.bigRect.w = width; this.bigRect.h = height;
    }

    // Small text is authored as a full-screen layout; map it directly to the canvas.
    // (This stays stable even if you change the source PNG dimensions.)
    this.smallRect.x = 0; this.smallRect.y = 0; this.smallRect.w = width; this.smallRect.h = height;
  },

  _buildTemplates() {
    // BIG: preserve template aspect for crisp big letters.
    this._big = this._buildTemplate(this.bigImg, {
      w: 820,
      h: 0,
      drawMode: "FIT",
      alphaThr: 48,
      step: 2,
      gamma: 1.15,
    });

    // SMALL: build mask in canvas aspect so it maps 1:1 to screen space.
    this._small = this._buildTemplate(this.smallImg, {
      w: 960,
      h: max(1, floor(960 * (height / max(1, width)))),
      drawMode: "STRETCH",
      // Small text is thin; include anti-aliased edges but avoid band leakage.
      alphaThr: 14,
      step: 1,
      gamma: 1.0,
      noTangent: true,
      dilate: 0,
    });

    this._stats.bigCand = this._big ? this._big.n : 0;
    this._stats.smallCand = this._small ? this._small.n : 0;
  },

  _buildTemplate(img, opt) {
    if (!img || !img.width) return null;
    const W = max(64, opt.w | 0);
    const H = opt.h ? max(64, opt.h | 0) : max(64, floor(W * (img.height / img.width)));
    const step = max(1, opt.step | 0);
    const thr = constrain(opt.alphaThr == null ? 48 : opt.alphaThr | 0, 1, 254);
    const gamma = opt.gamma == null ? 1.15 : opt.gamma;

    const g = createGraphics(W, H);
    g.pixelDensity(1);
    g.noSmooth();
    // We read pixels from this buffer to build a mask; optimize for frequent readback.
    try {
      const ctx = g.canvas && g.canvas.getContext ? g.canvas.getContext("2d", { willReadFrequently: true }) : null;
      if (ctx) {
        g.drawingContext = ctx;
        if (g._renderer) g._renderer.drawingContext = ctx;
        try {
          ctx.imageSmoothingEnabled = false;
          if (ctx.setTransform) ctx.setTransform(1, 0, 0, 1, 0, 0);
          ctx.globalCompositeOperation = "source-over";
          ctx.globalAlpha = 1;
        } catch (_) {}
      }
    } catch (_) {}
    g.clear();
    g.imageMode(CORNER);
    if (opt.drawMode === "STRETCH") {
      g.image(img, 0, 0, W, H);
    } else {
      // FIT into buffer (preserve aspect)
      const r = { x: 0, y: 0, w: 0, h: 0 };
      fitRect(img.width, img.height, W, H, r);
      // Quantize to avoid subpixel differences between hosts/GPUs.
      r.x = round(r.x);
      r.y = round(r.y);
      r.w = max(1, round(r.w));
      r.h = max(1, round(r.h));
      if (r.x < 0) r.x = 0;
      if (r.y < 0) r.y = 0;
      if (r.x + r.w > W) r.w = max(1, W - r.x);
      if (r.y + r.h > H) r.h = max(1, H - r.y);
      g.image(img, r.x, r.y, r.w, r.h);
    }
    g.loadPixels();

    const size = W * H;
    let mask = new Uint8Array(size);
    let maxA = 0;
    for (let i = 0; i < size; i++) {
      const a = g.pixels[4 * i + 3] | 0;
      mask[i] = a;
      if (a > maxA) maxA = a;
    }

    // Fallback if alpha isn't useful: treat dark pixels as ink on a solid background.
    let useAlpha = maxA > 8;
    if (!useAlpha) {
      for (let i = 0; i < size; i++) {
        const p = 4 * i;
        const br = (g.pixels[p] + g.pixels[p + 1] + g.pixels[p + 2]) / 765;
        const a = br < 0.70 ? 255 : 0;
        mask[i] = a;
      }
      useAlpha = true;
    }

    // Optional dilation (helps thin glyphs stay readable at low resolution).
    const dil = opt.dilate | 0;
    if (dil > 0) {
      for (let pass = 0; pass < dil; pass++) {
        const out = mask.slice();
        for (let y = 1; y < H - 1; y++) {
          const row = y * W;
          for (let x = 1; x < W - 1; x++) {
            const i = row + x;
            const a = mask[i] | 0;
            if (a < thr) continue;
            if (a > out[i - 1]) out[i - 1] = a;
            if (a > out[i + 1]) out[i + 1] = a;
            if (a > out[i - W]) out[i - W] = a;
            if (a > out[i + W]) out[i + W] = a;
          }
        }
        mask = out;
      }
    }

    const cx = [];
    const cy = [];
    const tx = [];
    const ty = [];
    const wgt = [];

    for (let y = 1; y < H - 1; y += step) {
      const row = y * W;
      for (let x = 1; x < W - 1; x += step) {
        const i = row + x;
        const a = mask[i] | 0;
        if (a < thr) continue;
        const a01 = a / 255;
        const ww = pow(a01, gamma);
        cx.push((x + 0.5) / W);
        cy.push((y + 0.5) / H);
        wgt.push(ww);

        if (opt.noTangent) {
          tx.push(0);
          ty.push(0);
        } else {
          const axp = mask[i + 1] | 0, axm = mask[i - 1] | 0;
          const ayp = mask[i + W] | 0, aym = mask[i - W] | 0;
          const gx = axp - axm;
          const gy = ayp - aym;
          let txx = -gy, tyy = gx;
          const l = sqrt(txx * txx + tyy * tyy) + 1e-9;
          txx /= l; tyy /= l;
          tx.push(txx);
          ty.push(tyy);
        }
      }
    }

    const n = cx.length | 0;
    if (n <= 0) return { W, H, mask, thr, n: 0 };

    const cX = new Float32Array(n);
    const cY = new Float32Array(n);
    const cTX = new Float32Array(n);
    const cTY = new Float32Array(n);
    const cW = new Float32Array(n);
    const cCDF = new Float32Array(n);
    let total = 0;
    for (let i = 0; i < n; i++) {
      cX[i] = cx[i];
      cY[i] = cy[i];
      cTX[i] = tx[i];
      cTY[i] = ty[i];
      const ww = wgt[i];
      cW[i] = ww;
      total += ww;
      cCDF[i] = total;
    }

    return { W, H, mask, thr, n, cX, cY, cTX, cTY, cW, cCDF, total };
  },

  _pick(template, u) {
    if (!template || template.n <= 0 || template.total <= 1e-9) return 0;
    const r = u * template.total;
    const cdf = template.cCDF;
    let lo = 0, hi = template.n - 1;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (cdf[mid] >= r) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  },

  _maskAt(template, x01, y01) {
    if (!template) return 0;
    const x = constrain((x01 * template.W) | 0, 0, template.W - 1);
    const y = constrain((y01 * template.H) | 0, 0, template.H - 1);
    return template.mask[y * template.W + x] | 0;
  },

  _seedParticles() {
    if (!this._big || !this._small) return;
    const bigN = min(this.n, this.BIG_N | 0);
    const smallN = this.n - bigN;

    // Stratified sampling so points spread nicely across the mask.
    const seedGroup = (gid, start, count, tpl) => {
      if (!tpl || tpl.n <= 0 || tpl.total <= 1e-9 || count <= 0) return;
      const stepC = tpl.total / count;
      let u0 = random(stepC);
      for (let i = 0; i < count; i++) {
        const pid = start + i;
        const u = (u0 + i * stepC) / tpl.total;
        const idx = this._pick(tpl, u);
        this.anchor[pid] = ((gid & 0xff) << 24) | (idx & 0xffffff);
        const ax = tpl.cX[idx];
        const ay = tpl.cY[idx];
        const jjx = (hash01((pid + 1) * 2654435761) - 0.5) * 0.0018;
        const jjy = (hash01((pid + 1) * 1597334677) - 0.5) * 0.0018;
        this.x[pid] = constrain(ax + jjx, 0, 1);
        this.y[pid] = constrain(ay + jjy, 0, 1);
        this.vx[pid] = (hash01((pid + 1) * 2246822519) - 0.5) * 0.0011;
        this.vy[pid] = (hash01((pid + 1) * 3266489917) - 0.5) * 0.0011;
        this.seed[pid] = hash01((pid + 1) * 1103515245) * 1000;
        this.kind[pid] = hash01((pid + 1) * 374761393) < (gid === 0 ? 0.70 : 0.62) ? 0 : 1;
      }
    };

    seedGroup(0, 0, bigN, this._big);
    seedGroup(1, bigN, smallN, this._small);

    // Initialize cached twinkle/jitter and 2nd-stamp offsets (avoids per-frame noise/hash for small text).
    for (let i = 0; i < this.n; i++) {
      this._tw[i] = 1;
      this._sx[i] = 0;
      this._sy[i] = 0;
      this._stampOx[i] = (hash01((i + 1) * 2654435761) - 0.5) * 0.42;
      this._stampOy[i] = (hash01((i + 1) * 1597334677) - 0.5) * 0.42;
    }
  },

  update(level) {
    if (!this.ready || !tAdvanced) return;
    const n = this.n | 0;
    if (n <= 0) return;

    const bigN = min(n, this.BIG_N | 0);
    const smallStart = bigN;
    const smallN = n - smallStart;

    const dt = min(0.08, max(0.004, abs(tDelta)));
    const kdt = constrain(dt / 0.03, 0.75, 1.8);

    const damp = 0.90;
    const springBig = (0.052 + level * 0.070) * kdt;
    const flowBig = (0.0014 + level * 0.0042) * kdt;
    const alongBig = (0.0015 + level * 0.0042) * kdt;
    const maxBig = (0.008 + level * 0.014) * kdt;

    // Small text: keep motion micro to preserve readability.
    const wobSmall = (0.0012 + level * 0.0018) * kdt; // ~sub-pixel to ~1px in screen space

    // BIG particles: full update (heavier motion and texture).
    for (let i = 0; i < bigN; i++) {
      const packed = this.anchor[i] >>> 0;
      const idx = packed & 0xffffff;
      const tpl = this._big;
      if (!tpl || idx >= tpl.n) continue;

      const ax = tpl.cX[idx];
      const ay = tpl.cY[idx];
      const txx = tpl.cTX[idx];
      const tyy = tpl.cTY[idx];

      let x = this.x[i], y = this.y[i];
      let vx = this.vx[i], vy = this.vy[i];

      const spring = springBig;
      const flow = flowBig;
      const along = alongBig;
      const maxSpd = maxBig;

      // Calm flow around anchor (keeps motion "in" the typography).
      const n0 = noise(ax * 4.2, ay * 4.2, t * 0.16 + this.seed[i]);
      const ang = n0 * TWO_PI * 2.0;
      vx += cos(ang) * flow;
      vy += sin(ang) * flow;

      // Along-stroke drift for big letters only (adds life without smearing small text).
      const sgn = (hash01((i + 1) * 747796405) < 0.5 ? -1 : 1);
      vx += txx * along * sgn;
      vy += tyy * along * sgn;

      // Spring to anchor.
      vx += (ax - x) * spring;
      vy += (ay - y) * spring;

      vx *= damp;
      vy *= damp;

      const sp2 = vx * vx + vy * vy;
      const ms2 = maxSpd * maxSpd;
      if (sp2 > ms2) {
        const s = maxSpd / sqrt(sp2);
        vx *= s; vy *= s;
      }

      let nx = x + vx;
      let ny = y + vy;

      // Constrain to letter alpha.
      const ok = this._maskAt(tpl, nx, ny) >= tpl.thr;
      if (nx < 0 || nx > 1 || ny < 0 || ny > 1 || !ok) {
        nx = lerp(nx, ax, 0.58);
        ny = lerp(ny, ay, 0.58);
        vx *= 0.25;
        vy *= 0.25;
      }

      this.x[i] = nx;
      this.y[i] = ny;
      this.vx[i] = vx;
      this.vy[i] = vy;
    }

    // SMALL particles: incremental updates (keeps the look but avoids thousands of noise calls per frame).
    if (smallN > 0) {
      const tpl = this._small;
      const step = max(64, min(this._chunk | 0, smallN));
      let cur = this._cursor | 0;
      const tTw = t * 0.05;
      const tWb = t * 0.10;

      for (let k = 0; k < step; k++) {
        const i = smallStart + cur;
        const packed = this.anchor[i] >>> 0;
        const idx = packed & 0xffffff;
        if (!tpl || idx >= tpl.n) { cur++; if (cur >= smallN) cur = 0; continue; }

        const ax = tpl.cX[idx];
        const ay = tpl.cY[idx];

        // Cached twinkle (slow) + cached wobble offsets (micro).
        this._tw[i] = 0.62 + 0.38 * noise(this.seed[i], tTw);
        this._sx[i] = (noise(this.seed[i] * 0.0123, tWb) - 0.5) * 2.0 * wobSmall;
        this._sy[i] = (noise(this.seed[i] * 0.0191 + 7.7, tWb) - 0.5) * 2.0 * wobSmall;

        let nx = ax + this._sx[i];
        let ny = ay + this._sy[i];
        if (this._maskAt(tpl, nx, ny) < tpl.thr) { nx = ax; ny = ay; }
        this.x[i] = nx;
        this.y[i] = ny;
        this.vx[i] *= 0.25;
        this.vy[i] *= 0.25;

        cur++;
        if (cur >= smallN) cur = 0;
      }
      this._cursor = cur;
    }
  },

  draw(level) {
    if (!this.ready) return;
    if (typeof sceneA === "number" && sceneA <= 0.004) return;
    if (!this._big || !this._small) return;

    this.update(level);
    const aScene = typeof sceneA === "number" ? sceneA : 1;
    const aBaseDim = this.BASE_ALPHA * 0.85 * aScene;
    const aBaseHi = this.BASE_ALPHA * 1.25 * aScene;
    const aMoveDim = this.BASE_ALPHA * 0.40 * aScene;
    const aMoveHi = this.BASE_ALPHA * 0.60 * aScene;

    push();
    blendMode(ADD);
    strokeCap(ROUND);

    // Base pass for readability (anchored).
    const sizeScale = (typeof cellW === "number" && typeof CELL_SIZE === "number" && CELL_SIZE > 0)
      ? (cellW / CELL_SIZE)
      : 1;
    let swBig = max(1.0, round(this.BIG_SIZE * sizeScale));
    if ((swBig & 1) && swBig > 1) swBig--;
    let swSmall = max(1.0, round(this.SMALL_SIZE * sizeScale));
    if ((swSmall & 1) && swSmall > 1) swSmall--;
    strokeWeight(swBig);
    const bigN = min(this.n, this.BIG_N | 0);
    // Batch big anchors by kind to reduce state changes (same visual output).
    stroke(255, 255, 255, aBaseDim);
    for (let i = 0; i < bigN; i++) {
      if (this.kind[i]) continue;
      const packed = this.anchor[i] >>> 0;
      const idx = packed & 0xffffff;
      const tpl = this._big;
      if (idx >= tpl.n) continue;
      point(((this.bigRect.x + tpl.cX[idx] * this.bigRect.w + 0.5) | 0), ((this.bigRect.y + tpl.cY[idx] * this.bigRect.h + 0.5) | 0));
    }
    stroke(255, 255, 255, aBaseHi);
    for (let i = 0; i < bigN; i++) {
      if (!this.kind[i]) continue;
      const packed = this.anchor[i] >>> 0;
      const idx = packed & 0xffffff;
      const tpl = this._big;
      if (idx >= tpl.n) continue;
      point(((this.bigRect.x + tpl.cX[idx] * this.bigRect.w + 0.5) | 0), ((this.bigRect.y + tpl.cY[idx] * this.bigRect.h + 0.5) | 0));
    }

    strokeWeight(swSmall);
    const smallStart = bigN;
    for (let i = smallStart; i < this.n; i++) {
      const packed = this.anchor[i] >>> 0;
      const idx = packed & 0xffffff;
      const tpl = this._small;
      if (idx >= tpl.n) continue;

      // Cached twinkle keeps the look but avoids per-frame noise() on every particle.
      const tw = this._tw[i] || 1;
      stroke(255, 255, 255, (this.kind[i] ? aBaseHi : aBaseDim) * tw);
      const ax = this.smallRect.x + tpl.cX[idx] * this.smallRect.w;
      const ay = this.smallRect.y + tpl.cY[idx] * this.smallRect.h;
      point(((ax + 0.5) | 0), ((ay + 0.5) | 0));
      if ((i & 1) === 0) point((((ax + this._stampOx[i]) + 0.5) | 0), (((ay + this._stampOy[i]) + 0.5) | 0));
    }

    // Moving texture for big letters only (half-rate for perf).
    let swMove = max(1.0, round(max(1.2, this.BIG_SIZE * 0.72) * sizeScale));
    if ((swMove & 1) && swMove > 1) swMove--;
    strokeWeight(swMove);
    for (let i = 0; i < this.n; i++) {
      const packed = this.anchor[i] >>> 0;
      const gid = (packed >>> 24) & 0xff;
      if (gid !== 0) continue;
      if ((i & 1) === 1) continue;
      const tw = 0.5 + 0.5 * sin(t * 0.12 + this.seed[i]);
      stroke(255, 255, 255, (this.kind[i] ? aMoveHi : aMoveDim) * (0.75 + 0.25 * tw));
      point(((this.bigRect.x + this.x[i] * this.bigRect.w + 0.5) | 0), ((this.bigRect.y + this.y[i] * this.bigRect.h + 0.5) | 0));
    }

    blendMode(BLEND);
    pop();
  },
};
