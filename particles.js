const Particles = {
  x: new Float32Array(N), y: new Float32Array(N), vx: new Float32Array(N), vy: new Float32Array(N),
  kind: new Uint8Array(N), cell: new Uint16Array(N), u: new Float32Array(N), v: new Float32Array(N),
  order: new Uint16Array(N),
  // Transition dance roam targets (full-canvas flight)
  roamX: new Float32Array(N),
  roamY: new Float32Array(N),
  cellCounts: new Uint16Array(1), desired: new Uint16Array(1),
  // walkable region = union of src0 and src1 (prevents routing through holes that are empty in both frames)
  walk: new Uint8Array(1),
  // transport lanes (for clean motion across empty space)
  trA: new Uint8Array(N),
  trP: new Float32Array(N),
  trSign: new Int8Array(N),
  trOff: new Float32Array(N),
  trSX: new Float32Array(N), trSY: new Float32Array(N),
  trEX: new Float32Array(N), trEY: new Float32Array(N),
  trNX: new Float32Array(N), trNY: new Float32Array(N),
  trAmp: new Float32Array(N),
  // spatial hash
  binSize: 6, binsX: 1, binsY: 1, binHead: new Int32Array(1), binNext: new Int32Array(N),
  // coarse hash for "bird flock" dance (larger bins => cheap neighbor queries)
  fBinSize: 24, fBinsX: 1, fBinsY: 1, fHead: new Int32Array(1), fNext: new Int32Array(N),

  realloc() {
    this.cellCounts = new Uint16Array(CELLS);
    this.desired = new Uint16Array(CELLS);
    this.walk = new Uint8Array(CELLS);
  },
  resizeBins() {
    this.binsX = max(1, ceil(width / this.binSize));
    this.binsY = max(1, ceil(height / this.binSize));
    this.binHead = new Int32Array(this.binsX * this.binsY);

    this.fBinsX = max(1, ceil(width / this.fBinSize));
    this.fBinsY = max(1, ceil(height / this.fBinSize));
    this.fHead = new Int32Array(this.fBinsX * this.fBinsY);
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
      this.trP[i] = 0;
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
  setWalkFromCounts(counts, off0, off1) {
    // Union mask: any cell with particles in either frame is considered walkable.
    // This is cheaper than using image alpha directly and matches our density-field intent.
    for (let i = 0; i < CELLS; i++) {
      this.walk[i] = (counts[off0 + i] | 0) > 0 || (counts[off1 + i] | 0) > 0 ? 1 : 0;
    }
  },
  _bestNeighborDeficit(cell, pid) {
    const r = (cell / COLS) | 0, c = cell - r * COLS;
    let best = -1, bestNeed = -1e9;
    const jitter = hash01(((pid + 1) * 1103515245) ^ (((t / 0.08) | 0) * 12345));

    // Ring 1 (8-neighborhood)
    for (let dr = -1; dr <= 1; dr++) {
      const rr = r + dr; if (rr < 0 || rr >= ROWS) continue;
      for (let dc = -1; dc <= 1; dc++) {
        if (dr === 0 && dc === 0) continue;
        const cc = c + dc; if (cc < 0 || cc >= COLS) continue;
        const nb = rr * COLS + cc;
        const need = (this.desired[nb] | 0) - (this.cellCounts[nb] | 0);
        if (need > bestNeed || (need === bestNeed && need > 0 && hash01(nb ^ (pid * 2654435761)) < jitter)) {
          bestNeed = need;
          best = nb;
        }
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
          if (need > bestNeed || (need === bestNeed && need > 0 && hash01((nb + 7) ^ (pid * 1597334677)) < jitter)) {
            bestNeed = need;
            best = nb;
          }
        }
      }
    }

    // If no deficit nearby, take a guided 1-cell hop toward the strongest deficit in a wider radius.
    // This fixes "left behind" islands without teleporting across the screen.
    if (bestNeed <= 0) {
      if (_hotCount <= 0) return -1;
      // IMPORTANT: don't always pick the same global best hotspot (creates traveling "lumps").
      // Sample a few candidates deterministically per particle so flow spreads calmly.
      let target = -1;
      let bestScore = 0;
      const maxK = min(_hotCount, 24);
      const h0 = ((pid + 1) * 374761393) ^ (((t / 0.08) | 0) * 668265263);
      const tries = min(6, maxK);
      for (let s = 0; s < tries; s++) {
        const k = (h0 + s * 1640531527) % maxK;
        const kk = (k < 0 ? k + maxK : k) | 0;
        const cell2 = _hotCells[kk];
        const need2 = _hotNeeds[kk];
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
      // Prefer to stay within the union silhouette when possible.
      if ((this.walk[nb] | 0) === 0) {
        // Try a perpendicular alternate step that still progresses.
        const altR = nr !== r ? r : nr;
        const altC = nc !== c ? c : nc;
        // Two alternates: swap one axis
        const nb1 = altR * COLS + nc;
        const nb2 = nr * COLS + altC;
        const ok1 = nb1 !== cell && (this.walk[nb1] | 0) === 1;
        const ok2 = nb2 !== cell && (this.walk[nb2] | 0) === 1;
        if (ok1) return nb1;
        if (ok2) return nb2;
      }
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
        const to = this._bestNeighborDeficit(from, i); if (to < 0) continue;
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
  softRetargetToDesired() {
    // Same deterministic reassignment as hardRetargetToDesired(), but WITHOUT teleporting positions.
    // This keeps motion continuous at act boundaries.
    this.cellCounts.set(this.desired);
    let cell = 0, cum = this.desired[0] | 0;
    for (let oi = 0; oi < N; oi++) {
      const i = this.order[oi], target = oi + 0.5;
      while (cum <= target && cell < CELLS - 1) { cell++; cum += this.desired[cell] | 0; }
      this.cell[i] = cell;
      const h = (i + 1) * 73856093 ^ (cell + 1) * 19349663;
      this.u[i] = hash01(h); this.v[i] = hash01(h ^ 0x68bc21eb);
      this.trA[i] = 0;
      this.vx[i] *= 0.35;
      this.vy[i] *= 0.35;
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
  _buildFlockBins() {
    this.fHead.fill(-1);
    for (let i = 0; i < N; i++) {
      const bx = constrain((this.x[i] / this.fBinSize) | 0, 0, this.fBinsX - 1);
      const by = constrain((this.y[i] / this.fBinSize) | 0, 0, this.fBinsY - 1);
      const b = by * this.fBinsX + bx;
      this.fNext[i] = this.fHead[b]; this.fHead[b] = i;
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
    const dtN = constrain(tDelta / LANE_DT_REF, 0, 2);

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
      const posEmpty = (this.walk[posCell] | 0) === 0;
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
          if ((this.walk[sc] | 0) === 0) emptyHits++;
        }
        useLane = emptyHits >= VOID_EMPTY_NEED;
      }

      // If lane is engaged, route via 1 of 3 deterministic trajectories: straight / arc+ / arc-
      if (useLane) {
        const lane = i % 3; // 0=straight, 1=arc+, 2=arc-
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

          // Many thin "lanes" per trajectory: stable per-particle offset so we don't clump into 3 visible bands.
          const h = hash01((i + 1) * 2654435761);
          const off01 = (h - 0.5) * 2; // -1..1
          this.trSign[i] = sign;
          this.trOff[i] = off01 * cellW * (sign === 0 ? 0.35 : 0.95);
          this.trAmp[i] = sign === 0 ? 0 : min(cellW * 7.5, len * 0.22) * (0.80 + 0.25 * abs(off01));

          // Initialize progress from projection (prevents snapping when lane starts)
          const len2 = dx * dx + dy * dy + 1e-6;
          let p0 = ((this.x[i] - this.trSX[i]) * dx + (this.y[i] - this.trSY[i]) * dy) / len2;
          this.trP[i] = constrain(p0, 0, 1);
        } else {
          // Follow destination smoothly so the lane doesn't aim at stale positions.
          this.trEX[i] = lerp(this.trEX[i], baseTx, 0.18);
          this.trEY[i] = lerp(this.trEY[i], baseTy, 0.18);
          const dx = this.trEX[i] - this.trSX[i];
          const dy = this.trEY[i] - this.trSY[i];
          const len = sqrt(dx * dx + dy * dy) + 1e-6;
          this.trNX[i] = -dy / len;
          this.trNY[i] = dx / len;
        }

        const sx = this.trSX[i], sy = this.trSY[i];
        const ex = this.trEX[i], ey = this.trEY[i];
        const dx = ex - sx, dy = ey - sy;
        const len = sqrt(dx * dx + dy * dy) + 1e-6;

        // Explicit progress (fixes stalling): advance progress based on a px-per-sound-tick speed.
        // This keeps trajectories clean and prevents particles getting "stuck" on curved arcs.
        const stepPx = (cellW * (0.55 + cu * 1.10)) * dtN;
        this.trP[i] = min(1, this.trP[i] + stepPx / len + 0.002 * dtN);
        const prog = this.trP[i];
        const prog2 = min(1, prog + max(LANE_LOOKAHEAD, (stepPx / len) * 1.35));

        const sgn = this.trSign[i] | 0;
        const sin1 = sgn === 0 ? 0 : sin(prog * PI);
        const sin2 = sgn === 0 ? 0 : sin(prog2 * PI);

        const amp = this.trAmp[i] * sgn;
        const nx = this.trNX[i], ny = this.trNY[i];
        const off = this.trOff[i];

        const px = sx + dx * prog + nx * (off + amp * sin1);
        const py = sy + dy * prog + ny * (off + amp * sin1);
        const ax = sx + dx * prog2 + nx * (off + amp * sin2);
        const ay = sy + dy * prog2 + ny * (off + amp * sin2);

        // Pull toward lane and push forward along it (orderly trajectories)
        const pull = 0.032 + cu * 0.010;
        const push = 0.42 + cu * 0.34;
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
        if (prog >= 0.999 || !posEmpty) this.trA[i] = 0;
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
    // Shooting-star flock: a single directed formation (pack) that travels with a clear heading.
    // We build an anisotropic "comet" formation around a moving head point.
    const cx = width * 0.5, cy = height * 0.5;
    const minD = min(width, height);

    // Keep the head roaming far for most of Stage A, then start homing toward the next focus late.
    const headBlend = map(s, 0.68, 1.0, 0, 0.90, true);
    const baseT = t * (0.22 + 0.06 * s);
    const hx0 = cx + sin(baseT * 0.9) * width * 0.46 + sin(baseT * 0.31 + 2.2) * width * 0.28;
    const hy0 = cy + sin(baseT * 0.7 + 1.1) * height * 0.40 + sin(baseT * 0.27 - 0.7) * height * 0.26;
    const hx = lerp(hx0, Render.fx || cx, headBlend);
    const hy = lerp(hy0, Render.fy || cy, headBlend);

    const baseT2 = (t + 0.03) * (0.22 + 0.06 * s);
    const hx02 = cx + sin(baseT2 * 0.9) * width * 0.46 + sin(baseT2 * 0.31 + 2.2) * width * 0.28;
    const hy02 = cy + sin(baseT2 * 0.7 + 1.1) * height * 0.40 + sin(baseT2 * 0.27 - 0.7) * height * 0.26;
    const hx2 = lerp(hx02, Render.fx || cx, headBlend);
    const hy2 = lerp(hy02, Render.fy || cy, headBlend);

    let dirX = hx2 - hx;
    let dirY = hy2 - hy;
    const dirL = sqrt(dirX * dirX + dirY * dirY) + 1e-6;
    dirX /= dirL; dirY /= dirL;
    const perpX = -dirY;
    const perpY = dirX;

    // Formation geometry
    const tailLen = minD * 0.78;
    const tailW = minD * 0.10;

    // Flight dynamics
    const maxSpd = 14.5 + map(level, 0, 0.2, 0, 4.0, true);
    const damp = 0.93;
    const steer = 0.34 + s * 0.20;
    const bankAmp = 0.55; // radians

    // Early break-away from last silhouette region
    const fx = Render.fx || cx;
    const fy = Render.fy || cy;
    const repR = minD * lerp(0.46, 0.22, s);
    const repR2 = repR * repR;
    const repW = (1 - s) * 1.35;

    const margin = minD * 0.06;
    const edgeW = 0.30;

    strokeCap(ROUND); // keep points round during flight

    for (let i = 0; i < N; i++) {
      const x = this.x[i], y = this.y[i];
      const vx = this.vx[i], vy = this.vy[i];

      // Stable slot along the comet tail (bias toward the head so it reads as a pack)
      const u = hash01((i + 1) * 2654435761);
      const tpos = u * u; // more near head
      const long = -tpos * tailLen;
      const w = tailW * (0.25 + 0.75 * tpos);
      const latBase = (hash01((i + 17) * 2246822519) - 0.5) * 2.0;
      const lat = latBase * w;

      // Bank/curve the whole formation a bit per-particle (bird-like arcs)
      const ca = (noise(i * 0.019, t * 0.18) - 0.5) * 2.0 * bankAmp;
      const cc = cos(ca), ss = sin(ca);
      const bDirX = dirX * cc - dirY * ss;
      const bDirY = dirX * ss + dirY * cc;
      const bPerpX = -bDirY;
      const bPerpY = bDirX;

      // Desired position in the moving formation
      const wig = (1 - tpos) * 10.0;
      const wx = (noise(i * 0.013, 9.1, t * 0.42) - 0.5) * wig;
      const wy = (noise(i * 0.017, 3.7, t * 0.42) - 0.5) * wig;
      const tx = hx + bDirX * long + bPerpX * lat + wx;
      const ty = hy + bDirY * long + bPerpY * lat + wy;

      // Seek steering toward formation slot
      let ax = (tx - x) * 0.010;
      let ay = (ty - y) * 0.010;

      // Early repulsion from old silhouette area
      if (repW > 0) {
        const dxF = x - fx, dyF = y - fy;
        const dF2 = dxF * dxF + dyF * dyF;
        if (dF2 < repR2) {
          const dF = sqrt(dF2) + 1e-6;
          const k = (1 - dF / repR) * repW * 0.9;
          ax += (dxF / dF) * k;
          ay += (dyF / dF) * k;
        }
      }

      // Boundary steering
      if (x < margin) ax += ((margin - x) / margin) * edgeW;
      else if (x > width - margin) ax -= ((x - (width - margin)) / margin) * edgeW;
      if (y < margin) ay += ((margin - y) / margin) * edgeW;
      else if (y > height - margin) ay -= ((y - (height - margin)) / margin) * edgeW;

      // Velocity-based steering (keeps smooth, controlled flight)
      let nvx = vx + ax * steer;
      let nvy = vy + ay * steer;
      nvx *= damp; nvy *= damp;

      const sp2 = nvx * nvx + nvy * nvy;
      const ms2 = maxSpd * maxSpd;
      if (sp2 > ms2) {
        const k = maxSpd / sqrt(sp2);
        nvx *= k; nvy *= k;
      }

      this.vx[i] = nvx;
      this.vy[i] = nvy;
      this.x[i] = x + nvx;
      this.y[i] = y + nvy;
    }
  },
  draw() {
    noStroke();
    // Base: always white lights (additive)
    blendMode(ADD);
    strokeCap(ROUND);
    strokeWeight(PARTICLE_SIZE);
    stroke(255, 255, 255, LIGHT_ALPHA);
    for (let i = 0; i < N; i++) point(this.x[i], this.y[i]);

    // Color overlay is *localized* (grouped) around a focus point.
    // Use normal alpha blending so tint shows over white (ADD would clamp to white).
    if (Render.tintA > 0.5) {
      blendMode(BLEND);
      const fx = Render.fx, fy = Render.fy;
      const outer2 = Render.outerR2;
      const sx = Render.stainX, sy = Render.stainY;
      const stainR = Render.stainR;
      const maxD = max(1, Render.maxD | 0);
      const edgeFreq = 0.010;

      // Majority tint: one calm "wash" over the focus region (no contouring by density)
      const aWash = Render.tintA | 0;
      stroke(red(Render.hi), green(Render.hi), blue(Render.hi), aWash);
      for (let i = 0; i < N; i++) {
        const dx = this.x[i] - fx, dy = this.y[i] - fy;
        const d2 = dx * dx + dy * dy;
        if (d2 <= outer2) point(this.x[i], this.y[i]);
      }

      // Minority tint: a single "stain" blob sitting on top of the wash (irregular edge via noise).
      const aStain = Render.tintB | 0;
      if (aStain > 0) {
        stroke(red(Render.lo), green(Render.lo), blue(Render.lo), aStain);
        const stainMax2 = (stainR * 1.35) * (stainR * 1.35);
        for (let i = 0; i < N; i++) {
          const dx0 = this.x[i] - fx, dy0 = this.y[i] - fy;
          const d0 = dx0 * dx0 + dy0 * dy0;
          if (d0 > outer2) continue;

          const dx = this.x[i] - sx, dy = this.y[i] - sy;
          const d2 = dx * dx + dy * dy;
          if (d2 > stainMax2) continue;

          // Rough edge: low-frequency noise in world space
          const n = noise(this.x[i] * edgeFreq, this.y[i] * edgeFreq, (Acts.act || 1) * 10.0);

          // Density influence (subtle): lower-density areas accept the stain slightly more
          const den = (this.desired[this.cell[i]] | 0) / maxD; // 0..1
          const denK = lerp(1.15, 0.85, den);

          const edge = stainR * denK * (0.82 + 0.30 * n);
          if (d2 <= edge * edge) point(this.x[i], this.y[i]);
        }
      }
    }

    blendMode(BLEND);
  },
};

Particles.onTransitionStart = function () {
  // Burst + de-correlate velocities so the flock immediately leaves the old silhouette.
  const cx = width * 0.5, cy = height * 0.5;
  const kick = min(width, height) * 0.0052;
  const margin = min(width, height) * 0.06;
  const minFar2 = (min(width, height) * 0.30) ** 2;
  for (let i = 0; i < N; i++) {
    const x = this.x[i], y = this.y[i];
    const dx = x - cx, dy = y - cy;
    const d = sqrt(dx * dx + dy * dy) + 1e-6;
    const h = hash01((i + 1) * 2654435761);
    const ang = TWO_PI * h;
    // Mostly outward + a bit of sideways "flap"
    this.vx[i] = this.vx[i] * 0.35 + (dx / d) * kick + cos(ang) * kick * 0.55;
    this.vy[i] = this.vy[i] * 0.35 + (dy / d) * kick + sin(ang) * kick * 0.55;
    this.trA[i] = 0;

    // Pick a far roam target so the flock travels the whole space.
    let rx = 0, ry = 0;
    for (let k = 0; k < 6; k++) {
      rx = random(margin, width - margin);
      ry = random(margin, height - margin);
      const ddx = rx - x, ddy = ry - y;
      if (ddx * ddx + ddy * ddy >= minFar2) break;
    }
    this.roamX[i] = rx;
    this.roamY[i] = ry;

    // Nudge initial velocity toward the roam target so the flock immediately traverses space.
    const rdx = rx - x, rdy = ry - y;
    const rl = sqrt(rdx * rdx + rdy * rdy) + 1e-6;
    this.vx[i] += (rdx / rl) * kick * 0.65;
    this.vy[i] += (rdy / rl) * kick * 0.65;
  }
};
