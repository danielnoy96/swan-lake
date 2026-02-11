const Particles = {
  x: new Float32Array(N), y: new Float32Array(N), vx: new Float32Array(N), vy: new Float32Array(N),
  kind: new Uint8Array(N), cell: new Uint16Array(N), u: new Float32Array(N), v: new Float32Array(N),
  order: new Uint16Array(N),
  // Stable per-particle ribbon slots (computed once; no per-frame randomness)
  slotU: new Float32Array(N),     // 0..1 (0=head, 1=tail)
  slotLat: new Float32Array(N),   // -1..1 (perpendicular offset)
  slotPhase: new Float32Array(N), // noise seed
  // Transition dance roam targets (full-canvas flight)
  roamX: new Float32Array(N),
  roamY: new Float32Array(N),
  // Transition ribbon path (global, keeps the flock formed and directional)
  pathAx: 0, pathAy: 0,
  pathBx: 0, pathBy: 0,
  pathCx: 0, pathCy: 0,
  pathDx: 0, pathDy: 0,
  pathSeed: 0,
  // Transition flock leader (head) + 2-turn waypoint navigation
  headX: 0, headY: 0,
  headVX: 0, headVY: 0,
  wp1X: 0, wp1Y: 0,
  wp2X: 0, wp2Y: 0,
  leaderSeed: 0,
  // Lemniscate (∞) ribbon path params (reference path)
  loopCX0: 0, loopCY0: 0,
  loopCX1: 0, loopCY1: 0,
  loopScale: 1,
  loopRot: 0,
  loopTheta0: 0,
  loopThetaSpan: 0,
  transitionMode: 0, // 0=lemniscate ribbon, 1=circle (special case)
  circleR: 0,
  circleThick: 0,
  circleTheta0: 0,
  circleThetaSpan: 0,
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
      const h0 = hash01((i + 1) * 2654435761);
      const h1 = hash01((i + 17) * 2246822519);
      const h2 = hash01((i + 99) * 1597334677);
      this.slotU[i] = pow(h0, 2.10);
      this.slotLat[i] = (h1 - 0.5) * 2;
      this.slotPhase[i] = h2 * 1000;
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

    // Safe default ribbon path so any fallback `dance()` calls (e.g., while sampler warms up)
    // never collapse to (0,0) and scatter particles into the top-left corner.
    const minD = min(width, height);
    this.pathSeed = random(10000);
    this.pathAx = cx; this.pathAy = cy;
    this.pathDx = cx + minD * 0.38; this.pathDy = cy - minD * 0.06;
    this.pathBx = lerp(this.pathAx, this.pathDx, 0.25); this.pathBy = this.pathAy - minD * 0.20;
    this.pathCx = lerp(this.pathAx, this.pathDx, 0.75); this.pathCy = this.pathDy + minD * 0.22;

    // Safe default leader state (used by the new "bird pack" rules).
    this.leaderSeed = random(10000);
    this.headX = cx; this.headY = cy;
    this.headVX = minD * 0.10; this.headVY = -minD * 0.02;
    this.wp1X = constrain(cx + minD * 0.18, minD * 0.08, width - minD * 0.08);
    this.wp1Y = constrain(cy - minD * 0.22, minD * 0.08, height - minD * 0.08);
    this.wp2X = constrain(cx + minD * 0.34, minD * 0.08, width - minD * 0.08);
    this.wp2Y = constrain(cy + minD * 0.20, minD * 0.08, height - minD * 0.08);

    // Safe default ∞ loop (prevents fallback `dance()` from collapsing to a corner).
    this.loopCX0 = this.loopCX1 = cx;
    this.loopCY0 = this.loopCY1 = cy;
    this.loopScale = minD * 0.42;
    this.loopRot = 0;
    this.loopTheta0 = 0;
    this.loopThetaSpan = TWO_PI * 1.12;
    this.transitionMode = 0;
    this.circleR = minD * 0.22;
    this.circleThick = minD * 0.020;
    this.circleTheta0 = 0;
    this.circleThetaSpan = TWO_PI * 1.25;
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
    // Large tDelta (e.g., when AUTO_SPEED is high) makes lane progress jumpy.
    // Clamp so motion stays smooth even if time advances faster.
    const dtN = constrain(tDelta / LANE_DT_REF, 0, 1);

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
    if ((this.transitionMode | 0) === 1) return this._danceCircle(level, s);
    // Reference ribbon path (∞): a deterministic lemniscate that reads like a true ribbon.
    // Tail is less "saturated" (lower alpha), head is more saturated (higher alpha).
    const minD = min(width, height);
    const margin = minD * 0.06;
    const uT = constrain(s, 0, 1);

    // dt in "sound-time" units (silence => 0 => freeze)
    const dt = constrain(abs(tDelta), 0, 0.06);
    if (dt <= 0) return;
    const dtK = dt / 0.02;

    // Smooth progress along the loop
    const p = uT * uT * (3 - 2 * uT);

    // Move the whole loop across the canvas (gives travel / direction).
    // We also pull toward next focus near the end to avoid freeze-then-zip.
    const fx = Render.fx || width * 0.5;
    const fy = Render.fy || height * 0.5;
    const homeBlend = map(uT, 0.78, 1.0, 0, 1, true);

    let cx = lerp(this.loopCX0, this.loopCX1, p);
    let cy = lerp(this.loopCY0, this.loopCY1, p);
    cx = lerp(cx, fx, homeBlend * 0.55);
    cy = lerp(cy, fy, homeBlend * 0.55);

    const scale = this.loopScale;
    const cr = cos(this.loopRot), sr = sin(this.loopRot);

    // Head parameter (theta) along the lemniscate
    const thetaH = this.loopTheta0 + this.loopThetaSpan * p;

    // Lemniscate of Gerono: x=sin(t), y=sin(2t)/2
    const sxH = sin(thetaH);
    const syH = sin(2 * thetaH) * 0.5;
    const rxH = sxH * cr - syH * sr;
    const ryH = sxH * sr + syH * cr;
    this.headX = constrain(cx + rxH * scale, margin, width - margin);
    this.headY = constrain(cy + ryH * scale * 0.92, margin, height - margin);

    // Dynamics (dt-scaled)
    const seekK = 24 + 14 * uT; // per sound-second^2
    const alignK = 16 + 8 * uT; // toward tangent velocity
    const maxSpd = (minD * 0.95) * (0.90 + map(level, 0, 0.2, 0, 0.20, true)); // px/sound-sec
    const dampP = pow(0.90, dtK);
    const flutter = (1 - uT) * 0.30 + map(level, 0, 0.2, 0, 0.55, true);

    // Ribbon width (narrow so it doesn't become a cloud)
    const headW = minD * 0.020;
    const tailW = minD * 0.010;

    // Tail span along the curve (in radians)
    const tailThetaSpan = TWO_PI * 0.95;

    for (let i = 0; i < N; i++) {
      const x = this.x[i], y = this.y[i];
      let vx = this.vx[i], vy = this.vy[i];

      const u = this.slotU[i];
      const w = lerp(headW, tailW, u);
      const lat = this.slotLat[i] * w;

      const th = thetaH - u * tailThetaSpan;
      const sx0 = sin(th);
      const sy0 = sin(2 * th) * 0.5;
      const rx0 = sx0 * cr - sy0 * sr;
      const ry0 = sx0 * sr + sy0 * cr;

      // Tangent for direction: dx=cos(t), dy=cos(2t)
      let tx = cos(th);
      let ty = cos(2 * th);
      const trx = tx * cr - ty * sr;
      const try_ = tx * sr + ty * cr;
      const tl = sqrt(trx * trx + try_ * try_) + 1e-6;
      tx = trx / tl;
      ty = try_ / tl;
      const px = -ty, py = tx;

      const n0 = noise(this.slotPhase[i], t * 0.18);
      const f = (n0 - 0.5) * flutter;

      const tx0 = cx + rx0 * scale + px * (lat + f);
      const ty0 = cy + ry0 * scale * 0.92 + py * (lat + f);

      const dx0 = tx0 - x;
      const dy0 = ty0 - y;

      // Spring toward ribbon slot
      vx += dx0 * seekK * dt;
      vy += dy0 * seekK * dt;

      // Match a forward velocity along the tangent (gives clear leading direction)
      const forward = lerp(maxSpd * 0.95, maxSpd * 0.35, u);
      vx += (tx * forward - vx) * alignK * dt;
      vy += (ty * forward - vy) * alignK * dt;

      vx *= dampP;
      vy *= dampP;

      const sp2 = vx * vx + vy * vy;
      const ms2 = maxSpd * maxSpd;
      if (sp2 > ms2) {
        const k = maxSpd / sqrt(sp2);
        vx *= k; vy *= k;
      }

      this.vx[i] = vx;
      this.vy[i] = vy;
      this.x[i] = x + vx * dt;
      this.y[i] = y + vy * dt;

      if (this.x[i] < 0) { this.x[i] = 0; this.vx[i] *= -0.25; }
      else if (this.x[i] > width) { this.x[i] = width; this.vx[i] *= -0.25; }
      if (this.y[i] < 0) { this.y[i] = 0; this.vy[i] *= -0.25; }
      else if (this.y[i] > height) { this.y[i] = height; this.vy[i] *= -0.25; }
    }
  },
  _danceCircle(level, s) {
    const minD = min(width, height);
    const margin = minD * 0.06;
    const uT = constrain(s, 0, 1);

    const dt = constrain(abs(tDelta), 0, 0.06);
    if (dt <= 0) return;
    const dtK = dt / 0.02;
    const p = uT * uT * (3 - 2 * uT);

    // Center stays fixed in the middle (simple, iconic circle).
    const cx = width * 0.5;
    const cy = height * 0.5;

    const thetaH = this.circleTheta0 + this.circleThetaSpan * p;
    const baseR = max(minD * 0.10, this.circleR);
    const thick = max(minD * 0.006, this.circleThick);

    this.headX = constrain(cx + cos(thetaH) * baseR, margin, width - margin);
    this.headY = constrain(cy + sin(thetaH) * baseR, margin, height - margin);

    const seekK = 26 + 18 * uT;
    const alignK = 16 + 10 * uT;
    const maxSpd = (minD * 0.85) * (0.90 + map(level, 0, 0.2, 0, 0.18, true));
    const dampP = pow(0.90, dtK);
    const flutter = (1 - uT) * 0.22 + map(level, 0, 0.2, 0, 0.45, true);

    // Tail span around the circle
    const tailSpan = TWO_PI * 0.95;

    for (let i = 0; i < N; i++) {
      const x = this.x[i], y = this.y[i];
      let vx = this.vx[i], vy = this.vy[i];

      const u = this.slotU[i];
      const th = thetaH - u * tailSpan;

      // Tangent (direction) around the circle
      const tx = -sin(th);
      const ty = cos(th);
      const px = -ty, py = tx;

      const n0 = noise(this.slotPhase[i], t * 0.18);
      const f = (n0 - 0.5) * flutter;

      const r = baseR + this.slotLat[i] * thick + f * thick;
      const tx0 = cx + cos(th) * r;
      const ty0 = cy + sin(th) * r;

      const dx0 = tx0 - x;
      const dy0 = ty0 - y;

      vx += dx0 * seekK * dt;
      vy += dy0 * seekK * dt;

      const forward = lerp(maxSpd * 0.95, maxSpd * 0.40, u);
      vx += (tx * forward - vx) * alignK * dt;
      vy += (ty * forward - vy) * alignK * dt;

      vx *= dampP;
      vy *= dampP;

      const sp2 = vx * vx + vy * vy;
      const ms2 = maxSpd * maxSpd;
      if (sp2 > ms2) {
        const k = maxSpd / sqrt(sp2);
        vx *= k; vy *= k;
      }

      this.vx[i] = vx;
      this.vy[i] = vy;
      this.x[i] = x + vx * dt;
      this.y[i] = y + vy * dt;

      if (this.x[i] < 0) { this.x[i] = 0; this.vx[i] *= -0.25; }
      else if (this.x[i] > width) { this.x[i] = width; this.vx[i] *= -0.25; }
      if (this.y[i] < 0) { this.y[i] = 0; this.vy[i] *= -0.25; }
      else if (this.y[i] > height) { this.y[i] = height; this.vy[i] *= -0.25; }
    }
  },
  draw() {
    noStroke();
    const aScene = typeof sceneA === "number" ? sceneA : 1;
    if (aScene <= 0.004) return;
    // Base: white lights (additive). During Stage A we render a head->tail alpha gradient
    // so the ribbon reads as "less saturated start" -> "more saturated end" like the reference.
    blendMode(ADD);
    strokeCap(ROUND);

    const p = Acts && Acts.mode === "TRANSITION"
      ? constrain((t - Acts.transitionStartT) / max(1e-6, TRANSITION_DURATION), 0, 1)
      : 1;
    const stageA = Acts && Acts.mode === "TRANSITION" && p < DANCE_PORTION;

    const sizeScale = (typeof cellW === "number" && typeof CELL_SIZE === "number" && CELL_SIZE > 0)
      ? (cellW / CELL_SIZE)
      : 1;
    const swBase = max(1.0, PARTICLE_SIZE * sizeScale);

    if (!stageA) {
      strokeWeight(swBase);
      stroke(255, 255, 255, LIGHT_ALPHA * aScene);
      for (let i = 0; i < N; i++) point(this.x[i], this.y[i]);
    } else {
      // Tail (faded)
      strokeWeight(max(1.0, swBase - 0.2 * sizeScale));
      stroke(255, 255, 255, (LIGHT_ALPHA * 0.28) * aScene);
      for (let i = 0; i < N; i++) if (this.slotU[i] > 0.58) point(this.x[i], this.y[i]);

      // Mid
      strokeWeight(swBase);
      stroke(255, 255, 255, (LIGHT_ALPHA * 0.62) * aScene);
      for (let i = 0; i < N; i++) if (this.slotU[i] > 0.22 && this.slotU[i] <= 0.58) point(this.x[i], this.y[i]);

      // Head (more saturated)
      strokeWeight(swBase + 1.2 * sizeScale);
      stroke(255, 255, 255, min(255, (LIGHT_ALPHA * 1.05) * aScene));
      for (let i = 0; i < N; i++) if (this.slotU[i] <= 0.22) point(this.x[i], this.y[i]);
    }

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
      const aWash = (Render.tintA * aScene) | 0;
      stroke(red(Render.hi), green(Render.hi), blue(Render.hi), aWash);
      for (let i = 0; i < N; i++) {
        const dx = this.x[i] - fx, dy = this.y[i] - fy;
        const d2 = dx * dx + dy * dy;
        if (d2 <= outer2) point(this.x[i], this.y[i]);
      }

      // Minority tint: a single "stain" blob sitting on top of the wash (irregular edge via noise).
      const aStain = (Render.tintB * aScene) | 0;
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

Particles.onTransitionStart = function (fromAct, toAct) {
  // Stage A: set leader + 2-turn waypoint navigation (S-curve).
  const minD = min(width, height);
  const margin = minD * 0.08;

  // Centroid becomes the leader start.
  let sx = 0, sy = 0;
  for (let i = 0; i < N; i++) { sx += this.x[i]; sy += this.y[i]; }
  sx /= N; sy /= N;
  sx = constrain(sx, margin, width - margin);
  sy = constrain(sy, margin, height - margin);

  // Special case: act4 -> act1 forms a centered circle.
  if ((fromAct | 0) === 4 && (toAct | 0) === 1) {
    this.transitionMode = 1;
    this.circleR = minD * 0.22;
    this.circleThick = minD * 0.020;
    this.circleTheta0 = 0;
    this.circleThetaSpan = TWO_PI * 1.25;

    // Takeoff: nudge velocities tangentially so it becomes a circle fast.
    const kick = minD * 0.48;
    for (let i = 0; i < N; i++) {
      const u = this.slotU[i];
      const th = -u * TWO_PI * 0.95;
      const tx = -sin(th);
      const ty = cos(th);
      this.vx[i] = this.vx[i] * 0.25 + tx * kick;
      this.vy[i] = this.vy[i] * 0.25 + ty * kick;
      this.trA[i] = 0;
    }
    return;
  }

  this.transitionMode = 0;

  // Pick a far end point in the opposite-ish quadrant so we travel across the screen.
  const minFar2 = (minD * 0.72) * (minD * 0.72);
  let ex = sx, ey = sy;
  for (let k = 0; k < 18; k++) {
    const wantRight = sx < width * 0.5;
    const wantDown = sy < height * 0.5;
    const rx = wantRight ? random(width * 0.62, width - margin) : random(margin, width * 0.38);
    const ry = wantDown ? random(height * 0.62, height - margin) : random(margin, height * 0.38);
    const dx = rx - sx, dy = ry - sy;
    if (dx * dx + dy * dy >= minFar2) { ex = rx; ey = ry; break; }
    ex = rx; ey = ry;
  }

  const vx0 = ex - sx, vy0 = ey - sy;
  const vl = sqrt(vx0 * vx0 + vy0 * vy0) + 1e-6;
  const ux = vx0 / vl, uy = vy0 / vl;

  // Lemniscate ribbon: start near the current silhouette, travel across the canvas,
  // and end near a far point (with late pull toward next focus handled in `dance()`).
  this.loopCX0 = sx; this.loopCY0 = sy;
  this.loopCX1 = ex; this.loopCY1 = ey;
  this.loopScale = min(width, height) * 0.42;
  // Keep the ∞ mostly horizontal like the reference, with a slight travel-aligned bias.
  this.loopRot = constrain(atan2(vy0, vx0) * 0.18 + (random() - 0.5) * 0.10, -0.45, 0.45);
  this.loopTheta0 = 0.0;                // start at center (no jump at transition start)
  this.loopThetaSpan = TWO_PI * 1.12;   // draw a full ∞ and a bit more

  // Takeoff: nudge velocities along the initial tangent so the ribbon "forms" immediately.
  const cr = cos(this.loopRot), sr = sin(this.loopRot);
  // Tangent at theta0=0: dx=cos(0)=1, dy=cos(0)=1
  let tx = 1 * cr - 1 * sr;
  let ty = 1 * sr + 1 * cr;
  const tl = sqrt(tx * tx + ty * ty) + 1e-6;
  tx /= tl; ty /= tl;
  const px = -ty, py = tx;

  const kick = minD * 0.42;
  const spread = minD * 0.08;
  for (let i = 0; i < N; i++) {
    const lat = this.slotLat[i];
    this.vx[i] = this.vx[i] * 0.25 + tx * kick + px * lat * spread;
    this.vy[i] = this.vy[i] * 0.25 + ty * kick + py * lat * spread;
    this.trA[i] = 0;
  }
};
