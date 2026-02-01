let mic, amp;
let audioStarted = false;

// זמן פנימי – מתקדם רק כשיש קול
let t = 0;

// רגישות מיקרופון
let threshold = 0.03;

// מערכות
let act = 1;          // 1–4
let actDuration = 6;  // "שניות קול" לכל מערכה
let actStartT = 0;
let framesPerSecond = 24;

// STEP 7: image sequences (already loaded in preload)
let actFrames = { 1: [], 2: [], 3: [], 4: [] };
let actFrameCounts = { 1: 192, 2: 172, 3: 0, 4: 0 };

// STEP 8: state machine + pooled material particles (mobile-friendly)
let mode = "ACT"; // "ACT" | "TRANSITION"
let nextAct = 2;
let transitionStartT = 0;
let transitionDuration = 2.0;
let transitionFromImg = null;
let transitionToImg = null;

// Small buffer for one-time target sampling (no per-frame loadPixels)
let targetG;
let targetW = 160;
let targetH = 284;

// Target points (preallocated, reused)
let maxTargets = 2600;
let targetX = new Float32Array(maxTargets);
let targetY = new Float32Array(maxTargets);
let targetCount = 0;
let _cdf = new Float32Array(targetW * targetH);
let _cdfTotal = 0;
let _usedPix = new Uint8Array(targetW * targetH);

// Transition particles (preallocated pool, reused)
let baseParticleCount = 2200;
let maxParticles = 2600;
let particleCount = baseParticleCount;
let pX = new Float32Array(maxParticles);
let pY = new Float32Array(maxParticles);
let pVX = new Float32Array(maxParticles);
let pVY = new Float32Array(maxParticles);
let pTarget = new Uint16Array(maxParticles);
let pKind = new Uint8Array(maxParticles); // 0=white, 1=black
let pSeedX = new Float32Array(maxParticles);
let pSeedY = new Float32Array(maxParticles);
let pUx = new Float32Array(maxParticles);
let pUy = new Float32Array(maxParticles);
let pJx = new Float32Array(maxParticles);
let pJy = new Float32Array(maxParticles);
let pCurTX = new Float32Array(maxParticles);
let pCurTY = new Float32Array(maxParticles);
let pNextTX = new Float32Array(maxParticles);
let pNextTY = new Float32Array(maxParticles);

// ACT mode: targets sampled only when frameIndex changes
let currentFrameIndex = 0;
let desiredFrameIndex = 0;
let lastSampledAct = 1;
let lastSampledFrameIndex = -1;

// Debug flag (must remain false): no PNG frames drawn to main canvas
let imagesDrawnToCanvas = false;

// Freeze all motion when t is frozen (silence)
let _prevT = 0;
let tAdvanced = false;

// Settling metrics (computed each draw; lightweight sampling)
let settleMeanDist = 0;
let settlePct = 0;
let settleIsReady = false;
let settleSampleCount = 200;
let settledRadius = 10;
let settleMeanThreshold = 14;
let settlePctThreshold = 0.7;

// Reusable rect objects to avoid per-frame allocations
let _rectA = { x: 0, y: 0, w: 0, h: 0 };
let _rectB = { x: 0, y: 0, w: 0, h: 0 };

function preload() {
  for (let a = 1; a <= 4; a++) {
    const count = actFrameCounts[a] || 0;
    actFrames[a] = new Array(count);
    for (let i = 0; i < count; i++) {
      const path = `assets/act${a}/act${a}_${nf(i, 5)}.png`;
      actFrames[a][i] = loadImage(path);
    }
  }
}

function setup() {
  createCanvas(windowWidth, windowHeight);
  background(0);

  mic = new p5.AudioIn();
  amp = new p5.Amplitude();

  initMaterial();
  initParticles();
  updateActTiming();
  lastSampledAct = act;
  lastSampledFrameIndex = -1;
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
  initMaterial();
  initParticles();
  updateActTiming();
  lastSampledAct = act;
  lastSampledFrameIndex = -1;
}

function updateActTiming() {
  const n = (actFrames[act] || []).length;
  if (n > 0 && framesPerSecond > 0) actDuration = n / framesPerSecond;
}

function initMaterial() {
  targetG = createGraphics(targetW, targetH);
  targetG.pixelDensity(1);
  targetG.noSmooth();
  particleCount = constrain(baseParticleCount, 200, maxParticles);
}

function initParticles() {
  const n = min(particleCount, maxParticles);
  const cx = width * 0.5;
  const cy = height * 0.5;
  const R = min(width, height) * 0.22;

  for (let i = 0; i < n; i++) {
    const ang = random(TWO_PI);
    const rad = sqrt(random()) * R;
    pX[i] = cx + cos(ang) * rad;
    pY[i] = cy + sin(ang) * rad;
    pVX[i] = random(-0.6, 0.6);
    pVY[i] = random(-0.6, 0.6);
    pKind[i] = random() < 0.1 ? 1 : 0; // always include black particles
    pTarget[i] = 0;
    pSeedX[i] = random(1000);
    pSeedY[i] = random(1000);
    pUx[i] = random(targetW);
    pUy[i] = random(targetH);
    pJx[i] = (noise(pSeedX[i], pSeedY[i], 0) - 0.5) * 4.0;
    pJy[i] = (noise(pSeedX[i] + 9.1, pSeedY[i] + 3.7, 0) - 0.5) * 4.0;
    pCurTX[i] = pX[i];
    pCurTY[i] = pY[i];
    pNextTX[i] = pX[i];
    pNextTY[i] = pY[i];
  }
}

function draw() {
  background(0);

  if (!audioStarted) {
    drawStartScreen();
    return;
  }

  const level = amp.getLevel();

  // אם יש קול — הזמן מתקדם
  if (level > threshold) {
    const speed = map(level, threshold, 0.2, 0.01, 0.05, true);
    t += speed;
  }
  tAdvanced = t !== _prevT;
  _prevT = t;

  updateAct();
  drawAct(level);
  drawDebugOverlay();
}

function mousePressed() {
  if (!audioStarted) {
    // חייב להיות בתוך gesture (לחיצה)
    userStartAudio();

    mic.start(() => {
      amp.setInput(mic);
      audioStarted = true;
      actStartT = t;
    });
  }
}

function drawStartScreen() {
  fill(255);
  noStroke();
  textAlign(CENTER, CENTER);
  textSize(18);
  text("Tap to activate microphone\nSound drives the animation", width / 2, height / 2);
}

function computeFitRectInto(img, w, h, out) {
  const s = min(w / img.width, h / img.height);
  const dw = img.width * s;
  const dh = img.height * s;
  out.x = (w - dw) * 0.5;
  out.y = (h - dh) * 0.5;
  out.w = dw;
  out.h = dh;
  return out;
}

function getNextActWithFrames(fromAct) {
  let a = fromAct;
  for (let k = 0; k < 4; k++) {
    a++;
    if (a > 4) a = 1;
    const frames = actFrames[a] || [];
    if (frames.length > 0) return a;
  }
  return fromAct;
}

function getActFrameAt(actNum, localT) {
  const frames = actFrames[actNum] || [];
  const n = frames.length;
  if (n === 0) return null;
  const progress = actDuration > 0 ? localT / actDuration : 0;
  const idx = constrain(floor(progress * n), 0, n - 1);
  return frames[idx] || null;
}

function prepareTargets(img) {
  targetCount = 0;
  if (!img || !targetG) return;

  // Draw next image into small buffer
  const gRect = computeFitRectInto(img, targetW, targetH, _rectA);
  targetG.clear();
  targetG.background(0);
  targetG.imageMode(CORNER);
  targetG.image(img, gRect.x, gRect.y, gRect.w, gRect.h);

  // One-time sampling at transition start (never in draw loop)
  targetG.loadPixels();

  // Map sampled points into the same fit rect on the canvas
  const cRect = computeFitRectInto(img, width, height, _rectB);

  _usedPix.fill(0);
  _cdfTotal = 0;

  // Build weighted CDF over the whole buffer area (weights proportional to brightness * alpha)
  for (let y = 0; y < targetH; y++) {
    for (let x = 0; x < targetW; x++) {
      const p = y * targetW + x;

      let w = 0;
      if (x >= gRect.x && x < gRect.x + gRect.w && y >= gRect.y && y < gRect.y + gRect.h) {
        const i = 4 * p;
        const r = targetG.pixels[i + 0];
        const g = targetG.pixels[i + 1];
        const b = targetG.pixels[i + 2];
        const a = targetG.pixels[i + 3];

        const an = a / 255;
        if (an > 0.05) {
          const br = (r + g + b) / 765; // 0..1
          if (br > 0.08) {
            // Bias toward brighter pixels (more particles) while keeping full coverage
            w = br * br * an;
          }
        }
      }

      _cdfTotal += w;
      _cdf[p] = _cdfTotal;
    }
  }

  const nTargets = min(particleCount, maxTargets);

  // Fallback: center cluster if image is too dark / empty
  if (_cdfTotal <= 0.000001) {
    const cx = width * 0.5;
    const cy = height * 0.5;
    for (let i = 0; i < nTargets; i++) {
      targetX[i] = cx + randomGaussian() * min(width, height) * 0.12;
      targetY[i] = cy + randomGaussian() * min(width, height) * 0.12;
    }
    targetCount = nTargets;
    return;
  }

  // Sample exactly nTargets target points (one per particle), with jitter to avoid stacking
  const jitter = 2.0;
  const gwInv = 1 / max(1, gRect.w);
  const ghInv = 1 / max(1, gRect.h);

  for (let i = 0; i < nTargets; i++) {
    let p = 0;
    let tries = 0;

    while (tries < 12) {
      const r = random(_cdfTotal) + 1e-9;
      let lo = 0;
      let hi = _cdf.length - 1;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (_cdf[mid] >= r) hi = mid;
        else lo = mid + 1;
      }
      p = lo;
      if (_usedPix[p] === 0) break;
      tries++;
    }

    _usedPix[p] = 1;

    const x = p % targetW;
    const y = (p / targetW) | 0;

    const ix = (x - gRect.x) * gwInv;
    const iy = (y - gRect.y) * ghInv;

    targetX[i] = cRect.x + ix * cRect.w + random(-jitter, jitter);
    targetY[i] = cRect.y + iy * cRect.h + random(-jitter, jitter);
  }

  targetCount = nTargets;
}

function prepareTargetsFromImgForTransition(img) {
  prepareTargets(img);
  // Scramble is allowed during transitions
  const n = min(particleCount, maxParticles);
  const tc = max(1, targetCount);
  for (let i = 0; i < n; i++) {
    pTarget[i] = floor(random(tc));
  }
}

function _bufWeightAt(x, y) {
  const xi = x | 0;
  const yi = y | 0;
  if (xi < 0 || xi >= targetW || yi < 0 || yi >= targetH) return 0;
  const p = yi * targetW + xi;
  const i = 4 * p;
  const a = targetG.pixels[i + 3] / 255;
  if (a <= 0.05) return 0;
  const r = targetG.pixels[i + 0];
  const g = targetG.pixels[i + 1];
  const b = targetG.pixels[i + 2];
  const br = (r + g + b) / 765;
  if (br <= 0.08) return 0;
  return br * br * a;
}

let _nbX = 0;
let _nbY = 0;

function _nearestBrightAround(ux, uy, gRect) {
  let bestX = constrain(ux, gRect.x, gRect.x + gRect.w - 1);
  let bestY = constrain(uy, gRect.y, gRect.y + gRect.h - 1);
  let bestW = _bufWeightAt(bestX, bestY);
  if (bestW > 0) {
    _nbX = bestX;
    _nbY = bestY;
    return;
  }

  // Sample a small neighborhood without heavy full scans
  const dirsX = [1, 0, -1, 0, 1, -1, -1, 1, 2, 0, -2, 0, 2, -2, -2, 2];
  const dirsY = [0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 2, 2, -2, -2];

  for (let r = 2; r <= 6; r++) {
    for (let k = 0; k < dirsX.length; k++) {
      const x = constrain(ux + dirsX[k] * r, gRect.x, gRect.x + gRect.w - 1);
      const y = constrain(uy + dirsY[k] * r, gRect.y, gRect.y + gRect.h - 1);
      const w = _bufWeightAt(x, y);
      if (w > bestW) {
        bestW = w;
        bestX = x;
        bestY = y;
      }
    }
    if (bestW > 0) break;
  }

  _nbX = bestX;
  _nbY = bestY;
}

function sampleActFrameIntoTargets(img, isActStart) {
  if (!img || !targetG) return;

  const gRect = computeFitRectInto(img, targetW, targetH, _rectA);
  targetG.clear();
  targetG.background(0);
  targetG.imageMode(CORNER);
  targetG.image(img, gRect.x, gRect.y, gRect.w, gRect.h);
  targetG.loadPixels();

  const cRect = computeFitRectInto(img, width, height, _rectB);

  // If starting a new act, choose persistent (ux,uy) biased to the object area
  if (isActStart) {
    _usedPix.fill(0);
    _cdfTotal = 0;

    for (let y = 0; y < targetH; y++) {
      for (let x = 0; x < targetW; x++) {
        const p = y * targetW + x;
        let w = 0;
        if (x >= gRect.x && x < gRect.x + gRect.w && y >= gRect.y && y < gRect.y + gRect.h) {
          const i = 4 * p;
          const a = targetG.pixels[i + 3] / 255;
          if (a > 0.05) {
            const r = targetG.pixels[i + 0];
            const g = targetG.pixels[i + 1];
            const b = targetG.pixels[i + 2];
            const br = (r + g + b) / 765;
            if (br > 0.08) w = br * br * a;
          }
        }
        _cdfTotal += w;
        _cdf[p] = _cdfTotal;
      }
    }

    const n = min(particleCount, maxParticles);
    const total = _cdfTotal;

    if (total > 0.000001) {
      for (let i = 0; i < n; i++) {
        const r = random(total) + 1e-9;
        let lo = 0;
        let hi = _cdf.length - 1;
        while (lo < hi) {
          const mid = (lo + hi) >> 1;
          if (_cdf[mid] >= r) hi = mid;
          else lo = mid + 1;
        }
        const p = lo;
        pUx[i] = p % targetW;
        pUy[i] = (p / targetW) | 0;
      }
    }
  }

  // For every particle: compute its next target from its own (ux,uy) + local bright search
  const n = min(particleCount, maxParticles);
  const gwInv = 1 / max(1, gRect.w);
  const ghInv = 1 / max(1, gRect.h);

  for (let i = 0; i < n; i++) {
    _nearestBrightAround(pUx[i], pUy[i], gRect);
    const ix = (_nbX - gRect.x) * gwInv;
    const iy = (_nbY - gRect.y) * ghInv;

    pNextTX[i] = cRect.x + ix * cRect.w + pJx[i];
    pNextTY[i] = cRect.y + iy * cRect.h + pJy[i];

    if (isActStart) {
      pCurTX[i] = pNextTX[i];
      pCurTY[i] = pNextTY[i];
    }
  }
}

function startTransition() {
  mode = "TRANSITION";
  transitionStartT = t;
  nextAct = getNextActWithFrames(act);

  transitionToImg = (actFrames[nextAct] || [])[0] || null;
  prepareTargetsFromImgForTransition(transitionToImg);
}

function endTransition() {
  act = nextAct;
  actStartT = t;
  mode = "ACT";
  updateActTiming();
  currentFrameIndex = 0;
  desiredFrameIndex = 0;
  lastSampledAct = act;
  lastSampledFrameIndex = -1;
}

function updateAct() {
  if (mode === "ACT") {
    const frames = actFrames[act] || [];
    const n = frames.length;
    if (n === 0) return;

    const actElapsed = t - actStartT;
    const desiredIdx = constrain(floor(max(0, actElapsed) * framesPerSecond), 0, n - 1);

    // Start act-to-act transition only after we reached the last frame AND the particles are settled.
    if (desiredIdx >= n - 1 && currentFrameIndex >= n - 1 && settleIsReady) startTransition();
    return;
  }

  if (mode === "TRANSITION") {
    const elapsed = t - transitionStartT;
    if (elapsed >= transitionDuration) endTransition();
  }
}

function updateActParticles(level) {
  if (!tAdvanced) return;
  const n = min(particleCount, maxParticles);
  const cx = width * 0.5;
  const cy = height * 0.5;

  const flow = 0.03 + map(level, 0, 0.2, 0, 0.10, true);
  const swirlBase = 0.45 + map(level, 0, 0.2, 0, 0.55, true);
  const driftX = 0.015;
  const driftY = -0.008;

  const damp = 0.78;
  const maxSpd = 3.8;
  const att = 0.28 + map(level, 0, 0.2, 0, 0.10, true);
  const targetLerp = 0.18;
  const snapRadius = 6;
  const wobbleAmpBase = 0.25 + map(level, 0, 0.2, 0, 0.55, true);
  const wobbleFreq = 0.35;

  for (let i = 0; i < n; i++) {
    pCurTX[i] += (pNextTX[i] - pCurTX[i]) * targetLerp;
    pCurTY[i] += (pNextTY[i] - pCurTY[i]) * targetLerp;

    const nx = pX[i] * 0.0020;
    const ny = pY[i] * 0.0020;
    const ang = noise(nx, ny, t * 0.28) * TWO_PI * 2.0;

    const dx = pX[i] - cx;
    const dy = pY[i] - cy;
    const d = sqrt(dx * dx + dy * dy) + 0.001;
    const centerBoost = constrain(1.0 - d / (min(width, height) * 0.65), 0, 1);
    const swirl = swirlBase * centerBoost;

    pVX[i] += (cos(ang) * flow + driftX + (-dy / d) * swirl) * 0.22;
    pVY[i] += (sin(ang) * flow + driftY + (dx / d) * swirl) * 0.22;

    const nx1 = noise(pSeedX[i], pSeedY[i], t * wobbleFreq);
    const ny1 = noise(pSeedX[i] + 17.3, pSeedY[i] + 91.7, t * wobbleFreq);
    const tx0 = pCurTX[i];
    const ty0 = pCurTY[i];
    const dx0 = tx0 - pX[i];
    const dy0 = ty0 - pY[i];
    const d0 = sqrt(dx0 * dx0 + dy0 * dy0);
    const wobbleAmp = wobbleAmpBase * constrain(d0 / 12, 0, 1);
    const wx = (nx1 - 0.5) * 2.0 * wobbleAmp;
    const wy = (ny1 - 0.5) * 2.0 * wobbleAmp;

    const desiredX = tx0 + wx;
    const desiredY = ty0 + wy;
    const dxT = desiredX - pX[i];
    const dyT = desiredY - pY[i];
    pVX[i] += dxT * att;
    pVY[i] += dyT * att;

    pVX[i] *= damp;
    pVY[i] *= damp;

    const spd = sqrt(pVX[i] * pVX[i] + pVY[i] * pVY[i]);
    if (spd > maxSpd) {
      const k = maxSpd / max(0.001, spd);
      pVX[i] *= k;
      pVY[i] *= k;
    }

    pX[i] += pVX[i];
    pY[i] += pVY[i];

    // Snap when very close, to fully form a crisp silhouette
    const d2 = dxT * dxT + dyT * dyT;
    if (d2 < snapRadius * snapRadius) {
      pX[i] = lerp(pX[i], tx0, 0.35);
      pY[i] = lerp(pY[i], ty0, 0.35);
      pVX[i] *= 0.6;
      pVY[i] *= 0.6;
    }

    // Soft containment (avoid wrap artifacts in ACT)
    if (pX[i] < 0) {
      pX[i] = 0;
      pVX[i] *= -0.25;
    } else if (pX[i] > width) {
      pX[i] = width;
      pVX[i] *= -0.25;
    }

    if (pY[i] < 0) {
      pY[i] = 0;
      pVY[i] *= -0.25;
    } else if (pY[i] > height) {
      pY[i] = height;
      pVY[i] *= -0.25;
    }
  }
}

function drawParticles() {
  const n = min(particleCount, maxParticles);
  noStroke();
  for (let i = 0; i < n; i++) {
    if (pKind[i] === 1) fill(0, 210);
    else fill(255, 235);
    circle(pX[i], pY[i], 2.1);
  }
}

function computeSettlingMetrics() {
  const n = min(particleCount, maxParticles);
  if (n <= 0) {
    settleMeanDist = 0;
    settlePct = 1;
    settleIsReady = true;
    return;
  }

  const samples = min(settleSampleCount, n);
  const stride = max(1, floor(n / samples));
  let sum = 0;
  let settled = 0;
  const r2 = settledRadius * settledRadius;

  for (let s = 0; s < samples; s++) {
    const i = s * stride;
    const dx = pX[i] - pNextTX[i];
    const dy = pY[i] - pNextTY[i];
    const d2 = dx * dx + dy * dy;
    sum += sqrt(d2);
    if (d2 < r2) settled++;
  }

  settleMeanDist = sum / samples;
  settlePct = settled / samples;
  settleIsReady = settleMeanDist < settleMeanThreshold || settlePct > settlePctThreshold;
}

function updateTransitionParticles(level, transP) {
  if (!tAdvanced) return;
  const n = min(particleCount, maxParticles);
  const cx = width * 0.5;
  const cy = height * 0.5;

  // Stage A: 0..0.6 DANCE (loose, choreographic)
  // Stage B: 0.6..1 COALESCE (attract to next-shape targets)
  const danceEnd = 0.6;
  const stageA = transP < danceEnd;
  const stageT = stageA ? transP / max(0.0001, danceEnd) : (transP - danceEnd) / max(0.0001, 1.0 - danceEnd);

  const flow = 0.22 + map(level, 0, 0.2, 0, 0.65, true);
  const swirlBase = 0.9 + map(level, 0, 0.2, 0, 1.2, true);

  // Drift guides the overall motion; swirl keeps it centered
  const driftX = lerp(0.25, 0.05, stageT);
  const driftY = lerp(-0.05, -0.15, stageT);

  const damp = stageA ? 0.92 : 0.82;
  const maxSpd = stageA ? 5.5 : 4.0;

  const att = stageA ? 0.0 : lerp(0.03, 0.22, stageT);
  const snap = stageA ? 0.0 : map(transP, 0.9, 1.0, 0, 1, true);
  const wobbleAmpBase = 0.6 + map(level, 0, 0.2, 0, 1.2, true);
  const wobbleAmp = stageA ? 0.0 : wobbleAmpBase * (1.0 - snap * 0.85);
  const wobbleFreq = 0.35;

  for (let i = 0; i < n; i++) {
    const nx = pX[i] * 0.0021;
    const ny = pY[i] * 0.0021;
    const ang = noise(nx, ny, t * 0.35) * TWO_PI * 2.4;

    const dx = pX[i] - cx;
    const dy = pY[i] - cy;
    const d = sqrt(dx * dx + dy * dy) + 0.001;
    const centerBoost = constrain(1.0 - d / (min(width, height) * 0.55), 0, 1);

    const swirl = swirlBase * centerBoost;
    const ax = cos(ang) * flow + driftX + (-dy / d) * swirl;
    const ay = sin(ang) * flow + driftY + (dx / d) * swirl;

    pVX[i] += ax * 0.18;
    pVY[i] += ay * 0.18;

    if (!stageA && targetCount > 0) {
      const ti = pTarget[i] % targetCount;
      const tx = targetX[ti];
      const ty = targetY[ti];

      const nx1 = noise(pSeedX[i], pSeedY[i], t * wobbleFreq);
      const ny1 = noise(pSeedX[i] + 17.3, pSeedY[i] + 91.7, t * wobbleFreq);
      const wx = (nx1 - 0.5) * 2.0 * wobbleAmp;
      const wy = (ny1 - 0.5) * 2.0 * wobbleAmp;

      const dxT = (tx + wx) - pX[i];
      const dyT = (ty + wy) - pY[i];
      pVX[i] += dxT * att;
      pVY[i] += dyT * att;

      if (snap > 0) {
        pX[i] = lerp(pX[i], tx, snap * 0.35);
        pY[i] = lerp(pY[i], ty, snap * 0.35);
        pVX[i] *= 1.0 - snap * 0.65;
        pVY[i] *= 1.0 - snap * 0.65;
      }
    }

    pVX[i] *= damp;
    pVY[i] *= damp;

    const spd = sqrt(pVX[i] * pVX[i] + pVY[i] * pVY[i]);
    if (spd > maxSpd) {
      const k = maxSpd / max(0.001, spd);
      pVX[i] *= k;
      pVY[i] *= k;
    }

    pX[i] += pVX[i];
    pY[i] += pVY[i];

    // Soft containment (avoid wrap rectangle artifacts)
    if (pX[i] < 0) {
      pX[i] = 0;
      pVX[i] *= -0.25;
    } else if (pX[i] > width) {
      pX[i] = width;
      pVX[i] *= -0.25;
    }

    if (pY[i] < 0) {
      pY[i] = 0;
      pVY[i] *= -0.25;
    } else if (pY[i] > height) {
      pY[i] = height;
      pVY[i] *= -0.25;
    }
  }
}

function drawTransitionParticles(level, transP) {
  const n = min(particleCount, maxParticles);
  const danceEnd = 0.6;
  const stageA = transP < danceEnd;
  const stageT = stageA ? transP / max(0.0001, danceEnd) : (transP - danceEnd) / max(0.0001, 1.0 - danceEnd);

  const aWhite = stageA ? (200 + level * 40) : (230 + stageT * 15);
  const aBlack = stageA ? (160 + level * 35) : (210 + stageT * 15);
  const sizeBase = stageA ? (2.2 + level * 1.2) : (2.0 + stageT * 0.6);

  noStroke();
  for (let i = 0; i < n; i++) {
    if (pKind[i] === 1) fill(0, aBlack);
    else fill(255, aWhite);

    const j = (noise(pX[i] * 0.004, pY[i] * 0.004, t * 0.25) - 0.5) * 0.25;
    circle(pX[i], pY[i], max(1.4, sizeBase + j));
  }
}

function drawAct(level) {
  if (mode === "ACT") {
    const frames = actFrames[act] || [];
    const n = frames.length;

    if (n === 0) {
      fill(255);
      noStroke();
      textAlign(CENTER, CENTER);
      textSize(24);
      text(`ACT ${act}`, width / 2, height / 2);
      return;
    }

    const actElapsed = max(0, t - actStartT);
    desiredFrameIndex = constrain(floor(actElapsed * framesPerSecond), 0, n - 1);

    // Gate progression: only advance one frame at a time when sufficiently settled.
    if (currentFrameIndex < 0 || currentFrameIndex >= n) currentFrameIndex = 0;

    // Ensure current frame targets are loaded (act start or after transitions)
    if (act !== lastSampledAct || currentFrameIndex !== lastSampledFrameIndex) {
      const isActStart = act !== lastSampledAct || lastSampledFrameIndex === -1;
      const curImg = frames[currentFrameIndex];
      if (curImg) sampleActFrameIntoTargets(curImg, isActStart);
      lastSampledAct = act;
      lastSampledFrameIndex = currentFrameIndex;
    }

    updateActParticles(level);
    computeSettlingMetrics();

    if (settleIsReady && currentFrameIndex < desiredFrameIndex) {
      currentFrameIndex = min(currentFrameIndex + 1, n - 1);
      const nextImg = frames[currentFrameIndex];
      if (nextImg) sampleActFrameIntoTargets(nextImg, false);
      lastSampledAct = act;
      lastSampledFrameIndex = currentFrameIndex;
      computeSettlingMetrics();
    }

    drawParticles();
    return;
  }

  if (mode === "TRANSITION") {
    const transP = constrain((t - transitionStartT) / max(0.0001, transitionDuration), 0, 1);

    updateTransitionParticles(level, transP);

    drawTransitionParticles(level, transP);
  }
}

function drawDebugOverlay() {
  push();
  noStroke();
  fill(0, 160);
  rect(10, 10, 330, 134, 8);
  fill(255);
  textAlign(LEFT, TOP);
  textSize(12);
  text(`mode: ${mode}`, 18, 16);
  text(`act: ${act}  next: ${mode === "TRANSITION" ? nextAct : getNextActWithFrames(act)}`, 18, 32);
  text(`desiredFrame: ${desiredFrameIndex}`, 18, 48);
  text(`currentFrame: ${currentFrameIndex}`, 18, 64);
  text(`meanDist: ${settleMeanDist.toFixed(2)}`, 18, 80);
  text(`settledPct: ${(settlePct * 100).toFixed(1)}%`, 18, 96);
  text(`actElapsed: ${(t - actStartT).toFixed(3)}`, 18, 112);
  text(`actStartT: ${actStartT.toFixed(3)}`, 175, 112);
  text(`particles: ${particleCount}`, 18, 128);
  text(`imagesDrawn: ${imagesDrawnToCanvas}`, 175, 128);
  pop();
}
