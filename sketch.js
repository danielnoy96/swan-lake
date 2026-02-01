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

// Sampling scratch buffers (CDF + used pixels) for the small offscreen buffer
let _cdf = new Float32Array(targetW * targetH);
let _cdfTotal = 0;
let _usedPix = new Uint8Array(targetW * targetH);
let _pickCount = new Uint16Array(targetW * targetH);
let lastGRect = { x: 0, y: 0, w: 0, h: 0 };
let lastCRect = { x: 0, y: 0, w: 0, h: 0 };

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
let pJx = new Float32Array(maxParticles);
let pJy = new Float32Array(maxParticles);
let pOrder = new Uint16Array(maxParticles); // stable correspondence across frames
let pUx = new Float32Array(maxParticles);
let pUy = new Float32Array(maxParticles);
let pTX = new Float32Array(maxParticles); // smoothed target
let pTY = new Float32Array(maxParticles);
let pDTX = new Float32Array(maxParticles); // desired target (on frame change)
let pDTY = new Float32Array(maxParticles);

// ACT mode: targets sampled only when frameIndex changes
let currentFrameIndex = 0;
let lastSampledAct = 1;
let lastSampledFrameIndex = -1;

// Debug flag (must remain false): no PNG frames drawn to main canvas
let imagesDrawnToCanvas = false;

// Freeze all motion when t is frozen (silence)
let _prevT = 0;
let tAdvanced = false;

// Debug (ACT settling + target sampling)
let debugMeanDist = 0;
let debugSettledPct = 0;
let debugLocalFound = 0;

// Local search
const _dirsX = [1, 0, -1, 0, 1, -1, -1, 1, 2, 0, -2, 0, 2, -2, -2, 2];
const _dirsY = [0, 1, 0, -1, 1, 1, -1, -1, 0, 2, 0, -2, 2, 2, -2, -2];

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

function _bufWeightAtXY(x, y) {
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
  if (br <= 0.06) return 0;
  return br * br * a;
}

function loadFrameToBuffer(img) {
  const gRect = computeFitRectInto(img, targetW, targetH, _rectA);
  const cRect = computeFitRectInto(img, width, height, _rectB);

  lastGRect.x = gRect.x;
  lastGRect.y = gRect.y;
  lastGRect.w = gRect.w;
  lastGRect.h = gRect.h;

  lastCRect.x = cRect.x;
  lastCRect.y = cRect.y;
  lastCRect.w = cRect.w;
  lastCRect.h = cRect.h;

  targetG.clear();
  targetG.background(0);
  targetG.imageMode(CORNER);
  targetG.image(img, gRect.x, gRect.y, gRect.w, gRect.h);
  targetG.loadPixels();
}

function assignPersistentUxUyFromCurrentBuffer(scramble) {
  const n = min(particleCount, maxParticles);

  _cdfTotal = 0;
  for (let y = 0; y < targetH; y++) {
    for (let x = 0; x < targetW; x++) {
      const p = y * targetW + x;
      let w = 0;
      if (
        x >= lastGRect.x &&
        x < lastGRect.x + lastGRect.w &&
        y >= lastGRect.y &&
        y < lastGRect.y + lastGRect.h
      ) {
        w = _bufWeightAtXY(x, y);
      }
      _cdfTotal += w;
      _cdf[p] = _cdfTotal;
    }
  }

  if (_cdfTotal <= 0.000001) {
    for (let i = 0; i < n; i++) {
      pUx[i] = lastGRect.x + random(lastGRect.w);
      pUy[i] = lastGRect.y + random(lastGRect.h);
    }
    return;
  }

  for (let i = 0; i < n; i++) {
    const r = scramble ? (random(_cdfTotal) + 1e-9) : (((pOrder[i] + 0.5) / n) * _cdfTotal + 1e-9);
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

function updateDesiredTargetsFromBuffer() {
  const n = min(particleCount, maxParticles);
  _pickCount.fill(0);
  debugLocalFound = 0;

  const gwInv = 1 / max(1, lastGRect.w);
  const ghInv = 1 / max(1, lastGRect.h);

  for (let i = 0; i < n; i++) {
    const prevX = pDTX[i];
    const prevY = pDTY[i];
    let bestX = constrain(pUx[i], lastGRect.x, lastGRect.x + lastGRect.w - 1);
    let bestY = constrain(pUy[i], lastGRect.y, lastGRect.y + lastGRect.h - 1);
    let bestW = _bufWeightAtXY(bestX, bestY);

    // Local search (2..6 buffer pixels), prefer less-used pixels to reduce clumping
    if (bestW <= 0) {
      for (let r = 2; r <= 6; r++) {
        for (let k = 0; k < _dirsX.length; k++) {
          const x = constrain(pUx[i] + _dirsX[k] * r, lastGRect.x, lastGRect.x + lastGRect.w - 1);
          const y = constrain(pUy[i] + _dirsY[k] * r, lastGRect.y, lastGRect.y + lastGRect.h - 1);
          let w = _bufWeightAtXY(x, y);
          if (w > 0) {
            const p = (y | 0) * targetW + (x | 0);
            const c = _pickCount[p];
            w = w / (1 + c * 0.65);
          }
          if (w > bestW) {
            bestW = w;
            bestX = x;
            bestY = y;
          }
        }
        if (bestW > 0) break;
      }
    }

    if (bestW > 0) debugLocalFound++;
    else {
      // Keep previous desired target if nothing bright was found near (ux,uy)
      pDTX[i] = prevX;
      pDTY[i] = prevY;
      continue;
    }

    const p = (bestY | 0) * targetW + (bestX | 0);
    _pickCount[p] = min(65535, _pickCount[p] + 1);

    const ix = (bestX + 0.5 - lastGRect.x) * gwInv;
    const iy = (bestY + 0.5 - lastGRect.y) * ghInv;

    pDTX[i] = lastCRect.x + ix * lastCRect.w + pJx[i];
    pDTY[i] = lastCRect.y + iy * lastCRect.h + pJy[i];
  }
}

function computeActSettlingMetrics() {
  const n = min(particleCount, maxParticles);
  const samples = min(200, n);
  const stride = max(1, floor(n / samples));
  let sum = 0;
  let settled = 0;
  const r = 10;
  const r2 = r * r;

  for (let s = 0; s < samples; s++) {
    const i = s * stride;
    const dx = pX[i] - pTX[i];
    const dy = pY[i] - pTY[i];
    const d2 = dx * dx + dy * dy;
    sum += sqrt(d2);
    if (d2 < r2) settled++;
  }

  debugMeanDist = sum / samples;
  debugSettledPct = settled / samples;
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
    pOrder[i] = i;
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
    pJx[i] = (noise(pSeedX[i], pSeedY[i], 0) - 0.5) * 4.0;
    pJy[i] = (noise(pSeedX[i] + 9.1, pSeedY[i] + 3.7, 0) - 0.5) * 4.0;
    pUx[i] = random(targetW);
    pUy[i] = random(targetH);
    pTX[i] = pX[i];
    pTY[i] = pY[i];
    pDTX[i] = pX[i];
    pDTY[i] = pY[i];
  }

  // Stable correspondence across frames: shuffle pOrder once (per (re)init)
  for (let i = n - 1; i > 0; i--) {
    const j = floor(random(i + 1));
    const tmp = pOrder[i];
    pOrder[i] = pOrder[j];
    pOrder[j] = tmp;
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

function startTransition() {
  mode = "TRANSITION";
  transitionStartT = t;
  nextAct = getNextActWithFrames(act);

  // TRANSITION: it's OK to scramble; prime coalesce targets from next act frame 0
  const img = (actFrames[nextAct] || [])[0] || null;
  if (img) {
    loadFrameToBuffer(img);
    assignPersistentUxUyFromCurrentBuffer(true);
    updateDesiredTargetsFromBuffer();
  }
}

function endTransition() {
  act = nextAct;
  actStartT = t;
  mode = "ACT";
  updateActTiming();
  currentFrameIndex = 0;
  lastSampledAct = act;
  lastSampledFrameIndex = -1;
}

function updateAct() {
  if (mode === "ACT") {
    const actElapsed = t - actStartT;
    if (actElapsed >= actDuration) startTransition();
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

  // Controlled motion: strong attraction to target, subtle flow + wobble
  const flow = 0.03 + map(level, 0, 0.2, 0, 0.10, true);
  const swirlBase = 0.40 + map(level, 0, 0.2, 0, 0.55, true);
  const driftX = 0.012;
  const driftY = -0.008;

  const damp = 0.78;
  const maxSpd = 3.6;
  const att = 0.30 + map(level, 0, 0.2, 0, 0.14, true);
  const targetLerp = 0.20;
  const snapRadius = 6;

  const wobbleAmpBase = 0.35 + map(level, 0, 0.2, 0, 1.0, true);
  const wobbleFreq = 0.35;

  for (let i = 0; i < n; i++) {
    // Smooth target updates toward desired target set on frame change
    pTX[i] += (pDTX[i] - pTX[i]) * targetLerp;
    pTY[i] += (pDTY[i] - pTY[i]) * targetLerp;

    const nx = pX[i] * 0.0020;
    const ny = pY[i] * 0.0020;
    const ang = noise(nx, ny, t * 0.28) * TWO_PI * 2.0;

    const dxC = pX[i] - cx;
    const dyC = pY[i] - cy;
    const dC = sqrt(dxC * dxC + dyC * dyC) + 0.001;
    const centerBoost = constrain(1.0 - dC / (min(width, height) * 0.65), 0, 1);
    const swirl = swirlBase * centerBoost;

    pVX[i] += (cos(ang) * flow + driftX + (-dyC / dC) * swirl) * 0.22;
    pVY[i] += (sin(ang) * flow + driftY + (dxC / dC) * swirl) * 0.22;

    // Subtle "living material" wobble around the interpolated target
    const nx1 = noise(pSeedX[i], pSeedY[i], t * wobbleFreq);
    const ny1 = noise(pSeedX[i] + 17.3, pSeedY[i] + 91.7, t * wobbleFreq);
    const dx0 = pTX[i] - pX[i];
    const dy0 = pTY[i] - pY[i];
    const d0 = sqrt(dx0 * dx0 + dy0 * dy0);
    const wobbleAmp = min(1.5, wobbleAmpBase) * constrain(d0 / 12, 0, 1);
    const wx = (nx1 - 0.5) * 2.0 * wobbleAmp;
    const wy = (ny1 - 0.5) * 2.0 * wobbleAmp;

    const desiredX = pTX[i] + wx;
    const desiredY = pTY[i] + wy;
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

    // Snap near target for crisp silhouette
    const d2 = dxT * dxT + dyT * dyT;
    if (d2 < snapRadius * snapRadius) {
      pX[i] = lerp(pX[i], pTX[i], 0.35);
      pY[i] = lerp(pY[i], pTY[i], 0.35);
      pVX[i] *= 0.6;
      pVY[i] *= 0.6;
    }

    // Containment
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

  const att = stageA ? 0.0 : lerp(0.04, 0.24, stageT);
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

    if (!stageA) {
      // Coalesce into desired targets computed once at transition start
      const targetLerp = 0.25 + stageT * 0.15;
      pTX[i] += (pDTX[i] - pTX[i]) * targetLerp;
      pTY[i] += (pDTY[i] - pTY[i]) * targetLerp;

      const nx1 = noise(pSeedX[i], pSeedY[i], t * wobbleFreq);
      const ny1 = noise(pSeedX[i] + 17.3, pSeedY[i] + 91.7, t * wobbleFreq);
      const wx = (nx1 - 0.5) * 2.0 * wobbleAmp;
      const wy = (ny1 - 0.5) * 2.0 * wobbleAmp;

      const dxT = (pTX[i] + wx) - pX[i];
      const dyT = (pTY[i] + wy) - pY[i];
      pVX[i] += dxT * att;
      pVY[i] += dyT * att;

      if (snap > 0) {
        pX[i] = lerp(pX[i], pTX[i], snap * 0.35);
        pY[i] = lerp(pY[i], pTY[i], snap * 0.35);
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
    currentFrameIndex = constrain(floor(actElapsed * framesPerSecond), 0, n - 1);

    const img = frames[currentFrameIndex];
    if (!img) return;

    if (act !== lastSampledAct || currentFrameIndex !== lastSampledFrameIndex) {
      loadFrameToBuffer(img);
      const isActStart = act !== lastSampledAct || lastSampledFrameIndex === -1;
      if (isActStart) assignPersistentUxUyFromCurrentBuffer(false);
      updateDesiredTargetsFromBuffer();

      if (isActStart) {
        const nn = min(particleCount, maxParticles);
        for (let i = 0; i < nn; i++) {
          pTX[i] = pDTX[i];
          pTY[i] = pDTY[i];
        }
      }

      lastSampledAct = act;
      lastSampledFrameIndex = currentFrameIndex;
    }

    updateActParticles(level);
    computeActSettlingMetrics();
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
  rect(10, 10, 330, 122, 8);
  fill(255);
  textAlign(LEFT, TOP);
  textSize(12);
  text(`mode: ${mode}`, 18, 16);
  text(`act: ${act}  next: ${mode === "TRANSITION" ? nextAct : getNextActWithFrames(act)}`, 18, 32);
  text(`frameIndex: ${currentFrameIndex}`, 18, 48);
  text(`meanDist: ${debugMeanDist.toFixed(2)}`, 18, 64);
  text(`settledPct: ${(debugSettledPct * 100).toFixed(1)}%`, 18, 80);
  text(`localFound: ${debugLocalFound}`, 18, 96);
  text(`particles: ${particleCount}  imagesDrawn: ${imagesDrawnToCanvas}`, 18, 112);
  pop();
}
