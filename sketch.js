function preload() {
  try {
    Typography.preload();
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
  }
  // Loading screen image (revealed by particles during boot)
  try {
    window.__loadingImg = loadImage("assets/loading/loadingscreen.jpg");
  } catch (_) {
    window.__loadingImg = null;
  }
}

// Keep last valid ACT silhouette when sampler is catching up (prevents "break" on claps).
let _actHasDesired = false;
let _actDesiredFor = 0;
const _prefOff = { off: 0 };
let _debugDtSec = 0;
let _debugRate = 0;
let _soundOK = false;
let compareMode = false;
let compareAct = 1;
let compareSrc0 = 0;
let compareAlpha = 0.5;
let _mainCanvasElt = null;
let _cssTick = 0;
let _warmupI = 0;
const WARMUP_ACT1_FRAMES = 24;
const WARMUP_OTHER_FRAMES = 2;
const PREFETCH_AHEAD = 12;
let _prefetchK = 0;

// Startup loading screen: precompute sampler frames so the first run is smooth.
// Per request: load ALL act frames (counts cache) before allowing the sketch to start.
let booting = true;
let _bootActs = [];
let _bootTotal = 0;
let _bootActI = 0;
let _bootIdx = 0;
let _bootReqTick = 0;
const BOOT_REQ_PER_FRAME = 6;   // request this many unique frames per draw (I/O), compute remains capped
const BOOT_SCAN_LIMIT = 600;    // max skips when searching for an unrequested frame
const BOOT_QUEUE_MAX = 220;     // avoid unbounded queue growth
let _loadingCounts = null;
let _loadingCountsOK = false;
let _bootShapeApplied = false;
let _bootSeedsInited = false;
function _ensureMainCanvasCss() {
  try {
    if (!_mainCanvasElt) _mainCanvasElt = document.getElementById ? document.getElementById("mainCanvas") : null;
    if (!_mainCanvasElt) return;
    _mainCanvasElt.style.position = "fixed";
    _mainCanvasElt.style.left = "0";
    _mainCanvasElt.style.top = "0";
    _mainCanvasElt.style.width = "100vw";
    _mainCanvasElt.style.height = "100vh";
    _mainCanvasElt.style.display = "block";
  } catch (_) {}
}
function _initLoadingCountsIfNeeded() {
  if (_loadingCountsOK) return;
  try {
    const img = window.__loadingImg;
    if (!img || !img.width) return;

    // Build the loading silhouette at higher sampling resolution so thin top/bottom text doesn't vanish.
    const maxDim = 720;
    let gW = 560;
    let gH = max(64, round((img.height / max(1, img.width)) * gW));
    if (max(gW, gH) > maxDim) {
      const s = maxDim / max(gW, gH);
      gW = max(64, round(gW * s));
      gH = max(64, round(gH * s));
    }
    const g = createGraphics(gW, gH);
    g.pixelDensity(1);
    g.noSmooth();
    g.clear();
    g.background(0);
    g.imageMode(CORNER);
    try { if (g.drawingContext) g.drawingContext.imageSmoothingEnabled = false; } catch (_) {}

    // Draw stretched so the loading image always maps to the full grid (no letterboxing / perceived crop).
    g.image(img, 0, 0, gW, gH);
    g.loadPixels();

    // Detect "white background + dark ink" and invert so the ink becomes the silhouette.
    let meanBr = 0;
    const sampleN = 96;
    for (let s = 0; s < sampleN; s++) {
      const x = (s * 37) % gW;
      const y = (s * 91) % gH;
      const p = 4 * (y * gW + x);
      const br = (g.pixels[p] + g.pixels[p + 1] + g.pixels[p + 2]) / 765;
      meanBr += br;
    }
    meanBr /= sampleN;
    const invert = meanBr > 0.62;

    const weights = new Float32Array(CELLS);
    const pixC = new Uint16Array(CELLS);
    const out = new Uint16Array(CELLS);

    // Map full buffer -> full grid (stretch). Deterministic and avoids "top cut" on thin text.
    for (let y = 0; y < gH; y++) {
      const rr = ((y + 0.5) * ROWS / gH) | 0;
      if (rr < 0 || rr >= ROWS) continue;
      const row = y * gW;
      for (let x = 0; x < gW; x++) {
        const cc = ((x + 0.5) * COLS / gW) | 0;
        if (cc < 0 || cc >= COLS) continue;
        const cell = rr * COLS + cc;
        const p = 4 * (row + x);
        const a = g.pixels[p + 3] / 255;
        if (a <= 0.02) continue;
        const br0 = (g.pixels[p] + g.pixels[p + 1] + g.pixels[p + 2]) / 765;
        const br = invert ? (1 - br0) : br0;
        if (br <= 0.02) continue;
        weights[cell] += br * a;
        if (pixC[cell] < 65535) pixC[cell]++;
      }
    }

    let maxBr = 0;
    for (let i = 0; i < CELLS; i++) {
      const pc = pixC[i];
      if (pc === 0) { weights[i] = 0; continue; }
      const br = weights[i] / pc;
      weights[i] = br;
      if (br > maxBr) maxBr = br;
    }
    if (maxBr <= 1e-6) {
      out.fill(0);
      out[((ROWS >> 1) * COLS + (COLS >> 1)) | 0] = N;
    } else {
      let total = 0;
      // Higher contrast for the loader so the image reads clearly.
      // (Keep some lift so thin strokes don't disappear.)
      const thr = invert ? 0.06 : 0.04;
      const gamma = 1.65;
      const lift = 0.02;
      const cdf = new Float32Array(CELLS);
      for (let i = 0; i < CELLS; i++) {
        const brN = weights[i] / maxBr;
        const x = (brN - thr) / max(1e-6, 1 - thr);
        if (x <= 0) { weights[i] = 0; continue; }
        const w = pow(x, gamma) + lift * x;
        weights[i] = w;
        total += w;
        cdf[i] = total;
      }
      out.fill(0);
      if (total <= 1e-9) {
        out[((ROWS >> 1) * COLS + (COLS >> 1)) | 0] = N;
      } else {
        const step = total / N;
        let u = 0.37 * step;
        let i = 0;
        for (let k = 0; k < N; k++) {
          const th = u + k * step;
          while (i < CELLS - 1 && cdf[i] < th) i++;
          out[i]++;
        }
      }
    }

    _loadingCounts = out;
    _loadingCountsOK = true;
  } catch (_) {}
}
function _applyLoadingShapeToParticles() {
  if (_bootShapeApplied) return;
  if (!_loadingCountsOK || !_loadingCounts) return;
  try {
    Particles.setDesired(_loadingCounts, 0, 0, 0);
    Particles.setWalkFromCounts(_loadingCounts, 0, 0);
    Particles.softRetargetToDesired();
    _bootShapeApplied = true;
  } catch (_) {}
}
function drawLoadingScreen(pct, ready, total) {
  // Pretty boot: particle-reveal of the loading image (organic build).
  // (We keep the old text-based loader below as fallback; early-return from the pretty path.)
  const p = constrain(pct, 0, 1);

  // Ensure we have a target silhouette to reveal.
  _initLoadingCountsIfNeeded();
  _applyLoadingShapeToParticles();

  background(0);
  sceneVis = 1;
  sceneA = 1;

  try {
    Render.tintA = 0;
    Render.tintB = 0;
    Render.clipOn = false;
    Render.revealOn = true;
    Render.revealP = p;
    Render.revealFreq = 0.0028;
    Render.revealSeed = 42.0;
    if (!_bootSeedsInited) {
      _bootSeedsInited = true;
      // 3-seed blob growth (triangle) so it doesn't read as left->right.
      Render.revealSx0 = width * 0.50;
      Render.revealSy0 = height * 0.52;
      Render.revealSx1 = width * 0.36;
      Render.revealSy1 = height * 0.40;
      Render.revealSx2 = width * 0.64;
      Render.revealSy2 = height * 0.68;
    }
  } catch (_) {}

  // Drive particles toward the loading silhouette during boot even though `t` is not advancing.
  // We do NOT advance `t`; we only provide a small dt for the integrator.
  tAdvanced = true;
  tDelta = 0.02;
  try {
    computeDeficitHotspots();
    Particles.rebalancePasses(6);
    Particles.updateInsideCells(0, 0.10, 0.12, 0);
    Particles.separate(0.018, 3.2);
  } catch (_) {}

  try { Particles.draw(); } catch (_) {}
  try {
    Render.revealOn = false;
    Render.clipOn = false;
  } catch (_) {}

  push();
  fill(255, 220);
  noStroke();
  textAlign(CENTER, TOP);
  textSize(14);
  text(`Loading... ${Math.round(p * 100)}%`, width / 2, 14);
  pop();
  return;

  background(0);
  fill(255);
  noStroke();
  textAlign(CENTER, CENTER);
  textSize(16);
  const msg = `Loading... ${Math.round(pct * 100)}%  (${ready}/${total})`;
  text(msg, width / 2, height / 2);
  textSize(12);
  text("Preparing all frames for smooth playback", width / 2, height / 2 + 22);
  try {
    const q = Sampler && Sampler.queueLen ? Sampler.queueLen() : 0;
    const cq = Sampler && Sampler.computeQueueLen ? Sampler.computeQueueLen() : 0;
    const inf = Sampler && Sampler.inFlight ? Sampler.inFlight : 0;
    text(`q:${q}  compute:${cq}  inFlight:${inf}`, width / 2, height / 2 + 40);
  } catch (_) {}
}
function _bootInitIfNeeded() {
  if (_bootActs.length > 0) return;
  _bootActs = [];
  _bootTotal = 0;
  for (let a of ACTS) {
    const c = SRC_COUNT[a] || 0;
    if (c > 0) { _bootActs.push(a); _bootTotal += c; }
  }
  _bootActI = 0;
  _bootIdx = 0;
  _bootReqTick = 0;
  _initLoadingCountsIfNeeded();
  _bootShapeApplied = false;
  _applyLoadingShapeToParticles();
  _bootSeedsInited = false;
}
function _bootStepAll() {
  if (!booting) return true;
  if (compareMode) { booting = false; return true; }
  _bootInitIfNeeded();
  if (_bootActs.length === 0 || _bootTotal <= 0) { booting = false; return true; }

  // Request frames steadily (I/O side). Keep queue bounded.
  try {
    const q = Sampler && Sampler.queueLen ? Sampler.queueLen() : 0;
    const cq = Sampler && Sampler.computeQueueLen ? Sampler.computeQueueLen() : 0;
    if (q + cq < BOOT_QUEUE_MAX) {
      let req = 0;
      let scans = 0;
      while (req < BOOT_REQ_PER_FRAME && scans < BOOT_SCAN_LIMIT) {
        const a = _bootActs[_bootActI] | 0;
        const cycle = SRC_COUNT[a] || 0;
        if (cycle > 0) {
          // Skip frames that are already done/loading/failed.
          const c = Sampler.cache ? Sampler.cache[a] : null;
          if (!c || c.cycle !== cycle) Sampler.ensureActCache(a, cycle);
          const cc = Sampler.cache ? Sampler.cache[a] : null;
          const idx = _bootIdx | 0;
          const ok = cc && (cc.done[idx] || cc.loading[idx] || cc.failed[idx]);
          if (!ok) {
            Sampler.ensure(a, idx, cycle, _prefOff);
            req++;
          }
        }
        // Advance cursor
        _bootIdx++;
        if (_bootIdx >= (cycle || 1)) {
          _bootIdx = 0;
          _bootActI++;
          if (_bootActI >= _bootActs.length) _bootActI = 0;
        }
        scans++;
      }
    }
  } catch (_) {}

  // Progress: done or failed counts across all acts.
  let ready = 0;
  for (let a of _bootActs) ready += (Sampler.readyOrFailed ? Sampler.readyOrFailed(a) : (Sampler.ready ? Sampler.ready(a) : 0));
  const pct = _bootTotal > 0 ? ready / _bootTotal : 1;

  const typOK = (typeof Typography !== "undefined" && Typography && Typography.ready);
  const q = Sampler && Sampler.queueLen ? Sampler.queueLen() : 0;
  const cq = Sampler && Sampler.computeQueueLen ? Sampler.computeQueueLen() : 0;
  const inf = Sampler && Sampler.inFlight ? Sampler.inFlight : 0;

  // Consider boot complete only when everything is computed (or failed) and queues are drained.
  if (ready >= _bootTotal && q === 0 && cq === 0 && inf === 0 && typOK) {
    booting = false;
    try {
      Render.clipOn = false;
      Render.revealOn = false;
    } catch (_) {}
    // Restore letterboxed simulation grid for the real sketch.
    try {
      resizeGrid();
      // `resizeGrid()` reallocates `Particles.cellCounts`, so rebuild counts from current positions
      // or redistribution will never happen and we'll stay stuck in the loading silhouette.
      if (typeof Particles !== "undefined" && Particles && typeof Particles.remapCellsFromPositions === "function") {
        Particles.remapCellsFromPositions();
      }
    } catch (_) {}
    return true;
  }

  drawLoadingScreen(pct, ready, _bootTotal);
  return false;
}
function _applyUrlParams() {
  try {
    if (typeof location === "undefined") return;
    const qs = new URLSearchParams(location.search || "");
    const cmp = qs.get("compare") || qs.get("cmp");
    if (cmp === "1" || cmp === "true") {
      compareMode = true;
      audioStarted = true;
    }
    const actQ = qs.get("act");
    if (actQ != null) compareAct = constrain(parseInt(actQ, 10) || 1, 1, 4);
    const srcQ = qs.get("src") || qs.get("src0");
    if (srcQ != null) compareSrc0 = parseInt(srcQ, 10) || 0;
    const aQ = qs.get("alpha") || qs.get("a");
    if (aQ != null) compareAlpha = constrain(parseFloat(aQ) || 0, 0, 1);
    const dbgQ = qs.get("debug") || qs.get("d");
    if (dbgQ === "1" || dbgQ === "true") debugOn = true;
  } catch (_) {}
}

function setup() {
  // Use a physical-pixel canvas size so the simulation/rendering look is consistent even when
  // Chrome uses different devicePixelRatio/zoom per-site (localhost vs GitHub Pages).
  pixelDensity(1);
  const cssW = (typeof window !== "undefined" ? (window.innerWidth || windowWidth) : windowWidth) | 0;
  const cssH = (typeof window !== "undefined" ? (window.innerHeight || windowHeight) : windowHeight) | 0;
  const dpr = (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1) || 1;
  const physW = max(1, round(cssW * dpr));
  const physH = max(1, round(cssH * dpr));
  const cnv = createCanvas(physW, physH);
  try { _mainCanvasElt = cnv?.elt || null; } catch (_) {}
  // Hard-set DOM id + CSS sizing here (more reliable than relying on p5.Element helpers across builds/hosts).
  try {
    if (_mainCanvasElt) {
      _mainCanvasElt.id = "mainCanvas";
      _ensureMainCanvasCss();
    }
  } catch (_) {}
  applySimSeed();
  _applyUrlParams();
  try { cnv?.style?.("display", "block"); } catch (_) {}
  try {
    _soundOK = (typeof p5 !== "undefined" && typeof p5.AudioIn === "function" && typeof p5.Amplitude === "function");
    if (_soundOK) {
      mic = new p5.AudioIn();
      amp = new p5.Amplitude();
    } else {
      // Fallback: allow Auto mode + visuals even if p5.sound failed to load on the host.
      mic = { start: () => {}, stop: () => {} };
      amp = { getLevel: () => 0, setInput: () => {} };
      console.warn("[Audio] p5.sound not available; microphone disabled. Press 'A' for auto.");
    }

    Style.init();
    updatePageZoom();
    resizeGrid();
    // During boot, use a full-bleed grid so the loading image isn't "cut" on tall viewports.
    if (booting) {
      gridX0 = 0;
      gridY0 = 0;
      cellW = width / COLS;
      cellH = height / ROWS;
      invCellW = 1 / max(1e-6, cellW);
      invCellH = 1 / max(1e-6, cellH);
    }
    Sampler.init();
    Particles.resizeBins();
    Particles.init();
    _initLoadingCountsIfNeeded();
    _bootShapeApplied = false;
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
    return;
  }

  try {
    Typography.init();
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
  }
}

function windowResized() {
  pixelDensity(1);
  const cssW = (typeof window !== "undefined" ? (window.innerWidth || windowWidth) : windowWidth) | 0;
  const cssH = (typeof window !== "undefined" ? (window.innerHeight || windowHeight) : windowHeight) | 0;
  const dpr = (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1) || 1;
  const physW = max(1, round(cssW * dpr));
  const physH = max(1, round(cssH * dpr));
  resizeCanvas(physW, physH);
  try {
    // Keep the main canvas visually full-screen in CSS pixels (robust even if other canvases exist).
    if (!_mainCanvasElt) _mainCanvasElt = document.getElementById ? document.getElementById("mainCanvas") : null;
    _ensureMainCanvasCss();
  } catch (_) {}
  applySimSeed();
  try {
    Style.init();
    updatePageZoom();
    resizeGrid();
    if (booting) {
      gridX0 = 0;
      gridY0 = 0;
      cellW = width / COLS;
      cellH = height / ROWS;
      invCellW = 1 / max(1e-6, cellW);
      invCellH = 1 / max(1e-6, cellH);
    }
    Particles.resizeBins();
    Particles.init();
    Typography.resize();
    _bootShapeApplied = false;
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
  }
}

function draw() {
  // Some hosts/extensions can mutate the canvas element styles; periodically re-assert full-screen CSS.
  if (((_cssTick++) & 31) === 0) _ensureMainCanvasCss();
  // Cap sampler compute work per frame; during boot allow more so loading finishes sooner.
  try { if (Sampler && Sampler.stepCompute) Sampler.stepCompute(booting ? 3 : 1); } catch (_) {}
  if (window.__fatalError) {
    background(0);
    fill(255);
    noStroke();
    textAlign(LEFT, TOP);
    textSize(12);
    text(String(window.__fatalError), 12, 12, width - 24, height - 24);
    noLoop();
    return;
  }

  // Always finish the boot loader first (precompute all frames), even if Auto is toggled.
  if (booting) {
    try {
      if (!_bootStepAll()) return;
    } catch (_) {}
  }

  // If neither Auto nor Mic is running, keep showing the start screen (even if `audioStarted`
  // was previously true due to Auto mode). This prevents a "dead" state where nothing renders.
  if (!autoRun && !micRunning) {
    background(0);
    drawStartScreen();
    return;
  }

  const level = autoRun ? 0.09 : (micRunning ? amp.getLevel() : 0);

  // Real-time scaled step (keeps loop length stable across devices/FPS).
  // NOTE: still freezes when `t` is not advanced (silence -> no motion).
  // Keep time stable across FPS, but avoid very large jumps that can overwhelm sampling/redistribution.
  // 0.10 means: down to ~10fps, timing stays correct; below that it will slow rather than spike.
  const dtSec = constrain((typeof deltaTime === "number" ? deltaTime : 16.7) / 1000, 0, 0.10);
  _debugDtSec = dtSec;

  if (!compareMode && autoRun) {
    t += AUTO_SPEED * dtSec;
    _debugRate = AUTO_SPEED;
  } else if (!compareMode && micRunning && level > MIC_THRESHOLD) {
    // Mic speed in sound-seconds per real second.
    // Quickly reaches 1.0 so "always loud" finishes a full loop in ~60s like Auto.
    const rate = map(level, MIC_THRESHOLD, MIC_THRESHOLD + 0.02, 0.35, 1.0, true);
    t += rate * dtSec;
    _debugRate = rate;
  } else {
    _debugRate = 0;
  }

  tDelta = t - _prevT;
  tAdvanced = tDelta !== 0;
  _prevT = t;
  if (compareMode) {
    // Drive motion even while time is frozen for deterministic cross-host comparisons.
    tAdvanced = true;
    tDelta = 0.016;
  }

  // Fade scene visibility based on mic level (objects/text appear only with noise).
  // NOT tied to `t`, so it can fade out when `t` is frozen.
  const targetVis = autoRun
    ? 1
    : micRunning
      ? pow(constrain(map(level, VIS_LEVEL_ON, VIS_LEVEL_FULL, 0, 1, true), 0, 1), 0.65)
      : 0;
  const k = targetVis > sceneVis ? VIS_IN : VIS_OUT;
  sceneVis += (targetVis - sceneVis) * k;
  sceneVis = constrain(sceneVis, 0, 1);
  sceneA = sceneVis * sceneVis;

  try {
    Sampler.resetFrameStats();
    if (compareMode) {
      Acts.mode = "ACT";
      Acts.act = constrain(compareAct | 0, 1, 4);
      Acts.next = nextActWithFrames(Acts.act);
      Acts.cycle = SRC_COUNT[Acts.act] || 0;
      Acts.fps = FPS_EFFECTIVE[Acts.act] || 24;
      Acts.elapsed = 0;
      Acts.dur = Acts.cycle > 0 && Acts.fps > 0 ? Acts.cycle / Acts.fps : 0;
      if (Acts.cycle > 0) {
        Acts.src0 = ((compareSrc0 | 0) % Acts.cycle + Acts.cycle) % Acts.cycle;
        Acts.src1 = (Acts.src0 + 1) % Acts.cycle;
        Acts.alpha = constrain(compareAlpha, 0, 1);
        Acts.cycles = 0;
      } else {
        Acts.src0 = Acts.src1 = 0;
        Acts.alpha = 0;
        Acts.cycles = 0;
      }
    } else {
      Acts.update();
    }

    // Streaming prefetch: keep a small window of upcoming frames requested so loads/compute happen
    // ahead of time (prevents "first loop" hitching when new frames are needed).
    try {
      if (Acts.mode === "ACT") {
        const a = Acts.act;
        const cycle = SRC_COUNT[a] || 0;
        if (cycle > 0) {
          const step = max(1, (FRAME_STEP && FRAME_STEP[a]) | 0);
          const k = ((_prefetchK++ % PREFETCH_AHEAD) + 2) | 0;
          const idx = (Acts.src0 + k * step) % cycle;
          Sampler.ensure(a, idx, cycle, _prefOff);
        }
      }
    } catch (_) {}
    Style.update(Acts, level);
    background(Render.bg);

    let desiredSum = 0;
    let moved = 0;
    let mismatch = 0;

    if (Acts.mode === "ACT") {
      if (Acts.cycle > 0) {
        const ok0 = Sampler.ensure(Acts.act, Acts.src0, Acts.cycle, _off0);
        const ok1 = Sampler.ensure(Acts.act, Acts.src1, Acts.cycle, _off1);

        // Prefetch a small window ahead so bursty inputs (claps) don't cause missing-frame fallbacks.
        if (tAdvanced && !compareMode) {
          Sampler.ensure(Acts.act, (Acts.src0 + 1) % Acts.cycle, Acts.cycle, _prefOff);
          Sampler.ensure(Acts.act, (Acts.src0 + 2) % Acts.cycle, Acts.cycle, _prefOff);
          Sampler.ensure(Acts.act, (Acts.src1 + 1) % Acts.cycle, Acts.cycle, _prefOff);
        }

        if (ok0 || ok1) {
          const counts = Sampler.counts(Acts.act);
          if (ok0 && ok1) {
            Particles.setDesired(counts, _off0.off, _off1.off, Acts.alpha);
            Particles.setWalkFromCounts(counts, _off0.off, _off1.off);
          } else if (ok0) {
            Particles.setDesired(counts, _off0.off, _off0.off, 0);
            Particles.setWalkFromCounts(counts, _off0.off, _off0.off);
          } else {
            Particles.setDesired(counts, _off1.off, _off1.off, 0);
            Particles.setWalkFromCounts(counts, _off1.off, _off1.off);
          }
          _actHasDesired = true;
          _actDesiredFor = Acts.act;

          for (let i = 0; i < CELLS; i++) desiredSum += Particles.desired[i];

          let maxD = 1;
          for (let i = 0; i < CELLS; i++) { const d = Particles.desired[i] | 0; if (d > maxD) maxD = d; }
          Render.maxD = maxD;

          if (tAdvanced) {
            // If sound-time advances in bursts (claps), do a few small sub-steps so the shape has time
            // to redistribute and doesn't "break" while the sampler catches up.
            const td0 = tDelta;
            const steps = constrain(ceil(abs(td0) / 0.02), 1, 4);
            const td = td0 / steps;
            for (let s = 0; s < steps; s++) {
              tDelta = td;
              mismatch = estimateMismatch();
              debugMeanCellDist = estimateMeanCellDist();

              const posCU = constrain((debugMeanCellDist - cellW * 0.85) / (cellW * 3.2), 0, 1);
              const targetCU = constrain((mismatch - N * 0.10) / (N * 0.45), 0, 1);
              catchUp = max(catchUp * 0.94, targetCU, posCU * 0.35);

              const passesTotal =
                (mismatch > 900 ? 10 : mismatch > 600 ? 8 : mismatch > 350 ? 6 : mismatch > 180 ? 4 : 3) + floor(catchUp * 2);
              const passes = min(8, max(2, ceil(passesTotal / steps)));

              computeDeficitHotspots();
              moved += Particles.rebalancePasses(passes);
              if (catchUp > 0.32) catchUpNudge(catchUp);

              Particles.updateInsideCells(level, 0.18 + catchUp * 0.26, 0.15, catchUp * 0.7);
              Particles.separate(0.018, 3.2);
              catchUp *= 0.97;
            }
            tDelta = td0;
          }
        } else {
          // If we already have a valid silhouette for this act, keep it until the sampler catches up.
          if (_actHasDesired && _actDesiredFor === Acts.act) {
            if (tAdvanced) {
              const td0 = tDelta;
              const steps = constrain(ceil(abs(td0) / 0.02), 1, 4);
              const td = td0 / steps;
              for (let s = 0; s < steps; s++) {
                tDelta = td;
                mismatch = estimateMismatch();
                debugMeanCellDist = estimateMeanCellDist();
                const posCU = constrain((debugMeanCellDist - cellW * 0.85) / (cellW * 3.2), 0, 1);
                const targetCU = constrain((mismatch - N * 0.10) / (N * 0.45), 0, 1);
                catchUp = max(catchUp * 0.94, targetCU, posCU * 0.35);
                computeDeficitHotspots();
                moved += Particles.rebalancePasses(3);
                Particles.updateInsideCells(level, 0.18 + catchUp * 0.26, 0.15, catchUp * 0.7);
                Particles.separate(0.018, 3.2);
                catchUp *= 0.97;
              }
              tDelta = td0;
            }
          } else {
            if (tAdvanced) Particles.dance(level, 0.20);
          }
        }
      } else {
        if (tAdvanced) Particles.dance(level, 0.35);
      }

      Particles.draw();
    } else {
      const p = constrain((t - Acts.transitionStartT) / max(1e-6, TRANSITION_DURATION), 0, 1);
      const stageA = p < DANCE_PORTION;
      const s = stageA ? p / max(1e-6, DANCE_PORTION) : (p - DANCE_PORTION) / max(1e-6, 1 - DANCE_PORTION);

      let maxD = 1;
      for (let i = 0; i < CELLS; i++) { const d = Particles.desired[i] | 0; if (d > maxD) maxD = d; }
      Render.maxD = maxD;

      if (tAdvanced) {
        if (stageA) {
          if (!Acts.prepDone) {
            const a = Acts.next;
            const cycle = SRC_COUNT[a] || 0;
            if (cycle > 0) {
              const step = max(1, (FRAME_STEP && FRAME_STEP[a]) | 0);
              const ok0 = Sampler.ensure(a, 0, cycle, _off0);
              const ok1 = Sampler.ensure(a, (0 + step) % cycle, cycle, _off1);
              if (ok0 && ok1) Acts.prepDone = true;
            }
          }
          Particles.dance(level, s);
        } else {
          const a = Acts.next;
          const cycle = SRC_COUNT[a] || 0;
          if (cycle > 0) {
            const step = max(1, (FRAME_STEP && FRAME_STEP[a]) | 0);
            const ok0 = Sampler.ensure(a, 0, cycle, _off0);
            const ok1 = Sampler.ensure(a, (0 + step) % cycle, cycle, _off1);
            if (ok0 || ok1) {
              const counts = Sampler.counts(a);
              if (ok0 && ok1) {
                Particles.setDesired(counts, _off0.off, _off1.off, s);
                Particles.setWalkFromCounts(counts, _off0.off, _off1.off);
              } else if (ok0) {
                Particles.setDesired(counts, _off0.off, _off0.off, 0);
                Particles.setWalkFromCounts(counts, _off0.off, _off0.off);
              } else {
                Particles.setDesired(counts, _off1.off, _off1.off, 0);
                Particles.setWalkFromCounts(counts, _off1.off, _off1.off);
              }
              for (let i = 0; i < CELLS; i++) desiredSum += Particles.desired[i];

              debugMeanCellDist = estimateMeanCellDist();
              computeDeficitHotspots();
              moved += Particles.rebalancePasses(4);
              if (p > (DANCE_PORTION + 0.10)) catchUpNudge(0.12);

              const snap01 = map(s, 0.70, 1.0, 0, 0.25, true);
              Particles.updateInsideCells(level, snap01, 0.0, 0.18);
              Particles.separate(0.020, 3.2);
            } else {
              Particles.dance(level, 1.0);
            }
          } else {
            Particles.dance(level, 1.0);
          }
        }
      }

      Particles.draw();
    }

    if (debugOn) drawDebug(level, desiredSum, moved, mismatch);
    if (showGrid) drawGridOverlay();
    Typography.draw(level);
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
  }
}

function drawStartScreen() {
  fill(255);
  noStroke();
  textAlign(CENTER, CENTER);
  textSize(18);
  if (_soundOK) {
    text("Tap to activate microphone\nSound drives the animation", width / 2, height / 2);
  } else {
    text("Audio unavailable (p5.sound failed to load)\nPress 'A' for Auto mode", width / 2, height / 2);
  }
}

function mousePressed() {
  // Allow starting the mic at any time (even during boot / even if Auto was toggled),
  // so "sound mode" is always available on mobile (no keyboard).
  if (compareMode) return;
  if (micRunning) return;
  if (!_soundOK || !mic || typeof mic.start !== "function") {
    // Don't crash on hosts where p5.sound didn't load; keep the sketch running.
    console.warn("[Audio] Mic start requested but p5.sound/mic is unavailable. Press 'A' for auto.");
    return;
  }
  if (typeof userStartAudio === "function") userStartAudio();
  mic.start(() => {
    amp.setInput(mic);
    audioStarted = true;
    micRunning = true;
    // Prefer mic-driven time once the mic is running (prevents "only Auto works").
    autoRun = false;
    Acts.actStartT = t;
  });
}

function touchStarted() {
  // Mobile browsers require a user gesture for audio; treat touch like mouse.
  try { mousePressed(); } catch (_) {}
  return false;
}

function keyPressed() {
  if (key === "d" || key === "D") debugOn = !debugOn;
  if (key === "g" || key === "G") showGrid = !showGrid;
  if (key === "a" || key === "A") {
    autoRun = !autoRun;
    if (autoRun) {
      audioStarted = true;
    } else {
      if (!micRunning) audioStarted = false;
    }
  }
  if (key === "c" || key === "C") {
    compareMode = !compareMode;
    audioStarted = true;
    // Helpful when sharing: reflect current compare state in the URL (no reload).
    try {
      if (typeof history !== "undefined" && typeof location !== "undefined") {
        const qs = new URLSearchParams(location.search || "");
        if (compareMode) qs.set("compare", "1");
        else qs.delete("compare");
        qs.set("act", String(compareAct));
        qs.set("src", String(compareSrc0));
        qs.set("alpha", String(compareAlpha.toFixed(3)));
        const url = location.pathname + (qs.toString() ? "?" + qs.toString() : "");
        history.replaceState(null, "", url);
      }
    } catch (_) {}
  }
  if (compareMode) {
    if (key === "1") compareAct = 1;
    if (key === "2") compareAct = 2;
    if (key === "3") compareAct = 3;
    if (key === "4") compareAct = 4;
    if (keyCode === LEFT_ARROW) compareSrc0 -= (keyIsDown(SHIFT) ? 10 : 1);
    if (keyCode === RIGHT_ARROW) compareSrc0 += (keyIsDown(SHIFT) ? 10 : 1);
    if (keyCode === UP_ARROW) compareAlpha = constrain(compareAlpha + (keyIsDown(SHIFT) ? 0.10 : 0.02), 0, 1);
    if (keyCode === DOWN_ARROW) compareAlpha = constrain(compareAlpha - (keyIsDown(SHIFT) ? 0.10 : 0.02), 0, 1);
    // Update URL for easy sharing/debugging (no reload).
    try {
      if (typeof history !== "undefined" && typeof location !== "undefined") {
        const qs = new URLSearchParams(location.search || "");
        qs.set("compare", "1");
        qs.set("act", String(compareAct));
        qs.set("src", String(compareSrc0));
        qs.set("alpha", String(compareAlpha.toFixed(3)));
        const url = location.pathname + "?" + qs.toString();
        history.replaceState(null, "", url);
      }
    } catch (_) {}
  }
  if (key === "p" || key === "P") {
    const dpr = (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1) || 1;
    const z = (typeof pageZoom === "number" ? pageZoom : 1) || 1;
    const pd = (typeof pixelDensity === "function" ? pixelDensity() : 1) || 1;
    const snap = {
      href: (typeof location !== "undefined" ? location.href : ""),
      dpr: +dpr.toFixed(3),
      zoom: +z.toFixed(3),
      // Canvas logical coords are in physical pixels (we size the canvas to innerWidth*dpr and set pixelDensity(1)).
      w: width, h: height,
      cssW: (typeof window !== "undefined" ? (window.innerWidth || 0) : 0),
      cssH: (typeof window !== "undefined" ? (window.innerHeight || 0) : 0),
      // Actual backing-store size (what the canvas really renders at in pixels)
      bufW: Math.round(width * pd),
      bufH: Math.round(height * pd),
      grid: `${COLS}x${ROWS}`,
      cell: +cellW.toFixed(3),
      soundOK: !!_soundOK,
      micRunning: !!micRunning,
      autoRun: !!autoRun,
      compareMode: !!compareMode,
      mode: Acts.mode,
      act: Acts.act,
      src0: Acts.src0,
      src1: Acts.src1,
      a: +Acts.alpha.toFixed(3),
      ready: (Sampler.ready ? Sampler.ready(Acts.act) : 0),
      hits: Sampler.hits, misses: Sampler.misses,
      inFlight: Sampler.inFlight || 0,
      q: Sampler.queueLen ? Sampler.queueLen() : 0,
    };
    try {
      if (typeof document !== "undefined") {
        const cs = document.querySelectorAll ? document.querySelectorAll("canvas") : null;
        snap.canvasCount = cs ? cs.length : 0;
        const main = document.getElementById ? document.getElementById("mainCanvas") : null;
        if (main && main.getBoundingClientRect) {
          const r = main.getBoundingClientRect();
          snap.mainRect = { x: +r.x.toFixed(1), y: +r.y.toFixed(1), w: +r.width.toFixed(1), h: +r.height.toFixed(1) };
        }
      }
    } catch (_) {}
    try {
      if (Sampler && typeof Sampler.framePixHash === "function") {
        snap.pixHash0 = Sampler.framePixHash(Acts.act, Acts.src0) >>> 0;
        snap.pixHash1 = Sampler.framePixHash(Acts.act, Acts.src1) >>> 0;
      }
      if (Sampler && typeof Sampler.frameCountsHash === "function") {
        snap.countsHash0 = Sampler.frameCountsHash(Acts.act, Acts.src0) >>> 0;
        snap.countsHash1 = Sampler.frameCountsHash(Acts.act, Acts.src1) >>> 0;
      }
      if (Sampler && typeof Sampler.frameImgWH === "function") {
        snap.img0 = Sampler.frameImgWH(Acts.act, Acts.src0);
        snap.img1 = Sampler.frameImgWH(Acts.act, Acts.src1);
      }
    } catch (_) {}
    try {
      if (typeof Particles !== "undefined" && Particles && Particles.desired && Particles.desired.length) {
        let h = 2166136261 >>> 0;
        // Hash a subset of cells for speed (still stable + comparable).
        const arr = Particles.desired;
        const step = max(1, floor(arr.length / 512));
        for (let i = 0; i < arr.length; i += step) {
          h ^= arr[i] & 0xffff;
          h = Math.imul(h, 16777619) >>> 0;
        }
        snap.desiredHash = h >>> 0;
      }
    } catch (_) {}
    console.log("[SNAP]", snap);
  }
}

function drawDebug(level, desiredSum, moved, mismatch) {
  const ready = Sampler.ready ? Sampler.ready(Acts.act) : 0;
  const q = Sampler.queueLen ? Sampler.queueLen() : 0;
  const typ = (Typography && Typography.debugString) ? Typography.debugString() : "";
  const dpr = (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1);
  const pd = (typeof pixelDensity === "function" ? pixelDensity() : 1) || 1;
  const z = (typeof pageZoom === "number" ? pageZoom : 1) || 1;
  const bufW = round(width * pd);
  const bufH = round(height * pd);
  const cssW = (typeof window !== "undefined" ? (window.innerWidth || 0) : 0);
  const cssH = (typeof window !== "undefined" ? (window.innerHeight || 0) : 0);
  let pix0 = 0, pix1 = 0, ch0 = 0, ch1 = 0;
  try {
    if (Sampler && typeof Sampler.framePixHash === "function") {
      pix0 = Sampler.framePixHash(Acts.act, Acts.src0) >>> 0;
      pix1 = Sampler.framePixHash(Acts.act, Acts.src1) >>> 0;
    }
    if (Sampler && typeof Sampler.frameCountsHash === "function") {
      ch0 = Sampler.frameCountsHash(Acts.act, Acts.src0) >>> 0;
      ch1 = Sampler.frameCountsHash(Acts.act, Acts.src1) >>> 0;
    }
  } catch (_) {}

  push();
  noStroke();
  fill(0, 160);
  rect(10, 10, 520, 122, 8);
  fill(255);
  textSize(12);
  textAlign(LEFT, TOP);
  text(`mode:${Acts.mode} act:${Acts.act} next:${Acts.next}  t:${t.toFixed(2)} el:${Acts.elapsed.toFixed(2)} dur:${Acts.dur.toFixed(2)} fps:${Acts.fps} dt:${_debugDtSec.toFixed(3)} rate:${_debugRate.toFixed(2)} vis:${sceneVis.toFixed(2)}`, 18, 16);
  const step = max(1, (FRAME_STEP && FRAME_STEP[Acts.act]) | 0);
  const framesPerCycle = Acts.cycle > 0 ? ceil(Acts.cycle / step) : 0;
  text(`cycle:${Acts.cycle} SRC:${SRC_COUNT[Acts.act] || 0} step:${step} fCycle:${framesPerCycle} ready:${ready}  src0/src1:${Acts.src0}/${Acts.src1} a:${Acts.alpha.toFixed(2)} cycles:${Acts.cycles}`, 18, 32);
  text(`grid:${COLS}x${ROWS} cell:${cellW.toFixed(2)} zoom:${z.toFixed(2)} dpr:${dpr.toFixed(2)} pd:${pd.toFixed(2)} css:${cssW}x${cssH} buf:${bufW}x${bufH} desiredSum:${desiredSum} moved:${moved} mismatch:${mismatch} meanDist:${debugMeanCellDist.toFixed(1)} catchUp:${catchUp.toFixed(2)} hot:${_hotCount} lane:${debugTransport}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} inFlight:${Sampler.inFlight || 0} q:${q} level:${level.toFixed(3)}  ${typ}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  if (pix0 || pix1 || ch0 || ch1) {
    text(`pixHash0/1:${pix0}/${pix1}  countsHash0/1:${ch0}/${ch1}`, 18, 96);
  }
  if (compareMode) {
    text(`COMPARE: act=${compareAct} src0=${compareSrc0} alpha=${compareAlpha.toFixed(2)}  (1-4 act, ←/→ src, ↑/↓ alpha, C toggle)`, 18, 80);
  }
  pop();
}
