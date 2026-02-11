function preload() {
  try {
    Typography.preload();
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
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
  const cnv = createCanvas(windowWidth, windowHeight);
  applyPixelDensity();
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
    Sampler.init();
    Particles.resizeBins();
    Particles.init();
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
  resizeCanvas(windowWidth, windowHeight);
  applyPixelDensity();
  applySimSeed();
  try {
    Style.init();
    updatePageZoom();
    resizeGrid();
    Particles.resizeBins();
    Particles.init();
    Typography.resize();
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
  }
}

function draw() {
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

  if (!audioStarted && !autoRun) {
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
        if (tAdvanced) {
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
  // Allow starting the mic even if we previously ran in Auto mode.
  if (autoRun) return;
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
    Acts.actStartT = t;
  });
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
      w: width, h: height,
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

  push();
  noStroke();
  fill(0, 160);
  rect(10, 10, 520, 106, 8);
  fill(255);
  textSize(12);
  textAlign(LEFT, TOP);
  text(`mode:${Acts.mode} act:${Acts.act} next:${Acts.next}  t:${t.toFixed(2)} el:${Acts.elapsed.toFixed(2)} dur:${Acts.dur.toFixed(2)} fps:${Acts.fps} dt:${_debugDtSec.toFixed(3)} rate:${_debugRate.toFixed(2)} vis:${sceneVis.toFixed(2)}`, 18, 16);
  const step = max(1, (FRAME_STEP && FRAME_STEP[Acts.act]) | 0);
  const framesPerCycle = Acts.cycle > 0 ? ceil(Acts.cycle / step) : 0;
  text(`cycle:${Acts.cycle} SRC:${SRC_COUNT[Acts.act] || 0} step:${step} fCycle:${framesPerCycle} ready:${ready}  src0/src1:${Acts.src0}/${Acts.src1} a:${Acts.alpha.toFixed(2)} cycles:${Acts.cycles}`, 18, 32);
  text(`grid:${COLS}x${ROWS} cell:${cellW.toFixed(2)} zoom:${z.toFixed(2)} dpr:${dpr.toFixed(2)} pd:${pd.toFixed(2)} buf:${bufW}x${bufH} desiredSum:${desiredSum} moved:${moved} mismatch:${mismatch} meanDist:${debugMeanCellDist.toFixed(1)} catchUp:${catchUp.toFixed(2)} hot:${_hotCount} lane:${debugTransport}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} inFlight:${Sampler.inFlight || 0} q:${q} level:${level.toFixed(3)}  ${typ}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  if (compareMode) {
    text(`COMPARE: act=${compareAct} src0=${compareSrc0} alpha=${compareAlpha.toFixed(2)}  (1-4 act, ←/→ src, ↑/↓ alpha, C toggle)`, 18, 80);
  }
  pop();
}
