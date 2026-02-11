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

function setup() {
  const cnv = createCanvas(windowWidth, windowHeight);
  applyPixelDensity();
  try { cnv?.style?.("display", "block"); } catch (_) {}
  try {
    mic = new p5.AudioIn();
    amp = new p5.Amplitude();

    Style.init();
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
  try {
    Style.init();
    resizeGrid();
    Particles.resizeBins();
    Particles.init();
    Sampler.resetAllCaches();
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

  if (autoRun) {
    t += AUTO_SPEED * dtSec;
    _debugRate = AUTO_SPEED;
  } else if (micRunning && level > MIC_THRESHOLD) {
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
    Acts.update();
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
  text("Tap to activate microphone\nSound drives the animation", width / 2, height / 2);
}

function mousePressed() {
  // Allow starting the mic even if we previously ran in Auto mode.
  if (autoRun) return;
  if (micRunning) return;
  userStartAudio();
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
}

function drawDebug(level, desiredSum, moved, mismatch) {
  const ready = Sampler.ready ? Sampler.ready(Acts.act) : 0;
  const q = Sampler.queueLen ? Sampler.queueLen() : 0;
  const typ = (Typography && Typography.debugString) ? Typography.debugString() : "";
  const dpr = (typeof window !== "undefined" ? (window.devicePixelRatio || 1) : 1);
  const pd = (typeof pixelDensity === "function" ? pixelDensity() : 1) || 1;

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
  text(`grid:${COLS}x${ROWS} cell:${cellW.toFixed(2)} dpr:${dpr.toFixed(2)} pd:${pd.toFixed(2)} desiredSum:${desiredSum} moved:${moved} mismatch:${mismatch} meanDist:${debugMeanCellDist.toFixed(1)} catchUp:${catchUp.toFixed(2)} hot:${_hotCount} lane:${debugTransport}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} inFlight:${Sampler.inFlight || 0} q:${q} level:${level.toFixed(3)}  ${typ}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  pop();
}
