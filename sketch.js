function preload() {
  Typography.preload();
}

function setup() {
  createCanvas(windowWidth, windowHeight);

  mic = new p5.AudioIn();
  amp = new p5.Amplitude();

  Style.init();
  resizeGrid();
  Sampler.init();
  Particles.resizeBins();
  Particles.init();

  const styleBtn = (btn) => {
    btn.style("position", "fixed");
    btn.style("z-index", "10");
    btn.style("padding", "8px 10px");
    btn.style("font-size", "14px");
    btn.style("border", "1px solid rgba(255,255,255,0.25)");
    btn.style("background", "rgba(0,0,0,0.55)");
    btn.style("color", "#fff");
    btn.style("border-radius", "8px");
  };

  const btnGrid = createButton("Grid: OFF");
  btnGrid.position(12, 12);
  styleBtn(btnGrid);
  btnGrid.mousePressed(() => {
    showGrid = !showGrid;
    btnGrid.html(showGrid ? "Grid: ON" : "Grid: OFF");
  });

  const btnAuto = createButton("Auto: OFF");
  btnAuto.position(110, 12);
  styleBtn(btnAuto);
  btnAuto.mousePressed(() => {
    autoRun = !autoRun;
    if (autoRun) audioStarted = true;
    btnAuto.html(autoRun ? "Auto: ON" : "Auto: OFF");
  });

  Typography.init();
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
  Style.init();
  resizeGrid();
  Particles.resizeBins();
  Particles.init();
  Sampler.resetAllCaches();
  Typography.resize();
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
  if (autoRun) t += AUTO_SPEED;
  else if (micRunning && level > MIC_THRESHOLD) t += map(level, MIC_THRESHOLD, 0.2, 0.01, 0.05, true);

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

          for (let i = 0; i < CELLS; i++) desiredSum += Particles.desired[i];

          let maxD = 1;
          for (let i = 0; i < CELLS; i++) { const d = Particles.desired[i] | 0; if (d > maxD) maxD = d; }
          Render.maxD = maxD;

          if (tAdvanced) {
            mismatch = estimateMismatch();
            debugMeanCellDist = estimateMeanCellDist();

            const posCU = constrain((debugMeanCellDist - cellW * 0.85) / (cellW * 3.2), 0, 1);
            const targetCU = constrain((mismatch - N * 0.10) / (N * 0.45), 0, 1);
            catchUp = max(catchUp * 0.94, targetCU, posCU * 0.35);

            const passes = (mismatch > 900 ? 10 : mismatch > 600 ? 8 : mismatch > 350 ? 6 : mismatch > 180 ? 4 : 3) + floor(catchUp * 2);
            computeDeficitHotspots();
            moved += Particles.rebalancePasses(passes);
            if (catchUp > 0.32) catchUpNudge(catchUp);

            Particles.updateInsideCells(level, 0.18 + catchUp * 0.26, 0.15, catchUp * 0.7);
            Particles.separate(0.018, 3.2);
            catchUp *= 0.97;
          }
        } else {
          if (tAdvanced) Particles.dance(level, 0.20);
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
              const ok0 = Sampler.ensure(a, 0, cycle, _off0);
              const ok1 = Sampler.ensure(a, (cycle > 1 ? 1 : 0), cycle, _off1);
              if (ok0 && ok1) Acts.prepDone = true;
            }
          }
          Particles.dance(level, s);
        } else {
          const a = Acts.next;
          const cycle = SRC_COUNT[a] || 0;
          if (cycle > 0) {
            const ok0 = Sampler.ensure(a, 0, cycle, _off0);
            const ok1 = Sampler.ensure(a, 1 % cycle, cycle, _off1);
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
  if (audioStarted) return;
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
}

function drawDebug(level, desiredSum, moved, mismatch) {
  const ready = Sampler.ready ? Sampler.ready(Acts.act) : 0;
  const q = Sampler.queueLen ? Sampler.queueLen() : 0;
  const typ = (Typography && Typography.debugString) ? Typography.debugString() : "";

  push();
  noStroke();
  fill(0, 160);
  rect(10, 10, 520, 106, 8);
  fill(255);
  textSize(12);
  textAlign(LEFT, TOP);
  text(`mode:${Acts.mode} act:${Acts.act} next:${Acts.next}  t:${t.toFixed(2)} el:${Acts.elapsed.toFixed(2)} dur:${Acts.dur.toFixed(2)} fps:${Acts.fps} vis:${sceneVis.toFixed(2)}`, 18, 16);
  text(`cycle:${Acts.cycle} SRC:${SRC_COUNT[Acts.act] || 0} ready:${ready}  src0/src1:${Acts.src0}/${Acts.src1} a:${Acts.alpha.toFixed(2)} cycles:${Acts.cycles}`, 18, 32);
  text(`grid:${COLS}x${ROWS} desiredSum:${desiredSum} moved:${moved} mismatch:${mismatch} meanDist:${debugMeanCellDist.toFixed(1)} catchUp:${catchUp.toFixed(2)} hot:${_hotCount} lane:${debugTransport}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} inFlight:${Sampler.inFlight || 0} q:${q} level:${level.toFixed(3)}  ${typ}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  pop();
}

