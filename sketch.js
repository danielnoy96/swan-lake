// ---- p5 lifecycle ----
// Typography overlay (top layer)
let typographyImg = null;
const typographyRect = { x: 0, y: 0, w: 0, h: 0 };

function preload() {
  // Only typography is preloaded (tiny). Act frames are loaded on-demand in Sampler.ensure().
  typographyImg = loadImage("assets/typography/swan_lake.png");
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

  const btnGrid = createButton("Grid: OFF");
  btnGrid.position(12, 12);
  btnGrid.style("position", "fixed");
  btnGrid.style("z-index", "10");
  btnGrid.style("padding", "8px 10px");
  btnGrid.style("font-size", "14px");
  btnGrid.style("border", "1px solid rgba(255,255,255,0.25)");
  btnGrid.style("background", "rgba(0,0,0,0.55)");
  btnGrid.style("color", "#fff");
  btnGrid.style("border-radius", "8px");
  btnGrid.mousePressed(() => {
    showGrid = !showGrid;
    btnGrid.html(showGrid ? "Grid: ON" : "Grid: OFF");
  });

  const btnAuto = createButton("Auto: OFF");
  btnAuto.position(110, 12);
  btnAuto.style("position", "fixed");
  btnAuto.style("z-index", "10");
  btnAuto.style("padding", "8px 10px");
  btnAuto.style("font-size", "14px");
  btnAuto.style("border", "1px solid rgba(255,255,255,0.25)");
  btnAuto.style("background", "rgba(0,0,0,0.55)");
  btnAuto.style("color", "#fff");
  btnAuto.style("border-radius", "8px");
  btnAuto.mousePressed(() => {
    autoRun = !autoRun;
    if (autoRun) audioStarted = true;
    btnAuto.html(autoRun ? "Auto: ON" : "Auto: OFF");
  });

  updateTypographyRect();
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
  Style.init();
  resizeGrid();
  Particles.resizeBins();
  Particles.init();
  Sampler.resetAllCaches();
  updateTypographyRect();
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
    drawTypographyOverlay();
    return;
  }

  const level = autoRun ? 0.09 : (micRunning ? amp.getLevel() : 0);
  if (autoRun) t += AUTO_SPEED;
  else if (micRunning && level > MIC_THRESHOLD) t += map(level, MIC_THRESHOLD, 0.2, 0.01, 0.05, true);
  tDelta = t - _prevT;
  tAdvanced = tDelta !== 0;
  _prevT = t;

  try {
    Sampler.resetFrameStats();
    Acts.update();

  Style.update(Acts, level);
  background(Render.bg);

  let desiredSum = 0, moved = 0;
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

        // Keep max density for subtle stain density influence.
        let maxD = 1;
        for (let i = 0; i < CELLS; i++) { const d = Particles.desired[i] | 0; if (d > maxD) maxD = d; }
        Render.maxD = maxD;

        if (tAdvanced) {
          mismatch = estimateMismatch();
          debugMeanCellDist = estimateMeanCellDist();
          // Position-based catch-up helps prevent "stuck in space", but keep it mild.
          const posCU = constrain((debugMeanCellDist - cellW * 0.85) / (cellW * 3.2), 0, 1);
          const targetCU = constrain((mismatch - N * 0.10) / (N * 0.45), 0, 1);
          catchUp = max(catchUp * 0.94, targetCU, posCU * 0.35);
          const passes = (mismatch > 900 ? 10 : mismatch > 600 ? 8 : mismatch > 350 ? 6 : mismatch > 180 ? 4 : 3) + floor(catchUp * 2);
          // IMPORTANT: compute hotspots BEFORE rebalancing so surplus "islands" can find long-range deficits
          // in the same update tick (otherwise particles may look "stuck" until the next t advance).
          computeDeficitHotspots();
          moved += Particles.rebalancePasses(passes);
          if (catchUp > 0.32) catchUpNudge(catchUp);
          // Catch-up increases snap/attraction but stays smooth (no teleport)
          Particles.updateInsideCells(level, 0.18 + catchUp * 0.26, 0.15, catchUp * 0.7);
          Particles.separate(0.018, 3.2);
          catchUp *= 0.97;
        }
      }
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
        // Prep next-act cache early so Stage B doesn't "stall" on first loadPixels() call.
        if (!Acts.prepDone) {
          const a = Acts.next;
          const cycle = SRC_COUNT[a] || 0;
          if (cycle > 0) {
            const idx0 = 0;
            const idx1 = cycle > 1 ? 1 : 0;
            const ok0 = Sampler.ensure(a, idx0, cycle, _off0);
            const ok1 = Sampler.ensure(a, idx1, cycle, _off1);
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
            // Coalesce with the SAME pace as the flight: fewer reassign passes, softer catch-up, no hard snap.
            moved += Particles.rebalancePasses(4);
            if (p > (DANCE_PORTION + 0.10)) catchUpNudge(0.12);
            const snap01 = map(s, 0.70, 1.0, 0, 0.25, true);
            Particles.updateInsideCells(level, snap01, 0.0, 0.18);
            Particles.separate(0.020, 3.2);
          } else {
            // If frames aren't ready yet, keep moving instead of freezing.
            Particles.dance(level, 1.0);
          }
        }
      }
    }
    Particles.draw();
  }

    if (debugOn) drawDebug(level, desiredSum, moved, mismatch);
    if (showGrid) drawGridOverlay();
    drawTypographyOverlay();
  } catch (e) {
    window.__fatalError = e?.stack || String(e);
    background(0);
    fill(255);
    noStroke();
    textAlign(LEFT, TOP);
    textSize(12);
    text(String(window.__fatalError), 12, 12, width - 24, height - 24);
    noLoop();
  }
}

function drawStartScreen() {
  fill(255); noStroke(); textAlign(CENTER, CENTER); textSize(18);
  text("Tap to activate microphone\nSound drives the animation", width / 2, height / 2);
}

function updateTypographyRect() {
  if (!typographyImg || !typographyImg.width) return;
  fitRect(typographyImg.width, typographyImg.height, width, height, typographyRect);
}

function drawTypographyOverlay() {
  if (!typographyImg || !typographyImg.width) return;
  push();
  // SCREEN makes the black areas invisible while keeping the white typography.
  blendMode(SCREEN);
  imageMode(CORNER);
  tint(255, 215);
  image(typographyImg, typographyRect.x, typographyRect.y, typographyRect.w, typographyRect.h);
  blendMode(BLEND);
  pop();
}

function mousePressed() {
  if (audioStarted) return;
  userStartAudio();
  mic.start(() => { amp.setInput(mic); audioStarted = true; micRunning = true; Acts.actStartT = t; });
}

function keyPressed() { if (key === "d" || key === "D") debugOn = !debugOn; }

function drawDebug(level, desiredSum, moved, mismatch) {
  const ready = Sampler.ready ? Sampler.ready(Acts.act) : 0;
  const q = Sampler.queueLen ? Sampler.queueLen() : 0;
  push(); noStroke(); fill(0, 160); rect(10, 10, 470, 106, 8); fill(255); textSize(12); textAlign(LEFT, TOP);
  text(`mode:${Acts.mode} act:${Acts.act} next:${Acts.next}  t:${t.toFixed(2)} el:${Acts.elapsed.toFixed(2)} dur:${Acts.dur.toFixed(2)} fps:${Acts.fps}`, 18, 16);
  text(`cycle:${Acts.cycle} SRC:${SRC_COUNT[Acts.act] || 0} ready:${ready}  src0/src1:${Acts.src0}/${Acts.src1} a:${Acts.alpha.toFixed(2)} cycles:${Acts.cycles}`, 18, 32);
  text(`grid:${COLS}x${ROWS} desiredSum:${desiredSum} moved:${moved} mismatch:${mismatch} meanDist:${debugMeanCellDist.toFixed(1)} catchUp:${catchUp.toFixed(2)} hot:${_hotCount} lane:${debugTransport}  cache(h/m):${Sampler.hitsF}/${Sampler.missesF} total:${Sampler.hits}/${Sampler.misses} inFlight:${Sampler.inFlight || 0} q:${q} level:${level.toFixed(3)}`, 18, 48);
  text(`imagesDrawnToCanvas:${imagesDrawnToCanvas}`, 18, 64);
  pop();
}
