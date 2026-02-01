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

let actFrames = { 1: [], 2: [], 3: [], 4: [] };
let actFrameCounts = { 1: 192, 2: 172, 3: 0, 4: 0 };

// STEP 8: state machine + fluid material
let mode = "ACT"; // "ACT" | "TRANSITION"
let nextAct = 2;
let transitionStartT = 0;
let transitionDuration = 1.5;

// Offscreen buffers (2D, mobile-friendly)
let srcG, outG, frameAG, frameBG;
let dispG, maskG;
let dispW = 120;
let dispH = 200;

// Coarse tiling for fast displacement (no per-pixel full-res work)
let tileCols = 40;
let tileRows = 70;

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

  initBuffers();
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
  initBuffers();
}

function initBuffers() {
  srcG = createGraphics(width, height);
  outG = createGraphics(width, height);
  frameAG = createGraphics(width, height);
  frameBG = createGraphics(width, height);

  dispG = createGraphics(dispW, dispH);
  maskG = createGraphics(dispW, dispH);

  srcG.pixelDensity(1);
  outG.pixelDensity(1);
  frameAG.pixelDensity(1);
  frameBG.pixelDensity(1);
  dispG.pixelDensity(1);
  maskG.pixelDensity(1);

  srcG.imageMode(CORNER);
  outG.imageMode(CORNER);
  frameAG.imageMode(CORNER);
  frameBG.imageMode(CORNER);

  const aspect = height / max(1, width);
  tileCols = 40;
  tileRows = constrain(floor(tileCols * aspect * 1.75), 50, 90);
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

  updateAct();
  drawAct(level);
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

function updateAct() {
  if (mode === "ACT") {
    const elapsed = t - actStartT;

    if (elapsed >= actDuration) {
      mode = "TRANSITION";
      transitionStartT = t;
      nextAct = act + 1;
      if (nextAct > 4) nextAct = 1;

      for (let k = 0; k < 4; k++) {
        const frames = actFrames[nextAct] || [];
        if (frames.length > 0) break;
        nextAct++;
        if (nextAct > 4) nextAct = 1;
      }
    }

    return;
  }

  if (mode === "TRANSITION") {
    const elapsed = t - transitionStartT;

    if (elapsed > transitionDuration) {
      act = nextAct;
      actStartT = t;
      mode = "ACT";
    }
  }
}

function smoothstep(a, b, x) {
  const t01 = constrain((x - a) / max(0.000001, b - a), 0, 1);
  return t01 * t01 * (3 - 2 * t01);
}

function drawFrameFit(g, img) {
  g.clear();
  g.background(0);

  const scale = min(width / img.width, height / img.height);
  const dw = img.width * scale;
  const dh = img.height * scale;

  g.imageMode(CENTER);
  g.image(img, width / 2, height / 2, dw, dh);
  g.imageMode(CORNER);
}

function updateFields() {
  // Low-res displacement field (RG) + localized transition mask (R)
  dispG.loadPixels();
  maskG.loadPixels();

  const df = 2.2;
  const mf = 2.0;
  const dz = t * 0.9;
  const mz = t * 0.7;

  for (let y = 0; y < dispH; y++) {
    const v = y / max(1, dispH - 1);
    for (let x = 0; x < dispW; x++) {
      const u = x / max(1, dispW - 1);

      const n1 = noise(u * df, v * df, dz);
      const n2 = noise(u * df + 10.0, v * df + 10.0, dz + 5.0);

      const dx = u - 0.5;
      const dy = v - 0.5;
      const r = sqrt(dx * dx + dy * dy);
      const radial = smoothstep(0.75, 0.1, r); // 1 at center, 0 at edges
      const nm = noise(u * mf + 20.0, v * mf + 20.0, mz);
      const mask = constrain(radial * (0.75 + 0.25 * nm), 0, 1);

      const di = 4 * (y * dispW + x);
      dispG.pixels[di + 0] = floor(n1 * 255);
      dispG.pixels[di + 1] = floor(n2 * 255);
      dispG.pixels[di + 2] = 0;
      dispG.pixels[di + 3] = 255;

      maskG.pixels[di + 0] = floor(mask * 255);
      maskG.pixels[di + 1] = 0;
      maskG.pixels[di + 2] = 0;
      maskG.pixels[di + 3] = 255;
    }
  }

  dispG.updatePixels();
  maskG.updatePixels();
}

function blendFramesToSrc(transP) {
  srcG.clear();
  srcG.background(0);

  const tileW = width / tileCols;
  const tileH = height / tileRows;

  maskG.loadPixels();

  for (let row = 0; row < tileRows; row++) {
    const y = row * tileH;
    const cy = y + tileH * 0.5;
    const my = floor((cy / max(1, height)) * (dispH - 1));

    for (let col = 0; col < tileCols; col++) {
      const x = col * tileW;
      const cx = x + tileW * 0.5;
      const mx = floor((cx / max(1, width)) * (dispW - 1));

      const mi = 4 * (my * dispW + mx);
      const maskVal = (maskG.pixels[mi] || 0) / 255;

      // Center transitions first, edges last (localized reveal, but completes by the end)
      const threshold = 1.0 - maskVal;
      const b = smoothstep(threshold - 0.12, threshold + 0.12, transP);

      srcG.image(frameAG, x, y, tileW, tileH, x, y, tileW, tileH);

      if (b > 0.001) {
        srcG.tint(255, 255 * b);
        srcG.image(frameBG, x, y, tileW, tileH, x, y, tileW, tileH);
        srcG.noTint();
      }
    }
  }
}

function applyMaterial(strengthPx) {
  outG.clear();
  outG.background(0);

  const tileW = width / tileCols;
  const tileH = height / tileRows;

  dispG.loadPixels();

  for (let row = 0; row < tileRows; row++) {
    const y = row * tileH;
    const cy = y + tileH * 0.5;
    const vy = cy / max(1, height);
    const sy = (cy - height * 0.5) / max(1, height);

    for (let col = 0; col < tileCols; col++) {
      const x = col * tileW;
      const cx = x + tileW * 0.5;
      const ux = cx / max(1, width);
      const sx = (cx - width * 0.5) / max(1, width);

      const px = floor(ux * (dispW - 1));
      const py = floor(vy * (dispH - 1));
      const di = 4 * (py * dispW + px);
      const r = ((dispG.pixels[di + 0] || 128) / 255 - 0.5) * 2.0;
      const g = ((dispG.pixels[di + 1] || 128) / 255 - 0.5) * 2.0;

      // Swirl-like flow around center (stronger near subject area)
      const centerBoost = constrain(1.0 - (sx * sx + sy * sy) * 3.0, 0, 1);
      const swirl = strengthPx * 0.9 * centerBoost;

      const dx = r * strengthPx + (-sy) * swirl;
      const dy = g * strengthPx + (sx) * swirl;

      outG.image(srcG, x + dx, y + dy, tileW, tileH, x, y, tileW, tileH);
    }
  }
}

function drawAct(level) {
  // כותרת מערכה
  // Material strength: always subtle, increases with sound
  const matStrengthPx = 1.2 + map(level, 0, 0.2, 0, 14, true);

  updateFields();

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

    const elapsed = max(0, t - actStartT);
    const progress = actDuration > 0 ? elapsed / actDuration : 0;
    const idx = constrain(floor(progress * n), 0, n - 1);
    const img = frames[idx];

    if (!img) return;

    drawFrameFit(srcG, img);
    applyMaterial(matStrengthPx);
    image(outG, 0, 0, width, height);
    return;
  }

  if (mode === "TRANSITION") {
    const fromFrames = actFrames[act] || [];
    const toFrames = actFrames[nextAct] || [];

    if (fromFrames.length === 0 || toFrames.length === 0) {
      actStartT = t;
      mode = "ACT";
      return;
    }

    const fromImg = fromFrames[fromFrames.length - 1] || fromFrames[0];
    const transP = constrain((t - transitionStartT) / max(0.0001, transitionDuration), 0, 1);
    const toIdx = constrain(floor(transP * toFrames.length), 0, toFrames.length - 1);
    const toImg = toFrames[toIdx] || toFrames[0];
    if (!fromImg || !toImg) return;

    drawFrameFit(frameAG, fromImg);
    drawFrameFit(frameBG, toImg);

    blendFramesToSrc(transP);
    applyMaterial(matStrengthPx);
    image(outG, 0, 0, width, height);
  }

  // ויזואליזציה של קול
}
