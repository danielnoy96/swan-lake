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
}

function windowResized() {
  resizeCanvas(windowWidth, windowHeight);
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
  const elapsed = t - actStartT;

  if (elapsed > actDuration) {
    act++;
    if (act > 4) act = 1;
    actStartT = t;
  }
}

function drawAct(level) {
  // כותרת מערכה
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

  const scale = min(width / img.width, height / img.height);
  const dw = img.width * scale;
  const dh = img.height * scale;
  imageMode(CENTER);
  image(img, width / 2, height / 2, dw, dh);

  // ויזואליזציה של קול
}
