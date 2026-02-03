// Color styling (per act)
const ACT_COLORS = {
  1: { bg: "#5F6E7C", hi: "#EFE9DD", lo: "#B8A76A" },
  2: { bg: "#6B665E", hi: "#F5F6F4", lo: "#C7CCD2" },
  3: { bg: "#121317", hi: "#4A1F2D", lo: "#7A2734" },
  4: { bg: "#3E4A57", hi: "#2A2D33", lo: "#BFC4CB" },
};
const ACT_FOCUS = {
  1: { x: 0.52, y: 0.55 },
  2: { x: 0.50, y: 0.54 },
  3: { x: 0.52, y: 0.52 },
  4: { x: 0.52, y: 0.54 },
};
const ACT_STAIN = {
  // Relative stain size (minor color) per act; act2 has more secondary color.
  1: { size: 0.40, ox: -0.10, oy: 0.04, a: 0.55 },
  2: { size: 0.72, ox: 0.12, oy: 0.06, a: 1.00 },
  3: { size: 0.62, ox: 0.08, oy: -0.06, a: 1.00 },
  4: { size: 0.26, ox: 0.10, oy: -0.04, a: 0.32 },
};

const Render = {
  bg: null,
  hi: null,
  lo: null,
  fx: 0,
  fy: 0,
  outerR2: 0,
  stainX: 0,
  stainY: 0,
  stainR: 0,
  tintA: 0,
  tintB: 0,
  maxD: 1,
  colorAmt: 0,
};

const Style = {
  scheme: { 1: null, 2: null, 3: null, 4: null },
  white: null,
  init() {
    this.white = color(255, 255, 255);
    for (let a of ACTS) {
      const base = ACT_COLORS[a] || ACT_COLORS[1];
      this.scheme[a] = { bg: color(base.bg), hi: color(base.hi), lo: color(base.lo) };
    }
  },
  schemeFor(a) { return this.scheme[a] || this.scheme[1]; },
  focusFor(a) { return ACT_FOCUS[a] || ACT_FOCUS[1]; },
  stainFor(a) { return ACT_STAIN[a] || ACT_STAIN[1]; },
  update(acts, level) {
    const s0 = this.schemeFor(acts.act);
    const s1 = this.schemeFor(acts.next);
    const p = acts.mode === "TRANSITION" ? constrain((t - acts.transitionStartT) / max(1e-6, TRANSITION_DURATION), 0, 1) : 0;

    const bgTarget = acts.mode === "TRANSITION" ? lerpColor(s0.bg, s1.bg, p) : s0.bg;

    const targetAmt = pow(constrain(map(level, 0.004, 0.045, 0, 1, true), 0, 1), 0.85);
    const k = targetAmt > Render.colorAmt ? 0.06 : 0.02;
    Render.colorAmt += (targetAmt - Render.colorAmt) * k;
    const colorAmt = Render.colorAmt;

    Render.bg = lerpColor(color(0), bgTarget, colorAmt);

    if (acts.mode === "TRANSITION") {
      const stageA = p < DANCE_PORTION;
      const s = stageA ? p / max(1e-6, DANCE_PORTION) : (p - DANCE_PORTION) / max(1e-6, 1 - DANCE_PORTION);
      const hiBlend = lerpColor(s0.hi, s1.hi, s);
      const loBlend = lerpColor(s0.lo, s1.lo, s);
      Render.hi = lerpColor(this.white, hiBlend, colorAmt);
      Render.lo = lerpColor(this.white, loBlend, colorAmt);
    } else {
      Render.hi = lerpColor(this.white, s0.hi, colorAmt);
      Render.lo = lerpColor(this.white, s0.lo, colorAmt);
    }

    // Grouped tint region parameters (blended by act during transitions)
    const f0 = this.focusFor(acts.act);
    const f1 = this.focusFor(acts.next);
    Render.fx = lerp(width * f0.x, width * f1.x, p);
    Render.fy = lerp(height * f0.y, height * f1.y, p);

    const minD = min(width, height);
    const coreR = minD * (0.08 + 0.16 * colorAmt);
    const outerR = coreR * 1.65;
    Render.outerR2 = outerR * outerR;

    // Single stain params (blended by act during transitions)
    const st0 = this.stainFor(acts.act);
    const st1 = this.stainFor(acts.next);
    const stA0 = st0.a == null ? 1 : st0.a;
    const stA1 = st1.a == null ? 1 : st1.a;
    const stainAlphaScale = lerp(stA0, stA1, p);
    Render.tintA = 110 * colorAmt;
    Render.tintB = 185 * colorAmt * stainAlphaScale;

    const stSize = lerp(st0.size, st1.size, p);
    const stOx = lerp(st0.ox, st1.ox, p);
    const stOy = lerp(st0.oy, st1.oy, p);
    const drift = outerR * 0.08 * colorAmt;
    const ndx = (noise(100 + acts.act * 7.1, t * 0.06) - 0.5) * 2;
    const ndy = (noise(200 + acts.act * 9.3, t * 0.06) - 0.5) * 2;
    Render.stainX = Render.fx + outerR * stOx + ndx * drift;
    Render.stainY = Render.fy + outerR * stOy + ndy * drift;
    Render.stainR = outerR * stSize;
  },
};

