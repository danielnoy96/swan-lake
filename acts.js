const Acts = {
  mode: "ACT", act: 1, next: 2, actStartT: 0, transitionStartT: 0, prepDone: false,
  elapsed: 0, dur: 0, fps: 24, cycle: 0, frameFloat: 0, idx: 0, src0: 0, src1: 0, alpha: 0, cycles: 0,
  compute() {
    const src = SRC_COUNT[this.act] || 0;
    this.cycle = src;
    this.fps = FPS_EFFECTIVE[this.act] || 24;
    this.elapsed = max(0, t - this.actStartT);
    const step = max(1, (FRAME_STEP && FRAME_STEP[this.act]) | 0);
    const framesPerCycle = this.cycle > 0 ? ceil(this.cycle / step) : 0;
    // Advance stepped-frame index at (fps/step): fewer density updates per second -> smoother / less CPU.
    this.frameFloat = step > 1 ? (this.elapsed * this.fps) / step : (this.elapsed * this.fps);
    this.dur = this.cycle > 0 && this.fps > 0 ? this.cycle / this.fps : 0;
    if (this.cycle <= 0) { this.idx = this.src0 = this.src1 = 0; this.alpha = 0; this.cycles = 0; return; }
    const base = floor(this.frameFloat);
    this.idx = framesPerCycle > 0 ? (base % framesPerCycle) : 0;
    this.src0 = (this.idx * step) % this.cycle;
    this.src1 = (this.src0 + step) % this.cycle;
    this.alpha = constrain(this.frameFloat - base, 0, 1);
    this.cycles = framesPerCycle > 0 ? floor(this.frameFloat / framesPerCycle) : 0;
  },
  startTransition() {
    this.mode = "TRANSITION";
    this.transitionStartT = t;
    this.next = nextActWithFrames(this.act);
    this.prepDone = false;
    // Give the flock an immediate "take off" so Stage A fills the whole space.
    if (typeof Particles !== "undefined" && Particles && Particles.onTransitionStart) Particles.onTransitionStart(this.act, this.next);
  },
  endTransition() {
    this.act = this.next;
    this.actStartT = t;
    this.mode = "ACT";
    this.compute();
    if (this.cycle > 0) {
      Sampler.resetFrameStats();
      if (Sampler.ensure(this.act, 0, this.cycle, _off0)) {
        const counts = Sampler.counts(this.act);
        Particles.setDesired(counts, _off0.off, _off0.off, 0);
        Particles.setWalkFromCounts(counts, _off0.off, _off0.off);
        Particles.softRetargetToDesired();
      }
    }
  },
  update() {
    this.compute();
    if (this.mode === "ACT") {
      if (this.cycle <= 0 || this.dur <= 0) return this.startTransition();
      if (this.cycles >= 1 || this.elapsed >= this.dur) return this.startTransition();
    } else {
      if (t - this.transitionStartT >= TRANSITION_DURATION) return this.endTransition();
    }
  },
};

const _off0 = { off: 0 }, _off1 = { off: 0 };

