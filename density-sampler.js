// ---- Runtime sampler facade: precomputed SWD data by default, PNGs on explicit request ----
const LEGACY_SAMPLER_MODE = (() => {
  try { return new URLSearchParams(location.search || "").get("legacySampler") === "1"; }
  catch (_) { return false; }
})();

const DensitySampler = {
  cache: { 1: null, 2: null, 3: null, 4: null },
  states: { 1: "idle", 2: "idle", 3: "idle", 4: "idle" },
  errors: { 1: null, 2: null, 3: null, 4: null },
  promises: { 1: null, 2: null, 3: null, 4: null },
  inFlight: 0,
  hits: 0, misses: 0, hitsF: 0, missesF: 0,

  init() {},
  realloc() {},
  stepCompute() {},
  resetFrameStats() { this.hitsF = 0; this.missesF = 0; },
  queueLen() {
    let n = 0;
    for (const a of ACTS) if (this.states[a] === "loading") n++;
    return n;
  },
  computeQueueLen() { return 0; },
  _url(a) { return `assets/density/act${a}.swd?v=${encodeURIComponent(DENSITY_REVISION)}`; },
  _delay(ms) { return new Promise((resolve) => setTimeout(resolve, ms)); },

  async _fetchNetwork(a, cacheMode, allowRetries) {
    const url = this._url(a);
    const delays = allowRetries === false ? [0] : [0, 400, 1200];
    let lastError = null;
    for (let attempt = 0; attempt < delays.length; attempt++) {
      if (delays[attempt]) await this._delay(delays[attempt]);
      try {
        const response = await fetch(url, { cache: cacheMode || "default" });
        if (!response.ok) {
          const error = new Error(`HTTP ${response.status} while loading ${url}`);
          if (response.status >= 400 && response.status < 500) error.noRetry = true;
          throw error;
        }
        return await response.arrayBuffer();
      } catch (error) {
        lastError = error;
        if (error && error.noRetry) break;
      }
    }
    throw lastError || new Error(`Unable to load ${url}`);
  },

  _decode(a, buffer) {
    return DensityCodec.decodeAct(buffer, {
      act: a,
      frameCount: SRC_COUNT[a] || 0,
      cols: COLS,
      rows: ROWS,
      particleCount: N,
    });
  },

  loadAct(a) {
    a |= 0;
    if (!ACTS.includes(a) || !(SRC_COUNT[a] > 0)) return Promise.reject(new Error(`Unknown act ${a}`));
    if (this.states[a] === "ready") return Promise.resolve(this.cache[a]);
    if (this.promises[a]) return this.promises[a];
    this.states[a] = "loading";
    this.errors[a] = null;
    this.inFlight++;
    this.promises[a] = (async () => {
      let buffer = await this._fetchNetwork(a, "default", true);
      let decoded;
      try {
        decoded = this._decode(a, buffer);
      } catch (firstError) {
        try {
          buffer = await this._fetchNetwork(a, "reload", false);
          decoded = this._decode(a, buffer);
        } catch (secondError) {
          throw new Error(`${firstError.message}; cache-bypass retry failed: ${secondError.message}`);
        }
      }
      this.cache[a] = decoded;
      this.states[a] = "ready";
      return decoded;
    })().catch((error) => {
      this.states[a] = "error";
      this.errors[a] = error;
      const message = `[DENSITY ERROR]\nAct ${a}: ${error && error.message ? error.message : String(error)}\nURL: ${this._url(a)}\nUse ?legacySampler=1 for the PNG recovery path.`;
      window.__fatalError = message;
      throw error;
    }).finally(() => {
      this.inFlight = Math.max(0, this.inFlight - 1);
      this.promises[a] = null;
    });
    return this.promises[a];
  },

  actState(a) { return this.states[a | 0] || "error"; },
  progress() {
    let ready = 0, loading = 0, error = 0;
    for (const a of ACTS) {
      if (this.states[a] === "ready") ready++;
      else if (this.states[a] === "loading") loading++;
      else if (this.states[a] === "error") error++;
    }
    return { ready, loading, error, total: ACTS.length };
  },
  ensureActCache(a) { if (this.states[a] === "idle") this.loadAct(a).catch(() => {}); },
  ensure(a, idx, cycle, outOff) {
    a |= 0; idx |= 0;
    const off = idx * CELLS;
    if (outOff) outOff.off = off;
    const c = this.cache[a];
    if (this.states[a] === "ready" && c && c.cycle === cycle && idx >= 0 && idx < c.cycle) {
      this.hits++; this.hitsF++;
      return true;
    }
    this.misses++; this.missesF++;
    if (this.states[a] === "idle") this.loadAct(a).catch(() => {});
    return false;
  },
  counts(a) { return this.cache[a] ? this.cache[a].counts : null; },
  ready(a) { const c = this.cache[a]; return this.states[a] === "ready" && c ? c.cycle : 0; },
  readyOrFailed(a) { return this.ready(a); },
  isDone(a, idx) { const c = this.cache[a]; return !!(c && c.done && c.done[idx]); },
  framePixHash(a, idx) { const c = this.cache[a]; return c && c.pixHash ? (c.pixHash[idx] >>> 0) : 0; },
  frameCountsHash(a, idx) { const c = this.cache[a]; return c && c.countsHash ? (c.countsHash[idx] >>> 0) : 0; },
  frameImgWH(a, idx) {
    const c = this.cache[a];
    return c && c.imgW && c.imgH ? { w: c.imgW[idx] | 0, h: c.imgH[idx] | 0 } : { w: 0, h: 0 };
  },
  resetAllCaches() {
    for (const a of ACTS) {
      this.cache[a] = null;
      this.states[a] = "idle";
      this.errors[a] = null;
      this.promises[a] = null;
    }
  },
};

const Sampler = LEGACY_SAMPLER_MODE ? PngSampler : DensitySampler;
