const DensityGenerator = {
  running: false,
  generated: {},
  urls: [],
  totalFrames: ACTS.reduce((sum, a) => sum + (SRC_COUNT[a] || 0), 0),
  processedBeforeAct: 0,

  get logEl() { return document.getElementById("log"); },
  get progressEl() { return document.getElementById("progress"); },
  log(message) {
    this.logEl.textContent += `${message}\n`;
    this.logEl.scrollTop = this.logEl.scrollHeight;
  },
  setBusy(busy) {
    this.running = busy;
    document.getElementById("generate").disabled = busy;
    document.getElementById("verify").disabled = busy;
  },
  nextFrame() { return new Promise((resolve) => requestAnimationFrame(resolve)); },

  async computeAct(a) {
    const cycle = SRC_COUNT[a] || 0;
    PngSampler.ensureActCache(a, cycle);
    const dummy = { off: 0 };
    for (let idx = 0; idx < cycle; idx++) PngSampler.ensure(a, idx, cycle, dummy);
    while (PngSampler.readyOrFailed(a) < cycle || PngSampler.queueLen() || PngSampler.computeQueueLen() || PngSampler.inFlight) {
      PngSampler.stepCompute(8);
      const ready = PngSampler.readyOrFailed(a);
      this.progressEl.value = this.processedBeforeAct + ready;
      if ((ready & 15) === 0) this.logEl.dataset.progress = `${a}:${ready}`;
      await this.nextFrame();
    }
    const cache = PngSampler.cache[a];
    let failed = 0;
    for (let i = 0; i < cycle; i++) failed += cache.failed[i] ? 1 : 0;
    if (failed) throw new Error(`Act ${a} has ${failed} missing or failed source frames`);
    this.processedBeforeAct += cycle;
    this.progressEl.value = this.processedBeforeAct;
    return cache;
  },

  encodeAndCheck(a, cache) {
    const encoded = DensityCodec.encodeAct({
      act: a,
      frameCount: cache.cycle,
      cols: COLS,
      rows: ROWS,
      particleCount: N,
      counts: cache.counts,
      pixHash: cache.pixHash,
      countsHash: cache.countsHash,
      imgW: cache.imgW,
      imgH: cache.imgH,
    });
    const decoded = DensityCodec.decodeAct(encoded, {
      act: a, frameCount: cache.cycle, cols: COLS, rows: ROWS, particleCount: N,
    });
    if (decoded.counts.length !== cache.counts.length) throw new Error(`Act ${a} round-trip length mismatch`);
    for (let i = 0; i < cache.counts.length; i++) {
      if (decoded.counts[i] !== cache.counts[i]) throw new Error(`Act ${a} round-trip mismatch at count ${i}`);
    }
    for (let i = 0; i < cache.cycle; i++) {
      if (decoded.pixHash[i] !== cache.pixHash[i] || decoded.countsHash[i] !== cache.countsHash[i]) {
        throw new Error(`Act ${a} metadata hash mismatch at frame ${i}`);
      }
    }
    return encoded;
  },

  addDownload(a, bytes) {
    let binary = "";
    const chunk = 0x8000;
    for (let i = 0; i < bytes.length; i += chunk) binary += String.fromCharCode(...bytes.subarray(i, i + chunk));
    const url = `data:application/octet-stream;base64,${btoa(binary)}`;
    const link = document.createElement("a");
    link.className = "download";
    link.href = url;
    link.download = `act${a}.swd`;
    link.textContent = `Download act${a}.swd (${(bytes.length / 1048576).toFixed(2)} MiB)`;
    document.getElementById("downloads").appendChild(link);
  },

  clearDownloads() {
    this.urls.length = 0;
    this.generated = {};
    document.getElementById("downloads").replaceChildren();
  },

  async prepareSources() {
    this.processedBeforeAct = 0;
    this.progressEl.max = this.totalFrames;
    this.progressEl.value = 0;
    PngSampler.resetAllCaches();
    PngSampler.maxInFlight = 6;
    PngSampler.maxPendingCompute = 16;
    for (const a of ACTS) {
      this.log(`Act ${a}: sampling ${SRC_COUNT[a]} PNG frames…`);
      const cache = await this.computeAct(a);
      const bytes = this.encodeAndCheck(a, cache);
      this.generated[a] = { cache, bytes };
      this.log(`Act ${a}: exact round trip passed, ${(bytes.length / 1048576).toFixed(2)} MiB.`);
    }
  },

  async generate() {
    if (this.running) return;
    this.setBusy(true);
    this.clearDownloads();
    this.logEl.textContent = "";
    try {
      await this.prepareSources();
      for (const a of ACTS) this.addDownload(a, this.generated[a].bytes);
      this.log("All acts passed. Download all four files and place them in assets/density/.");
    } catch (error) {
      this.log(`ERROR: ${error && error.stack ? error.stack : error}`);
    } finally {
      this.setBusy(false);
    }
  },

  async verify() {
    if (this.running) return;
    this.setBusy(true);
    this.clearDownloads();
    this.logEl.textContent = "";
    try {
      await this.prepareSources();
      for (const a of ACTS) {
        const response = await fetch(`../assets/density/act${a}.swd?v=${Date.now()}`, { cache: "no-store" });
        if (!response.ok) throw new Error(`Act ${a}: committed file returned HTTP ${response.status}`);
        const decoded = DensityCodec.decodeAct(await response.arrayBuffer(), {
          act: a, frameCount: SRC_COUNT[a], cols: COLS, rows: ROWS, particleCount: N,
        });
        const source = this.generated[a].cache;
        let rasterHashDifferences = 0;
        for (let frame = 0; frame < source.cycle; frame++) {
          if (decoded.pixHash[frame] !== source.pixHash[frame]) rasterHashDifferences++;
          if (decoded.countsHash[frame] !== source.countsHash[frame]) {
            throw new Error(`Act ${a}: count hash mismatch at frame ${frame} (file ${decoded.countsHash[frame]}, current ${source.countsHash[frame]})`);
          }
        }
        for (let i = 0; i < source.counts.length; i++) {
          if (decoded.counts[i] !== source.counts[i]) {
            const frame = Math.floor(i / CELLS);
            const cell = i % CELLS;
            throw new Error(`Act ${a}: committed count mismatch at frame ${frame}, cell ${cell} (file ${decoded.counts[i]}, current ${source.counts[i]})`);
          }
        }
        this.log(`Act ${a}: committed file exactly matches all PNG-derived counts (${rasterHashDifferences} diagnostic raster hash differences).`);
      }
      this.log("Verification passed for all four committed files.");
    } catch (error) {
      this.log(`ERROR: ${error && error.stack ? error.stack : error}`);
    } finally {
      this.setBusy(false);
    }
  },
};

function setup() {
  noCanvas();
  COLS = GRID_COLS;
  ROWS = GRID_ROWS;
  CELLS = COLS * ROWS;
  PngSampler.init();
  PngSampler.realloc();
  PngSampler.assetBase = "../";
  document.getElementById("generate").disabled = false;
  document.getElementById("verify").disabled = false;
  document.getElementById("generate").addEventListener("click", () => DensityGenerator.generate());
  document.getElementById("verify").addEventListener("click", () => DensityGenerator.verify());
  DensityGenerator.logEl.textContent = `Ready. ${DensityGenerator.totalFrames} source frames found in configuration.\n`;
}

function draw() {}
