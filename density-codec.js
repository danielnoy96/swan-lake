// ---- SWD1 codec: deterministic sparse particle-density data ----
const DensityCodec = {
  MAGIC: [0x53, 0x57, 0x44, 0x31], // "SWD1"
  FORMAT_VERSION: 1,
  HEADER_SIZE: 32,
  RECORD_SIZE: 16,
  _crcTable: null,

  _asBytes(input) {
    if (input instanceof Uint8Array) return input;
    if (input instanceof ArrayBuffer) return new Uint8Array(input);
    if (ArrayBuffer.isView(input)) return new Uint8Array(input.buffer, input.byteOffset, input.byteLength);
    throw new Error("SWD input must be an ArrayBuffer or Uint8Array");
  },

  _crc32(bytes) {
    if (!this._crcTable) {
      const table = new Uint32Array(256);
      for (let n = 0; n < 256; n++) {
        let c = n;
        for (let k = 0; k < 8; k++) c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
        table[n] = c >>> 0;
      }
      this._crcTable = table;
    }
    let crc = 0xffffffff;
    const table = this._crcTable;
    for (let i = 0; i < bytes.length; i++) crc = table[(crc ^ bytes[i]) & 0xff] ^ (crc >>> 8);
    return (crc ^ 0xffffffff) >>> 0;
  },

  _countsHash(arr, cells) {
    let h = 2166136261 >>> 0;
    const stride = Math.max(1, Math.floor(cells / 512));
    for (let i = 0; i < arr.length; i += stride) {
      h ^= arr[i] & 0xffff;
      h = Math.imul(h, 16777619) >>> 0;
    }
    return h >>> 0;
  },

  _pushUleb(out, value) {
    let v = value >>> 0;
    do {
      let b = v & 0x7f;
      v >>>= 7;
      if (v) b |= 0x80;
      out.push(b);
    } while (v);
  },

  _readUleb(bytes, cursor, end) {
    let value = 0;
    let shift = 0;
    while (cursor.i < end && shift <= 28) {
      const b = bytes[cursor.i++];
      value |= (b & 0x7f) << shift;
      if ((b & 0x80) === 0) return value >>> 0;
      shift += 7;
    }
    throw new Error("Invalid or truncated ULEB128 value");
  },

  encodeAct(source) {
    const act = source.act | 0;
    const frameCount = source.frameCount | 0;
    const cols = source.cols | 0;
    const rows = source.rows | 0;
    const particleCount = source.particleCount | 0;
    const cells = cols * rows;
    const counts = source.counts;
    if (act < 1 || frameCount <= 0 || cols <= 0 || rows <= 0 || particleCount <= 0) {
      throw new Error("Invalid SWD act metadata");
    }
    if (!(counts instanceof Uint16Array) || counts.length !== frameCount * cells) {
      throw new Error("SWD counts length does not match frameCount * grid cells");
    }

    const payload = [];
    const payloadOffsets = new Uint32Array(frameCount);
    for (let frame = 0; frame < frameCount; frame++) {
      payloadOffsets[frame] = payload.length >>> 0;
      const off = frame * cells;
      let previous = -1;
      let sum = 0;
      for (let cell = 0; cell < cells; cell++) {
        const count = counts[off + cell] | 0;
        if (count <= 0) continue;
        this._pushUleb(payload, cell - previous);
        this._pushUleb(payload, count);
        previous = cell;
        sum += count;
      }
      if (sum !== particleCount) throw new Error(`Act ${act} frame ${frame} sums to ${sum}, expected ${particleCount}`);
    }

    const tableOffset = this.HEADER_SIZE;
    const payloadOffset = tableOffset + frameCount * this.RECORD_SIZE;
    const out = new Uint8Array(payloadOffset + payload.length);
    const view = new DataView(out.buffer);
    for (let i = 0; i < 4; i++) out[i] = this.MAGIC[i];
    view.setUint16(4, this.FORMAT_VERSION, true);
    view.setUint16(6, act, true);
    view.setUint16(8, frameCount, true);
    view.setUint16(10, cols, true);
    view.setUint16(12, rows, true);
    view.setUint16(14, particleCount, true);
    view.setUint32(16, tableOffset, true);
    view.setUint32(20, payloadOffset, true);
    view.setUint32(24, payload.length, true);

    const pixHash = source.pixHash || new Uint32Array(frameCount);
    const countsHash = source.countsHash || new Uint32Array(frameCount);
    const imgW = source.imgW || new Uint16Array(frameCount);
    const imgH = source.imgH || new Uint16Array(frameCount);
    for (let frame = 0; frame < frameCount; frame++) {
      const rec = tableOffset + frame * this.RECORD_SIZE;
      const frameCounts = counts.subarray(frame * cells, (frame + 1) * cells);
      view.setUint32(rec, payloadOffsets[frame], true);
      view.setUint32(rec + 4, (pixHash[frame] || 0) >>> 0, true);
      view.setUint32(rec + 8, (countsHash[frame] || this._countsHash(frameCounts, cells)) >>> 0, true);
      view.setUint16(rec + 12, (imgW[frame] || 0) & 0xffff, true);
      view.setUint16(rec + 14, (imgH[frame] || 0) & 0xffff, true);
    }
    out.set(payload, payloadOffset);
    view.setUint32(28, this._crc32(out.subarray(this.HEADER_SIZE)), true);
    return out;
  },

  decodeAct(input, expected) {
    const bytes = this._asBytes(input);
    if (bytes.byteLength < this.HEADER_SIZE) throw new Error("SWD file is shorter than its header");
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    for (let i = 0; i < 4; i++) if (bytes[i] !== this.MAGIC[i]) throw new Error("Invalid SWD magic");
    const version = view.getUint16(4, true);
    const act = view.getUint16(6, true);
    const frameCount = view.getUint16(8, true);
    const cols = view.getUint16(10, true);
    const rows = view.getUint16(12, true);
    const particleCount = view.getUint16(14, true);
    const tableOffset = view.getUint32(16, true);
    const payloadOffset = view.getUint32(20, true);
    const payloadLength = view.getUint32(24, true);
    const storedCrc = view.getUint32(28, true) >>> 0;
    if (version !== this.FORMAT_VERSION) throw new Error(`Unsupported SWD version ${version}`);
    if (frameCount <= 0 || cols <= 0 || rows <= 0 || particleCount <= 0) throw new Error("Invalid SWD dimensions");
    if (tableOffset !== this.HEADER_SIZE || payloadOffset !== tableOffset + frameCount * this.RECORD_SIZE) {
      throw new Error("Invalid SWD table layout");
    }
    if (payloadOffset + payloadLength !== bytes.byteLength) throw new Error("SWD payload length does not match file size");
    if (this._crc32(bytes.subarray(this.HEADER_SIZE)) !== storedCrc) throw new Error("SWD CRC32 mismatch");
    if (expected) {
      if ((expected.act | 0) !== act) throw new Error(`SWD act ${act} does not match requested act ${expected.act}`);
      if ((expected.frameCount | 0) !== frameCount) throw new Error(`SWD frame count ${frameCount} does not match ${expected.frameCount}`);
      if ((expected.cols | 0) !== cols || (expected.rows | 0) !== rows) throw new Error("SWD grid dimensions do not match runtime constants");
      if ((expected.particleCount | 0) !== particleCount) throw new Error("SWD particle count does not match runtime constant");
    }

    const cells = cols * rows;
    const counts = new Uint16Array(frameCount * cells);
    const pixHash = new Uint32Array(frameCount);
    const countsHash = new Uint32Array(frameCount);
    const imgW = new Uint16Array(frameCount);
    const imgH = new Uint16Array(frameCount);
    const starts = new Uint32Array(frameCount);
    for (let frame = 0; frame < frameCount; frame++) {
      const rec = tableOffset + frame * this.RECORD_SIZE;
      starts[frame] = view.getUint32(rec, true);
      pixHash[frame] = view.getUint32(rec + 4, true);
      countsHash[frame] = view.getUint32(rec + 8, true);
      imgW[frame] = view.getUint16(rec + 12, true);
      imgH[frame] = view.getUint16(rec + 14, true);
      if (starts[frame] > payloadLength || (frame > 0 && starts[frame] < starts[frame - 1])) {
        throw new Error(`Invalid payload offset for frame ${frame}`);
      }
    }
    if (starts[0] !== 0) throw new Error("First SWD frame does not begin at payload offset zero");

    for (let frame = 0; frame < frameCount; frame++) {
      const start = payloadOffset + starts[frame];
      const end = payloadOffset + (frame + 1 < frameCount ? starts[frame + 1] : payloadLength);
      const cursor = { i: start };
      const outOff = frame * cells;
      let previous = -1;
      let sum = 0;
      while (cursor.i < end) {
        const delta = this._readUleb(bytes, cursor, end);
        const count = this._readUleb(bytes, cursor, end);
        if (delta === 0 || count === 0) throw new Error(`Invalid sparse pair in frame ${frame}`);
        const cell = previous + delta;
        if (cell <= previous || cell >= cells) throw new Error(`Cell index outside grid in frame ${frame}`);
        if (count > particleCount || sum + count > particleCount) throw new Error(`Particle overflow in frame ${frame}`);
        counts[outOff + cell] = count;
        previous = cell;
        sum += count;
      }
      if (cursor.i !== end) throw new Error(`Frame ${frame} ends inside a sparse pair`);
      if (sum !== particleCount) throw new Error(`Frame ${frame} sums to ${sum}, expected ${particleCount}`);
      const hash = this._countsHash(counts.subarray(outOff, outOff + cells), cells);
      if (countsHash[frame] && hash !== countsHash[frame]) throw new Error(`Counts hash mismatch in frame ${frame}`);
      countsHash[frame] = hash;
    }

    const done = new Uint8Array(frameCount); done.fill(1);
    return {
      act, cycle: frameCount, cols, rows, particleCount, counts,
      pixHash, countsHash, imgW, imgH, done,
      loading: new Uint8Array(frameCount), failed: new Uint8Array(frameCount),
    };
  },
};
