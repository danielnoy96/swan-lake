function estimateMismatch() {
  // Sum of surplus particles over all cells; higher means more redistribution needed
  let surplus = 0;
  for (let i = 0; i < CELLS; i++) {
    const d = (Particles.cellCounts[i] | 0) - (Particles.desired[i] | 0);
    if (d > 0) surplus += d;
  }
  return surplus;
}

function estimateMeanCellDist() {
  const samples = min(220, N);
  const stride = max(1, (N / samples) | 0);
  let sum = 0;
  for (let s = 0; s < samples; s++) {
    const i = s * stride;
    const cell = Particles.cell[i];
    const ox = cellOriginX(cell), oy = cellOriginY(cell);
    const tx = ox + Particles.u[i] * cellW;
    const ty = oy + Particles.v[i] * cellH;
    const dx = Particles.x[i] - tx, dy = Particles.y[i] - ty;
    sum += sqrt(dx * dx + dy * dy);
  }
  return sum / samples;
}

const _hotCells = new Int32Array(48);
const _hotNeeds = new Int32Array(48);
let _hotCount = 0;

function computeDeficitHotspots() {
  // Top-k deficit cells by (desired - current). O(CELLS * K) with small K.
  const K = _hotCells.length;
  let count = 0;
  for (let i = 0; i < K; i++) { _hotCells[i] = -1; _hotNeeds[i] = 0; }
  for (let cell = 0; cell < CELLS; cell++) {
    const need = (Particles.desired[cell] | 0) - (Particles.cellCounts[cell] | 0);
    if (need <= 0) continue;
    let j = count < K ? count++ : K - 1;
    if (count === K && need <= _hotNeeds[j]) continue;
    while (j > 0 && need > _hotNeeds[j - 1]) {
      _hotNeeds[j] = _hotNeeds[j - 1];
      _hotCells[j] = _hotCells[j - 1];
      j--;
    }
    _hotNeeds[j] = need;
    _hotCells[j] = cell;
  }
  _hotCount = count;
}

function catchUpNudge(strength01) {
  if (_hotCount <= 0) return;
  const cu = constrain(strength01, 0, 1);
  // Gentle "hard catch-up" (no teleport): keep this subtle to avoid feeling too fast.
  const kick = 0.00025 * (0.20 + cu * 0.60);
  const stride = 10; // much fewer nudged particles per update
  const seed = (t * 60) | 0;
  for (let i = 0; i < N; i += stride) {
    const from = Particles.cell[i];
    const surplus = (Particles.cellCounts[from] | 0) - (Particles.desired[from] | 0);
    if (surplus <= 0) continue;
    const idx = (i + seed) % _hotCount;
    const targetCell = _hotCells[idx];
    const tx = cellOriginX(targetCell) + cellW * 0.5;
    const ty = cellOriginY(targetCell) + cellH * 0.5;
    const dx = tx - Particles.x[i];
    const dy = ty - Particles.y[i];
    Particles.vx[i] += dx * kick;
    Particles.vy[i] += dy * kick;
  }
}

