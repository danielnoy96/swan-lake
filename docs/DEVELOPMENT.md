# Development guide

## Requirements

The artwork itself requires only a modern browser. p5.js and p5.sound 1.9.0 are loaded from cdnjs by `index.html`.

For local development, use any static HTTP server. For example:

```powershell
python -m http.server 8000
```

Open <http://localhost:8000>. Do not rely on opening `index.html` through `file://`; browser asset and microphone policies can differ.

## Editing model

There is no compiler or bundler. Refresh the browser after changing JavaScript. Scripts communicate through shared global names, so their load order in `index.html` matters.

When deploying JavaScript changes, update the `?v=` query strings in `index.html` to prevent stale cached files.

## Important configuration

The main tuning constants live near the top of `core.js`:

| Constant | Current value | Meaning |
| --- | ---: | --- |
| `SRC_COUNT` | 192 / 172 / 192 / 136 | Number of source frames in each act |
| `FPS_EFFECTIVE` | 14 for every act | Frames per internal sound-second |
| `FRAME_STEP` | 1 for every act | Uses every source frame |
| `GRID_COLS` × `GRID_ROWS` | 137 × 74 | Fixed density simulation grid |
| `N` | 4,200 | Main particle population |
| `TRANSITION_DURATION` | 2.6 | Internal seconds per act transition |
| `MIC_THRESHOLD` | 0.03 | Sound threshold used by playback logic |
| `AUTO_SPEED` | 1.0 | Internal seconds per real second |
| `SIM_SEED` | 1337 | Deterministic random and noise seed |
| `DENSITY_REVISION` | `20260825-1` | Cache key for generated `.swd` files |

Typography counts and sizes are configured near the top of `typography.js`. Act colors and tint placement are configured in `style.js`.

## Debugging

- Press `D` for runtime timing, cache, grid, frame, and particle statistics.
- Press `G` to see grid alignment.
- Press `C` to inspect and interpolate individual source frames.
- Press `P` to log a reproducible runtime snapshot and hashes to the console.
- Add `?legacySampler=1` to exercise the original PNG-backed sampler.
- Runtime errors and unhandled promise rejections are displayed in an on-screen overlay.

See [CONTROLS.md](CONTROLS.md) for all parameters and keys.

## Manual verification

After behavioral or rendering changes:

1. Confirm that the loading image resolves through particles.
2. Start automatic mode and observe a complete act and transition.
3. Test microphone input after an explicit click or tap.
4. Resize the browser and test a narrow/mobile aspect ratio.
5. Confirm that only one visible canvas exists.
6. Check comparison mode on each act.
7. Check the console and on-screen error overlay.

## Performance notes

- Startup decodes only Act 1 (or the requested comparison act) before showing the start screen. Remaining acts load in the background.
- The four SWD files total approximately 1.25 MB; normal playback does not request the roughly 75 MB PNG sequence.
- The fixed grid keeps cached density data stable when the viewport changes.
- Source images and typography masks are sampled off-screen; act frames are not directly composited onto the visible canvas.
- High-resolution displays increase canvas cost. Test changes on the actual presentation hardware.

## Density generation

Open `/tools/generate-density.html` through the local HTTP server whenever source frames, grid dimensions, particle count, or sampling rules change. Generate all four files, place them in `assets/density/`, then run the page's committed-file verification before updating `DENSITY_REVISION`.
