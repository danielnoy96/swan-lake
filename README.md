# Swan Lake Card

Swan Lake Card is a full-screen, sound-reactive browser artwork made with p5.js. A fixed population of particles reconstructs animated image sequences inspired by *Swan Lake*, then travels across the screen while changing from one act into the next.

The work can run from microphone input or in an automatic presentation mode. It is designed for desktop and mobile browsers and does not require a build step.

## Run locally

The project loads p5.js from a CDN, so an internet connection is required unless those libraries are hosted locally. Serve the folder through a local web server rather than opening `index.html` directly.

```powershell
python -m http.server 8000
```

Then open <http://localhost:8000>.

## Interaction

- Tap or click to grant microphone access and begin sound-driven playback.
- Press `A` to toggle automatic playback.
- Louder input advances the animation and reveals more color.
- The four acts loop continuously, with particle-flight transitions between them.

See [docs/CONTROLS.md](docs/CONTROLS.md) for development and comparison controls.

## Project structure

| Path | Responsibility |
| --- | --- |
| `index.html` | Loads p5.js, p5.sound, and the project scripts in dependency order. |
| `core.js` | Global configuration, runtime state, grid sizing, and shared helpers. |
| `style.js` | Per-act palettes and sound-responsive color rendering. |
| `sampler.js` | PNG source sampler used by the density generator and explicit legacy mode. |
| `density-codec.js` / `density-sampler.js` | Decode compact precomputed density files for normal playback. |
| `particles.js` | Particle storage, movement, transitions, and drawing. |
| `redistribute.js` | Density balancing and deficit-hotspot helpers. |
| `acts.js` | Four-act timeline and transition state machine. |
| `typography.js` | Separate particle system for typography masks. |
| `sketch.js` | p5 lifecycle, preload sequence, input, playback, and debug UI. |
| `assets/` | Authored PNG sources, runtime density files, loading image, typography masks, and debug images. |
| `tools/generate-density.html` | Regenerates and verifies runtime density files from the authored PNGs. |

## Key technical facts

- p5.js and p5.sound 1.9.0
- 4,200 main particles
- 9,000 typography particles
- Fixed 137 × 74 simulation grid
- 692 act frames across four sequences
- Approximately 1.25 MB of precomputed runtime density data across four act files
- 14 source frames per internal sound-second
- Approximately one minute per complete loop at automatic speed

## Documentation

- [Concept](docs/CONCEPT.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Controls](docs/CONTROLS.md)
- [Asset guide](docs/ASSETS.md)
- [Development guide](docs/DEVELOPMENT.md)
- [Deployment guide](docs/DEPLOYMENT.md)
- [Project tasks](TODO.md)
- [Change history](CHANGELOG.md)

## Browser permissions

Microphone input requires user permission and a secure context. `localhost` is accepted during development; a deployed version should use HTTPS. If p5.sound is unavailable or permission is not granted, use automatic mode with `A`.

Normal playback does not request the 692 authored act PNGs. Add `?legacySampler=1` to a URL only when testing or recovering through the original PNG sampler.
