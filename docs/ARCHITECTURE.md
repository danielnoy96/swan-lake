# Architecture

## Runtime flow

```text
Precomputed SWD act data
      ↓
Sampler: sparse data → density per grid cell
      ↓
Particles: desired counts → movement and redistribution
      ↓
Style + Typography
      ↓
p5 canvas
```

`sketch.js` coordinates this flow through p5's `preload`, `setup`, `draw`, resize, and input callbacks.

## Script loading order

The order in `index.html` is significant because the project uses shared browser globals:

1. `core.js`
2. `style.js`
3. `sampler.js`
4. `density-codec.js`
5. `density-sampler.js`
6. `particles.js`
7. `redistribute.js`
8. `acts.js`
9. `typography.js`
10. `sketch.js`

The code is intentionally split by responsibility but is not packaged as ES modules.

## Major systems

### Core configuration

`core.js` defines act IDs, source counts, playback rates, particle count, the fixed simulation grid, microphone thresholds, visual constants, shared runtime state, and error reporting. It also compensates for device-pixel ratio and browser zoom differences.

### Frame sampler

Normal playback fetches one compact `assets/density/actN.swd` file per act. `density-sampler.js` validates and decodes its sparse cell/count pairs into the same dense 137 × 74 frame cache used by the particle system. Act 1 (or a comparison deep-link act) gates startup; later acts load sequentially in the background.

`sampler.js` retains the PNG source sampler. The developer generator uses it to convert brightness and alpha into exactly 4,200 deterministic particle assignments per frame. It is also available at runtime only through `?legacySampler=1`.

### Particle system

`particles.js` owns particle positions, velocity, grid occupancy, retargeting, transition paths, separation, and rendering. During an act, particles follow interpolated density targets from adjacent source frames. During a transition, they leave the silhouette and travel through full-canvas paths before settling into the next act.

### Redistribution

`redistribute.js` identifies cells that contain too few particles and supports balancing particles toward density deficits. This helps silhouettes resolve without obvious holes or long-lived stacks.

### Act state machine

`acts.js` alternates between `ACT` and `TRANSITION` states. Each act plays one source cycle. A transition lasts 2.6 internal seconds, after which the next act begins. Act four loops back to act one.

### Styling

`style.js` blends each act's background, highlight, and secondary colors. Microphone level controls how strongly those colors appear. It also manages focused tint regions and organic stains.

### Typography

`typography.js` samples two authored masks and renders them with a separate 9,000-particle system. The large and small typography groups have independent sizes and particle allocations.

### Application coordinator

`sketch.js` manages:

- primary-act density loading, background act loading, and the particle-reveal loading screen;
- canvas creation and resize handling;
- microphone and automatic timing;
- act sampling, interpolation, and particle updates;
- typography rendering;
- keyboard, pointer, and URL-based debug controls.

## Timing model

The animation uses internal time `t`, measured in sound-seconds. Automatic mode advances it at one internal second per real second. Microphone mode derives its rate from the current input level. Silence can therefore freeze act progression without stopping the browser render loop.

At the configured 14 frames per internal second, the source cycles and four 2.6-second transitions produce a complete loop of roughly one minute.
