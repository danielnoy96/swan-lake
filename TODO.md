# Project tasks

This list is intentionally limited to work that can improve presentation, maintainability, or confidence without changing the artistic concept.

## Before presenting or submitting

- [ ] Test microphone permission on the final HTTPS URL.
- [ ] Test automatic mode with `A` when microphone access is unavailable.
- [ ] Watch a complete four-act loop and check every transition.
- [ ] Test at least one desktop and one mobile viewport.
- [ ] Confirm that all 692 frame images load without 404 errors.
- [ ] Check the browser console for runtime errors.
- [ ] Record the final deployment URL in `README.md`.
- [ ] Add author, course, year, and asset credits to `README.md` when ready.

## Possible improvements

- [ ] Decide whether p5.js should be stored locally for offline exhibition use.
- [ ] Measure startup time and memory use on the intended exhibition device.
- [x] Replace the runtime PNG sequence with compact, precomputed per-act density files.
- [ ] Add a small automated asset audit for frame names, counts, and missing files.
- [ ] Add a favicon and social-preview metadata.
- [ ] Document the source and license of every photographic asset.

## Completed

- [x] Organize the project documentation.
- [x] Provide microphone and automatic playback modes.
- [x] Provide debug and frame-comparison controls.
- [x] Add a visual loading sequence while act frames are prepared.
