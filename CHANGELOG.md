# Changelog

Meaningful project changes are recorded here. Dates use `YYYY-MM-DD`.

## Unreleased

### Added

- Project documentation covering the concept, architecture, controls, assets, development, and deployment.
- Precomputed SWD density assets and a reusable browser generator/verifier.

### Changed

- Normal playback now loads compact per-act density data instead of sampling all 692 PNG frames at startup.
- Later acts load in the background, with an explicit `?legacySampler=1` PNG recovery mode.

## 2026-02-12

### Changed

- Final project revision.

## 2026-02-11

### Added

- Particle-based loading screen.

### Changed

- Adjustments for consistent appearance across local and hosted environments.

## 2026-02-10

### Changed

- Final visual and behavioral refinements.
