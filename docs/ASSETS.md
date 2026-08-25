# Asset guide

## Inventory

| Folder | Files | Approximate size | Purpose |
| --- | ---: | ---: | --- |
| `assets/act1/` | 192 | 13.9 MB | Act 1 animation frames |
| `assets/act2/` | 172 | 16.9 MB | Act 2 animation frames |
| `assets/act3/` | 192 | 25.5 MB | Act 3 animation frames |
| `assets/act4/` | 136 | 18.0 MB | Act 4 animation frames |
| `assets/loading/` | 1 | 0.2 MB | Particle-reveal loading image |
| `assets/typography/` | 2 | 0.4 MB | Large and small typography masks |
| `assets/_debug/` | 2 | less than 0.1 MB | Typography mask diagnostics |
| `assets/density/` | 4 | 1.25 MB | Compact runtime particle-density data |

The four act folders contain 692 PNG frames in total.

The PNG sequences are authored sources and are not requested during normal playback. The browser loads `act1.swd` first, then fetches the remaining act files in the background. The PNG source path remains available through `?legacySampler=1`.

## Regenerating runtime density data

1. Serve the repository over HTTP.
2. Open `/tools/generate-density.html`.
3. Select **Generate and validate** and wait for all 692 frames to pass.
4. Download all four `.swd` files into `assets/density/`.
5. Select **Verify committed files against PNGs** and require an exact match for every frame.
6. Bump `DENSITY_REVISION` in `core.js` and the JavaScript cache version in `index.html`.

The generator validates the SWD header, CRC32, grid and particle constants, frame sums, hashes, and every decoded cell count before enabling downloads.

## Frame naming

Act frames must follow this exact pattern:

```text
assets/act{act}/act{act}_{frame}.png
```

The frame index is zero-based and padded to five digits. Examples:

```text
assets/act1/act1_00000.png
assets/act1/act1_00191.png
assets/act4/act4_00135.png
```

Frames for an act must be contiguous. If frames are added or removed, update `SRC_COUNT` in `core.js`.

## Typography

The typography masks are:

- `assets/typography/bigtypography.png`
- `assets/typography/smalltypography.png`

They are sampled as masks rather than drawn directly. Preserve clear contrast and transparency when replacing them. Very large sources are downscaled at runtime before mask generation.

## Loading image

`assets/loading/loadingscreen.jpg` is converted into a temporary particle-density target during startup. The loader detects whether the image uses dark marks on a light background and can invert it automatically.

## Replacing an act sequence

1. Export every frame as a PNG with the same canvas dimensions.
2. Start numbering at `00000` and avoid gaps.
3. Place the sequence in the corresponding act folder.
4. Update `SRC_COUNT` in `core.js`.
5. If the intended playback rate changes, update `FPS_EFFECTIVE`.
6. Run a local server and watch the complete sequence in comparison and automatic modes.
7. Update the cache-busting version in `index.html` before deployment.
8. Regenerate and verify all four density files before deployment.

## Source and licensing

Asset authorship and license information is not currently recorded in the repository. Add those details before public distribution or exhibition if any source material is not wholly original.
