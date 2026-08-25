# Controls

## Audience controls

| Input | Result |
| --- | --- |
| Click or tap | Requests microphone access and starts sound-driven playback. |
| `A` | Toggles automatic playback. |

On mobile, the initial tap is required by the browser before microphone audio can begin.

## Development controls

| Key | Result |
| --- | --- |
| `D` | Toggles the debug information overlay. |
| `G` | Toggles the simulation grid. |
| `C` | Toggles frame-comparison mode and updates the URL. |
| `P` | Writes a detailed runtime snapshot to the browser console. |

## Comparison mode

While comparison mode is active:

| Input | Result |
| --- | --- |
| `1`–`4` | Selects an act. |
| Left / Right | Moves backward or forward by one source frame. |
| Shift + Left / Right | Moves backward or forward by ten source frames. |
| Up / Down | Changes interpolation alpha by 0.02. |
| Shift + Up / Down | Changes interpolation alpha by 0.10. |

Comparison mode stores its current values in the address bar, making a state easy to share or reproduce.

## URL parameters

| Parameter | Accepted values | Purpose |
| --- | --- | --- |
| `compare` or `cmp` | `1`, `true` | Opens comparison mode. |
| `act` | `1`–`4` | Selects the comparison act. |
| `src` or `src0` | Integer | Selects the source frame. |
| `alpha` or `a` | `0`–`1` | Sets interpolation between adjacent frames. |
| `debug` or `d` | `1`, `true` | Opens the debug overlay. |
| `legacySampler` | `1` | Uses the original PNG sampler for verification or recovery. |

Example:

```text
http://localhost:8000/?compare=1&act=3&src=96&alpha=0.5&debug=1
```
