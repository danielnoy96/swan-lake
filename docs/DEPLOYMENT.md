# Deployment guide

## Static hosting

The project can be published on any static HTTPS host, including GitHub Pages. Upload the repository contents with `index.html` at the site root. No build output is required.

## GitHub Pages

1. Push the repository to GitHub.
2. Open the repository's **Settings → Pages**.
3. Select the branch and root folder that contain `index.html`.
4. Save and wait for the published HTTPS URL.
5. Open the site, allow microphone access, and run the checks below.

## Cache busting

`index.html` appends a version query to every local JavaScript file:

```html
<script src="core.js?v=20260211-25"></script>
```

When deploying JavaScript changes, replace the version value on all local script tags with a new shared value. A date plus revision number works well. This helps browsers and GitHub Pages avoid serving stale JavaScript.

When density files change, also update `DENSITY_REVISION` in `core.js`. Runtime requests include that value as a query parameter so a previously cached act cannot survive a data update.

## Deployment checks

- The page is served over HTTPS.
- p5.js and p5.sound load successfully from cdnjs.
- All typography requests and all four `assets/density/actN.swd` requests succeed.
- Normal mode makes no `assets/actN/*.png` requests.
- The loading sequence completes.
- Click or tap triggers the browser microphone prompt.
- Sound-driven playback responds after permission is granted.
- Automatic mode works with `A`.
- A complete four-act loop plays without an error overlay.
- The layout fills desktop and mobile screens.
- A hard refresh loads the latest JavaScript version.
- `?legacySampler=1` still reaches the PNG-backed recovery path.

## Offline exhibition

The current project depends on cdnjs for p5.js and p5.sound. For an offline or unreliable-network installation, download the exact library versions, serve them locally, update the two CDN script tags, and test microphone behavior in the final browser environment.

## Microphone security

Browsers generally allow microphone access only on HTTPS origins or `localhost`, and only after a user gesture. Do not attempt to start the microphone automatically on page load.
