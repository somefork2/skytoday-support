# Sky Today — support site and data feed

GitHub Pages: support (`index.html`), privacy policy (`privacy.html`) and `feed.json`, the app's daily data:
bright comets (JPL SBDB/Horizons elements), asteroid close approaches (JPL CNEOS), satellite TLEs (CelesTrak:
ISS, Tiangong, Hubble), aurora Kp forecast (NOAA SWPC) and curated `notes.json`.
`.github/workflows/feed.yml` rebuilds the feed every 6 hours with `tools/build_feed.py`.
Add a note: edit `notes.json` (`[{"id","start","end","title":{"en","ru"},"body":{"en","ru"}}]`) and push.
