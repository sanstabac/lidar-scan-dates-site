# When Was This Spot Scanned?

A static web page that shows when each spot in a Washington DNR lidar survey was captured.
Click anywhere on the map to see:

- which day or days the spot was scanned
- every flight pass over that 6 × 6 ft square, with the time it passed over
- the nearest individual laser point, with its exact timestamp, what it hit, and its elevation

Current data: a block of about 6 × 6 miles (15 × 15 tiles of 2,250 ft) from the
**Cascades North Wali 2023** project on the [WA DNR Lidar Portal](https://lidarportal.dnr.wa.gov/),
centered near Punkin Center on SR 20. Every point in every tile was decoded.
Despite the project name, this area was flown in 2022.

## Publish on GitHub Pages

1. Create a new repository and push the contents of this folder to it (`index.html` at the root).
   The data is about 0.7 GB in roughly 700 files, all under GitHub's 100 MB per-file limit.
2. In the repository, open **Settings → Pages**.
3. Under **Build and deployment**, set **Source** to *Deploy from a branch*. Pick `main` and `/ (root)`, then save.
4. After a few minutes the site is live at `https://<user-or-org>.github.io/<repo>/`.

GitHub Pages allows sites up to 1 GB. This one fits, but a much larger area won't; see "Bigger areas" below.

To preview locally, serve the folder over HTTP. Opening `index.html` straight from disk won't work,
because browsers block `fetch()` of local files.

```bash
python -m http.server 8000
```

Then open http://localhost:8000.

## How the page works

- **Zoomed out**, it shows `data/overview_hillshade.png` (15 ft pixels) with colors from
  `data/overview.bin.gz` (main scan day, main flight pass and number of days, per 30 ft square).
- **Zoomed in**, it loads each tile in view from `data/tiles/<tile>/`:
  - `hillshade.png`: shaded relief at 3 ft pixels
  - `cells.bin.gz`: for each 6 ft square and flight pass, the point count and typical capture time
  - `points.bin.gz`: a random sample of 150,000 individual points with exact timestamps
- **Clicking** loads the tile under the click (if it isn't loaded yet) and fills the side panel.
- `data/index.json` lists the tiles, flight passes, scan days, the lat/lon conversion and the highlighted pin.

External resources: Leaflet 1.9.4 from cdnjs and IBM Plex fonts from Google Fonts.
If your network blocks those, download `leaflet.js` into the repo and change the `<script src>`.
The fonts fall back to system fonts on their own.

## Rebuild or change the area

Requires Python 3.10+:

```bash
pip install -r tools/requirements.txt
```

1. Edit the settings at the top of `tools/build_area.py`:
   - `PIN_LATLON`: the center and the highlighted spot
   - `RADIUS_FT`: how far out to go (15,840 ft = 3 miles)
   - `SRC_CRS`: must match the tiles (EPSG:2927 for this project)
2. Get the point cloud. On the DNR portal, use the rectangle tool to draw the area, tick only
   **Point Cloud** in the popup, and download. Put the zip or zips in a folder, for example `incoming/`.
   Zips holding several tiles are fine.
3. Decode the tiles. This deletes each zip once it's processed, so keep a copy if you need one:

   ```bash
   python tools/build_area.py process incoming work data
   ```

4. Build the overview and index:

   ```bash
   python tools/build_area.py assemble work data
   ```

5. Update the project name in the header text of `index.html` if you changed projects.

`python tools/build_area.py tiles work` lists the tile names and download boxes for the configured area.

Notes:

- **Timestamp format.** Tiles must use "Adjusted Standard GPS time" (LAS global encoding bit 0 set).
  Modern DNR projects do. Some older projects (around 2006) store only seconds-of-week,
  which can't be turned into a date. The script skips those tiles and says so in `work/process.log`.
- **Leap seconds.** `LEAP = 18` is the GPS-to-UTC offset, which is right for any data after 2017.
- **Speed.** Each tile takes about 6–15 seconds. The script runs 6 tiles at a time.

## Bigger areas

Each tile adds about 3 MB of web data, so GitHub Pages' 1 GB limit allows roughly 300 tiles
(about a 7 × 7 mile block). For more, either host the files elsewhere (any static web server or
cloud storage bucket works) or lower `N_SAMPLE` and raise `CELL` in `build_area.py` to make each
tile lighter.

## How the times are worked out

Every lidar point stores the GPS time of its laser pulse, in seconds since 1980-01-06 minus 10⁹.
The page adds 10⁹ back, subtracts the 18 leap seconds between GPS time and UTC, and shows the result
in Pacific time (PDT or PST, whichever applied that day). Times for a square are the average time for that
pass within the square (0.1 s resolution). Times for an individual point are exact to the millisecond.
