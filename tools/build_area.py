"""Build the multi-tile data for index.html from WA DNR lidar tiles.

Three steps:

  1. python build_area.py tiles  <work_dir>
       Lists the 2,250 ft tiles within RADIUS_FT of PIN_LATLON and writes work_dir/tiles.json
       (each tile's name plus a tiny lat/lon box at its center, for the portal's download URL).

  2. python build_area.py process <incoming_dir> <work_dir> <site_data_dir> [--watch]
       Decodes every tile zip in incoming_dir (in parallel), writes the per-tile web files to
       site_data_dir/tiles/<name>/, keeps a small summary in work_dir, and deletes the zip.
       With --watch it keeps polling until a file named DONE appears in incoming_dir.

  3. python build_area.py assemble <work_dir> <site_data_dir>
       Combines the tile summaries into the overview layers and data/index.json.

Per tile (site_data_dir/tiles/<name>/):
  hillshade.png   shaded relief, 3 ft pixels (750 x 750)
  cells.bin.gz    per 6 ft square and flight pass: point count (u8) then typical capture time (u16, 0.1 s)
  points.bin.gz   random sample of individual points, 13-byte records (see SAMPLE_DTYPE)
"""
import gzip, io, json, sys, time, zipfile, zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta, timezone
from pathlib import Path
from zoneinfo import ZoneInfo

import numpy as np

SRC_CRS = "EPSG:2927"          # NAD83(HARN) / Washington South (US ft), per the tile headers
PIN_LATLON = (48.53100, -122.00675)
RADIUS_FT = 15_840             # 3 miles
TILE = 2250.0                  # DNR tile size, ft
PIX = 3.0                      # detail hillshade pixel, ft
CELL = 6.0                     # detail "when" square, ft
OV_PIX = 15.0                  # overview hillshade pixel, ft
OV_CELL = 30.0                 # overview date square, ft
N_SAMPLE = 150_000
WORKERS = 6
GPS_EPOCH = datetime(1980, 1, 6, tzinfo=timezone.utc)
LEAP = 18
PACIFIC = ZoneInfo("America/Los_Angeles")
SAMPLE_DTYPE = np.dtype([("x", "<u2"), ("y", "<u2"), ("z", "<i2"), ("ms", "<u4"),
                         ("p", "u1"), ("c", "u1"), ("r", "u1")])  # 13 bytes


def gps_to_utc(t):
    return GPS_EPOCH + timedelta(seconds=float(t) + 1e9 - LEAP)


def pacific_day(t):
    return int(gps_to_utc(t).astimezone(PACIFIC).strftime("%Y%m%d"))


def tile_origin(name):
    return float(name[1:8]), float(name[9:15]) * 10


# ---------------------------------------------------------------- step 1
def cmd_tiles(work):
    from pyproj import Transformer
    work = Path(work); work.mkdir(parents=True, exist_ok=True)
    inv = Transformer.from_crs("EPSG:4326", SRC_CRS, always_xy=True)
    fwd = Transformer.from_crs(SRC_CRS, "EPSG:4326", always_xy=True)
    px, py = inv.transform(PIN_LATLON[1], PIN_LATLON[0])
    ix0, ix1 = int((px - RADIUS_FT) // TILE), int((px + RADIUS_FT) // TILE)
    iy0, iy1 = int((py - RADIUS_FT) // TILE), int((py + RADIUS_FT) // TILE)
    tiles = []
    for iy in range(iy1, iy0 - 1, -1):
        for ix in range(ix0, ix1 + 1):
            x0, y0 = ix * TILE, iy * TILE
            lon, lat = fwd.transform(x0 + TILE / 2, y0 + TILE / 2)
            d = 0.0002
            box = [[round(lon - d, 5), round(lat - d, 5)], [round(lon + d, 5), round(lat - d, 5)],
                   [round(lon + d, 5), round(lat + d, 5)], [round(lon - d, 5), round(lat + d, 5)],
                   [round(lon - d, 5), round(lat - d, 5)]]
            tiles.append({"name": f"w{int(x0):07d}n{int(y0 // 10):06d}", "geojson": json.dumps({"type": "Polygon", "coordinates": [box]}, separators=(",", ":"))})
    (work / "tiles.json").write_text(json.dumps(tiles))
    print(f"{len(tiles)} tiles ({ix1 - ix0 + 1} x {iy1 - iy0 + 1}), pin at {px:.0f}, {py:.0f}")


# ---------------------------------------------------------------- step 2
def process_zip(zpath, work, site):
    import laspy
    from PIL import Image
    from scipy import ndimage

    zpath, work, site = Path(zpath), Path(work), Path(site)
    with zipfile.ZipFile(zpath) as z:
        lazs = [n for n in z.namelist() if n.lower().endswith(".laz")]
        if not lazs:
            return {"zip": zpath.name, "skipped": "no point cloud in this square"}
        results = []
        for entry in lazs:
            name = Path(entry).stem
            las = laspy.read(io.BytesIO(z.read(entry)))
            if not las.header.global_encoding.value & 1:
                results.append({"name": name, "skipped": "timestamps are GPS week time, not dates"})
                continue
            results.append(_process_las(name, las, work, site, ndimage, Image))
    return results


def _process_las(name, las, work, site, ndimage, Image):
    x0, y0 = tile_origin(name)
    x1, y1 = x0 + TILE, y0 + TILE
    x = np.asarray(las.x); y = np.asarray(las.y); z = np.asarray(las.z)
    t = np.asarray(las.gps_time)
    src = np.asarray(las.point_source_id)
    cls = np.asarray(las.classification)
    ret = np.asarray(las.return_number); nret = np.asarray(las.number_of_returns)
    good = (cls != 7) & (cls != 18) & (x >= x0) & (x < x1) & (y >= y0) & (y < y1)

    pass_ids = [int(s) for s in np.unique(src[good])]
    if not pass_ids:
        return {"name": name, "skipped": "no usable points"}
    pidx = np.minimum(np.searchsorted(pass_ids, src), len(pass_ids) - 1).astype(np.int64)
    passes = []
    for s in pass_ids:
        m = (src == s) & good
        ta, tb = float(t[m].min()), float(t[m].max())
        passes.append({"id": s, "gps_start": ta, "duration_s": round(tb - ta, 3), "points": int(m.sum()), "day": pacific_day(ta)})
    NP = len(passes)
    pstart = np.array([p["gps_start"] for p in passes])

    out = site / "tiles" / name
    out.mkdir(parents=True, exist_ok=True)

    # hillshade of the highest point per 3 ft pixel
    W = H = int(TILE / PIX)
    col = np.clip(((x - x0) / PIX).astype(int), 0, W - 1)
    row = np.clip(((y1 - y) / PIX).astype(int), 0, H - 1)
    dsm = np.full(H * W, -np.inf)
    np.maximum.at(dsm, (row * W + col)[good], z[good])
    dsm = dsm.reshape(H, W)
    holes = ~np.isfinite(dsm)
    if holes.all():
        dsm[:] = 0
    elif holes.any():
        _, (ri, ci) = ndimage.distance_transform_edt(holes, return_indices=True)
        dsm = dsm[ri, ci]
    dzdy, dzdx = np.gradient(dsm, PIX)
    az, alt = np.radians(315), np.radians(45)
    slope = np.arctan(np.hypot(dzdx, dzdy) * 1.5)
    aspect = np.arctan2(-dzdx, dzdy)
    shade = np.clip(np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect), 0, 1)
    shade8 = (40 + shade * 205).astype(np.uint8)
    shade8[ndimage.binary_erosion(holes, iterations=3)] = 228  # large gaps (water, no data) drawn flat
    Image.fromarray(shade8, "L").save(out / "hillshade.png", optimize=True)

    # per 6 ft square, per pass: count and typical time offset from the pass's start in this tile
    CW = CH = int(TILE / CELL)
    cc = np.clip(((x - x0) / CELL).astype(int), 0, CW - 1)
    cr = np.clip(((y1 - y) / CELL).astype(int), 0, CH - 1)
    key = ((cr * CW + cc) * NP + pidx)[good]
    counts = np.bincount(key, minlength=CH * CW * NP)
    offs = (t - pstart[pidx])[good]
    sums = np.bincount(key, weights=offs, minlength=CH * CW * NP)
    mean = np.divide(sums, counts, out=np.zeros_like(sums), where=counts > 0)
    cnt = np.minimum(counts, 255).astype(np.uint8)
    off = np.where(counts > 0, np.clip(np.round(mean * 10), 0, 65534), 65535).astype("<u2")
    with gzip.open(out / "cells.bin.gz", "wb", compresslevel=6) as f:
        f.write(cnt.tobytes()); f.write(off.tobytes())

    # random sample of individual points
    gi = np.flatnonzero(good)
    rng = np.random.default_rng(zlib.crc32(name.encode()))
    pick = np.sort(rng.choice(gi, size=min(N_SAMPLE, gi.size), replace=False))
    rec = np.zeros(pick.size, dtype=SAMPLE_DTYPE)
    rec["x"] = np.round((x[pick] - x0) / 0.04)
    rec["y"] = np.round((y[pick] - y0) / 0.04)
    rec["z"] = np.clip(np.round(z[pick] * 10), -32768, 32767)
    rec["ms"] = np.round((t[pick] - pstart[pidx[pick]]) * 1000)
    rec["p"] = pidx[pick]
    rec["c"] = cls[pick]
    rec["r"] = (ret[pick] << 4) | nret[pick]
    with gzip.open(out / "points.bin.gz", "wb", compresslevel=6) as f:
        f.write(rec.tobytes())

    # overview summaries
    hs_ov = shade8.reshape(int(TILE / OV_PIX), 5, int(TILE / OV_PIX), 5).mean(axis=(1, 3)).astype(np.uint8)
    OW = int(TILE / OV_CELL)
    oc = np.clip(((x - x0) / OV_CELL).astype(int), 0, OW - 1)
    orr = np.clip(((y1 - y) / OV_CELL).astype(int), 0, OW - 1)
    okey = ((orr * OW + oc) * NP + pidx)[good]
    ocount = np.bincount(okey, minlength=OW * OW * NP).reshape(OW * OW, NP)
    pass_day = np.array([p["day"] for p in passes])
    days = sorted(set(pass_day.tolist()))
    dcount = np.stack([ocount[:, pass_day == d].sum(axis=1) for d in days], axis=1)
    has = ocount.sum(axis=1) > 0
    ov_pass = np.where(has, np.array(pass_ids)[ocount.argmax(axis=1)], 0).astype(np.int32).reshape(OW, OW)
    ov_day = np.where(has, np.array(days)[dcount.argmax(axis=1)], 0).astype(np.int32).reshape(OW, OW)
    ov_ndays = (dcount >= 5).sum(axis=1).astype(np.uint8).reshape(OW, OW)
    np.savez_compressed(work / f"{name}.npz", hs=hs_ov, ov_pass=ov_pass, ov_day=ov_day, ov_ndays=ov_ndays)

    meta = {"name": name, "x0": x0, "y0": y0, "passes": passes, "n_sample": int(pick.size),
            "points": int(x.size), "noise": int(((cls == 7) | (cls == 18)).sum())}
    (work / f"{name}.json").write_text(json.dumps(meta))
    return meta


def cmd_process(incoming, work, site, watch=False):
    incoming, work, site = Path(incoming), Path(work), Path(site)
    work.mkdir(parents=True, exist_ok=True); site.mkdir(parents=True, exist_ok=True)
    log = open(work / "process.log", "a", buffering=1)
    pending = {}
    with ProcessPoolExecutor(max_workers=WORKERS) as pool:
        while True:
            for zp in sorted(incoming.glob("*.zip")):
                if zp not in pending:
                    pending[zp] = pool.submit(process_zip, zp, work, site)
            for zp, fut in list(pending.items()):
                if fut.done():
                    try:
                        res = fut.result()
                        res = res if isinstance(res, list) else [res]
                        for r in res:
                            msg = r.get("skipped") or f"{len(r['passes'])} passes, {r['points']:,} pts"
                            log.write(f"{time.strftime('%H:%M:%S')} {r.get('name', zp.stem)}: {msg}\n")
                        zp.unlink()
                    except Exception as e:  # keep the zip for a retry
                        log.write(f"{time.strftime('%H:%M:%S')} {zp.stem}: ERROR {e!r}\n")
                        zp.rename(zp.with_suffix(".failed"))
                    del pending[zp]
            if not pending and (not watch or (incoming / "DONE").exists()) and not list(incoming.glob("*.zip")):
                break
            time.sleep(1)
    log.write("finished\n")


# ---------------------------------------------------------------- step 3
def cmd_assemble(work, site):
    from PIL import Image
    from pyproj import Transformer
    work, site = Path(work), Path(site)
    metas = [json.loads(p.read_text()) for p in sorted(work.glob("w*.json"))]
    if not metas:
        sys.exit("no processed tiles")
    xs = [m["x0"] for m in metas]; ys = [m["y0"] for m in metas]
    X0, Y0, X1, Y1 = min(xs), min(ys), max(xs) + TILE, max(ys) + TILE

    # global passes and days
    gp = {}
    for m in metas:
        for p in m["passes"]:
            g = gp.setdefault(p["id"], {"id": p["id"], "gps_start": p["gps_start"], "points": 0, "tiles": 0})
            g["gps_start"] = min(g["gps_start"], p["gps_start"]); g["points"] += p["points"]; g["tiles"] += 1
    passes = sorted(gp.values(), key=lambda g: g["gps_start"])
    for g in passes:
        g["day"] = pacific_day(g["gps_start"])
        g["utc_start"] = gps_to_utc(g["gps_start"]).isoformat()
    pid_index = {g["id"]: i for i, g in enumerate(passes)}
    days = sorted({g["day"] for g in passes})
    day_index = {d: i for i, d in enumerate(days)}
    for g in passes:
        g["day"] = day_index[g["day"]]

    # overview rasters
    nx, ny = int((X1 - X0) / TILE), int((Y1 - Y0) / TILE)
    hs_n, ov_n = int(TILE / OV_PIX), int(TILE / OV_CELL)
    hs = np.full((ny * hs_n, nx * hs_n), 228, np.uint8)
    ovd = np.full((ny * ov_n, nx * ov_n), 255, np.uint8)
    ovp = np.full((ny * ov_n, nx * ov_n), 65535, "<u2")
    ovk = np.zeros((ny * ov_n, nx * ov_n), np.uint8)
    for m in metas:
        a = np.load(work / f"{m['name']}.npz")
        cx, ry = int((m["x0"] - X0) / TILE), int((Y1 - TILE - m["y0"]) / TILE)
        hs[ry * hs_n:(ry + 1) * hs_n, cx * hs_n:(cx + 1) * hs_n] = a["hs"]
        d = a["ov_day"]; p = a["ov_pass"]
        ovd[ry * ov_n:(ry + 1) * ov_n, cx * ov_n:(cx + 1) * ov_n] = np.where(d > 0, np.vectorize(lambda v: day_index.get(v, 255))(d), 255)
        ovp[ry * ov_n:(ry + 1) * ov_n, cx * ov_n:(cx + 1) * ov_n] = np.where(p > 0, np.vectorize(lambda v: pid_index.get(v, 65535))(p), 65535)
        ovk[ry * ov_n:(ry + 1) * ov_n, cx * ov_n:(cx + 1) * ov_n] = a["ov_ndays"]
    Image.fromarray(hs, "L").save(site / "overview_hillshade.png", optimize=True)
    with gzip.open(site / "overview.bin.gz", "wb") as f:
        f.write(ovd.tobytes()); f.write(ovp.tobytes()); f.write(ovk.tobytes())

    # lat/lon conversion: quadratic fit over the area (feet from SW corner -> degrees)
    fwd = Transformer.from_crs(SRC_CRS, "EPSG:4326", always_xy=True)
    inv = Transformer.from_crs("EPSG:4326", SRC_CRS, always_xy=True)
    gx, gy = np.meshgrid(np.linspace(X0, X1, 25), np.linspace(Y0, Y1, 25))
    lon, lat = fwd.transform(gx.ravel(), gy.ravel())
    dx, dy = (gx.ravel() - X0) / 10000, (gy.ravel() - Y0) / 10000
    A = np.column_stack([np.ones_like(dx), dx, dy, dx * dx, dx * dy, dy * dy])
    clon, *_ = np.linalg.lstsq(A, lon, rcond=None); clat, *_ = np.linalg.lstsq(A, lat, rcond=None)
    err = np.hypot((A @ clon - lon) * 111320 * np.cos(np.radians(lat.mean())), (A @ clat - lat) * 110540).max()
    px, py = inv.transform(PIN_LATLON[1], PIN_LATLON[0])

    tiles = [{"name": m["name"], "x0": m["x0"], "y0": m["y0"], "n": m["n_sample"],
              "passes": [[pid_index[p["id"]], p["gps_start"]] for p in m["passes"]]} for m in metas]
    index = {
        "project": "Cascades North Wali 2023", "crs": SRC_CRS, "bounds": [X0, Y0, X1, Y1], "tile": TILE,
        "pix": PIX, "cell": CELL, "ov_pix": OV_PIX, "ov_cell": OV_CELL,
        "sample": {"rec": SAMPLE_DTYPE.itemsize, "xy_unit": 0.04, "z_unit": 0.1},
        "days": [{"key": d, "label": datetime.strptime(str(d), "%Y%m%d").strftime("%a, %b %-d, %Y") if sys.platform != "win32" else datetime.strptime(str(d), "%Y%m%d").strftime("%a, %b %#d, %Y")} for d in days],
        "passes": passes, "tiles": tiles,
        "total_points": sum(m["points"] for m in metas), "noise_removed": sum(m["noise"] for m in metas),
        "lonlat_quad": {"lon": clon.tolist(), "lat": clat.tolist(), "scale_ft": 10000, "max_err_m": float(err)},
        "pin": {"lat": PIN_LATLON[0], "lon": PIN_LATLON[1], "x": px, "y": py},
    }
    (site / "index.json").write_text(json.dumps(index, separators=(",", ":")))
    print(f"{len(metas)} tiles, {len(passes)} passes, {len(days)} days, lat/lon fit error {err:.3f} m")
    print("days:", [d["label"] for d in index["days"]])


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "tiles":
        cmd_tiles(sys.argv[2])
    elif cmd == "process":
        cmd_process(sys.argv[2], sys.argv[3], sys.argv[4], watch="--watch" in sys.argv)
    elif cmd == "assemble":
        cmd_assemble(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
