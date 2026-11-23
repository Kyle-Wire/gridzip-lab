"""Surprise score on a real block: FloodSim's finished South Street runs (Climate Home Inspection, 16 Nov 2026).

Unit = a building footprint (Philadelphia L&I footprints). Reading = the 90th-percentile max water depth in the 2 m ring
of open ground around the footprint, for each design storm (10, 100, 500 yr; 1 h and 24 h), averaged over the run
variants FloodSim made for that storm; the variants' spread is the model noise. Coding scale: log10(depth + 3 cm), at
0.02 (about 5 % in depth). Neighbours: footprints within 30 m (weights fall off with sigma 15 m). Descriptor: the ring's median ground height (3DEP 1 m
DEM): "surprise given height" compares a house only with neighbours at about the same height, so the score can say how
much of a house's oddity its height explains.

Also two planted tests on the same block, because the real block has no known answer:
  plant   30 random houses get +X cm of water in every storm's ring (a local low spot or concentrated roof runoff);
          detection = share ranked above the block's 95th percentile, X = 2, 5, 10, 20 cm.
  shuffle readings are shuffled between houses within each storm (location carries no information): the score
          should then find nothing beyond chance (no structure in where the surprises sit).
Usage: python surprise/flood_block.py <CHI folder with runs/geotiff, data/raw> <out dir> [--private]
--private adds the two sample homes (their keys are read from private/ and never written to the output).
"""
import json, sys, os, glob, re
import numpy as np
import rasterio
from rasterio.features import rasterize
from rasterio.warp import transform_geom
from scipy import ndimage
from pyproj import Transformer
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from surprise import core

RES, LO, HI, FLOOR, RADIUS = 0.02, -2.0, 1.0, 0.03, 30.0
STORMS = ["10-yr-1-h", "100-yr-1-h", "500-yr-1-h", "10-yr-24-h", "100-yr-24-h", "500-yr-24-h"]


def load(chi):
    runs = sorted(glob.glob(f"{chi}/runs/geotiff/POC South Street/*/max_depth.tif"))
    by = {s: [] for s in STORMS}
    for p in runs:
        name = os.path.basename(os.path.dirname(p))
        for s in STORMS:
            if re.search(r"_" + s + r"-Exist", name): by[s].append(p)
    ref = rasterio.open(by[STORMS[0]][0])
    fp = json.load(open(f"{chi}/data/raw/south-street/li_footprints.p0.geojson"))["features"]
    geoms = [transform_geom("EPSG:4326", ref.crs, f["geometry"]) for f in fp]
    ids = rasterize(((g, i + 1) for i, g in enumerate(geoms)), out_shape=ref.shape, transform=ref.transform, fill=0, dtype="int32")
    bld = ids > 0
    dem = rasterio.open(f"{chi}/data/raw/south-street/dem_3dep_1m_block.tif")
    rows, cols = np.mgrid[0:ref.shape[0], 0:ref.shape[1]]
    xs, ys = rasterio.transform.xy(ref.transform, rows.ravel(), cols.ravel())
    t = Transformer.from_crs(ref.crs, dem.crs, always_xy=True); X, Y = t.transform(np.array(xs), np.array(ys))
    z = np.array([v[0] for v in dem.sample(zip(X, Y))], float).reshape(ref.shape); z[z < -100] = np.nan
    depth = {s: [np.nan_to_num(rasterio.open(p).read(1), nan=0.0) for p in by[s]] for s in STORMS}
    return ref, geoms, fp, ids, bld, z, depth, {s: len(by[s]) for s in STORMS}


def rings(ids, bld, n, width=2):
    out = []
    for i in range(1, n + 1):
        sl = ndimage.find_objects((ids == i).astype(np.int8))
        if not sl or sl[0] is None: out.append(None); continue
        r0, c0 = max(sl[0][0].start - width - 1, 0), max(sl[0][1].start - width - 1, 0)
        sub = ids[r0:sl[0][0].stop + width + 1, c0:sl[0][1].stop + width + 1] == i
        ring = ndimage.binary_dilation(sub, iterations=width) & ~bld[r0:r0 + sub.shape[0], c0:c0 + sub.shape[1]]
        rr, cc = np.nonzero(ring); out.append((rr + r0, cc + c0) if len(rr) >= 3 else None)
    return out


def readings(ring_cells, depth, z):
    U = len(ring_cells); Y = np.full((U, len(STORMS)), np.nan); N = np.zeros_like(Y); G = np.full(U, np.nan)
    for i, rc in enumerate(ring_cells):
        if rc is None: continue
        G[i] = np.nanmedian(z[rc])
        for s, k in zip(STORMS, range(len(STORMS))):
            v = np.array([np.log10(np.percentile(d[rc], 90) + FLOOR) for d in depth[s]])
            Y[i, k] = v.mean(); N[i, k] = v.std() if len(v) > 1 else 0.0
    return Y, N, G


def main(chi, out, private=False):
    os.makedirs(out, exist_ok=True)
    ref, geoms, fp, ids, bld, z, depth, nvar = load(chi)
    rc = rings(ids, bld, len(geoms))
    Y, N, G = readings(rc, depth, z)
    ok = np.isfinite(Y).all(1) & np.isfinite(G)
    cent = np.array([[np.mean(np.array(g["coordinates"][0])[:, 0]), np.mean(np.array(g["coordinates"][0])[:, 1])] if g["type"] == "Polygon" else [np.nan, np.nan] for g in geoms])
    ok &= np.isfinite(cent).all(1)
    idx = np.where(ok)[0]; y, nz, xy, g = Y[idx], N[idx], cent[idx], G[idx]
    nb = core.neighbours(xy, RADIUS); gm = np.array([w @ g[j] for j, w in nb]); dz = (g - gm)[:, None]
    r = core.score(y, xy, RADIUS, RES, LO, HI, noise=nz, descriptor=g, verify=True)
    res = {"block": "south-street", "units": int(len(idx)), "storms": STORMS, "variants": nvar, "radius_m": RADIUS,
           "coding": {"scale": "log10(depth_m + 0.03)", "res": RES, "lo": LO, "hi": HI}, "h": r["h"], "tau_m": r["tau"],
           "roundtrip_ok": r["roundtrip_ok"], "total_bits": r["total_bits"],
           "percentiles": {p: float(np.percentile(r["surprise"], p)) for p in (5, 25, 50, 75, 90, 95, 99)},
           "percentiles_left": {p: float(np.percentile(r["surprise_given"], p)) for p in (50, 90, 95, 99)}}
    p95 = np.percentile(r["surprise"], 95)
    # planted tests
    rng = np.random.default_rng(7); plants = {}
    d_lin = 10 ** y - FLOOR
    for X in (0.02, 0.05, 0.10, 0.20):
        hits = []; gains = []
        for rep in range(3):
            pick = rng.choice(len(idx), 30, replace=False); y2 = y.copy()
            y2[pick] = np.log10(np.maximum(d_lin[pick], 0) + X + FLOOR)
            r2 = core.score(y2, xy, RADIUS, RES, LO, HI, noise=nz, hs=[r["h"]], verify=False)
            hits.append(float(np.mean(r2["surprise"][pick] > np.percentile(r2["surprise"], 95))))
            gains.append(float(np.median(r2["surprise"][pick] - r["surprise"][pick])))
        plants[f"+{int(X*100)} cm"] = {"detected_share_mean": float(np.mean(hits)), "reps": hits, "median_extra_bits": float(np.mean(gains))}
    res["planted"] = plants
    sh = y.copy()
    for k in range(y.shape[1]): sh[:, k] = rng.permutation(sh[:, k])
    rs = core.score(sh, xy, RADIUS, RES, LO, HI, noise=nz, verify=False)
    # spatial structure of surprise: Moran's I-like neighbour correlation, real vs shuffled
    def nbcorr(s):
        s = (s - s.mean()) / s.std(); return float(np.mean([s[i] * (w @ s[j]) for i, (j, w) in enumerate(nb)]))
    res["shuffle"] = {"percentiles": {p: float(np.percentile(rs["surprise"], p)) for p in (50, 95, 99)},
                      "neighbour_corr_real": nbcorr(r["surprise"]), "neighbour_corr_shuffled": nbcorr(rs["surprise"])}
    homes = {}
    if private:
        for key in sorted(glob.glob(f"{chi}/private/property-key*.json")):
            k = json.load(open(key)); t = Transformer.from_crs("EPSG:4326", ref.crs, always_xy=True)
            x, yy = t.transform(k["lon"], k["lat"]); row, col = rasterio.transform.rowcol(ref.transform, x, yy)
            cand = ids[max(row - 3, 0):row + 4, max(col - 3, 0):col + 4]; cand = cand[cand > 0]
            if not len(cand): continue
            b = np.bincount(cand).argmax() - 1; m = np.where(idx == b)[0]
            if not len(m): continue
            m = int(m[0]); s = r["surprise"]
            homes[k["id"]] = {"surprise_bits": float(s[m]), "percentile": float((s < s[m]).mean() * 100), "sign": int(r["sign"][m]),
                              "left_after_height_bits": float(r["surprise_given"][m]), "height_vs_neighbours_m": float(dz[m, 0]),
                              "ring_p90_depth_cm": [float(100 * (10 ** v - FLOOR)) for v in y[m]],
                              "neighbour_p90_depth_cm": [float(100 * (10 ** v - FLOOR)) for v in r["mu"][m]]}
        res["homes"] = homes
    json.dump(res, open(f"{out}/flood_block_south_street.json", "w"), indent=1)
    np.savez_compressed(f"{out}/flood_block_south_street.npz", idx=idx, surprise=r["surprise"], left=r["surprise_given"], sign=r["sign"], y=y, mu=r["mu"], dz=dz[:, 0], xy=xy, shuffled=rs["surprise"])
    return res, r, (ref, geoms, ids, idx, homes)


if __name__ == "__main__":
    res, *_ = main(sys.argv[1], sys.argv[2], "--private" in sys.argv)
    print(json.dumps({k: v for k, v in res.items()}, indent=1)[:4000])
