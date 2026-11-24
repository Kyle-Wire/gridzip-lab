"""Protocol 6, step 1a: storm terms for every thunderstorm pair, from the station's own ISD hourly record.

python physics/storm_terms.py [--windsim ../WindSim] [--workers 8]
Reads WindSim's F2 pair files; streams each station-year from the ISD bucket (AWS Open Data, public domain), keeps
only the rows near each gust and deletes nothing it did not write. Writes results/physics/storm_terms.json.
"""
import argparse, json, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
URL = "https://noaa-global-hourly-pds.s3.amazonaws.com/{year}/{sid}.csv"
WIN = pd.Timedelta(hours=3)
OK_Q = set("014579")


def pairs(windsim: Path):
    out = []
    for grp in ("cal", "blind"):
        for f in sorted((windsim / "docs/validation/layer1/f2" / grp).glob("*.json")):
            d = json.loads(f.read_text()); ids = {s["icao"]: s["id"] for s in d.get("stations", [])}
            for p in d.get("pairs", []):
                if p["type"] == "thunderstorm":
                    out.append(dict(place=d["region"], group=grp, station=p["station"], sid=ids[p["station"]],
                                    storm=p["storm_time"], obs=p["observed_ms"], pred=p["predicted_ms"]))
    return out


def _val(s, scale, bad, qi=1):
    """ISD 'value,quality' fields -> float (scaled) or nan."""
    v = np.full(len(s), np.nan)
    for i, x in enumerate(s.fillna("").to_numpy()):
        parts = str(x).split(",")
        if len(parts) > qi and parts[0].strip() and parts[qi].strip() in OK_Q:
            try:
                n = int(parts[0])
            except ValueError:
                continue
            if abs(n) < bad: v[i] = n / scale
    return v


def met(raw: str) -> pd.DataFrame:
    cols = ["DATE", "REPORT_TYPE", "TMP", "DEW", "SLP", "MA1"]
    import io
    df = pd.read_csv(io.StringIO(raw), dtype=str, usecols=lambda c: c in cols)
    for c in cols:
        if c not in df: df[c] = np.nan
    df = df[df["REPORT_TYPE"].fillna("").str.strip().isin(["FM-15", "FM-16"])]
    t = pd.to_datetime(df["DATE"], utc=True)
    ma1_stn = _val(df["MA1"].fillna("").str.split(",").str[2:4].str.join(","), 10, 99999)   # station pressure
    out = pd.DataFrame({"time": t.to_numpy(), "T": _val(df["TMP"], 10, 9999), "Td": _val(df["DEW"], 10, 9999),
                        "p_stn": ma1_stn, "slp": _val(df["SLP"], 10, 99999)})
    out["time"] = pd.to_datetime(out["time"], utc=True)
    return out.sort_values("time").reset_index(drop=True)


def terms_for(p, gust: pd.DataFrame, m: pd.DataFrame) -> dict:
    t0 = pd.Timestamp(p["storm"]); w = gust[(gust.time >= t0 - WIN) & (gust.time <= t0 + WIN) & gust.gust_ms.notna()]
    r = {"matched": False}
    if w.empty:
        return r | {"no_gust_rows": True}
    c = w[(w.gust_ms - p["obs"]).abs() <= 0.3]
    if len(c):
        row = c.iloc[(c.time - t0).abs().argmin()]; r["matched"] = True
    else:
        row = w.loc[w.gust_ms.idxmax()]
    tg = pd.Timestamp(row.gust_time) if pd.notna(row.gust_time) else pd.Timestamp(row.time)
    r["gust_time"] = tg.isoformat(); r["gust_ms"] = float(row.gust_ms)
    before = m[(m.time >= tg - pd.Timedelta(minutes=120)) & (m.time <= tg)]
    after = m[(m.time > tg) & (m.time <= tg + pd.Timedelta(minutes=90))]
    ref = m[m.time <= tg - pd.Timedelta(minutes=30)].tail(1)
    T0 = before["T"].max(); Tmin = after["T"].min()
    r["T0"] = None if np.isnan(T0) else float(T0)
    r["dT"] = None if (np.isnan(T0) or np.isnan(Tmin)) else float(T0 - Tmin)
    pcol = "p_stn" if m["p_stn"].notna().any() else "slp"
    p_ref = ref[pcol].iloc[0] if len(ref) else np.nan; p_max = after[pcol].max()
    r["dp"] = None if (np.isnan(p_ref) or np.isnan(p_max)) else float(p_max - p_ref); r["p_source"] = pcol
    td_ref = ref["Td"].iloc[0] if len(ref) else np.nan; td_min = after["Td"].min()
    r["dTd"] = None if (np.isnan(td_ref) or np.isnan(td_min)) else float(td_ref - td_min)
    t_ref = ref["T"].iloc[0] if len(ref) else np.nan
    r["dryness"] = None if (np.isnan(t_ref) or np.isnan(td_ref)) else float(t_ref - td_ref)
    return r


def fetch(sid, year, tries=4):
    for k in range(tries):
        try:
            with urllib.request.urlopen(URL.format(year=year, sid=sid), timeout=120) as f:
                return f.read().decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001
            if k == tries - 1: raise
            time.sleep(2 ** k)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--windsim", default=str(ROOT.parent / "WindSim"))
    ap.add_argument("--workers", type=int, default=8); a = ap.parse_args()
    sys.path.insert(0, a.windsim)
    from windsim.layer1.isd import parse_year
    P = pairs(Path(a.windsim)); groups = {}
    for i, p in enumerate(P): groups.setdefault((p["sid"], p["storm"][:4]), []).append(i)
    out = ROOT / "results/physics/storm_terms.json"
    done = json.loads(out.read_text()) if out.exists() else {}
    res = {int(k): v for k, v in done.get("pairs", {}).items()}

    def job(key):
        sid, yr = key
        raw = fetch(sid, yr); g = parse_year(raw); m = met(raw)
        return key, {i: terms_for(P[i], g, m) for i in groups[key]}
    todo = [k for k in groups if not all(i in res for i in groups[k])]
    print(f"{len(P)} pairs, {len(groups)} station-years, {len(todo)} to fetch", flush=True)
    with ThreadPoolExecutor(a.workers) as ex:
        futs = [ex.submit(job, k) for k in todo]
        for n, f in enumerate(as_completed(futs), 1):
            try:
                key, r = f.result(); res.update(r)
            except Exception as e:  # noqa: BLE001
                print("fail", e, flush=True)
            if n % 25 == 0 or n == len(futs):
                out.write_text(json.dumps({"pairs": {str(k): v for k, v in res.items()}, "meta": P}, default=float))
                print(f"{n}/{len(futs)}", flush=True)
    out.write_text(json.dumps({"pairs": {str(k): v for k, v in res.items()}, "meta": P}, default=float))


if __name__ == "__main__":
    main()
