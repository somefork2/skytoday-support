#!/usr/bin/env python3
"""Build the Sky Today data feed (tools/feed/out/feed.json).

Sources (all public, no keys):
  * JPL SBDB Query API  — every comet with q < 5 AU, M1/K1 and orbit, screened for brightness
  * JPL Horizons API    — fresh osculating elements (heliocentric, J2000 ecliptic) for what we publish
  * JPL SBDB API        — designation, SPK-ID, H/G (and fallback elements) for close-approach asteroids
  * JPL CNEOS CAD API   — asteroid close approaches within 0.05 AU over the next 60 days
  * CelesTrak GP API    — TLEs for ISS, Tiangong, Hubble (one request each)
  * NOAA SWPC           — planetary Kp: last 24 h observed/estimated + 3-day forecast
  * tools/feed/notes.json (optional) — curated notes, merged in

Usage:  .venv/bin/python tools/feed/build_feed.py [--now 2026-10-01T00:00:00Z] [--out PATH]
Exits non-zero (and writes nothing) when a source fails or returns nothing.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import re
import sys
import time
from pathlib import Path

import requests

HERE = Path(__file__).resolve().parent
UA = "SkyTodayFeed/1.0 (+https://github.com/; daily astronomy feed)"
SBDB_QUERY = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"
CAD = "https://ssd-api.jpl.nasa.gov/cad.api"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
CELESTRAK = "https://celestrak.org/NORAD/elements/gp.php"
KP_FORECAST = "https://services.swpc.noaa.gov/products/noaa-planetary-k-index-forecast.json"

COMET_MAG_LIMIT = 12.0      # publish comets predicted at least this bright…
COMET_WINDOW_DAYS = 180     # …at some point in the next ~6 months
SCREEN_SLACK = 1.5          # screening with catalog elements keeps comets up to limit + slack
MAX_COMETS = 15
CA_DAYS = 60
CA_DIST_AU = 0.05
CA_BODIES = 5               # closest approaches whose elements go into "bodies"
SATELLITES = [("ISS", "ISS", 25544), ("CSS", "Tiangong", 48274), ("HST", "Hubble", 20580)]

AU_KM = 149_597_870.7
LD_KM = 384_400.0
K_GAUSS = 0.01720209895
TT_MINUS_UTC = 69.184       # s (37 leap seconds + 32.184); TDB−TT < 2 ms


class FeedError(RuntimeError):
    pass


def log(*a):
    print(*a, file=sys.stderr, flush=True)


# ---------------------------------------------------------------- HTTP

_session = requests.Session()
_session.headers["User-Agent"] = UA


def http_get(url: str, params: dict | None = None, *, kind: str = "json", attempts: int = 4, timeout: float = 60):
    last: Exception | None = None
    for n in range(attempts):
        if n:
            time.sleep(2 ** n)
        try:
            r = _session.get(url, params=params, timeout=timeout)
            if r.status_code == 429 or r.status_code >= 500:
                last = FeedError(f"HTTP {r.status_code} from {url}")
                continue
            if r.status_code != 200:
                raise FeedError(f"HTTP {r.status_code} from {url}: {r.text[:200]}")
            return r.json() if kind == "json" else r.text
        except (requests.ConnectionError, requests.Timeout, ValueError) as e:
            last = e
    raise FeedError(f"giving up on {url}: {last}")


# ---------------------------------------------------------------- time

def jd_of(t: dt.datetime) -> float:
    return t.timestamp() / 86400 + 2440587.5


def iso(t: dt.datetime) -> str:
    return t.astimezone(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def utc_of_jd_tdb(jd: float) -> dt.datetime:
    return dt.datetime.fromtimestamp((jd - 2440587.5) * 86400 - TT_MINUS_UTC, dt.timezone.utc)


def parse_iso(s: str) -> dt.datetime:
    t = dt.datetime.fromisoformat(s.replace("Z", "+00:00"))
    return t if t.tzinfo else t.replace(tzinfo=dt.timezone.utc)


# ---------------------------------------------------------------- two-body (mirrors SkyKit SmallBody)

def helio(el: dict, jd: float) -> tuple[float, float, float]:
    q, e = el["q"], el["e"]
    t = jd - el["tp"]
    if abs(e - 1) < 1e-4:
        w = 3 * K_GAUSS / math.sqrt(2 * q) / q * t
        yb = (w / 2 + math.sqrt(w * w / 4 + 1)) ** (1 / 3)
        s = yb - 1 / yb
        nu = 2 * math.atan(s)
        r = q * (1 + s * s)
        x, y = r * math.cos(nu), r * math.sin(nu)
    elif e < 1:
        a = q / (1 - e)
        M = math.remainder(K_GAUSS / a ** 1.5 * t, 2 * math.pi)
        E = M if e < 0.8 else math.copysign(math.pi, M)
        for _ in range(60):
            d = (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
            E -= d
            if abs(d) < 1e-13:
                break
        x, y = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E)
    else:
        a = q / (e - 1)
        M = K_GAUSS / a ** 1.5 * t
        H = math.asinh(M / e)
        for _ in range(60):
            d = (e * math.sinh(H) - H - M) / (e * math.cosh(H) - 1)
            H -= d
            if abs(d) < 1e-13:
                break
        x, y = a * (e - math.cosh(H)), a * math.sqrt(e * e - 1) * math.sinh(H)
    r = math.radians
    cw, sw, cn, sn, ci, si = (math.cos(r(el["peri"])), math.sin(r(el["peri"])), math.cos(r(el["node"])),
                              math.sin(r(el["node"])), math.cos(r(el["i"])), math.sin(r(el["i"])))
    return ((cw * cn - sw * sn * ci) * x + (-sw * cn - cw * sn * ci) * y,
            (cw * sn + sw * cn * ci) * x + (-sw * sn + cw * cn * ci) * y,
            (sw * si) * x + (cw * si) * y)


# Earth–Moon barycentre, Standish (1992) mean elements, J2000 ecliptic — ~0.0001 AU, plenty for magnitudes
_EMB = dict(a=(1.00000261, 0.00000562), e=(0.01671123, -0.00004392), I=(-0.00001531, -0.01294668),
            L=(100.46457166, 35999.37244981), w=(102.93768193, 0.32327364), O=(0.0, 0.0))


def earth(jd: float) -> tuple[float, float, float]:
    T = (jd - 2451545.0) / 36525
    p = {k: v[0] + v[1] * T for k, v in _EMB.items()}
    e = p["e"]
    a = p["a"]
    M = math.radians((p["L"] - p["w"]) % 360)
    E = M
    for _ in range(20):
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    x, y = a * (math.cos(E) - e), a * math.sqrt(1 - e * e) * math.sin(E)
    peri = math.radians(p["w"] - p["O"])
    node, inc = math.radians(p["O"]), math.radians(p["I"])
    cw, sw, cn, sn, ci, si = math.cos(peri), math.sin(peri), math.cos(node), math.sin(node), math.cos(inc), math.sin(inc)
    return ((cw * cn - sw * sn * ci) * x + (-sw * cn - cw * sn * ci) * y,
            (cw * sn + sw * cn * ci) * x + (-sw * sn + cw * cn * ci) * y,
            (sw * si) * x + (cw * si) * y)


def comet_mag(el: dict, jd: float) -> tuple[float, float, float]:
    p = helio(el, jd)
    g = earth(jd)
    r = math.dist(p, (0, 0, 0))
    delta = math.dist(p, g)
    k1 = el["k1"] if el.get("k1") is not None else 10.0
    return el["m1"] + 5 * math.log10(delta) + k1 * math.log10(r), r, delta


def peak(el: dict, jd0: float, days: int = COMET_WINDOW_DAYS, step: float = 2.0) -> tuple[float, float]:
    """Brightest predicted total magnitude in [jd0, jd0+days] and its JD."""
    best = (99.0, jd0)
    t = jd0
    while t <= jd0 + days:
        m = comet_mag(el, t)[0]
        if m < best[0]:
            best = (m, t)
        t += step
    return best


# ---------------------------------------------------------------- Horizons elements

def horizons_elements(spkid: str | int, epoch_jd: float, comet: bool) -> dict:
    cmd = f"DES={spkid};" + ("CAP;" if comet else "")
    data = http_get(HORIZONS, {
        "format": "json", "COMMAND": f"'{cmd}'", "EPHEM_TYPE": "ELEMENTS", "CENTER": "'500@10'",
        "REF_PLANE": "ECLIPTIC", "REF_SYSTEM": "J2000", "OUT_UNITS": "AU-D", "TLIST": f"'{epoch_jd:.6f}'",
        "TLIST_TYPE": "JD", "TIME_TYPE": "TDB", "CSV_FORMAT": "YES", "OBJ_DATA": "NO", "MAKE_EPHEM": "YES",
    })
    if "error" in data:
        raise FeedError(f"Horizons {spkid}: {data['error']}")
    res = data["result"]
    if "$$SOE" not in res:
        raise FeedError(f"Horizons {spkid}: no ephemeris ({res[-300:].strip()})")
    head = res[:res.index("$$SOE")].strip().splitlines()
    cols = [c.strip() for c in head[-2].split(",")]
    row = [c.strip() for c in res[res.index("$$SOE") + 5:res.index("$$EOE")].strip().split(",")]
    v = dict(zip(cols, row))
    return dict(epoch=float(v["JDTDB"]), e=float(v["EC"]), q=float(v["QR"]), i=float(v["IN"]),
                node=float(v["OM"]), peri=float(v["W"]), tp=float(v["Tp"]))


def with_fresh_elements(body: dict, spkid, epoch_jd: float) -> dict:
    try:
        body.update(horizons_elements(spkid, epoch_jd, body["kind"] == "comet"))
        body["src"] = "horizons"
    except FeedError as e:
        log(f"  ! {body['id']}: {e} — keeping SBDB elements (epoch {body['epoch']})")
        body["src"] = "sbdb"
    return body


# ---------------------------------------------------------------- comets

_NUMBERED = re.compile(r"^\d+[PDI](-[A-Z]+)?$")


def comet_id(pdes: str, prefix: str | None) -> str:
    return pdes if _NUMBERED.match(pdes) or not prefix else f"{prefix}/{pdes}"


def fetch_comets(now: dt.datetime) -> list[dict]:
    log("comets: SBDB query…")
    data = http_get(SBDB_QUERY, {
        "fields": "spkid,pdes,prefix,name,epoch,e,q,i,om,w,tp,M1,K1,last_obs,data_arc",
        "sb-kind": "c", "full-prec": "true", "sb-cdata": json.dumps({"AND": ["q|LT|5"]}),
    })
    rows = data.get("data") or []
    if not rows:
        raise FeedError("SBDB query returned no comets")
    jd0 = jd_of(now)
    cands = []
    for spkid, pdes, prefix, name, epoch, e, q, i, om, w, tp, m1, k1, last_obs, arc in rows:
        if m1 is None or tp is None or prefix in ("D", "X"):
            continue
        # Drop orbits nobody can trust to return: a few days of arc long ago (e.g. SOHO sungrazer fragments).
        # Numbered comets have several apparitions; anything seen in the last two months is a live comet.
        y, mo, d = (re.findall(r"\d+", (last_obs or "").replace("??", "01")) + ["1", "1", "1"])[:3]
        seen_days = (now - dt.datetime(int(y), int(mo), int(d), tzinfo=dt.timezone.utc)).days if last_obs else 1e9
        if not _NUMBERED.match(pdes) and seen_days > 60 and int(arc or 0) < 30:
            continue
        el = dict(id=comet_id(pdes, prefix), name=(name or comet_id(pdes, prefix)).replace("-", "–"),
                  kind="comet", q=float(q), e=float(e), i=float(i), node=float(om), peri=float(w), tp=float(tp),
                  m1=float(m1), k1=float(k1) if k1 is not None else None, epoch=float(epoch))
        m, _ = peak(el, jd0, step=4)
        if m <= COMET_MAG_LIMIT + SCREEN_SLACK:
            cands.append((m, int(spkid), el))
    log(f"comets: {len(rows)} in SBDB, {len(cands)} pass screening; fresh elements from Horizons…")
    epoch = math.floor(jd0 - 0.5) + 0.5
    out = []
    for _, spkid, el in sorted(cands, key=lambda c: (c[0], c[2]["id"])):
        el = with_fresh_elements(el, spkid, epoch)
        m, when = peak(el, jd0)
        if m <= COMET_MAG_LIMIT:
            el["peakMag"], el["peakDate"] = round(m, 1), iso(utc_of_jd_tdb(when))[:10]
            el["mag"] = round(comet_mag(el, jd0)[0], 1)
            out.append(el)
        log(f"  {el['id']:<14} {el['name']:<24} peak {m:5.1f} {'✓' if m <= COMET_MAG_LIMIT else ''}")
    out.sort(key=lambda b: (b["peakMag"], b["id"]))
    return out[:MAX_COMETS]


# ---------------------------------------------------------------- asteroids

def diameter_m(h: float) -> list[int]:
    d = lambda p: 1329e3 / math.sqrt(p) * 10 ** (-h / 5)
    return [round(d(0.25)), round(d(0.05))]


def asteroid_name(fullname: str, des: str) -> str:
    f = fullname.strip()
    m = re.match(r"^(\d+)\s+(.+?)\s+\((.+)\)$", f)          # 99942 Apophis (2004 MN4)
    if m:
        return m.group(2)
    m = re.match(r"^(\d+)\s+\((.+)\)$", f)                  # 523934 (1998 FF14)
    if m:
        return f"{m.group(1)} ({m.group(2)})"
    m = re.match(r"^(\d+)\s+(\S.*)$", f)                    # 433 Eros
    if m:
        return m.group(2)
    return f.strip("()") or des


def fetch_approaches(now: dt.datetime) -> tuple[list[dict], list[dict]]:
    log("approaches: CNEOS CAD…")
    data = http_get(CAD, {"date-min": now.strftime("%Y-%m-%dT%H:%M:%S"), "date-max": f"+{CA_DAYS}",
                          "dist-max": str(CA_DIST_AU), "sort": "date", "fullname": "true", "diameter": "true"})
    f = data.get("fields") or []
    rows = [dict(zip(f, r)) for r in data.get("data") or []]
    if not rows:
        raise FeedError("CAD API returned no close approaches")
    approaches = []
    for r in rows:
        jd = float(r["jd"])
        km = float(r["dist"]) * AU_KM
        h = float(r["h"]) if r.get("h") else None
        if r.get("diameter"):
            dm = round(float(r["diameter"]) * 1000)
            dia = [dm, dm]
        else:
            dia = diameter_m(h) if h is not None else None
        approaches.append({
            "id": r["des"], "name": asteroid_name(r.get("fullname") or "", r["des"]),
            "date": iso(utc_of_jd_tdb(jd)), "distanceKm": round(km), "lunarDistances": round(km / LD_KM, 2),
            "diameterM": dia, "speedKms": round(float(r["v_rel"]), 2), "h": h, "_jd": jd,
        })
    approaches.sort(key=lambda a: (a["date"], a["id"]))

    log(f"approaches: {len(approaches)}; elements for the {CA_BODIES} closest…")
    bodies = []
    for a in sorted(approaches, key=lambda a: (a["distanceKm"], a["id"]))[:CA_BODIES]:
        d = http_get(SBDB, {"sstr": a["id"], "full-prec": "true", "phys-par": "true"})
        if "object" not in d:
            raise FeedError(f"SBDB has no object {a['id']}: {d}")
        el = {x["name"]: x["value"] for x in d["orbit"]["elements"]}
        phys = {p["name"]: p["value"] for p in d.get("phys_par") or []}
        h = float(phys["H"]) if phys.get("H") else a["h"]
        g = float(phys["G"]) if phys.get("G") else 0.15
        b = dict(id=a["id"], name=a["name"], kind="asteroid", q=float(el["q"]), e=float(el["e"]), i=float(el["i"]),
                 node=float(el["om"]), peri=float(el["w"]), tp=float(el["tp"]), h=h, g=g,
                 epoch=float(d["orbit"]["epoch"]))
        # elements osculating at the moment of the flyby: exact where pointing matters most (a close pass bends
        # the heliocentric orbit, so these drift away from the flyby — fine with a feed refreshed daily)
        bodies.append(with_fresh_elements(b, d["object"]["spkid"], round(a["_jd"], 6)))
        log(f"  {a['id']:<12} {a['date']} {a['distanceKm']:>10,} km  H={h}")
    bodies.sort(key=lambda b: (next(x["date"] for x in approaches if x["id"] == b["id"]), b["id"]))
    for a in approaches:
        del a["_jd"]
    return approaches, bodies


# ---------------------------------------------------------------- satellites

def tle_checksum(line: str) -> int:
    return sum(int(c) if c.isdigit() else c == "-" for c in line[:68]) % 10


def fetch_satellites(now: dt.datetime) -> list[dict]:
    out = []
    for sid, name, norad in SATELLITES:
        txt = http_get(CELESTRAK, {"CATNR": norad, "FORMAT": "TLE"}, kind="text")
        lines = [l.rstrip() for l in txt.splitlines() if l.strip()]
        l1 = next((l for l in lines if l.startswith("1 ")), None)
        l2 = next((l for l in lines if l.startswith("2 ")), None)
        ok = (l1 and l2 and len(l1) == 69 and len(l2) == 69 and int(l1[2:7]) == norad and int(l2[2:7]) == norad
              and tle_checksum(l1) == int(l1[68]) and tle_checksum(l2) == int(l2[68]))
        if not ok:
            raise FeedError(f"CelesTrak: bad TLE for {norad}: {txt[:200]!r}")
        out.append({"id": sid, "name": name, "norad": norad, "tle1": l1, "tle2": l2, "fetched": iso(now)})
        log(f"satellite {sid}: epoch {l1[18:32]}")
    return out


# ---------------------------------------------------------------- aurora

def fetch_aurora(now: dt.datetime) -> dict:
    data = http_get(KP_FORECAST)
    if data and isinstance(data[0], list):                    # older array-of-arrays format
        head, data = data[0], [dict(zip(data[0], r)) for r in data[1:]]
    rows = []
    for r in data or []:
        t = parse_iso(r["time_tag"])
        status = (r.get("observed") or "").lower()
        # "estimated" covers the whole current UT day, including 3-h slots still ahead of us
        observed = status == "observed" or (status == "estimated" and t <= now)
        if r.get("kp") is None or (observed and t < now - dt.timedelta(hours=24)):
            continue
        rows.append({"time": iso(t), "kp": round(float(r["kp"]), 2), "observed": observed})
    rows.sort(key=lambda r: r["time"])
    past = [r for r in rows if r["observed"]]
    if not rows or not past or not any(not r["observed"] for r in rows):
        raise FeedError("SWPC Kp forecast returned no usable rows")
    cur = past[-1]
    log(f"aurora: {len(rows)} Kp rows, current {cur['kp']} at {cur['time']}")
    return {"current": {"time": cur["time"], "kp": cur["kp"]}, "kp": rows, "fetched": iso(now)}


# ---------------------------------------------------------------- notes

def load_notes(path: Path, now: dt.datetime) -> list[dict]:
    if not path.exists():
        return []
    notes = json.loads(path.read_text())
    if not isinstance(notes, list):
        raise FeedError(f"{path}: expected a JSON array")
    out = []
    for n in notes:
        for k in ("id", "start", "end", "title", "body"):
            if k not in n:
                raise FeedError(f"{path}: note {n.get('id')!r} lacks {k!r}")
        if "en" not in n["title"] or "en" not in n["body"]:
            raise FeedError(f"{path}: note {n['id']!r} needs English title and body")
        start, end = parse_iso(n["start"]), parse_iso(n["end"])
        if end < now:
            continue                                       # expired
        out.append({"id": n["id"], "start": iso(start), "end": iso(end), "title": n["title"], "body": n["body"]})
    ids = [n["id"] for n in out]
    if len(ids) != len(set(ids)):
        raise FeedError(f"{path}: duplicate note ids")
    return sorted(out, key=lambda n: (n["start"], n["id"]))


# ---------------------------------------------------------------- output

def rounded(b: dict) -> dict:
    """Stable key order and enough digits for arc-second positions."""
    out = {k: b[k] for k in ("id", "name", "kind")}
    for k, nd in (("q", 10), ("e", 10), ("i", 7), ("node", 7), ("peri", 7), ("tp", 6), ("epoch", 6)):
        out[k] = round(b[k], nd)
    for k in ("m1", "k1", "h", "g", "mag", "peakMag", "peakDate"):
        if b.get(k) is not None:
            out[k] = b[k]
    return out


def build(now: dt.datetime, notes_path: Path) -> dict:
    comets = fetch_comets(now)
    approaches, asteroids = fetch_approaches(now)
    return {
        "version": 1,
        "generated": iso(now),
        "bodies": [rounded(b) for b in comets + asteroids],
        "approaches": approaches,
        "satellites": fetch_satellites(now),
        "aurora": fetch_aurora(now),
        "notes": load_notes(notes_path, now),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--now", help="ISO time to build for (default: current UTC)")
    ap.add_argument("--out", default=str(HERE / "out" / "feed.json"))
    ap.add_argument("--notes", default=str(HERE / "notes.json"))
    args = ap.parse_args()
    now = parse_iso(args.now) if args.now else dt.datetime.now(dt.timezone.utc)
    try:
        feed = build(now, Path(args.notes))
    except FeedError as e:
        log(f"FEED BUILD FAILED: {e}")
        return 1
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(feed, ensure_ascii=False, separators=(",", ":")) + "\n")
    tmp.replace(out)
    log(f"wrote {out} ({out.stat().st_size:,} bytes): {sum(b['kind'] == 'comet' for b in feed['bodies'])} comets, "
        f"{sum(b['kind'] == 'asteroid' for b in feed['bodies'])} asteroids, {len(feed['approaches'])} approaches, "
        f"{len(feed['satellites'])} satellites, {len(feed['aurora']['kp'])} Kp rows, {len(feed['notes'])} notes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
