"""
Data-driven temperature × feedstock property model.

Maps pyrolysis temperature (and optionally feedstock) to biochar composition
targets (H/C, O/C, aromaticity) and a set of reference properties, fit from
public characterization data.

Two stages:

* :func:`build_model` (offline / dev) parses the source CSVs, computes molar
  H/C and O/C from elemental mass %, robustly smooths each property versus
  temperature, fits ``aromaticity = f(H/C)``, and writes a compact JSON
  artifact ``biochar/data/temperature_model.json``.
* :class:`TemperatureModel` (runtime) loads that JSON and answers queries with
  ``numpy.interp`` only — no scikit-learn/scipy/pandas needed.

Primary data source (definitive): **UC Davis Biochar Database**
(https://biochar.ucdavis.edu/), ``davis-biochar-db.csv``.
Aromaticity calibration only: ``biochar_data.csv`` (NMR aromaticity vs H/C).
Methodological parent: Wood, Mašek & Erastova, *Cell Reports Physical Science*
5(7), 2024, DOI 10.1016/j.xcrp.2024.102036.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).parent / "data"
_MODEL_PATH = _DATA_DIR / "temperature_model.json"
_DAVIS_CSV = _DATA_DIR / "davis-biochar-db.csv"
_AROM_CSV = _DATA_DIR / "biochar_data.csv"

# Molar masses (g/mol) for converting elemental mass % to molar ratios.
_M_H, _M_C, _M_O = 1.008, 12.011, 15.999

# Smoothing grid and kernel.
_GRID = list(range(100, 1001, 25))     # °C
_BANDWIDTH = 50.0                       # °C, Gaussian kernel
_MAD_K = 3.0                            # windowed outlier trim threshold

# Feedstock taxonomy: davis ``Feedstock Composition`` label → normalized name.
# Exact (case-insensitive) match — "wood" is its OWN category, not soft/hard.
_FEEDSTOCK_MAP = {
    "soft wood": "softwood",
    "hard wood": "hardwood",
    "grass": "grass",
    "manure": "manure",
    "corn stover": "corn_stover",
    "wood": "wood",
}
#: First-class feedstocks eligible for per-group H/C and O/C overrides.
VALID_FEEDSTOCKS = tuple(sorted(set(_FEEDSTOCK_MAP.values())))

# Per-feedstock override is emitted only for these properties and only if a
# group passes the sufficiency gate below.
_OVERRIDE_PROPS = ("H_C_ratio", "O_C_ratio")
_GATE_MIN_N = 30
_GATE_MIN_TSPAN = 300.0

# Davis column indices (0-based). See module docstring / source header.
_DAVIS = {
    "feedstock": 1,
    "temp": 10,
    "ash_pct": 11,
    "C_pct": 13,        # "Total C (%) Used in Plot" — per project decision
    "N_pct": 14,
    "H_pct": 15,
    "O_pct": 16,
    "pH": 17,
    "ec_dsm": 18,
    "VM_pct": 41,
    "surface_area_m2g": 42,
}
# Reference properties read directly from davis columns (besides computed ratios).
_DAVIS_DIRECT = ("C_pct", "H_pct", "O_pct", "N_pct", "ash_pct", "pH",
                 "ec_dsm", "VM_pct", "surface_area_m2g")


# ---------------------------------------------------------------------------
# Small numeric helpers (numpy only)
# ---------------------------------------------------------------------------

def _to_float(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    x = x.strip().rstrip("%").strip()
    if x == "":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def _weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cw = np.cumsum(w)
    if cw[-1] <= 0:
        return float(np.median(values))
    return float(v[np.searchsorted(cw, 0.5 * cw[-1])])


def _mad_trim(t: np.ndarray, y: np.ndarray,
              window: float = 75.0, k: float = _MAD_K) -> np.ndarray:
    """Boolean mask keeping points within k·MAD of the local (windowed) median."""
    keep = np.ones(len(y), dtype=bool)
    for i in range(len(y)):
        local = np.abs(t - t[i]) <= window
        yl = y[local]
        if yl.size < 4:
            continue
        med = np.median(yl)
        mad = np.median(np.abs(yl - med))
        if mad > 0 and abs(y[i] - med) > k * mad:
            keep[i] = False
    return keep


def _smooth(t: np.ndarray, y: np.ndarray) -> dict:
    """Gaussian-kernel weighted-median smoothing of (t, y) onto :data:`_GRID`."""
    keep = _mad_trim(t, y)
    t, y = t[keep], y[keep]
    grid = np.array(_GRID, dtype=float)
    mean, spread, n = [], [], []
    for tc in grid:
        w = np.exp(-0.5 * ((t - tc) / _BANDWIDTH) ** 2)
        wsum = w.sum()
        if wsum < 1e-6 or t.size == 0:
            mean.append(float("nan")); spread.append(float("nan")); n.append(0)
            continue
        m = _weighted_median(y, w)
        sp = _weighted_median(np.abs(y - m), w)  # weighted MAD
        mean.append(m); spread.append(sp)
        n.append(int(np.count_nonzero(np.abs(t - tc) <= 2 * _BANDWIDTH)))
    mean = _fill_and_clamp(np.array(mean))
    return {
        "mean": [round(float(v), 5) for v in mean],
        "spread": [None if np.isnan(s) else round(float(s), 5) for s in spread],
        "n": n,
        "t_min": float(t.min()) if t.size else None,
        "t_max": float(t.max()) if t.size else None,
    }


def _fill_and_clamp(arr: np.ndarray) -> np.ndarray:
    """Replace NaN grid points by nearest finite neighbour (endpoint clamp)."""
    arr = arr.copy()
    finite = np.where(~np.isnan(arr))[0]
    if finite.size == 0:
        return np.zeros_like(arr)
    for i in range(len(arr)):
        if np.isnan(arr[i]):
            j = finite[np.argmin(np.abs(finite - i))]
            arr[i] = arr[j]
    return arr


# ---------------------------------------------------------------------------
# Offline build
# ---------------------------------------------------------------------------

def _read_csv(path: Path) -> List[List[str]]:
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.reader(fh))[1:]  # drop header


def _davis_records(rows: List[List[str]]) -> Dict[str, List[Tuple[float, float, Optional[str]]]]:
    """Return {property: [(temp, value, feedstock_norm), ...]} from davis rows."""
    out: Dict[str, list] = {p: [] for p in
                            ("H_C_ratio", "O_C_ratio", *_DAVIS_DIRECT)}
    ci = _DAVIS
    for r in rows:
        if len(r) <= ci["surface_area_m2g"]:
            continue
        temp = _to_float(r[ci["temp"]])
        if temp is None:
            continue
        feed = _classify_feedstock(r[ci["feedstock"]])
        c = _to_float(r[ci["C_pct"]])
        h = _to_float(r[ci["H_pct"]])
        o = _to_float(r[ci["O_pct"]])
        if c and h and c > 0:
            out["H_C_ratio"].append((temp, (h / _M_H) / (c / _M_C), feed))
        if c and o and c > 0:
            out["O_C_ratio"].append((temp, (o / _M_O) / (c / _M_C), feed))
        for prop in _DAVIS_DIRECT:
            v = _to_float(r[ci[prop]])
            if v is not None:
                out[prop].append((temp, v, feed))
    return out


def _fit_aromaticity(arom_rows: List[List[str]]) -> dict:
    """Linear (clamped) fit aromaticity% = a + b·(H/C) from biochar_data.csv."""
    HC, AR = 4, 16
    pts = []
    for r in arom_rows:
        if len(r) <= AR:
            continue
        hc, ar = _to_float(r[HC]), _to_float(r[AR])
        if hc is not None and ar is not None:
            pts.append((hc, ar))
    x = np.array([p[0] for p in pts]); y = np.array([p[1] for p in pts])
    b, a = np.polyfit(x, y, 1)  # slope, intercept
    pred = a + b * x
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return {"a": round(float(a), 5), "b": round(float(b), 5),
            "n": len(pts), "r2": round(r2, 4),
            "hc_min": round(float(x.min()), 4), "hc_max": round(float(x.max()), 4)}


def _dedup_report(davis_rows, arom_rows) -> int:
    """Count biochar_data rows whose (Temp±10, H/C±0.02) match a davis row."""
    ci = _DAVIS
    keys = set()
    for r in davis_rows:
        if len(r) <= ci["H_pct"]:
            continue
        t = _to_float(r[ci["temp"]]); c = _to_float(r[ci["C_pct"]]); h = _to_float(r[ci["H_pct"]])
        if t and c and h and c > 0:
            keys.add((round(t / 10) * 10, round((h / _M_H) / (c / _M_C), 2)))
    dups = 0
    for r in arom_rows:
        if len(r) <= 5:
            continue
        t, hc = _to_float(r[5]), _to_float(r[4])
        if t is not None and hc is not None and (round(t / 10) * 10, round(hc, 2)) in keys:
            dups += 1
    return dups


def build_model(output_path: Optional[Path] = None) -> Path:
    """
    Build the temperature-model JSON artifact from the bundled CSVs.

    Dev/offline entry point::

        python -c "from biochar.temperature_model import build_model; build_model()"

    Returns the path written. Prints a short QC report (dedup count, per-property
    point counts, aromaticity fit R²).
    """
    davis_rows = _read_csv(_DAVIS_CSV)
    arom_rows = _read_csv(_AROM_CSV)

    dups = _dedup_report(davis_rows, arom_rows)
    print(f"[build] davis rows={len(davis_rows)}  aromaticity rows={len(arom_rows)}")
    print(f"[build] davis is definitive; {dups} biochar_data rows duplicate davis "
          f"(Temp±10, H/C±0.02) — not merged (aromaticity sourced separately)")

    recs = _davis_records(davis_rows)

    properties: Dict[str, dict] = {}
    for prop, triples in recs.items():
        if not triples:
            continue
        t = np.array([x[0] for x in triples], dtype=float)
        y = np.array([x[1] for x in triples], dtype=float)
        properties[prop] = _smooth(t, y)
        print(f"[build]   {prop:18s} n={len(triples):4d}")

    # Per-feedstock overrides (H/C, O/C only) for gate-passing first-class groups.
    overrides: Dict[str, dict] = {}
    for fs in VALID_FEEDSTOCKS:
        for prop in _OVERRIDE_PROPS:
            triples = [x for x in recs[prop] if x[2] == fs]
            if len(triples) < _GATE_MIN_N:
                continue
            t = np.array([x[0] for x in triples], dtype=float)
            if float(t.max() - t.min()) < _GATE_MIN_TSPAN:
                continue
            y = np.array([x[1] for x in triples], dtype=float)
            overrides.setdefault(fs, {})[prop] = _smooth(t, y)
        if fs in overrides:
            print(f"[build]   override[{fs}] props={list(overrides[fs])}")

    arom = _fit_aromaticity(arom_rows)
    print(f"[build] aromaticity = {arom['a']} + {arom['b']}·(H/C)  "
          f"(n={arom['n']}, R²={arom['r2']})")

    artifact = {
        "schema": 1,
        "grid_celsius": _GRID,
        "bandwidth_celsius": _BANDWIDTH,
        "properties": properties,
        "feedstock_overrides": overrides,
        "aromaticity_fit": arom,
        "provenance": {
            "primary_source": "UC Davis Biochar Database",
            "primary_url": "https://biochar.ucdavis.edu/",
            "aromaticity_source": "biochar_data.csv (NMR aromaticity vs H/C)",
            "method_ref": "Wood, Mašek & Erastova, Cell Rep. Phys. Sci. 5(7) 2024, "
                          "DOI 10.1016/j.xcrp.2024.102036",
            "davis_rows": len(davis_rows),
            "duplicate_rows_vs_biochar_data": dups,
            "feedstock_note": "woody/lignocellulosic-dominated; unmapped labels pooled",
        },
    }
    path = output_path or _MODEL_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(artifact, fh, indent=2)
    print(f"[build] wrote {path}")
    return path


def _classify_feedstock(raw: str) -> Optional[str]:
    """Normalize a davis feedstock label; return None if not a first-class group."""
    return _FEEDSTOCK_MAP.get(raw.strip().lower())


# ---------------------------------------------------------------------------
# Runtime
# ---------------------------------------------------------------------------

class TemperatureModel:
    """
    Query the bundled temperature × feedstock property model.

    Args:
        model_path: Path to the JSON artifact (defaults to the bundled file).
    """

    def __init__(self, model_path: Optional[Path] = None):
        path = Path(model_path) if model_path else _MODEL_PATH
        with open(path, encoding="utf-8") as fh:
            self._m = json.load(fh)
        self._grid = np.array(self._m["grid_celsius"], dtype=float)

    # -- property access ----------------------------------------------------

    def _series(self, prop: str, feedstock: Optional[str]):
        """Return (mean_array, t_min, t_max) honouring a feedstock override."""
        if (feedstock and prop in _OVERRIDE_PROPS
                and feedstock in self._m["feedstock_overrides"]
                and prop in self._m["feedstock_overrides"][feedstock]):
            s = self._m["feedstock_overrides"][feedstock][prop]
            return np.array(s["mean"], dtype=float), s.get("t_min"), s.get("t_max")
        s = self._m["properties"].get(prop)
        if s is None:
            return None, None, None
        return np.array(s["mean"], dtype=float), s.get("t_min"), s.get("t_max")

    def predict(self, temperature: float, prop: str,
                feedstock: Optional[str] = None) -> float:
        """
        Predict *prop* at *temperature* (°C). Uses the feedstock-specific curve
        only for H/C and O/C when that group has one and *temperature* is within
        its support; otherwise the pooled curve. Clamped to the grid ends.
        """
        use_fs = feedstock
        # Fall back to pooled if the override doesn't cover this temperature.
        if feedstock and prop in _OVERRIDE_PROPS:
            ov = self._m["feedstock_overrides"].get(feedstock, {}).get(prop)
            if ov is None or not (_in_range(temperature, ov.get("t_min"), ov.get("t_max"))):
                use_fs = None
        mean, _, _ = self._series(prop, use_fs)
        if mean is None:
            raise KeyError(f"Unknown property: {prop!r}")
        return float(np.interp(temperature, self._grid, mean))

    def aromaticity_from_hc(self, hc: float) -> float:
        f = self._m["aromaticity_fit"]
        return float(min(100.0, max(0.0, f["a"] + f["b"] * hc)))

    def composition(self, temperature: float,
                    feedstock: Optional[str] = None) -> Dict[str, float]:
        """The three generator targets at *temperature* (and *feedstock*)."""
        hc = self.predict(temperature, "H_C_ratio", feedstock)
        oc = self.predict(temperature, "O_C_ratio", feedstock)
        return {
            "H_C_ratio": hc,
            "O_C_ratio": oc,
            "aromaticity_percent": self.aromaticity_from_hc(hc),
        }

    def predict_all(self, temperature: float,
                    feedstock: Optional[str] = None) -> Dict[str, float]:
        """Every modeled property at *temperature* (reference query)."""
        out = self.composition(temperature, feedstock)
        for prop in self._m["properties"]:
            if prop not in out:
                out[prop] = self.predict(temperature, prop, feedstock)
        return out

    def get_valid_range(
        self, feedstock: Optional[str] = None
    ) -> Optional[Tuple[float, float]]:
        """
        Return ``(T_min, T_max)`` (°C) of the data used to fit the model curves.

        Uses the feedstock-specific override range when one exists for H/C or O/C;
        otherwise returns the pooled range.  Returns ``None`` if the model artifact
        contains no range information (older schema).
        """
        for prop in _OVERRIDE_PROPS:
            if (feedstock and feedstock in self._m["feedstock_overrides"]
                    and prop in self._m["feedstock_overrides"][feedstock]):
                s = self._m["feedstock_overrides"][feedstock][prop]
            else:
                s = self._m["properties"].get(prop, {})
            t_min = s.get("t_min")
            t_max = s.get("t_max")
            if t_min is not None and t_max is not None:
                return float(t_min), float(t_max)
        return None

    @property
    def feedstocks(self) -> Tuple[str, ...]:
        return VALID_FEEDSTOCKS

    @property
    def provenance(self) -> dict:
        return dict(self._m.get("provenance", {}))


def _in_range(t: float, lo, hi) -> bool:
    return lo is not None and hi is not None and lo <= t <= hi


# Cached default model (loaded once; runtime needs only numpy + json).
_DEFAULT_MODEL: Optional[TemperatureModel] = None


def get_default_model() -> TemperatureModel:
    """Return a process-wide cached :class:`TemperatureModel` for the bundled artifact."""
    global _DEFAULT_MODEL
    if _DEFAULT_MODEL is None:
        _DEFAULT_MODEL = TemperatureModel()
    return _DEFAULT_MODEL


def get_valid_range(feedstock: Optional[str] = None) -> Optional[Tuple[float, float]]:
    """
    Return ``(T_min, T_max)`` data range for the given *feedstock* (or pooled).

    Example::

        from biochar.temperature_model import get_valid_range
        lo, hi = get_valid_range("softwood")  # e.g. (200.0, 900.0)
    """
    if feedstock is not None and feedstock not in VALID_FEEDSTOCKS:
        raise ValueError(
            f"feedstock must be one of {VALID_FEEDSTOCKS} or None, got {feedstock!r}"
        )
    return get_default_model().get_valid_range(feedstock)


def properties(temperature: float, feedstock: Optional[str] = None) -> Dict[str, float]:
    """
    Convenience: all modeled biochar properties at a pyrolysis *temperature* (°C),
    optionally for a *feedstock* (one of :data:`VALID_FEEDSTOCKS`).

    Example::

        import biochar
        biochar.properties(600, feedstock="softwood")
    """
    if feedstock is not None and feedstock not in VALID_FEEDSTOCKS:
        raise ValueError(
            f"feedstock must be one of {VALID_FEEDSTOCKS} or None, got {feedstock!r}"
        )
    return get_default_model().predict_all(temperature, feedstock)


# ---------------------------------------------------------------------------
# Refit-and-compare helper (for when Charchive DB arrives)
# ---------------------------------------------------------------------------

def compare_models(path_a: Path, path_b: Path) -> Dict[str, dict]:
    """
    Compare two model artifacts on their shared grid; report per-property mean
    absolute and max prediction delta. Used to evaluate a refit (e.g. davis-only
    vs davis+Charchive).
    """
    a = json.load(open(path_a)); b = json.load(open(path_b))
    out = {}
    shared = set(a["properties"]) & set(b["properties"])
    for prop in sorted(shared):
        ma = np.array(a["properties"][prop]["mean"], dtype=float)
        mb = np.array(b["properties"][prop]["mean"], dtype=float)
        n = min(len(ma), len(mb))
        d = np.abs(ma[:n] - mb[:n])
        out[prop] = {"mean_abs_delta": round(float(d.mean()), 5),
                     "max_abs_delta": round(float(d.max()), 5)}
    return out
