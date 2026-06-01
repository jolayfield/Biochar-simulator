# Plan: Temperature × feedstock composition model (davis dataset primary)

> Status: **REVISED** with `data/davis-biochar-db.csv` as the definitive source.
> Supersedes the earlier wood-only draft. Feedstock splitting is now first-class.

## Context

Goal: drive the generator's composition targets (H/C, O/C, aromaticity) from
real data as a continuous function of **pyrolysis temperature** and
**feedstock**, and expose all other measured properties as a reference query.

Two datasets are involved:

| Dataset | Rows | Role |
|---|---|---|
| `data/davis-biochar-db.csv` | 1176 | **PRIMARY / definitive.** Clean `Feedstock Composition` column; Temp 100–1000 °C; C/H/O/N %, ash, pH, surface area, VM, EC, CEC, P, metals, PAH content, … |
| `data/biochar_data.csv` | 454 | **Supplement only** — provides the one thing davis lacks: NMR aromaticity, used to calibrate `aromaticity = f(H/C)`. |

### Deduplication (measured)

- All **50/50** source papers in the old CSV also appear in davis (davis is a
  citation-level superset).
- **71 %** of old H/C rows are near-duplicates of davis rows (same Temp ±10 °C,
  H/C ±0.02).
- **Rule (per user): davis is authoritative on any duplicate.** Architecturally
  guaranteed: each property is sourced from **exactly one** file — H/C, O/C, and
  all reference properties from **davis**; aromaticity calibration from the old
  CSV only. So duplicate rows cannot double-count or corrupt any model. The
  build still runs an explicit duplicate detector (match on author+year+Temp, or
  Temp+H/C) and **reports** the count for QC.

## Key data facts (drive the design)

- **H/C and O/C are not stored** in davis — compute molar ratios from mass %:
  `H/C = (H%/1.008)/(C%/12.011)`, `O/C = (O%/15.999)/(C%/12.011)`, using
  `Total C (%) Used in Plot` (col 13), `H (%)` (15), `O (%)` (16).
  → **750 rows** give computable H/C (100–1000 °C); **633** give O/C.
- **Clean feedstock taxonomy** (col 1) — no ID parsing. Per-feedstock H/C counts:

  | feedstock | n (H/C) | n (O/C) | T-range | gate ≥30 & ≥300 °C |
  |---|---|---|---|---|
  | Grass | 150 | 126 | 100–900 | ✅ first-class |
  | Other/ Mixed | 153 | 134 | 200–900 | ✅ (pooled-ish bucket) |
  | Soft Wood | 106 | 88 | 100–950 | ✅ |
  | Hard Wood | 87 | 82 | 200–1000 | ✅ |
  | Manure | 77 | 62 | 200–800 | ✅ (non-lignocellulosic — warn) |
  | Corn stover | 39 | 30 | 300–850 | ✅ |
  | Nutshell | 25 | 23 | 250–900 | ⚠ borderline → low-confidence |
  | Hull | 21 | 21 | 300–800 | ⚠ borderline |
  | Wood | 19 | 17 | 300–900 | ⚠ borderline |
  | Sludge | 17 | 11 | 300–900 | ⚠ non-lignocellulosic (H/C≈1.0) — warn |
  | Algae | 16 | 1 | 300–750 | ⚠ non-lignocellulosic |
  | Pomace | 13 | 13 | 150–700 | ⚠ borderline |

- **Aromaticity ← H/C** (old CSV, 20 paired pts, **r = −0.956**, monotonic
  H/C 0.03→100 % … 1.6→~25 %). Fit `aromaticity% = clamp(f(H/C), 0, 100)` once,
  then derive aromaticity from the davis H/C — making it feedstock/temperature
  aware indirectly.
- **Lignocellulosic appropriateness:** the generator builds aromatic PAH sheets.
  Wood/grass/nutshell/hull/corn-stover fit well; **manure/sludge/algae** are
  high-H/C, low-aromaticity — supported but the generator **warns** when the
  derived `aromaticity_percent < 50` that a PAH model is a poor representation.

## Design

**Offline fit → compact JSON artifact → runtime `numpy.interp`** (runtime needs
only numpy; core functionality, not an optional extra).

### Stage 1 — build (`biochar/temperature_model.py::build_model()`)

1. Parse both CSVs with stdlib `csv` (quoted citations contain commas).
2. **davis**: compute molar H/C, O/C; collect every numeric property with Temp.
3. Run the duplicate detector vs old CSV; report counts (davis wins).
4. Per (property, feedstock-or-pooled): robust **MAD outlier trim**, then
   **Gaussian-kernel weighted-median** smoothing onto a 100–1000 °C / 25 °C grid
   (bandwidth ≈ 50 °C); store `mean`, `spread`, `n`, `[t_min, t_max]`; clamp
   outside support.
5. Feedstock overrides for H/C & O/C only, for groups passing the gate
   (≥30 pts, ≥300 °C span). Others omitted (→ pooled fallback) or flagged low-n.
6. Fit `aromaticity = f(H/C)` from old CSV; store coefficients.
7. Write `biochar/data/temperature_model.json`.

### Stage 2 — runtime (`TemperatureModel`)

- `predict(temperature, prop, feedstock=None)` — interp, clamped, feedstock
  override iff present+in-range else pooled.
- `composition(temperature, feedstock=None)` → `{H_C_ratio, O_C_ratio,
  aromaticity_percent}` (aromaticity via the H/C calibration).
- `predict_all(...)` → every reference property.

### Integration (`biochar/biochar_generator.py`, `constants.py`, `__init__.py`)

- `GeneratorConfig`: add `temperature` and `feedstock` (validated against the
  normalized taxonomy); composition fields default `None`; `__post_init__`
  precedence **explicit > (temperature,feedstock)-derived > old default**.
- Feedstock label normalization: `"Soft Wood"→softwood`, `"Hard Wood"→hardwood`,
  `"Grass"→grass`, `"Corn stover"→corn_stover`, etc.
- `PRESETS` (constants) become temperature labels (`BC600→600`).
- `biochar.properties(temperature, feedstock=None)` query; warn when derived
  aromaticity < 50 %.

## Critical files

| File | Change |
|---|---|
| `biochar/temperature_model.py` | **New** — dual-source `build_model()`, dedup, smoother, `aromaticity=f(H/C)` fit, `TemperatureModel` |
| `biochar/data/temperature_model.json` | **New** — committed artifact |
| `data/davis-biochar-db.csv`, `data/biochar_data.csv` | Source data (repo-root, not shipped) |
| `biochar/biochar_generator.py` | `GeneratorConfig.temperature`/`.feedstock`; resolution; `properties()`; low-aromaticity warning |
| `biochar/constants.py` | `PRESETS` labels; feedstock taxonomy map |
| `biochar/__init__.py` | export `TemperatureModel`, `properties` |
| `pyproject.toml` | `data/*.json` package-data |
| `tests/test_temperature_model.py` | **New** |

## Verification

```bash
python3 -c "from biochar.temperature_model import build_model; build_model()"  # writes JSON, prints dedup count + LOO-RMSE per property/feedstock
python3 -m pytest tests/test_temperature_model.py -v
python3 -m pytest tests/ -v   # full suite stays green
```

Tests: monotonic H/C↓ & aromaticity↑ with T; H/C(600)∈[0.2,0.6]; clamping;
`composition(450,"softwood")` ≠ `composition(450,"grass")`; out-of-range
feedstock → pooled; unknown feedstock raises; `GeneratorConfig(temperature=600,
feedstock="hardwood")` fills 3 fields, explicit override wins, bare config
unchanged (0.5/0.1/90); aromaticity derived from H/C tracks the r=−0.956 curve;
low-aromaticity warning fires for `feedstock="sludge"`; `properties(600)` returns
ash/pH/surface-area/etc.

## Extensibility — refitting when new datasets arrive

`build_model()` takes an **ordered list of source specs** (path + column map +
priority); highest priority wins on dedup (davis first, future sources slotted
in, old CSV last for aromaticity). Adding a source = one list entry + rerun the
build; the JSON regenerates and **no runtime/integration code changes**. A
`compare_models(json_a, json_b)` helper reports per-property mean/max prediction
delta and RMSE change, so refit-and-compare is one command. The `f(H/C)`
aromaticity fit is a tiny regression — trivial to refit as points are added.
Net: building now is **not** wasted work; the Charchive DB drops in later.

## Decisions confirmed (2026-05-31)

- ✅ **Davis data may be redistributed/bundled** with the package (user-approved).
  Still record the formal citation; redistribution itself is cleared.
- ✅ **Aromaticity derived from H/C** (`aromaticity = f(H/C)`, r=−0.956) — approved.
- ✅ **Non-lignocellulosic feedstocks** (manure/sludge/algae): **support with
  warnings** — expose them, but warn when derived `aromaticity_percent < 50` that
  the PAH model is a poor fit. (Not excluded.)
- ✅ **Wood et al. citation** resolved: Wood, Mašek & Erastova, *Cell Rep. Phys.
  Sci.* 5(7), 2024, DOI 10.1016/j.xcrp.2024.102036. The `data/*.pptx` are its
  figures, built from `biochar_data.csv` (no new data).
- ✅ **Package version** reconciled to **0.1.4** across all files.
- ✅ **Feedstock taxonomy** — first-class `feedstock=` options:
  `softwood` (Soft Wood), `hardwood` (Hard Wood), `grass` (Grass),
  `manure` (Manure), `corn_stover` (Corn stover), **`wood` (its own category,
  NOT merged into soft/hard — per user)**. Remaining labels (Other, Other/Mixed,
  Nutshell, Hull, Sludge, Algae, Pomace) **feed the pooled baseline** and are
  selectable but fall back to the **pooled curve with a low-confidence warning**.
- ✅ **H/C and O/C computed from `Total C (%) Used in Plot`** (col 13), not
  `C Organic (%)` (user-confirmed).
- ✅ **Low predicted aromaticity → clamp + warn**: when `f(H/C)` returns below the
  PAH engine's buildable floor, **clamp `aromaticity_percent` to that floor and
  emit a warning** (the generator builds near-fully-aromatic sheets).
- ✅ **Davis DB citation/URL:** **https://biochar.ucdavis.edu/** (UC Davis Biochar
  Database) — record in `constants.py` provenance + README/docs.

## Open items (remaining)

1. **Charchive DB (pending access)** — slots in as a higher-priority source;
   refit + `compare_models` against the davis-only artifact when available. If it
   carries direct NMR aromaticity, upgrade the `f(H/C)` proxy to real (possibly
   per-feedstock) aromaticity data.
2. **Aromaticity fit form** — linear-clamped vs logistic on ~20 pts (pick by R²
   during implementation; approach already approved above). *(implementation detail)*
