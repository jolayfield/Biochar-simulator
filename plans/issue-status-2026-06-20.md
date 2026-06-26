# Biochar Simulator — Issue Status & Outstanding Work Plan

**Generated:** 2026-06-20 (automated)
**Version:** 0.2.0
**Open GitHub issues:** 0 (all #1–#5 closed)

---

## Summary

No open GitHub issues. All five tracked features (#1 amorphous packing, #2 ring-N
doping, #3 S-doping, #4 ML charges, #5 named BiocharResult) are implemented and
passing in the 389-test suite.

---

## Changes Since 2026-06-13 Status Report

**Completed:**
- `defect_fraction` README docs section added (lines 95–107) — was listed as
  outstanding in the previous plan.

**Still outstanding:**
- `SurfaceConfig` `amorphous_fallback` graceful degradation param (see Issue A below).

**New item:**
- `examples/hardwood_400C_series.ipynb` is untracked in git.

---

## Proposed GitHub Issues

### Issue A — `SurfaceConfig`: add `amorphous_fallback` parameter

**Label:** `enhancement` | **Effort:** 2–3 hr

When `pore_type="amorphous"` cannot place a sheet within `max_attempts` it raises
`RuntimeError`. An optional `amorphous_fallback: str = None` on `SurfaceConfig`
(accepting `"slit"`) would degrade gracefully to slit-pore geometry and emit a
warning instead of crashing.

**Touchpoints:**
- `biochar/surface_builder.py` — `SurfaceConfig` dataclass + `_pack_amorphous` catch block
- `tests/test_surface_builder.py` — new test for the fallback path

---

### Issue B — Commit example notebook `hardwood_400C_series.ipynb`

**Label:** `docs` | **Effort:** 15 min

`examples/hardwood_400C_series.ipynb` demonstrates a 5-molecule series conditioned on
the UC Davis temperature-feedstock model (400 °C, hardwood). The notebook is untracked
in git and should be committed so it is discoverable in the repo.

**Action:**
```bash
git add examples/hardwood_400C_series.ipynb
git commit -m "docs: add hardwood 400°C series example notebook"
```

---

### Issue C — Stale module docstring in `surface_builder.py`

**Label:** `docs` | **Effort:** 5 min

`biochar/surface_builder.py:4` still reads "amorphous (Phase 2, deferred — see #1)"
but GitHub issue #1 is closed and the implementation is complete. Update the docstring
to remove the stale "deferred" note.
