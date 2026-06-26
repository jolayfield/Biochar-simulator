# Biochar Simulator — Issue Status & Outstanding Work Plan

**Generated:** 2026-06-26 (automated)
**Version:** 0.3.0 (bumped this run)
**Open GitHub issues:** 0 (all #1–#5 closed)

---

## Summary

No open GitHub issues. All previously tracked items are resolved. This run found 4 issues
related to unreleased work that had accumulated since the v0.2.0 git tag; all 4 were
addressed within this automated run.

---

## Changes Since 2026-06-20 Status Report

**Completed (committed before this run):**
- `amorphous_fallback="slit"` on `SurfaceConfig` — committed in `069b86c`
- Stale module docstring in `surface_builder.py` — fixed in `069b86c`
- Example notebook `examples/hardwood_400C_series.ipynb` — committed in `069b86c`

**Addressed this run:**
- Version bumped `0.2.0` → `0.3.0` in `pyproject.toml` and `biochar/__init__.py`
- `CHANGELOG.md` updated with entries for v0.1.5, v0.2.0, and v0.3.0
- `QMChargeError` now exported from `biochar.__init__` (`__all__`)
- README updated: `charge_method` row added to GeneratorConfig parameter table; new
  "Partial charge methods" section documents all three backends with requirements

---

## Proposed GitHub Issues

None outstanding. All work from the June 9 backlog (ISSUE-A through ISSUE-K), the
June 13 plan, and the June 20 plan is now implemented and the repo is tagged and
documented through v0.3.0.

### Future candidates (not urgent)

- **LBCC bond-charge corrections** — `qm_charges.py` notes that 1.14*CM1A-LBCC
  corrections are not applied. Low priority; the base 1.14*CM1A model already matches
  LigParGen's default.
- **N-doped amorphous structures** — combining `pore_type="amorphous"` with
  pyridinic/graphitic N has not been explicitly tested. May surface edge cases in
  `carbon_skeleton.py` when N substitution interacts with ring parity constraints.
