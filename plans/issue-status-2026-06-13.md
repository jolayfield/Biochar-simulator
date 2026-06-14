# Biochar Simulator — Issue Status & Outstanding Work Plan

**Generated:** 2026-06-13 (automated)  
**Version:** 0.2.0 (Alpha)  
**Open GitHub issues:** 0 (all #1–#5 closed)

---

## Summary

No open GitHub issues exist. The internal backlog (`plans/issue-backlog-2026-06-09.md`) proposed 11 issues (ISSUE-A through ISSUE-K); all were addressed in the working diff committed since that backlog was written.

---

## Backlog Implementation Status

| Issue | Title | Status |
|-------|-------|--------|
| ISSUE-A | Global random seed isolation | ✅ Done — instance-level `self._rng = random.Random(seed)` throughout; global `random.seed()` removed |
| ISSUE-B | Amorphous pore packing diagnostics | ✅ Done — `_is_clash_free` → `_min_separation`; error includes `"Best achieved separation: X.XX Å"` |
| ISSUE-C | `max_ether_span` default inconsistency | ✅ Done — convenience wrapper default changed from `5` → `None` (defers to `GeneratorConfig` default of 3) |
| ISSUE-D | Temperature model validity range | ✅ Done — `TemperatureModel.get_valid_range()` + module-level function; warning when temp is outside data range |
| ISSUE-E | `functional_groups` count validation | ✅ Done — warning when requested count > 1.5× feasible edge-site estimate |
| ISSUE-F | Valence errors buried in validation | ✅ Done — dedicated `"Valence Issues: N"` line in `print_results()`, always shown |
| ISSUE-G | `max_ether_span < 3` not enforced | ✅ Done — `ValueError` guard in `OxygenAssigner.__init__` and `GeneratorConfig.__post_init__` |
| ISSUE-H | Box padding validation | ✅ Done — `box_padding_xy/z ≤ 0` raises; values > 10 nm emit unit-error warning |
| ISSUE-I | Batch progress callbacks | ✅ Done — `generate_biochar_series()` gains `progress_callback` and `on_error` params |
| ISSUE-J | `defect_fraction` → ring composition | ✅ Done — `ring_composition` on `CarbonSkeleton`, `BiocharGenerator`, and `BiocharResult` |
| ISSUE-K | Static OPLS charge neutrality | ✅ Done — debug log for residual > 0.01e; `test_charges_sum_to_zero` asserts neutrality |

---

## Remaining Items

### 1. ISSUE-J follow-up: README section for `defect_fraction`

No documentation exists explaining that `defect_fraction` is a request probability, not a guaranteed pentagon fraction — parity constraints reject some placements.

**Suggested addition to `README.md`:**

> ### Understanding `defect_fraction`
>
> `defect_fraction=0.15` does **not** mean exactly 15% of rings will be pentagons.  
> During skeleton growth each ring-addition event has a 15% chance of *attempting*  
> a pentagon; parity and adjacency constraints reject many of these attempts.  
> Typical observed pentagon fractions:
>
> | `defect_fraction` | ~50 C skeleton | ~100 C skeleton |
> |-------------------|---------------|-----------------|
> | 0.05              | 0–1 pentagons | 1–3 pentagons   |
> | 0.15              | 1–3 pentagons | 2–6 pentagons   |
> | 0.30              | 2–5 pentagons | 4–10 pentagons  |
>
> Use `result.ring_composition` to inspect the actual counts after generation.

**Effort:** ~30 minutes.

---

### 2. ISSUE-B follow-up: Auto-relax fallback for amorphous packing (optional)

The improved diagnostic error message is in place. The optional `amorphous_fallback="slit"` param on `SurfaceConfig` (silently degrade to slit pore on failure) was not implemented. This is a nice-to-have and can be filed as a future GitHub issue.

---

## Suggested GitHub Issues to File

| Title | Label | Effort |
|-------|-------|--------|
| `defect_fraction` docs: README section explaining pentagon fraction vs. request probability | `docs` | 30 min |
| `SurfaceConfig`: add `amorphous_fallback` param to degrade gracefully on packing failure | `enhancement` | 2–3 hr |
