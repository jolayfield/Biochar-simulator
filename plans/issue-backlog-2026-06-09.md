# Biochar Simulator — Issue Backlog & Implementation Plans

**Generated:** 2026-06-09  
**Current version:** 0.2.0 (Alpha)  
**Status:** All existing GitHub issues (#1–#5) are closed. This document proposes new issues based on codebase audit.

---

## Priority 1 — Medium severity

### ISSUE-A: Global random seed isolation (reproducibility bug)

**Problem:** Multiple pipeline classes (`PAHAssembler`, `OxygenAssigner`, `HydrogenAssigner`, `NitrogenSubstitutor`) call `random.seed()` on the *global* `random` module in their `__init__`. When `generate_biochar()` is called multiple times in the same process (e.g., `generate_biochar_series()`), each call resets the global seed, but any other random draws between instantiations break determinism. Two calls with `seed=42` can produce different molecules in the same process.

**Affected files:**
- `biochar/heteroatom_assignment.py` — `OxygenAssigner.__init__`, `HydrogenAssigner.__init__`, `NitrogenSubstitutor.__init__`
- `biochar/carbon_skeleton.py` — `PAHAssembler.__init__`

**Implementation plan:**
1. Replace every `random.seed(seed)` call with an instance-level `self._rng = random.Random(seed)`.
2. Replace bare `random.choice()`, `random.shuffle()`, `random.sample()` calls in those classes with `self._rng.*`.
3. For numpy usage, replace `np.random.*` calls with a `np.random.default_rng(seed)` generator stored as `self._np_rng`.
4. Pass the seed down explicitly from `GeneratorConfig` through `BiocharGenerator` to each sub-component.
5. Add a test in `tests/test_generator.py`: call `generate_biochar(seed=42)` twice in the same process and assert SMILES equality.

**Effort:** ~2 days. No API changes needed.

---

### ISSUE-B: Amorphous pore packing failure diagnostics

**Problem:** When `pore_type="amorphous"` exhausts `max_attempts`, the `RuntimeError` message offers general advice but no diagnostics — users don't know how close the packing came or which parameter is the bottleneck.

**Affected files:**
- `biochar/surface_builder.py` — `_generate_amorphous_pore()` (lines ~553–601)

**Implementation plan:**
1. Track the closest approach distance and last-rejected reason inside the placement loop.
2. Include these in the error message: `"Best achieved separation: {min_sep:.2f} Å (needed {target:.2f} Å)"`.
3. Add optional auto-relax: if packing fails at `min_separation=s`, retry once with `s * 0.95` and warn the user.
4. Add an optional `amorphous_fallback="slit"` param to `SurfaceConfig` that silently falls back to slit pore.
5. Update tests with a tight-box scenario to cover the new error message path.

**Effort:** ~1 day.

---

## Priority 2 — Low severity / high value

### ISSUE-C: `max_ether_span` default inconsistency between entry points

**Problem:** `GeneratorConfig.max_ether_span` defaults to `3` but the `generate_biochar()` convenience wrapper hardcodes `max_ether_span=5` in its signature. Callers using the convenience API get a different default than callers using the config object. The README documents 3 as the safe value.

**Affected files:**
- `biochar/biochar_generator.py` — `generate_biochar()` signature (~line 677) and `generate_biochar_series()` (~line 746)

**Implementation plan:**
1. Change the `generate_biochar()` signature default to `max_ether_span=None`.
2. When `None`, use `GeneratorConfig`'s own default (3) rather than specifying it explicitly.
3. Add a `GeneratorConfig.__post_init__()` check: warn if `max_ether_span > 4`.
4. Update docstrings on both functions to document `None` → "uses GeneratorConfig default (3)".

**Effort:** ~2 hours.

---

### ISSUE-D: Temperature model validity range — missing validation and documentation

**Problem:** The `temperature` parameter (derived from UC Davis Biochar Database model) has a valid range (~300–900 °C based on training data), but out-of-range inputs (e.g., `temperature=50`) are silently accepted and produce extrapolated, potentially nonsensical composition ratios. No warning is emitted.

**Affected files:**
- `biochar/biochar_generator.py` — `GeneratorConfig.__post_init__()`
- `biochar/temperature_model.py`

**Implementation plan:**
1. Add a `get_valid_range(feedstock)` function to `temperature_model.py` that returns the `(T_min, T_max)` range for each feedstock.
2. In `GeneratorConfig.__post_init__()`, if `temperature` is set and falls outside the valid range, emit a `logging.warning()`.
3. Add the valid range to the `--temperature` CLI help text.
4. Update the README temperature table to include a "Valid range" row.

**Effort:** ~3 hours.

---

### ISSUE-E: Input validation for `functional_groups` counts

**Problem:** `functional_groups={"phenolic": 100}` is silently accepted for a 50-carbon molecule. The placer will try to place 100 groups, succeed for far fewer, and only log warnings. Batch users may silently get under-oxygenated structures.

**Affected files:**
- `biochar/heteroatom_assignment.py` — `OxygenAssigner.__init__`
- `biochar/biochar_generator.py` — `BiocharGenerator.generate()`

**Implementation plan:**
1. After skeleton generation, compute an upper-bound estimate for each group type based on edge-atom count.
2. Before running `OxygenAssigner`, check that requested totals are ≤150% of the estimate.
3. In strict mode (`strict=True`), raise `ValidationError`; otherwise log a warning.
4. Add a helper method `estimate_max_functional_groups(mol)` to `heteroatom_assignment.py`.
5. Add tests for over-requested groups: verify warning emitted, structure still valid.

**Effort:** ~4 hours.

---

### ISSUE-F: Valence errors buried in validation summary

**Problem:** `ValidationEngine.validate_complete()` collects valence errors but they can be invisible in the printed summary if the molecule otherwise passes. Users in non-strict mode may miss them entirely.

**Affected files:**
- `biochar/validation.py` — `ValidationEngine.print_summary()`
- `biochar/biochar_generator.py` — post-generation summary call

**Implementation plan:**
1. Add a dedicated `"Valence Issues"` section to `print_summary()` that always prints, even if the count is zero.
2. Ensure any valence error always logs at `ERROR` level (not just stored in the results dict).
3. In strict mode, any valence error should raise `ValidationError` before returning the result.
4. Add a test asserting that a deliberately malformed molecule triggers the new log output.

**Effort:** ~2 hours.

---

### ISSUE-G: `max_ether_span < 3` not enforced

**Problem:** `OxygenAssigner` docstring states "Minimum enforced = 3" but the `__init__` has no actual guard. Passing `max_ether_span=1` or `max_ether_span=2` will not raise an error and may produce invalid ring geometries or runtime errors during placement.

**Affected files:**
- `biochar/heteroatom_assignment.py` — `OxygenAssigner.__init__`

**Implementation plan:**
1. Add `if max_ether_span < 3: raise ValueError(f"max_ether_span must be ≥ 3, got {max_ether_span}")` in `OxygenAssigner.__init__`.
2. Add a matching guard in `GeneratorConfig._validate()`.
3. Add a unit test covering the `ValueError`.

**Effort:** 30 minutes.

---

### ISSUE-H: Box padding values not validated in `SurfaceConfig`

**Problem:** Negative or extreme padding values (e.g., `box_padding_xy=-1`) are accepted silently and will produce invalid `.gro` boxes at export time.

**Affected files:**
- `biochar/surface_builder.py` — `SurfaceConfig._validate()`

**Implementation plan:**
1. Add guards: `box_padding_xy > 0`, `box_padding_z > 0`.
2. Warn if either exceeds 10 nm (likely a unit error by the user).
3. Add minimal tests for the new validation.

**Effort:** 30 minutes.

---

## Priority 3 — Nice-to-have / longer-term

### ISSUE-I: Batch generation progress reporting

**Problem:** `generate_biochar_series()` prints progress inline but has no support for progress bars, callbacks, or cancellation — making it unsuitable for GUI/web frontends or long overnight runs.

**Affected files:**
- `biochar/biochar_generator.py` — `generate_biochar_series()`

**Implementation plan:**
1. Add optional `progress_callback: Callable[[int, int, str], None]` parameter (called with `(completed, total, name)` after each molecule).
2. Add optional `progress_bar: bool = False` that activates a `tqdm` bar if tqdm is installed (soft dependency with try-import).
3. Add `on_error: Literal["raise", "skip", "warn"] = "raise"` to allow partial-batch success.
4. Add tests covering the callback and skip-on-error paths.

**Effort:** ~1 day.

---

### ISSUE-J: Document `defect_fraction` → actual pentagon fraction relationship

**Problem:** `defect_fraction=0.15` doesn't mean 15% of rings will be pentagons — parity constraints reject some pentagonal placements. Users have no way to know what actual defect density they'll get.

**Affected files:**
- `biochar/carbon_skeleton.py` — `PAHAssembler`
- `README.md`

**Implementation plan:**
1. Add a `ring_composition` property to the assembled skeleton returning `{"hexagons": N, "pentagons": M}`.
2. Include this in the `BiocharResult` named tuple or as a field on `BiocharGenerator`.
3. Add a README section "Understanding defect_fraction" with a table of `defect_fraction` vs. observed pentagon % for 50 C, 100 C structures.

**Effort:** ~4 hours.

---

### ISSUE-K: Static OPLS charge neutrality not guaranteed

**Problem:** The ML charge refinement in `ml_charges.py` enforces total charge = 0, but the static OPLS assignment in `ChargeAssigner` does not. For large molecules the total charge can drift by ±0.01e, which may cause GROMACS electrostatics warnings.

**Affected files:**
- `biochar/opls_typing.py` — `ChargeAssigner`

**Implementation plan:**
1. After assigning all charges, compute `total = sum(charges)` and distribute the residual evenly across all non-hydrogen atoms (or all atoms — document the choice).
2. Log a debug message if the residual before correction exceeds 0.01e.
3. Add a test asserting `abs(sum(charges)) < 1e-6` for a generated structure.

**Effort:** ~2 hours.

---

## Recommended GitHub issue creation order

| Priority | Proposed Issue | Effort | Suggests |
|----------|---------------|--------|----------|
| 1 | ISSUE-A: Seed isolation | 2 days | File as `bug` + `reproducibility` |
| 1 | ISSUE-B: Amorphous diagnostics | 1 day | File as `enhancement` + `ux` |
| 2 | ISSUE-G: max_ether_span guard | 30 min | File as `bug` (easy win) |
| 2 | ISSUE-H: box padding validation | 30 min | File as `bug` (easy win) |
| 2 | ISSUE-C: max_ether_span default | 2 hr | File as `bug` + `api` |
| 2 | ISSUE-F: Valence error visibility | 2 hr | File as `enhancement` |
| 2 | ISSUE-D: Temperature model range | 3 hr | File as `enhancement` + `docs` |
| 2 | ISSUE-E: functional_groups validation | 4 hr | File as `enhancement` |
| 3 | ISSUE-J: defect_fraction docs | 4 hr | File as `docs` |
| 3 | ISSUE-K: charge neutrality | 2 hr | File as `enhancement` |
| 3 | ISSUE-I: batch progress reporting | 1 day | File as `enhancement` |
