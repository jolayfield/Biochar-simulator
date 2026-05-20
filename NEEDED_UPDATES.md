# Needed Updates

Analysis of the biggest gaps in the Biochar Simulator codebase, ordered by scientific/engineering impact.

---

## 1. N-Doping via Amino Groups

**Problem:** The OPLS atom types for nitrogen (`N`, `NT`) are defined in `constants.py` and `opls_typing.py` with placeholder comments, but there is no way to actually add nitrogen atoms to structures. Biochar from protein-rich feedstocks contains 1–5% nitrogen by weight, and N-doped biochar is an active research area.

**Files changed:** `biochar/constants.py`, `biochar/heteroatom_assignment.py`, `biochar/opls_typing.py`

**Implementation:** Add an `amino` functional group (`Ar–NH₂`, analogous to `phenolic` but nitrogen) with full OPLS-AA parameters (aniline type: `NA`/`HNA`). Extend `OxygenAssigner._place_group` to handle the new type; add `num_nitrogens` and `N_C_ratio` tracking to `CompositionResult`.

---

## 2. `molecular_formula` and `molecular_weight` on `CompositionResult`

**Problem:** `CompositionResult` stores raw atom counts but provides no molecular formula string or molecular weight. Every user must compute these themselves — error-prone and unnecessary boilerplate.

**Files changed:** `biochar/heteroatom_assignment.py`

**Implementation:** Add `@property molecular_formula` (returns Hill-order string, e.g. `"C24H12O2N1"`) and `@property molecular_weight` (returns float in g/mol) to `CompositionResult`. Also add `num_nitrogens: int = 0` and `N_C_ratio: float = 0.0` fields (needed for item 1 above and the formula/weight computations).

---

## 3. Integration Tests for `generate_biochar_series()`

**Problem:** `generate_biochar_series()` is a key public API that orchestrates batch generation and creates a combined GROMACS topology. It has **zero test coverage** in the pytest suite. A regression here could silently break batch workflows.

**Files changed:** `tests/test_generator.py`

**Implementation:** Add a `TestBiocharSeries` class with tests covering: successful 2-structure generation, correct returned keys, file existence on disk, combined topology creation, error handling for missing `molecule_name`, and the new amino-group functional_groups path.

---

## 4. Replace `print()` with `logging` Module

**Problem:** The entire pipeline uses bare `print()` calls for progress and warnings. As a library, this is incorrect: print output cannot be silenced, filtered, or redirected. Library code must use `logging` so the calling application controls output.

**Files changed:** `biochar/__init__.py`, `biochar/biochar_generator.py`, `biochar/heteroatom_assignment.py`

**Implementation:**
- Add `logging.getLogger(__name__).addHandler(logging.NullHandler())` to `biochar/__init__.py` (standard library best practice).
- Add `logger = logging.getLogger(__name__)` to each module and replace `print()` with `logger.info()` / `logger.warning()` throughout the generation pipeline.
- Keep `BiocharGenerator.print_summary()` as-is — it is explicitly a display function.

---

## 5. `GeneratorConfig.to_dict()` / `from_dict()` for Config Serialization

**Problem:** There is no way to save or reload a `GeneratorConfig`. This makes reproducibility harder: to share exact generation parameters, users must manually reconstruct Python code. A round-trippable dict representation enables JSON storage and programmatic config management.

**Files changed:** `biochar/biochar_generator.py`

**Implementation:** Add `to_dict() -> dict` (returns a JSON-serializable dict, converting `numpy.ndarray` box_size to a list) and `from_dict(cls, d: dict) -> GeneratorConfig` classmethod to `GeneratorConfig`.

---

## 6. CLI Entry Point (`biochar-gen`)

**Problem:** There is no command-line interface. Every user must write a Python script to generate even a single structure. A simple CLI drastically lowers the barrier to use for wet-lab collaborators and quick one-off generations.

**Files changed (new):** `biochar/cli.py`  
**Files changed:** `pyproject.toml`

**Implementation:** Add `biochar/cli.py` with an argparse-based `main()` function exposing the most common `generate_biochar` parameters (`--carbons`, `--hc-ratio`, `--oc-ratio`, `--defects`, `--name`, `--seed`, `--output-dir`, `--phenolic`, `--carboxyl`, `--ether`, `--amino`). Register it as `biochar-gen` in `pyproject.toml` `[project.scripts]`.

---

## 7. Expose `generate_biochar_series` in Public API

**Problem:** `generate_biochar_series` is exported from `biochar_generator.py` but not re-exported in `biochar/__init__.py`, so `from biochar import generate_biochar_series` fails with an `ImportError` even though the function exists.

**Files changed:** `biochar/__init__.py`

**Implementation:** Add `generate_biochar_series` to the imports and `__all__` list in `__init__.py`.

---

## Out of Scope (Future Work)

The following are genuine gaps but require substantially more work and are left for future phases:

- **Pyridinic / pyrrolic / graphitic N**: Ring-substituting nitrogen requires modifying the carbon skeleton topology, not just pendant group attachment. Phase 3 work.
- **Amorphous porous packing**: Already flagged as Phase 2 in the codebase (`surface_builder.py:106`). Requires a 3D packing algorithm.
- **S-doping (thiol, thioether)**: OPLS types exist as placeholders; requires the same pipeline extension as amino groups but for S chemistry.
- **ML-based charge refinement**: Mentioned in roadmap; out of scope for structural generation.
