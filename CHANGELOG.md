# Biochar Simulator — Changelog

All significant changes to the Biochar Simulator project are documented here.

## [0.3.0] — June 26, 2026

### Added

- **`amorphous_fallback` on `SurfaceConfig`** — optional `amorphous_fallback="slit"` degrades gracefully to slit-pore geometry when amorphous packing exhausts `max_attempts`, emitting a warning instead of raising `RuntimeError`. `None` (default) preserves the original error-raising behaviour.
- **`QMChargeError` exported from `biochar`** — `from biochar import QMChargeError` now works without reaching into the submodule.
- **Example notebook** — `examples/hardwood_400C_series.ipynb` demonstrating a 5-molecule 400 °C hardwood series.

### Fixed (ISSUE-A through ISSUE-K)

- **Seed isolation** (ISSUE-A) — all pipeline classes now use instance-level `self._rng = random.Random(seed)`; global `random.seed()` calls removed. Two `generate_biochar(seed=42)` calls in the same process now always produce identical SMILES.
- **Amorphous packing diagnostics** (ISSUE-B) — `RuntimeError` on packing failure now includes `"Best achieved separation: X.XX Å (needed Y.YY Å)"`.
- **`max_ether_span` default inconsistency** (ISSUE-C) — `generate_biochar()` convenience wrapper default changed from `5` → `None`, deferring to `GeneratorConfig` default of 3.
- **Temperature model validity range** (ISSUE-D) — `TemperatureModel.get_valid_range(feedstock)` added; `GeneratorConfig` warns when `temperature` is outside the UC Davis training data range.
- **`functional_groups` count validation** (ISSUE-E) — warning emitted when requested count exceeds 1.5× feasible edge-site estimate.
- **Valence error visibility** (ISSUE-F) — dedicated `"Valence Issues: N"` line in `print_results()`, always shown even when count is zero.
- **`max_ether_span < 3` not enforced** (ISSUE-G) — `ValueError` guard added in `OxygenAssigner.__init__` and `GeneratorConfig.__post_init__`.
- **Box padding validation** (ISSUE-H) — `box_padding_xy/z ≤ 0` raises `ValueError`; values > 10 nm emit a unit-error warning.
- **Batch progress callbacks** (ISSUE-I) — `generate_biochar_series()` gains `progress_callback` and `on_error` parameters.
- **`ring_composition` on results** (ISSUE-J) — `ring_composition` dict (`{"hexagons": N, "pentagons": M}`) exposed on `CarbonSkeleton`, `BiocharGenerator`, and `BiocharResult`.
- **Static OPLS charge neutrality** (ISSUE-K) — residual charge redistributed after static assignment; debug log when residual > 0.01 e before correction.

### Changed

- **Version bump** — `0.2.0` → `0.3.0` to reflect new public API surface.

---

## [0.2.0] — June 7, 2026

### Added

- **Named `BiocharResult` dataclass** — `generate_biochar()` now returns a `BiocharResult` with named fields (`mol`, `coords`, `gro_path`, `top_path`, `itp_path`, `ring_composition`). Positional unpacking still works via `__iter__` (issue #5).
- **`write_files=False` parameter** — pass `write_files=False` to `generate_biochar()` to skip all GROMACS file I/O; `result.gro_path` / `.top_path` / `.itp_path` will be `None` (issue #5).
- **QM 1.14*CM1A charge backend** — `charge_method="qm"` runs a single-point AM1 calculation via an external MOPAC binary, maps Mulliken charges through the CM1A correction, and scales by 1.14 (LigParGen methodology). Requires `conda install -c conda-forge mopac`. See `docs/qm-charge-backend.md`.

---

## [0.1.5] — June 1, 2026

### Added

- **Temperature × feedstock composition model** — `temperature` and `feedstock` parameters on `GeneratorConfig` / `generate_biochar()` derive H/C and O/C ratios from a regression model trained on the UC Davis Biochar Database (Davis *et al.*, 2024, DOI 10.1016/j.xcrp.2024.102036).
- **CLI `--temperature` and `--feedstock` flags** — expose the composition model in `biochar-gen`.
- **`TemperatureModel` public class** — `from biochar import TemperatureModel` for direct access to the underlying model.

---

## [0.1.4] — May 31, 2026

### Added

- **ML-based partial charge refinement** — opt-in `charge_method="ml"` using a bundled Gaussian-process model trained on OPLS reference charges (issue #4).

### Changed

- CI and Read the Docs now install scikit-learn for the `ml` extra.

## [0.1.3] — May 29, 2026

### Added

- **Amorphous porous packing** — `pore_type="amorphous"` for disordered sheet packing (issue #1).
- **S-doping** — thiol and thioether functional groups (issue #3).
- **Ring-substituting nitrogen** — pyridinic / pyrrolic / graphitic (issue #2).
- Expanded test coverage (~83%).

> **Note:** the `[1.x]` entries below predate the current `0.1.x` version scheme
> and are retained for history only. The authoritative package version is
> **0.1.4** (`pyproject.toml`, `biochar/__init__.py`, git tag `v0.1.4`).

## [1.2.0] — April 16, 2026

### Added

- **Pentagon Ring Defects** — New `defect_fraction` parameter (0.0–1.0) for inserting 5-membered rings during PAH graph growth
  - Parity-aware pentagon/hexagon selection for valid Kekulé structures
  - Up to 5 retries with different sub-seeds for kekulization of non-bipartite graphs
  - Produces topologically disordered graphitic structures mimicking amorphous biochar
  - Available in `GeneratorConfig`, `generate_biochar()`, `SurfaceConfig`, and `generate_surface()`

- **Porous Surface Generation** — New `generate_surface()` convenience function for slit-pore systems
  - Parallel graphene-like sheets with user-controlled pore diameter (Ångströms)
  - Per-sheet chemistry control via `sheet_overrides` 
  - Identical sheet optimization (single `.itp` + count) and distinct sheet support
  - Automatic SVD-based flattening and z-positioning
  - Multi-residue `.gro` and `.top` export via `MultiSheetGROWriter` and `SurfaceTopologyWriter`

- **Test Coverage** — 9 new tests for pentagon defect ring growth (`TestDefectRings` in `test_generator.py`)

### Changed

- Updated `README.md` with defect_fraction and porous surface examples
- Updated `BEST_PRACTICES.md` with pentagon ring defect section
- Updated `BATCH_GENERATION_GUIDE.md` with defect_fraction parameter and examples
- Updated `RELEASE_SUMMARY.md` with current feature set and v1.2.0 date

### Internal Changes

- Added `_fuse_pentagon()` and `_fuse_hexagon()` helper functions in `carbon_skeleton.py`
- Rewrote `_grow_graph()` with parity-aware pentagon/hexagon selection
- Updated `PAHAssembler.generate()` and `._build_from_seed()` to handle defect mode
- Added `CARBON_VDW_DIAMETER` constant (3.4 Å) in `constants.py`
- New `surface_builder.py` module with `SurfaceBuilder`, `SurfaceConfig`, `SheetResult` classes
- New `MultiSheetGROWriter` and `SurfaceTopologyWriter` in `gromacs_export.py`

### Tests

All 59 tests pass (50 existing + 9 new defect ring tests)

---

## [1.1.0] — March 31, 2026

### Added

- **Valence Validation System** — Comprehensive atom valence checking
  - `ValenceValidator` class with `validate_molecule()` and `print_valence_report()`
  - `SafeBondAdder` for safe bond addition with valence constraints
  - `ValenceReport` for generating valence statistics
  - Enforcement of standard valences (H:1, C:4, O:2, N:3, S:2-6, halogens:1)
  - Integration throughout heteroatom assignment pipeline

- **VALENCE_SYSTEM.md** — Comprehensive guide to the valence validation system

### Changed

- **Fixed coordinate units**: `.gro` files now export in nanometers (nm), not Ångströms
  - Critical fix for GROMACS compatibility (1 Å = 0.1 nm conversion in `GROFileWriter`)
- Updated `README.md` with valence documentation
- Improved error messages in heteroatom assignment and validation

### Internal Changes

- New `src/valence.py` module with valence checking classes
- Integrated `ValenceValidator` into `src/validation.py`
- Used `SafeBondAdder` in `OxygenAssigner` and `HydrogenAssigner`
- Fixed coordinate unit conversion in `GROFileWriter.write()`

### Bug Fixes

- ✅ Coordinate units (Ångströms → nanometers)
- ✅ Missing valence checks during bond addition
- ✅ Generic error messages → detailed valence error reporting

---

## [1.0.0] — March 2026

### Initial Release

- **PAH Assembly** — Polycyclic aromatic hydrocarbon skeleton generation
  - PAH library with 18 validated entries (benzene → coronene)
  - Hex-lattice seed builder for >24 carbon structures
  - Ring-growth engine with fusable edge detection

- **Structure Generation Pipeline**
  - Carbon skeleton assembly (100% aromatic)
  - Oxygen functional group placement (hydroxyl, carboxyl, ether, etc.)
  - Hydrogen saturation to match H/C and O/C ratios
  - 3D coordinate generation with RDKit + MMFF94 force field minimization
  - Aromatic planarity enforcement (≤80 heavy atoms: ETKDGv3; >80: flat graphene + FF)

- **OPLS-AA Integration**
  - Automatic atom type assignment (CA, HA, CT, HC, OH, OS, OC, etc.)
  - Partial charge calculation from OPLS-AA parameters
  - Bond/angle/dihedral definitions for GROMACS

- **GROMACS Export**
  - `.gro` structure files (3D coordinates, box vectors)
  - `.top` topology files (forcefield includes, molecule definitions)
  - `.itp` molecule definition files (reusable)
  - Proper GROMACS format compliance (residue naming ≤5 chars)

- **Batch Generation**
  - `generate_biochar_series()` for multiple structures
  - Automatic combined `.top` file for mixed simulations
  - Perfect for temperature series (BC400, BC600, BC800) and composition studies

- **Validation**
  - Composition validation (H/C, O/C ratios within tolerance)
  - Chemical feasibility checks (aromaticity, charge distribution)
  - Structure quality checks (connectivity, ring planarity)
  - 3-stage validation engine with detailed error reporting

- **Documentation**
  - Comprehensive README with examples
  - BATCH_GENERATION_GUIDE with naming conventions and use cases
  - BEST_PRACTICES guide for users
  - GROMACS_WORKFLOW for integration with molecular dynamics

### Capabilities

- **Size**: 6–200+ carbons (exact for small, ~80–95% for medium, ~70–90% for large)
- **H/C ratio**: 0.3–1.0 with configurable tolerance
- **O/C ratio**: 0.0–0.5 with configurable tolerance
- **Aromaticity**: 100% aromatic (default), configurable % target
- **Functional groups**: 7 types (phenolic, hydroxyl, carboxyl, ether, carbonyl, quinone, lactone)

---

## Historical Notes

### Valence System (v1.1)

Before v1.1, atoms could have incorrect valence counts (e.g., carbon with 3 bonds instead of 4). The comprehensive valence validation system in v1.1 ensures all atoms have correct minimum and maximum bonds:

- **Carbon**: exactly 4 bonds
- **Hydrogen**: exactly 1 bond  
- **Oxygen**: exactly 2 bonds
- **Nitrogen** (future): exactly 3 bonds

### Coordinate Units (v1.1)

Before v1.1, `.gro` files were in Ångströms (Å), which GROMACS doesn't use. The fix in v1.1 converts to nanometers (nm), the standard SI unit for GROMACS simulations.

### Pentagon Ring Defects (v1.2)

The v1.2 defect system introduces pentagon (5-membered) rings during PAH growth, producing topologically disordered structures. This is motivated by:

1. **Amorphous biochar** has both hexagonal and pentagonal rings
2. **Defects introduce curvature** in otherwise flat graphitic structures
3. **Parity constraints** ensure all structures can be kekulized (valid electron distributions)

Key insight: Pentagon addition (+3 nodes, odd) fixes parity when the current graph has odd nodes, whereas hexagon (+4 nodes, even) doesn't. This is automatically handled by `_grow_graph()`.

### Surface Generation (v1.2)

Porous surface generation was added to support:

1. **Slit-pore models** — Parallel sheets mimicking activated carbon pores
2. **Controllable chemistry** — Per-sheet functional groups and composition
3. **GROMACS compatibility** — Proper multi-residue `.gro` and `.top` format

The implementation reuses the existing `BiocharGenerator` pipeline for each sheet, then:
1. Flattens sheets to xy plane via SVD
2. Positions sheets along z at `i * (pore_diameter + 3.4 Å)` spacing
3. Centers the system in a periodic box
4. Exports as a single multi-residue system

---

## Future Enhancements

- [ ] Amorphous packing (Phase 2 of porous surfaces)
- [ ] Nitrogen-doped biochar
- [ ] Sulfur-containing biochar
- [ ] Machine learning-based charge refinement
- [ ] Direct GROMACS simulation validation
- [ ] GUI/Web interface
- [ ] Pre-computed structure database

---

**Current Version**: 0.1.4  
**Last Updated**: May 31, 2026  
**Status**: Production Ready
