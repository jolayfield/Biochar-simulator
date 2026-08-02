# Biochar Simulator — Changelog

All significant changes to the Biochar Simulator project are documented here.

## [Unreleased]

### Changed


- **`TemperatureModel.get_valid_range` now takes the property it is being asked
  about.** *Breaking, for direct callers of the model.* The signature is
  `get_valid_range(prop, feedstock=None)`; `prop` is required and comes first.

  It previously ignored its argument in effect — it looped over H/C and O/C and
  returned the first range it found — so pH (fitted 200–900 °C) and electrical
  conductivity (220–900 °C) were both reported as spanning 100–1000 °C. Code
  that called `get_valid_range("softwood")` will now raise `KeyError` for an
  unknown property rather than silently returning the wrong range. It is not
  part of the public API declared in `biochar.__all__`.

- **A prediction that falls back from a feedstock curve to the pooled curve says
  so.** A feedstock's own curve covers only the temperature range its data
  spans; outside that the pooled curve answers, which was previously
  indistinguishable from a feedstock-backed answer. Reported at `INFO`, so a
  feedstock sweep can tell the two apart.

- **An aromaticity derived from an H/C outside the fitted range warns.**
  `aromaticity_from_hc` clamps into 0–100%, which kept the output physical and
  hid how far outside the fit the input was — an H/C of 5.0 against a relation
  fitted over 0.03–1.64 returned exactly 0.0%, reading as a confident prediction
  rather than a refusal. The clamp is unchanged; it is now accompanied by a
  `UserWarning`. Reachable through ordinary configuration, since H/C up to 2.0
  is allowed.


- **Strict mode now fails on a shortfall, not only on total failure.** A request
  for explicit `functional_groups` that the skeleton cannot host was previously
  accepted as long as at least one group of each kind was placed — a request for
  40 ether bridges answered with 6 returned successfully. `strict=True` now
  raises `ValidationError` naming each group as `placed/requested`.

  This is scoped to groups the caller named explicitly. Counts derived from
  `O_C_ratio` are unaffected: there the ratio is the target and the count is an
  implementation detail of reaching it, judged by `O_C_tolerance` as before.

  **Sweeps will report more `fallback` and fewer `strict_pass` rows.** A point
  whose group request was only partly met now exhausts its seed retries and
  falls back, which is what those statuses were always meant to mean. Nothing
  crashes and no structure is lost — only the label changes, and it changes to
  the honest one. Manifest counts are not comparable with earlier runs.

- **Sheets copied by the identical-sheet optimisation no longer share a
  composition record.** The copy duplicated the molecule, coordinates, atom
  types and charges but passed `composition` through by reference, so annotating
  one sheet's composition changed every other sheet on the surface. Latent
  rather than observed — nothing in the package mutated it — but it made the
  copies not copies.

- **A surface that fell back from amorphous to slit geometry is no longer
  labelled amorphous.** `amorphous_fallback="slit"` degrades gracefully when
  packing cannot converge, but `pore_type` was never updated and the `.gro`
  title was chosen from it — so a slit stack was written to disk titled
  `"Amorphous surface"`. The title now states the geometry actually built and
  notes the substitution.

- **`generate_surface` exposes the options that decide the geometry it returns**
  — `amorphous_fallback`, `aromaticity_percent`, `box_padding_xy` and
  `box_padding_z`. `amorphous_fallback` matters most: without it, an amorphous
  request that cannot be packed raised instead of degrading, and the graceful
  path was unreachable from the entry point the documentation points at.
  Defaults preserve the previous behaviour.

### Added


- **`TemperatureModel.predict_with_evidence(temperature, prop, feedstock=None)`**
  returns a prediction together with what it rests on: the observation count and
  spread at the nearest grid point, whether that point had no observations and
  the value was carried in from a neighbour, whether the temperature is outside
  the curve's support, and which curve answered. The model already recorded all
  of it; nothing surfaced it. Electrical conductivity, for instance, has zero
  observations at 100 °C and predicts 1.0 there.

## [0.5.0] — July 31, 2026

> **Exported topologies changed. Re-run anything whose results matter.**
> Structures exported before this release were missing every 1–4 non-bonded
> interaction and had no term holding aromatic rings planar. Trajectories from
> those topologies are not comparable with trajectories from new ones. The
> structures themselves are unaffected — only the force field applied to them.

### Fixed — force field

- **1–4 interactions are no longer silently dropped.** `[ moleculetype ]`
  declares `nrexcl = 3`, which correctly excludes 1–4 pairs from the plain
  non-bonded loop, but OPLS-AA restores them at half strength via `gen-pairs` /
  `fudgeLJ` / `fudgeQQ` — and `gen-pairs` supplies the pair *parameters*, never
  the pair *list*. No `[ pairs ]` section was written, so every 1–4 interaction
  was excluded and none restored. Exported topologies now carry the full 1–4
  pair list.
- **Aromatic rings now carry improper torsions.** Proper dihedrals restrain
  rotation about a bond and do nothing to stop a substituted ring carbon
  pyramidalising, so sheets were planar only because they were written that way.
  Impropers are emitted for every three-connected aromatic ring carbon as
  `improper_Z_CA_X_Y` (function type 1, per `oplsaa.ff`'s own residue templates).
- **Atom names above five characters no longer collide or diverge.** The `.gro`
  atom-name field is five characters; truncating made `C10000` and `C1000` the
  same name inside one file, and disagree with the untruncated `.itp`. Indices
  that do not fit are now base-36, and all four writers share one helper.
  Affects systems above ~10,000 atoms of one element — stacked surfaces mainly.

### Fixed — MD run setup

- **The local run pipeline could not pass solvation.** `gmx solvate -p` updates
  a topology in place, and the generated script passed it a `wet.top` that
  nothing created. Under `set -euo pipefail` the run aborted there. (The SLURM
  path was unaffected.)
- **Multi-species ion placement discarded all but the last species.** One
  `genion.tpr` was built from pre-ion coordinates and reused for every cation,
  so each was placed into a structure containing none of the previous ones while
  the topology accumulated all of them. The run input is now rebuilt between
  species, in both the local and SLURM paths.
- **Annealing follows the pyrolysis temperature.** The schedule applied the
  400 °C Wood protocol (1000 K peak, 1 fs) to every structure regardless of
  temperature; an 800 °C char now anneals to 3000 K at 0.5 fs. Rendered from
  `workflows.condensation.anneal_spec_for_htt`, so the protocol has one
  implementation rather than two that had drifted apart. When a manifest records
  no temperature the mildest anchor applies, and the run directory says so.
- **Include files are resolved from the sweep manifest** rather than guessed as
  `<topology stem>.itp`. The guess held for a single molecule and failed for a
  surface, whose topology is `<base>.top` while its include is
  `<base>_sheet.itp` — so the include was silently not copied and `grompp` could
  not find the moleculetype.
- **`grompp` warnings are no longer blanket-suppressed.** Every call passed
  `-maxwarn 2`, absorbing whichever two warnings arrived — including the two
  that matter most here, a net-charged system under PME and a structure/topology
  atom mismatch. `-maxwarn 1` is now spent only where the topology declares a
  non-zero charge and only before `genion -neutral`; a neutral structure
  suppresses nothing.

### Changed

- **Validation reports every geometry error, not the first three.** This does
  not change any pass/fail decision — `is_valid` was already `len(errors) == 0`
  — but `n_validation_errors` in sweep manifests is now the true count rather
  than one capped at 3, so numbers will not match manifests generated earlier.
  Console output is unchanged: `print_summary` still shows three and states the
  total.
- **Run directories carry `run_provenance.json`** recording label, status,
  source files, net charge, and the annealing schedule used. A `fallback`
  structure is no longer indistinguishable from a `strict_pass` one at the point
  results are collected.

### Added — project

- `rqm/gromacs-export.md` and `rqm/md-setup.md`, plus a `docs/reviews/`
  category for dated, commit-pinned assessments. 64 requirement scenarios, each
  with a defending test.
- CI now fails on any specified scenario with no defending test, and on a stale
  traceability reference.
- The `pull_request` CI trigger no longer filters on `branches: [main]`, which
  had let a PR based on any other branch run zero checks while reporting
  `MERGEABLE / CLEAN`.

## [0.4.0] — July 11, 2026

### Added

- **Declarative parameter-sweep driver** (`biochar.sweep`) — describe a factorial
  grid (e.g. pyrolysis temperature × feedstock) in YAML/JSON; the driver expands
  the cartesian product, generates every structure strict-first with seed-retry →
  fallback, and writes a manifest (CSV + JSON). New `biochar-sweep` CLI and
  `examples/sweeps/temperature_grid.yaml` / `oxygen_group_grid.yaml`.
- **GROMACS MD run-setup pipeline** (`biochar.md_setup`) — turns each sweep
  manifest row into a ready-to-submit run directory (dry-anneal → solvate → ion →
  wet-equilibrate), writing files only (no `gmx` invoked). Includes
  Minnesota-groundwater `IonProfile`/`ION_PROFILES` presets and a
  `biochar-md-setup` CLI.
- **Generic pre-solvation insertion seam** — `PreSolvationStage` +
  `MoleculeInsertion` on `MDSetupConfig.pre_solvation_stage` let a downstream
  workflow inject an `insert-molecules` stage (after dry-anneal, before
  solvation) and switch the topology used from solvation onward, without
  `md_setup` knowing what the molecules are. (The PFAS-sorption workflow now
  lives in the separate `biochar-pfas` project and consumes this seam.)
- **`allow_aliphatic` on `GeneratorConfig`** — allow pendant sp3 (aliphatic)
  carbon so H/C targets above the pure-aromatic ceiling are reachable.

### Fixed

- **H/C ratio is now actually reached** — previously every structure landed well
  below the requested H/C (~0.14 low at 100 C) because H/C was bounded by the
  aromatic perimeter and the target never fed back into the skeleton. Now: honest
  ceiling reporting (`CompositionResult.h_c_ceiling` / `.h_c_target_unreachable`),
  H/C-aware skeleton growth (elongation), and aliphatic decoration reach targets
  0.4–0.8 within tolerance for 30/50/100 C, and strict composition validation
  passes where it previously could not.
- **`test_ml_charges` skips instead of erroring** when the optional `scikit-learn`
  extra is absent (`pytest.importorskip`).

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

**Released version**: 0.5.0 (see `pyproject.toml`, which is authoritative)  
**Status**: Production Ready
