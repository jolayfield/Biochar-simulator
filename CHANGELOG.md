# Biochar Simulator — Changelog

All significant changes to the Biochar Simulator project are documented here.

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

**Current Version**: 1.2.0  
**Last Updated**: April 16, 2026  
**Status**: Production Ready
