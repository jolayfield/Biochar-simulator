# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project knowledge

- `docs/solutions/` — documented solutions to past problems (bugs, conventions, patterns), organized by category with YAML frontmatter (`module`, `tags`, `problem_type`). Relevant when implementing or debugging in a documented area.
- `CONCEPTS.md` — shared domain vocabulary. Relevant when orienting to the codebase or discussing domain concepts.

## Commands

```bash
# Install from conda-forge (recommended)
conda install -c conda-forge biochar

# or from PyPI
pip install biochar

# Development install (editable, with extras)
pip install -e ".[dev,viz]"

# Run unit tests
python -m pytest tests/test_generator.py -v
python -m pytest tests/test_surface_builder.py -v

# Run a single test
python -m pytest tests/test_generator.py::TestCarbonSkeleton::test_pah_assembler_benzene -v

# PAH quality suite (not pytest — run directly)
python tests/test_pah_quality.py

# Build docs
cd docs && make html
```

The test discovery root is `tests/` (configured in `pyproject.toml`). The two loose test files at the project root (`test_valence_comprehensive.py`, `test_valence_update.py`) are older scratch tests, not part of the suite.

## Architecture

The package lives entirely in `biochar/`. The public API is in `biochar_generator.py`; all other modules are internal.

### Generation pipeline (single molecule)

`BiocharGenerator.generate()` runs five sequential steps, each in its own module:

1. **`carbon_skeleton.py`** — `PAHAssembler` builds a PAH graph. For targets ≤40 C it picks an exact SMILES from `PAH_LIBRARY` (in `constants.py`). For larger targets it selects the closest library seed then iteratively fuses rings (hexagons or pentagons controlled by `defect_fraction`) until the carbon count is within tolerance. A custom hex-lattice position tracker keeps all coordinates consistent so `geometry_3d` receives a flat, strain-free sheet.

2. **`heteroatom_assignment.py`** — `OxygenAssigner` places functional groups (phenolic, carboxyl, ether, etc.) either from an explicit dict or derived from `O_C_ratio`. `HydrogenAssigner` then fills remaining free valences. A critical function `_fix_heteroatom_bond_types` (also called in `biochar_generator.py` after geometry and after validation) corrects RDKit's tendency to mark ether C–O bonds as `AROMATIC` during sanitisation.

3. **`geometry_3d.py`** — `CoordinateGenerator` embeds in 3D. Molecules ≤80 heavy atoms use RDKit ETKDGv3/v2 + MMFF94; larger molecules use the pre-computed hex-lattice 2D coords. When `generator.used_hex_lattice` is `True`, clash resolution is **skipped** because peri-H contacts in large fused PAHs are real geometry, not errors — displacing atoms would shatter the ring lattice.

4. **`opls_typing.py`** — `AtomTyper` maps each atom to an internal OPLS-AA type (CA, HA, OH, OS, …). At export time `GROMACS_OPLS_TYPE_MAP` in `constants.py` translates these to GROMACS `opls_XXX` names.

5. **`validation.py`** — `ValidationEngine` checks composition ratios and geometry. Steric clash warnings from the flat hex-lattice path are expected and do not indicate a bug.

### Surface pipeline

`generate_surface()` → `SurfaceBuilder` (`surface_builder.py`). Generates each sheet through the single-molecule pipeline, flattens to the xy plane (SVD best-fit rotation), then stacks along z with spacing `pore_diameter + 3.4 Å` (graphene interlayer distance). Identical sheets are generated once and deep-copied. The combined `.gro` has all sheets as separate residues; `.top` references one `.itp` (identical sheets) or one `.itp` per unique sheet type.

### GROMACS export

`gromacs_export.py` contains three writers: `GROFileWriter`, `ITPFileWriter`, `TOPFileWriter`, orchestrated by `GromacsExporter`. Coordinates are converted from RDKit Å → GROMACS nm (× 0.1). Residue names are hard-limited to 5 characters by GROMACS `.gro` format.

### Key constraints

- `molecule_name` / residue name: **≤5 characters** (GROMACS hard limit).
- Carbonyl, quinone, and lactone functional groups fall back silently to phenolic — pure aromatic PAH edges have no free valence for these.
- `_fix_heteroatom_bond_types` must be called after any RDKit `SanitizeMol` pass that touches a molecule containing ether oxygens.
- The hex-lattice path (large molecules) produces geometrically perfect structures; clash warnings from `GeometryValidator` on these are artefacts and should be ignored — GROMACS energy minimisation resolves them.
