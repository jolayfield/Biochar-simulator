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

`tests/` holds 16 files; the three named above are examples, not the whole suite. A `PreToolUse` hook in
`.claude/settings.json` runs the full suite before any `git commit` and blocks the commit if it fails.

## Architecture

The package lives entirely in `biochar/`.

**`pyproject.toml` is the source of truth for what is public.** `[project.scripts]` declares the four console
entry points; `biochar/__init__.py`'s `__all__` declares the public Python API. Do not infer the public surface
from this file — check those two.

| Console script | Module | Purpose |
|---|---|---|
| `biochar-gen` | `cli.py` | Generate a single structure or series |
| `biochar-sweep` | `sweep_cli.py` | Run a declarative factorial parameter sweep |
| `biochar-md-setup` | `md_setup_cli.py` | Generate GROMACS run directories from a sweep manifest |
| `biochar-condense` | `condensation_cli.py` | Set up a Wood et al. 2024 condensation-annealing run |

### Module map

`biochar_generator.py` orchestrates the single-molecule pipeline (five steps, below). Beyond that pipeline:

- **`sweep.py`** — declarative parameter-sweep driver; `run_sweep`, `expand_grid`, `build_point`, and the
  sweep manifest that `md_setup` consumes.
- **`md_setup.py`** — writes GROMACS run inputs (`.mdp` templates plus a driver script) per structure;
  `setup_md_from_manifest`, ion profiles, pre-solvation stages.
- **`condensation.py`** — Wood et al. 2024 condensation-annealing setup (parallel construction mode).
- **`temperature_model.py`** — data-driven temperature × feedstock property model; `TemperatureModel`,
  `properties`, `VALID_FEEDSTOCKS`.
- **`qm_charges.py`** — LigParGen-style QM partial charges (1.14*CM1A) via an external MOPAC binary.
  Requires a MOPAC install; raises `QMChargeError` when unavailable.
- **`ml_charges.py`** — ML-based partial charge refinement.
- **`valence.py`** — valence validation system (see `VALENCE_SYSTEM.md`).
- **`cli.py`**, **`sweep_cli.py`**, **`md_setup_cli.py`**, **`condensation_cli.py`** — argument parsing only;
  the logic lives in the module each one wraps.

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
