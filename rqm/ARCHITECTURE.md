# Architecture <!-- rq-7efb46e4 -->

System design overview for `biochar`. This is the entry point for the requirements set in
`rqm/`; individual feature documents specify behaviour and carry the traceability IDs that link
requirements to the tests defending them.

Domain vocabulary lives in `CONCEPTS.md` at the repository root and is not repeated here. The
distinction between an **internal atom type**, an **OPLS type name**, and a **bonded type** is
load-bearing throughout the export path — read it there before working on typing or export.

## What this project is <!-- rq-bb785e6d -->

A Python library that generates biochar polycyclic aromatic hydrocarbon structures with
realistic composition and exports them as GROMACS topologies parameterised against OPLS-AA.
Inputs are composition targets (carbon count, O/C and H/C ratios, aromaticity, functional group
distribution, pH) or a pyrolysis temperature and feedstock from which those targets are derived.
Outputs are `.gro`, `.itp`, and `.top` files ready for molecular dynamics.

## Public surface <!-- rq-12b97d18 -->

`pyproject.toml` is the source of truth for what is public — `[project.scripts]` declares the
console entry points and `src/biochar/__init__.py`'s `__all__` declares the Python API. Do not infer
the public surface from this document.

| Console script | Module | Purpose |
|---|---|---|
| `biochar-gen` | `cli.py` | Generate a single structure or series |
| `biochar-sweep` | `sweep_cli.py` | Run a declarative factorial parameter sweep |
| `biochar-md-setup` | `md_setup_cli.py` | Generate GROMACS run directories from a sweep manifest |
| `biochar-condense` | `condensation_cli.py` | Set up a Wood et al. 2024 condensation-annealing run |

The four CLI modules parse arguments only; logic lives in the module each one wraps.

## Layering <!-- rq-18307b81 -->

The package is a single flat module directory under `src/biochar/`, but its imports form an
almost-clean DAG. The layers below are the real structure and the intended seams for any future
decomposition.

The src layout is deliberate: because there is no importable `biochar/` at the repository root, the
package must be installed to be imported. A test run therefore cannot pass by exercising the working
tree while the installed package is broken — the failure mode that lets a packaging bug reach users
with CI green.

- **Layer 0** — no intra-package imports: `constants`, `valence`, `qm_charges`
- **Layer 1** — depend on `constants`: `carbon_skeleton`, `geometry_3d`, `opls_typing`,
  `temperature_model`
- **Layer 2** — `heteroatom_assignment`, `gromacs_export`, `ml_charges`, `protonation`,
  `validation`
- **Layer 3** — `biochar_generator`, which orchestrates the single-molecule pipeline
- **Layer 4** — `surface_builder`, `condensation`, `sweep`, `md_setup`
- **Layer 5** — the four CLI modules

One layering violation remains: `surface_builder` imports `biochar_generator` at module scope
while `biochar_generator` imports `SurfaceBuilder` inside a function to break the resulting
cycle. `generate_surface` is a surface-level entry point living in the single-molecule module.

## Single-molecule pipeline <!-- rq-60a0144a -->

`BiocharGenerator.generate()` runs five sequential stages, each in its own module.

1. **Carbon skeleton** (`carbon_skeleton.py`) — `PAHAssembler` builds the PAH graph. Targets of
   40 carbons or fewer take an exact SMILES from `PAH_LIBRARY`; larger targets select the
   closest seed and fuse rings until the count is within tolerance, with pentagon frequency
   controlled by `defect_fraction`. A hex-lattice position tracker keeps coordinates consistent
   so the geometry stage receives a flat, strain-free sheet.

2. **Heteroatom assignment** (`heteroatom_assignment.py`) — `OxygenAssigner` places functional
   groups; `HydrogenAssigner` fills the remaining valences. Specified in
   `heteroatom-assignment.md`.

3. **Geometry** (`geometry_3d.py`) — `CoordinateGenerator` embeds in 3D; `GeometryValidator`
   judges the result. Specified in `geometry-embedding.md`.

4. **Typing** (`opls_typing.py`) — `AtomTyper` assigns internal types; `GROMACS_OPLS_TYPE_MAP`
   translates them at export. Specified in `opls-typing.md`.

5. **Validation** (`validation.py`) — `ValidationEngine` checks composition ratios and geometry.

## Surface pipeline <!-- rq-e13b2348 -->

`generate_surface()` delegates to `SurfaceBuilder`. Each sheet goes through the single-molecule
pipeline, is flattened to the xy plane by an SVD best-fit rotation, then stacked along z with
spacing `pore_diameter + 3.4 Å` — the graphene interlayer distance. Identical sheets are
generated once and deep-copied. The combined `.gro` carries each sheet as a separate residue;
the `.top` references one `.itp` for identical sheets, or one per unique sheet type.

## Export <!-- rq-92c69795 -->

`gromacs_export.py` holds `GROFileWriter`, `ITPFileWriter`, and `TOPFileWriter`, orchestrated by
`GromacsExporter`. Coordinates convert from RDKit ångström to GROMACS nanometres. Residue names
are hard-limited to 5 characters by the `.gro` format, which constrains `molecule_name`.

## Standing constraints <!-- rq-faf5d18d -->

These hold across the whole system and are easy to violate locally.

- `molecule_name` and residue name must be **5 characters or fewer** — a GROMACS format limit.
- `_fix_heteroatom_bond_types` must be called after **any** `SanitizeMol` pass touching a
  molecule with ether oxygens.
- Clash warnings from `GeometryValidator` on hex-lattice structures are artefacts of a
  geometrically exact lattice, not defects.
- Requested composition is a target, not a promise. Measure the realised structure.
- Only the `opls_XXX` name reaches GROMACS. A name identifying the wrong element silently
  simulates the wrong chemistry.

## Requirements set <!-- rq-4f998061 -->

| Document | Covers |
|---|---|
| `geometry-embedding.md` | Embedding path selection, clash detection, bond-length validation |
| `heteroatom-assignment.md` | Functional group placement, fallbacks, oxygen spill, hydrogen saturation |
| `opls-typing.md` | Internal typing, GROMACS mapping, forcefield verification depths |

Not yet specified, in rough priority order: carbon skeleton growth, valence validation, GROMACS
export, parameter sweep, pH protonation, temperature/feedstock model, surface stacking, MD
setup, condensation annealing, QM and ML charge backends.

An absent document is not an absent requirement — it means the behaviour is currently defended
by tests that no specification points at. `tools/rqm/rq check` reports the converse: specified
requirements that no test defends.
