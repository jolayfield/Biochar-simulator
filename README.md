# Biochar Simulator — Structure Generator

A Python package for generating realistic biochar molecular structures for GROMACS molecular dynamics simulations. Supports single molecules, temperature/composition series, and **porous slit-pore surfaces**.

---

## Overview

The Biochar Simulator builds polycyclic aromatic hydrocarbon (PAH) structures with user-specified compositional parameters and exports GROMACS-ready force field files. Supports both pure hexagonal aromatic skeletons and topologically disordered structures with pentagon ring defects.

**Key capabilities:**

- **Size**: 6 to 200+ carbons (exact count for most targets)
- **Composition**: H/C and O/C ratios with configurable tolerances
- **Aromaticity**: 100% aromatic skeleton (graphene-nanoflake topology, or defective with pentagons)
- **Ring Defects**: Optional pentagon insertion during graph growth (`defect_fraction` parameter)
- **Functional Groups**: Exact counts via dict API — phenolic, carboxyl, ether, carbonyl, quinone, lactone, hydroxyl
- **Porous Surfaces**: Slit-pore systems of stacked PAH sheets with user-controlled pore diameter
- **Output**: GROMACS `.gro`, `.top`, `.itp` files with OPLS-AA force field

---

## Installation

### conda (recommended)

```bash
conda install -c conda-forge biochar
```

### PyPI

```bash
pip install biochar
```

### Requirements

- Python 3.9+
- RDKit ≥ 2023.9
- NumPy ≥ 1.24, SciPy ≥ 1.10, NetworkX ≥ 3.1

---

## Quick Start

### Single molecule

```python
from biochar.biochar_generator import generate_biochar

mol, coords, gro_path, top_path, itp_path = generate_biochar(
    target_num_carbons=100,
    H_C_ratio=0.5,
    O_C_ratio=0.1,
    output_directory="output",
    basename="biochar_100C",
    seed=42,
)
```

### With specific functional groups

Use a dict to place **exact counts** of each group type:

```python
mol, coords, gro, top, itp = generate_biochar(
    target_num_carbons=50,
    functional_groups={"phenolic": 3, "carboxyl": 1, "ether": 2},
    output_directory="output",
    basename="biochar_fg",
    seed=42,
)
```

When `functional_groups` is `None` (default), total oxygen is controlled by `O_C_ratio` and placed as phenolic groups.

### With pentagon ring defects

Add topological disorder by inserting 5-membered rings during growth:

```python
# ~15% of rings will be pentagons instead of hexagons
mol, coords, gro, top, itp = generate_biochar(
    target_num_carbons=60,
    defect_fraction=0.15,  # probability per ring addition
    output_directory="output",
    basename="biochar_defects",
    seed=42,
)
```

`defect_fraction` ranges from 0.0 (pure hexagonal PAH) to 1.0 (all pentagons). Typical values: 0.1–0.2.

#### Understanding `defect_fraction`

`defect_fraction` is a *request probability*, not a guaranteed pentagon fraction. During skeleton growth each ring-addition event has a `defect_fraction` chance of *attempting* a pentagon; parity and adjacency constraints reject many of these attempts. Observed pentagon counts are lower than the requested fraction:

| `defect_fraction` | ~50 C skeleton | ~100 C skeleton |
|:-----------------:|:--------------:|:---------------:|
| 0.05              | 0–1 pentagons  | 1–3 pentagons   |
| 0.15              | 1–3 pentagons  | 2–6 pentagons   |
| 0.30              | 2–5 pentagons  | 4–10 pentagons  |

Use `result.ring_composition` to inspect the actual counts after generation:

```python
result = generate_biochar(target_num_carbons=60, defect_fraction=0.15, seed=42)
print(result.ring_composition)  # e.g. {"hexagons": 12, "pentagons": 2}
```

### Slit-pore surface

```python
from biochar.biochar_generator import generate_surface

# Two identical sheets, 10 Å pore
sheets, gro, top, itps = generate_surface(
    target_num_carbons=50,
    functional_groups={"phenolic": 2, "ether": 1},
    pore_diameter=10.0,
    output_directory="output",
    basename="slit_pore",
    seed=42,
)

# Asymmetric pore — different chemistry on each wall
sheets, gro, top, itps = generate_surface(
    pore_diameter=8.0,
    sheet_overrides=[
        {"functional_groups": {"phenolic": 3}, "target_num_carbons": 40},
        {"functional_groups": {"carboxyl": 2}, "target_num_carbons": 50},
    ],
    output_directory="output",
    basename="asymmetric_pore",
)
```

### Batch generation (temperature/composition series)

```python
from biochar.biochar_generator import generate_biochar_series

configs = [
    {"molecule_name": "BC400", "target_num_carbons": 80,  "H_C_ratio": 0.65, "O_C_ratio": 0.20, "seed": 1},
    {"molecule_name": "BC600", "target_num_carbons": 100, "H_C_ratio": 0.55, "O_C_ratio": 0.12, "seed": 2},
    {"molecule_name": "BC800", "target_num_carbons": 120, "H_C_ratio": 0.40, "O_C_ratio": 0.05, "seed": 3},
]

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/temperature_series",
    create_combined_top=True,   # writes combined.top
)
```

Then run in GROMACS:

```bash
cd output/temperature_series
gmx grompp -f md.mdp -p combined.top -o topol.tpr
gmx mdrun -deffnm topol
```

---

## Parameter sweeps (declarative grids)

`generate_biochar_series` takes a hand-written list of configs. For a full
**factorial grid** — every combination of pyrolysis temperature × feedstock ×
oxygen content — use the sweep driver instead. You describe the axes once in a
YAML/JSON file and the pipeline expands the grid, generates every structure,
and writes a manifest table for downstream analysis.

### CLI

```bash
# Emit a starter config you can edit
biochar-sweep template > my_sweep.yaml

# Run a sweep
biochar-sweep run examples/sweeps/temperature_grid.yaml
```

### Config format

```yaml
name: temperature_grid
output_directory: sweep_out/temperature_grid
seed: 1000              # base RNG seed; each grid point gets a distinct offset
max_retries: 8          # strict-validation seed retries before fallback
on_validation_fail: fallback   # fallback | skip | strict
name_template: "T{temperature}_{feedstock}"

axes:                   # cartesian product of every list below
  temperature: [300, 400, 500, 600, 700]
  feedstock:   [softwood, hardwood]

fixed:                  # applied to every grid point
  target_num_carbons: 100
```

Temperature and feedstock are mapped to `H_C_ratio` / `O_C_ratio` targets via
`biochar.temperature_model`, so a single temperature axis drives the oxygen
chemistry automatically. The example
[`examples/sweeps/oxygen_group_grid.yaml`](examples/sweeps/oxygen_group_grid.yaml)
sweeps explicit `functional_groups` counts instead.

### Python API

```python
from biochar.sweep import load_sweep_config, run_sweep

cfg = load_sweep_config("examples/sweeps/temperature_grid.yaml")
summary = run_sweep(cfg, quiet=True)

print(summary["n_built"], "structures")   # -> 10
import pandas as pd
df = pd.read_csv(summary["manifest_csv"])
```

### Output

```
<output_directory>/
├── manifest.csv          # one row per grid point (see below)
├── manifest.json         # same data + run metadata
└── structures/
    ├── 000_T300_softwood/   # .gro / .top / .itp for this point
    ├── 001_T300_hardwood/
    └── ...
```

The **manifest** carries one row per structure with the axis values
(`axis_temperature`, `axis_feedstock`, …), the achieved composition
(`molecular_formula`, `H_C_ratio`, `O_C_ratio`, `functional_groups`), the
build `status`, `seed_used`, `n_attempts`, validation counts, and the paths to
the three GROMACS files — ready to join against downstream analysis.

### Validation and the `status` column

Each grid point is generated **strict-first**: the driver retries with
successive seeds (`max_retries`) seeking a structure that passes full
composition + geometry validation. If none passes, `on_validation_fail`
decides what happens:

| `status` | Meaning |
|---|---|
| `strict_pass` | Passed strict validation on one of the retry seeds |
| `fallback` | Strict never passed; built with `strict=False` (files still written) |
| `skipped` | Strict never passed and `on_validation_fail: skip` |
| `failed` | A non-validation error occurred (captured in the `error` column) |

The seed-retry → fallback design exists so a sweep completes deterministically
rather than stalling when a particular grid point cannot pass strict validation
(e.g. minor flat-sheet steric clashes that resolve under GROMACS energy
minimisation) — the `.gro/.top/.itp` files are still well-formed.

---

## Configuration Parameters

### Single molecule (`GeneratorConfig` / `generate_biochar`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `target_num_carbons` | int | 50 | Target carbon count |
| `H_C_ratio` | float | 0.5 | Target H/C molar ratio |
| `H_C_tolerance` | float | 0.10 | Allowed H/C error (fraction) |
| `O_C_ratio` | float | 0.1 | Target O/C molar ratio |
| `O_C_tolerance` | float | 0.10 | Allowed O/C error (fraction) |
| `aromaticity_percent` | float | 90.0 | Target % aromatic carbons |
| `functional_groups` | dict\|None | `None` | Exact group counts, e.g. `{"phenolic": 2}` |
| `defect_fraction` | float | 0.0 | Probability [0, 1) each ring is a pentagon |
| `charge_method` | str | `"opls"` | Partial charge source: `"opls"`, `"ml"`, or `"qm"` |
| `molecule_name` | str | `"BC"` | Residue name (≤5 chars for GROMACS) |
| `periodic_box` | bool | False | Add periodic box vectors to `.gro` |
| `seed` | int\|None | None | RNG seed for reproducibility |

### Slit-pore surface (`SurfaceConfig` / `generate_surface`)

| Parameter | Type | Default | Description |
|---|---|---|---|
| `target_num_carbons` | int | 50 | Carbons per sheet |
| `H_C_ratio` | float | 0.3 | Target H/C per sheet |
| `O_C_ratio` | float | 0.05 | Target O/C per sheet |
| `functional_groups` | dict\|None | `None` | Groups applied to all sheets |
| `defect_fraction` | float | 0.0 | Pentagon probability per ring per sheet |
| `pore_diameter` | float | 10.0 | Gap between sheet surfaces (Å) |
| `num_sheets` | int | 2 | Number of parallel sheets |
| `sheet_overrides` | list\|None | `None` | Per-sheet config dicts (length = `num_sheets`) |
| `box_padding_xy` | float | 1.0 | Box padding in x/y (nm) |
| `box_padding_z` | float | 1.0 | Box padding in z (nm) |
| `system_name` | str | `"SLIT"` | Name in `.top [ system ]` section |
| `sheet_base_name` | str | `"SHT"` | Residue name base (≤3 chars) |
| `seed` | int\|None | None | RNG seed |

### Available functional groups

| Group | O added | Notes |
|---|---|---|
| `"phenolic"` | 1 | Ar–OH; always works |
| `"hydroxyl"` | 1 | Same as phenolic for pure PAH |
| `"carboxyl"` | 2 | Ar–C(=O)(OH); adds extra C |
| `"ether"` | 1 | Ar–O–Ar bridge across two edge sites |
| `"carbonyl"` | 1 | Falls back to phenolic with warning |
| `"quinone"` | 2 | Falls back to phenolic with warning |
| `"lactone"` | 2 | Falls back to phenolic with warning |

Carbonyl, quinone, and lactone require ≥2 free valence on one carbon, which is unavailable on pure aromatic PAH edge sites — they warn and substitute phenolic automatically.

### Partial charge methods (`charge_method`)

| `charge_method` | Description | Requirement |
|---|---|---|
| `"opls"` (default) | Static OPLS-AA lookup table | None |
| `"ml"` | GP model trained on OPLS reference charges | `pip install biochar[ml]` (scikit-learn) |
| `"qm"` | 1.14×CM1A via MOPAC AM1 semiempirical | `conda install -c conda-forge mopac` |

The `"qm"` backend reproduces the LigParGen 1.14*CM1A methodology using an external MOPAC binary. It runs a single-point AM1 calculation on the generated 3D geometry and maps the resulting Mulliken charges through the CM1A correction. See `docs/qm-charge-backend.md` for full details. Raises `biochar.QMChargeError` if MOPAC is not on PATH or exits with an error.

---

## Generation Pipeline

### Single molecule

```
target_num_carbons
      │
      ▼
Carbon skeleton (PAH seed + ring-growth)
      │  Parity-aware hex-lattice builder; 100% aromatic for any size
      ▼
Oxygen assignment (functional groups dict or O/C-ratio-driven)
      │
      ▼
Hydrogen assignment (fill valences → target H/C ratio)
      │
      ▼
3D geometry
      │  ≤80 heavy atoms: ETKDGv3 / ETKDGv2 embedding + MMFF94
      │  >80 heavy atoms: 2D-first embedding (flat graphene sheet) + FF minimisation
      ▼
OPLS-AA atom typing & partial charges
      │
      ▼
Validation (composition, geometry, steric clashes)
      │
      ▼
GROMACS export (.gro / .top / .itp)
```

### Slit-pore surface

```
SurfaceConfig
      │
      ▼
Generate N sheets (each via single-molecule pipeline above)
      │  Identical sheets: generate once, deep-copy remainder
      │  Distinct sheets:  generate each independently
      ▼
Flatten each sheet to xy plane (SVD best-fit plane rotation)
      │
      ▼
Stack along z: sheet_i centroid at z = i × (pore_diameter + 3.4 Å)
      │
      ▼
Compute periodic box (bounding box + padding)
      │
      ▼
Centre system in box
      │
      ▼
GROMACS export
      │  .gro  — all N sheets as separate residues, single file
      │  .itp  — one file (identical sheets) or one per sheet (distinct)
      └─ .top  — includes forcefield + itp(s); [ molecules ] count = N
```

---

## Supported Sizes

| Range | Strategy | Aromaticity | Count accuracy |
|---|---|---|---|
| 6–40 C | Exact PAH library match | 100% | Exact |
| 41–200+ C | Library seed + 4-node ring growth | 100% | ≤5% error |

### PAH library (18 validated entries)

| Molecule | Carbons | Type |
|---|---|---|
| benzene | 6 | Classic |
| naphthalene | 10 | Classic |
| anthracene / phenanthrene | 14 | Linear / angular |
| pyrene | 16 | Pericondensed |
| chrysene / tetracene / triphenylene | 18 | Various |
| pentacene / picene / hex_lattice_22 | 22 | Various |
| coronene | 24 | 7-ring pericondensed |
| hexacene / dibenzo_bc_ef_coronene | 26 | Various |
| hex_lattice_28 | 28 | Compact nanoflake |
| hex_lattice_30 | 30 | Compact nanoflake |
| hex_lattice_38 | 38 | Compact nanoflake |
| hex_lattice_40 | 40 | Compact nanoflake |

---

## Typical H/C ratios by pyrolysis temperature

| Temperature | H/C | O/C | Notes |
|---|---|---|---|
| 300–400 °C | 0.6–0.8 | 0.15–0.25 | Partially carbonised, high O |
| 500–600 °C | 0.4–0.6 | 0.08–0.15 | Moderately graphitic |
| 700–800 °C | 0.2–0.4 | 0.02–0.08 | Highly graphitic, low O |

---

## Output Files

| File | Format | Contents |
|---|---|---|
| `.gro` | GROMACS structure | Atom positions in **nm**, box vectors |
| `.top` | GROMACS topology | Force field include, molecule definitions |
| `.itp` | Include topology | Atoms, bonds, angles, dihedrals for one molecule type |

Coordinates are in **nanometers** (GROMACS convention; RDKit Å × 0.1).

For surfaces, a single `.gro` contains all sheets as separate residues, and the `.top` references one `.itp` with a molecule count (identical sheets) or one `.itp` per unique sheet type.

---

## Source Modules

| Module | Responsibility |
|---|---|
| `carbon_skeleton.py` | PAH library lookup and ring-growth engine |
| `heteroatom_assignment.py` | Oxygen (functional groups) and hydrogen placement |
| `geometry_3d.py` | 3D coordinate generation and clash resolution |
| `opls_typing.py` | OPLS-AA atom types and partial charges |
| `gromacs_export.py` | `.gro` / `.top` / `.itp` writers (single and multi-sheet) |
| `surface_builder.py` | Slit-pore surface assembly (`SurfaceBuilder`, `SurfaceConfig`) |
| `validation.py` | Composition, chemistry, and geometry checks |
| `constants.py` | OPLS-AA parameters, PAH library, VdW radii |
| `biochar_generator.py` | Public API: `generate_biochar`, `generate_surface`, `generate_biochar_series` |

---

## Testing

```bash
# Unit tests (constants, skeleton, heteroatoms, geometry, OPLS, validation, generator)
python3 -m pytest tests/test_generator.py -v

# Surface builder tests (config, geometry, GROMACS export, convenience function)
python3 -m pytest tests/test_surface_builder.py -v

# PAH quality suite (sizes 6–200C, compositions, seeds)
python3 tests/test_pah_quality.py
```

The PAH quality suite reports:
- PAH library SMILES validity and kekulizability
- Atom count accuracy (target vs actual)
- Bond errors and steric clashes per size
- Ring planarity (Å deviation from best-fit plane)
- H/C and O/C ratio accuracy across compositions

---

## Known Limitations

| Issue | Notes |
|---|---|
| H/C ratio loose for very small structures (6–14 C) | Edge-to-interior carbon ratio limits control; use ≥30 C for tight H/C |
| Steric clash count increases with size | Large flat aromatics have H···H near-contacts; use GROMACS energy minimisation after generation for production runs |
| Geometry validation thresholds | The built-in validator uses strict VdW radii; some reported "clashes" are artefacts of the flat starting structure and resolve under MD |
| Amorphous porous surfaces not yet implemented | Only slit pores (parallel sheets) are supported; `pore_type="amorphous"` is reserved for a future release |

---

## References

- Jorgensen, W. L. et al. "Development and Testing of the OPLS All-Atom Force Field." *J. Am. Chem. Soc.* 118.45 (1996): 11225–11236.
- RDKit: Open-source cheminformatics. https://www.rdkit.org
