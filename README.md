# Biochar Simulator — Structure Generator

A Python package for generating realistic biochar molecular structures for GROMACS molecular dynamics simulations. Supports structures from 6 to 200+ carbons with 100% aromatic carbon content and accurate H/C and O/C ratios.

---

## Overview

The Biochar Simulator builds polycyclic aromatic hydrocarbon (PAH) structures with user-specified compositional parameters and exports GROMACS-ready force field files.

**Key capabilities:**

- **Size**: 6 to 200+ carbons (exact count for most targets)
- **Composition**: H/C and O/C ratios with configurable tolerances
- **Aromaticity**: 100% aromatic skeleton (graphene-nanoflake topology)
- **Functional Groups**: Hydroxyl, carboxyl, ether, carbonyl, phenolic, lactone, quinone
- **Output**: GROMACS `.gro`, `.top`, `.itp` files with OPLS-AA force field

---

## Installation

### Requirements

- Python 3.8+
- RDKit
- NumPy, SciPy
- NetworkX

```bash
pip install -r requirements.txt
```

---

## Quick Start

### Single structure

```python
from src.biochar_generator import generate_biochar

mol, coords, gro_path, top_path, itp_path = generate_biochar(
    target_num_carbons=100,
    H_C_ratio=0.5,
    O_C_ratio=0.1,
    output_directory="output",
    basename="biochar_100C",
    molecule_name="BC",   # residue name, max 5 chars (GROMACS requirement)
    seed=42,
)
```

### Full configuration

```python
from src.biochar_generator import BiocharGenerator, GeneratorConfig

config = GeneratorConfig(
    target_num_carbons=200,
    H_C_ratio=0.5,
    H_C_tolerance=0.10,      # ±10%
    O_C_ratio=0.1,
    O_C_tolerance=0.10,
    aromaticity_percent=90.0,
    functional_groups=["hydroxyl", "carboxyl", "ether"],
    seed=42,
)

gen = BiocharGenerator(config)
mol, coords, composition = gen.generate()
gen.print_summary()
gro, top, itp = gen.export_gromacs(output_directory="output", basename="my_biochar")
```

### Batch generation (temperature/composition series)

```python
from src.biochar_generator import generate_biochar_series

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

## Configuration Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `target_num_carbons` | int | 50 | Target carbon count |
| `H_C_ratio` | float | 0.5 | Target H/C molar ratio |
| `H_C_tolerance` | float | 0.10 | Allowed H/C error (fraction) |
| `O_C_ratio` | float | 0.1 | Target O/C molar ratio |
| `O_C_tolerance` | float | 0.10 | Allowed O/C error (fraction) |
| `aromaticity_percent` | float | 90.0 | Target % aromatic carbons |
| `functional_groups` | list | `["hydroxyl","carboxyl","ether"]` | Groups to place |
| `molecule_name` | str | `"BC"` | Residue name (≤5 chars) |
| `periodic_box` | bool | False | Add periodic box vectors |
| `seed` | int | None | RNG seed for reproducibility |

### Available functional groups

`"hydroxyl"`, `"carboxyl"`, `"phenolic"`, `"ether"`, `"carbonyl"`, `"lactone"`, `"quinone"`

---

## Generation Pipeline

```
target_num_carbons
      │
      ▼
Carbon skeleton (PAH seed + ring-growth)
      │  Parity-aware hex-lattice builder; 100% aromatic for any size
      ▼
Oxygen assignment (hydroxyl → target O/C ratio)
      │
      ▼
Hydrogen assignment (fill valences → target H/C ratio)
      │
      ▼
3D geometry
      │  ≤80 heavy atoms: ETKDGv3 / ETKDGv2 embedding + MMFF94
      │  >80 heavy atoms: 2D-first embedding (flat graphene sheet) + FF minimization
      ▼
OPLS-AA atom typing & partial charges
      │
      ▼
Validation (composition, geometry, steric clashes)
      │
      ▼
GROMACS export (.gro / .top / .itp)
```

---

## Supported Sizes

The skeleton builder uses a pre-validated PAH library for exact matches and parity-aware ring growth for all other sizes:

| Range | Strategy | Aromaticity | Count accuracy |
|---|---|---|---|
| 6, 10, 14, 16, 18, 24 C | Exact library match | 100% | Exact |
| 17–50 C | Library seed + 4-node ring growth | 100% | Exact |
| 51–200+ C | Hex-lattice seed + ring growth | 100% | ≤5% error |

### PAH library

| Molecule | Carbons |
|---|---|
| benzene | 6 |
| naphthalene | 10 |
| anthracene / phenanthrene | 14 |
| pyrene | 16 |
| chrysene / triphenylene | 18 |
| coronene | 24 |

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
| `.gro` | GROMACS structure | Atom positions in **nm** |
| `.top` | GROMACS topology | Atoms, bonds, angles, dihedrals |
| `.itp` | Include topology | Reusable molecule definition |

Coordinates are in **nanometers** (GROMACS convention; RDKit Å values × 0.1).

---

## Source Modules

| Module | Responsibility |
|---|---|
| `carbon_skeleton.py` | PAH library lookup and ring-growth engine |
| `heteroatom_assignment.py` | Oxygen and hydrogen placement |
| `geometry_3d.py` | 3D coordinate generation and clash resolution |
| `opls_typing.py` | OPLS-AA atom types and partial charges |
| `gromacs_export.py` | `.gro` / `.top` / `.itp` file writer |
| `validation.py` | Composition, chemistry, and geometry checks |
| `constants.py` | OPLS-AA parameters, PAH library, VdW radii |
| `biochar_generator.py` | Public API (`BiocharGenerator`, `generate_biochar`) |

---

## Testing

```bash
# Quality test suite (PAH sizes 6–200C, compositions, seeds, large structures)
python3 tests/test_pah_quality.py
```

The suite reports:
- PAH library SMILES validity
- Kekulizability (100% aromatic required)
- Atom count accuracy (target vs actual)
- Bond errors and steric clashes per size
- Ring planarity (Å deviation)
- H/C and O/C ratio accuracy across compositions

---

## Known Limitations

| Issue | Notes |
|---|---|
| Only hydroxyl oxygens are placed | Other functional groups (carboxyl, ether…) are defined in constants but not yet wired into the oxygen assigner |
| H/C ratio loose for very small structures (6–14 C) | Edge-to-interior carbon ratio limits control; use ≥30 C for tight H/C |
| Steric clash count increases with size | Large flat aromatics have many H···H near-contacts; use GROMACS energy minimisation after generation for production runs |
| Geometry validation thresholds | The built-in validator uses strict VdW radii; some reported "clashes" are artefacts of the 2D-flat starting structure and resolve under MD |

---

## Examples

```bash
# Single structure examples
python3 examples/example_usage.py

# Batch generation (5 series: temperature, composition, oxygen, size, mixed)
python3 examples/batch_generation.py
```

---

## References

- Jorgensen, W. L. et al. "Development and Testing of the OPLS All-Atom Force Field." *J. Am. Chem. Soc.* 118.45 (1996): 11225–11236.
- RDKit: Open-source cheminformatics. https://www.rdkit.org
