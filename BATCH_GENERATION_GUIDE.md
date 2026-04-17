# Batch Generation Guide - Mixed Biochar Simulations

## Overview

The `generate_biochar_series()` function generates multiple biochar structures in a single call, perfect for:

- **Temperature series simulations** (BC400, BC600, BC800)
- **Composition studies** (varying H/C or O/C ratios)
- **Size-dependent studies** (10, 20, 50+ carbon atoms)
- **Mixed biochar systems** (multiple types in one simulation)

## Quick Example

```python
from biochar_generator import generate_biochar_series

# Define configurations
configs = [
    {"molecule_name": "BC400", "H_C_ratio": 0.65, "O_C_ratio": 0.20},
    {"molecule_name": "BC600", "H_C_ratio": 0.55, "O_C_ratio": 0.12},
    {"molecule_name": "BC800", "H_C_ratio": 0.40, "O_C_ratio": 0.05},
]

# Generate all structures
results = generate_biochar_series(
    configurations=configs,
    output_directory="my_simulation",
    create_combined_top=True,
)

# Ready for GROMACS!
# gmx grompp -p my_simulation/combined.top -o topol.tpr
```

## Configuration Format

Each configuration dictionary should contain:

### Required
- `molecule_name` (str): Residue name, **max 5 characters**

### Optional (defaults provided)
- `target_num_carbons` (int): Default 50
- `H_C_ratio` (float): Default 0.5
- `O_C_ratio` (float): Default 0.1
- `aromaticity_percent` (float): Default 90.0
- `functional_groups` (list): Default ["hydroxyl", "carboxyl", "ether"]
- `seed` (int): For reproducibility (no default - random if not set)

## Naming Conventions

### Temperature-Based (Recommended)
```python
configs = [
    {"molecule_name": "BC400", ...},  # 400°C pyrolysis
    {"molecule_name": "BC600", ...},  # 600°C pyrolysis
    {"molecule_name": "BC800", ...},  # 800°C pyrolysis
]
```

Advantages:
- Matches literature conventions
- Aligns with research standards
- Intuitive for temperature-dependent studies

### Composition-Based
```python
# H/C ratio series
configs = [
    {"molecule_name": "BCH04", "H_C_ratio": 0.4, ...},
    {"molecule_name": "BCH06", "H_C_ratio": 0.6, ...},
    {"molecule_name": "BCH08", "H_C_ratio": 0.8, ...},
]

# O/C ratio series
configs = [
    {"molecule_name": "BCO05", "O_C_ratio": 0.05, ...},
    {"molecule_name": "BCO12", "O_C_ratio": 0.12, ...},
    {"molecule_name": "BCO20", "O_C_ratio": 0.20, ...},
]
```

### Sequential
```python
configs = [
    {"molecule_name": "BC001", ...},
    {"molecule_name": "BC002", ...},
    {"molecule_name": "BC999", ...},
]
```

### Custom/Combined
```python
configs = [
    {"molecule_name": "BC4OL", ...},  # BC, ~400°C, Oxygen-Low
    {"molecule_name": "BC6OM", ...},  # BC, ~600°C, Oxygen-Medium
    {"molecule_name": "BC8OH", ...},  # BC, ~800°C, Oxygen-High
]
```

**Constraint**: All names must be ≤5 characters (GROMACS .gro format requirement)

## Output Files

For each configuration, generates:

```
output_directory/
├── bc400.gro          # Structure file (3D coordinates)
├── bc400.top          # Full topology with forcefield includes
├── bc400.itp          # Reusable molecule definition
├── bc600.gro
├── bc600.top
├── bc600.itp
├── bc800.gro
├── bc800.top
├── bc800.itp
└── combined.top       # Combined topology for mixed simulation
```

### combined.top Format

```
; Combined topology for mixed biochar simulation
#include "oplsaa.ff/forcefield.itp"

; Molecule definitions
#include "bc400.itp"
#include "bc600.itp"
#include "bc800.itp"

[ system ]
Mixed Biochar System

[ molecules ]
BC400  1
BC600  1
BC800  1
```

## Return Value

Returns a dictionary mapping molecule names to file paths:

```python
results = generate_biochar_series(configs)

# Access results
for mol_name, (gro_path, top_path, itp_path) in results.items():
    print(f"{mol_name}: {gro_path}")
```

## Example 1: Temperature Series

```python
from biochar_generator import generate_biochar_series

# Typical composition trends with pyrolysis temperature
configs = [
    {
        "molecule_name": "BC400",
        "target_num_carbons": 80,
        "H_C_ratio": 0.65,      # High H at low temp
        "O_C_ratio": 0.20,      # High O at low temp
        "seed": 401,
    },
    {
        "molecule_name": "BC600",
        "target_num_carbons": 85,
        "H_C_ratio": 0.55,      # Intermediate
        "O_C_ratio": 0.12,      # Intermediate
        "seed": 402,
    },
    {
        "molecule_name": "BC800",
        "target_num_carbons": 90,
        "H_C_ratio": 0.40,      # Low H at high temp
        "O_C_ratio": 0.05,      # Low O at high temp
        "seed": 403,
    },
]

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/temperature_series",
    create_combined_top=True,
)

print("Temperature series generated!")
print("Run with: gmx grompp -p output/temperature_series/combined.top")
```

## Example 2: Composition Study

```python
# Study effect of H/C ratio
configs = []
h_c_values = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

for i, h_c in enumerate(h_c_values):
    mol_name = f"BCH{int(h_c*10):02d}"  # BCH03, BCH04, BCH05, etc.
    configs.append({
        "molecule_name": mol_name,
        "target_num_carbons": 70,
        "H_C_ratio": h_c,
        "O_C_ratio": 0.1,
        "seed": 500 + i,
    })

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/h_c_study",
    create_combined_top=True,
)
```

## Example 3: With Ring Defects

```python
# Compare pure hexagonal vs defective structures
configs = [
    {
        "molecule_name": "BC800",
        "target_num_carbons": 80,
        "H_C_ratio": 0.40,
        "O_C_ratio": 0.05,
        "defect_fraction": 0.0,    # Pure hexagonal
        "seed": 501,
    },
    {
        "molecule_name": "BD805",  # BD = Biochar Disordered
        "target_num_carbons": 80,
        "H_C_ratio": 0.40,
        "O_C_ratio": 0.05,
        "defect_fraction": 0.05,   # 5% pentagons
        "seed": 502,
    },
    {
        "molecule_name": "BD815",
        "target_num_carbons": 80,
        "H_C_ratio": 0.40,
        "O_C_ratio": 0.05,
        "defect_fraction": 0.15,   # 15% pentagons
        "seed": 503,
    },
]

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/defect_study",
    create_combined_top=True,
)
```

## Example 4: Mixed System

```python
# Combine different biochar types in one simulation
configs = [
    {
        "molecule_name": "BC400",
        "target_num_carbons": 60,
        "H_C_ratio": 0.65,
        "O_C_ratio": 0.25,  # Oxidized biochar
        "seed": 601,
    },
    {
        "molecule_name": "BC800",
        "target_num_carbons": 100,
        "H_C_ratio": 0.35,
        "O_C_ratio": 0.03,  # Reduced biochar
        "seed": 602,
    },
    {
        "molecule_name": "BCHWP",  # Activated/water-treated
        "target_num_carbons": 80,
        "H_C_ratio": 0.5,
        "O_C_ratio": 0.30,
        "seed": 603,
    },
]

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/mixed_biochar",
    create_combined_top=True,
    verbose=True,
)
```

## Using in GROMACS

### Basic Workflow

```bash
# 1. Generate structures
python3 generate_biochar_series(configs)

# 2. Prepare system (using combined.top)
gmx grompp -f md.mdp -p combined.top -c combined.gro -o topol.tpr

# 3. Run simulation
gmx mdrun -deffnm topol -v
```

### Combined Structure File (Optional)

If you need a single `.gro` file with all structures:

```bash
# Use genconf to create periodic/replicated structure
gmx genconf -f bc400.gro -nbox 1 1 1 -o system.gro
```

Then edit `combined.top` to reference the combined structure.

## API Reference

```python
generate_biochar_series(
    configurations: List[Dict],           # List of config dicts
    output_directory: str = ".",          # Output path
    create_combined_top: bool = True,     # Generate combined.top
    verbose: bool = True,                 # Print progress
) -> Dict[str, Tuple[Path, Path, Path]]:
    """
    Returns: {
        "BC400": (Path(bc400.gro), Path(bc400.top), Path(bc400.itp)),
        "BC600": (Path(bc600.gro), Path(bc600.top), Path(bc600.itp)),
        ...
    }
    """
```

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `configurations` | List[Dict] | Required | List of config dicts |
| `output_directory` | str | "." | Output directory path |
| `create_combined_top` | bool | True | Generate combined topology |
| `verbose` | bool | True | Print progress information |

## Configuration Dict Keys

| Key | Type | Default | Notes |
|-----|------|---------|-------|
| `molecule_name` | str | Required | Max 5 chars for .gro format |
| `target_num_carbons` | int | 50 | Target structure size |
| `H_C_ratio` | float | 0.5 | Target H/C ratio |
| `O_C_ratio` | float | 0.1 | Target O/C ratio |
| `aromaticity_percent` | float | 90.0 | Target aromaticity % |
| `functional_groups` | dict | default | Custom functional groups |
| `defect_fraction` | float | 0.0 | Probability each ring is a pentagon |
| `seed` | int | random | Random seed (reproducible if set) |

## Best Practices

### 1. Naming
- Use descriptive ≤5 character names
- Match naming convention to your study type
- Include seed for reproducibility

### 2. Parameters
- Set realistic H/C and O/C ratios
- Use different seeds for each structure (or omit for random)
- Keep target_num_carbons reasonable (10-1000 atoms)

### 3. Output
- Always set `create_combined_top=True` for mixed simulations
- Check generated files before using in GROMACS
- Review validation reports (check for warnings)

### 4. GROMACS Integration
```bash
# Check topology before simulation
gmx grompp -pp processed.top -p combined.top -f md.mdp -o topol.tpr
```

## Troubleshooting

### "molecule_name exceeds 5 character limit"
Solution: Use shorter names (BC400, BCH05, etc.)

### Validation warnings for small structures
This is expected for small molecules. Either:
- Increase `target_num_carbons`
- Increase tolerance in validation settings
- Accept as-is (still valid for MD)

### File count mismatch
Check that all configs have unique `molecule_name` values.

## Running the Examples

```bash
# Run all batch generation examples
python3 examples/batch_generation.py

# This generates:
# - output/temperature_series/ (BC400, BC600, BC800)
# - output/composition_series/ (BCH04, BCH06, BCH08)
# - output/oxygen_series/ (BCO05, BCO12, BCO20)
# - output/size_series/ (BC10, BC20, BC30)
# - output/mixed_system/ (BC400, BC800, BCHWP)
```

Each directory contains individual `.gro`, `.top`, `.itp` files plus a `combined.top` ready for GROMACS simulations.
