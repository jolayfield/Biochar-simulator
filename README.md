# Biochar Simulator - Structure Generator

A Python package for generating realistic biochar molecular structures for GROMACS molecular dynamics simulations.

## Overview

The Biochar Simulator creates building blocks for biochar structures with user-specified compositional parameters:

- **Size**: Target number of carbon atoms (adjustable from 10 to 10,000+)
- **Composition**: H/C and O/C ratios (with configurable tolerances)
- **Aromaticity**: Target aromatic carbon percentage
- **Functional Groups**: Carboxyl, hydroxyl, ether, carbonyl, phenolic, lactone, quinone
- **Output**: GROMACS-ready `.gro`, `.top`, and `.itp` files using OPLS-AA forcefield

## Features

- **Modular Architecture**: Independent components for each generation stage
- **PAH-Based Generation**: Uses Polycyclic Aromatic Hydrocarbons as building blocks
- **3D Geometry**: Generates realistic 3D coordinates with planarity constraints
- **OPLS-AA Typing**: Automatic atom type and partial charge assignment
- **Comprehensive Validation**: Multi-stage validation of composition, chemistry, and geometry
- **Reproducible**: Supports seeded RNG for reproducible structure generation
- **Flexible**: Supports isolated molecules and periodic boxes

## Installation

### Requirements

- Python 3.8+
- RDKit (for molecular structure manipulation)
- NumPy, SciPy (numerical computing)
- NetworkX (graph operations)
- Pandas (data handling)

### Setup

```bash
# Install dependencies
pip install -r requirements.txt

# Install package (from project root)
pip install -e .
```

## Quick Start

### Basic Usage (Single Structure)

```python
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.biochar_generator import generate_biochar

# Generate a 16-carbon biochar (exact PAH match)
mol, coords, gro_path, top_path, itp_path = generate_biochar(
    target_num_carbons=16,
    H_C_ratio=0.5,
    O_C_ratio=0.1,
    aromaticity_percent=90.0,
    output_directory="output",
    basename="my_biochar",
    molecule_name="BC",  # Residue name (max 5 chars for GROMACS)
)

print(f"Generated GROMACS files:")
print(f"  Structure: {gro_path}")
print(f"  Topology: {top_path}")
print(f"  Include:  {itp_path}")
```

> **Tip**: Use target sizes of 6, 10, 14, or 16 carbons for 100% carbon count accuracy.
> See [BEST_PRACTICES.md](BEST_PRACTICES.md) for size recommendations.

### Batch Generation (Multiple Structures)

For mixed simulations (temperature series, composition studies, etc.):

```python
from src.biochar_generator import generate_biochar_series

# Generate temperature series
configs = [
    {"molecule_name": "BC400", "H_C_ratio": 0.65, "O_C_ratio": 0.20},
    {"molecule_name": "BC600", "H_C_ratio": 0.55, "O_C_ratio": 0.12},
    {"molecule_name": "BC800", "H_C_ratio": 0.40, "O_C_ratio": 0.05},
]

results = generate_biochar_series(
    configurations=configs,
    output_directory="output/temperature_series",
    create_combined_top=True,  # Auto-generate combined.top
)

# Use in GROMACS:
# gmx grompp -f md.mdp -c combined.gro -p combined.top -o topol.tpr
```

### Advanced Configuration

```python
from src.biochar_generator import BiocharGenerator, GeneratorConfig

config = GeneratorConfig(
    target_num_carbons=200,
    H_C_ratio=0.55,
    H_C_tolerance=0.15,
    O_C_ratio=0.15,
    O_C_tolerance=0.12,
    aromaticity_percent=92.0,
    aromaticity_tolerance=5.0,
    functional_groups=["hydroxyl", "carboxyl", "ether", "phenolic"],
    periodic_box=False,
    seed=42,  # For reproducibility
)

generator = BiocharGenerator(config)
mol, coords, composition = generator.generate()
generator.print_summary()

# Export to GROMACS
gro_path, top_path, itp_path = generator.export_gromacs(
    output_directory="output",
    basename="large_biochar",
)
```

## Architecture

### Component Modules

1. **`constants.py`** - OPLS-AA parameters, atom types, functional group definitions
2. **`carbon_skeleton.py`** - Aromatic carbon framework generation (PAHs)
3. **`heteroatom_assignment.py`** - Hydrogen and oxygen placement
4. **`geometry_3d.py`** - 3D coordinate generation with RDKit
5. **`opls_typing.py`** - OPLS-AA atom typing and charge assignment
6. **`gromacs_export.py`** - Writing `.gro`, `.top`, `.itp` files
7. **`validation.py`** - Composition, chemistry, and structure validation
8. **`biochar_generator.py`** - Main orchestrator and user API

### Generation Pipeline

```
Input Parameters
    ↓
Carbon Skeleton Generation (PAHs)
    ↓
Oxygen Assignment (functional groups)
    ↓
Hydrogen Assignment (fill valences)
    ↓
3D Geometry Generation (RDKit + constraints)
    ↓
OPLS-AA Typing & Charges
    ↓
Validation (multi-stage)
    ↓
GROMACS Export (.gro, .top, .itp)
```

## Configuration Parameters

### Size & Structure

- `target_num_carbons` (int): Target number of carbon atoms (10-10000)
- `aromaticity_percent` (float): Target % aromatic carbons (0-100)
- `aromaticity_tolerance` (float): Aromaticity tolerance in %

### Composition

- `H_C_ratio` (float): Target H/C molar ratio
- `H_C_tolerance` (float): H/C tolerance (e.g., 0.10 = ±10%)
- `O_C_ratio` (float): Target O/C molar ratio
- `O_C_tolerance` (float): O/C tolerance (e.g., 0.10 = ±10%)

### Functional Groups

- `functional_groups` (list): Types to include. Options:
  - `"hydroxyl"` - Hydroxyl groups (-OH)
  - `"carboxyl"` - Carboxylic acid (-COOH)
  - `"ether"` - Ether linkages (-O-)
  - `"carbonyl"` - Carbonyl (C=O)
  - `"phenolic"` - Aromatic hydroxyl
  - `"lactone"` - Cyclic ester
  - `"quinone"` - Quinone (C=O on aromatic)

### System Setup

- `periodic_box` (bool): Generate periodic boundary box
- `box_size` (array): Custom box dimensions [Lx, Ly, Lz]
- `seed` (int): Random seed for reproducibility

## Unit System

All GROMACS output files use **SI units** (the GROMACS standard):

| Quantity | Unit | Notes |
|----------|------|-------|
| Coordinates | nanometers (nm) | Converted from RDKit Ångströms |
| Box vectors | nanometers (nm) | 1 Å = 0.1 nm |
| Masses | atomic mass units (u) | Standard |
| Charges | elementary charge (e) | From OPLS-AA |
| Distances | nanometers (nm) | Bond lengths, etc. |

The generator automatically converts from RDKit's Ångströms (Å) to GROMACS's nanometers (nm).

## Output Files

### `.gro` File (Structure)
GROMACS structure format containing:
- Atom positions (in **nanometers**, nm)
- Box vectors (in **nanometers**, nm)
- Velocity information (optional, zeros)

### `.top` File (Topology)
Full GROMACS topology including:
- Moleculetype definition
- Atom definitions with types and charges
- Bond/angle/dihedral parameters
- System and molecules sections

### `.itp` File (Include Topology)
Reusable molecule definition for inclusion in larger simulations.

## Batch Generation for Mixed Simulations

The `generate_biochar_series()` function is designed for generating multiple biochar structures for mixed simulations. This is ideal for:

- **Temperature Series**: BC400, BC600, BC800 (different pyrolysis temperatures)
- **Composition Studies**: BCH04, BCH06, BCH08 (varying H/C ratios)
- **Oxygen Content**: BCO05, BCO12, BCO20 (varying O/C ratios)
- **Size Studies**: BC10, BC20, BC30 (different carbon counts)
- **Mixed Systems**: Combined biochar types in single simulation

### Naming Convention (GROMACS .gro Format)

Residue names must be **≤5 characters** (GROMACS requirement):

| Name | Purpose | Example |
|------|---------|---------|
| Temperature | Pyrolysis temperature | BC400, BC600, BC800 |
| Composition | H/C or O/C ratio | BCH05, BCO10 |
| Sequential | Batch numbering | BC001, BC002, BC003 |
| Custom | Mixed types | BC4OL, BC6OM, BC8OH |

### Example: Mixed Simulation

```python
from biochar_generator import generate_biochar_series

# Define biochar configurations
configs = [
    {
        "molecule_name": "BC400",
        "target_num_carbons": 80,
        "H_C_ratio": 0.65,
        "O_C_ratio": 0.20,
        "seed": 101,
    },
    {
        "molecule_name": "BC600",
        "target_num_carbons": 85,
        "H_C_ratio": 0.55,
        "O_C_ratio": 0.12,
        "seed": 102,
    },
    {
        "molecule_name": "BC800",
        "target_num_carbons": 90,
        "H_C_ratio": 0.40,
        "O_C_ratio": 0.05,
        "seed": 103,
    },
]

# Generate all structures and combined topology
results = generate_biochar_series(
    configurations=configs,
    output_directory="output/mixed_biochar",
    create_combined_top=True,
    verbose=True,
)

# The function automatically creates:
# - Individual: bc400.gro, bc400.top, bc400.itp
# - Individual: bc600.gro, bc600.top, bc600.itp
# - Individual: bc800.gro, bc800.top, bc800.itp
# - Combined: combined.top (references all .itp files)
```

### Using Combined Topology in GROMACS

```bash
cd output/mixed_biochar

# Prepare system for simulation
gmx grompp -f ../md.mdp -p combined.top -o topol.tpr

# Run molecular dynamics
gmx mdrun -deffnm topol
```

## Examples

See examples directory:

- **`examples/example_usage.py`** - Basic structure generation with different sizes/compositions
- **`examples/batch_generation.py`** - Batch generation for mixed simulations (5 complete examples)

```bash
cd examples
python3 batch_generation.py
```

Generates 5 example batches:
1. Temperature series (BC400, BC600, BC800)
2. H/C composition series (BCH04, BCH06, BCH08)
3. Oxygen content series (BCO05, BCO12, BCO20)
4. Size series (BC10, BC20, BC30)
5. Custom mixed system (BC400, BC800, BCHWP)

## Testing

Run unit tests:

```bash
pytest tests/test_generator.py -v
```

Tests cover:
- OPLS constants and PAH library
- Carbon skeleton generation
- Heteroatom assignment
- 3D geometry generation
- OPLS typing and charges
- Validation engine
- Complete generator workflow

## Output Structure

```
output/
├── biochar_50C.gro          # Structure file
├── biochar_50C.top          # Topology file
├── biochar_50C.itp          # Include file
├── biochar_large_500C.gro
├── biochar_large_500C.top
└── biochar_large_500C.itp
```

## Validation

The generator performs comprehensive multi-stage validation:

### Three-Stage Validation

**Stage 1: Valence Checking** ⭐ NEW
- Ensures all atoms have proper bond counts
- H: max 1 bond, C: max 4, O: max 2, N: max 3, S: max 2-6
- Prevents chemically invalid structures before geometry generation
- See [VALENCE_SYSTEM.md](VALENCE_SYSTEM.md) for details

**Stage 2: Composition Validation**
- Checks H/C and O/C ratios against targets

1. **Composition Validation**: Checks H/C and O/C ratios against targets
2. **Chemical Feasibility**: Validates valences, atom types, charges
3. **Structure Quality**: Checks connectivity, geometry, aromaticity

All validation issues are reported as errors or warnings in the summary.

## Known Limitations

| Issue | Status | Workaround |
|-------|--------|------------|
| Carbon count accuracy for large targets | ⚠ 70–90% for targets > 50C | Request slightly more than needed |
| H/C ratio accuracy for small structures | ⚠ Limited by PAH edge geometry | Use targets ≥ 30C for better control |
| Functional groups (only –OH placed) | ⚠ Others defined but not wired in | Use hydroxyl for oxygen content |
| Very large structures (> 500C) | ⚠ May fail with conformer error | Split into multiple smaller molecules |
| Chrysene/Coronene SMILES | ⚠ Currently unverified | Generator falls back to pyrene (16C) |

> For detailed guidance, see [BEST_PRACTICES.md](BEST_PRACTICES.md).

## OPLS-AA Support

The package includes comprehensive OPLS-AA atom type definitions:

- **Aromatic**: CA (aromatic C), HA (aromatic H)
- **Aliphatic**: CT (sp³ C), HC (H on aliphatic C)
- **Oxygens**: OH (hydroxyl), OC (carbonyl), OS (ether), O (carboxylic)
- **Other**: N, NT, S, SH (for future extensions)

Partial charges and bond/angle parameters are automatically assigned from OPLS-AA tables.

## Known Limitations

1. **PAH Combination**: Current implementation uses simple PAH assembly. More sophisticated fusion algorithms could improve size coverage.

2. **Functional Group Placement**: Groups are added greedily. More intelligent placement considering spatial constraints could improve results.

3. **Large Molecules**: Structures > 2000 atoms may be slower to generate due to 3D coordinate optimization.

4. **Geometry Relaxation**: Uses RDKit MMFF94/UFF force fields. For production, GROMACS energy minimization is recommended.

## Future Enhancements

- [ ] Better PAH fusion algorithms
- [ ] Machine learning-based atom typing refinement
- [ ] Interface with scikit-nano for additional structure generation
- [ ] Support for heteroatom-doped biochar (N, S)
- [ ] Batch structure generation
- [ ] Integration with GROMACS for direct validation
- [ ] Database of pre-computed structures

## References

The biochar simulator is inspired by research on molecular modeling of biochar and related carbon materials:

- Jorgensen, W. L., et al. "Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids." *J. Am. Chem. Soc.* 118.45 (1996): 11225-11236.
- Wood, R., Masek, O., & Erastova, V. (Referenced in parent project)

## License

This project is provided as-is for research purposes.

## Contributing

To contribute improvements:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Submit a pull request

## Support

For issues, questions, or contributions, please check:
- README files in each module
- Docstrings in source code
- Example scripts in `examples/`
- Unit tests in `tests/`
