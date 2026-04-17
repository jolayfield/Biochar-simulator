# Biochar Simulator — Best Practices & Known Limitations

## Quick Reference

| Feature | Status | Notes |
|---------|--------|-------|
| Carbon count matching (6–16C) | ✅ Exact | 100% — uses verified PAH templates |
| Carbon count matching (17–50C) | ✅ Good | 80–94% of requested size |
| Carbon count matching (50–200C) | ⚠ Acceptable | 70–90% of requested size |
| Hydrogen assignment | ✅ Fixed | All atoms properly saturated |
| Valence enforcement (min + max) | ✅ Complete | Detects under- and over-saturation |
| Oxygen assignment (–OH groups) | ✅ Working | Added to aromatic carbons with available valence |
| GROMACS file output | ✅ Working | .gro, .top, .itp in nm units |
| Relative imports | ✅ Fixed | Use `from src.module import ...` |

---

## Correct Import Pattern

Always run scripts from the **project root** and import from the `src` package:

```python
import sys
from pathlib import Path

# Add project root (NOT src/) to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import from the src package
from src.biochar_generator import BiocharGenerator, GeneratorConfig, generate_biochar
from src.biochar_generator import generate_biochar_series
```

> ⚠ **Common mistake**: Adding `src/` to the path and importing `from biochar_generator import ...` will
> fail with `ImportError: attempted relative import with no known parent package`.

---

## Carbon Count Recommendations

### How Carbon Count Works

The generator uses two methods depending on target size:

| Target Size | Method | Expected Match |
|------------|--------|----------------|
| ≤ 6C | Benzene template | **100%** |
| 7–10C | Naphthalene template | **100%** |
| 11–14C | Anthracene template | **100%** |
| 15–16C | Pyrene template | **100%** |
| 17–50C | Pyrene + benzene rings | **80–94%** |
| 50C+ | Modular PAH assembly | **70–90%** |

### Best Practice: Target Sizes

For the most predictable results, use targets that match a PAH template exactly:

```python
RELIABLE_TARGETS = [6, 10, 14, 16]  # Exact matches guaranteed

# These work well (80–95% match):
GOOD_TARGETS = list(range(20, 60, 6))  # 20, 26, 32, 38, 44, 50...

# For exact sizes, set seed and verify the actual count
mol, _, _, _, _ = generate_biochar(target_num_carbons=50, seed=42, ...)
actual_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
print(f"Requested 50C, got {actual_c}C")
```

### Compensating for Size Discrepancy

Because the generator may produce slightly fewer carbons than requested, you can overshoot the target slightly:

```python
# If you want ~30 actual carbons, request ~35
mol, _, _, _, _ = generate_biochar(
    target_num_carbons=35,  # Request slightly more
    ...
)

actual_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
print(f"Actual carbons: {actual_c}")
```

---

## Valence and Chemical Validity

### Guaranteed Valid Structures

As of v1.2, all generated structures satisfy both minimum and maximum valence:

- **Carbon**: exactly 4 bonds
- **Hydrogen**: exactly 1 bond
- **Oxygen**: exactly 2 bonds
- **Nitrogen**: exactly 3 bonds (if present)

### How to Verify

```python
from src.valence import ValenceValidator

# Check after generation
is_valid, errors = ValenceValidator.validate_molecule(mol)

if is_valid:
    print("✓ All atoms have valid valence")
else:
    for error in errors:
        print(f"  ✗ {error}")
```

### Interpreting the Valence Report

```python
ValenceValidator.print_valence_report(mol)
```

Output:
```
Idx  Sym Min  Max  Bonds  Need  Avail Valid  Charge
0    C   4    4    4      0     0     ✓ OK   0
6    O   2    2    2      0     0     ✓ OK   0
12   H   1    1    1      0     0     ✓ OK   0
```

| Column | Meaning |
|--------|---------|
| Min | Minimum bonds required |
| Max | Maximum bonds allowed |
| Bonds | Current bond count |
| Need | Additional bonds needed to reach Min |
| Avail | Bonds available before reaching Max |
| Valid | ✓ OK if Min ≤ Bonds ≤ Max |

---

## Composition Ratios

### Why Ratios May Drift from Target

Three factors affect composition accuracy:

1. **Carbon skeleton size** — If fewer carbons are generated than requested, the absolute number of H/O atoms is also smaller.
2. **Aromatic saturation** — Aromatic carbons naturally attract fewer hydrogens than aliphatic ones.
3. **Oxygen placement** — Only aromatic carbons with available valence receive –OH groups.

### Expected Composition Behavior

| Scenario | Effect |
|----------|--------|
| Small skeleton (≤16C) | H/C tends to be higher due to edge-heavy aromatic rings |
| Large skeleton (50C+) | H/C closer to target; more interior carbons |
| High O/C target | May not be reached if few aromatic carbons have available valence |
| Low O/C target (< 0.05) | Usually well-matched |

### Practical H/C Ranges

The physical chemistry of biochar imposes natural constraints:

| Pyrolysis Temperature | H/C Range | O/C Range |
|-----------------------|-----------|-----------|
| 300–400 °C (BC400) | 0.60–0.75 | 0.15–0.25 |
| 500–600 °C (BC600) | 0.40–0.60 | 0.08–0.15 |
| 700–800 °C (BC800) | 0.20–0.45 | 0.02–0.08 |

---

## GROMACS Output

### Residue Name Rules

GROMACS limits residue names to **≤5 characters**. The generator enforces this:

```python
# ✓ Valid residue names
config = GeneratorConfig(molecule_name="BC")
config = GeneratorConfig(molecule_name="BC400")
config = GeneratorConfig(molecule_name="BC6OL")

# ✗ Invalid — will raise an error
config = GeneratorConfig(molecule_name="BIOCHAR")  # 7 chars
```

### Recommended Naming Convention

| Series | Template | Examples |
|--------|----------|---------|
| Pyrolysis temperature | BC{T} | BC400, BC600, BC800 |
| H/C ratio | BCH{H} | BCH04, BCH06, BCH08 |
| O/C ratio | BCO{O} | BCO05, BCO10, BCO20 |
| Sequential batch | BC{N} | BC001, BC002, BC003 |
| Mixed descriptor | BC{X} | BC4OL (400°C, low O) |

### Unit System

All output is in GROMACS standard SI units:

| Quantity | Unit |
|---------|------|
| Coordinates | nanometers (nm) |
| Box vectors | nanometers (nm) |
| Masses | u (amu) |
| Charges | e (elementary charge) |

---

## Running the Examples

Always run from the **project root**:

```bash
# From the project root directory:
cd /path/to/Biochar-simulator

# Single structure examples
python examples/example_usage.py

# Batch generation examples
python examples/batch_generation.py
```

---

## Validation Messages Explained

### "Validation FAILED" vs Errors

The validator produces two types of messages:

| Type | Meaning | Action |
|------|---------|--------|
| `Warning: Composition validation issues` | Ratio slightly outside tolerance | Normal for small structures; check actual ratios |
| `Validation FAILED` | Ratio or valence outside limits | Review structure; may still be usable |
| `Atom X (C): Valence 3 below minimum 4` | Under-saturated carbon | Bug — report it; fixed in v1.2 |
| `Atom X (C): Valence 5 exceeds maximum 4` | Over-saturated carbon | Bug — report it |

### Common Validation Warnings

```
H/C ratio 0.750 exceeds tolerance (target: 0.500, error: 50.0%)
```
**Cause**: Small aromatic skeleton (e.g., pyrene) has high H/C due to peripheral H atoms.
**Solution**: Use larger target size, or accept the actual composition for simulation.

```
O/C ratio 0.000 exceeds tolerance (target: 0.100, error: 100.0%)
```
**Cause**: No aromatic carbons with available valence for –OH attachment.
**Solution**: Use a less aromatic starting structure, or reduce O/C target.

---

## Performance Notes

| Target Size | Generation Time | Notes |
|------------|-----------------|-------|
| ≤ 22C | < 1 second | Single PAH template |
| 22–100C | 1–5 seconds | Template + RDKit 3D embedding |
| 100–500C | 5–30 seconds | Large 3D embedding, slow sanitization |
| 500C+ | Minutes | May hit RDKit memory/conformer limits |

### For Large Structures

```python
# Use a seed for reproducibility
# Set periodic_box=True for large simulations
config = GeneratorConfig(
    target_num_carbons=200,
    periodic_box=True,
    box_size=[5.0, 5.0, 5.0],  # nm
    seed=42,
)
```

---

## Troubleshooting

### `Bad Conformer Id` Error

**Cause**: RDKit cannot generate a 3D conformer for a very large or unusual molecule.
**Fix**: Reduce target size, or use `seed` parameter to try different structures.

### `ImportError: attempted relative import with no known parent package`

**Cause**: Script is directly importing from inside `src/` directory.
**Fix**: Add project root to `sys.path`, not `src/`:

```python
sys.path.insert(0, str(Path(__file__).parent.parent))
from src.biochar_generator import generate_biochar
```

### `Can't kekulize mol`

**Cause**: RDKit cannot assign alternating single/double bonds to the aromatic system.
**Fix**: This is an internal error; the generator falls back to a simpler structure. The output will be valid, but may be smaller than requested.

### Atoms with 3 bonds instead of 4

**Cause**: This was a known bug in hydrogen assignment. **Fixed in v1.2**.
**Fix**: Update to v1.2 — the new `HydrogenAssigner._saturate_valences()` method ensures all atoms reach minimum valence.

---

## Pentagon Ring Defects (defect_fraction)

### What Are Ring Defects?

By default, all ring additions during graph growth are hexagons (6 carbons). Setting `defect_fraction > 0` inserts occasional pentagons (5 carbons), producing topologically disordered structures that mimic amorphous graphitic biochar.

### Usage

```python
# 10% of rings will be pentagons
mol, coords, gro, top, itp = generate_biochar(
    target_num_carbons=60,
    defect_fraction=0.10,
    seed=42,
)

# In batch generation
configs = [
    {"molecule_name": "BC800", "defect_fraction": 0.0},  # pure hexagonal
    {"molecule_name": "BC8D5", "defect_fraction": 0.05}, # 5% defects
    {"molecule_name": "BC8D15", "defect_fraction": 0.15}, # 15% defects
]

results = generate_biochar_series(configs)
```

### Key Points

- **Range**: 0.0 (pure hexagon) to 1.0 (all pentagons)
- **Typical values**: 0.05–0.20 (5–20% pentagons)
- **Parity handling**: Pentagon additions (+3 nodes, odd) automatically fix parity constraints
- **Kekulization**: Non-bipartite graphs may fail kekulization; up to 5 retries with different sub-seeds
- **Reproducibility**: Use same `seed` to get identical defect patterns

### Effects

- Introduces curvature and topological disorder
- Maintains aromaticity and valence validity
- Produces slightly different atom counts than pure hexagon mode
- More realistic model of amorphous biochar

---

## Known Limitations (v1.2)

1. **Carbon count matching for large targets**: Structures with >50 requested carbons may be smaller than requested (70–90% accuracy). The generator still produces chemically valid structures.

2. **Chrysene/coronene SMILES**: The SMILES strings for chrysene (C18) and coronene (C24) in the PAH library are not validated. The generator avoids these and falls back to smaller templates.

3. **Functional groups beyond –OH**: The oxygen assigner currently only places hydroxyl groups (–OH). Carboxyl, ether, carbonyl, lactone, and quinone groups are defined but not yet fully wired in.

4. **Very large structures (>500C)**: May fail with `Bad Conformer Id` from RDKit's 3D embedding. Use periodic boundary boxes and split into multiple smaller molecules.

5. **H/C ratio control**: For small structures (≤22C), the H/C ratio is largely determined by the PAH template. The hydrogen trimming may not always reach the exact target.

---

## Version History

| Version | Key Changes |
|---------|-------------|
| v1.0 | Initial release — PAH assembly, GROMACS output |
| v1.1 | Valence system, coordinate unit fix (Å → nm) |
| v1.2 | **Pentagon ring defects** (`defect_fraction`), **porous surface generation** (`generate_surface`), min+max valence enforcement, hydrogen assignment fix |
