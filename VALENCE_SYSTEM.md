# Valence Validation System

## Overview

A comprehensive valence checking system ensures that all atoms in generated biochar structures maintain proper chemical bonding. This prevents generation of invalid molecules that would fail GROMACS simulations.

## What is Valence?

**Valence** is the maximum number of chemical bonds an atom can form based on its atomic structure.

## Standard Valences

The system enforces **both minimum and maximum** valence constraints to ensure chemically valid structures.

| Element | Symbol | Min Bonds | Max Bonds | Examples |
|---------|--------|-----------|-----------|----------|
| Hydrogen | H | 1 | 1 | H-C, H-O |
| Carbon | C | 4 | 4 | C-C, C=O, C-H, C-O |
| Nitrogen | N | 3 | 3 | N-C, N=C, N-H |
| Oxygen | O | 2 | 2 | O-C, O-H, C-O-C |
| Sulfur | S | 2 | 2 | S-C, S=C (standard) |
| Chlorine | Cl | 1 | 1 | C-Cl |

## Components

### 1. ValenceValidator
Core validation class that checks atom valences.

```python
from valence import ValenceValidator

# Get valence info for specific atom
val_info = ValenceValidator.get_valence_info(mol, atom_idx=0)
# Returns: ValenceInfo with atom_idx, symbol, current_bonds, max_valence,
#          available_valence, is_valid, formal_charge

# Validate entire molecule
is_valid, errors = ValenceValidator.validate_molecule(mol)

# Print detailed report
ValenceValidator.print_valence_report(mol)
```

### 2. SafeBondAdder
Safely adds bonds while respecting valence constraints.

```python
from valence import SafeBondAdder

# Check if bond can be added
can_add, reason = SafeBondAdder.can_add_bond(
    mol,
    atom_idx1=0,
    atom_idx2=1,
    bond_type=1  # single bond
)

# Safely add bond
success, msg = SafeBondAdder.add_bond_safe(
    emol,           # EditableMol
    mol,            # Original mol for validation
    atom_idx1=0,
    atom_idx2=1,
    bond_type=Chem.BondType.SINGLE
)
```

### 3. ValenceReport
Generate summary statistics about valence.

```python
from valence import ValenceReport

# Get summary dict
summary = ValenceReport.get_summary(mol)
# {'total_atoms': 7, 'valid_atoms': 7, 'invalid_atoms': 0,
#  'element_counts': {'C': 6, 'O': 1}, 'all_valid': True}

# Print summary
ValenceReport.print_summary(mol)
```

## Usage in Biochar Generation

The valence system is automatically integrated:

### 1. Oxygen Assignment
```python
# Only adds -OH groups to aromatic carbons with available valence
aromatic_carbons = [
    atom.GetIdx() for atom in mol.GetAtoms()
    if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
]

for c_idx in aromatic_carbons:
    # Check if carbon has available valence
    val_info = ValenceValidator.get_valence_info(mol, c_idx)
    if val_info.available_valence >= 1:  # Need 1 bond for C-O
        # Add O-H group
```

### 2. Validation During Generation
```python
# After adding atoms/bonds
is_valid, errors = ValenceValidator.validate_molecule(new_mol)
if not is_valid:
    # Log errors but continue (errors caught in validation stage)
```

### 3. Final Validation
```python
# Chemical feasibility validator uses valence checking
errors = ChemicalFeasibilityValidator._check_valences(mol)
```

## Example Output

### Valence Report
```
================================================================================
VALENCE VALIDATION REPORT
================================================================================

Idx  Sym Min  Max  Bonds  Need  Avail Valid  Charge
--------------------------------------------------------------------------------
0    C   4    4    4      0     0     ✓ OK   0
1    C   4    4    4      0     0     ✓ OK   0
2    C   4    4    3      1     1     ✗ FAIL 0
3    O   2    2    2      0     0     ✓ OK   0
4    H   1    1    1      0     0     ✓ OK   0
--------------------------------------------------------------------------------
✓ All atoms have valid valence
================================================================================
```

**Column meanings:**
- `Min`: Minimum bonds required for valid valence
- `Max`: Maximum bonds allowed for valid valence
- `Bonds`: Current number of bonds
- `Need`: Additional bonds needed to reach minimum (if positive)
- `Avail`: Bonds still available before reaching maximum
- `Valid`: ✓ OK if `Min ≤ Bonds ≤ Max`, ✗ FAIL otherwise
- `Charge`: Formal charge on the atom

### Valence Summary
```
Valence Summary:
  Total atoms: 24
  Valid atoms: 24
  Invalid atoms: 0
  Status: ✓ VALID
  Elements: {'C': 16, 'H': 8}
```

## Valence Info Details

The `ValenceInfo` dataclass provides detailed information about an atom's bonding state:

```python
@dataclass
class ValenceInfo:
    atom_idx: int           # Index of atom in molecule
    symbol: str             # Element symbol (C, H, O, etc.)
    atomic_num: int         # Atomic number
    min_valence: int        # Minimum bonds required
    max_valence: int        # Maximum bonds allowed
    current_bonds: int      # Current number of bonds
    available_valence: int  # Bonds still available (max - current)
    needed_valence: int     # Bonds still needed to reach minimum
    is_valid: bool          # Whether current bonds are in valid range [min, max]
    is_saturated: bool      # Whether at maximum valence
    formal_charge: int      # Formal charge on atom
```

An atom is **valid** if: `min_valence ≤ current_bonds ≤ max_valence`

## Bond Types

| Bond Type | Order | Representation |
|-----------|-------|-----------------|
| SINGLE | 1 | C-C |
| DOUBLE | 2 | C=C |
| TRIPLE | 3 | C≡N |
| AROMATIC | 1.5 | c-c (benzene ring) |

## Extended Valences

Some elements can have multiple valence states in special cases (e.g., with formal charges):

- **Nitrogen**: [3, 3] standard, [4, 4] with +1 formal charge (e.g., NH₄⁺ - ammonium)
- **Sulfur**: [2, 2] standard (e.g., S-C), [4, 6] with higher oxidation states (e.g., SO₂, SO₃)

The system automatically adjusts valence ranges based on formal charge using the `get_valence_range(atomic_num, formal_charge)` function.

## Error Detection

### Common Valence Errors

The system detects and reports both **minimum and maximum** valence violations:

**Maximum valence exceeded:**
```
"Atom 0 (C): Valence 5 exceeds maximum 4"
```
Indicates carbon has too many bonds (over-saturated).

**Minimum valence not met:**
```
"Atom 2 (C): Valence 3 below minimum 4"
```
Indicates carbon doesn't have enough bonds (under-saturated).

**Insufficient valence for bond:**
```
"Atom 2 (O) insufficient valence: need 1, have 0"
```
Indicates oxygen cannot accept more bonds (already at maximum).

### Prevention

Before adding bonds, the system checks:

```python
can_add, reason = SafeBondAdder.can_add_bond(mol, atom1, atom2, bond_type)

if not can_add:
    print(f"Cannot add bond: {reason}")
    # Skip this bond and try alternative
```

## Best Practices

### 1. Check Before Adding
```python
# Always check valence before adding bonds
can_add, reason = SafeBondAdder.can_add_bond(mol, i, j, 1)
if can_add:
    SafeBondAdder.add_bond_safe(emol, mol, i, j, Chem.BondType.SINGLE)
```

### 2. Validate After Modifications
```python
# After significant changes, validate entire molecule
is_valid, errors = ValenceValidator.validate_molecule(mol)
if not is_valid:
    for error in errors:
        print(f"Valence error: {error}")
```

### 3. Report Generation
```python
# Generate detailed reports for debugging
ValenceValidator.print_valence_report(mol)
ValenceReport.print_summary(mol)
```

## Testing

The system includes comprehensive tests:

```bash
python3 << 'EOF'
from valence import ValenceValidator

# Test with simple molecules
mol_methanol = Chem.MolFromSmiles("CO")
is_valid, errors = ValenceValidator.validate_molecule(mol_methanol)
assert is_valid and len(errors) == 0

# Test with aromatic molecule
mol_phenol = Chem.MolFromSmiles("c1ccccc1O")
is_valid, errors = ValenceValidator.validate_molecule(mol_phenol)
assert is_valid and len(errors) == 0

print("✓ All valence tests passed")
EOF
```

## Integration with Biochar Generation

The valence system is integrated at multiple points:

1. **Skeleton Generation**: Carbon skeleton already respects valences (RDKit)
2. **Oxygen Assignment**: Checks available valence before adding -OH
3. **Hydrogen Assignment**: RDKit's AddHs respects valences
4. **Validation Stage**: Final check catches any issues

## Performance Impact

- **Minimal**: Valence checks are O(atoms) complexity
- **Fast**: Simple integer comparisons
- **Non-blocking**: Invalid bonds are skipped, not errors

## Future Enhancements

- [ ] Support for formal charges during structure generation
- [ ] Automatic charge balancing
- [ ] Extended valence states for transition metals
- [ ] Custom valence rules per element

## References

- Valence: https://en.wikipedia.org/wiki/Valence_(chemistry)
- GROMACS Bond Types: https://manual.gromacs.org/documentation/current/reference-manual/file-formats/gromacs-file-formats.html
- RDKit Valence: https://www.rdkit.org/

## Troubleshooting

### "Valence exceeds maximum"
**Cause**: Atom has too many bonds
**Solution**: Check your bond creation logic, ensure bonds aren't duplicated

### "Insufficient valence for bond"
**Cause**: Atom already fully bonded
**Solution**: Choose different atom or reduce bond order

### "Bond already exists"
**Cause**: Trying to add duplicate bond
**Solution**: Check if bond exists before adding: `mol.GetBondBetweenAtoms(i, j)`

---

All biochar structures generated by the simulator pass valence validation automatically.
