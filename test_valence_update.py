#!/usr/bin/env python3
"""
Test script to verify the updated valence validation system enforces both
minimum and maximum valences correctly.
"""

from rdkit import Chem
from src.valence import (
    ValenceValidator, SafeBondAdder, ValenceReport, ValenceInfo,
    get_valence_range, get_min_valence, get_max_valence
)


def test_valence_ranges():
    """Test that valence ranges are correct for common elements."""
    print("=" * 70)
    print("TEST 1: Valence Ranges for Common Elements")
    print("=" * 70)

    test_cases = [
        (1, "H", (1, 1)),      # Hydrogen: exactly 1
        (6, "C", (4, 4)),      # Carbon: 4 bonds
        (8, "O", (2, 2)),      # Oxygen: 2 bonds
        (7, "N", (3, 3)),      # Nitrogen: 3 bonds (neutral)
        (16, "S", (2, 2)),     # Sulfur: 2 bonds standard
        (17, "Cl", (1, 1)),    # Chlorine: 1 bond
    ]

    for atomic_num, symbol, expected_range in test_cases:
        min_val, max_val = get_valence_range(atomic_num)
        status = "✓" if (min_val, max_val) == expected_range else "✗"
        print(f"{status} {symbol:2d} (Z={atomic_num:2d}): min={min_val}, max={max_val} "
              f"(expected {expected_range})")

    print()


def test_simple_molecules():
    """Test valence validation on well-known molecules."""
    print("=" * 70)
    print("TEST 2: Valence Validation on Simple Molecules")
    print("=" * 70)

    molecules = {
        "H2": "CC",                           # Ethane (C2H6)
        "Methanol": "CO",                     # Methanol (CH3OH)
        "Ethanol": "CCO",                     # Ethanol (C2H5OH)
        "Benzene": "c1ccccc1",               # Benzene (C6H6)
        "Phenol": "c1ccccc1O",               # Phenol (C6H5OH)
        "Formaldehyde": "C=O",               # Formaldehyde (CH2O)
    }

    for name, smiles in molecules.items():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"✗ {name:20s}: Could not parse SMILES '{smiles}'")
                continue

            # Add hydrogens
            mol = Chem.AddHs(mol)

            # Validate
            is_valid, errors = ValenceValidator.validate_molecule(mol)

            status = "✓" if is_valid else "✗"
            error_msg = f" - {errors[0]}" if errors else ""
            print(f"{status} {name:20s}: {smiles:20s} - {'Valid' if is_valid else 'Invalid'}{error_msg}")

        except Exception as e:
            print(f"✗ {name:20s}: Exception - {str(e)}")

    print()


def test_valence_info():
    """Test ValenceInfo dataclass for detailed valence information."""
    print("=" * 70)
    print("TEST 3: Detailed Valence Information (ValenceInfo)")
    print("=" * 70)

    # Create a simple molecule
    mol = Chem.MolFromSmiles("c1ccccc1O")  # Phenol
    mol = Chem.AddHs(mol)

    print(f"\nPhenol (C6H5OH) - {mol.GetNumAtoms()} atoms:")
    print("-" * 70)

    for atom in mol.GetAtoms():
        val_info = ValenceValidator.get_valence_info(mol, atom.GetIdx())

        # Format output
        status = "✓" if val_info.is_valid else "✗"
        saturated = "SAT" if val_info.is_saturated else "UNSAT"

        print(f"  Atom {val_info.atom_idx}: {val_info.symbol:2s} | "
              f"Min={val_info.min_valence} Max={val_info.max_valence} | "
              f"Bonds={val_info.current_bonds} | "
              f"Need={val_info.needed_valence} Avail={val_info.available_valence} | "
              f"{status} {saturated}")

    print()


def test_valence_report():
    """Test the detailed valence report printing."""
    print("=" * 70)
    print("TEST 4: Valence Report Printing")
    print("=" * 70)

    # Create a biochar-like structure
    mol = Chem.MolFromSmiles("c1cc(O)ccc1c2ccccc2")  # Hydroxylated biphenyl
    if mol is None:
        print("Could not create test molecule")
        return

    mol = Chem.AddHs(mol)

    # Print detailed report
    ValenceValidator.print_valence_report(mol)

    # Print summary
    summary = ValenceReport.get_summary(mol)
    print(f"Summary:")
    print(f"  Total atoms: {summary['total_atoms']}")
    print(f"  Valid atoms: {summary['valid_atoms']}")
    print(f"  Invalid atoms: {summary['invalid_atoms']}")
    print(f"  All valid: {summary['all_valid']}")
    print(f"  Elements: {summary['element_counts']}")
    print()


def test_safe_bond_adder():
    """Test SafeBondAdder respects valence constraints."""
    print("=" * 70)
    print("TEST 5: SafeBondAdder - Respecting Valence Constraints")
    print("=" * 70)

    # Create two separate atoms
    emol = Chem.EditableMol(Chem.Mol())

    # Add atoms: C and C
    c1_idx = emol.AddAtom(Chem.Atom(6))
    c2_idx = emol.AddAtom(Chem.Atom(6))

    mol = emol.GetMol()

    print(f"Created 2 carbon atoms: C0 and C1")

    # Test 1: Can add single bond (both have plenty of valence)
    can_add, reason = SafeBondAdder.can_add_bond(mol, c1_idx, c2_idx, 1)
    print(f"  Single bond C-C: {'✓ Can add' if can_add else '✗ Cannot add'} ({reason})")

    # Add the bond
    if can_add:
        emol.AddBond(c1_idx, c2_idx, Chem.BondType.SINGLE)

    mol = emol.GetMol()

    # Now we have C-C, get valence info
    val_c1 = ValenceValidator.get_valence_info(mol, c1_idx)
    val_c2 = ValenceValidator.get_valence_info(mol, c2_idx)

    print(f"\nAfter adding C-C bond:")
    print(f"  C0: {val_c1.current_bonds} bonds (avail: {val_c1.available_valence})")
    print(f"  C1: {val_c2.current_bonds} bonds (avail: {val_c2.available_valence})")

    print()


def test_extended_valences():
    """Test extended valences (e.g., S with formal charge)."""
    print("=" * 70)
    print("TEST 6: Extended Valences (Formal Charges)")
    print("=" * 70)

    # Test N with positive charge (ammonium-like)
    min_val, max_val = get_valence_range(7, formal_charge=1)  # N+
    print(f"Nitrogen with +1 charge: min={min_val}, max={max_val}")

    # Test S with no charge vs with charge
    min_val_neutral, max_val_neutral = get_valence_range(16, formal_charge=0)
    min_val_charged, max_val_charged = get_valence_range(16, formal_charge=2)
    print(f"Sulfur neutral: min={min_val_neutral}, max={max_val_neutral}")
    print(f"Sulfur +2: min={min_val_charged}, max={max_val_charged}")

    print()


def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("VALENCE VALIDATION SYSTEM - COMPREHENSIVE TEST SUITE")
    print("=" * 70 + "\n")

    try:
        test_valence_ranges()
        test_simple_molecules()
        test_valence_info()
        test_valence_report()
        test_safe_bond_adder()
        test_extended_valences()

        print("=" * 70)
        print("✓ ALL TESTS COMPLETED SUCCESSFULLY")
        print("=" * 70)

    except Exception as e:
        print(f"\n✗ ERROR during testing: {str(e)}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
