"""
Valence Validation System

Ensures all atoms maintain proper valence (maximum bonds per atom).
Based on standard chemistry valence rules and periodic table properties.
"""

from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import BondType


@dataclass
class ValenceInfo:
    """Information about an atom's valence."""

    atom_idx: int
    symbol: str
    atomic_num: int
    min_valence: int           # Minimum required bonds
    max_valence: int           # Maximum allowed bonds
    current_bonds: int         # Current number of bonds
    available_valence: int     # Max bonds still available
    needed_valence: int        # Min bonds still needed
    is_valid: bool             # Whether current bonds are in valid range
    formal_charge: int         # Formal charge on atom
    is_saturated: bool         # Whether fully saturated (at max valence)


# Standard valence ranges for common elements
# Format: (min_valence, max_valence, element_name)
STANDARD_VALENCES = {
    1: (1, 1, "H"),        # Hydrogen: exactly 1 bond
    6: (4, 4, "C"),        # Carbon: 4 bonds (can be lower in radicals, but we avoid those)
    7: (3, 3, "N"),        # Nitrogen: 3 bonds (neutral)
    8: (2, 2, "O"),        # Oxygen: 2 bonds
    9: (1, 1, "F"),        # Fluorine: 1 bond
    16: (2, 2, "S"),       # Sulfur: 2 bonds standard
    17: (1, 1, "Cl"),      # Chlorine: 1 bond
    35: (1, 1, "Br"),      # Bromine: 1 bond
}

# Extended valences for some elements in special cases
# Format: (min_valence, max_valence)
EXTENDED_VALENCES = {
    7: (4, 4),             # N can have 4 bonds (formal +1, ammonium-like)
    16: (4, 6),            # S can have 4 or 6 bonds (sulfone, sulfate)
}


def get_valence_range(atomic_num: int, formal_charge: int = 0) -> Tuple[int, int]:
    """
    Get valence range (min, max) for an atom.

    Args:
        atomic_num: Atomic number
        formal_charge: Formal charge on atom

    Returns:
        (min_valence, max_valence) tuple
    """
    if atomic_num not in STANDARD_VALENCES:
        # Default for unknown elements (allow 1-4 bonds)
        return (1, 4)

    min_val, max_val, _ = STANDARD_VALENCES[atomic_num]

    # Adjust for formal charge
    if formal_charge != 0 and atomic_num in EXTENDED_VALENCES:
        # Can use extended valence states
        ext_min, ext_max = EXTENDED_VALENCES[atomic_num]
        if formal_charge > 0:
            # Positive charge: use extended max
            return (ext_min, ext_max)
        elif formal_charge < 0:
            # Negative charge: minimum might decrease
            return (max(0, min_val - 1), max_val)

    return (min_val, max_val)


def get_max_valence(atomic_num: int, formal_charge: int = 0) -> int:
    """
    Get maximum valence for an atom. (Deprecated: use get_valence_range)

    Args:
        atomic_num: Atomic number
        formal_charge: Formal charge on atom

    Returns:
        Maximum number of bonds (valence)
    """
    _, max_val = get_valence_range(atomic_num, formal_charge)
    return max_val


def get_min_valence(atomic_num: int, formal_charge: int = 0) -> int:
    """
    Get minimum valence for an atom.

    Args:
        atomic_num: Atomic number
        formal_charge: Formal charge on atom

    Returns:
        Minimum number of bonds (valence)
    """
    min_val, _ = get_valence_range(atomic_num, formal_charge)
    return min_val


class ValenceValidator:
    """Validates and tracks atom valences in molecules."""

    @staticmethod
    def get_valence_info(mol: Chem.Mol, atom_idx: int) -> ValenceInfo:
        """
        Get detailed valence information for an atom.

        Args:
            mol: RDKit molecule
            atom_idx: Atom index

        Returns:
            ValenceInfo with current and valence range
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        atomic_num = atom.GetAtomicNum()
        symbol = atom.GetSymbol()
        formal_charge = atom.GetFormalCharge()

        # Count current bonds
        current_bonds = sum(
            bond.GetBondTypeAsDouble()
            for bond in atom.GetBonds()
        )
        current_bonds = int(current_bonds)

        # Get valence range
        min_val, max_val = get_valence_range(atomic_num, formal_charge)

        # Calculate available and needed
        available = max_val - current_bonds
        needed = max(0, min_val - current_bonds)

        # Check if valence is valid (within range)
        is_valid = min_val <= current_bonds <= max_val
        is_saturated = current_bonds >= max_val

        return ValenceInfo(
            atom_idx=atom_idx,
            symbol=symbol,
            atomic_num=atomic_num,
            min_valence=min_val,
            max_valence=max_val,
            current_bonds=current_bonds,
            available_valence=available,
            needed_valence=needed,
            is_valid=is_valid,
            formal_charge=formal_charge,
            is_saturated=is_saturated,
        )

    @staticmethod
    def validate_molecule(mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Validate all atoms in molecule have proper valence.

        Args:
            mol: RDKit molecule

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        for atom in mol.GetAtoms():
            val_info = ValenceValidator.get_valence_info(mol, atom.GetIdx())

            if not val_info.is_valid:
                # Report both minimum and maximum violations
                if val_info.current_bonds < val_info.min_valence:
                    errors.append(
                        f"Atom {val_info.atom_idx} ({val_info.symbol}): "
                        f"Valence {val_info.current_bonds} below minimum {val_info.min_valence}"
                    )
                elif val_info.current_bonds > val_info.max_valence:
                    errors.append(
                        f"Atom {val_info.atom_idx} ({val_info.symbol}): "
                        f"Valence {val_info.current_bonds} exceeds maximum {val_info.max_valence}"
                    )

        return len(errors) == 0, errors

    @staticmethod
    def get_all_valences(mol: Chem.Mol) -> List[ValenceInfo]:
        """
        Get valence information for all atoms.

        Args:
            mol: RDKit molecule

        Returns:
            List of ValenceInfo for all atoms
        """
        return [
            ValenceValidator.get_valence_info(mol, atom.GetIdx())
            for atom in mol.GetAtoms()
        ]

    @staticmethod
    def print_valence_report(mol: Chem.Mol):
        """Print detailed valence report for molecule."""
        print("\n" + "=" * 80)
        print("VALENCE VALIDATION REPORT")
        print("=" * 80)

        all_valid = True
        print(f"\n{'Idx':<4} {'Sym':<3} {'Min':<4} {'Max':<4} {'Bonds':<6} {'Need':<5} {'Avail':<5} {'Valid':<6} {'Charge':<6}")
        print("-" * 80)

        for val_info in ValenceValidator.get_all_valences(mol):
            symbol = val_info.symbol.ljust(3)
            valid_str = "✓ OK" if val_info.is_valid else "✗ FAIL"
            charge_str = f"{val_info.formal_charge:+d}" if val_info.formal_charge != 0 else "0"

            print(
                f"{val_info.atom_idx:<4} {symbol} {val_info.min_valence:<4} "
                f"{val_info.max_valence:<4} {val_info.current_bonds:<6} "
                f"{val_info.needed_valence:<5} {val_info.available_valence:<5} "
                f"{valid_str:<6} {charge_str:<6}"
            )

            if not val_info.is_valid:
                all_valid = False

        print("-" * 80)
        if all_valid:
            print("✓ All atoms have valid valence")
        else:
            print("✗ Some atoms have invalid valence")

        print("=" * 80 + "\n")


class SafeBondAdder:
    """Safely adds bonds while respecting valence constraints."""

    @staticmethod
    def can_add_bond(
        mol: Chem.Mol,
        atom_idx1: int,
        atom_idx2: int,
        bond_type: int = 1,  # 1 = single, 2 = double, 3 = triple
    ) -> Tuple[bool, str]:
        """
        Check if bond can be added without violating valence.

        Args:
            mol: RDKit molecule
            atom_idx1: First atom index
            atom_idx2: Second atom index
            bond_type: Bond type (1=single, 2=double, 3=triple)

        Returns:
            (can_add, reason)
        """
        val1 = ValenceValidator.get_valence_info(mol, atom_idx1)
        val2 = ValenceValidator.get_valence_info(mol, atom_idx2)

        # Check if atoms already bonded
        if mol.GetBondBetweenAtoms(atom_idx1, atom_idx2) is not None:
            return False, f"Bond already exists between atoms {atom_idx1} and {atom_idx2}"

        # Check maximum valence (atoms not over-saturated)
        if val1.available_valence < bond_type:
            return (
                False,
                f"Atom {atom_idx1} ({val1.symbol}) insufficient valence: "
                f"need {bond_type}, have {val1.available_valence}",
            )

        if val2.available_valence < bond_type:
            return (
                False,
                f"Atom {atom_idx2} ({val2.symbol}) insufficient valence: "
                f"need {bond_type}, have {val2.available_valence}",
            )

        # Check minimum valence (adding bond helps reach minimum if needed)
        # Note: We allow bonds that help satisfy minimum requirements
        # The final validation checks if minimums are met after all bonds added

        return True, "OK"

    @staticmethod
    def add_bond_safe(
        emol: Chem.EditableMol,
        mol: Chem.Mol,
        atom_idx1: int,
        atom_idx2: int,
        bond_type: Chem.BondType = Chem.BondType.SINGLE,
        verbose: bool = False,
    ) -> Tuple[bool, str]:
        """
        Safely add a bond if valence allows.

        Args:
            emol: RDKit EditableMol
            mol: Original RDKit molecule (for validation)
            atom_idx1: First atom index
            atom_idx2: Second atom index
            bond_type: Bond type (SINGLE, DOUBLE, TRIPLE)
            verbose: Print debug info

        Returns:
            (success, message)
        """
        bond_order = int(bond_type)

        can_add, reason = SafeBondAdder.can_add_bond(
            mol, atom_idx1, atom_idx2, bond_order
        )

        if not can_add:
            if verbose:
                print(f"  Cannot add bond {atom_idx1}-{atom_idx2}: {reason}")
            return False, reason

        try:
            emol.AddBond(atom_idx1, atom_idx2, bond_type)
            if verbose:
                print(f"  ✓ Added {bond_type.name} bond {atom_idx1}-{atom_idx2}")
            return True, "Bond added"
        except Exception as e:
            return False, f"Error adding bond: {str(e)}"

    @staticmethod
    def add_atom_safe(
        emol: Chem.EditableMol,
        atomic_num: int,
        formal_charge: int = 0,
    ) -> int:
        """
        Safely add an atom to molecule.

        Args:
            emol: RDKit EditableMol
            atomic_num: Atomic number
            formal_charge: Formal charge on atom

        Returns:
            New atom index
        """
        atom = Chem.Atom(atomic_num)
        if formal_charge != 0:
            atom.SetFormalCharge(formal_charge)

        atom_idx = emol.AddAtom(atom)
        return atom_idx


class ValenceReport:
    """Generate detailed valence reports."""

    @staticmethod
    def get_summary(mol: Chem.Mol) -> Dict:
        """
        Get summary of molecule valence.

        Returns:
            Dictionary with valence statistics
        """
        all_valences = ValenceValidator.get_all_valences(mol)

        valid_atoms = sum(1 for v in all_valences if v.is_valid)
        invalid_atoms = len(all_valences) - valid_atoms

        element_counts = {}
        for val_info in all_valences:
            if val_info.symbol not in element_counts:
                element_counts[val_info.symbol] = 0
            element_counts[val_info.symbol] += 1

        return {
            "total_atoms": len(all_valences),
            "valid_atoms": valid_atoms,
            "invalid_atoms": invalid_atoms,
            "element_counts": element_counts,
            "all_valid": invalid_atoms == 0,
        }

    @staticmethod
    def print_summary(mol: Chem.Mol):
        """Print summary report."""
        summary = ValenceReport.get_summary(mol)

        print(f"\nValence Summary:")
        print(f"  Total atoms: {summary['total_atoms']}")
        print(f"  Valid atoms: {summary['valid_atoms']}")
        print(f"  Invalid atoms: {summary['invalid_atoms']}")
        print(f"  Status: {'✓ VALID' if summary['all_valid'] else '✗ INVALID'}")
        print(f"  Elements: {summary['element_counts']}")
