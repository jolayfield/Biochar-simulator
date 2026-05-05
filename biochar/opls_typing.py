"""
OPLS-AA Atom Typing and Charge Assignment

Assign OPLS-AA atom types and partial charges to atoms.
"""

from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass

from rdkit import Chem

from .constants import OPLS_ATOM_TYPES, OPLS_BOND_PARAMS, OPLS_ANGLE_PARAMS


@dataclass
class AtomProperty:
    """Atom type and charge information."""

    atom_idx: int
    symbol: str
    opls_type: str
    charge: float
    mass: float
    aromatic: bool


class AtomTyper:
    """Assign OPLS-AA atom types based on chemical environment."""

    def __init__(self):
        self.opls_types = OPLS_ATOM_TYPES
        self.bond_params = OPLS_BOND_PARAMS
        self.angle_params = OPLS_ANGLE_PARAMS

    def assign_atom_types(self, mol: Chem.Mol) -> Dict[int, str]:
        """
        Assign OPLS-AA atom types to all atoms.

        Args:
            mol: RDKit molecule

        Returns:
            Dictionary of {atom_idx: opls_type}
        """
        atom_types = {}

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            atom_type = self._determine_atom_type(mol, atom)
            atom_types[idx] = atom_type

        return atom_types

    def _determine_atom_type(self, mol: Chem.Mol, atom: Chem.Atom) -> str:
        """
        Determine OPLS atom type for a single atom.

        Logic:
        - Aromatic C -> CA, aromatic H -> HA
        - Aliphatic sp3 C -> CT, H on sp3 C -> HC
        - Carbonyl O -> OC, hydroxyl O -> OH, ether O -> OS
        """
        atomic_num = atom.GetAtomicNum()
        is_aromatic = atom.GetIsAromatic()

        # Carbon
        if atomic_num == 6:
            # For large fused-ring systems RDKit Kekulization can fail, leaving
            # is_aromatic=False on all atoms.  Fall back to ring membership +
            # degree-3 connectivity (2 ring C neighbours + 1 H or O, or 3 ring C
            # neighbours for interior junction atoms) as a reliable proxy.
            ring_info = mol.GetRingInfo()
            in_ring = ring_info.NumAtomRings(atom.GetIdx()) > 0
            if is_aromatic or (in_ring and atom.GetDegree() == 3):
                return "CA"
            else:
                # Check if connected to C=O
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            return "C"  # Carboxylic acid carbon
                return "CT"  # Aliphatic carbon

        # Hydrogen
        elif atomic_num == 1:
            # Check what it's bonded to
            neighbors = atom.GetNeighbors()
            if neighbors:
                neighbor = neighbors[0]
                neighbor_type = self._determine_atom_type(mol, neighbor)
                if neighbor_type in ["CA"]:
                    return "HA"
                elif neighbor_type in ["OH"]:
                    return "HO"
                elif neighbor_type in ["OH2"]:
                    return "HO2"
                else:
                    return "HC"
            return "HC"

        # Oxygen
        elif atomic_num == 8:
            # Check bonding environment
            num_bonds = len(atom.GetBonds())
            neighbors = list(atom.GetNeighbors())

            if num_bonds == 1:
                # Terminal oxygen (likely carbonyl)
                neighbor = neighbors[0]
                if neighbor.GetAtomicNum() == 6:
                    # Check if double bond
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        return "OC"  # Carbonyl
                    else:
                        return "O"  # Carboxylic acid O
            elif num_bonds == 2:
                # Connected to two atoms
                if any(n.GetAtomicNum() == 1 for n in neighbors):
                    return "OH"  # Hydroxyl
                else:
                    return "OS"  # Ether
            else:
                # Special case for carboxylic acid OH
                return "OH2"

        # Nitrogen (placeholder)
        elif atomic_num == 7:
            num_bonds = len(atom.GetBonds())
            if num_bonds == 3:
                return "N"
            else:
                return "NT"

        # Sulfur (placeholder)
        elif atomic_num == 16:
            if any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()):
                return "SH"
            else:
                return "S"

        # Default
        return f"X{atomic_num}"

    def get_bond_parameters(
        self, atom_type_1: str, atom_type_2: str
    ) -> Optional[Tuple[float, float]]:
        """
        Get bond parameters for atom type pair.

        Returns:
            (force_constant_kcal_mol_A2, equilibrium_length_A)
        """
        key = tuple(sorted([atom_type_1, atom_type_2]))
        return self.bond_params.get(key)

    def get_angle_parameters(
        self, atom_type_1: str, atom_type_2: str, atom_type_3: str
    ) -> Optional[Tuple[float, float]]:
        """
        Get angle parameters for atom type triple.

        Returns:
            (force_constant_kcal_mol_rad2, equilibrium_angle_deg)
        """
        return self.angle_params.get((atom_type_1, atom_type_2, atom_type_3))


class ChargeAssigner:
    """
    Assign partial charges to atoms using OPLS-AA parameters.

    Strategy: Use predefined OPLS charges, with adjustments for unusual groups.
    """

    def __init__(self):
        self.opls_types = OPLS_ATOM_TYPES

    def assign_charges(
        self, mol: Chem.Mol, atom_types: Dict[int, str]
    ) -> Dict[int, float]:
        """
        Assign partial charges to all atoms.

        Args:
            mol: RDKit molecule
            atom_types: Dictionary of {atom_idx: opls_type}

        Returns:
            Dictionary of {atom_idx: charge}
        """
        charges = {}

        for idx, opls_type in atom_types.items():
            if opls_type in self.opls_types:
                # Use predefined charge
                charges[idx] = self.opls_types[opls_type][2]
            else:
                # Estimate charge based on atom type
                charges[idx] = self._estimate_charge(mol, idx, opls_type)

        # Apply charge equilibration to ensure overall neutrality
        charges = self._equilibrate_charges(mol, charges)

        return charges

    def _estimate_charge(self, mol: Chem.Mol, atom_idx: int, opls_type: str) -> float:
        """
        Estimate charge for atom if not in standard OPLS table.

        Args:
            mol: RDKit molecule
            atom_idx: Atom index
            opls_type: OPLS atom type string

        Returns:
            Estimated partial charge
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        atomic_num = atom.GetAtomicNum()

        # Default charge based on element and hybridization
        if atomic_num == 6:  # Carbon
            return -0.15  # Typical carbon charge
        elif atomic_num == 8:  # Oxygen
            return -0.50  # Typical oxygen charge
        elif atomic_num == 1:  # Hydrogen
            return 0.10  # Typical hydrogen charge
        elif atomic_num == 7:  # Nitrogen
            return -0.30
        elif atomic_num == 16:  # Sulfur
            return -0.20
        else:
            return 0.0

    def _equilibrate_charges(
        self, mol: Chem.Mol, charges: Dict[int, float]
    ) -> Dict[int, float]:
        """
        Adjust charges to ensure overall neutrality.

        Scales all charges proportionally to sum to zero.

        Args:
            mol: RDKit molecule
            charges: Dictionary of {atom_idx: charge}

        Returns:
            Equilibrated charges
        """
        total_charge = sum(charges.values())

        if abs(total_charge) < 1e-6:
            return charges

        # Find polar atoms to adjust
        # Prefer to adjust oxygen and nitrogen
        adjustable_atoms = [
            idx
            for idx, charge in charges.items()
            if mol.GetAtomWithIdx(idx).GetAtomicNum() in [8, 7]
        ]

        if not adjustable_atoms:
            adjustable_atoms = list(charges.keys())

        # Scale adjustment
        scale = total_charge / len(adjustable_atoms)
        adjusted_charges = charges.copy()

        for idx in adjustable_atoms:
            adjusted_charges[idx] -= scale

        return adjusted_charges


class OPLSPropertyTable:
    """
    Generate OPLS property table for GROMACS topology.

    Contains atom type, mass, charge information for all atoms.
    """

    def __init__(self, mol: Chem.Mol, atom_types: Dict[int, str], charges: Dict[int, float]):
        self.mol = mol
        self.atom_types = atom_types
        self.charges = charges
        self.properties: List[AtomProperty] = []
        self._build_properties()

    def _build_properties(self):
        """Build atom property table."""
        for atom in self.mol.GetAtoms():
            idx = atom.GetIdx()
            opls_type = self.atom_types.get(idx, "CT")

            if opls_type in OPLS_ATOM_TYPES:
                desc, mass, default_charge = OPLS_ATOM_TYPES[opls_type]
            else:
                mass = atom.GetMass()
                default_charge = 0.0

            charge = self.charges.get(idx, default_charge)

            prop = AtomProperty(
                atom_idx=idx,
                symbol=atom.GetSymbol(),
                opls_type=opls_type,
                charge=charge,
                mass=mass,
                aromatic=atom.GetIsAromatic(),
            )
            self.properties.append(prop)

    def get_properties(self) -> List[AtomProperty]:
        """Get atom properties."""
        return self.properties

    def get_total_charge(self) -> float:
        """Calculate total molecular charge."""
        return sum(p.charge for p in self.properties)

    def get_mass(self) -> float:
        """Calculate total molecular mass."""
        return sum(p.mass for p in self.properties)

    def validate(self) -> Tuple[bool, List[str]]:
        """
        Validate OPLS properties.

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        # Check for unassigned atom types
        for prop in self.properties:
            if prop.opls_type.startswith("X"):
                errors.append(f"Atom {prop.atom_idx} has unrecognized type")

        # Check for extreme charges
        for prop in self.properties:
            if abs(prop.charge) > 2.0:
                errors.append(
                    f"Atom {prop.atom_idx} has extreme charge: {prop.charge:.2f}"
                )

        # Check total charge is reasonable
        total_charge = self.get_total_charge()
        if abs(total_charge) > 1.0:
            # Warning, not error
            pass

        return len(errors) == 0, errors
