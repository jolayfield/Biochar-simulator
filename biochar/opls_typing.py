"""
OPLS-AA Atom Typing and Charge Assignment

Assign OPLS-AA atom types and partial charges to atoms.
"""

import logging
from typing import Dict, Tuple, List, Optional
from dataclasses import dataclass

from rdkit import Chem

logger = logging.getLogger(__name__)

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

    @staticmethod
    def _is_carboxyl_carbon(atom: Chem.Atom) -> bool:
        """
        True when *atom* is the carbon of a carboxyl or carboxylate group.

        A carboxyl carbon carries two oxygens (C(=O)-OH or C(=O)-O-); a ketone
        or quinone carbon carries one, and an ether carbon's oxygen is not
        doubly bonded. Two oxygens on one carbon is therefore sufficient here --
        the generator builds no esters (the 'lactone' group falls back to
        phenolic and is never actually placed).

        This is what separates a carboxyl -OH from a phenolic -OH; without it
        both look like "an oxygen with one hydrogen".
        """
        if atom.GetAtomicNum() != 6:
            return False
        return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 8) >= 2

    @classmethod
    def _is_carboxylate_carbon(cls, atom: Chem.Atom) -> bool:
        """
        True when *atom* is the carbon of a *deprotonated* carboxylate.

        Detected by any of its oxygens carrying a negative formal charge. The
        two oxygens of a carboxylate are equivalent by resonance, so both take
        the carboxylate type -- RDKit's C(=O)[O-] Kekule form assigns the charge
        to only one of them, and typing the other from its double bond alone
        would leave one oxygen on the neutral carboxylic-acid type.
        """
        if not cls._is_carboxyl_carbon(atom):
            return False
        return any(
            n.GetAtomicNum() == 8 and n.GetFormalCharge() < 0
            for n in atom.GetNeighbors()
        )

    def _determine_atom_type(self, mol: Chem.Mol, atom: Chem.Atom) -> str:
        """
        Determine OPLS atom type for a single atom.

        Logic:
        - Aromatic C -> CA, aromatic H -> HA
        - Aliphatic sp3 C -> CT, H on sp3 C -> HC
        - Carbonyl O -> OC, hydroxyl O -> OH, ether O -> OS
        - Formal charge selects the ionized types: an anionic O is carboxylate
          (O2M) or phenolate (OM) by its neighbour; a cationic N is graphitic
          (NGR), pyridinium (NPYP), or anilinium (NAP) by ring membership and
          hydrogen count.
        """
        atomic_num = atom.GetAtomicNum()
        is_aromatic = atom.GetIsAromatic()
        formal_charge = atom.GetFormalCharge()

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
                if neighbor_type == "CA":
                    return "HA"
                elif neighbor_type == "OH":
                    return "HO"
                elif neighbor_type == "OH2":
                    return "HO2"
                elif neighbor_type == "NA":
                    return "HNA"
                elif neighbor_type == "SH_":
                    return "HSH"
                elif neighbor_type == "NPR":
                    return "HNPR"
                elif neighbor_type == "NPYP":
                    return "HPYP"
                elif neighbor_type == "NAP":
                    return "HNAP"
                else:
                    return "HC"
            return "HC"

        # Oxygen
        elif atomic_num == 8:
            # Check bonding environment
            num_bonds = len(atom.GetBonds())
            neighbors = list(atom.GetNeighbors())

            # Deprotonated oxygen: carboxylate or phenolate, by its neighbour.
            if formal_charge < 0:
                heavy = [n for n in neighbors if n.GetAtomicNum() != 1]
                if heavy and self._is_carboxyl_carbon(heavy[0]):
                    return "O2M"  # Carboxylate  Ar-COO-
                return "OM"       # Phenolate    Ar-O-

            if num_bonds == 1:
                # Terminal oxygen: a carbonyl of some kind.
                neighbor = neighbors[0]
                if neighbor.GetAtomicNum() == 6:
                    # Check if double bond
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        # Both oxygens of a carboxylate are equivalent by
                        # resonance, so the C=O of a deprotonated group takes
                        # the carboxylate type too -- not the neutral acid's.
                        if self._is_carboxylate_carbon(neighbor):
                            return "O2M"  # Carboxylate  Ar-COO-
                        # A carboxylic acid C=O is not a ketone C=O -- it has
                        # its own OPLS type and its own charge.
                        if self._is_carboxyl_carbon(neighbor):
                            return "O"   # Carboxylic acid C=O
                        return "OC"      # Ketone / aldehyde C=O
                    else:
                        return "O"  # Carboxylic acid O
            elif num_bonds == 2:
                # Connected to two atoms
                if any(n.GetAtomicNum() == 1 for n in neighbors):
                    # Hydroxyl -- but a carboxyl -OH and a phenolic -OH are
                    # different types. Both have exactly one H and one heavy
                    # neighbour, so only the neighbour tells them apart.
                    heavy = [n for n in neighbors if n.GetAtomicNum() != 1]
                    if heavy and self._is_carboxyl_carbon(heavy[0]):
                        return "OH2"  # Carboxylic acid -OH
                    return "OH"       # Phenolic / alcohol -OH
                else:
                    return "OS"  # Ether

            # Oxygen with 3+ bonds is not a species this pipeline builds.
            return "OH"

        # Nitrogen
        elif atomic_num == 7:
            neighbors = list(atom.GetNeighbors())
            h_count = sum(1 for n in neighbors if n.GetAtomicNum() == 1)
            heavy_neighbors = [n for n in neighbors if n.GetAtomicNum() != 1]
            has_aromatic_c = any(
                n.GetAtomicNum() == 6 and n.GetIsAromatic() for n in neighbors
            )

            # Ring-substituting nitrogen (pyridinic / pyrrolic / graphitic).
            # A ring N replaced INTO the skeleton is bonded only to ring carbons
            # (plus, for pyrrolic, one H).  Aniline N (pendant Ar-NH2) is not in
            # a ring, so it is excluded here.
            ring_info = mol.GetRingInfo()
            in_ring = ring_info.NumAtomRings(atom.GetIdx()) > 0

            # Cationic nitrogen. Graphitic N is quaternary and carries a formal
            # +1 by construction, so it must be matched before pyridinium --
            # both are cationic ring nitrogens, and only the hydrogen count and
            # heavy-neighbour count separate them.
            if formal_charge > 0:
                if in_ring:
                    if h_count == 0 and len(heavy_neighbors) >= 3:
                        return "NGR"   # Graphitic: quaternary, no H
                    if h_count >= 1:
                        return "NPYP"  # Pyridinium: protonated ring N
                elif h_count >= 2:
                    return "NAP"       # Anilinium Ar-NH3+

            if in_ring:
                ring_sizes = [
                    len(r) for r in ring_info.AtomRings() if atom.GetIdx() in r
                ]
                if h_count >= 1 and 5 in ring_sizes:
                    return "NPR"  # Pyrrolic N (5-ring, N-H)
                if h_count == 0 and 6 in ring_sizes:
                    if len(heavy_neighbors) >= 3:
                        return "NGR"  # Graphitic / quaternary N (interior)
                    return "NPY"  # Pyridinic N (edge 6-ring, no H)

            if has_aromatic_c and h_count >= 1:
                return "NA"  # Aniline-type aromatic primary amine (Ar-NH2)
            elif len(atom.GetBonds()) == 3:
                return "N"
            else:
                return "NT"

        # Sulfur
        elif atomic_num == 16:
            if formal_charge < 0:
                return "SM"   # Thiophenolate sulfur (Ar-S-)
            if any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()):
                return "SH_"  # Thiol sulfur (Ar-SH)
            else:
                return "SS"   # Thioether sulfur (Ar-S-Ar)

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
        Adjust charges so the molecule sums to its **total formal charge**.

        The target is the sum of the atoms' formal charges, not zero. A neutral
        molecule has a formal charge of 0, so this reduces exactly to the old
        sum-to-zero behaviour and leaves the non-pH path untouched.

        Targeting zero unconditionally -- as this did previously -- made an
        ionized structure impossible to express: any charge placed on a
        carboxylate was redistributed away until the molecule was neutral
        again, with no warning. It also quietly erased the formal +1 that
        graphitic nitrogen has carried since it was introduced.

        Neutralising a genuinely charged system is the job of `genion -neutral`
        at solvation time (see md_setup), not of the molecule definition.

        Args:
            mol: RDKit molecule
            charges: Dictionary of {atom_idx: charge}

        Returns:
            Equilibrated charges summing to the molecule's formal charge
        """
        target = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        residual = sum(charges.values()) - target

        if abs(residual) < 1e-6:
            return charges

        if abs(residual) > 0.01:
            logger.debug(
                "Charge residual before correction: %.4f e (target %+d e, "
                "distributing across %d heteroatoms)",
                residual,
                target,
                sum(1 for idx in charges if mol.GetAtomWithIdx(idx).GetAtomicNum() in [7, 8]),
            )

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
        scale = residual / len(adjustable_atoms)
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

        # The partial charges must sum to the molecule's formal charge.
        #
        # This replaces an `abs(total_charge) > 1.0` check that could never fire
        # -- charges were forced to zero before it ran, and it did nothing but
        # `pass` regardless. Comparing against the formal charge is a real
        # check: it catches an atom typed to the wrong ionization state, which
        # otherwise produces a topology that grompp accepts and simulates wrong.
        formal_charge = sum(a.GetFormalCharge() for a in self.mol.GetAtoms())
        total_charge = self.get_total_charge()
        if abs(total_charge - formal_charge) > 0.01:
            errors.append(
                f"Net charge {total_charge:+.4f} e does not match total formal "
                f"charge {formal_charge:+d} e (difference "
                f"{total_charge - formal_charge:+.4f} e)"
            )

        return len(errors) == 0, errors
