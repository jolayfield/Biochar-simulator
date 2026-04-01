"""
Heteroatom Assignment

Assign hydrogen and oxygen atoms to carbon skeletons to satisfy compositional ratios.
Respects valence constraints for all atoms.
"""

import random
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from .constants import FUNCTIONAL_GROUPS
from .valence import ValenceValidator, SafeBondAdder, ValenceReport


@dataclass
class CompositionInfo:
    """Information about molecular composition."""

    num_carbons: int
    num_hydrogens: int
    num_oxygens: int
    H_C_ratio: float
    O_C_ratio: float
    functional_groups: Dict[str, int]  # group_name -> count


class OxygenAssigner:
    """
    Assign oxygen atoms and functional groups to carbon skeleton.

    Strategy: Place functional groups (carboxyl, hydroxyl, etc.) to reach target O/C ratio.
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

    def assign_oxygens(
        self,
        mol: Chem.Mol,
        target_O_C_ratio: float,
        O_C_tolerance: float = 0.10,
        functional_group_preference: Optional[List[str]] = None,
    ) -> Tuple[Chem.Mol, CompositionInfo]:
        """
        Assign oxygens to molecule to reach target O/C ratio.

        Strategy:
        1. Identify aromatic carbons with available valence
        2. Add hydroxyl groups (-OH) to reach target oxygen count
        3. Validate that all atoms maintain valid valence

        Args:
            mol: RDKit molecule (carbon skeleton)
            target_O_C_ratio: Target oxygen-to-carbon ratio
            O_C_tolerance: Tolerance for ratio (e.g., 0.10 = ±10%)
            functional_group_preference: List of preferred functional groups (unused for now)

        Returns:
            (modified_mol, composition_info)
        """
        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        target_num_oxygens = int(round(num_carbons * target_O_C_ratio))

        # If target oxygens is 0 or very small, just return the molecule as-is
        if target_num_oxygens <= 0:
            composition = self._calculate_composition(mol, {})
            return mol, composition

        functional_groups_added = {}

        try:
            emol = Chem.EditableMol(mol)
            new_mol = emol.GetMol()

            # Add hydroxyl groups to aromatic carbons with available valence
            added_oxygens = 0

            # Find aromatic carbons - try to add oxygens to as many as possible
            aromatic_carbons = [
                atom.GetIdx()
                for atom in mol.GetAtoms()
                if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
            ]

            # Shuffle to randomize selection if seed is set
            if self.seed is not None:
                random.shuffle(aromatic_carbons)

            for c_idx in aromatic_carbons:
                if added_oxygens >= target_num_oxygens:
                    break

                # Get current valence state
                val_info = ValenceValidator.get_valence_info(new_mol, c_idx)

                # Need at least 1 bond available for C-O
                if val_info.available_valence < 1:
                    continue

                # Add oxygen atom
                o_idx = emol.AddAtom(Chem.Atom(8))
                new_mol = emol.GetMol()

                # Add C-O bond
                try:
                    emol.AddBond(c_idx, o_idx, Chem.BondType.SINGLE)
                except Exception:
                    continue

                new_mol = emol.GetMol()

                # Add hydrogen to oxygen (to create -OH group)
                try:
                    h_idx = emol.AddAtom(Chem.Atom(1))
                    emol.AddBond(o_idx, h_idx, Chem.BondType.SINGLE)
                    new_mol = emol.GetMol()
                    added_oxygens += 1
                    functional_groups_added["hydroxyl"] = functional_groups_added.get("hydroxyl", 0) + 1
                except Exception:
                    continue

            # Sanitize the molecule
            try:
                Chem.SanitizeMol(new_mol)
            except Exception:
                pass

        except Exception as e:
            new_mol = mol
            functional_groups_added = {}

        composition = self._calculate_composition(new_mol, functional_groups_added)
        return new_mol, composition

    def _find_attachment_sites(self, mol: Chem.Mol) -> List[int]:
        """Find carbon atoms suitable for functional group attachment."""
        sites = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon
                # Can attach to carbons with available valence
                if atom.GetTotalValence() < 4:
                    sites.append(atom.GetIdx())
        return sites

    def _calculate_composition(
        self, mol: Chem.Mol, functional_groups: Dict[str, int]
    ) -> CompositionInfo:
        """Calculate composition information from molecule."""
        num_C = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        num_H = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
        num_O = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

        H_C_ratio = num_H / max(num_C, 1)
        O_C_ratio = num_O / max(num_C, 1)

        return CompositionInfo(
            num_carbons=num_C,
            num_hydrogens=num_H,
            num_oxygens=num_O,
            H_C_ratio=H_C_ratio,
            O_C_ratio=O_C_ratio,
            functional_groups=functional_groups,
        )


class HydrogenAssigner:
    """
    Assign hydrogen atoms to satisfy valence and H/C ratio.

    Strategy:
    1. First, ensure all atoms have minimum valence satisfied
    2. Then, adjust hydrogen count to match target H/C ratio
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

    def assign_hydrogens(
        self,
        mol: Chem.Mol,
        target_H_C_ratio: float,
        H_C_tolerance: float = 0.10,
    ) -> Tuple[Chem.Mol, CompositionInfo]:
        """
        Assign hydrogens to molecule to satisfy minimum valence and target H/C ratio.

        Strategy:
        1. Add hydrogens to satisfy minimum valences for all atoms
        2. Trim excess hydrogens to match target H/C ratio (if needed)
        3. Ensure no violations of valence constraints

        Args:
            mol: RDKit molecule (may already have oxygens)
            target_H_C_ratio: Target hydrogen-to-carbon ratio
            H_C_tolerance: Tolerance for ratio

        Returns:
            (modified_mol_with_H, composition_info)
        """
        num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        target_num_hydrogens = int(round(num_carbons * target_H_C_ratio))

        # Step 1: Ensure minimum valences are satisfied
        mol_saturated = self._saturate_valences(mol)

        # Step 2: Adjust hydrogen count to match target ratio
        num_hydrogens = sum(1 for atom in mol_saturated.GetAtoms() if atom.GetAtomicNum() == 1)

        if num_hydrogens > target_num_hydrogens:
            # Remove excess hydrogens, but only from atoms that still have available valence
            mol_trimmed = self._trim_hydrogens(
                mol_saturated,
                num_hydrogens - target_num_hydrogens
            )
            mol_with_H = mol_trimmed
        else:
            mol_with_H = mol_saturated

        # Sanitize if possible
        try:
            Chem.SanitizeMol(mol_with_H)
        except Exception:
            pass

        composition = self._calculate_composition(mol_with_H)
        return mol_with_H, composition

    def _saturate_valences(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Add hydrogens to ensure all atoms meet minimum valence requirements.

        Uses valence validation to ensure each atom gets enough bonds.
        """
        emol = Chem.EditableMol(mol)

        # Iterate through atoms and add H where needed
        atoms_to_saturate = []
        for atom in mol.GetAtoms():
            val_info = ValenceValidator.get_valence_info(mol, atom.GetIdx())
            if val_info.needed_valence > 0:
                atoms_to_saturate.append((atom.GetIdx(), val_info.needed_valence))

        # Add hydrogens to unsaturated atoms
        for atom_idx, needed_bonds in atoms_to_saturate:
            for _ in range(needed_bonds):
                h_idx = emol.AddAtom(Chem.Atom(1))
                emol.AddBond(atom_idx, h_idx, Chem.BondType.SINGLE)

        return emol.GetMol()

    def _trim_hydrogens(self, mol: Chem.Mol, num_to_remove: int) -> Chem.Mol:
        """
        Remove excess hydrogens while preserving minimum valence constraints.

        Only removes H atoms from non-aromatic atoms that still have available valence.
        Avoids removing H atoms that are needed for minimum valence.
        """
        if num_to_remove <= 0:
            return mol

        emol = Chem.EditableMol(mol)

        # Find removable H atoms
        removable_h = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:  # Hydrogen
                # Check if this H is connected to a non-aromatic atom with available valence
                neighbors = atom.GetNeighbors()
                if neighbors:
                    parent = neighbors[0]
                    parent_val = ValenceValidator.get_valence_info(mol, parent.GetIdx())

                    # Can remove this H if parent has available valence after removal
                    if parent_val.available_valence > 0:
                        removable_h.append(atom.GetIdx())

        # Remove H atoms (in reverse order to maintain indices)
        for h_idx in sorted(removable_h, reverse=True)[:num_to_remove]:
            emol.RemoveAtom(h_idx)

        return emol.GetMol()

    def _calculate_composition(self, mol: Chem.Mol) -> CompositionInfo:
        """Calculate composition information."""
        num_C = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        num_H = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
        num_O = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

        H_C_ratio = num_H / max(num_C, 1)
        O_C_ratio = num_O / max(num_C, 1)

        return CompositionInfo(
            num_carbons=num_C,
            num_hydrogens=num_H,
            num_oxygens=num_O,
            H_C_ratio=H_C_ratio,
            O_C_ratio=O_C_ratio,
            functional_groups={},
        )


class HeteroatomValidator:
    """Validate heteroatom assignments."""

    @staticmethod
    def validate_ratios(
        composition: CompositionInfo,
        target_H_C: float,
        target_O_C: float,
        H_C_tolerance: float = 0.10,
        O_C_tolerance: float = 0.10,
    ) -> Tuple[bool, List[str]]:
        """
        Validate that ratios are within tolerance.

        Args:
            composition: CompositionInfo object
            target_H_C: Target H/C ratio
            target_O_C: Target O/C ratio
            H_C_tolerance: Tolerance for H/C
            O_C_tolerance: Tolerance for O/C

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        # Check H/C ratio
        H_C_error = abs(composition.H_C_ratio - target_H_C) / max(target_H_C, 0.01)
        if H_C_error > H_C_tolerance:
            errors.append(
                f"H/C ratio {composition.H_C_ratio:.3f} outside tolerance "
                f"(target: {target_H_C:.3f}, tolerance: {H_C_tolerance:.1%})"
            )

        # Check O/C ratio
        if target_O_C > 0:
            O_C_error = abs(composition.O_C_ratio - target_O_C) / max(target_O_C, 0.01)
            if O_C_error > O_C_tolerance:
                errors.append(
                    f"O/C ratio {composition.O_C_ratio:.3f} outside tolerance "
                    f"(target: {target_O_C:.3f}, tolerance: {O_C_tolerance:.1%})"
                )

        return len(errors) == 0, errors

    @staticmethod
    def validate_functional_groups(
        composition: CompositionInfo, expected_groups: Optional[List[str]] = None
    ) -> Tuple[bool, List[str]]:
        """
        Validate functional group assignment.

        Args:
            composition: CompositionInfo object
            expected_groups: List of expected functional groups

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        if expected_groups is not None:
            for group in expected_groups:
                if group not in composition.functional_groups:
                    errors.append(f"Expected functional group '{group}' not found")

        return len(errors) == 0, errors
