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


# Flags for sanitisation that skip SANITIZE_KEKULIZE and SANITIZE_SETAROMATICITY.
#
# SANITIZE_KEKULIZE:      fails on large graphene-like systems and partially
#                         corrupts aromatic bond types before raising.
#
# SANITIZE_SETAROMATICITY: uses RDKit's own aromaticity model (RDKIT, not MDL),
#                          which marks ether C-O bonds as AROMATIC when the bridge
#                          closes a ring that satisfies Hückel's rule (furan-like).
#                          We re-perceive aromaticity separately with
#                          SetAromaticity(mol, AROMATICITY_MDL) which correctly
#                          leaves C-O bonds as SINGLE while still marking ring C-C
#                          bonds as AROMATIC.
_SAFE_FLAGS = (
    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS
    | Chem.SanitizeFlags.SANITIZE_ADJUSTHS
)


def _safe_sanitize(mol: Chem.Mol) -> None:
    """Sanitize mol in-place without Kekulization, then re-perceive aromaticity."""
    try:
        Chem.SanitizeMol(mol, _SAFE_FLAGS)
    except Exception:
        pass
    try:
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    except Exception:
        pass


def _fix_heteroatom_bond_types(mol: Chem.Mol) -> Chem.Mol:
    """
    Correct bonds that involve heteroatoms (non-C) and were spuriously marked
    AROMATIC by aromaticity perception.

    This can happen when an ether oxygen bridges two aromatic edge carbons and
    the resulting large ring satisfies Hückel's rule (e.g., furan-like 5-ring).
    Oxygen, hydrogen, and carboxyl carbons outside the aromatic ring system should
    never have AROMATIC bonds.

    Returns a new mol with corrected bond types.
    """
    rwmol = Chem.RWMol(mol)
    changed = False
    for bond in rwmol.GetBonds():
        if bond.GetBondTypeAsDouble() == 1.5:  # currently AROMATIC
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only fix bonds where at least one atom is a heteroatom (not C)
            if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
                bond.SetBondType(Chem.BondType.SINGLE)
                if a1.GetAtomicNum() != 6:
                    a1.SetIsAromatic(False)
                if a2.GetAtomicNum() != 6:
                    a2.SetIsAromatic(False)
                changed = True
    return rwmol.GetMol() if changed else mol


# Groups that are not implementable on pure all-aromatic PAH edge carbons.
# They require sp2 C with ≥2 free valence (absent in our skeletons).
_FALLBACK_GROUPS: Dict[str, str] = {
    "carbonyl": "phenolic",
    "quinone":  "phenolic",
    "lactone":  "phenolic",
}


class OxygenAssigner:
    """
    Assign oxygen atoms and functional groups to carbon skeleton.

    Supported groups and oxygen count per placement:
        phenolic  — Ar-OH           (1 O)
        hydroxyl  — Ar-OH (= phenolic for pure PAH)  (1 O)
        carboxyl  — Ar-C(=O)(OH)    (2 O)
        ether     — Ar-O-Ar bridge  (1 O, needs 2 edge sites)
        carbonyl  — not supported for pure aromatic PAH; falls back to phenolic
        quinone   — not supported for pure aromatic PAH; falls back to phenolic
        lactone   — not supported for pure aromatic PAH; falls back to phenolic

    Usage — exact-count dict:
        assigner.assign_oxygens(mol, ...,
            functional_group_preference={"phenolic": 3, "carboxyl": 1})

    Usage — O_C_ratio-driven (default, phenolic only):
        assigner.assign_oxygens(mol, target_O_C_ratio=0.1)
    """

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            random.seed(seed)

    # ------------------------------------------------------------------ #
    #  Public interface                                                    #
    # ------------------------------------------------------------------ #

    def assign_oxygens(
        self,
        mol: Chem.Mol,
        target_O_C_ratio: float,
        O_C_tolerance: float = 0.10,
        functional_group_preference: Optional[Dict[str, int]] = None,
    ) -> Tuple[Chem.Mol, CompositionInfo]:
        """
        Assign oxygen-containing functional groups to the molecule.

        Args:
            mol: RDKit molecule (carbon skeleton, all-aromatic)
            target_O_C_ratio: Target O/C ratio (used only when
                functional_group_preference is None).
            O_C_tolerance: Tolerance for ratio check (informational only).
            functional_group_preference: Dict mapping group name → exact count,
                e.g. {"phenolic": 3, "carboxyl": 1}.
                If None, places phenolic groups to reach target_O_C_ratio.

        Returns:
            (modified_mol, composition_info)
        """
        num_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)

        # ── Mode 1: exact-count dict ───────────────────────────────────────
        if isinstance(functional_group_preference, dict) and functional_group_preference:
            groups_spec = self._validate_and_normalise_spec(functional_group_preference)
        # ── Mode 2: O_C_ratio-driven (legacy / default) ────────────────────
        else:
            target_num_oxygens = int(round(num_carbons * target_O_C_ratio))
            if target_num_oxygens <= 0:
                return mol, self._calculate_composition(mol, {})
            groups_spec = {"phenolic": target_num_oxygens}

        if not groups_spec:
            return mol, self._calculate_composition(mol, {})

        # ── Build edge-site pool ───────────────────────────────────────────
        site_pool = self._get_edge_sites(mol)
        if self.seed is not None:
            random.seed(self.seed)
        random.shuffle(site_pool)
        used_sites: Set[int] = set()

        emol = Chem.EditableMol(mol)
        new_mol = mol
        functional_groups_added: Dict[str, int] = {}

        # ── Place each group type the requested number of times ────────────
        for group, count in groups_spec.items():
            placed = 0
            while placed < count:
                free = [s for s in site_pool if s not in used_sites]
                ok, new_mol, emol, n_O, sites_used = self._place_group(
                    group, emol, new_mol, free
                )
                if not ok:
                    print(
                        f"Warning: Could not place all '{group}' groups "
                        f"({placed}/{count} placed — not enough edge sites)"
                    )
                    break
                placed += 1
                used_sites.update(sites_used)
                functional_groups_added[group] = (
                    functional_groups_added.get(group, 0) + 1
                )

        _safe_sanitize(new_mol)
        # Fix any heteroatom bonds spuriously marked AROMATIC (e.g., ether O
        # closing a ring that satisfies Hückel's rule through the PAH system).
        new_mol = _fix_heteroatom_bond_types(new_mol)
        return new_mol, self._calculate_composition(new_mol, functional_groups_added)

    # ------------------------------------------------------------------ #
    #  Private helpers                                                    #
    # ------------------------------------------------------------------ #

    def _validate_and_normalise_spec(
        self, spec: Dict[str, int]
    ) -> Dict[str, int]:
        """
        Validate group names and substitute unsupported ones with their fallback.
        Returns the cleaned spec dict.
        """
        valid_keys = set(FUNCTIONAL_GROUPS.keys())

        # Remove unknown keys
        unknown = [k for k in spec if k not in valid_keys]
        if unknown:
            print(f"Warning: Unknown functional groups ignored: {unknown}")
        cleaned: Dict[str, int] = {k: v for k, v in spec.items() if k in valid_keys}

        # Substitute unsupported groups (warn once each)
        warned: Set[str] = set()
        for g in list(cleaned):
            if g in _FALLBACK_GROUPS:
                fallback = _FALLBACK_GROUPS[g]
                if g not in warned:
                    print(
                        f"Warning: '{g}' is not supported for pure aromatic PAH "
                        f"(needs ≥2 free valence on one carbon); "
                        f"substituting with '{fallback}'"
                    )
                    warned.add(g)
                cleaned[fallback] = cleaned.get(fallback, 0) + cleaned.pop(g)

        return cleaned

    def _get_edge_sites(self, mol: Chem.Mol) -> List[int]:
        """Return indices of edge aromatic carbons with ≥1 free valence."""
        return [
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6
            and a.GetIsAromatic()
            and ValenceValidator.get_valence_info(mol, a.GetIdx()).available_valence >= 1
        ]

    def _place_group(
        self,
        group: str,
        emol: Chem.EditableMol,
        new_mol: Chem.Mol,
        free_sites: List[int],
    ) -> Tuple[bool, Chem.Mol, Chem.EditableMol, int, Set[int]]:
        """
        Attempt to attach one instance of *group* to the molecule.

        Returns:
            (success, updated_mol, updated_emol, oxygens_added, sites_consumed)
        """
        BT = Chem.BondType

        if group in ("phenolic", "hydroxyl"):
            # Ar-OH: single-bond O then H
            if not free_sites:
                return False, new_mol, emol, 0, set()
            c = free_sites[0]
            o = emol.AddAtom(Chem.Atom(8))
            h = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(c, o, BT.SINGLE)
            emol.AddBond(o, h, BT.SINGLE)
            return True, emol.GetMol(), emol, 1, {c}

        elif group == "carboxyl":
            # Ar-C(=O)(OH): one ring C → new sp2 C → =O and -OH
            if not free_sites:
                return False, new_mol, emol, 0, set()
            c_ring = free_sites[0]
            c_new  = emol.AddAtom(Chem.Atom(6))  # carboxyl carbon
            o1     = emol.AddAtom(Chem.Atom(8))   # C=O
            o2     = emol.AddAtom(Chem.Atom(8))   # C-OH
            h      = emol.AddAtom(Chem.Atom(1))
            emol.AddBond(c_ring, c_new, BT.SINGLE)
            emol.AddBond(c_new,  o1,   BT.DOUBLE)
            emol.AddBond(c_new,  o2,   BT.SINGLE)
            emol.AddBond(o2,     h,    BT.SINGLE)
            return True, emol.GetMol(), emol, 2, {c_ring}

        elif group == "ether":
            # Ar-O-Ar bridge: two NON-ADJACENT edge C atoms connected through one O.
            # Adjacent C atoms would create a 3-membered ring whose C-O bonds can be
            # miscategorised as aromatic by RDKit's aromaticity model.
            c1, c2 = None, None
            for i in range(len(free_sites)):
                for j in range(i + 1, len(free_sites)):
                    s1, s2 = free_sites[i], free_sites[j]
                    if new_mol.GetBondBetweenAtoms(s1, s2) is None:
                        c1, c2 = s1, s2
                        break
                if c1 is not None:
                    break
            if c1 is None:
                return False, new_mol, emol, 0, set()
            o = emol.AddAtom(Chem.Atom(8))
            emol.AddBond(c1, o, BT.SINGLE)
            emol.AddBond(c2, o, BT.SINGLE)
            return True, emol.GetMol(), emol, 1, {c1, c2}

        # Unknown / unsupported (should have been filtered before reaching here)
        return False, new_mol, emol, 0, set()

    def _find_attachment_sites(self, mol: Chem.Mol) -> List[int]:
        """Find carbon atoms suitable for functional group attachment."""
        sites = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
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
            mol_trimmed = self._trim_hydrogens(
                mol_saturated,
                num_hydrogens - target_num_hydrogens
            )
            mol_with_H = mol_trimmed
        else:
            mol_with_H = mol_saturated

        _safe_sanitize(mol_with_H)
        # After aromaticity re-perception, ether oxygens that bridge two ring
        # carbons can be spuriously marked AROMATIC.  Fix them back to SINGLE.
        mol_with_H = _fix_heteroatom_bond_types(mol_with_H)

        composition = self._calculate_composition(mol_with_H)
        return mol_with_H, composition

    def _saturate_valences(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Add hydrogens to ensure all atoms meet minimum valence requirements.

        Uses valence validation to ensure each atom gets enough bonds.
        """
        emol = Chem.EditableMol(mol)

        atoms_to_saturate = []
        for atom in mol.GetAtoms():
            val_info = ValenceValidator.get_valence_info(mol, atom.GetIdx())
            if val_info.needed_valence > 0:
                atoms_to_saturate.append((atom.GetIdx(), val_info.needed_valence))

        for atom_idx, needed_bonds in atoms_to_saturate:
            for _ in range(needed_bonds):
                h_idx = emol.AddAtom(Chem.Atom(1))
                emol.AddBond(atom_idx, h_idx, Chem.BondType.SINGLE)

        return emol.GetMol()

    def _trim_hydrogens(self, mol: Chem.Mol, num_to_remove: int) -> Chem.Mol:
        """
        Remove excess hydrogens while preserving minimum valence constraints.

        Only removes H atoms from atoms that still have available valence.
        """
        if num_to_remove <= 0:
            return mol

        emol = Chem.EditableMol(mol)

        removable_h = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                neighbors = atom.GetNeighbors()
                if neighbors:
                    parent = neighbors[0]
                    parent_val = ValenceValidator.get_valence_info(mol, parent.GetIdx())
                    if parent_val.available_valence > 0:
                        removable_h.append(atom.GetIdx())

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

        H_C_error = abs(composition.H_C_ratio - target_H_C) / max(target_H_C, 0.01)
        if H_C_error > H_C_tolerance:
            errors.append(
                f"H/C ratio {composition.H_C_ratio:.3f} outside tolerance "
                f"(target: {target_H_C:.3f}, tolerance: {H_C_tolerance:.1%})"
            )

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
