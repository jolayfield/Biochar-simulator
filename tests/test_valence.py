"""
Tests for biochar/valence.py — valence ranges, validator, report, safe bond adder.
"""

import pytest
from rdkit import Chem

from biochar.valence import (
    get_valence_range,
    get_max_valence,
    get_min_valence,
    ValenceValidator,
    ValenceReport,
    SafeBondAdder,
    STANDARD_VALENCES,
)


class TestGetValenceRange:
    """Test get_valence_range for standard elements."""

    def test_hydrogen_is_exactly_one(self):
        lo, hi = get_valence_range(1)
        assert lo == 1
        assert hi == 1

    def test_carbon_is_exactly_four(self):
        lo, hi = get_valence_range(6)
        assert lo == 4
        assert hi == 4

    def test_oxygen_is_exactly_two(self):
        lo, hi = get_valence_range(8)
        assert lo == 2
        assert hi == 2

    def test_nitrogen_neutral_is_three(self):
        lo, hi = get_valence_range(7)
        assert lo == 3
        assert hi == 3

    def test_nitrogen_positive_charge_allows_four(self):
        lo, hi = get_valence_range(7, formal_charge=1)
        assert hi == 4

    def test_unknown_element_defaults_to_one_four(self):
        lo, hi = get_valence_range(99)  # Einsteinium — not in table
        assert lo >= 1
        assert hi >= 1

    def test_fluorine_is_exactly_one(self):
        lo, hi = get_valence_range(9)
        assert lo == 1
        assert hi == 1

    def test_chlorine_is_exactly_one(self):
        lo, hi = get_valence_range(17)
        assert lo == 1
        assert hi == 1


class TestGetMaxMinValence:
    def test_max_valence_carbon(self):
        assert get_max_valence(6) == 4

    def test_min_valence_carbon(self):
        assert get_min_valence(6) == 4

    def test_max_valence_hydrogen(self):
        assert get_max_valence(1) == 1

    def test_min_valence_oxygen(self):
        assert get_min_valence(8) == 2


class TestValenceValidatorGetInfo:
    def test_aromatic_carbon_in_benzene(self):
        # Need explicit H so aromatic C counts 2×1.5 + 1×1.0 = 4 bonds
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        info = ValenceValidator.get_valence_info(mol, 0)
        assert info.symbol == "C"
        assert info.is_valid

    def test_hydrogen_in_water_explicit_h(self):
        mol = Chem.AddHs(Chem.MolFromSmiles("O"))
        h_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
        assert len(h_atoms) > 0
        info = ValenceValidator.get_valence_info(mol, h_atoms[0].GetIdx())
        assert info.symbol == "H"
        assert info.current_bonds == 1
        assert info.is_valid

    def test_available_valence_is_max_minus_current(self):
        mol = Chem.MolFromSmiles("CC")
        # Ethane: each C has 4 bonds total — but we need explicit H
        mol = Chem.AddHs(mol)
        info = ValenceValidator.get_valence_info(mol, 0)
        assert info.available_valence == info.max_valence - info.current_bonds

    def test_saturated_flag_on_fully_bonded_carbon(self):
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        info = ValenceValidator.get_valence_info(mol, 0)
        assert info.is_saturated


class TestValenceValidatorMolecule:
    def test_benzene_is_valid(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        is_valid, errors = ValenceValidator.validate_molecule(mol)
        assert is_valid
        assert errors == []

    def test_ethanol_is_valid(self):
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        is_valid, errors = ValenceValidator.validate_molecule(mol)
        assert is_valid

    def test_all_valences_returns_correct_length(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        all_vals = ValenceValidator.get_all_valences(mol)
        assert len(all_vals) == mol.GetNumAtoms()

    def test_returns_list_of_valence_info_objects(self):
        from biochar.valence import ValenceInfo
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        all_vals = ValenceValidator.get_all_valences(mol)
        for vi in all_vals:
            assert isinstance(vi, ValenceInfo)


class TestValenceReport:
    def test_get_summary_returns_dict(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert isinstance(summary, dict)

    def test_summary_has_required_keys(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert "total_atoms" in summary
        assert "valid_atoms" in summary
        assert "invalid_atoms" in summary
        assert "all_valid" in summary
        assert "element_counts" in summary

    def test_summary_total_atoms(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert summary["total_atoms"] == mol.GetNumAtoms()

    def test_summary_all_valid_for_benzene(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert summary["all_valid"] is True
        assert summary["invalid_atoms"] == 0

    def test_summary_element_counts(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert summary["element_counts"]["C"] == 6
        assert summary["element_counts"]["H"] == 6


class TestSafeBondAdder:
    def _editable_mol(self, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        return mol

    def test_can_add_bond_between_unsaturated_atoms(self):
        # [CH2] radicals would have open valence, but safer: use editmol trick
        # Build CH4 + CH4 manually to test
        mol = Chem.RWMol()
        mol.AddAtom(Chem.Atom(6))  # C at 0
        mol.AddAtom(Chem.Atom(6))  # C at 1
        mol = mol.GetMol()
        can_add, reason = SafeBondAdder.can_add_bond(mol, 0, 1, bond_type=1)
        # Both carbons have 0 bonds so far → should be able to add
        assert can_add, reason

    def test_cannot_add_existing_bond(self):
        mol = Chem.MolFromSmiles("CC")
        # Atoms 0 and 1 are already bonded
        can_add, reason = SafeBondAdder.can_add_bond(mol, 0, 1, bond_type=1)
        assert not can_add
        assert "already exists" in reason.lower() or "bond" in reason.lower()

    def test_cannot_add_bond_to_saturated_hydrogen(self):
        mol = Chem.MolFromSmiles("[H][H]")  # H2 — H is saturated
        # Both H atoms are saturated (1 bond max, 1 bond current)
        h_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
        assert len(h_idxs) == 2
        can_add, reason = SafeBondAdder.can_add_bond(mol, h_idxs[0], h_idxs[1], bond_type=1)
        # Already bonded → cannot add
        assert not can_add
