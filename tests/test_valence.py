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

    # -- Characterization: charged-atom ranges relied on by the protonation
    # -- pipeline.  Pinned before get_valence_range was extended for anionic O
    # -- so that the extension is provably additive.

    def test_nitrogen_cation_is_exactly_four(self):
        """Pyridinium / anilinium N: 4 sigma bonds, no room for more."""
        assert get_valence_range(7, formal_charge=1) == (4, 4)

    def test_sulfur_cation_extends_beyond_two(self):
        """Hypervalent S (sulfone/sulfate) keeps its extended range."""
        lo, hi = get_valence_range(16, formal_charge=1)
        assert hi > 2

    def test_neutral_ranges_are_unaffected_by_charge_zero(self):
        """Every neutral range must be identical whether or not charge is passed."""
        for atomic_num in (1, 6, 7, 8, 9, 16, 17):
            assert get_valence_range(atomic_num) == get_valence_range(
                atomic_num, formal_charge=0
            ), f"neutral range for Z={atomic_num} shifted"

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


class TestAnionicValenceRanges:
    """
    Anionic heteroatoms carry one fewer bond than their neutral form.

    Without this, a deprotonated oxygen reports ``needed_valence == 1`` and
    :meth:`HydrogenAssigner._saturate_valences` silently re-protonates it,
    undoing the deprotonation with no error.
    """

    def test_oxygen_anion_is_exactly_one(self):
        """Phenolate / carboxylate O-: one bond, and no room for another."""
        assert get_valence_range(8, formal_charge=-1) == (1, 1)

    def test_sulfur_anion_is_exactly_one(self):
        """Thiophenolate S-: one bond, and no room for another."""
        assert get_valence_range(16, formal_charge=-1) == (1, 1)

    def test_nitrogen_anion_is_exactly_two(self):
        """Amide-like N-: two bonds."""
        assert get_valence_range(7, formal_charge=-1) == (2, 2)

    def test_oxygen_anion_needs_no_hydrogen(self):
        """The regression this whole unit exists to prevent."""
        min_val, _ = get_valence_range(8, formal_charge=-1)
        bonds_on_phenolate_o = 1
        assert max(0, min_val - bonds_on_phenolate_o) == 0

    def test_neutral_oxygen_still_needs_a_second_bond(self):
        """Guard: the neutral path must keep its existing behaviour."""
        min_val, _ = get_valence_range(8, formal_charge=0)
        assert max(0, min_val - 1) == 1


class TestAnionSurvivesHydrogenAssignment:
    """
    Integration coverage through the real HydrogenAssigner.

    The range tests above prove get_valence_range returns the right numbers.
    These prove the assigner actually honours them -- which is the behaviour the
    protonation feature depends on, and which no unit-level assertion can show.
    """

    @staticmethod
    def _phenolate() -> Chem.Mol:
        """Benzene bearing a single deprotonated O- (phenolate)."""
        mol = Chem.MolFromSmiles("[O-]c1ccccc1")
        assert mol is not None
        return Chem.AddHs(mol)

    def _anionic_oxygens(self, mol):
        return [
            a for a in mol.GetAtoms()
            if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1
        ]

    def test_phenolate_oxygen_keeps_its_charge_and_gains_no_hydrogen(self):
        from biochar.heteroatom_assignment import HydrogenAssigner

        mol = self._phenolate()
        before = self._anionic_oxygens(mol)
        assert len(before) == 1, "fixture should carry exactly one O-"
        assert before[0].GetTotalNumHs() == 0

        out, _ = HydrogenAssigner(seed=0).assign_hydrogens(mol, target_H_C_ratio=1.0)

        after = self._anionic_oxygens(out)
        assert len(after) == 1, "the O- lost its formal charge during H assignment"
        o = after[0]
        h_neighbours = [n for n in o.GetNeighbors() if n.GetAtomicNum() == 1]
        assert h_neighbours == [], (
            "HydrogenAssigner re-protonated the phenolate oxygen -- "
            "get_valence_range is not reporting needed_valence == 0 for O-"
        )

    def test_neutral_phenol_oxygen_still_keeps_its_hydrogen(self):
        """Regression guard: the neutral path must be unchanged."""
        from biochar.heteroatom_assignment import HydrogenAssigner

        mol = Chem.AddHs(Chem.MolFromSmiles("Oc1ccccc1"))
        out, _ = HydrogenAssigner(seed=0).assign_hydrogens(mol, target_H_C_ratio=1.0)

        oxygens = [a for a in out.GetAtoms() if a.GetAtomicNum() == 8]
        assert len(oxygens) == 1
        h_neighbours = [
            n for n in oxygens[0].GetNeighbors() if n.GetAtomicNum() == 1
        ]
        assert len(h_neighbours) == 1, "neutral phenol lost its hydroxyl H"


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
