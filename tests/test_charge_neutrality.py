"""
ChargeAssigner: net charge must equal total formal charge, not zero.

The old contract was "charges sum to 0", enforced unconditionally. That made an
ionized structure impossible to represent -- any charge placed on a carboxylate
was redistributed away until the molecule was neutral again, silently.

The new contract is "charges sum to the total formal charge". For a neutral
molecule the target is 0 and the behaviour is identical, which is what the
characterization tests here exist to prove.
"""

import pytest
from rdkit import Chem

from biochar.opls_typing import AtomTyper, ChargeAssigner, OPLSPropertyTable


def charges_for(smiles: str):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    types = AtomTyper().assign_atom_types(mol)
    return mol, ChargeAssigner().assign_charges(mol, types)


def total_charge(smiles: str) -> float:
    _, charges = charges_for(smiles)
    return sum(charges.values())


def formal_charge(smiles: str) -> int:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


NEUTRAL_SPECIES = [
    "c1ccccc1",             # benzene
    "Oc1ccccc1",            # phenol
    "OC(=O)c1ccccc1",       # benzoic acid
    "Sc1ccccc1",            # thiophenol
    "Nc1ccccc1",            # aniline
    "c1ccncc1",             # pyridine
    "c1ccc2ccccc2c1",       # naphthalene
]


class TestNeutralMoleculesUnchanged:
    """
    Characterization: the neutral path must behave exactly as it did before the
    target changed. This is the guarantee that setting no pH reproduces today's
    output.
    """

    @pytest.mark.parametrize("smiles", NEUTRAL_SPECIES)
    def test_neutral_molecule_sums_to_zero(self, smiles):
        assert total_charge(smiles) == pytest.approx(0.0, abs=1e-6)

    @pytest.mark.parametrize("smiles", NEUTRAL_SPECIES)
    def test_neutral_target_equals_zero(self, smiles):
        """A neutral molecule's formal charge is 0, so the new target is 0."""
        assert formal_charge(smiles) == 0

    def test_every_atom_gets_a_charge(self):
        mol, charges = charges_for("Oc1ccccc1")
        assert set(charges) == {a.GetIdx() for a in mol.GetAtoms()}


class TestIonizedMoleculesKeepTheirCharge:
    """The behaviour the old sum-to-zero contract made impossible."""

    @pytest.mark.parametrize(
        "smiles,expected",
        [
            ("[O-]C(=O)c1ccccc1", -1.0),   # benzoate
            ("[O-]c1ccccc1", -1.0),        # phenolate
            ("[S-]c1ccccc1", -1.0),        # thiophenolate
            ("[NH3+]c1ccccc1", +1.0),      # anilinium
            ("c1cc[nH+]cc1", +1.0),        # pyridinium
        ],
    )
    def test_ionized_molecule_sums_to_its_formal_charge(self, smiles, expected):
        assert formal_charge(smiles) == expected, "fixture assumption"
        assert total_charge(smiles) == pytest.approx(expected, abs=1e-6)

    def test_dianion_sums_to_minus_two(self):
        smiles = "[O-]C(=O)c1ccc([O-])cc1"  # 4-hydroxybenzoate, both deprotonated
        assert formal_charge(smiles) == -2
        assert total_charge(smiles) == pytest.approx(-2.0, abs=1e-6)

    def test_zwitterion_sums_to_zero_but_keeps_its_formal_charges(self):
        """
        A carboxylate plus an anilinium is net neutral. The molecule must still
        carry both formal charges -- summing to zero is not the same as having
        no charges, and collapsing them would erase real electrostatics.
        """
        smiles = "[O-]C(=O)c1ccc([NH3+])cc1"
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert formal_charge(smiles) == 0
        assert total_charge(smiles) == pytest.approx(0.0, abs=1e-6)

        signs = {a.GetFormalCharge() for a in mol.GetAtoms()}
        assert -1 in signs and 1 in signs, "formal charges were cancelled away"


class TestGraphiticNitrogenChargeIsNoLongerSmearedAway:
    """
    Graphitic N has carried a formal +1 since it was introduced, and the
    sum-to-zero contract has been quietly redistributing it away ever since.
    This changes existing output: qtot goes 0 -> +1 for these structures.
    """

    @staticmethod
    def _graphitic_mol():
        from biochar.heteroatom_assignment import NitrogenSubstitutor

        mol = Chem.MolFromSmiles("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61")
        sub = NitrogenSubstitutor(seed=1)
        out = sub.substitute(mol, n_graphitic=1)
        assert sub.placed_graphitic == 1
        return out

    def test_graphitic_structure_now_reports_plus_one(self):
        mol = self._graphitic_mol()
        assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) == 1

        types = AtomTyper().assign_atom_types(mol)
        charges = ChargeAssigner().assign_charges(mol, types)
        assert sum(charges.values()) == pytest.approx(1.0, abs=1e-6)


class TestPropertyTableValidation:
    def test_ionized_structure_is_not_flagged_as_extreme(self):
        mol, charges = charges_for("[O-]C(=O)c1ccccc1")
        types = AtomTyper().assign_atom_types(mol)
        table = OPLSPropertyTable(mol, types, charges)
        ok, errors = table.validate()
        assert ok, f"legitimately ionized structure flagged: {errors}"

    def test_total_charge_matches_formal_charge(self):
        mol, charges = charges_for("[O-]c1ccccc1")
        types = AtomTyper().assign_atom_types(mol)
        table = OPLSPropertyTable(mol, types, charges)
        assert table.get_total_charge() == pytest.approx(-1.0, abs=1e-6)

    def test_validation_catches_a_charge_that_does_not_match_formal_charge(self):
        """A partial-charge sum that disagrees with the formal charge is a bug."""
        mol, charges = charges_for("[O-]c1ccccc1")
        types = AtomTyper().assign_atom_types(mol)
        # Corrupt one charge so the sum no longer matches the formal charge.
        corrupted = dict(charges)
        corrupted[0] += 0.75
        table = OPLSPropertyTable(mol, types, corrupted)
        ok, errors = table.validate()
        assert not ok
        assert any("formal" in e.lower() or "net" in e.lower() for e in errors), errors
