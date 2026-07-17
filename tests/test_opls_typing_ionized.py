"""
AtomTyper coverage for ionized species and for the carboxyl/phenolic split.

Every fixture here is hand-built rather than generated, so the typing rules are
pinned independently of whatever the generator happens to produce.
"""

import pytest
from rdkit import Chem

from biochar.opls_typing import AtomTyper


def types_of(smiles: str) -> dict:
    """Map each heavy-atom symbol+index to its OPLS type for a SMILES string."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assigned = AtomTyper().assign_atom_types(mol)
    return mol, assigned


def types_for_element(smiles: str, atomic_num: int) -> list:
    mol, assigned = types_of(smiles)
    return [
        assigned[a.GetIdx()]
        for a in mol.GetAtoms()
        if a.GetAtomicNum() == atomic_num
    ]


def h_types_on(smiles: str, atomic_num: int) -> list:
    """OPLS types of the hydrogens attached to atoms of the given element."""
    mol, assigned = types_of(smiles)
    out = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != atomic_num:
            continue
        out += [
            assigned[n.GetIdx()] for n in a.GetNeighbors() if n.GetAtomicNum() == 1
        ]
    return out


class TestCarboxylVersusPhenolic:
    """
    A carboxyl -OH and a phenolic -OH both look like 'an O with one H'. Only the
    neighbouring carbon separates them, and before this split every carboxyl was
    typed as a phenol.
    """

    def test_carboxyl_hydroxyl_is_not_typed_as_phenol(self):
        assert "OH2" in types_for_element("OC(=O)c1ccccc1", 8)

    def test_carboxyl_carbonyl_is_not_typed_as_ketone(self):
        """A carboxylic C=O has its own type; it is not a ketone C=O (OC)."""
        oxygens = types_for_element("OC(=O)c1ccccc1", 8)
        assert set(oxygens) == {"OH2", "O"}, oxygens
        assert "OC" not in oxygens

    def test_carboxyl_hydroxyl_hydrogen_follows_its_oxygen(self):
        assert h_types_on("OC(=O)c1ccccc1", 8) == ["HO2"]

    def test_carboxyl_carbon_types_as_carboxyl_carbon(self):
        assert "C" in types_for_element("OC(=O)c1ccccc1", 6)

    def test_phenol_still_types_as_phenol(self):
        assert types_for_element("Oc1ccccc1", 8) == ["OH"]

    def test_phenol_hydrogen_still_types_as_ho(self):
        assert h_types_on("Oc1ccccc1", 8) == ["HO"]

    def test_ether_oxygen_is_untouched(self):
        assert types_for_element("c1ccc2c(c1)Oc1ccccc1-2", 8) == ["OS"]

    def test_quinone_carbonyl_is_still_a_ketone_oxygen(self):
        """A real ketone C=O must not be caught by the carboxyl branch."""
        assert types_for_element("O=C1C=CC(=O)C=C1", 8) == ["OC", "OC"]


class TestAnionicOxygen:
    def test_phenolate_types_as_phenolate(self):
        assert types_for_element("[O-]c1ccccc1", 8) == ["OM"]

    def test_phenolate_is_not_typed_as_neutral_phenol(self):
        assert "OH" not in types_for_element("[O-]c1ccccc1", 8)

    def test_both_carboxylate_oxygens_type_as_carboxylate(self):
        """
        The two oxygens of a carboxylate are equivalent by resonance. RDKit's
        Kekule form puts the charge on only one of them, so typing the other
        from its double bond alone would strand it on the neutral acid's type.
        """
        oxygens = types_for_element("[O-]C(=O)c1ccccc1", 8)
        assert oxygens == ["O2M", "O2M"], (
            f"carboxylate oxygens are not equivalent: {oxygens}"
        )

    def test_carboxylate_does_not_retain_the_neutral_acid_carbonyl_type(self):
        assert "O" not in types_for_element("[O-]C(=O)c1ccccc1", 8)

    def test_neutral_acid_carbonyl_is_unaffected_by_the_carboxylate_branch(self):
        """Guard: a protonated COOH must keep its distinct C=O type."""
        assert types_for_element("OC(=O)c1ccccc1", 8) == ["OH2", "O"]

    def test_carboxylate_is_not_typed_as_phenolate(self):
        """The neighbour is what separates the two anionic oxygens."""
        oxygens = types_for_element("[O-]C(=O)c1ccccc1", 8)
        assert "OM" not in oxygens

    def test_carboxylate_carbon_still_types_as_carboxyl_carbon(self):
        assert "C" in types_for_element("[O-]C(=O)c1ccccc1", 6)

    def test_phenolate_has_no_hydrogen(self):
        assert h_types_on("[O-]c1ccccc1", 8) == []


class TestAnionicSulfur:
    def test_thiophenolate_types_as_thiolate(self):
        assert types_for_element("[S-]c1ccccc1", 16) == ["SM"]

    def test_neutral_thiol_still_types_as_thiol(self):
        assert types_for_element("Sc1ccccc1", 16) == ["SH_"]

    def test_thioether_still_types_as_thioether(self):
        assert types_for_element("c1ccc2c(c1)Sc1ccccc1-2", 16) == ["SS"]


class TestCationicNitrogen:
    def test_pyridinium_types_as_pyridinium(self):
        assert types_for_element("c1cc[nH+]cc1", 7) == ["NPYP"]

    def test_pyridinium_is_not_typed_as_neutral_pyridinic(self):
        assert "NPY" not in types_for_element("c1cc[nH+]cc1", 7)

    def test_pyridinium_hydrogen_types_as_pyridinium_h(self):
        assert h_types_on("c1cc[nH+]cc1", 7) == ["HPYP"]

    def test_neutral_pyridine_still_types_as_pyridinic(self):
        assert types_for_element("c1ccncc1", 7) == ["NPY"]

    def test_anilinium_types_as_anilinium(self):
        assert types_for_element("[NH3+]c1ccccc1", 7) == ["NAP"]

    def test_anilinium_is_not_typed_as_neutral_aniline(self):
        assert "NA" not in types_for_element("[NH3+]c1ccccc1", 7)

    def test_anilinium_hydrogens_type_as_anilinium_h(self):
        assert h_types_on("[NH3+]c1ccccc1", 7) == ["HNAP"] * 3

    def test_neutral_aniline_still_types_as_aniline(self):
        assert types_for_element("Nc1ccccc1", 7) == ["NA"]

    def test_pyrrolic_nitrogen_is_untouched(self):
        assert types_for_element("c1cc[nH]c1", 7) == ["NPR"]


class TestGraphiticNitrogenNotCaughtByCationicBranch:
    """
    Graphitic N already carries a formal +1, so the new cationic branch could
    easily swallow it. It is quaternary and hydrogen-free, which is what keeps
    it distinct from pyridinium.
    """

    @staticmethod
    def _graphitic():
        from biochar.heteroatom_assignment import NitrogenSubstitutor

        # Coronene has genuine interior junction carbons to substitute.
        mol = Chem.MolFromSmiles("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61")
        sub = NitrogenSubstitutor(seed=1)
        out = sub.substitute(mol, n_graphitic=1)
        assert sub.placed_graphitic == 1, "fixture failed to place graphitic N"
        return out

    def test_graphitic_nitrogen_still_types_as_ngr(self):
        mol = self._graphitic()
        assigned = AtomTyper().assign_atom_types(mol)
        n_types = [
            assigned[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() == 7
        ]
        assert n_types == ["NGR"], f"graphitic N typed as {n_types}"

    def test_graphitic_nitrogen_is_cationic_but_not_pyridinium(self):
        mol = self._graphitic()
        n = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
        assert n.GetFormalCharge() == 1, "fixture assumption: graphitic N is +1"
        assigned = AtomTyper().assign_atom_types(mol)
        assert assigned[n.GetIdx()] != "NPYP"


class TestNoUnrecognisedTypes:
    @pytest.mark.parametrize(
        "smiles",
        [
            "[O-]C(=O)c1ccccc1",   # carboxylate
            "[O-]c1ccccc1",        # phenolate
            "[S-]c1ccccc1",        # thiophenolate
            "c1cc[nH+]cc1",        # pyridinium
            "[NH3+]c1ccccc1",      # anilinium
            "OC(=O)c1ccccc1",      # carboxylic acid
        ],
    )
    def test_every_atom_gets_a_real_type(self, smiles):
        mol, assigned = types_of(smiles)
        unknown = {
            a.GetSymbol(): assigned[a.GetIdx()]
            for a in mol.GetAtoms()
            if assigned[a.GetIdx()].startswith("X")
        }
        assert not unknown, f"{smiles} produced unrecognised types: {unknown}"

    @pytest.mark.parametrize(
        "smiles",
        [
            "[O-]C(=O)c1ccccc1",
            "[O-]c1ccccc1",
            "[S-]c1ccccc1",
            "c1cc[nH+]cc1",
            "[NH3+]c1ccccc1",
            "OC(=O)c1ccccc1",
        ],
    )
    def test_every_type_is_exportable_to_gromacs(self, smiles):
        """A type with no GROMACS mapping would be written raw into the .top."""
        from biochar.constants import GROMACS_OPLS_TYPE_MAP

        mol, assigned = types_of(smiles)
        missing = {t for t in assigned.values() if t not in GROMACS_OPLS_TYPE_MAP}
        assert not missing, f"{smiles} produced unmappable types: {missing}"
