"""
Tests for the rqm/valence-validation.md scenarios -- whether the bond count the
rest of the pipeline is built on is the right one.

Distinct from tests/test_valence.py, which checks the tables and the happy path
on saturated aliphatic fragments. These use the ring chemistries this package
actually builds: a furan-type bridge, a pyrrolic nitrogen, a phenolate. Those
are where summing aromatic bond orders and where an anion's traded lone pair
both stop agreeing with the chemistry.

Seven scenarios carry xfail(strict=True); each names the defect it defers.
"""

import io
import contextlib
import warnings

import pytest
from rdkit import Chem

from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig
from biochar.pipeline.valence import (
    SafeBondAdder,
    ValenceValidator,
)


def _mol(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.SanitizeMol(mol)
    return mol


def _violations(mol, symbol=None):
    _, errors = ValenceValidator.validate_molecule(mol)
    if symbol is None:
        return errors
    return [e for e in errors if f"({symbol})" in e]


# --------------------------------------------------------------------------- #
# Aromatic ring members
# --------------------------------------------------------------------------- #
class TestAnAromaticRingMemberIsCountedByItsBonds:
    # rq-c70dd417
    @pytest.mark.xfail(
        strict=True,
        reason="aromatic bond orders are summed at 1.5 each, so a ring oxygen "
               "donating its lone pair reports three bonds and is refused for "
               "exceeding two -- the furan bridge is a core structure here",
    )
    def test_a_furan_type_ring_oxygen_is_valid(self):
        # Dibenzofuran: the oxygen is an aromatic ring member with two sigma bonds.
        assert not _violations(_mol("c1ccc2c(c1)Oc1ccccc1-2"), "O")

    # rq-302b177e
    @pytest.mark.xfail(
        strict=True,
        reason="same cause: a pyrrolic nitrogen holds two aromatic bonds and an "
               "N-H, which sums to four against a maximum of three",
    )
    def test_a_pyrrolic_ring_nitrogen_is_valid(self):
        assert not _violations(_mol("c1cc[nH]c1"), "N")

    # rq-e9202195
    @pytest.mark.xfail(
        strict=True,
        reason="same cause again: sulfur's extended range is consulted only "
               "for a positive formal charge, so a neutral aromatic ring "
               "sulfur is judged against a maximum of two",
    )
    def test_a_thiophene_ring_sulfur_is_valid(self):
        assert not _violations(_mol("c1ccsc1"), "S")

    # rq-c179a6e7
    def test_a_genuinely_over_bonded_nitrogen_is_still_reported(self):
        """The counting fix must not answer the two above by not looking."""
        # Neutral N with four single bonds: not a molecule.
        rw = Chem.RWMol(Chem.MolFromSmiles("N"))
        n = rw.GetAtomWithIdx(0)
        n.SetNoImplicit(True)
        n.SetNumExplicitHs(0)
        for _ in range(4):
            h = rw.AddAtom(Chem.Atom(1))
            rw.AddBond(0, h, Chem.BondType.SINGLE)
        mol = rw.GetMol()

        errors = _violations(mol, "N")
        assert errors, "a neutral four-bonded nitrogen was reported as valid"
        assert "Atom 0" in errors[0] and "3" in errors[0], errors


# --------------------------------------------------------------------------- #
# Anions
# --------------------------------------------------------------------------- #
class TestAnAnionHasOneFewerBondAndNoRoom:
    @staticmethod
    def _phenolate_oxygen():
        mol = _mol("[O-]c1ccccc1")
        idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() == -1)
        return ValenceValidator.get_valence_info(mol, idx)

    # rq-eb28d0a1
    def test_a_deprotonated_oxygen_needs_no_further_bond(self):
        info = self._phenolate_oxygen()
        assert info.needed_valence == 0, (
            "a phenolate looks like it still needs a hydrogen; saturation would "
            "put one back and undo the deprotonation"
        )

    # rq-a467c647
    def test_a_deprotonated_oxygen_has_no_free_valence(self):
        assert self._phenolate_oxygen().available_valence == 0


# --------------------------------------------------------------------------- #
# Nitrogen doping
# --------------------------------------------------------------------------- #
def _generate(**kw):
    cfg = dict(target_num_carbons=60, seed=2, strict=False, defect_fraction=0.3,
               H_C_tolerance=1.0, O_C_tolerance=1.0)
    cfg.update(kw)
    gen = BiocharGenerator(GeneratorConfig(**cfg))
    with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("ignore")
        mol, _, _ = gen.generate()
    return gen, mol


class TestPyrrolicDopingPicksASiteItCanUse:
    # rq-ee235774
    @pytest.mark.xfail(
        strict=True,
        reason="pyrrolic substitution takes any five-ring carbon, including one "
               "already carrying a functional-group oxygen, and then adds the "
               "N-H on top -- four bonds on a neutral nitrogen",
    )
    def test_a_doped_structure_has_no_over_bonded_nitrogen(self):
        _, mol = _generate(num_pyrrolic=2, O_C_ratio=0.15)

        nitrogens = [a for a in mol.GetAtoms() if a.GetSymbol() == "N"]
        assert nitrogens, "fixture is void: no nitrogen was placed"
        assert not _violations(mol, "N"), (
            "a doped nitrogen carries more bonds than nitrogen can: "
            + "; ".join(_violations(mol, "N"))
        )

    # rq-070a2c21
    def test_pyrrolic_nitrogen_is_still_placed(self):
        """The site check must not answer the above by placing nothing."""
        gen, mol = _generate(num_pyrrolic=2, O_C_ratio=0.15)
        assert sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N") >= 1, (
            "no nitrogen was placed at all; a stricter site rule that finds no "
            "site is not a fix"
        )


# --------------------------------------------------------------------------- #
# SafeBondAdder
# --------------------------------------------------------------------------- #
def _carbon_with_two_free_valences():
    """CH2 as an explicit-H fragment: one carbon, two hydrogens, room for two."""
    rw = Chem.RWMol()
    c = rw.AddAtom(Chem.Atom(6))
    rw.GetAtomWithIdx(c).SetNoImplicit(True)
    for _ in range(2):
        h = rw.AddAtom(Chem.Atom(1))
        rw.AddBond(c, h, Chem.BondType.SINGLE)
    return rw.GetMol()


class TestABondIsJudgedAgainstTheMoleculeBeingEdited:
    # rq-e0d6b415
    @pytest.mark.xfail(
        strict=True,
        reason="add_bond_safe validates against the snapshot it was handed, so "
               "an atom appended to the editor does not exist there and the "
               "lookup raises -- which is the workflow VALENCE_SYSTEM.md shows",
    )
    def test_a_newly_added_atom_can_be_bonded(self):
        mol = _carbon_with_two_free_valences()
        emol = Chem.EditableMol(mol)
        new_idx = SafeBondAdder.add_atom_safe(emol, 1)

        added, message = SafeBondAdder.add_bond_safe(emol, mol, 0, new_idx)
        assert added, message

    # rq-998f84fb
    @pytest.mark.xfail(
        strict=True,
        reason="each addition is weighed against the state before any of them, "
               "so a carbon with room for two accepts four",
    )
    def test_a_sequence_of_bonds_is_bounded_by_the_atom_s_valence(self):
        mol = _carbon_with_two_free_valences()
        emol = Chem.EditableMol(mol)
        partners = [SafeBondAdder.add_atom_safe(emol, 1) for _ in range(4)]

        accepted = sum(
            SafeBondAdder.add_bond_safe(emol, mol, 0, idx)[0] for idx in partners
        )
        assert accepted == 2, (
            f"a carbon with room for two more bonds accepted {accepted}"
        )
        assert emol.GetMol().GetAtomWithIdx(0).GetDegree() == 4

    # rq-95aaba31
    @pytest.mark.xfail(
        strict=True,
        reason="Chem.BondType.AROMATIC is the enumeration value 12, so an "
               "aromatic bond asked for by name demands twelve free valences "
               "and is refused on every atom of every molecule",
    )
    def test_an_aromatic_bond_is_weighed_as_one_bond(self):
        # Two biphenyl-like carbons, each with a free valence after removing an H.
        mol = _mol("c1ccccc1")
        rw = Chem.RWMol(mol)
        hydrogens = sorted(
            (a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1),
            reverse=True,
        )[:2]
        for h in hydrogens:
            rw.RemoveAtom(h)
        stripped = rw.GetMol()
        free = [
            a.GetIdx()
            for a in stripped.GetAtoms()
            if ValenceValidator.get_valence_info(stripped, a.GetIdx()).available_valence
        ]
        assert len(free) >= 2, "fixture is void: no carbon has a free valence"

        can_add, reason = SafeBondAdder.can_add_bond(
            stripped, free[0], free[1], bond_type=Chem.BondType.AROMATIC
        )
        assert can_add, f"an aromatic bond was refused: {reason}"


# --------------------------------------------------------------------------- #
# Reporting
# --------------------------------------------------------------------------- #
class TestEveryViolationIsNamed:
    # rq-08130b9d
    def test_a_correct_molecule_validates_with_no_messages(self):
        is_valid, errors = ValenceValidator.validate_molecule(_mol("c1ccccc1"))
        assert is_valid and errors == []

    # rq-b8611a4d
    def test_an_atom_below_its_minimum_is_named(self):
        rw = Chem.RWMol(Chem.MolFromSmiles("C"))
        rw.GetAtomWithIdx(0).SetNoImplicit(True)
        rw.GetAtomWithIdx(0).SetNumExplicitHs(0)
        mol = rw.GetMol()

        is_valid, errors = ValenceValidator.validate_molecule(mol)
        assert not is_valid
        assert "Atom 0" in errors[0] and "below minimum 4" in errors[0], errors
