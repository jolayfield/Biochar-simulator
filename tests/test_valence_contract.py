"""
Tests for the rqm/valence-validation.md scenarios -- whether the bond count the
rest of the pipeline is built on is the right one.

Distinct from tests/test_valence.py, which checks the tables and the happy path
on saturated aliphatic fragments. These use the ring chemistries this package
actually builds: a furan-type bridge, a pyrrolic nitrogen, a phenolate. Those
are where summing aromatic bond orders and where an anion's traded lone pair
both stop agreeing with the chemistry.

Six of the seven scenarios that carried xfail(strict=True) were retired by
XPASS when their fixes landed. The seventh -- pyrrolic doping -- had three
causes; two are fixed and the marker stays, naming the one that is not.
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
    def test_a_furan_type_ring_oxygen_is_valid(self):
        # Dibenzofuran: the oxygen is an aromatic ring member with two sigma bonds.
        assert not _violations(_mol("c1ccc2c(c1)Oc1ccccc1-2"), "O")

    # rq-302b177e
    def test_a_pyrrolic_ring_nitrogen_is_valid(self):
        assert not _violations(_mol("c1cc[nH]c1"), "N")

    # rq-e9202195
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
    # Several seeds, because the outcome depends on which carbons the shuffle
    # offers. A single seed here passes or fails by luck and says nothing.
    SEEDS = [0, 1, 2, 3, 4, 6, 8]

    # rq-ee235774
    @pytest.mark.xfail(
        strict=True,
        reason="two of the three causes are fixed -- the site is no longer a "
               "decorated carbon, and no longer a second nitrogen in the same "
               "ring, and the RingInfo crash that used to abort some seeds "
               "outright is fixed. The last cause remains: on some skeletons "
               "the pentagon loses "
               "aromaticity downstream and kekulises, leaving the N-H on a C=N "
               "and the nitrogen at four bonds. Chasing that means working "
               "through the H/C-shaping and sanitisation interaction",
    )
    def test_a_doped_structure_has_no_over_bonded_nitrogen(self):
        offenders = {}
        for seed in self.SEEDS:
            _, mol = _generate(seed=seed, num_pyrrolic=2, O_C_ratio=0.15)
            violations = _violations(mol, "N")
            if violations:
                offenders[seed] = violations[0]

        assert not offenders, (
            f"{len(offenders)} of {len(self.SEEDS)} seeds produced a nitrogen "
            f"carrying more bonds than nitrogen can: {offenders}"
        )

    # rq-070a2c21
    def test_pyrrolic_nitrogen_is_still_placed(self):
        """The site check must not answer the above by placing nothing.

        A stricter rule that finds no qualifying site would make every seed
        pass the scenario above while making the feature useless, so the
        placement rate is pinned across the same seeds.
        """
        placed = 0
        for seed in self.SEEDS:
            _, mol = _generate(seed=seed, num_pyrrolic=2, O_C_ratio=0.15)
            placed += sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")

        assert placed >= len(self.SEEDS), (
            f"only {placed} nitrogens placed across {len(self.SEEDS)} seeds "
            f"asking for two each; the site rule has become too strict to use"
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
    def test_a_newly_added_atom_can_be_bonded(self):
        mol = _carbon_with_two_free_valences()
        emol = Chem.EditableMol(mol)
        new_idx = SafeBondAdder.add_atom_safe(emol, 1)

        added, message = SafeBondAdder.add_bond_safe(emol, mol, 0, new_idx)
        assert added, message

    # rq-998f84fb
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
    def test_an_aromatic_bond_is_weighed_as_one_bond(self):
        """Two ring carbons stripped to a genuine free valence.

        Removing the hydrogen is not enough on its own -- RDKit puts an
        implicit one back -- so each is pinned with SetNoImplicit.
        """
        rw = Chem.RWMol(_mol("c1ccccc1"))
        # Two carbons on opposite sides of the ring, so they are not already
        # bonded -- an existing bond is refused for a different reason.
        hydrogens = [
            h.GetIdx()
            for h in rw.GetAtoms()
            if h.GetAtomicNum() == 1
            and h.GetNeighbors()[0].GetIdx() in (0, 3)
        ]
        parents = [rw.GetAtomWithIdx(h).GetNeighbors()[0].GetIdx() for h in hydrogens]
        hydrogens = sorted(hydrogens, reverse=True)
        for h in hydrogens:
            rw.RemoveAtom(h)
        for c in parents:
            atom = rw.GetAtomWithIdx(c)
            atom.SetNoImplicit(True)
            atom.SetNumExplicitHs(0)
        stripped = rw.GetMol()

        free = [
            c
            for c in parents
            if ValenceValidator.get_valence_info(stripped, c).available_valence >= 1
        ]
        assert len(free) == 2, (
            f"fixture is void: {len(free)} of 2 stripped carbons have a free "
            f"valence"
        )

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


# --------------------------------------------------------------------------- #
# Typing a molecule that has not been sanitised
# --------------------------------------------------------------------------- #
class TestTypingAnUnperceivedMolecule:
    # rq-c6ab7cbe
    def test_a_molecule_with_no_ring_information_is_typed(self):
        """RDKit raises rather than answering "is this in a ring?" when nothing
        has perceived the rings. Molecules arrive here that way from the
        nitrogen-doping paths, and the whole generation used to abort on it.
        """
        from biochar.pipeline.opls_typing import AtomTyper

        rw = Chem.RWMol()
        idx = [rw.AddAtom(Chem.Atom(6)) for _ in range(6)]
        for a, b in zip(idx, idx[1:] + idx[:1]):
            rw.AddBond(a, b, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        mol.UpdatePropertyCache(strict=False)
        with pytest.raises(RuntimeError):
            mol.GetRingInfo().NumRings()

        types = AtomTyper().assign_atom_types(mol)
        assert set(types) == set(idx)
