"""
Tests for the rqm/heteroatom-assignment.md scenarios not already defended
elsewhere: the carbonyl/quinone/lactone-to-phenolic fallback, the oxygen count
surviving that substitution, the ether bond-type repair surviving every
re-sanitising stage, full valence saturation, and the small-molecule H/C
perimeter ceiling.

The fallback group is the sharpest gap this file closes. CLAUDE.md documents
it as a known surprise -- carbonyl, quinone, and lactone are accepted by the
API and silently placed as phenolic instead -- but until now nothing asserted
it. Deleting the fallback entirely would have kept every other test green.
"""

import pytest
from rdkit import Chem

from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig
from biochar.pipeline.heteroatom_assignment import _fix_heteroatom_bond_types
from biochar.pipeline.valence import ValenceValidator


def _generate(**kw):
    cfg = GeneratorConfig(**kw)
    mol, coords, comp = BiocharGenerator(cfg).generate()
    return mol, comp


def _num_co_double_bonds(mol):
    return sum(
        1 for b in mol.GetBonds()
        if b.GetBondType() == Chem.BondType.DOUBLE
        and {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()} == {6, 8}
    )


def _num_hydroxyls(mol):
    return len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))


def _aromatic_fraction(mol):
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    return sum(1 for a in carbons if a.GetIsAromatic()) / len(carbons)


class TestUnsupportedGroupsFallBackToPhenolic:
    """A pure aromatic PAH edge carbon has one free valence: not enough for a
    carbonyl's second C=O bond, a quinone's paired carbonyls, or a lactone's
    ring closure. All three are silently placed as phenolic instead."""

    @pytest.mark.parametrize(
        "group,rq",
        [
            ("carbonyl", "rq-d3ef932d"),
            ("quinone", "rq-7d3e9432"),
            ("lactone", "rq-2f5ab3f8"),
        ],
    )
    def test_request_places_a_phenolic_group(self, group, rq):
        # rq id above documents which scenario this parametrised case defends.
        phenolic_mol, _ = _generate(
            target_num_carbons=24, functional_groups={"phenolic": 2},
            molecule_name="fb", seed=3, strict=False,
        )
        fallback_mol, _ = _generate(
            target_num_carbons=24, functional_groups={group: 2},
            molecule_name="fb", seed=3, strict=False,
        )

        # Same seed and skeleton target, so a real substitution (not a
        # different structure entirely) produces an identical result: no C=O
        # double bond anywhere (a phenolic ring carbon has none to give), and
        # the same hydroxyl count as an explicit phenolic request.
        assert _num_co_double_bonds(fallback_mol) == 0
        assert _num_hydroxyls(fallback_mol) == _num_hydroxyls(phenolic_mol)

    # rq-fbb6a9c0
    def test_oxygen_count_still_matches_the_request(self):
        for group in ("carbonyl", "quinone", "lactone"):
            _, comp = _generate(
                target_num_carbons=24, functional_groups={group: 2},
                molecule_name="fb", seed=3, strict=False,
            )
            assert comp.num_oxygens == 2, (
                f"'{group}' requested 2 O via 2 groups; substitution must not "
                f"change the requested count, only which group supplies it"
            )


class TestEtherRepairSurvivesEveryStage:
    # rq-b8fce9a6
    def test_repair_runs_more_than_once_and_bonds_stay_single(self, monkeypatch):
        calls = {"n": 0}
        orig = _fix_heteroatom_bond_types

        def spy(mol):
            calls["n"] += 1
            return orig(mol)

        # biochar_generator imported the name directly (`from
        # .heteroatom_assignment import _fix_heteroatom_bond_types`), so it
        # holds its own binding distinct from the defining module's -- both
        # must be patched to see every call site.
        monkeypatch.setattr("biochar.pipeline.biochar_generator._fix_heteroatom_bond_types", spy)
        monkeypatch.setattr("biochar.pipeline.heteroatom_assignment._fix_heteroatom_bond_types", spy)

        mol, _ = _generate(
            target_num_carbons=24, functional_groups={"ether": 2},
            molecule_name="et", seed=1, strict=False,
        )

        # "after assignment, after geometry, and after validation" -- at
        # least the two call sites in biochar_generator.py fire in addition to
        # whatever heteroatom_assignment.py does internally during placement.
        assert calls["n"] >= 2

        ether_bonds = [
            b for b in mol.GetBonds()
            if {b.GetBeginAtom().GetAtomicNum(), b.GetEndAtom().GetAtomicNum()} == {6, 8}
        ]
        assert ether_bonds, "expected at least one ether C-O bond"
        for b in ether_bonds:
            assert b.GetBondType() == Chem.BondType.SINGLE, (
                "an ether C-O bond survived re-sanitisation as AROMATIC"
            )


class TestNoFreeValenceAfterAssignment:
    # rq-85c279ce
    @pytest.mark.parametrize(
        "kw",
        [
            dict(target_num_carbons=24, H_C_ratio=0.5, O_C_ratio=0.1),
            dict(
                target_num_carbons=60, H_C_ratio=0.6, O_C_ratio=0.15,
                functional_groups={"carboxyl": 1, "ether": 1, "amino": 1, "thiol": 1},
            ),
        ],
    )
    def test_every_atom_is_fully_saturated(self, kw):
        mol, _ = _generate(molecule_name="sv", seed=2, strict=False, **kw)
        bad = [
            info for info in ValenceValidator.get_all_valences(mol)
            if info.available_valence != 0
        ]
        assert not bad, f"atoms with an unfilled valence: {bad}"


class TestSmallMoleculesAreBoundedByPerimeter:
    """A fully aromatic flake can only carry hydrogen on its perimeter: one H
    per aromatic carbon is the absolute ceiling, since a second H would break
    aromaticity.

    allow_aliphatic=False forces the genuinely fully-aromatic case the
    requirement describes. Without it, a high enough H/C request makes the
    generator reach for aliphatic decoration even at 10-14 carbons (verified:
    at H_C_ratio=0.9 with the default allow_aliphatic=True, a 10- or 14-carbon
    structure comes back partially sp3, aromatic_fraction 0.875 not 1.0) --
    a different, already-covered mechanism (see test_hc_control.py's Step 3),
    not the perimeter ceiling this scenario is about.
    """

    # rq-8176d084
    @pytest.mark.parametrize("num_carbons", [6, 10, 14])
    def test_achieved_h_c_never_exceeds_the_aromatic_ceiling(self, num_carbons):
        mol, comp = _generate(
            target_num_carbons=num_carbons, H_C_ratio=0.9, O_C_ratio=0.0,
            allow_aliphatic=False, molecule_name="sm", seed=1, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)
        assert comp.H_C_ratio <= 1.0 + 1e-9, (
            "a fully aromatic flake cannot carry more than one H per carbon"
        )
