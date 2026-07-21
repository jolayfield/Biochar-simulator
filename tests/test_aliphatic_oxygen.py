"""
Tests for oxygen placement on aliphatic (sp3) carbons.

Oxygen groups historically attached only to aromatic edge carbons
(`_get_edge_sites`).  On low-aromaticity chars the H/C-shaping stage fills those
edges with pendant methyls and H-saturation, leaving almost none for oxygen, so
high-O/C low-temperature points could not reach their O/C target and fell to the
sweep's fallback path (worst case: hardwood 300 °C placed 1 O against a target of
13).

Aliphatic hydroxyls (-CH2-OH on the sp3 carbons) lift that ceiling.  These tests
cover the new site pool, the placement, the O_C_ratio-mode spill that reaches the
target, the opt-out flag, and that the group types cleanly with no new OPLS type.
"""

import pytest

from rdkit import Chem

from biochar.constants import FUNCTIONAL_GROUPS, GROMACS_OPLS_TYPE_MAP
from biochar.heteroatom_assignment import (
    OxygenAssigner,
    attach_aliphatic_carbons,
)
from biochar.opls_typing import AtomTyper, ChargeAssigner
import random


def _aromatic_flake_with_methyls(n_methyl: int, seed: int = 0):
    """A small PAH (pyrene) decorated with *n_methyl* pendant sp3 carbons."""
    mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")  # pyrene, 16 C
    Chem.SanitizeMol(mol)
    return attach_aliphatic_carbons(mol, n_methyl, random.Random(seed))


def _count_oxygens(mol):
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)


def _aliphatic_hydroxyl_carbons(mol):
    """sp3 carbons that carry an -OH (i.e. a hydroxymethyl carbon)."""
    out = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 6 or a.GetIsAromatic():
            continue
        for nbr in a.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and any(
                x.GetAtomicNum() == 1 for x in nbr.GetNeighbors()
            ):
                out.append(a.GetIdx())
    return out


# ---------------------------------------------------------------------------
# Registry + site pool
# ---------------------------------------------------------------------------

class TestRegistryAndSites:
    def test_aliphatic_hydroxyl_is_registered(self):
        assert "aliphatic_hydroxyl" in FUNCTIONAL_GROUPS
        spec = FUNCTIONAL_GROUPS["aliphatic_hydroxyl"]
        # Attaches to a CT (sp3) carbon, not CA.
        assert ("O", "OH") in spec["atoms"]
        assert spec["connectivity"][0][1] == "CT"

    def test_aliphatic_hydroxyl_uses_no_new_opls_type(self):
        for _atom, code in FUNCTIONAL_GROUPS["aliphatic_hydroxyl"]["atoms"]:
            assert code in GROMACS_OPLS_TYPE_MAP

    def test_get_aliphatic_sites_finds_pendant_methyls(self):
        mol = _aromatic_flake_with_methyls(3)
        assigner = OxygenAssigner(seed=0)
        sites = assigner._get_aliphatic_sites(mol)
        assert len(sites) == 3
        for idx in sites:
            a = mol.GetAtomWithIdx(idx)
            assert a.GetAtomicNum() == 6 and not a.GetIsAromatic()

    def test_get_aliphatic_sites_empty_for_pure_aromatic(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        Chem.SanitizeMol(mol)
        assert OxygenAssigner(seed=0)._get_aliphatic_sites(mol) == []


# ---------------------------------------------------------------------------
# Placement in explicit-dict mode
# ---------------------------------------------------------------------------

class TestExplicitPlacement:
    def test_places_aliphatic_hydroxyl_on_sp3_carbon(self):
        mol = _aromatic_flake_with_methyls(2)
        assigner = OxygenAssigner(seed=1)
        out, comp = assigner.assign_oxygens(
            mol, 0.0, functional_group_preference={"aliphatic_hydroxyl": 2}
        )
        assert comp.placed_counts.get("aliphatic_hydroxyl", 0) == 2
        assert _count_oxygens(out) == 2
        assert len(_aliphatic_hydroxyl_carbons(out)) == 2

    def test_cannot_exceed_aliphatic_site_count(self):
        mol = _aromatic_flake_with_methyls(2)
        assigner = OxygenAssigner(seed=1)
        out, comp = assigner.assign_oxygens(
            mol, 0.0, functional_group_preference={"aliphatic_hydroxyl": 5}
        )
        # Only 2 sp3 sites exist; placement stops there rather than crashing.
        assert comp.placed_counts.get("aliphatic_hydroxyl", 0) == 2

    def test_typed_and_neutral(self):
        mol = _aromatic_flake_with_methyls(2)
        out, _ = OxygenAssigner(seed=1).assign_oxygens(
            mol, 0.0, functional_group_preference={"aliphatic_hydroxyl": 2}
        )
        out = Chem.AddHs(out)
        Chem.SanitizeMol(out)
        types = AtomTyper().assign_atom_types(out)
        assert all(t in GROMACS_OPLS_TYPE_MAP for t in types.values())
        charges = ChargeAssigner().assign_charges(out, types)
        assert abs(sum(charges.values())) < 1e-6


# ---------------------------------------------------------------------------
# O_C_ratio-mode spill: the ceiling break
# ---------------------------------------------------------------------------

class TestRatioModeSpill:
    def test_spills_to_aliphatic_when_aromatic_exhausted(self):
        """With few aromatic sites and a high O target, the shortfall lands on
        the sp3 carbons instead of being silently dropped."""
        mol = _aromatic_flake_with_methyls(8)  # pyrene (16C) + 8 sp3 = 24 C
        assigner = OxygenAssigner(seed=2)
        # Target ~10 O on 24 C; pyrene has far fewer than 10 free aromatic edges
        # once 8 are consumed by the methyls.
        out, comp = assigner.assign_oxygens(mol, 10 / 24)
        assert comp.placed_counts.get("aliphatic_hydroxyl", 0) > 0, \
            "expected spill onto aliphatic carbons"
        # The aliphatic hydroxyls are really there.
        assert len(_aliphatic_hydroxyl_carbons(out)) == \
            comp.placed_counts["aliphatic_hydroxyl"]

    def test_opt_out_reproduces_aromatic_only(self):
        mol = _aromatic_flake_with_methyls(8)
        assigner = OxygenAssigner(seed=2)
        out, comp = assigner.assign_oxygens(
            mol, 10 / 24, allow_aliphatic_oxygen=False
        )
        assert comp.placed_counts.get("aliphatic_hydroxyl", 0) == 0
        assert _aliphatic_hydroxyl_carbons(out) == []

    def test_pure_aromatic_unaffected_by_flag(self):
        """A structure with no sp3 carbons behaves identically either way."""
        mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")  # pyrene
        Chem.SanitizeMol(mol)
        on, _ = OxygenAssigner(seed=3).assign_oxygens(mol, 0.2)
        off, _ = OxygenAssigner(seed=3).assign_oxygens(
            mol, 0.2, allow_aliphatic_oxygen=False
        )
        assert _count_oxygens(on) == _count_oxygens(off)


# ---------------------------------------------------------------------------
# End-to-end: the sweep point that hit the composition ceiling
# ---------------------------------------------------------------------------

class TestHighOxygenCharComposition:
    """hardwood 300 °C (O/C 0.245) placed 1 oxygen before; must now reach target.

    NOTE: this asserts *composition*, not strict_pass.  These structures also
    form intramolecular O-H...O contacts that the geometry validator reports as
    steric clashes until the H-bond-aware clash fix (PR #30) is present; that is
    a separate change.  Here we prove only that the oxygen can now be placed.
    """

    @pytest.mark.parametrize("feedstock,T,target", [
        ("hardwood", 300, 0.245),
        ("grass", 300, 0.270),
        ("manure", 300, 0.307),
    ])
    def test_reaches_oc_target(self, feedstock, T, target):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        cfg = GeneratorConfig(
            temperature=T, feedstock=feedstock, molecule_name="x",
            seed=0, strict=False,
        )
        _, _, comp = BiocharGenerator(cfg).generate()
        # Within the 10% O/C tolerance the validator uses.
        assert comp.O_C_ratio == pytest.approx(cfg.O_C_ratio, rel=0.10)

    def test_opt_out_still_short(self):
        """With the spill disabled, hardwood 300 stays oxygen-starved -- proving
        the aliphatic placement is what closes the gap."""
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        cfg = GeneratorConfig(
            temperature=300, feedstock="hardwood", molecule_name="x",
            seed=0, strict=False, allow_aliphatic_oxygen=False,
        )
        _, _, comp = BiocharGenerator(cfg).generate()
        assert comp.O_C_ratio < 0.5 * cfg.O_C_ratio
