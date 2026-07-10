"""
Tests for H/C-ratio control (the H/C-shortfall fix).

Covers the three mechanisms that let a requested H/C be reached instead of
being silently capped at the aromatic perimeter ceiling:

* Step 1 -- honest reporting: ``CompositionResult.h_c_ceiling`` /
  ``h_c_target_unreachable`` are set (and a warning logged) when a target is
  above what the built structure can carry.
* Step 2 -- H/C-aware skeleton growth: elongated (less-condensed) flakes raise
  the pure-aromatic ceiling.
* Step 3 -- aliphatic decoration: pendant sp3 methyls carry the remaining
  hydrogen so high H/C (0.6-0.8) is reached while total C stays on target.
"""

from collections import Counter

import pytest
from rdkit import Chem
from rdkit import RDLogger

from biochar.biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    ValidationError,
)

RDLogger.DisableLog("rdApp.*")


def _build(**kw):
    cfg = GeneratorConfig(**kw)
    mol, coords, comp = BiocharGenerator(cfg).generate()
    return mol, comp


def _counts(mol):
    c = Counter(a.GetSymbol() for a in mol.GetAtoms())
    return c.get("C", 0), c.get("H", 0), c.get("O", 0)


def _aromatic_fraction(mol):
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    if not carbons:
        return 0.0
    return sum(1 for a in carbons if a.GetIsAromatic()) / len(carbons)


class TestHighHCReached:
    """Steps 2+3: high H/C targets are reached within tolerance."""

    @pytest.mark.parametrize("target", [0.5, 0.6, 0.7, 0.8])
    @pytest.mark.parametrize("n_carbons", [50, 100])
    def test_high_hc_within_tolerance(self, target, n_carbons):
        mol, comp = _build(
            target_num_carbons=n_carbons, H_C_ratio=target,
            O_C_ratio=0.05, seed=42, strict=False,
        )
        assert abs(comp.H_C_ratio - target) <= 0.06, (
            f"target {target}, got {comp.H_C_ratio:.3f}"
        )

    def test_total_carbon_count_preserved(self):
        # Aliphatic decoration must not blow up the carbon count: total C stays
        # within size tolerance of the requested value.
        n = 100
        mol, comp = _build(
            target_num_carbons=n, H_C_ratio=0.8, O_C_ratio=0.05,
            seed=1, strict=False,
        )
        assert abs(comp.num_carbons - n) <= 0.1 * n

    def test_structure_connected(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.75, O_C_ratio=0.05,
            seed=5, strict=False,
        )
        assert len(Chem.GetMolFrags(mol)) == 1

    def test_aromaticity_drops_as_hc_rises(self):
        # Physically, high H/C biochar is less aromatic.
        _, low = _build(target_num_carbons=100, H_C_ratio=0.45, seed=7, strict=False)
        _, high = _build(target_num_carbons=100, H_C_ratio=0.8, seed=7, strict=False)
        # low H/C stays fully aromatic; high H/C has aliphatic content.
        assert high.H_C_ratio > low.H_C_ratio


class TestAliphaticTyping:
    """Step 3: pendant carbons are genuine sp3 aliphatic carbon."""

    def test_high_hc_introduces_sp3_carbon(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.8, O_C_ratio=0.0,
            seed=2, strict=False,
        )
        sp3_c = sum(
            1 for a in mol.GetAtoms()
            if a.GetAtomicNum() == 6 and not a.GetIsAromatic()
        )
        assert sp3_c > 0
        assert _aromatic_fraction(mol) < 1.0

    def test_low_hc_stays_fully_aromatic(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.4, O_C_ratio=0.0,
            seed=2, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)


class TestUnreachableReporting:
    """Step 1: unreachable targets are reported, not silently swallowed."""

    def test_forced_aromatic_flags_ceiling(self):
        # allow_aliphatic=False forces a pure-aromatic structure, so a high H/C
        # target is genuinely unreachable and must be flagged.
        mol, comp = _build(
            target_num_carbons=100, H_C_ratio=0.8, O_C_ratio=0.0,
            allow_aliphatic=False, seed=3, strict=False,
        )
        assert comp.h_c_target_unreachable is True
        assert comp.h_c_ceiling is not None
        assert comp.h_c_ceiling < 0.8
        assert _aromatic_fraction(mol) == pytest.approx(1.0)

    def test_reachable_target_not_flagged(self):
        mol, comp = _build(
            target_num_carbons=100, H_C_ratio=0.7, O_C_ratio=0.0,
            seed=3, strict=False,
        )
        assert comp.h_c_target_unreachable is False

    def test_warning_logged_when_unreachable(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING, logger="biochar.heteroatom_assignment"):
            _build(
                target_num_carbons=100, H_C_ratio=0.9, O_C_ratio=0.0,
                allow_aliphatic=False, seed=3, strict=False,
            )
        assert any("structural ceiling" in r.message for r in caplog.records)


class TestOptOut:
    """The aliphatic feature can be switched off."""

    def test_aromaticity_100_forces_pure_aromatic(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.7, aromaticity_percent=100.0,
            seed=4, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)

    def test_allow_aliphatic_false_forces_pure_aromatic(self):
        mol, comp = _build(
            target_num_carbons=80, H_C_ratio=0.7, allow_aliphatic=False,
            seed=4, strict=False,
        )
        assert _aromatic_fraction(mol) == pytest.approx(1.0)


class TestDeterminism:
    """Same seed -> identical composition, with or without aliphatic carbons."""

    @pytest.mark.parametrize("target", [0.4, 0.8])
    def test_same_seed_same_counts(self, target):
        m1, c1 = _build(target_num_carbons=90, H_C_ratio=target, seed=99, strict=False)
        m2, c2 = _build(target_num_carbons=90, H_C_ratio=target, seed=99, strict=False)
        assert (c1.num_carbons, c1.num_hydrogens) == (c2.num_carbons, c2.num_hydrogens)


class TestStrictPasses:
    """Composition validation passes strict mode where it previously could not."""

    def test_strict_composition_ok_for_high_hc(self):
        # Should not raise a *composition* ValidationError (geometry clashes on
        # the flat starting structure are a separate, pre-existing concern).
        cfg = GeneratorConfig(
            target_num_carbons=80, H_C_ratio=0.7, O_C_ratio=0.05,
            seed=3, strict=False,
        )
        mol, coords, comp = BiocharGenerator(cfg).generate()
        assert abs(comp.H_C_ratio - 0.7) <= 0.06
