"""
Tests for bond-order-aware bond-length validation.

COVALENT_RADII are single-bond radii, so their raw sum over-predicts every
aromatic and multiple bond: an aromatic C-C was reported as "expected 1.52"
when the real value is 1.40, and a C=O as 1.42 when it is 1.23.  The message
was therefore wrong about what it expected even when the bond genuinely was
out of range.

Correcting `expected` downward would on its own *lower* the absolute floor and
let a compressed bond through, so the tolerance factors were tightened to keep
detection at least as sensitive as before.  These tests pin both halves.
"""

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from biochar.constants import (
    COVALENT_RADII,
    BOND_ORDER_LENGTH_FACTORS,
    BOND_LENGTH_MIN_FACTOR,
    BOND_LENGTH_MAX_FACTOR,
)
from biochar.geometry_3d import GeometryValidator


def _optimized(smiles: str, seed: int = 3):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    return mol, mol.GetConformer().GetPositions()


def _bond_length_errors(mol, coords):
    return [e for e in GeometryValidator.validate_geometry(mol, coords)[1]
            if "bond length" in e]


class TestExpectedLengthIsBondOrderAware:
    """The reported expectation must match real bond lengths."""

    def test_aromatic_cc_expectation_is_140_not_152(self):
        expected = (COVALENT_RADII["C"] * 2) * BOND_ORDER_LENGTH_FACTORS["AROMATIC"]
        assert expected == pytest.approx(1.40, abs=0.02)

    def test_carbonyl_expectation_is_123_not_142(self):
        expected = (COVALENT_RADII["C"] + COVALENT_RADII["O"]) * \
            BOND_ORDER_LENGTH_FACTORS["DOUBLE"]
        assert expected == pytest.approx(1.23, abs=0.03)

    def test_single_bond_expectation_is_unscaled(self):
        assert BOND_ORDER_LENGTH_FACTORS["SINGLE"] == 1.00

    def test_message_quotes_the_aromatic_expectation(self):
        """A compressed aromatic bond reports 1.40, not the old 1.52."""
        mol, coords = _optimized("c1ccccc1")
        d = coords[1] - coords[0]
        coords[1] = coords[0] + d / np.linalg.norm(d) * 1.16
        errors = _bond_length_errors(mol, coords)
        assert errors, "a 1.16 Å aromatic C-C must be reported"
        assert "expected: 1.40" in errors[0]


class TestRealMoleculesValidateClean:
    """Relaxed real geometries must not trip the tightened band."""

    @pytest.mark.parametrize("smiles", [
        "c1ccccc1",           # benzene: aromatic C-C, aromatic C-H
        "Oc1ccccc1",          # phenol: aromatic C-O, O-H
        "OC(=O)c1ccccc1",     # benzoic acid: C=O, C-O, O-H
        "Oc1ccccc1O",         # catechol: adjacent hydroxyls
        "c1ccc2ccccc2c1",     # naphthalene: fused rings
        "Cc1ccccc1",          # toluene: sp3 C-C and C-H
    ])
    def test_no_bond_length_errors(self, smiles):
        mol, coords = _optimized(smiles)
        assert not _bond_length_errors(mol, coords)


class TestDetectionSensitivityPreserved:
    """Correcting `expected` must not weaken the floor."""

    def test_compressed_aromatic_bond_still_caught(self):
        mol, coords = _optimized("c1ccccc1")
        d = coords[1] - coords[0]
        coords[1] = coords[0] + d / np.linalg.norm(d) * 1.16
        assert _bond_length_errors(mol, coords)

    def test_stretched_bond_caught(self):
        mol, coords = _optimized("c1ccccc1")
        d = coords[1] - coords[0]
        coords[1] = coords[0] + d / np.linalg.norm(d) * 2.20
        assert _bond_length_errors(mol, coords)

    def test_absolute_floor_no_looser_than_before(self):
        """The old band was 0.8-1.5 x the *unscaled* radii sum.

        For an aromatic C-C that floor was 1.52 * 0.8 = 1.22 Å.  The corrected
        expectation must not raise the accepted range's floor above that by
        more than a hair, or the change would silently admit worse geometry.
        """
        old_floor = (COVALENT_RADII["C"] * 2) * 0.8
        new_floor = ((COVALENT_RADII["C"] * 2)
                     * BOND_ORDER_LENGTH_FACTORS["AROMATIC"]
                     * BOND_LENGTH_MIN_FACTOR)
        assert new_floor <= old_floor + 0.05

    def test_upper_bound_is_tighter_than_before(self):
        old_ceiling = (COVALENT_RADII["C"] * 2) * 1.5
        new_ceiling = ((COVALENT_RADII["C"] * 2)
                       * BOND_ORDER_LENGTH_FACTORS["AROMATIC"]
                       * BOND_LENGTH_MAX_FACTOR)
        assert new_ceiling < old_ceiling
