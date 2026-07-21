"""
Tests for the clash-severity tolerance.

Neither clash floor (0.75 x vdW-sum, or the 1.5 A hydrogen-bond floor) is
meaningful to 0.01 A, but a hard `distance < floor` knife-edge reports a contact
a hundredth of an angstrom inside the floor as a clash.  On the densest chars
(400 C manure, O/C 0.24) every retried seed carried one such boundary contact --
an H-bond 0.01 A shorter than 1.5 A, or a peri H...C 0.005 A inside the vdW floor
-- and strict validation failed on all of them even though GROMACS EM relaxes the
atom in its first step.

CLASH_SEVERITY_TOLERANCE (0.05 A) makes a contact a clash only when it is deeper
than the tolerance below the floor.  These tests pin that a boundary contact is
ignored, a genuine overlap is still caught, and the tolerance sits well below the
depth of a real clash.
"""

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from biochar.constants import (
    VDW_RADII,
    CLASH_SEVERITY_TOLERANCE,
    HBOND_MIN_H_ACCEPTOR_DISTANCE,
)
from biochar.geometry_3d import GeometryValidator, ClashResolver


def _two_atoms(sym_a, sym_b, distance):
    m = Chem.RWMol()
    m.AddAtom(Chem.Atom(sym_a))
    m.AddAtom(Chem.Atom(sym_b))
    coords = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]], float)
    return m.GetMol(), coords


def _clashes(mol, coords):
    return [e for e in GeometryValidator.validate_geometry(mol, coords)[1]
            if "Steric clash" in e]


class TestBoundaryContactIgnored:
    def test_contact_just_inside_floor_is_not_a_clash(self):
        floor = 0.75 * (VDW_RADII["O"] + VDW_RADII["O"])  # 2.28
        mol, coords = _two_atoms(8, 8, floor - 0.03)  # 0.03 < 0.05 tolerance
        assert _clashes(mol, coords) == []

    def test_contact_at_the_floor_is_not_a_clash(self):
        floor = 0.75 * (VDW_RADII["C"] + VDW_RADII["H"])
        mol, coords = _two_atoms(6, 1, floor)  # exactly at floor
        assert _clashes(mol, coords) == []

    def test_tolerance_value_is_below_a_real_overlap(self):
        # Genuine overlaps (two groups embedded on each other) sit 0.2-0.5 A in;
        # the tolerance must stay well under that.
        assert 0.0 < CLASH_SEVERITY_TOLERANCE < 0.1


class TestGenuineOverlapStillCaught:
    def test_overlap_beyond_tolerance_is_a_clash(self):
        floor = 0.75 * (VDW_RADII["O"] + VDW_RADII["O"])
        mol, coords = _two_atoms(8, 8, floor - 0.10)  # 0.10 > 0.05 tolerance
        assert len(_clashes(mol, coords)) == 1

    def test_deep_overlap_is_a_clash(self):
        mol, coords = _two_atoms(8, 8, 1.20)  # far inside the 2.28 floor
        assert len(_clashes(mol, coords)) == 1

    def test_boundary_at_exactly_tolerance_depth(self):
        """A contact exactly CLASH_SEVERITY_TOLERANCE below the floor is the
        edge of the band: `distance < floor - tol` is strict, so it is NOT a
        clash, and one hair deeper IS."""
        floor = 0.75 * (VDW_RADII["O"] + VDW_RADII["O"])
        mol, at = _two_atoms(8, 8, floor - CLASH_SEVERITY_TOLERANCE)
        assert _clashes(mol, at) == []
        mol2, deeper = _two_atoms(8, 8, floor - CLASH_SEVERITY_TOLERANCE - 0.01)
        assert len(_clashes(mol2, deeper)) == 1


class TestResolverSharesTolerance:
    def test_resolver_ignores_boundary_contact(self):
        """The resolver must use the same definition of a clash, or it would
        nudge an atom the validator already accepts."""
        floor = 0.75 * (VDW_RADII["O"] + VDW_RADII["O"])
        mol, coords = _two_atoms(8, 8, floor - 0.03)
        assert ClashResolver._detect_clashes(mol, coords) == []

    def test_resolver_still_detects_real_overlap(self):
        mol, coords = _two_atoms(8, 8, 1.20)
        assert ClashResolver._detect_clashes(mol, coords) == [(0, 1)]


class TestManure400Regression:
    """400 C manure (O/C 0.236) failed strict on seeds 0-2, each on a single
    boundary contact (severity 0.00-0.03 A).  It must now pass on seed 0.

    Requires the H-bond clash fix and the aliphatic-oxygen fix already on main;
    this asserts the tolerance closes the last gap.
    """

    def test_manure_400_passes_strict(self):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        cfg = GeneratorConfig(
            temperature=400, feedstock="manure", molecule_name="mn400",
            seed=0, strict=True,
        )
        mol, coords, comp = BiocharGenerator(cfg).generate()
        assert mol is not None
        assert comp.O_C_ratio == pytest.approx(cfg.O_C_ratio, rel=0.10)
