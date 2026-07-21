"""
Tests for hydrogen-bond-aware steric clash detection.

The generic clash floor is 0.75 × vdW-sum, which for an O/H pair is 2.04 Å --
squarely inside the physical H...A range of a hydrogen bond (~1.6-2.2 Å).  Every
intramolecular O-H...O contact between adjacent -OH groups was therefore
reported as a steric clash, and on high-oxygen chars (O/C >= ~0.2, e.g. 400 °C
softwood) such contacts are unavoidable -- so strict validation failed on every
seed and the sweep degraded to its fallback path.

These tests pin the reduced H-bond floor, the angle gate that keeps it from
excusing a genuine overlap, and the end-to-end 400 °C softwood regression.
"""

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from biochar.constants import (
    HBOND_MIN_H_ACCEPTOR_DISTANCE,
    HBOND_MIN_DHA_ANGLE_DEG,
)
from biochar.geometry_3d import (
    GeometryValidator,
    _clash_floor,
    _get_hbond_pairs,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _catechol_3d(seed: int = 0xF00D):
    """Benzene-1,2-diol: two adjacent -OH groups that form an intramolecular
    O-H...O hydrogen bond once the geometry is minimised."""
    mol = Chem.AddHs(Chem.MolFromSmiles("Oc1ccccc1O"))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    return mol, mol.GetConformer().GetPositions()


def _hydroxyl_indices(mol):
    """Return (donor_O, polar_H) index pairs for every -OH in *mol*."""
    out = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                out.append((atom.GetIdx(), nbr.GetIdx()))
    return out


def _place_hbond(distance: float, angle_deg: float) -> tuple:
    """Build a minimal water-like O-H...O system with an exact H...A distance
    and D-H...A angle, so the geometric gate can be probed directly.

    Layout: donor O at the origin, H at +x (1.0 Å), acceptor placed at
    *distance* from H, rotated by (180 - angle) from the D->H direction.
    """
    mol = Chem.RWMol()
    donor = mol.AddAtom(Chem.Atom(8))       # 0: donor O
    h = mol.AddAtom(Chem.Atom(1))           # 1: polar H
    acceptor = mol.AddAtom(Chem.Atom(8))    # 2: acceptor O
    h2 = mol.AddAtom(Chem.Atom(1))          # 3: H on the acceptor
    mol.AddBond(donor, h, Chem.BondType.SINGLE)
    mol.AddBond(acceptor, h2, Chem.BondType.SINGLE)
    mol = mol.GetMol()

    # Angle at H between (H->D) and (H->A).  H->D points along -x, so an
    # acceptor at angle theta from -x sits at H + distance*(cos, sin) rotated
    # accordingly.
    theta = np.radians(angle_deg)
    h_pos = np.array([1.0, 0.0, 0.0])
    hd_dir = np.array([-1.0, 0.0, 0.0])          # H -> D
    ha_dir = np.array([np.cos(theta) * hd_dir[0],
                       np.sin(theta), 0.0])
    coords = np.zeros((4, 3))
    coords[donor] = [0.0, 0.0, 0.0]
    coords[h] = h_pos
    coords[acceptor] = h_pos + distance * ha_dir
    # Park the acceptor's H far away so it cannot itself form a contact.
    coords[h2] = coords[acceptor] + np.array([0.0, 0.0, 5.0])
    return mol, coords, (donor, h, acceptor)


# ---------------------------------------------------------------------------
# _get_hbond_pairs / _clash_floor unit behaviour
# ---------------------------------------------------------------------------

class TestHBondDetection:
    def test_catechol_intramolecular_hbond_is_detected(self):
        """The adjacent -OH pair in catechol registers as a hydrogen bond."""
        mol, coords = _catechol_3d()
        pairs = _get_hbond_pairs(mol, coords)
        assert pairs, "expected at least one O-H...O pair in catechol"

        # Every detected pair must be a polar H paired with an O acceptor.
        for i, j in pairs:
            symbols = {mol.GetAtomWithIdx(i).GetSymbol(),
                       mol.GetAtomWithIdx(j).GetSymbol()}
            assert symbols == {"H", "O"}

    def test_hbond_pair_gets_reduced_floor(self):
        mol, coords = _catechol_3d()
        pairs = _get_hbond_pairs(mol, coords)
        i, j = next(iter(pairs))
        assert _clash_floor(mol, i, j, pairs) == HBOND_MIN_H_ACCEPTOR_DISTANCE

    def test_non_hbond_pair_keeps_vdw_floor(self):
        """A plain aromatic C/H pair is untouched by the H-bond rule."""
        mol, coords = _catechol_3d()
        pairs = _get_hbond_pairs(mol, coords)
        carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
        i, j = carbons[0], carbons[2]
        # 0.75 * (1.70 + 1.70)
        assert _clash_floor(mol, i, j, pairs) == pytest.approx(2.55)

    def test_no_polar_hydrogens_yields_no_pairs(self):
        """Benzene has no donors, so the helper short-circuits to empty."""
        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(mol, randomSeed=7)
        coords = mol.GetConformer().GetPositions()
        assert _get_hbond_pairs(mol, coords) == set()


class TestHBondAngleGate:
    """The angle gate must admit near-linear contacts and reject sideways ones."""

    def test_linear_contact_is_an_hbond(self):
        mol, coords, (donor, h, acceptor) = _place_hbond(distance=1.9, angle_deg=175.0)
        pairs = _get_hbond_pairs(mol, coords)
        assert (min(h, acceptor), max(h, acceptor)) in pairs

    def test_sideways_contact_is_not_excused(self):
        """An acceptor jammed into the side of the D-H bond is a real clash."""
        mol, coords, (donor, h, acceptor) = _place_hbond(distance=1.9, angle_deg=30.0)
        pairs = _get_hbond_pairs(mol, coords)
        assert (min(h, acceptor), max(h, acceptor)) not in pairs
        # ...and it therefore keeps the generic 2.04 Å O/H floor.
        assert _clash_floor(mol, h, acceptor, pairs) == pytest.approx(2.04)

    @pytest.mark.parametrize("angle", [HBOND_MIN_DHA_ANGLE_DEG + 1.0, 120.0, 180.0])
    def test_angles_at_or_above_gate_are_hbonds(self, angle):
        mol, coords, (donor, h, acceptor) = _place_hbond(distance=1.9, angle_deg=angle)
        pairs = _get_hbond_pairs(mol, coords)
        assert (min(h, acceptor), max(h, acceptor)) in pairs


class TestGenuineOverlapStillCaught:
    """The reduced floor must not become a blanket exemption."""

    def test_overlapping_donor_acceptor_is_still_a_clash(self):
        mol, coords, (donor, h, acceptor) = _place_hbond(distance=1.0, angle_deg=175.0)
        errors = GeometryValidator._check_steric_clashes(mol, coords)
        clashes = [e for e in errors if "Steric clash" in e]
        # At this separation the two oxygens are 2.0 Å apart and clash too, so
        # pick out the H...acceptor pair specifically.
        h_a = [e for e in clashes if f"atoms {min(h, acceptor)} and {max(h, acceptor)}" in e]
        assert h_a, f"an H...O contact of 1.0 Å must still be reported; got {clashes}"
        assert "type: H-bond" in h_a[0]

    def test_normal_hbond_distance_is_not_a_clash(self):
        mol, coords, (donor, h, acceptor) = _place_hbond(distance=1.9, angle_deg=175.0)
        errors = GeometryValidator._check_steric_clashes(mol, coords)
        assert not [e for e in errors if "Steric clash" in e]

    def test_catechol_reports_no_steric_clash(self):
        mol, coords = _catechol_3d()
        errors = GeometryValidator._check_steric_clashes(mol, coords)
        assert not [e for e in errors if "Steric clash" in e]


# ---------------------------------------------------------------------------
# End-to-end regression: the sweep point that drove everything to fallback
# ---------------------------------------------------------------------------

class TestHighOxygenCharRegression:
    """400 °C softwood (O/C ~= 0.205) must pass strict validation.

    Before H-bond-aware clash detection every seed failed on O-H...O contacts
    reported as steric clashes, so `sweep.build_point` fell back on all points.
    """

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_400C_softwood_builds_clash_free(self, seed):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        cfg = GeneratorConfig(
            temperature=400, feedstock="softwood",
            molecule_name="sw400", seed=seed, strict=False,
        )
        gen = BiocharGenerator(cfg)
        mol, coords, comp = gen.generate()

        errors = GeometryValidator.validate_geometry(mol, coords)[1]
        clashes = [e for e in errors if "Steric clash" in e]
        assert not clashes, f"unexpected clashes at seed {seed}: {clashes}"

        # The composition was always on target; guard against regressing it
        # while chasing the geometry.
        assert comp.H_C_ratio == pytest.approx(cfg.H_C_ratio, rel=0.10)
        assert comp.O_C_ratio == pytest.approx(cfg.O_C_ratio, rel=0.10)

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_400C_softwood_bond_lengths_are_sane(self, seed):
        """Geometry refinement must not be skipped just because there is no clash.

        The FF pass in `_generate_geometry` used to sit behind `if
        steric_clashes:`.  Every high-oxygen structure carried at least one
        false H-bond "clash", so the pass always ran and the coupling was
        invisible -- until H-bond-aware detection removed the clashes and left
        a 1.16 Å aromatic bond unrelaxed.
        """
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        cfg = GeneratorConfig(
            temperature=400, feedstock="softwood",
            molecule_name="sw400", seed=seed, strict=False,
        )
        mol, coords, _ = BiocharGenerator(cfg).generate()

        errors = GeometryValidator.validate_geometry(mol, coords)[1]
        bad_bonds = [e for e in errors if "bond length" in e]
        assert not bad_bonds, f"unrelaxed geometry at seed {seed}: {bad_bonds}"

    def test_400C_softwood_passes_strict_mode(self):
        """The full strict path -- what the sweep actually invokes."""
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        cfg = GeneratorConfig(
            temperature=400, feedstock="softwood",
            molecule_name="sw400", seed=0, strict=True,
        )
        # strict=True raises ValidationError if the structure does not validate.
        mol, coords, comp = BiocharGenerator(cfg).generate()
        assert mol is not None
