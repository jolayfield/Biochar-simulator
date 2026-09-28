"""
Tests for hex-lattice crowding relief (rqm/geometry-embedding.md, rq-8e5f5dc1).

The hex-lattice path places ring carbons on an exact flat lattice and every
substituent radially. At a cove -- two edge carbons 2.46 A apart, four bonds
apart -- both hydrogens point into the same vacant lattice site. The fixture is
benzo[c]phenanthrene, the smallest PAH with a cove, built exactly the way the
flat path builds a sheet: carbons on the lattice, each H 1.09 A out along the
bisector of its carbon's two ring neighbours.
"""

import re

import numpy as np
import pytest

from rdkit import Chem

import biochar.pipeline.geometry_3d as geometry_3d
from biochar.constants import LATTICE_MAX_SHEET_RMS
from biochar.pipeline.biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    ValidationError,
)
from biochar.pipeline.geometry_3d import (
    CoordinateGenerator,
    GeometryValidator,
    _clash_pairs,
    _tilt_crowded_substituents,
)

LATTICE_BOND = 1.42


def _flat_pah_on_lattice(ring_centres):
    """A PAH whose rings are hexagons centred at *ring_centres* (2D), flat at
    z = 0, with radially placed hydrogens -- the flat path's placement rule."""
    vertices = []

    def vertex_index(p):
        for k, q in enumerate(vertices):
            if np.linalg.norm(p - q) < 1e-3:
                return k
        vertices.append(p)
        return len(vertices) - 1

    rings = []
    for cx, cy in ring_centres:
        ring = []
        for k in range(6):
            a = np.radians(30.0 + 60.0 * k)
            ring.append(vertex_index(np.array([cx, cy]) + LATTICE_BOND * np.array([np.cos(a), np.sin(a)])))
        rings.append(ring)

    mol = Chem.RWMol()
    for _ in vertices:
        atom = Chem.Atom(6)
        atom.SetIsAromatic(True)
        mol.AddAtom(atom)
    bonded = set()
    for ring in rings:
        for k in range(6):
            i, j = sorted((ring[k], ring[(k + 1) % 6]))
            if (i, j) not in bonded:
                bonded.add((i, j))
                mol.AddBond(i, j, Chem.BondType.AROMATIC)

    coords = [np.array([p[0], p[1], 0.0]) for p in vertices]
    for c in range(len(vertices)):
        nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors()]
        if len(nbrs) != 2:
            continue
        h = mol.AddAtom(Chem.Atom(1))
        mol.AddBond(c, h, Chem.BondType.SINGLE)
        out = coords[c] - np.mean([coords[n] for n in nbrs], axis=0)
        coords.append(coords[c] + 1.09 * out / np.linalg.norm(out))

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    coords = np.array(coords)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(coords):
        conformer.SetAtomPosition(i, xyz.tolist())
    mol.AddConformer(conformer, assignId=True)
    return mol, coords


def _benzo_c_phenanthrene():
    """[4]helicene: four rings, each fused one step further round the same way."""
    step = LATTICE_BOND * np.sqrt(3.0)
    centres, p = [np.zeros(2)], np.zeros(2)
    for turn in (0.0, 60.0, 120.0):
        p = p + step * np.array([np.cos(np.radians(turn)), np.sin(np.radians(turn))])
        centres.append(p)
    return _flat_pah_on_lattice(centres)


def _cove_carbons(mol, coords):
    """The two edge carbons 2.46 A apart that are four bonds apart."""
    dm = Chem.GetDistanceMatrix(mol)
    carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    for i in carbons:
        for j in carbons:
            if i < j and dm[i][j] == 4 and abs(np.linalg.norm(coords[i] - coords[j]) - 2.46) < 0.02:
                return i, j
    raise AssertionError("fixture has no cove")


def _sheet_rms(mol, coords):
    ring = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()]
    centred = coords[ring] - coords[ring].mean(axis=0)
    normal = np.linalg.svd(centred, full_matrices=False)[2][-1]
    return float(np.sqrt(np.mean((centred @ normal) ** 2)))


@pytest.fixture
def cove():
    mol, coords = _benzo_c_phenanthrene()
    assert mol.GetNumAtoms() == 30  # C18H12
    assert _clash_pairs(mol, coords), "the flat placement must crowd the cove"
    return mol, coords


class TestCoveRelief:
    # rq-57d3fc46
    def test_cove_substituents_are_pulled_apart(self, cove):
        mol, coords = cove
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert _clash_pairs(mol, relieved) == []

    # rq-da5ef2fc
    def test_cove_opens_rather_than_only_moving_its_hydrogens(self, cove):
        mol, coords = cove
        i, j = _cove_carbons(mol, coords)
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.linalg.norm(relieved[i] - relieved[j]) > 2.46 + 0.1

    # rq-13532295
    def test_relieved_sheet_stays_flat_and_keeps_its_bonds(self, cove):
        mol, coords = cove
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert _sheet_rms(mol, relieved) <= LATTICE_MAX_SHEET_RMS
        for bond in mol.GetBonds():
            if bond.GetIsAromatic():
                length = np.linalg.norm(
                    relieved[bond.GetBeginAtomIdx()] - relieved[bond.GetEndAtomIdx()]
                )
                assert length == pytest.approx(LATTICE_BOND, abs=0.05)

    def test_relief_is_deterministic(self, cove):
        mol, coords = cove
        a = CoordinateGenerator(seed=3).relieve_lattice_crowding(Chem.Mol(mol), coords)
        b = CoordinateGenerator(seed=3).relieve_lattice_crowding(Chem.Mol(mol), coords)
        assert np.array_equal(a, b)

    def test_relieved_coordinates_are_carried_by_the_molecule(self, cove):
        mol, coords = cove
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.allclose(mol.GetConformer().GetPositions(), relieved)


class TestReliefNeverMakesThingsWorse:
    # rq-ac64435a
    def test_a_worse_relaxation_is_discarded(self, cove, monkeypatch):
        mol, coords = cove

        def wrecked(mol_, coords_, seed=None):
            broken = coords_.copy()
            broken[0] += np.array([3.0, 0.0, 0.0])  # tear a ring bond
            return broken

        monkeypatch.setattr(geometry_3d, "_tethered_lattice_relax", wrecked)
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.array_equal(relieved, _tilt_crowded_substituents(mol, coords))

    def test_a_relaxation_that_cannot_be_set_up_falls_back_to_the_tilt(self, cove, monkeypatch):
        mol, coords = cove
        monkeypatch.setattr(geometry_3d, "_tethered_lattice_relax", lambda *a, **k: None)
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.array_equal(relieved, _tilt_crowded_substituents(mol, coords))

    # rq-37a3d53e
    def test_a_sheet_with_nothing_to_relieve_is_left_alone(self):
        # Naphthalene has no cove: its flat lattice placement is already clean.
        mol, coords = _flat_pah_on_lattice([(0.0, 0.0), (LATTICE_BOND * np.sqrt(3.0), 0.0)])
        assert _clash_pairs(mol, coords) == []
        relieved = CoordinateGenerator(seed=0).relieve_lattice_crowding(mol, coords)
        assert np.array_equal(relieved, coords)


class TestOneDefinitionOfAClash:
    # rq-5079b4df -- the validator and the relief share one clash definition
    # rq-4eb90e63
    def test_validator_reports_exactly_the_clash_finders_pairs(self, cove):
        mol, coords = cove
        reported = [
            tuple(int(x) for x in re.search(r"atoms (\d+) and (\d+)", e).groups())
            for e in GeometryValidator._check_steric_clashes(mol, coords)
        ]
        assert reported == [(i, j) for i, j, _d, _f in _clash_pairs(mol, coords)]
        assert reported, "fixture must have clashes for this to mean anything"


def _generate_on_hex_lattice(monkeypatch, **config):
    """Generate a structure, asserting it really took the relief path."""
    calls = {"n": 0}
    original = CoordinateGenerator.relieve_lattice_crowding

    def spy(self, mol, coords):
        calls["n"] += 1
        return original(self, mol, coords)

    monkeypatch.setattr(CoordinateGenerator, "relieve_lattice_crowding", spy)
    result = BiocharGenerator(GeneratorConfig(**config)).generate()
    assert calls["n"] == 1, "expected the hex-lattice relief path"
    return result


class TestHydroxylPlacement:
    # rq-7515b128 -- a hydroxyl hydrogen keeps its bond angle
    # rq-653a5fd2
    def test_hydroxyls_on_a_large_sheet_are_bent(self, monkeypatch):
        mol, coords, _comp = _generate_on_hex_lattice(
            monkeypatch, target_num_carbons=100, H_C_ratio=0.35, O_C_ratio=0.1,
            seed=2, strict=False,
        )
        angles = []
        for o in mol.GetAtoms():
            if o.GetAtomicNum() != 8:
                continue
            hs = [n.GetIdx() for n in o.GetNeighbors() if n.GetAtomicNum() == 1]
            heavy = [n.GetIdx() for n in o.GetNeighbors() if n.GetAtomicNum() != 1]
            if len(hs) == 1 and len(heavy) == 1:
                u = coords[hs[0]] - coords[o.GetIdx()]
                v = coords[heavy[0]] - coords[o.GetIdx()]
                cos = u @ v / (np.linalg.norm(u) * np.linalg.norm(v))
                angles.append(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))
        assert angles, "expected hydroxyl groups on this structure"
        assert all(100.0 <= a <= 120.0 for a in angles), sorted(angles)


class TestStrictValidationOnLargeSheets:
    # rq-9a59d62c
    @pytest.mark.slow
    @pytest.mark.parametrize("o_c", [0.0, 0.1, 0.2])
    def test_no_steric_clash_failures_at_100_carbons(self, monkeypatch, o_c):
        # H/C 0.35 is reachable at 100 C, so composition cannot be what fails;
        # before relief every one of these seeds failed on a clash.
        for seed in range(1, 6):
            try:
                _generate_on_hex_lattice(
                    monkeypatch, target_num_carbons=100, H_C_ratio=0.35,
                    O_C_ratio=o_c, seed=seed,
                )
            except ValidationError as e:
                assert "Steric clash" not in str(e), f"seed {seed}: {str(e)[:300]}"
