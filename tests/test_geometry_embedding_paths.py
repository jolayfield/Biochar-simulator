"""
Tests for embedding-path selection, the bounded embedding retry, clash-
resolution routing, and clash-severity reporting -- the rqm/geometry-
embedding.md scenarios not already defended in tests/test_hbond_clash.py and
tests/test_clash_tolerance.py.

Path selection (ETKDG vs. the flat 2D-first/hex-lattice path) has no public
flag of its own except `CoordinateGenerator.used_hex_lattice`, which is only
set inside the 2D-first branch. To tell the two paths apart from outside, most
tests here spy on `AllChem.EmbedMolecule` -- it is called only by the ETKDG
loop; the 2D-first path builds its conformer by hand and never calls it.
"""

import numpy as np
import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig
from biochar.pipeline.carbon_skeleton import PAHAssembler
from biochar.pipeline.geometry_3d import ClashResolver, CoordinateGenerator, GeometryValidator


def _naphthalene_with_hs():
    """10 heavy atoms -- well under the 80-heavy-atom threshold."""
    return Chem.AddHs(Chem.MolFromSmiles("c1ccc2ccccc2c1"))


def _large_hex_lattice_mol(target_num_carbons=100, seed=0):
    """A skeleton built through carbon_skeleton.py, which tags atoms with the
    init_x/init_y hex-lattice positions that geometry_3d.py reads. Targets
    above 40 carbons take the iterative-fusion path (see rqm/ARCHITECTURE.md),
    which always produces these tags -- so this is >80 heavy atoms *and* on
    the hex-lattice path, which is the realistic case for a structure this
    size.
    """
    return PAHAssembler(seed=seed).generate(target_num_carbons=target_num_carbons).mol


class TestEmbeddingPathSelection:
    # rq-bb1fbb12
    def test_small_molecule_uses_the_distance_geometry_embedder(self, monkeypatch):
        calls = {"n": 0}
        orig = AllChem.EmbedMolecule

        def spy(mol, params=None, **kw):
            calls["n"] += 1
            return orig(mol, params, **kw) if params is not None else orig(mol, **kw)

        monkeypatch.setattr(AllChem, "EmbedMolecule", spy)

        mol = _naphthalene_with_hs()
        assert mol.GetNumHeavyAtoms() <= 80

        gen = CoordinateGenerator(seed=0)
        mol2, coords = gen.generate_3d_coordinates(mol)

        assert calls["n"] > 0, "ETKDG's EmbedMolecule was never called"
        assert not gen.used_hex_lattice
        assert np.all(np.isfinite(coords))
        assert coords.shape[0] == mol2.GetNumAtoms()

    # rq-265ec385
    def test_large_flat_sheet_uses_the_2d_first_path(self, monkeypatch):
        calls = {"n": 0}

        def spy(mol, params=None, **kw):
            calls["n"] += 1
            return -1  # unused: EmbedMolecule must not be reached at all

        monkeypatch.setattr(AllChem, "EmbedMolecule", spy)

        mol = _large_hex_lattice_mol(target_num_carbons=100)
        assert mol.GetNumAtoms() > 80  # heavy-atom-only mol at this stage

        gen = CoordinateGenerator(seed=0)
        mol2, coords = gen.generate_3d_coordinates(mol)

        assert calls["n"] == 0, (
            "the 2D-first path builds its own conformer and must never call "
            "the ETKDG embedder"
        )
        # "Not folded": every heavy atom sits in a single global plane, not
        # merely locally planar per ring (which force_aromatic_planarity
        # enforces on both paths and so cannot distinguish them).
        heavy_idx = [a.GetIdx() for a in mol2.GetAtoms() if a.GetAtomicNum() != 1]
        z = coords[heavy_idx, 2]
        assert (z.max() - z.min()) < 0.5, "sheet is not flat -- looks folded"

    # rq-c4e683c4
    def test_embedding_retries_a_bounded_number_of_times(self, monkeypatch):
        calls = {"n": 0}

        def always_fail(mol, params=None, **kw):
            calls["n"] += 1
            return -1

        monkeypatch.setattr(AllChem, "EmbedMolecule", always_fail)

        mol = _naphthalene_with_hs()
        gen = CoordinateGenerator(seed=0)
        # Must not hang or raise even though every embedder fails outright.
        mol2, coords = gen.generate_3d_coordinates(mol, max_embedding_attempts=3)

        # 3 attempts x 3 embedders (ETKDGv3, ETKDGv2, EmbedParameters) is the
        # loosest possible bound; the current implementation stops re-trying
        # an embedder once it has failed once, so the real count is lower.
        assert calls["n"] <= 3 * 3
        assert np.all(np.isfinite(coords)), "must still produce a 2D fallback"
        assert coords.shape[0] == mol2.GetNumAtoms()


class TestClashResolutionRouting:
    # rq-49e2ed82
    def test_hex_lattice_structures_keep_their_coordinates(self, monkeypatch):
        calls = {"n": 0}
        orig = ClashResolver.resolve_clashes

        def spy(*a, **k):
            calls["n"] += 1
            return orig(*a, **k)

        monkeypatch.setattr(ClashResolver, "resolve_clashes", staticmethod(spy))

        # >40 carbons takes the hex-lattice growth path (rqm/ARCHITECTURE.md);
        # its edge/peri contacts routinely trip the generic clash floor (see
        # geometry-embedding.md), which is exactly the case this skip exists
        # for -- so this is not a quiet-molecule accident.
        cfg = GeneratorConfig(
            target_num_carbons=100, H_C_ratio=0.5, O_C_ratio=0.15,
            molecule_name="lg", seed=0, strict=False,
        )
        mol, coords, comp = BiocharGenerator(cfg).generate()

        assert calls["n"] == 0, "clash resolution must not run on the hex lattice"

        # "Ring lattice intact": aromatic C-C bonds sit at the hex-lattice
        # value (1.42 Angstrom), not displaced by a resolver that never ran.
        conf = mol.GetConformer()
        aromatic_cc = [
            b for b in mol.GetBonds()
            if b.GetIsAromatic()
            and b.GetBeginAtom().GetAtomicNum() == 6
            and b.GetEndAtom().GetAtomicNum() == 6
        ]
        assert aromatic_cc, "expected a fused aromatic sheet"
        for b in aromatic_cc[:20]:
            p1 = np.array(conf.GetAtomPosition(b.GetBeginAtomIdx()))
            p2 = np.array(conf.GetAtomPosition(b.GetEndAtomIdx()))
            length = np.linalg.norm(p1 - p2)
            assert length == pytest.approx(1.42, abs=0.05)

    # rq-3c400653
    def test_small_molecule_structures_still_get_clash_resolution(self, monkeypatch):
        calls = {"n": 0}
        orig_resolve = ClashResolver.resolve_clashes

        def spy_resolve(*a, **k):
            calls["n"] += 1
            return orig_resolve(*a, **k)

        monkeypatch.setattr(ClashResolver, "resolve_clashes", staticmethod(spy_resolve))

        # Force a clash deterministically rather than hoping a real small
        # structure happens to embed with one -- this isolates the routing
        # decision (used_hex_lattice False -> resolver is reachable) from
        # incidental embedding luck.
        orig_validate = GeometryValidator.validate_geometry

        def forced_clash(mol, coords):
            valid, errors = orig_validate(mol, coords)
            fake = (
                "Steric clash: atoms 0 and 1 (distance: 1.00, min: 2.00, "
                "severity: 1.00 Å, type: Other)"
            )
            return False, list(errors) + [fake]

        monkeypatch.setattr(GeometryValidator, "validate_geometry", staticmethod(forced_clash))

        cfg = GeneratorConfig(
            target_num_carbons=10, H_C_ratio=0.5, O_C_ratio=0.0,
            molecule_name="sm", seed=0, strict=False,
        )
        BiocharGenerator(cfg).generate()

        assert calls["n"] >= 1, (
            "a small (non-hex-lattice) structure with a reported clash must "
            "reach the resolver -- unlike the hex-lattice path, it is not "
            "gated off"
        )


class TestClashSeverityReporting:
    # rq-92fec87b
    def test_severity_distinguishes_a_marginal_contact_from_a_real_one(self):
        def _clash_message(distance):
            m = Chem.RWMol()
            m.AddAtom(Chem.Atom("C"))
            m.AddAtom(Chem.Atom("C"))
            coords = np.array([[0.0, 0.0, 0.0], [distance, 0.0, 0.0]], float)
            errors = GeometryValidator.validate_geometry(m.GetMol(), coords)[1]
            clashes = [e for e in errors if "Steric clash" in e]
            assert clashes, f"expected a clash at distance {distance}"
            return clashes[0]

        # 0.75 * (1.70 + 1.70) = 2.55 Angstrom is the generic C-C floor.
        marginal = _clash_message(2.55 - 0.10)   # just past the 0.05 Angstrom tolerance
        real = _clash_message(2.55 - 0.40)       # a genuine overlap

        for msg in (marginal, real):
            assert "severity:" in msg
            # The reported severity is a parseable float, not a placeholder.
            token = msg.split("severity:")[1].split("Å")[0].strip()
            float(token)

        marginal_severity = float(marginal.split("severity:")[1].split("Å")[0])
        real_severity = float(real.split("severity:")[1].split("Å")[0])
        assert real_severity > marginal_severity, (
            "a deeper overlap must report a larger severity, so a marginal "
            "contact can be told apart from a real one"
        )


class TestRefinementIsNotGatedOnClashes:
    """rqm/geometry-embedding.md: force-field refinement is unconditional on the
    non-hex-lattice path.

    Both scenarios read the same run. A structure that reports no clash gives the
    sharpest evidence available from outside: the clash resolver has nothing to do
    and is never invoked, so any refinement that still happens cannot have been
    gated on the clash report.
    """

    @staticmethod
    def _run_with_spies(monkeypatch):
        calls = {"relax": 0, "resolve": 0}

        orig_relax = CoordinateGenerator.validate_and_relax
        orig_resolve = ClashResolver.resolve_clashes

        def spy_relax(self, mol, coords, max_iterations=100):
            calls["relax"] += 1
            return orig_relax(self, mol, coords, max_iterations=max_iterations)

        def spy_resolve(mol, coords, **kw):
            calls["resolve"] += 1
            return orig_resolve(mol, coords, **kw)

        monkeypatch.setattr(CoordinateGenerator, "validate_and_relax", spy_relax)
        monkeypatch.setattr(ClashResolver, "resolve_clashes", staticmethod(spy_resolve))

        # 20 carbons with no oxygen: small enough for ETKDG (<= 80 heavy atoms),
        # and clean enough that no contact is reported.
        gen = BiocharGenerator(
            GeneratorConfig(target_num_carbons=20, O_C_ratio=0.0, strict=False, seed=1)
        )
        mol, coords, _ = gen.generate()

        errors = GeometryValidator.validate_geometry(mol, coords)[1]
        clashes = [e for e in errors if "Steric clash" in e]
        assert not clashes, (
            f"precondition failed -- this structure was chosen because it has no "
            f"clashes, but validation reported {len(clashes)}: {clashes[:3]}"
        )
        return calls

    # rq-acf97ed2
    def test_structure_with_no_reported_clashes_is_still_refined(self, monkeypatch):
        calls = self._run_with_spies(monkeypatch)
        assert calls["relax"] > 0, (
            "force-field refinement never ran on a clash-free structure -- "
            "refinement has been coupled to the clash report"
        )

    # rq-5658c4b6
    def test_clash_resolver_is_not_invoked_when_there_is_nothing_to_resolve(
        self, monkeypatch
    ):
        calls = self._run_with_spies(monkeypatch)
        assert calls["resolve"] == 0, (
            f"clash resolution ran {calls['resolve']} time(s) on a structure with "
            f"no reported clash"
        )


class TestEveryGeometryErrorIsReported:
    """rqm/geometry-embedding.md: a structure report names every geometry error."""

    @staticmethod
    def _mol_with_many_geometry_errors():
        """Naphthalene with six atoms collapsed onto one another.

        Produces far more than three geometry errors, so a report that truncates
        is distinguishable from one that does not.
        """
        mol = _naphthalene_with_hs()
        AllChem.EmbedMolecule(mol, randomSeed=7)
        coords = np.array(mol.GetConformer().GetPositions())
        for i in range(6):
            coords[i] = coords[0] + np.array([0.02 * i, 0.0, 0.0])
        return mol, coords

    # rq-1f2ce9fc
    def test_report_contains_every_geometry_error(self):
        from biochar.pipeline.validation import StructureValidator

        mol, coords = self._mol_with_many_geometry_errors()

        found = GeometryValidator.validate_geometry(mol, coords)[1]
        assert len(found) > 3, (
            "this fixture is meant to produce more than three geometry errors; "
            f"it produced {len(found)}"
        )

        report = StructureValidator.validate(mol, coords)
        reported = [e for e in report.errors if e in found]
        assert len(reported) == len(found), (
            f"geometry validation found {len(found)} errors but the structure "
            f"report carries {len(reported)}"
        )
