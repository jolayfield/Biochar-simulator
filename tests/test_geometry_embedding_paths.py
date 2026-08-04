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
from biochar.pipeline.geometry_3d import (
    ClashResolver,
    CoordinateGenerator,
    GeometryValidator,
    _kekulize_or_dearomatize,
)


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


class TestEmbeddingPreparationKeepsRingInformation:
    """rqm/geometry-embedding.md: the molecule handed on still knows its rings.

    `_kekulize_or_dearomatize` is allowed to fail -- a graph carrying one bad
    valence is still one the force field can be handed -- but `SanitizeMol`
    clears the computed properties, ring information among them, before it
    does any of its work. A failure partway through therefore used to hand
    back a molecule that knew less than the one that went in, and the crash
    landed three stages downstream in atom typing.
    """

    @staticmethod
    def _unkekulizable_with_a_bad_valence():
        """A non-kekulisable ring, and separately an atom over its valence.

        Both halves are needed and neither can be the other. Kekulisation has
        to fail so the de-aromatising branch runs at all, and the sanitisation
        in that branch has to fail so it raises after clearing the rings.

        They are kept apart because putting them on the same atom does not
        work: an over-bonded ring atom is a substituted one, which gives
        kekulisation the sp3 escape it needed and the branch is never reached.
        Kekulisation failing is also not enough on its own -- it perceives the
        rings on its own way to raising, and a sanitisation that then succeeds
        hands them back.

        The de-aromatising branch is reached by an odd ring of aromatic
        carbons each holding one pinned hydrogen, which has no kekule
        structure; its sanitisation is failed by a neutral nitrogen carrying
        four, which is the shape of the real failure -- a valence the pipeline
        wrote earlier and geometry is the first thing to sanitise.
        """
        rw = Chem.RWMol()
        carbons = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for a, b in zip(carbons, carbons[1:] + carbons[:1]):
            bond = rw.GetBondWithIdx(rw.AddBond(a, b, Chem.BondType.AROMATIC) - 1)
            bond.SetIsAromatic(True)
        for idx in carbons:
            atom = rw.GetAtomWithIdx(idx)
            atom.SetIsAromatic(True)
            atom.SetNoImplicit(True)
            rw.AddBond(idx, rw.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)

        nitrogen = rw.AddAtom(Chem.Atom(7))
        rw.GetAtomWithIdx(nitrogen).SetNoImplicit(True)
        for _ in range(4):
            rw.AddBond(nitrogen, rw.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)
        return rw.GetMol()

    # rq-a157df17
    def test_a_molecule_whose_sanitisation_failed_still_knows_its_rings(self):
        mol = self._unkekulizable_with_a_bad_valence()
        with pytest.raises(Chem.KekulizeException):
            Chem.Kekulize(Chem.RWMol(mol), clearAromaticFlags=False)
        with pytest.raises(Chem.AtomValenceException):
            Chem.SanitizeMol(Chem.RWMol(mol))

        prepared = _kekulize_or_dearomatize(mol)

        assert prepared.GetRingInfo().NumRings() == 1

    # rq-ac433b55
    def test_a_molecule_that_kekulises_also_knows_its_rings(self):
        prepared = _kekulize_or_dearomatize(_naphthalene_with_hs())

        assert prepared.GetRingInfo().NumRings() == 2


class TestEmbeddingReturnsTheMoleculeItWasGiven:
    """rqm/geometry-embedding.md: the bond-order-rewritten copy does not escape.

    Embedders and force fields need integer bond orders, and a fused biochar
    sheet frequently has no kekule structure, so the preparation step strips
    every aromatic bond to single. Those coordinates are a result; that
    chemistry is scaffolding.
    """

    @staticmethod
    def _unkekulizable_ring():
        """Five aromatic carbons, each holding one pinned explicit hydrogen.

        An odd number of π centres cannot be paired into double bonds, so
        kekulisation fails -- the same reason it fails on a real char sheet,
        where `HydrogenAssigner` pins the hydrogens that fix the parity. A bare
        `PAHAssembler` skeleton always kekulises, so it is no fixture for this:
        the failure only appears once the hydrogens are on.
        """
        rw = Chem.RWMol()
        carbons = [rw.AddAtom(Chem.Atom(6)) for _ in range(5)]
        for a, b in zip(carbons, carbons[1:] + carbons[:1]):
            bond = rw.GetBondWithIdx(rw.AddBond(a, b, Chem.BondType.AROMATIC) - 1)
            bond.SetIsAromatic(True)
        for idx in carbons:
            atom = rw.GetAtomWithIdx(idx)
            atom.SetIsAromatic(True)
            atom.SetNoImplicit(True)
            rw.AddBond(idx, rw.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)
        mol = rw.GetMol()
        with pytest.raises(Chem.KekulizeException):
            Chem.Kekulize(Chem.RWMol(mol), clearAromaticFlags=False)
        return mol

    # rq-d5b18b3a
    def test_a_structure_that_will_not_kekulise_comes_back_aromatic(self):
        mol = self._unkekulizable_ring()
        aromatic_bonds_before = sum(
            1 for b in mol.GetBonds()
            if b.GetBondType() == Chem.BondType.AROMATIC
        )
        assert aromatic_bonds_before > 0

        out, _ = CoordinateGenerator(seed=0).generate_3d_coordinates(mol)

        aromatic_bonds_after = sum(
            1 for b in out.GetBonds()
            if b.GetBondType() == Chem.BondType.AROMATIC
        )
        assert aromatic_bonds_after == aromatic_bonds_before, (
            "coordinate generation handed back the de-aromatised working copy: "
            f"{aromatic_bonds_before} aromatic bonds went in, "
            f"{aromatic_bonds_after} came out"
        )
        assert any(a.GetIsAromatic() for a in out.GetAtoms())

    # rq-79908020
    def test_the_returned_molecule_carries_the_coordinates(self):
        mol = _naphthalene_with_hs()

        out, coords = CoordinateGenerator(seed=0).generate_3d_coordinates(mol)

        conf = np.array(out.GetConformer(0).GetPositions())
        assert conf.shape == coords.shape
        assert np.allclose(conf, coords, atol=1e-6) or np.allclose(
            conf, coords, atol=0.6
        ), "the returned conformer is not the coordinates that were returned"

    # rq-592b8c66
    def test_refinement_relaxes_an_aromatic_molecule(self):
        """MMFF and UFF both refuse aromatic bond orders, and both refusals are
        swallowed -- so a refinement that parametrised the molecule it was
        handed would return the coordinates untouched and refined-looking.
        """
        mol = _naphthalene_with_hs()
        AllChem.EmbedMolecule(mol, randomSeed=11)
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
        assert any(
            b.GetBondType() == Chem.BondType.AROMATIC for b in mol.GetBonds()
        )

        coords = np.array(mol.GetConformer().GetPositions())
        strained = coords.copy()
        strained[0] += np.array([0.4, 0.25, 0.15])

        relaxed, _ = CoordinateGenerator(seed=0).validate_and_relax(
            mol, strained, max_iterations=200
        )

        assert not np.allclose(relaxed, strained), (
            "refinement left the coordinates exactly as it found them"
        )


class TestValidationDoesNotRewriteWhatItValidates:
    """rqm/geometry-embedding.md: structure validation is a read."""

    # rq-75c9724b
    def test_validating_an_unkekulizable_molecule_leaves_its_bonds_alone(self):
        from biochar.pipeline.validation import ChemicalFeasibilityValidator

        mol = TestEmbeddingReturnsTheMoleculeItWasGiven._unkekulizable_ring()
        before = [b.GetBondType() for b in mol.GetBonds()]
        aromatic_before = [a.GetIsAromatic() for a in mol.GetAtoms()]

        ChemicalFeasibilityValidator.validate(mol)

        assert [b.GetBondType() for b in mol.GetBonds()] == before, (
            "validation rewrote the bond types of the molecule it was checking"
        )
        assert [a.GetIsAromatic() for a in mol.GetAtoms()] == aromatic_before

    # rq-8ebaed02
    def test_a_molecule_rdkit_rejects_is_still_reported(self):
        from biochar.pipeline.validation import ChemicalFeasibilityValidator

        rw = Chem.RWMol()
        carbon = rw.AddAtom(Chem.Atom(6))
        for _ in range(5):
            rw.AddBond(carbon, rw.AddAtom(Chem.Atom(1)), Chem.BondType.SINGLE)
        report = ChemicalFeasibilityValidator.validate(rw.GetMol())

        assert any("sanitization" in w for w in report.warnings), report.warnings
