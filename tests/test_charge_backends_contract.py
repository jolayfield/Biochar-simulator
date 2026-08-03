"""
Tests for the rqm/charge-backends.md scenarios that no existing test defends.

The arithmetic half of that document -- charge conservation through the CM1A
mapping, the sum after scaling, the refusal of a missing MOPAC binary -- is
back-annotated onto tests/test_qm_charges.py and tests/test_ml_charges.py.

What is left is provenance: whether a structure's charges say where they came
from. Three ways they currently do not. A fallback model substitutes itself
with only a log line; a pickle written by a different scikit-learn is flagged
by scikit-learn rather than in terms of the charges; and an empirical factor
fitted on neutral liquids is applied to ions without comment, on the path the
pH-plus-ML refusal explicitly recommends.

Every scenario here passes. The five that carried xfail(strict=True) were
retired by XPASS when their fixes landed.
"""

import warnings

import pytest

from biochar.charges.ml_charges import MLChargeRefinement, _DEFAULT_MODEL_PATH
from biochar.charges.qm_charges import scale_and_neutralize

sklearn = pytest.importorskip("sklearn")


def _caught(fn, *args, **kwargs):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fn(*args, **kwargs)
    return result, [str(w.message) for w in caught]


# --------------------------------------------------------------------------- #
# The 1.14 factor on an ion
# --------------------------------------------------------------------------- #
class TestTheNeutralMoleculeFactor:
    CHARGES = [-0.6, -0.6, -0.6, 0.4, 0.4, 0.4, 0.4, 0.2]

    # rq-817541fc
    def test_scaling_a_charged_molecule_reports_the_extrapolation(self):
        charged = [q - 3.0 / len(self.CHARGES) for q in self.CHARGES]
        _, messages = _caught(scale_and_neutralize, charged, total_charge=-3.0)

        assert any("neutral" in m.lower() for m in messages), (
            f"scaling an anion said nothing about the factor being a "
            f"neutral-molecule parameterisation; messages were {messages}"
        )
        assert any("-3" in m for m in messages), (
            f"the warning does not name the formal charge: {messages}"
        )

    # rq-2d0042a6
    def test_scaling_a_neutral_molecule_reports_nothing(self):
        """The complement: the warning must not fire on the common case."""
        neutral = [q - sum(self.CHARGES) / len(self.CHARGES) for q in self.CHARGES]
        _, messages = _caught(scale_and_neutralize, neutral, total_charge=0.0)
        assert not messages, messages

    def test_the_correction_still_lands_on_the_formal_charge(self):
        """Pinned so a fix cannot answer the above by dropping the correction."""
        charged = [q - 3.0 / len(self.CHARGES) for q in self.CHARGES]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = scale_and_neutralize(charged, total_charge=-3.0)
        assert sum(out) == pytest.approx(-3.0, abs=1e-9)


# --------------------------------------------------------------------------- #
# Which model answered
# --------------------------------------------------------------------------- #
class TestTheRefinerSaysWhichModelAnswered:
    # rq-c202e1dc
    def test_the_refiner_reports_the_bundled_model(self):
        assert _DEFAULT_MODEL_PATH.exists(), "fixture is void: no bundled model"
        refiner = MLChargeRefinement()
        assert refiner.model_source == "bundled"

    # rq-3cf63133
    def test_a_substituted_fallback_announces_itself(self, tmp_path):
        refiner, messages = _caught(
            MLChargeRefinement, model_path=tmp_path / "absent.pkl"
        )
        assert refiner.model_source == "fallback"
        assert any("fallback" in m.lower() for m in messages), messages


# --------------------------------------------------------------------------- #
# Library version behind the pickle
# --------------------------------------------------------------------------- #
class TestTheLibraryIsCheckedAgainstTheArtifact:
    # rq-1e274fe6
    def test_a_model_from_a_different_sklearn_is_reported(self, monkeypatch):
        monkeypatch.setattr(
            MLChargeRefinement, "_recorded_sklearn_version",
            staticmethod(lambda path: "0.0.1-not-a-real-version"),
        )
        _, messages = _caught(MLChargeRefinement)

        ours = [m for m in messages if "charge" in m.lower()]
        assert ours, (
            f"the mismatch was not reported in terms of the charges it "
            f"affects; messages were {messages}"
        )
        assert any("0.0.1-not-a-real-version" in m for m in ours), ours
        assert any(sklearn.__version__ in m for m in ours), ours

    # rq-75d3624b
    def test_a_matching_version_is_not_reported(self, monkeypatch):
        monkeypatch.setattr(
            MLChargeRefinement, "_recorded_sklearn_version",
            staticmethod(lambda path: sklearn.__version__),
        )
        _, messages = _caught(MLChargeRefinement)
        assert not [m for m in messages if "scikit-learn" in m], messages


# --------------------------------------------------------------------------- #
# Geometry the charges belong to
# --------------------------------------------------------------------------- #
class TestChargesBelongToTheExportedGeometry:
    # rq-1e6d149c
    @pytest.mark.parametrize(
        "optimize,expect_single_point", [(False, True), (True, False)]
    )
    def test_the_default_is_a_single_point_calculation(
        self, optimize, expect_single_point
    ):
        """1SCF is MOPAC's keyword for one SCF pass and no optimisation.

        The default matters because it decides which geometry the charges
        describe: the one written to the .gro, or one AM1 re-optimised behind
        the caller's back.
        """
        from rdkit import Chem

        from biochar.charges.qm_charges import QMChargeAssigner

        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        Chem.SanitizeMol(mol)
        coords = [[0.0, 0.0, 0.0]] * mol.GetNumAtoms()

        assigner = QMChargeAssigner(optimize=optimize)
        deck = assigner._build_input(mol, coords, total_charge=0)

        assert ("1SCF" in deck) is expect_single_point, deck.splitlines()[0]
