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

import json
import warnings

import numpy as np
import pytest

from biochar.charges.ml_charges import (
    MLChargeRefinement,
    _DEFAULT_MODEL_PATH,
    _generate_training_data,
)
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

    # rq-712dd0d8
    def test_a_model_that_does_not_rebuild_to_its_charges_is_reported(
        self, tmp_path
    ):
        """The version string is provenance; this is the question it stood in for.

        A recorded version tells you which library fitted the numbers. It does
        not tell you whether the library in hand computes the same charges from
        them -- and that is the only thing a reader of those charges needs.
        """
        payload = json.loads(_DEFAULT_MODEL_PATH.read_text())
        payload["reference_charges"] = [
            q + 0.01 for q in payload["reference_charges"]
        ]
        tampered = tmp_path / "drifted.json"
        tampered.write_text(json.dumps(payload))

        _, messages = _caught(MLChargeRefinement, model_path=tampered)

        ours = [m for m in messages if "invalid" in m.lower()]
        assert ours, (
            f"a model that does not reproduce its own recorded charges was "
            f"loaded without saying so; messages were {messages}"
        )
        assert any("1e-06" in m for m in ours), (
            f"the warning does not name the tolerance it failed: {ours}"
        )
        assert any("1.00e-02" in m for m in ours), (
            f"the warning does not name the deviation it measured: {ours}"
        )

    # rq-fa93a423
    def test_a_model_that_rebuilds_exactly_is_not_reported(self):
        """The complement: silence has to be earned on the shipped artifact."""
        _, messages = _caught(MLChargeRefinement)
        assert not [m for m in messages if "invalid" in m.lower()], messages


# --------------------------------------------------------------------------- #
# The artifact against the reference data it claims
# --------------------------------------------------------------------------- #
class TestTheArtifactIsCheckedAgainstItsReferenceData:
    # rq-da875aee
    def test_the_bundled_training_data_is_the_current_reference_data(self):
        """The artifact is pinned data; the data it was pinned from is not.

        This has already happened once. The typer stopped calling benzoic
        acid's carboxyl a ketone plus a phenol, and the bundled model went on
        predicting from targets that said it was -- through every release after
        the fix, because a predicted charge carries no provenance and nothing
        compared the two.

        A failure here is not necessarily a bug. It is a decision to make in a
        commit: rebuild the artifact and take the force-field change, or record
        why it stays behind.
        """
        payload = json.loads(_DEFAULT_MODEL_PATH.read_text())
        X_now, y_now = _generate_training_data()

        stored_X = np.asarray(payload["X"], dtype=float)
        stored_y = np.asarray(payload["y"], dtype=float)

        assert stored_X.shape == X_now.shape, (
            f"the bundled artifact was fitted to {stored_X.shape[0]} reference "
            f"atoms; the package now generates {X_now.shape[0]}"
        )
        assert np.array_equal(stored_X, X_now), (
            "the bundled artifact's training features are not the ones the "
            "package generates now -- the featuriser or the reference "
            "structures moved without the artifact"
        )

        # Exact equality is the wrong test. The equilibration step sums the raw
        # table charges, and the summation order differs enough between RDKit
        # releases to move a target by one ULP -- 1e-16 e, on a third of the
        # reference atoms, meaning nothing. The drift worth catching is the
        # kind that has already happened: a retyped functional group, tenths of
        # an electron. Anything above 1e-12 e is a change somebody made.
        drifted = np.flatnonzero(np.abs(stored_y - y_now) > 1e-12)
        assert drifted.size == 0, (
            f"the bundled artifact's target charges disagree with the OPLS "
            f"table at {drifted.size} of {y_now.size} reference atoms "
            f"(max {np.abs(stored_y - y_now).max():.4f} e). Rebuild it with "
            f"biochar.charges.ml_charges.build_and_save_bundled_model() and "
            f"record the force-field change, or say here why it stays behind."
        )


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
