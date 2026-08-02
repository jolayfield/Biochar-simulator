"""
Tests for the rqm/generation-config.md scenarios -- how GeneratorConfig resolves
a request, and which requests it refuses rather than quietly adjusts.

The distinction these turn on: a ValueError at construction means the request was
impossible, a ValidationError from generate() means the attempt fell short, and a
warning means the request was changed on the caller's behalf. Confusing the three
is how a request gets missed by a wide margin and reported as success.

Two scenarios assert behaviour the module does not have yet and carry
xfail(strict=True) naming the gap.
"""

import logging

import pytest

from biochar.constants import MIN_BUILDABLE_AROMATICITY
from biochar.pipeline.biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    ValidationError,
)


class TestCompositionResolution:
    # rq-90382d04
    def test_explicit_ratio_wins_over_the_derived_one(self):
        explicit = GeneratorConfig(
            target_num_carbons=20, temperature=600, feedstock="softwood",
            H_C_ratio=0.77,
        )
        derived = GeneratorConfig(
            target_num_carbons=20, temperature=600, feedstock="softwood",
        )

        assert explicit.H_C_ratio == 0.77
        assert derived.H_C_ratio != 0.77, (
            "fixture is void: the model happens to predict the explicit value"
        )
        # The targets that were not given are still derived, not defaulted.
        assert explicit.O_C_ratio == derived.O_C_ratio
        assert explicit.O_C_ratio != 0.1, "O/C fell back to the default, not the model"

    # rq-c58a2ebe
    def test_no_temperature_falls_back_to_the_fixed_defaults(self):
        cfg = GeneratorConfig(target_num_carbons=20)

        assert cfg.H_C_ratio == 0.5
        assert cfg.O_C_ratio == 0.1
        assert cfg.aromaticity_percent == 90.0


class TestAromaticFloor:
    # rq-1762f896
    def test_derived_aromaticity_below_the_floor_is_clamped(self, caplog):
        # A low pyrolysis temperature predicts an aromaticity the sheet builder
        # cannot realise.
        with caplog.at_level(logging.WARNING, logger="biochar"):
            cfg = GeneratorConfig(
                target_num_carbons=20, temperature=250, feedstock="manure",
            )

        assert cfg.aromaticity_percent == MIN_BUILDABLE_AROMATICITY
        assert any("clamp" in r.message.lower() for r in caplog.records), (
            f"clamping happened without a warning; records: "
            f"{[r.message[:60] for r in caplog.records]}"
        )


class TestImpossibleRequestsRefused:
    # rq-9aa547f6
    def test_unattainable_hc_is_refused_at_construction(self):
        with pytest.raises(ValueError, match=r"H_C_ratio"):
            GeneratorConfig(target_num_carbons=20, H_C_ratio=2.0)

        # Refused at construction, not deferred to generate().
        try:
            GeneratorConfig(target_num_carbons=20, H_C_ratio=2.5)
        except ValueError as exc:
            assert "2.0" in str(exc), "the message should state the ceiling"
        else:
            pytest.fail("H_C_ratio=2.5 was accepted")

    # rq-a9c68e17
    @pytest.mark.parametrize("bad_ph", [-0.5, 14.5])
    def test_ph_outside_the_meaningful_range_is_refused(self, bad_ph):
        with pytest.raises(ValueError, match=r"pH"):
            GeneratorConfig(target_num_carbons=20, pH=bad_ph)


class TestChargeBackendCompatibility:
    # rq-9bace66f
    def test_ph_with_the_ml_backend_is_refused(self):
        with pytest.raises(ValueError) as exc:
            GeneratorConfig(target_num_carbons=20, pH=5.0, charge_method="ml")

        message = str(exc.value)
        assert "ml" in message
        # The caller is told what to use instead, not just what not to use.
        assert "opls" in message or "qm" in message, (
            f"refusal names no alternative backend: {message}"
        )

    def test_the_other_backends_accept_a_ph(self):
        for method in ("opls", "qm"):
            cfg = GeneratorConfig(target_num_carbons=20, pH=5.0, charge_method=method)
            assert cfg.pH == 5.0


class TestStrictnessCoversShortfall:
    """A request for N groups answered with far fewer is not a success.

    The fixture asks for 40 ether bridges on a 60-carbon skeleton. Ether consumes
    two aromatic sites per bridge, so the request cannot be met; the tolerances
    are opened so that composition validation does not fire first and mask which
    guard is under test.
    """

    @staticmethod
    def _shortfall_config(**extra):
        kwargs = dict(
            target_num_carbons=60, seed=1, functional_groups={"ether": 40},
            H_C_tolerance=1.0, O_C_tolerance=1.0,
        )
        kwargs.update(extra)
        return GeneratorConfig(**kwargs)

    def test_the_fixture_really_does_under_place(self):
        """Guard the guard: if placement ever succeeds, the scenario below is void."""
        _, _, comp = BiocharGenerator(self._shortfall_config(strict=False)).generate()

        requested = comp.requested_counts.get("ether", 0)
        placed = comp.placed_counts.get("ether", 0)
        assert requested == 40
        assert 0 < placed < requested, (
            f"fixture no longer under-places: requested {requested}, placed {placed}"
        )

    # rq-4fc714fe
    def test_partial_placement_does_not_satisfy_strict_validation(self):
        with pytest.raises(ValidationError) as exc:
            BiocharGenerator(self._shortfall_config(strict=True)).generate()

        message = str(exc.value)
        assert "ether" in message
        assert "40" in message, (
            f"the failure should state what was requested: {message}"
        )

    # rq-ab986122
    def test_a_ratio_derived_shortfall_is_judged_by_the_ratio(self):
        """O/C-driven placement populates requested_counts too.

        The shortfall check must not fire there: the caller asked for a ratio,
        and the count is an implementation detail of reaching it, judged by
        O_C_tolerance instead. Without this scoping the fix above would make
        every ratio-driven strict run fail on a one-group miss.
        """
        cfg = GeneratorConfig(
            target_num_carbons=40, seed=1, strict=True, O_C_ratio=0.15,
            functional_groups=None,
        )
        _, _, comp = BiocharGenerator(cfg).generate()

        assert comp.requested_counts, (
            "fixture is void: ratio-driven placement recorded no requested counts"
        )


class TestExtrapolationIsReported:
    """Scoped to the properties this config derives -- H/C and O/C.

    An earlier version of this scenario asked for a temperature inside H/C's
    range but outside "another derived property's". No such temperature exists:
    H/C and O/C are both fitted over the whole 100-1000 C grid. The properties
    with narrower support -- pH, conductivity -- are ones the generator never
    derives, and warning about them here would attach noise to a structure
    request that never consulted them.
    """

    @staticmethod
    def _derived_range(prop="H_C_ratio"):
        from biochar.models.temperature_model import get_default_model

        return get_default_model().get_valid_range(prop)

    # rq-a1da340c
    def test_a_temperature_beyond_a_derived_property_is_reported(self, caplog):
        _, t_max = self._derived_range()
        beyond = t_max + 500

        with caplog.at_level(logging.WARNING, logger="biochar"):
            GeneratorConfig(target_num_carbons=20, temperature=beyond)

        messages = [r.getMessage() for r in caplog.records]
        assert any("outside the data range" in m for m in messages), (
            f"T={beyond} is beyond H/C's support of {self._derived_range()} and "
            f"nothing warned; records: {[m[:60] for m in messages]}"
        )
        assert any("H_C_ratio" in m for m in messages), (
            f"the warning does not name the property: {messages}"
        )

    # rq-0d326d5b
    def test_a_property_this_config_never_derives_is_not_judged(self, caplog):
        from biochar.models.temperature_model import get_default_model

        model = get_default_model()
        derived_max = min(
            model.get_valid_range(p)[1] for p in ("H_C_ratio", "O_C_ratio")
        )
        narrower = [
            name for name, spec in model._m["properties"].items()
            if name not in ("H_C_ratio", "O_C_ratio")
            and spec.get("t_max") is not None
            and spec["t_max"] < derived_max
        ]
        assert narrower, "no property has narrower support than the derived ones"

        # In range for everything derived; out of range for pH and conductivity.
        with caplog.at_level(logging.WARNING, logger="biochar"):
            GeneratorConfig(target_num_carbons=20, temperature=derived_max)

        assert not [
            r for r in caplog.records if "outside the data range" in r.getMessage()
        ], (
            f"warned at {derived_max} C about a property this config never derives; "
            f"narrower properties are {sorted(narrower)}"
        )


class TestQuestionableCombinationsWarn:
    # rq-8155aecf
    def test_crossed_nitrogen_doping_warns_but_builds(self, caplog):
        with caplog.at_level(logging.WARNING, logger="biochar"):
            cfg = GeneratorConfig(
                target_num_carbons=40, seed=1, num_pyridinic=1, strict=False,
                functional_groups={"carboxyl": 1},
            )
            # The warning describes the structure being built, so it is emitted
            # by generate(), not by construction.
            mol, _, _ = BiocharGenerator(cfg).generate()

        assert mol is not None, "the structure should still be produced" 

        assert any(
            "doping" in r.message.lower() for r in caplog.records
        ), (
            f"crossed doping did not warn; records: "
            f"{[r.message[:60] for r in caplog.records]}"
        )
