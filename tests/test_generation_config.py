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
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Strict mode raises only when a group places zero "
            "(biochar_generator.py, `zero_placed`), so 6 of 40 counts as success. "
            "Retire this marker with the fix."
        ),
    )
    def test_partial_placement_does_not_satisfy_strict_validation(self):
        with pytest.raises(ValidationError) as exc:
            BiocharGenerator(self._shortfall_config(strict=True)).generate()

        message = str(exc.value)
        assert "ether" in message
        assert "40" in message, (
            f"the failure should state what was requested: {message}"
        )


class TestExtrapolationIsReported:
    # rq-a1da340c
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "TemperatureModel.get_valid_range takes no property argument and "
            "always returns H/C's range, so a temperature outside a narrower "
            "property's support warns about nothing. Retire this marker with the fix."
        ),
    )
    def test_temperature_beyond_a_narrower_property_is_reported(self, caplog):
        from biochar.models.temperature_model import get_default_model

        model = get_default_model()
        props = model._m["properties"]
        hc = props["H_C_ratio"]
        narrower = {
            name: spec for name, spec in props.items()
            if spec.get("t_max") is not None and spec["t_max"] < hc["t_max"]
        }
        assert narrower, "no property has a narrower range than H/C; scenario is void"

        # Inside H/C's range, outside the others'.
        temperature = hc["t_max"]

        with caplog.at_level(logging.WARNING, logger="biochar"):
            GeneratorConfig(target_num_carbons=20, temperature=temperature)

        assert any("outside the data range" in r.message for r in caplog.records), (
            f"T={temperature} is beyond the support of "
            f"{sorted(narrower)[:3]} and no warning was emitted"
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
