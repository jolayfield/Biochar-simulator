"""
Tests for the rqm/temperature-model.md scenarios -- what a prediction rests on,
and whether the caller can tell.

Distinct from tests/test_temperature_model.py, which checks that predicted
values land in plausible ranges. These check the reporting contract: the model
records an observation count, a spread, a per-property temperature support, and
the H/C range the aromaticity fit was taken over, and none of that currently
reaches a caller. A value carried in from a neighbouring grid point is returned
in the same shape as one fitted from two hundred observations.

Fixtures read the ranges and counts out of the artifact rather than hard-coding
them, so a refit moves the tests with the data instead of breaking them.
"""

import pytest

from biochar.models.temperature_model import (
    VALID_FEEDSTOCKS,
    TemperatureModel,
    get_default_model,
    properties,
)


@pytest.fixture(scope="module")
def model():
    return get_default_model()


def _override(model, feedstock, prop="H_C_ratio"):
    return model._m["feedstock_overrides"].get(feedstock, {}).get(prop)


def _feedstock_with_override(model, prop="H_C_ratio"):
    for fs in sorted(model._m["feedstock_overrides"]):
        if _override(model, fs, prop):
            return fs
    pytest.skip("no feedstock has its own curve in this artifact")


class TestFeedstockCurvesAreNarrow:
    # rq-39bd2478
    def test_only_hc_and_oc_vary_by_feedstock(self, model):
        with_curves = [
            fs for fs in VALID_FEEDSTOCKS if _override(model, fs) is not None
        ]
        assert len(with_curves) >= 2, "need two feedstocks with their own curves"
        a, b = with_curves[0], with_curves[1]

        table_a = properties(600, a)
        table_b = properties(600, b)
        assert set(table_a) == set(table_b)

        differing = {k for k in table_a if table_a[k] != table_b[k]}
        # aromaticity is derived from H/C, so it moves with it.
        assert differing <= {"H_C_ratio", "O_C_ratio", "aromaticity_percent"}, (
            f"properties other than H/C and O/C differ by feedstock: "
            f"{differing - {'H_C_ratio', 'O_C_ratio', 'aromaticity_percent'}}"
        )
        assert {"H_C_ratio", "O_C_ratio"} <= differing, (
            "H/C and O/C are identical across two feedstocks that each have a curve"
        )

    # rq-3b2ef55d
    def test_a_feedstock_without_a_curve_uses_the_pooled_one(self, model):
        without = [fs for fs in VALID_FEEDSTOCKS if _override(model, fs) is None]
        if not without:
            pytest.skip("every feedstock passed the sufficiency gate in this artifact")

        fs = without[0]
        assert model.predict(600, "H_C_ratio", fs) == model.predict(600, "H_C_ratio")


class TestFallbackIsVisible:
    # rq-b6522656
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "predict() drops to the pooled curve outside a feedstock's support and "
            "returns a bare float, identical to an override-backed answer. Retire "
            "this marker with the fix."
        ),
    )
    def test_a_fallback_to_the_pooled_curve_is_reported(self, model, caplog):
        import logging

        fs = _feedstock_with_override(model)
        support = _override(model, fs)
        beyond = float(support["t_max"]) + 50

        with caplog.at_level(logging.INFO, logger="biochar"):
            value = model.predict(beyond, "H_C_ratio", fs)

        # The precondition: it really did fall back.
        assert value == model.predict(beyond, "H_C_ratio"), (
            f"{fs} still has its own curve at {beyond} C; pick a further temperature"
        )
        assert caplog.records, (
            f"{fs} at {beyond} C silently returned the pooled value; nothing in the "
            f"call says the feedstock was ignored"
        )


class TestSupportRangeIsPerProperty:
    @staticmethod
    def _two_properties_with_different_ranges(model):
        props = model._m["properties"]
        ranges = {
            name: (spec["t_min"], spec["t_max"])
            for name, spec in props.items()
            if spec.get("t_min") is not None and spec.get("t_max") is not None
        }
        wide = max(ranges, key=lambda k: ranges[k][1])
        narrow = min(ranges, key=lambda k: ranges[k][1])
        if ranges[wide] == ranges[narrow]:
            pytest.skip("every property shares one temperature range in this artifact")
        return (wide, ranges[wide]), (narrow, ranges[narrow])

    # rq-876c3bc5
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "get_valid_range takes no property argument and returns H/C's range "
            "for everything. Retire this marker with the fix."
        ),
    )
    def test_the_reported_range_belongs_to_the_property(self, model):
        (wide, wide_range), (narrow, narrow_range) = (
            self._two_properties_with_different_ranges(model)
        )

        assert model.get_valid_range(wide) == wide_range
        assert model.get_valid_range(narrow) == narrow_range

    # rq-a234cdbe
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "Same root cause: one range is reported for every property, so a "
            "temperature cannot be in range for one and out for another. Retire "
            "this marker with the fix."
        ),
    )
    def test_a_temperature_can_be_in_range_for_one_property_and_not_another(self, model):
        (wide, wide_range), (narrow, narrow_range) = (
            self._two_properties_with_different_ranges(model)
        )
        between = narrow_range[1] + 1.0
        assert between <= wide_range[1], "fixture is void: no gap between the ranges"

        assert wide_range[0] <= between <= wide_range[1]
        reported_narrow = model.get_valid_range(narrow)
        assert not (reported_narrow[0] <= between <= reported_narrow[1]), (
            f"{between} C is beyond {narrow}'s support of {narrow_range} but the "
            f"reported range {reported_narrow} contains it"
        )


class TestEvidenceTravelsWithTheEstimate:
    @staticmethod
    def _property_with_an_empty_grid_point(model):
        for name, spec in model._m["properties"].items():
            counts = spec.get("n") or []
            for i, count in enumerate(counts):
                if count == 0:
                    return name, float(model._grid[i])
        pytest.skip("no property has an empty grid point in this artifact")

    # rq-c0b934a0
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "n and spread are recorded per grid point but no query surfaces them. "
            "Retire this marker with the fix."
        ),
    )
    def test_a_prediction_reports_the_observations_behind_it(self, model):
        evidence = model.predict_with_evidence(600, "H_C_ratio")

        assert "value" in evidence
        assert evidence["n"] > 0
        assert "spread" in evidence

    # rq-a687f7a1
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "A grid point with zero observations is filled from its nearest finite "
            "neighbour and returned indistinguishably from a fitted value. Retire "
            "this marker with the fix."
        ),
    )
    def test_an_empty_grid_point_is_reported_as_carried_in(self, model):
        prop, temperature = self._property_with_an_empty_grid_point(model)

        evidence = model.predict_with_evidence(temperature, prop)
        assert evidence["n"] == 0
        assert evidence.get("filled") is True, (
            f"{prop} at {temperature} C rests on no observations but is not "
            f"marked as carried in from a neighbour"
        )


class TestAromaticityExtrapolation:
    # rq-27f04167
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "aromaticity_from_hc clamps into 0-100 and never consults the fit's "
            "hc_min/hc_max, so an H/C far outside the fit returns exactly 0.0 as "
            "though it were a confident prediction. Retire this marker with the fix."
        ),
    )
    def test_an_hc_beyond_the_fit_is_reported_as_extrapolated(self, model):
        fit = model._m["aromaticity_fit"]
        beyond = float(fit["hc_max"]) * 3

        with pytest.warns(UserWarning, match="extrapolat"):
            model.aromaticity_from_hc(beyond)

    # rq-2bf0ee76
    def test_an_hc_inside_the_fit_is_answered_normally(self, model, recwarn):
        fit = model._m["aromaticity_fit"]
        inside = (float(fit["hc_min"]) + float(fit["hc_max"])) / 2

        value = model.aromaticity_from_hc(inside)
        assert 0.0 <= value <= 100.0
        assert not [w for w in recwarn if "extrapolat" in str(w.message).lower()], (
            "an in-range H/C reported extrapolation"
        )


class TestProvenance:
    # rq-f08462b9
    def test_the_model_reports_its_source_and_fit_statistics(self, model):
        prov = model.provenance
        assert prov.get("primary_url"), "provenance names no source dataset"
        assert "davis" in prov["primary_url"].lower(), (
            f"unexpected source: {prov['primary_url']}"
        )

        fit = model._m["aromaticity_fit"]
        for key in ("a", "b", "n", "r2", "hc_min", "hc_max"):
            assert key in fit, f"aromaticity fit does not report {key}"
        assert fit["n"] > 0 and 0.0 <= fit["r2"] <= 1.0


def test_the_default_model_is_shared():
    """get_default_model caches; a fresh TemperatureModel does not have to."""
    assert get_default_model() is get_default_model()
    assert isinstance(TemperatureModel(), TemperatureModel)
