# Feature: Temperature × Feedstock Composition Model <!-- rq-7e4059d9 -->

`src/biochar/models/temperature_model.py` turns a pyrolysis temperature, and optionally a
feedstock, into the composition targets the generator builds from. It is fitted offline from the
UC Davis biochar database and shipped as a JSON artifact; at run time it interpolates that
artifact and never touches the source data.

Every number it returns is an estimate from a finite, uneven dataset. Some grid points rest on
hundreds of observations, some on none. Some feedstocks have enough measurements for their own
curve over part of the temperature range and fall back to the pooled curve elsewhere. A
prediction at 1000 °C for a property fitted only to 900 °C is an extrapolation.

**The model already records what each prediction rests on** — the observation count, the spread,
the temperature support per property, the H/C range the aromaticity fit was taken over. None of
that reaches the caller. A number derived from zero observations is returned in the same shape,
and with the same apparent authority, as one derived from two hundred. The requirements here are
mostly about closing that gap: the estimate and its evidence travel together, or the caller
cannot tell a measurement from a guess.

## Feature API <!-- rq-766e9c83 -->

- `TemperatureModel(model_path: Optional[Path] = None)` <!-- rq-7ad3970f -->
  - Loads the fitted JSON artifact. `get_default_model()` returns a shared instance.
  - `predict(temperature: float, prop: str, feedstock: Optional[str] = None) -> float` —
    interpolates the fitted curve, clamping at the ends of the temperature grid. Raises `KeyError`
    for an unknown property.
  - `composition(temperature, feedstock=None) -> Dict[str, float]` — the three generator targets:
    `H_C_ratio`, `O_C_ratio`, and `aromaticity_percent` derived from the predicted H/C.
  - `predict_all(temperature, feedstock=None) -> Dict[str, float]` — every modelled property.
  - `aromaticity_from_hc(hc: float) -> float` — the linear aromaticity fit, clamped to 0–100.
  - `get_valid_range(prop: str, feedstock: Optional[str] = None) -> Optional[Tuple[float, float]]`
    — the temperature support of the curve that a `predict` call for `prop` would use.
  - `feedstocks -> Tuple[str, ...]` and `provenance -> dict` — what the artifact was built from.

- `properties(temperature: float, feedstock: Optional[str] = None) -> Dict[str, float]` <!-- rq-f1f301e0 -->
  - Module-level convenience over the default model, returning the full reference property table.

- `VALID_FEEDSTOCKS: Tuple[str, ...]` <!-- rq-5f4e2e68 -->
  - The feedstock names accepted anywhere in the package. Any other value is an error at the
    caller, not a silent fallback here.

- `build_model(output_path: Optional[Path] = None) -> Path` <!-- rq-c145dabd -->
  - Offline refit from the source CSVs. Not called at run time.

## Feedstock-Specific Curves Are Narrow and Gated <!-- rq-f01d97ff -->

A per-feedstock curve is emitted only for **`H_C_ratio` and `O_C_ratio`**, and only when that
feedstock's data passes a sufficiency gate: at least `_GATE_MIN_N` = 30 observations spanning at
least `_GATE_MIN_TSPAN` = 300 °C. Every other property — pH, ash, surface area, conductivity — is
pooled across feedstocks, so `properties(600, "manure")` and `properties(600, "softwood")` differ
in two entries and agree in the rest.

This is a deliberate limit rather than an oversight: a curve fitted to a handful of points
clustered in one temperature decade is worse than the pooled curve, not better.

```gherkin
Feature: Use a feedstock-specific curve only where the data supports one

  @rq-39bd2478
  Scenario: Only H/C and O/C vary by feedstock
    Given two feedstocks that both have their own curves
    When the full property table is requested at one temperature for each
    Then H/C and O/C differ between them
    And every other property is identical

  @rq-3b2ef55d
  Scenario: A feedstock without sufficient data has no curve of its own
    Given a feedstock whose observations fail the sufficiency gate
    When its H/C is requested
    Then the pooled curve answers
```

## A Prediction Says Which Curve Answered It <!-- rq-c5d8e160 -->

A feedstock override covers only the temperature range its own data spans. Outside that range the
pooled curve answers instead, which is the right numerical choice and an invisible one: the return
value is a bare float, identical in type and shape to an override-backed answer.

`corn_stover` at 900 °C returns exactly the pooled value, and nothing in the call says so.

The caller is told which curve was used, so that a sweep varying feedstock can distinguish points
where the feedstock mattered from points where it silently did not.

```gherkin
Feature: Make the fallback from a feedstock curve to the pooled curve visible

  @rq-b6522656
  Scenario: A temperature outside a feedstock's support falls back visibly
    Given a feedstock whose curve does not span the requested temperature
    When its H/C is predicted there
    Then the pooled value is returned
    And the fallback is reported rather than silent
```

## The Support Range Is Per Property <!-- rq-c2dd4be4 -->

Properties are fitted over different temperature ranges: H/C spans 100–1000 °C, pH 200–900 °C,
electrical conductivity 220–900 °C. A temperature can therefore be comfortably inside the support
of one property and well outside another's.

`get_valid_range` answers for the property asked about. A range that ignores its argument
describes H/C and mislabels everything else, so a caller checking before predicting is told the
request is in range when it is not — which is exactly how a 1000 °C request reaches pH and
conductivity with no warning at all.

```gherkin
Feature: Report the support of the property being asked about

  @rq-876c3bc5
  Scenario: The reported range is the one belonging to the property
    Given two properties fitted over different temperature ranges
    When the valid range is requested for each
    Then each answer matches that property's own support

  @rq-a234cdbe
  Scenario: A temperature inside one property's support but outside another's is distinguishable
    Given a temperature within the H/C range but beyond the pH range
    When the valid range is consulted for each property
    Then the H/C answer contains the temperature and the pH answer does not
```

## An Estimate Carries the Evidence Behind It <!-- rq-5d9abc46 -->

The artifact records, for every property at every grid point, the number of observations `n` and
their `spread`. Those are the difference between a fitted value and a plausible-looking one.

Grid points with no observations are filled from the nearest finite neighbour so the curve stays
continuous — a reasonable numerical choice that becomes a reporting problem when the filled value
is returned indistinguishably from a fitted one. Electrical conductivity has no observations at
all at 100 °C and predicts 1.0 there, carried in from a neighbouring point.

A caller can therefore ask what a prediction rests on, and a prediction resting on nothing says so.

```gherkin
Feature: Let a caller tell a fitted value from a filled one

  @rq-c0b934a0
  Scenario: A prediction reports the observations behind it
    Given a property and a temperature
    When the prediction is requested with its evidence
    Then the observation count and spread at that temperature are available

  @rq-a687f7a1
  Scenario: A grid point with no observations is reported as absent, not fitted
    Given a property with zero observations at some grid temperature
    When it is predicted there
    Then the result is marked as carried in from a neighbouring point
```

## Aromaticity Is Extrapolated Beyond Its Fit <!-- rq-7ef42231 -->

`aromaticity_percent` is not measured against temperature. It comes from a linear fit of
aromaticity against H/C, taken over H/C values between **0.03 and 1.64**, then clamped into
0–100%.

The clamp keeps the output physical and hides how far outside the fit an input was. An H/C of 5.0
is triple the largest value ever fitted, and the clamp turns the resulting negative number into
exactly 0.0% — a value that reads as a confident prediction of a fully aromatic char rather than
as a refusal.

Since `generation-config.md` allows any H/C below 2.0, inputs beyond the fit are reachable through
ordinary configuration rather than only through misuse.

```gherkin
Feature: Do not present an extrapolated aromaticity as a fitted one

  @rq-27f04167
  Scenario: An H/C beyond the fit range is reported rather than silently clamped
    Given an H/C larger than the maximum the aromaticity fit was taken over
    When the aromaticity is derived from it
    Then the caller can tell the value was extrapolated

  @rq-2bf0ee76
  Scenario: An H/C inside the fit range is answered normally
    Given an H/C within the fitted range
    When the aromaticity is derived from it
    Then the value is returned with no extrapolation reported
```

## Provenance Is Recorded <!-- rq-304bf22a -->

The artifact names the dataset it was built from, so a structure generated today can be traced to
the data that shaped its targets. The source is the UC Davis biochar database; the model records
its URL alongside the fit statistics.

One column choice is a judgement rather than a reading: total carbon is taken from the database's
"Total C (%) Used in Plot" column rather than the raw elemental column, and that decision is
recorded where the column indices are defined.

```gherkin
Feature: Record what the model was built from

  @rq-f08462b9
  Scenario: The model reports its source dataset
    Given a loaded model
    When its provenance is read
    Then it names the source database
    And it reports the fit statistics for the aromaticity relation
```

## Cross-references <!-- rq-6deffbcf -->

- The generator's use of these predictions — which targets are derived, which are explicit, and
  the buildable-aromaticity floor applied afterwards — is specified in `generation-config.md`.
  The out-of-range warning it emits depends on the per-property support described above.
- The dataset and its provenance are described in `docs/` and in the module docstring; this
  document specifies the model's behaviour, not the data's content.
