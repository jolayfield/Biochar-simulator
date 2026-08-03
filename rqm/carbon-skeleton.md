# Feature: Carbon Skeleton Growth <!-- rq-f2d55972 -->

`PAHAssembler` (`src/biochar/pipeline/carbon_skeleton.py`) builds the polycyclic aromatic
skeleton every structure starts from. It is stage one of the pipeline: everything downstream —
heteroatoms, geometry, typing, export — decorates what this produces.

A skeleton is assembled from **whole rings**. That single fact shapes every requirement here.
Carbon counts are reachable only in the increments ring fusion allows, a "10% pentagon" request is
a probability applied per ring rather than a share of the finished structure, and aromaticity is
not an input at all — it is whatever the ring topology turns out to be.

So most of what a caller passes is an approximation, and the module's job is to be honest about
which. The failure mode is a request that is quietly not honoured: a parameter accepted and
ignored, a count met by a different number, or a fallback that returns something unrelated to
what was asked for.

## Feature API <!-- rq-f61f4d50 -->

- `PAHAssembler(seed: int = None)` <!-- rq-5840d1b3 -->
  - `generate(target_num_carbons: int, prefer_larger_pahs: bool = True, defect_fraction: float = 0.0, target_h_c: Optional[float] = None, heptagon_fraction: float = 0.0) -> CarbonSkeleton`
    builds a skeleton of approximately `target_num_carbons` carbons.
  - Every parameter it accepts changes the result. A knob that cannot affect the output is not
    part of the signature.
  - Targets within the library's range are satisfied exactly from a pre-validated structure;
    larger targets take the closest library seed and grow by ring fusion.

- `CarbonSkeleton` <!-- rq-419ffa85 -->
  - The assembled `mol` plus its ring composition. The realised carbon count is read from the
    molecule, not assumed from the request.

- `SkeletonValidator.validate(skeleton) -> Tuple[bool, List[str]]` <!-- rq-83521d1c -->
  - Structural checks on an assembled skeleton.

## A Target Within the Library Is Met Exactly <!-- rq-f7dda304 -->

`PAH_LIBRARY` holds pre-validated structures from benzene at 6 carbons up to 40. A target matching
one of those sizes is satisfied from the library rather than grown, so it is met exactly and the
result is a structure whose aromaticity and ring closure are already known good.

```gherkin
Feature: Use the validated library where it covers the request

  @rq-4af71330
  Scenario: A target matching a library size is met exactly
    Given a target carbon count that a library structure provides
    When a skeleton is generated
    Then the skeleton has exactly that many carbons
```

## A Larger Target Is Reached From Above, Never From Below <!-- rq-3f012547 -->

Above the library, growth proceeds by fusing whole rings, so only certain carbon counts are
reachable. When the requested count is not one of them, growth continues to the next reachable
count rather than stopping short: **the realised count is greater than or equal to the request,
never less.**

The direction is the requirement. A skeleton one ring short of the request has fewer edge sites
than the caller planned for, which propagates into the achievable H/C and O/C; overshooting by a
ring leaves the composition targets reachable. Deviations of up to about five carbons occur at
awkward targets.

The realised count is a property of the returned molecule, so a caller who needs the exact figure
measures it rather than assuming the request was met.

```gherkin
Feature: Reach an unachievable count from above

  @rq-59e23c2e
  Scenario: A target above the library is met or exceeded
    Given a target carbon count larger than any library structure
    When a skeleton is generated
    Then the skeleton has at least that many carbons

  @rq-5435c9d8
  Scenario: The overshoot stays within a ring or two
    Given a range of awkward target carbon counts
    When skeletons are generated for each
    Then none exceeds its target by more than eight carbons
```

## Defect Fractions Are Per-Ring Probabilities <!-- rq-9c8fe79c -->

`defect_fraction` and `heptagon_fraction` are the probability that **any one ring added during
growth** is a pentagon or a heptagon. They are not the fraction of rings in the finished
structure.

The realised share lands below the requested probability, because the seed a skeleton grows from
is all hexagons and only the added rings are subject to the draw. Asking for 0.3 yields roughly
0.22 of the finished rings.

Naming them as probabilities rather than fractions matters because the parameter is how the Wood
et al. curvature model is reached, and a caller reproducing published ring statistics is comparing
against realised fractions.

```gherkin
Feature: Apply defect probabilities per ring addition

  @rq-c6b675bb
  Scenario: A zero defect fraction produces a pure hexagonal skeleton
    Given a defect fraction and a heptagon fraction of zero
    When a skeleton is generated
    Then every ring is a hexagon

  @rq-0c91f5d5
  Scenario: A positive defect fraction produces pentagons
    Given a positive defect fraction
    When a large skeleton is generated
    Then the skeleton contains at least one pentagon
    And the realised pentagon share does not exceed the requested probability
```

## Aromaticity Is an Outcome, Not a Request <!-- rq-cf3d7a56 -->

Aromaticity is decided by the ring topology the assembler builds. There is no way to ask for a
skeleton of a given aromaticity and no mechanism that could honour such a request, so the
assembler does not accept one.

A parameter that is accepted and documented as unused is worse than no parameter: it reads at the
call site as a request being made, and a caller who sets it has no way to discover that nothing
changed. The composition-level aromaticity target is handled in `generation-config.md`, which
clamps it to what the builder can realise — that is where the concept belongs.

```gherkin
Feature: Do not accept a request that cannot be honoured

  @rq-b6eec273
  Scenario: The assembler takes no aromaticity argument
    Given the skeleton generation entry point
    When its accepted parameters are inspected
    Then none of them is an aromaticity target
```

## A Skeleton Is Never Returned at an Unrelated Size <!-- rq-fb1a5f14 -->

If growth cannot produce a skeleton, that is a failure of the request and is raised. Returning a
fixed small structure instead — the library's pyrene, at 16 carbons — answers a request for 200
carbons with something an eighth of the size, and every downstream stage then works on it as
though it were what was asked for.

A caller can recover from an exception. A caller cannot recover from a plausible-looking molecule
of the wrong size, because nothing later in the pipeline re-checks the request.

```gherkin
Feature: Fail rather than substitute an unrelated skeleton

  @rq-f953b9b0
  Scenario: A skeleton that cannot be grown raises
    Given a request that graph growth cannot satisfy
    When a skeleton is generated
    Then the failure is raised rather than answered with a fixed structure
```

## A Library Entry That Cannot Be Loaded Is Reported <!-- rq-e7260a44 -->

The library is parsed once at construction. An entry whose SMILES will not parse, or will not
sanitise, is dropped from the working set.

Dropping it silently shrinks the library without changing anything a caller can see: targets that
were met exactly begin to be grown instead, and the structure that was pre-validated is quietly
replaced by one that was not. `PAH_LIBRARY` states that all its entries are validated; a dropped
entry means that statement is no longer true, and the discrepancy is worth surfacing rather than
absorbing.

```gherkin
Feature: Surface a library entry that does not load

  @rq-4dc5ce35
  Scenario: An unparseable library entry is reported
    Given a library containing an entry whose SMILES cannot be parsed
    When the assembler is constructed
    Then the entry is reported rather than silently dropped
```

## Cross-references <!-- rq-a60c7e3b -->

- The composition targets that drive `target_num_carbons` and `target_h_c`, and the aromaticity
  floor applied above this layer, are specified in `generation-config.md`.
- The hex-lattice coordinates this stage attaches, and why large sheets skip clash resolution,
  are specified in `geometry-embedding.md`.
- Functional-group placement onto the edge sites this stage produces is specified in
  `heteroatom-assignment.md`.
