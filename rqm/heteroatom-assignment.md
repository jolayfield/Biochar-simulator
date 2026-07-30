# Feature: Heteroatom and Functional Group Assignment <!-- rq-ffb45a6e -->

`OxygenAssigner` and `HydrogenAssigner` (`src/biochar/pipeline/heteroatom_assignment.py`) decorate a bare PAH
skeleton with functional groups and then saturate the remaining free valences.

Two things in this module surprise people. First, three of the functional groups a caller can
name are **silently substituted** for something else. Second, hitting a requested O/C ratio can
require placing oxygen somewhere the caller did not ask for. Both behaviours are correct, both
are load-bearing, and neither is discoverable from the call site — which is why they are
specified here rather than left to the docstring.

## Supported Functional Groups <!-- rq-50c742cf -->

The groups that place as named:

- `phenolic` — Ar-OH, 1 oxygen
- `hydroxyl` — Ar-OH, identical to phenolic on a pure PAH, 1 oxygen
- `carboxyl` — Ar-C(=O)(OH), 2 oxygens
- `ether` — Ar-O-Ar bridge, 1 oxygen, consumes 2 edge sites
- `aliphatic_hydroxyl` — an -OH on an sp3 carbon, 1 oxygen
- `amino` — Ar-NH2, 0 oxygens, 1 nitrogen
- `thiol` — Ar-SH, 0 oxygens, 1 sulfur
- `thioether` — Ar-S-Ar bridge, 0 oxygens, 1 sulfur, consumes 2 edge sites

## Groups That Fall Back to Phenolic <!-- rq-a9492b86 -->

`carbonyl`, `quinone`, and `lactone` are accepted by the API and **placed as phenolic instead**.

A pure aromatic PAH edge carbon has one free valence. A carbonyl needs the ring carbon to give
up two, a quinone needs two such carbons in a specific para relationship, and a lactone needs a
fused ring closure. None of these exist on the skeletons this generator builds, so the groups
cannot be realised as named.

The fallback is deliberate — refusing the request outright would break callers that pass a
literature-derived functional-group distribution containing carbonyl. But it means **the
composition a caller requests is not always the composition they receive**, and nothing in the
return value announces the substitution. A caller that asks for carbonyl and later assumes a C=O
is present in the topology will be wrong.

```gherkin
Feature: Substitute unrealisable functional groups rather than failing

  @rq-d3ef932d
  Scenario: A carbonyl request places a phenolic group
    Given a generator configured with functional_groups containing carbonyl
    When the structure is generated
    Then the structure carries a phenolic hydroxyl in its place
    And generation does not raise

  @rq-7d3e9432
  Scenario: A quinone request places a phenolic group
    Given a generator configured with functional_groups containing quinone
    When the structure is generated
    Then the structure carries a phenolic hydroxyl in its place

  @rq-2f5ab3f8
  Scenario: A lactone request places a phenolic group
    Given a generator configured with functional_groups containing lactone
    When the structure is generated
    Then the structure carries a phenolic hydroxyl in its place

  @rq-fbb6a9c0
  Scenario: The oxygen count still matches the request
    Given a request for a group that falls back to phenolic
    When the structure is generated
    Then the total oxygen count reflects the substituted group, not the requested one
```

## Two Placement Modes <!-- rq-be69350a -->

**Exact-count mode.** When `functional_group_preference` is a non-empty dict, it is taken as an
exact specification: `{"phenolic": 3, "carboxyl": 1}` means three phenolic groups and one
carboxyl. If the requested counts exceed the sites actually available, a warning is emitted
rather than silently placing fewer.

**Ratio mode.** Otherwise the target oxygen count is `round(num_carbons × target_O_C_ratio)`.
Aromatic edge sites are filled with phenolic groups first, and any shortfall spills onto sp3
carbons as aliphatic hydroxyls.

## Aliphatic Oxygen Spill <!-- rq-e9a1bf8b -->

The spill exists because of a specific, reproducible failure. A low-aromaticity char carries
aliphatic decoration and hydrogen saturation that consume most of its aromatic edge sites. Once
those are gone there is nowhere left to put a phenolic group, so the structure cannot reach its
oxygen target **at all** — not approximately, but not at all. Every point in that region of
parameter space failed strict validation on every seed.

Spilling the shortfall onto sp3 carbons as aliphatic hydroxyls reaches the target. The behaviour
is controlled by `allow_aliphatic_oxygen` (default true) and applies **only in ratio mode** —
exact-count mode places exactly what was named.

```gherkin
Feature: Reach the oxygen target when aromatic edges run out

  @rq-4fa5b8be
  Scenario: Shortfall spills onto sp3 carbons
    Given a low-aromaticity structure whose aromatic edge sites cannot hold the target oxygen count
    And allow_aliphatic_oxygen is enabled
    When oxygens are assigned by O/C ratio
    Then the remaining oxygens are placed as hydroxyls on sp3 carbons
    And the achieved O/C ratio reaches the target

  @rq-c2a9e7d9
  Scenario: The spill can be disabled
    Given a structure generated with allow_aliphatic_oxygen disabled
    When oxygens are assigned by O/C ratio
    Then no aliphatic hydroxyls are placed
    And the achieved O/C ratio may fall short of the target

  @rq-ad6d4b29
  Scenario: Exact-count mode does not spill
    Given a generator configured with an explicit functional_groups dict
    When oxygens are assigned
    Then only the named groups are placed
    And no aliphatic hydroxyl appears that was not requested

  @rq-f16a145c
  Scenario: A high oxygen target is reachable at low aromaticity
    Given a structure targeting an O/C ratio at or above 0.2 with low aromaticity
    When the structure is generated
    Then strict validation passes without falling back
```

## Ether Bond Types Must Be Repaired After Sanitisation <!-- rq-29741bdc -->

RDKit's `SanitizeMol` marks ether C–O bonds as `AROMATIC` when the ether bridges two aromatic
rings. That is wrong — the C–O bond is single — and it propagates into atom typing and the
exported topology.

`_fix_heteroatom_bond_types` corrects this, and **must be called after any `SanitizeMol` pass
that touches a molecule containing ether oxygens**. It is called in three places: after
assignment, after geometry, and after validation. Each call is load-bearing because each of
those stages can re-sanitise.

```gherkin
Feature: Keep ether bonds single through sanitisation

  @rq-44655f38
  Scenario: An ether C-O bond survives sanitisation as a single bond
    Given a structure containing an Ar-O-Ar ether bridge
    When the molecule is sanitised during generation
    Then both C-O bonds of the ether are single, not aromatic

  @rq-b8fce9a6
  Scenario: The repair runs after every sanitising stage
    Given a structure containing an ether bridge
    When generation completes
    Then the exported topology contains no aromatic C-O bond
```

## Hydrogen Saturation <!-- rq-fc4fac56 -->

`HydrogenAssigner` fills every remaining free valence after functional groups are placed. The
achievable H/C ratio is bounded by structure, not by request: a fully aromatic flake can only
carry hydrogen on its perimeter, so small molecules have a hard H/C ceiling that no
configuration can exceed. Tests must target achievable ratios rather than derived ones.

```gherkin
Feature: Saturate remaining valences

  @rq-85c279ce
  Scenario: No free valence remains after assignment
    Given a structure whose functional groups have been placed
    When hydrogens are assigned
    Then no atom has an unfilled valence

  @rq-8176d084
  Scenario: Small molecules are bounded by their perimeter
    Given a fully aromatic flake of 14 or fewer carbons
    When hydrogens are assigned
    Then the achieved H/C ratio does not exceed the perimeter ceiling
```

## Cross-references <!-- rq-157f1ec7 -->

- The hydrogen bonds between adjacent hydroxyls placed here are what
  `geometry-embedding.md` must not report as steric clashes.
- Atom types for the groups placed here are specified in `opls-typing.md`.
- The aliphatic-carbon population this module draws on is produced by the H/C-aware growth
  described in the carbon skeleton, documented in `docs/solutions/` pending its own `rqm/` file.
