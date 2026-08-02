# Feature: Surface Stacking <!-- rq-63c93813 -->

`src/biochar/workflows/surface_builder.py` builds a multi-sheet biochar surface: several PAH
sheets, each generated through the single-molecule pipeline, arranged into a pore geometry and
exported as one GROMACS system.

Two geometries exist. A **slit pore** stacks parallel sheets along z at a fixed separation. An
**amorphous** pore places randomly rotated sheets into a box until none are closer than a minimum
separation, which is a search and can fail.

The recurring hazard here is that a surface is a physical claim about pore structure, and the
files carry that claim in three separate places — the coordinates, the box, and the title written
into the `.gro`. When the builder cannot deliver the geometry it was asked for, those three can
disagree, and the file is the only thing a reader has weeks later. So the requirements below pin
what the geometry *is*, and pin that the exported files say the same thing about it.

## Feature API <!-- rq-c39502d1 -->

- `SurfaceConfig(target_num_carbons: int = ..., num_sheets: int = 2, pore_type: str = "slit", pore_diameter: float = 10.0, min_separation: float = 3.0, max_attempts: int = 500, amorphous_fallback: Optional[str] = None, box_padding_xy: float = 1.0, box_padding_z: float = 1.0, sheet_overrides: Optional[List[Dict]] = None, sheet_base_name: str = "SHT", seed: Optional[int] = None, ...)` <!-- rq-d0631f70 -->
  - `pore_diameter` is the gap between the **inner van der Waals surfaces** of adjacent sheets,
    not the distance between their centroids.
  - `num_sheets` is at least 2 — one sheet is a molecule, not a surface.
  - `sheet_base_name` is at most 3 characters, so `<name><index>` stays inside the `.gro`
    5-character residue field.
  - Validated on construction; an unusable combination raises rather than producing a surface that
    does not match the request.

- `SurfaceBuilder(config: SurfaceConfig)` <!-- rq-169cc55e -->
  - `build() -> Tuple[List[SheetResult], np.ndarray]` — generates and positions the sheets,
    returning them with the box vectors in nanometres.
  - `export_gromacs(output_directory, basename="surface") -> Tuple[Path, Path, List[Path]]` —
    writes the combined `.gro`, the `.top`, and one `.itp` per distinct sheet.

- `SheetResult` <!-- rq-464bbae1 -->
  - One sheet's `mol`, `coords` (ångström), `atom_types`, `charges`, `composition`, and
    `molecule_name`. Sheets that are copies of one another do not share mutable state.

- `generate_surface(...) -> Tuple[Path, Path, List[Path]]` <!-- rq-c7ad2878 -->
  - One-call convenience over the two above. Exposes the options that change the geometry a caller
    receives, including how a failed amorphous packing is handled.

## Sheet Separation Is Measured Between Surfaces, Not Centres <!-- rq-ad61986f -->

Adjacent sheets in a slit pore sit `pore_diameter + CARBON_VDW_DIAMETER` apart, centroid to
centroid — 3.4 Å being the graphite interlayer spacing, which is the thickness a flat aromatic
sheet occupies.

The distinction is the whole meaning of the parameter. A caller asking for a 10 Å pore wants 10 Å
of space for water and sorbates, not a 10 Å centroid separation that leaves 6.6 Å of usable gap.

```gherkin
Feature: Space sheets by the gap the caller asked for

  @rq-dbd8829a
  Scenario: Adjacent sheets are separated by the pore diameter plus the sheet thickness
    Given a slit pore with a requested pore diameter
    When the sheets are positioned
    Then each adjacent pair of centroids is separated by the pore diameter plus 3.4 Angstrom
```

## The Exported Files Describe the Geometry That Was Built <!-- rq-028b78d5 -->

Amorphous packing is a search: sheets are randomly rotated and placed until none violate the
minimum separation, and with enough sheets in too small a box it does not converge.
`amorphous_fallback="slit"` asks for a slit stack instead of an exception — a reasonable request,
because a slit surface is still useful and losing the whole run is worse.

What the caller must not receive is a slit stack that describes itself as amorphous. The `.gro`
title is the only human-readable record of what the geometry is, and it is read weeks later by
someone deciding whether a trajectory answers their question.

```gherkin
Feature: Describe the geometry actually produced

  @rq-0a48ccd8
  Scenario: A surface that fell back to slit geometry is not labelled amorphous
    Given an amorphous request that cannot be packed
    And a fallback to slit geometry
    When the surface is exported
    Then the structure file describes a slit pore

  @rq-4e6abc16
  Scenario: The fallback is reported to the caller
    Given an amorphous request that cannot be packed
    And a fallback to slit geometry
    When the surface is built
    Then a warning names the geometry that was substituted

  @rq-2c190d16
  Scenario: Without a fallback an unpackable request raises
    Given an amorphous request that cannot be packed
    And no fallback configured
    When the surface is built
    Then the build raises rather than returning another geometry
```

## Copied Sheets Share Nothing Mutable <!-- rq-9caef4cc -->

When every sheet is the same — no per-sheet overrides, no pH — the builder generates one and
copies it, which is a large saving on a many-sheet surface.

A copy is a copy. Two sheets that share a mutable object are one sheet wearing two names: a caller
that adjusts one sheet's composition record, or a later pipeline stage that annotates it, silently
changes the others.

```gherkin
Feature: Copy identical sheets without sharing state

  @rq-9cdd2239
  Scenario: Copied sheets do not share their composition record
    Given a surface whose sheets are identical
    When the sheets are built
    Then no two sheets refer to the same composition object

  @rq-1f914136
  Scenario: Copied sheets do not share coordinates
    Given a surface whose sheets are identical
    When one sheet's coordinates are modified
    Then the other sheets are unchanged
```

## Per-Sheet Protonation Consumes a Seed Range <!-- rq-0569cb45 -->

Setting `pH` makes each sheet an independent sample of a protonation state, so the identical-sheet
optimisation cannot apply and each sheet is generated separately. Each takes a seed derived from
the configured one by adding its index.

That makes a surface reproducible from its seed, and it means one surface consumes a **range** of
seeds rather than a single value. Two ensemble members at consecutive seeds would otherwise share
most of their sheets while appearing independent.

```gherkin
Feature: Derive per-sheet seeds without silently overlapping ensembles

  @rq-6ba78ac5
  Scenario: A protonated surface is reproducible from its seed
    Given a surface built at a pH with a fixed seed
    When it is built twice
    Then both builds produce identical coordinates

  @rq-682a21cd
  Scenario: The seeds a surface consumes are documented for the caller
    Given a surface of N sheets built at a pH from seed S
    Then it consumes seeds S through S plus N minus 1
```

## The Topology and the Include Files Agree on the Sheet Count <!-- rq-fc302462 -->

Identical sheets are written as one `.itp` referenced with `count = N` in `[ molecules ]`.
Distinct sheets are written as N `.itp` files, each referenced once.

The two halves are decided by the same predicate, and they must be: a `.top` that lists one
moleculetype `N` times while only one `.itp` was written is correct, but a `.top` written for
distinct sheets against a single `.itp` leaves `grompp` unable to resolve the moleculetypes it
names.

```gherkin
Feature: Keep the topology consistent with the include files written beside it

  @rq-17ffc86b
  Scenario: Identical sheets are written once and referenced many times
    Given a surface whose sheets are identical
    When it is exported
    Then exactly one include file is written
    And the topology references it with a count equal to the sheet count

  @rq-884c2220
  Scenario: Distinct sheets each get their own include file
    Given a surface whose sheets differ
    When it is exported
    Then one include file is written per sheet
    And every moleculetype the topology names is defined in one of them
```

## The Convenience Wrapper Reaches the Options That Change the Geometry <!-- rq-ff5fd763 -->

`generate_surface` is the entry point the documentation points at, and a caller who never
constructs a `SurfaceConfig` directly should still be able to reach the options that decide what
geometry they receive — `amorphous_fallback` above all, since without it an amorphous request that
cannot be packed raises instead of degrading.

An option that exists only on the configuration object is unavailable to most callers in practice.

```gherkin
Feature: Expose geometry-determining options through the one-call entry point

  @rq-bfb7e848
  Scenario: The fallback behaviour is reachable without building a config
    Given a caller using the convenience function
    When an amorphous surface is requested with a slit fallback
    Then the fallback applies
```

## Cross-references <!-- rq-6dd192a1 -->

- Each sheet is produced by the single-molecule pipeline; the composition targets, the aromatic
  floor, and strictness are specified in `generation-config.md`.
- The `.gro`, `.itp`, and `.top` writers, the ångström-to-nanometre conversion, and the
  5-character residue limit are specified in `gromacs-export.md`.
- Clash warnings on the flat sheets that make up a surface are specified in
  `geometry-embedding.md`.
