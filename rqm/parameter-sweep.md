# Feature: Parameter Sweep <!-- rq-80dfec16 -->

`workflows/sweep.py` turns a hand-run loop into a declarative pipeline. A sweep config names a
set of axes; the driver expands their Cartesian product, builds every point through the
single-molecule pipeline, writes GROMACS files into an organised tree, and emits a **manifest**
recording what each point actually became.

The manifest is the module's real product. The structures are on disk either way, but a sweep of
sixty points is read through its manifest — by `export/md_setup.py`, which turns rows into run
directories, and by a person deciding months later which trajectories answer their question.
Everything here is therefore a requirement about **the manifest telling the truth about the files
beside it**: that each row names a distinct structure, that the name in the row is the name inside
the topology, that a status means one thing, and that the row records enough to say which version
of this package produced it.

A sweep is also the layer where failure stops being exceptional. Individual points fail routinely —
a composition target that a given seed cannot reach — and the driver's job is to keep going while
recording precisely what happened. The failure mode is a point that is quietly recorded as
something it is not.

## Feature API <!-- rq-4a9e69bb -->

- `expand_grid(axes, fixed=None, name_template="BC{i:03d}") -> List[GridPoint]` <!-- rq-56adfdbe -->
  - Expands the axes into one `GridPoint` per combination. An axis value overrides the same key in
    `fixed`. A key that is not a `GeneratorConfig` field is refused here, before anything is built.
  - Every point in the returned grid is distinguishable from every other: no two share a molecule
    name.

- `build_point(point, output_root, base_seed=0, max_retries=8, on_validation_fail="fallback", ...) -> PointResult` <!-- rq-bbb23eaa -->
  - Builds one point, retrying seeds while validation fails and then applying the chosen failure
    mode. Returns a result rather than raising, except where the failure mode says otherwise.

- `run_sweep(config, output_directory=None, progress_callback=None, quiet=True) -> Dict` <!-- rq-50e10082 -->
  - Expands, builds every point, writes `manifest.csv` and `manifest.json`, and returns a summary.

- `load_sweep_config(path) -> Dict` <!-- rq-5450da17 -->
  - Reads a sweep config from YAML or JSON.

- `GridPoint` <!-- rq-0bf15b93 -->
  - One point: its index, its filesystem label, its molecule name, the axis values that vary at
    this point, and the full keyword set the generator is configured from.

- `PointResult` <!-- rq-7081e03f -->
  - What a point became: status, the seed that produced it, the achieved composition, the
    validation report, and the paths written. `to_row()` flattens it to a manifest row.

- `POINT_STATUSES` <!-- rq-7e5c5b4c -->
  - The complete set of statuses a manifest row can carry. A consumer filtering rows compares
    against this set rather than against a list copied from prose.

- `SweepError` <!-- rq-78db3aec -->
  - Raised for a malformed sweep configuration, and for a sweep that its failure mode says must
    not continue.

## Every Point in the Grid Is Distinguishable From Every Other <!-- rq-3cf07edf -->

The molecule name is a GROMACS residue name, so it is capped at five characters. Names are
generated from a template, and a template plus a cap is a hashing function: distinct points can
collide.

A collision is not cosmetic. The manifest carries one row per point, each naming its molecule; two
rows carrying the same name describe two different structures that nothing downstream can tell
apart. Points 100 and 1000 of a default-named grid are the clearest case — `BC{i:03d}` renders
`BC0100` and `BC1000`, and both truncate to `BC100`.

The grid is expanded before anything is built, which is where the collision is cheap to detect and
where the caller can still change the template.

```gherkin
Feature: Give every grid point its own name

  @rq-2a033e95
  Scenario: No two grid points share a molecule name
    Given a grid large enough that the name template collides under the five-character cap
    When the grid is expanded
    Then the expansion is refused rather than returning two points of the same name
```

## The Name a Point Is Recorded Under Is the Name in Its Topology <!-- rq-64a42a92 -->

`molecule_name` is a `GeneratorConfig` field, so a sweep config can set it in `fixed` — where it
applies to every point, overriding the per-point name the template produced.

The point still reports the templated name, so the manifest says `BC000` while the topology beside
it says something else entirely. Every consumer that joins a manifest row to its files by name is
then joining on a name that does not appear in them, and `md_setup` writes run directories for
structures whose residue names it has never seen.

A name that varies per point is not a parameter a sweep can fix. Setting it is a configuration
error, and it is refused where the rest of them are.

```gherkin
Feature: Keep the recorded name and the written name the same

  @rq-4156245d
  Scenario: A molecule name fixed across the sweep is refused
    Given a sweep config that sets a molecule name in its fixed parameters
    When the grid is expanded
    Then the configuration is refused as one the per-point name cannot survive

  @rq-a118d6d2
  Scenario: The name a point reports is the name its generator is configured with
    Given an expanded grid
    When each point's reported name is compared with the name in its generator configuration
    Then they agree for every point
```

## A Seed Retry Is a Variance Remedy, Not a Repair <!-- rq-1bdc6066 -->

Generation is stochastic: functional-group placement, defect draws and embedding all consume the
seed, so a composition target that one seed misses another may reach. That is what the retry loop
is for, and it is why the default sweep completes rather than stopping at the first miss.

It is only that. A request that no seed can satisfy fails identically at every seed, and running it
eight times produces eight identical failures more slowly. Two kinds of request are in that
position and both stop the loop:

- one that raises something other than a validation failure — an unbuildable skeleton, a malformed
  request — where the exception ends the attempt loop immediately and is recorded with its text;
- one that fails validation with **the same errors it failed with last time**. Two seeds producing
  a byte-identical report is the loop's evidence that the failure does not depend on the seed.
  Asking for five hundred carboxyl groups on a twenty-carbon skeleton reports
  `carboxyl 12/500` at every seed, word for word.

The distinction is visible in the manifest. `n_attempts` on a failed row says whether the driver
believed variance was worth exploring, and a deterministic failure that reported eight attempts
would be claiming a search it never made.

```gherkin
Feature: Retry only what a different seed could change

  @rq-bedc92d0
  Scenario: A validation failure is retried at a new seed
    Given a point whose first seeds fail validation and a later one passes
    When the point is built
    Then the passing seed is recorded and the point is marked as passing strictly

  @rq-d600f3fc
  Scenario: A failure that is not a validation failure is not retried
    Given a point whose configuration raises something other than a validation failure
    When the point is built
    Then only one attempt is recorded and the failure text is carried into the result

  @rq-a47e48e0
  Scenario: A validation failure that repeats unchanged stops the retry loop
    Given a point whose validation report is identical at every seed
    When the point is built with a large retry budget
    Then the loop stops once the report repeats rather than spending the budget
```

## Each Failure Mode Is a Different Answer <!-- rq-af33646e -->

`on_validation_fail` offers three modes, and a caller chooses between them because they want three
different things:

- `fallback` — keep the structure that strict validation rejected, record why, and carry on. The
  default, because a sweep exists to produce a grid and a partly-out-of-tolerance point is usually
  still worth having.
- `skip` — do not keep it. The point is recorded as **skipped**, carrying its validation errors,
  and no files are written.
- `strict` — do not accept a grid containing it. The sweep is refused, naming the point and what
  was wrong with it.

`strict` is the one that must not merely be a synonym. A caller who selects it over `skip` is
saying the grid is only meaningful if every point is in tolerance, and the answer they need is an
error rather than a manifest they must remember to re-check.

A skipped point is also not a failed one. `failed` says something went wrong that no configuration
asked for — a crash, an unbuildable request — and it is the row a reader investigates. A point the
sweep declined to keep because its own config said not to is an expected outcome of a healthy run,
and collapsing the two makes a working sweep indistinguishable from a broken one at a glance.

```gherkin
Feature: Give each failure mode its own behaviour

  @rq-0cf2e5d4
  Scenario: The default mode keeps a structure that strict validation rejected
    Given a point that no seed builds within tolerance
    When the point is built under the fallback mode
    Then the structure is written and recorded as a fallback carrying its validation errors

  @rq-086427f6
  Scenario: Skip records the point as skipped and writes nothing
    Given a point that no seed builds within tolerance
    When the point is built under the skip mode
    Then no structure files are written and the point is recorded as skipped with its errors

  @rq-626e9fb8
  Scenario: Strict refuses the sweep rather than recording a row
    Given a point that no seed builds within tolerance
    When the point is built under the strict mode
    Then the failure is raised, naming the point and its validation errors
```

## The Manifest Records What Produced It <!-- rq-826be977 -->

A manifest is read long after the run, and the statuses in it are not stable across releases of
this package. `strict_pass` in particular has meant two different things: a functional-group
shortfall used to satisfy strict validation and no longer does, so an old manifest and a new one
can describe the same grid with different counts in every column.

Nothing in the manifest says which meaning applies. The axes, the seed, the retry budget and the
failure mode are all recorded — everything about the request, and nothing about the code that
answered it. The package version is the one field that makes the rest interpretable.

```gherkin
Feature: Record the code that produced the manifest

  @rq-43924feb
  Scenario: The manifest names the package version that wrote it
    Given a completed sweep
    When its manifest is read
    Then the recorded version is the version of the package that ran
```

## A Status Means One Thing <!-- rq-1ae3a84a -->

Every manifest row carries a status, and consumers branch on it: `md_setup` writes a run directory
for a row that produced files and reports back any row that did not.

That branch is only as good as the vocabulary both sides share. A status named in the
documentation but never written is a filter that silently matches nothing, and a status written but
never documented is a row a consumer has no rule for. The set is small enough to state exactly, and
is stated once in the module rather than transcribed into prose that drifts from it.

```gherkin
Feature: Keep one vocabulary of statuses

  @rq-2da595b7
  Scenario: Every status a build can return is one the documented set names
    Given the statuses the point builder can assign
    When they are compared with the module's declared set
    Then the two agree exactly

  @rq-7f22837b
  Scenario: The documented status table names the statuses the module declares
    Given the sweep status table in the project README
    When the statuses it names are compared with the module's declared set
    Then the two agree exactly
```

## Cross-references <!-- rq-8d8ebbfe -->

- The strict-validation semantics a point is judged against, including the functional-group
  shortfall rule that changed what `strict_pass` means, are specified in `generation-config.md`.
- The consumption of a manifest into GROMACS run directories, and the provenance each run
  directory records, are specified in `md-setup.md`.
- The five-character residue-name limit the molecule name is capped by is specified in
  `gromacs-export.md`.
