# Feature: Command-Line Entry Points <!-- rq-b9751811 -->

The four console scripts declared in `pyproject.toml`'s `[project.scripts]` — `biochar-gen`,
`biochar-sweep`, `biochar-md-setup`, `biochar-condense` — are the only interface most users of this
package ever touch. Each one parses arguments and maps them onto a library call.

They were long described as holding no logic of their own. That is not quite true, and the gap is
where this document lives. Three decisions are made at this layer and nowhere else:

- **Precedence** between a config file and the command line. The library sees one resolved
  `GeneratorConfig` and cannot tell which of its fields the user asked for.
- **Exit status.** The library raises or returns; turning that into a number a shell can branch on
  is the entry point's job, and a pipeline chaining these commands has nothing else to read.
- **Cross-argument consistency.** Two flags that are individually valid can be jointly impossible.
  The entry point is the only place that sees both before any work starts.

The failure mode throughout is the same one `carbon-skeleton.md` names for parameters: a request
accepted and quietly not honoured. At this layer it has a second form — a command that reports
success for a run that produced nothing.

## Feature API <!-- rq-b31f2ffe -->

- `biochar.cli.cli.main(argv=None) -> int` <!-- rq-fa8c5dbc -->
  - Generates a single structure. Returns `0` on success, `1` on a configuration or generation
    failure.
  - `_build_parser() -> argparse.ArgumentParser` builds the parser separately so the accepted
    surface can be inspected without running anything.

- `biochar.cli.sweep_cli.main(argv=None) -> int` <!-- rq-46f0dae6 -->
  - `run` executes a declarative grid; `template` writes a starter config to stdout. Returns `0`
    when every point built, `2` when any point failed, `1` on a config error.

- `biochar.cli.md_setup_cli.main(argv=None) -> int` <!-- rq-a2c4ed0e -->
  - Turns a sweep manifest into GROMACS run directories. Returns `0` when at least one directory
    was written, `2` when none was, `1` on a setup error.

- `biochar.cli.condensation_cli.main(argv=None) -> int` <!-- rq-c46e42dc -->
  - `generate` builds a molecule then sets up its condensation run; `from-files` condenses an
    existing pair of files. Returns `0` on success, `1` on an inconsistent request.

Every `main` takes `argv` explicitly and returns a status rather than calling `sys.exit`, so the
behaviour above is reachable from a test without a subprocess.

## A Saved Config Reloads as the Config That Was Saved <!-- rq-a71b7f66 -->

`--save-config` and `--load-config` exist to make a run reproducible: save what was resolved, reload
it later, get the same structure. That promise requires the pair to round-trip.

Reloading must therefore preserve every field of the saved file. The field that matters most is
`seed` — it is the entire mechanism by which a run is reproducible, and it is exactly the kind of
field a user does not re-type on the reload because the config file is supposed to carry it.

`temperature` and `feedstock` matter for a second reason. The composition fields they derive —
`H_C_ratio`, `O_C_ratio`, `aromaticity_percent` — are written into the saved config as numbers. A
reload that keeps the numbers but drops the temperature they came from produces a config that still
generates the right structure while no longer recording what that structure rests on, which is the
provenance failure `temperature-model.md` was written to prevent.

```gherkin
Feature: Round-trip a saved configuration

  @rq-b4e2b5cc
  Scenario: A reloaded config matches the one that was saved
    Given a config saved from a run with a non-default carbon count, name, seed, temperature and feedstock
    When it is reloaded with no other arguments and saved again
    Then every field of the second config equals the first

  @rq-e2c50b04
  Scenario: Derived composition keeps the conditions it was derived from
    Given a config saved from a run driven by temperature and feedstock
    When it is reloaded with no other arguments
    Then the reloaded config still records that temperature and feedstock
```

## A Default Is Not a Request <!-- rq-e6faa000 -->

A flag the user did not type carries no instruction. When a config file is loaded, only the flags
actually present on the command line override it; the rest of the parser's defaults stand aside.

Without that rule every field with a default silently overwrites the loaded config with the
default's value, and the overwrite is invisible: the user sees a config file being loaded, no
warning, and a structure built from something else. It is the parameter-accepted-and-ignored failure
in its most costly form, because the ignored parameter is the entire file.

The rule replaces a three-way split — fields that always overrode, fields that overrode only when
non-`None`, and `pH` — with one that holds for every field. Where no config is loaded, the parser's
defaults apply exactly as before; a default is a fallback, not an override.

```gherkin
Feature: Distinguish a supplied flag from an unsupplied default

  @rq-e6099cd4
  Scenario: An unsupplied flag leaves the loaded value alone
    Given a config file specifying a non-default carbon count
    When the config is loaded without a carbon-count flag
    Then the loaded carbon count is used

  @rq-b8fb0aeb
  Scenario: A supplied flag overrides the loaded value
    Given a config file specifying a carbon count
    When the config is loaded together with a different explicit carbon count
    Then the command-line count is used

  @rq-1a4b6dae
  Scenario: A flag supplied at its own default value still overrides
    Given a config file specifying a non-default carbon count
    When the config is loaded together with a carbon count equal to the parser default
    Then the command-line count is used

  @rq-b81e79b5
  Scenario: Without a config file the parser defaults still apply
    Given no config file and no composition flags
    When a structure is generated
    Then the documented default carbon count, name and charge method are used
```

## Functional Groups Are Replaced as a Set <!-- rq-fd1e5f16 -->

The six group-count flags describe one dictionary between them. Supplying any of them replaces a
loaded config's group dictionary entirely rather than merging into it.

Merging would produce a request the user never made: asking for two phenolics against a config
holding one carboxyl would silently generate both. The counts interact — they compete for the same
edge sites — so a set assembled from two sources is not a composition anyone chose. Replacement
keeps the rule the help text already states, that explicit group counts displace the ratio-driven
placement, and extends it to the config file.

```gherkin
Feature: Treat the group flags as one request

  @rq-92f2fb37
  Scenario: Any group flag replaces the loaded group set
    Given a config file specifying one functional group
    When the config is loaded together with a different explicit group flag
    Then only the group named on the command line is requested

  @rq-c72f3fd1
  Scenario: No group flag leaves the loaded group set intact
    Given a config file specifying functional groups
    When the config is loaded with no group flags
    Then the loaded groups are requested
```

## A Command That Produced Nothing Does Not Report Success <!-- rq-1eae7bef -->

`biochar-md-setup` skips any manifest row without usable structure files, which is normal — a sweep
point that failed legitimately has none. A run in which *every* row was skipped is different: the
command wrote no run directories at all, and there is nothing to submit.

Exiting `0` there tells a pipeline the stage succeeded. The next stage then submits an empty
directory tree, and the failure surfaces at the scheduler rather than at the command that had the
information. `biochar-sweep` already distinguishes this case with `2`; the manifest consumer uses
the same code for the same meaning, so a shell can branch on it without knowing which command it
ran.

Partial skips stay at `0`. The signal is "produced nothing", not "produced less than requested" —
skips are the expected shape of a real sweep and an exit code that fires on them would be ignored
within a week.

```gherkin
Feature: Report an empty result as an empty result

  @rq-2b6f2c93
  Scenario: A manifest yielding no run directories exits non-zero
    Given a manifest whose every row lacks structure files
    When run directories are set up from it
    Then no run directory is written
    And the exit status is not success

  @rq-59fb4726
  Scenario: A manifest yielding some run directories succeeds
    Given a manifest with one usable row and one unusable row
    When run directories are set up from it
    Then the usable row is written
    And the exit status is success
```

## A Setup Run Does Not Point at a Repeat It Will Not Produce <!-- rq-a48afb56 -->

`biochar-condense` writes `run_condensation.sh`, which runs `--repeats` independent annealing
repeats, and `run_surface.sh`, which lifts one of them — `--which-repeat` — into a surface. The two
numbers are only ever seen together here.

An out-of-range repeat is accepted today and written into the surface script as a path that the
condensation run will never create. Nothing fails at setup time. It fails after the annealing run,
which is the expensive part, at the moment the user reaches for the result — the worst possible
place to learn that a flag was wrong, since setup is instant and re-running it costs nothing.

The check belongs at the entry point and before any generation: `generate` builds a molecule first,
so validating afterwards would still waste the build.

```gherkin
Feature: Refuse a surface request the condensation run cannot satisfy

  @rq-c53b3aeb
  Scenario: A repeat above the repeat count is refused
    Given a condensation setup requesting three repeats
    When the surface is requested from repeat seven
    Then the request is refused before any structure is generated
    And no run directory is written

  @rq-a4f5a7d7
  Scenario: A repeat below one is refused
    Given a condensation setup requesting three repeats
    When the surface is requested from repeat zero
    Then the request is refused

  @rq-2ed54e2c
  Scenario: A repeat within range is accepted
    Given a condensation setup requesting three repeats
    When the surface is requested from repeat three
    Then the setup is written and references that repeat
```

## The Starter Template Names Only Options That Exist <!-- rq-eb1c8214 -->

`biochar-sweep template` writes a config a user is expected to edit and run. Its comments are the
only documentation most users will read for the sweep format, so an option named there is one they
will type.

A comment naming a charge backend the config rejects turns the template into a source of errors
rather than a starting point, and the same comment omitting a backend that does exist hides it. The
template's own defaults are exercised whenever it is run; its comments are not, so they need a test
that reads them against the accepted values rather than a reviewer noticing.

```gherkin
Feature: Keep the template honest about what it offers

  @rq-cba5c290
  Scenario: Every charge method the template names is accepted
    Given the starter sweep config the template command writes
    When each charge method named in it is used to build a configuration
    Then none is rejected

  @rq-32e0e0e6
  Scenario: The template runs as written
    Given the starter sweep config the template command writes
    When it is parsed as a sweep configuration
    Then it expands to a grid without error
```

## A Failure Is Reported on stderr and in the Status <!-- rq-ac71d20f -->

Every entry point writes its diagnostics to stderr and its results to stdout, and reports failure in
both the message and the exit status.

Splitting the streams is what makes the commands composable: `biochar-sweep template > sweep.yaml`
depends on progress lines not landing in the file, and `biochar-gen`'s printed output paths are
meant to be read by the next command in a pipeline. A diagnostic on stdout corrupts both silently.

```gherkin
Feature: Keep diagnostics off stdout

  @rq-b697d2e6
  Scenario: A missing config file is reported and fails
    Given a config path that does not exist
    When a structure is requested from it
    Then the failure is reported on stderr
    And the exit status is a failure

  @rq-c8d1f9e0
  Scenario: The template command writes only the template to stdout
    Given the template subcommand
    When it is run
    Then stdout holds a parseable config and nothing else
```

## Cross-references <!-- rq-fb2e2ae7 -->

- The `GeneratorConfig` fields these flags resolve to, the composition defaults, and the strictness
  that decides whether a structure is accepted are specified in `generation-config.md`.
- The temperature and feedstock model behind `--temperature` / `--feedstock` is specified in
  `temperature-model.md`.
- The sweep config format, grid expansion and the manifest `biochar-md-setup` consumes are specified
  in `parameter-sweep.md`.
- The run directories, ion profiles and stage ordering `biochar-md-setup` writes are specified in
  `md-setup.md`.
- The condensation schedule, packing and surface stage `biochar-condense` sets up are specified in
  `condensation-annealing.md`.
