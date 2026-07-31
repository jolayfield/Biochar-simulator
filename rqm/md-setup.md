# Feature: GROMACS Run Directory Setup <!-- rq-1433518c -->

`src/biochar/export/md_setup.py` turns generated structures into runnable GROMACS jobs. Given a
sweep manifest it writes one directory per structure containing the `.mdp` files for every stage,
copies of the structure and topology, and a driver script that walks the eight-stage protocol:
dry energy minimisation, annealing under NVT then NPT, a final NPT, then box padding, solvation,
ion placement, and wet EM/NVT/NPT.

What this module emits is not read by Python again. It is `.mdp` text and shell script that a
scheduler runs hours later, often on another machine, and the first reader is usually `grompp` or
`mdrun`. A mistake here surfaces as a job that dies at 3 a.m., or worse as one that completes and
reports a plausible density having simulated something other than what was asked for.

That distance is the constraint the requirements below are written against. The generated script
runs under `set -euo pipefail`, so a stage that references a file an earlier stage did not create
aborts the run rather than continuing; and every physical quantity in the `.mdp` set is a claim
about the protocol being reproduced, not a default that GROMACS would sanity-check.

## Feature API <!-- rq-15d0c9f2 -->

- `setup_md_from_manifest(manifest_csv: str | Path, output_root: str | Path, ion_profile: str = "mn_calcareous_default", config: Optional[MDSetupConfig] = None) -> list[dict]` <!-- rq-86c5ecb3 -->
  - Reads a `workflows.sweep` manifest and writes one run directory per row that produced files.
  - Processes rows whose `status` is `strict_pass` or `fallback`; reports every other row with a
    reason and writes nothing for it.
  - Resolves the structure, topology, and include-file paths recorded in the manifest, trying the
    path as given and then relative to the manifest's own location.
  - Returns one dict per row: `label`, `status`, `run_dir`, `gro_path`, `top_path`,
    `skipped_reason`.

- `setup_one_structure(gro_path: str | Path, top_path: str | Path, output_dir: str | Path, label: str = "structure", config: Optional[MDSetupConfig] = None, pyrolysis_temperature_c: Optional[float] = None, status: Optional[str] = None, include_paths: Optional[list] = None) -> Path` <!-- rq-6f1c55c5 -->
  - Writes the `.mdp` set, the driver script, and copies of the structure, topology, and every
    include file the topology needs.
  - Raises `MDSetupError` when a required input is missing, rather than writing a directory that
    cannot run.
  - Renders the annealing stages for `pyrolysis_temperature_c`. Without one it uses the lowest
    Wood anchor, 400 °C, so an unknown structure gets the mildest schedule rather than one that
    could over-anneal it.
  - Writes `run_provenance.json` recording the label, status, source files, net charge, and the
    annealing schedule used — including whether the temperature came from the manifest or the
    default.

- `MDSetupConfig(solvent_pad_nm: float = 1.2, water_model: str = "spce", ion_profile: str = "mn_calcareous_default", ntomp: int = 8, gmx_bin: str = "gmx", cluster: Optional[str] = None, ...)` <!-- rq-00cf61e3 -->
  - Options for one run directory. `cluster="slurm"` additionally writes a dependent-job chain
    and requires the site-specific submit script, image, and partition fields.
  - `pre_solvation_stage` splices an `insert-molecules` stage in after the final anneal; every
    later stage then uses that stage's merged topology.

- `get_ion_profile(name_or_profile) -> IonProfile` <!-- rq-8d689517 -->
  - Resolves a named background electrolyte, or passes an `IonProfile` through.

- `IonProfile(name: str, ca_mM: float = 0.0, mg_mM: float = 0.0, na_mM: float = 0.0, k_mM: float = 0.0, counter_ion: str = "CL", description: str = "")` <!-- rq-2533ac5e -->
  - Target bulk concentrations. These are representative compositions, not measurements of any
    particular site, and the description says which.

## Every Stage Consumes What an Earlier Stage Produced <!-- rq-4f5a5569 -->

The driver script runs under `set -euo pipefail`. Each stage names its inputs explicitly, and a
stage's inputs are files an earlier stage has already written.

Solvation is where this is easy to get wrong, because `gmx solvate -p` **updates** a topology in
place rather than creating one. The solvation stage therefore needs the run's working topology to
exist before it runs, not after.

```gherkin
Feature: Order the pipeline so every input exists when it is read

  @rq-62e0732a
  Scenario: The topology the solvation step updates exists before that step runs
    Given a rendered pipeline script
    When the solvation stage is reached
    Then the topology it passes to solvate has already been created by an earlier line
```

## Ion Placement Accumulates <!-- rq-2fe208c2 -->

`gmx genion` replaces solvent molecules with ions, reading coordinates from a `.tpr` and writing
a new structure. Placing several species means running it once per species, and each run must see
what the previous one did — both the structure it writes and the `.tpr` the next run is built
from.

Building one `.tpr` from the pre-ion coordinates and reusing it for every species discards each
species as the next is placed, while the topology accumulates all of them. The result is a
structure and a topology that disagree about how many atoms exist, which is one of the warning
classes below.

```gherkin
Feature: Place each ion species on top of the last

  @rq-39ce3988
  Scenario: Every ion species is placed against the coordinates the previous species produced
    Given an ion profile naming more than one cation
    When the ion stages are rendered
    Then each genion call is preceded by a grompp that rebuilds the run input from the structure
         the previous genion wrote
```

## A Suppressed Warning Is Named <!-- rq-93d3a759 -->

`grompp` warnings are refusals to proceed silently, and the two this pipeline can provoke are
both fatal to the science: a net-charged system under PME, and a structure whose atom names or
count disagree with the topology.

`grompp` has no way to suppress a *named* warning — `-maxwarn` takes a count and absorbs whichever
warnings arrive first. That makes the count the whole of the control, and it makes a generous
count indiscriminate: a stage allowed two warnings accepts the two the author had in mind, or any
other two, including ones nobody has seen yet.

So `-maxwarn` appears only on stages where a specific warning is expected, with a comment naming
it, and the count is the number of warnings named. A stage with nothing to excuse passes no
`-maxwarn` at all and fails on the first surprise.

```gherkin
Feature: Suppress only the warnings that were understood

  @rq-dd80065b
  Scenario: A stage with no expected warning suppresses nothing
    Given a rendered pipeline script for a neutral structure
    When the dry stages are read
    Then none of them passes -maxwarn

  @rq-af56cae9
  Scenario: A suppressed warning is named where it is allowed
    Given a rendered pipeline script
    When a grompp stage passes -maxwarn
    Then a comment on that stage names the warning it absorbs
    And the count equals the number of warnings named
```

## The Annealing Schedule Follows the Pyrolysis Temperature <!-- rq-317d0c21 -->

The Wood et al. protocol anneals to a peak temperature that scales with the heat-treatment
temperature of the char being modelled: 1000 K at 400 °C, 2000 K at 600 °C, and 3000 K at 800 °C,
interpolated between those anchors and clamped outside them, with the timestep dropping to 0.5 fs
once the peak exceeds roughly 1500 K.

A fixed peak applies the 400 °C schedule to every structure. An 800 °C char annealed to 1000 K is
not driven through the graphitisation the protocol exists to reproduce, and the run completes and
reports a structure that looks converged.

The pyrolysis temperature reaches setup only when the sweep varied it, since it travels as a
sweep axis. When it is absent the 400 °C anchor applies — the mildest schedule, chosen so an
unknown structure is under-annealed rather than over-annealed — and the run directory records
that the value was a default rather than a measurement.

The scaling has one implementation, `workflows.condensation.anneal_spec_for_htt`. Two
implementations of one published protocol will disagree, and the disagreement will be silent
because both produce runnable input.

```gherkin
Feature: Anneal each structure on the schedule its pyrolysis temperature implies

  @rq-5b7a5ab8
  Scenario: The annealing peak temperature follows the pyrolysis temperature
    Given a manifest row for a structure pyrolysed at 800 degrees Celsius
    When its run directory is written
    Then the annealing peak temperature is the one the Wood scaling gives for 800 degrees
    And it is not the 400 degree value

  @rq-997c1640
  Scenario: The timestep drops with the schedule it belongs to
    Given a structure whose annealing peak exceeds 1500 K
    When its run directory is written
    Then the annealing timestep is 0.5 femtoseconds
```

## The Wet Stage Holds the Sheet and Relaxes the Solvent <!-- rq-6a530d39 -->

The wet NPT stage couples pressure semi-isotropically with the xy plane held and z free.

A biochar sheet is a slab. Its lateral dimensions are set by the carbon lattice, which is not
compressible in any physically meaningful sense on this timescale, while the solvent above and
below it must be free to find its own density. Coupling xy would squeeze the lattice to satisfy a
barostat.

GROMACS reads the two `compressibility` values in the order `xy z`, so the held axis is the first
value. The condensation workflow uses the opposite assignment for the opposite reason — it models
dry sheets separated by a vacuum gap along z, and freezing z is what preserves that gap. Both are
deliberate; the axis order is what makes them look alike.

```gherkin
Feature: Couple pressure on the axes that should relax

  @rq-51c8f769
  Scenario: The wet stage holds the lateral lattice and lets the box relax along z
    Given the wet NPT stage
    When its pressure coupling is read
    Then coupling is semi-isotropic
    And the xy compressibility is zero while the z compressibility is not
```

## A Run Directory Records What It Was Made From <!-- rq-c039505d -->

A structure that reached `fallback` did not meet the composition targets that were asked for. It
is still worth simulating, and the sweep records it as such, but a run directory that looks
identical to one built from a `strict_pass` structure loses that distinction at exactly the point
where results are collected and compared.

The run directory states the label, the status, and the structure it was built from.

```gherkin
Feature: Keep provenance with the run

  @rq-9e4f2d8e
  Scenario: A structure that only passed as a fallback is marked as such in its run directory
    Given a manifest row whose status is fallback
    When its run directory is written
    Then the directory records that status

  @rq-da9561c4
  Scenario: A row that produced no structure files gets no run directory
    Given a manifest row whose status is neither strict_pass nor fallback
    When the manifest is processed
    Then no run directory is written for it
    And the returned result names the reason
```

## Include Files Come From the Manifest <!-- rq-212c17a8 -->

A topology that `#include`s a molecule definition cannot run without it. The sweep records the
path of every file it wrote, including the `.itp`, and setup resolves include files from those
recorded paths.

Deriving the include's name from the topology's name assumes the two differ only by suffix. That
holds for a single molecule and fails for a surface, whose topology is `<basename>.top` while its
include is `<basename>_sheet.itp`. The failure is silent: the include is simply not copied, and
the run directory looks complete until `grompp` cannot find the moleculetype.

```gherkin
Feature: Copy the include files the topology actually names

  @rq-2fdbb62b
  Scenario: An include file the topology needs is resolved from the manifest, not a filename guess
    Given a manifest row whose include file does not share its topology's stem
    When its run directory is written
    Then the include file is present in the run directory
```

## Cross-references <!-- rq-68fe810a -->

- The topology contents these runs consume — and the atom-name agreement that avoids one of the
  suppressed warning classes — are specified in `gromacs-export.md`.
- `anneal_spec_for_htt` and the Wood annealing protocol it implements belong to the condensation
  workflow; see `condensation.md` when that document exists.
- The `status` values a manifest row can carry are set by the sweep driver; see
  `parameter-sweep.md` when that document exists.
