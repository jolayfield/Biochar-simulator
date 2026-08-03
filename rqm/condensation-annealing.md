# Feature: Condensation Annealing <!-- rq-881e08bd -->

`workflows/condensation.py` sets up the Wood et al. 2024 condensation protocol: island-type
building blocks are packed into a periodic box and condensed into an amorphous bulk solid by an
HTT-scaled heat/cool cycle, then expanded into an exposed surface.

**It is setup-only.** It writes `.mdp` files, a `.top`, and a shell script; it never invokes `gmx`.
Producing a finished model means running roughly 45 ns × 3 on the user's own build, and that
separation is deliberate — the script is meant to be read before it is run, and a module that ran
it would make reading it optional.

What this module is really writing down is somebody else's published protocol. The temperatures,
durations, barostats and repeat count all come from one paper, and their value is that they match
it. So the requirements here are mostly about **fidelity to a citation**: the anchors are the
anchors, the schedule is the schedule, and where the protocol is silent the choice is stated as a
choice rather than presented as theirs.

The second theme is that a script nobody watched must fail rather than drift. A packing step that
places fewer molecules than the topology declares, or a `grompp` allowed to absorb warnings nobody
named, turns a 45 ns run into a result that looks like an answer.

## Feature API <!-- rq-23755c47 -->

- `anneal_spec_for_htt(htt_c: float) -> AnnealSpec` <!-- rq-f32774a4 -->
  - Maps a pyrolysis heat-treatment temperature to the annealing peak temperature and timestep.

- `AnnealSpec` <!-- rq-55893f23 -->
  - The two HTT-dependent knobs: `peak_T_K` and `timestep_fs`. The final equilibration stage does
    not use the latter — it is always 2 fs.

- `render_mdp_set(spec) -> Dict[str, str]` <!-- rq-039cc5bf -->
  - The four stage `.mdp` files: energy minimisation, NVT, the NPT anneal, and the NPT final.

- `render_condensation_script(gro_name, top_name, spec, n_repeats=3, ...) -> str` <!-- rq-dd22e354 -->
  - The run script chaining the four stages across independent repeats.

- `setup_condensation(output_dir, molecule_gro, molecule_itp, n_copies, ...) -> Path` <!-- rq-bb992384 -->
  - Writes a complete run directory around one existing molecule: its files, a `.top` declaring
    `n_copies` of it, the `.mdp` set, and the script that packs and runs.

- `generate_and_condense(output_dir, n_copies, generator_config=None, ...) -> Path` <!-- rq-142d216f -->
  - The same, generating the molecule first.

- `setup_surface(...)`, `add_surface_and_validation(...)` <!-- rq-be66ed17 -->
  - The expansion of a condensed bulk into an exposed surface, and the pore-probe validation run
    beside it.

- `estimate_box_nm(gro_text, n_copies, loose_factor=1.6) -> float` <!-- rq-475641a0 -->
  - A packing box loose enough for insertion to succeed before condensation closes it.

- `CondensationError` <!-- rq-ee864495 -->
  - Raised for a setup that cannot be written as asked.

## The Published Anchors Are Reproduced Exactly <!-- rq-7dee31f6 -->

Wood anchored three heat-treatment temperatures: 400 °C to 1000 K, 600 °C to 2000 K, 800 °C to
3000 K. Those three pairs are the citation, and a request at one of them returns it unchanged.

Between them the peak temperature is interpolated linearly, and outside them it is clamped. Both
are extensions past what the paper states, which is why the anchors themselves must be exact: a
reader checking this module against Table 6 is checking the three rows that are actually there.

```gherkin
Feature: Return the published values at the published points

  @rq-9ba33460
  Scenario: Each anchored heat-treatment temperature returns its published pair
    Given one of the three heat-treatment temperatures Wood anchored
    When the annealing spec is derived
    Then the peak temperature and timestep are the published ones

  @rq-50d5f27c
  Scenario: A temperature between anchors is interpolated
    Given a heat-treatment temperature between two anchors
    When the annealing spec is derived
    Then the peak temperature lies between the two anchored peaks

  @rq-df4cf342
  Scenario: A temperature outside the anchors is clamped
    Given a heat-treatment temperature below the lowest anchor or above the highest
    When the annealing spec is derived
    Then the peak temperature is the nearest anchored peak
```

## The Timestep Is Chosen by the Anchor, and Says So <!-- rq-c1dbabdb -->

Wood used 1.0 fs at 400 °C and 0.5 fs at both 600 °C and 800 °C. Only three points are given, so
any rule for the temperatures in between is this module's, not theirs.

The rule is the conservative one: 1.0 fs at or below the lowest anchor, 0.5 fs above it. A request
at 450 °C therefore anneals at roughly 1250 K with a 0.5 fs step — twice the wall-clock of a 1.0 fs
step, chosen because the halving point is unknown and an unstable high-temperature run is worse
than a slow one.

That reasoning has to be the one recorded. Describing the rule as a temperature threshold when it
is keyed to the anchor sends a reader looking for a peak-temperature comparison that is not in the
code, and invites a later change that "restores" a threshold the module never had.

```gherkin
Feature: Say which rule the timestep follows

  @rq-0a4c0f30
  Scenario: The timestep halves above the lowest anchor
    Given a heat-treatment temperature just above the lowest anchor
    When the annealing spec is derived
    Then the timestep is the smaller of the two published values

  @rq-fa9763e9
  Scenario: The recorded rationale is the rule the code follows
    Given the documented explanation of the timestep choice
    When it is compared with the condition the code tests
    Then they describe the same rule
```

## The Annealing Schedule Is the Published One <!-- rq-7b79f60e -->

Twenty-five nanoseconds in three parts: hold the peak for ten, cool linearly to 300 K over the next
ten, hold 300 K for the last five. The NVT stage before it is ten nanoseconds at the peak, and the
final equilibration is ten more at 300 K and 1 bar with a 2 fs step.

Step counts follow from the durations and the stage's timestep, so a change to either is visible in
the `.mdp` rather than needing to be applied twice.

```gherkin
Feature: Render the protocol's own durations

  @rq-3d4eeaf7
  Scenario: The annealing schedule holds, cools, then holds
    Given a rendered annealing stage
    When its temperature schedule is read
    Then it holds the peak, cools to 300 K, and holds 300 K over the published intervals

  @rq-a16e0b49
  Scenario: Step counts follow from the duration and the timestep
    Given a rendered stage
    When its step count is compared with its duration and timestep
    Then they agree

  @rq-3a97bd4a
  Scenario: The final equilibration is always two femtoseconds
    Given a rendered final stage at any heat-treatment temperature
    When its timestep is read
    Then it is two femtoseconds
```

## Every Repeat Runs the Full Chain <!-- rq-53200bd1 -->

Wood ran three repeats from different starting configurations for ergodicity. Each repeat packs
with its own insertion seed, then runs energy minimisation, NVT, the anneal and the final
equilibration in that order, into its own directory.

```gherkin
Feature: Run each repeat independently and in order

  @rq-5084d61d
  Scenario: Each repeat runs the four stages in order
    Given a rendered run script
    When its stages are read
    Then each repeat runs minimisation, NVT, anneal and final in that order

  @rq-516fd1c6
  Scenario: Repeats are seeded differently
    Given a rendered run script that packs
    When the insertion seeds are read
    Then each repeat uses its own
```

## The Packed Box Contains What the Topology Declares <!-- rq-e151b7e1 -->

`gmx insert-molecules` is a placement search. Asked for two hundred copies in a box too tight for
them, it places what it can, reports the shortfall on standard output, and exits successfully.

The `.top` beside it declares the number that was asked for. Nothing reconciles the two, so the
first thing that notices is `grompp`, several commands later, complaining about a coordinate count
that does not match the topology — a message that describes the symptom and not the cause, at the
point furthest from it.

The packing step checks its own result against the number the topology declares and stops there if
they differ, naming both.

```gherkin
Feature: Stop where the packing shortfall happens

  @rq-9a027a37
  Scenario: The script verifies the packed count against the topology
    Given a rendered run script that packs
    When the packing stage is read
    Then it compares the molecules placed with the number the topology declares
    And stops with both numbers if they differ
```

## A Suppressed grompp Warning Is Named <!-- rq-d32ee60c -->

`grompp` warnings are refusals to proceed silently, and `-maxwarn` takes a count rather than a
name: a stage allowed two warnings accepts the two the author had in mind, or any other two.

Every stage in this script passes `-maxwarn 2`. Nothing in the four stages is expected to warn at
all — the topology is the one the exporter wrote, the coordinates are the ones packing produced —
so the allowance excuses nothing that was foreseen and absorbs anything that turns up.

The same rule applies here as to the run directories `md-setup.md` specifies: `-maxwarn` appears
only where a specific warning is expected, with a comment naming it and a count equal to the number
named. A stage with nothing to excuse passes none and fails on the first surprise.

```gherkin
Feature: Suppress only the warnings that were understood

  @rq-30534b7a
  Scenario: A stage with no expected warning suppresses nothing
    Given a rendered condensation script
    When its grompp stages are read
    Then a stage with no named expected warning passes no -maxwarn

  @rq-8410a38d
  Scenario: A suppressed warning is named where it is allowed
    Given a rendered condensation script
    When a grompp stage passes -maxwarn
    Then a comment on that stage names the warning it absorbs
    And the count equals the number of warnings named
```

## A Run Directory Records What Set It Up <!-- rq-c0d9fdf2 -->

These directories outlive the call that produced them by months, and the thing a reader needs from
one is what it is a simulation *of*: which heat-treatment temperature, which peak and timestep the
mapping produced, how many copies of what molecule, and which version of this package decided.

None of it is recoverable from the files. The `.mdp` files carry the peak temperature but not the
HTT it came from, the `.top` carries a molecule count but not the box it was packed into, and
neither says whether the mapping that produced them has since changed.

```gherkin
Feature: Record the setup beside the run

  @rq-75419e1a
  Scenario: A condensation run directory records its provenance
    Given a condensation run directory
    When its provenance record is read
    Then it names the heat-treatment temperature, the annealing spec, the copy count and the package version
```

## The Module Writes Files and Nothing Else <!-- rq-05ee02a8 -->

Setting up a run and running it are different acts with different costs, and this module only does
the first. A script that would be reviewed before running is worth writing; a module that ran it
would make the review optional, and the review is where a 45 ns mistake is cheap to catch.

```gherkin
Feature: Set up without running

  @rq-ec21c31b
  Scenario: Setting up a condensation run invokes no external program
    Given a request to set up a condensation run
    When the setup completes
    Then no external process was started
```

## Cross-references <!-- rq-94d2e77f -->

- The `-maxwarn` rule this module shares, and the run provenance the MD run directories record, are
  specified in `md-setup.md`. The annealing schedule in those directories is rendered from this
  module, so the two cannot drift.
- The single-molecule `.gro` and `.itp` that packing consumes are specified in `gromacs-export.md`.
- The generator configuration `generate_and_condense` builds its molecule from is specified in
  `generation-config.md`.
