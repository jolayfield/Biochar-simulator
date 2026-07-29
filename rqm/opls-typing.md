# Feature: OPLS-AA Atom Typing and Force-Field Mapping <!-- rq-c2c2dbbd -->

`AtomTyper` (`src/biochar/opls_typing.py`) assigns each atom an internal OPLS type, and
`GROMACS_OPLS_TYPE_MAP` (`src/biochar/constants.py`) translates those to the `opls_XXX` names that
GROMACS resolves against a real `oplsaa.ff`.

This module has the worst failure characteristics in the project. When the exported `.itp`
`#include`s a real forcefield, **only the `opls_XXX` name reaches GROMACS** — mass, charge,
Lennard-Jones, and every bonded parameter are looked up from that name. A name that identifies
the wrong element therefore simulates the wrong chemistry, produces plausible-looking output,
and nothing downstream complains. Every requirement below exists because some version of that
failure actually shipped.

## Typing Must Not Depend on Molecule Size <!-- rq-7709e5f5 -->

The single question "is this an aromatic carbon?" must be answered in exactly one place,
`_is_aromatic_carbon(mol, atom)`.

RDKit kekulization fails on large fused sheets — it emits `Can't kekulize mol` and leaves
`GetIsAromatic()` returning false on **every** atom. Any caller reading `GetIsAromatic()`
directly therefore gets an answer that changes with molecule size. The carbon branch carried a
ring-membership fallback for this; the nitrogen branch did not, so the same Ar-NH2 group typed
as `NA` on a 40-carbon flake and as `N` on a 150-carbon sheet.

The helper returns true when the atom is a carbon and either RDKit reports it aromatic, or it is
in a ring with degree 3 (two ring-carbon neighbours plus one H or O, or three ring carbons at an
interior junction).

```gherkin
Feature: Type atoms the same way regardless of molecule size

  @rq-3789a326
  Scenario: An amino nitrogen types as aniline N on a small flake
    Given a structure of 40 carbons carrying one amino group
    When atom types are assigned
    Then the nitrogen types as NA
    And its hydrogens type as HNA

  @rq-c7715c26
  Scenario: An amino nitrogen types as aniline N on a large sheet
    Given a structure of 150 carbons carrying one amino group
    When atom types are assigned
    Then the nitrogen types as NA
    And its hydrogens type as HNA

  @rq-4063a024
  Scenario: Kekulization failure does not change any atom type
    Given a structure large enough that RDKit cannot kekulize it
    When atom types are assigned
    Then every atom receives the type it would receive if kekulization had succeeded
```

## Every Assignable Type Must Be Mapped <!-- rq-bed2cb4c -->

`GromacsExporter` writes `GROMACS_OPLS_TYPE_MAP.get(t, t)`. An internal type with no mapping is
therefore written into the topology verbatim, and `grompp` fails with `Atomtype <name> not
found`. Nothing in the export path defines types locally, so there is no recovery.

Two invariants follow, and both are checked mechanically rather than by review:

- Every entry in `OPLS_ATOM_TYPES` has a `GROMACS_OPLS_TYPE_MAP` entry.
- Every entry in `GROMACS_OPLS_TYPE_MAP` corresponds to a real `OPLS_ATOM_TYPES` entry.

An internal type that the pipeline can no longer assign must be **removed**, not left mapped.
The tertiary `N` and quaternary `NT` types were unreachable-but-defined for exactly as long as
it took a kekulization failure to make them reachable again.

```gherkin
Feature: Keep the internal type table and the GROMACS map in agreement

  @rq-70e5f033
  Scenario: No assignable type lacks a GROMACS mapping
    Given the set of internal OPLS atom types
    When each is looked up in the GROMACS type map
    Then every one resolves

  @rq-3e169a31
  Scenario: The map holds no entries for types that do not exist
    Given the GROMACS type map
    When each key is looked up in the internal type table
    Then every one is present

  @rq-70f6b542
  Scenario: An unreachable nitrogen environment is flagged, not guessed
    Given a nitrogen the generator cannot produce, such as a trialkylamine
    When atom types are assigned
    Then it receives the X<Z> fallback type
    And validation reports it as an unrecognised type
```

## Three Depths of Force-Field Verification <!-- rq-794e1efe -->

Asserting that a mapped `opls_XXX` name exists is not sufficient. Verification has three depths,
each of which has caught a real bug that the previous depth missed:

1. **The name exists** in `oplsaa.ff/atomtypes.atp`.
2. **The element and mass are consistent** with the internal type's intent. `SS` once mapped to
   `opls_209`, which is a *carbon*; the name existed and resolved cleanly.
3. **Every bond and angle the topology emits resolves** through the bonded types into
   `ffbonded.itp`. `SS → opls_222` has the right element and the right mass, yet a thioether
   emits a `CA-S-CA` angle that stock OPLS does not define, and `grompp` refused the topology.

Depth 3 is the one that is easy to skip and expensive to skip.

```gherkin
Feature: Verify types against a real forcefield, not against a table

  @rq-03a40f15
  Scenario: Every mapped type names the right element
    Given each internal type and its mapped opls_XXX name
    When the forcefield atomtypes are consulted
    Then the element implied by the internal type matches the forcefield entry

  @rq-35a5407d
  Scenario: Every emitted angle resolves in the forcefield
    Given a structure exercising every functional group the generator can place
    When the topology is written
    Then every emitted angle resolves to an angletype in ffbonded.itp or the supplement

  @rq-b40b353d
  Scenario: Forcefield-backed tests fail loudly when the forcefield is absent
    Given no discoverable oplsaa.ff
    When the forcefield-backed tests run
    Then they skip with a message naming the missing forcefield
    And the skip is visible rather than reported as a pass
```

## The Supplementary Angle Table Is Narrowly Scoped <!-- rq-915e673b -->

Lennard-Jones, bond, and angle parameters live in `oplsaa.ff`, not in this repository. A
hand-copied table here is dead weight that can only drift — a previous one did, ending up a mix
of AMBER and OPLS values with no single provenance, and it was deleted.

`SUPPLEMENTARY_ANGLE_PARAMS` is the sole exception, and the rule that keeps it from becoming
that table again is strict: **it may hold only combinations stock `oplsaa.ff` does not define.**
A value that also exists in the forcefield is duplication and will drift; a value that exists
nowhere else cannot. Every entry carries a provenance comment naming its analog.

```gherkin
Feature: Supplement the forcefield only where it has a genuine gap

  @rq-78f96107
  Scenario: The supplement never shadows a stock entry
    Given each entry in the supplementary angle table
    When the stock forcefield is consulted
    Then the combination is absent from the stock forcefield

  @rq-8e4335b8
  Scenario: Every supplement entry states where its value came from
    Given each entry in the supplementary angle table
    Then it carries a provenance comment naming the analog it derives from
```

## Nearest-Analog Values Are Chosen, Not Validated <!-- rq-1729edcb -->

Several parameters are nearest-analog choices with provenance comments rather than
QM-validated values. This is a deliberate, recorded limitation, not an oversight:

- `CA-S-CA` — the thioether angle stock OPLS omits.
- `NPY-CA-OH` — hydroxyl on a ring carbon adjacent to pyridinic N, derived from the phenol
  angle `CA-CA-OH`; 3-hydroxypyridine is real chemistry OPLS simply does not cover.
- `SS → opls_222` — "S in thioanisoles", the nearest aryl-S.
- `NGR → opls_379` — "CytH+ N3", a protonated cationic aromatic ring N. **Chosen deliberately
  on 2026-07-17 and reviewed.** Do not "restore" `opls_520` on the assumption that this drifted.
- The phenolate, thiophenolate, pyridinium, and anilinium types — stock OPLS has no aryl oxide,
  aryl thiolate, or quaternary aromatic N.

These are chemically plausible and internally consistent. They are **not publication-grade**,
and a study whose conclusions turn on these specific parameters needs QM validation first.

```gherkin
Feature: Record analog-derived parameters as a known limitation

  @rq-1548650e
  Scenario: Analog-derived types remain stable across edits
    Given the analog-derived types NGR, SS, NPY-CA-OH and CA-S-CA
    When the constants tables are read
    Then each holds the value recorded in its provenance comment
```

## Carboxyl and Phenolic Oxygens Are Different Types <!-- rq-a7480504 -->

A carboxyl -OH and a phenolic -OH both present as "an oxygen with one hydrogen and one heavy
neighbour". Only the neighbour distinguishes them: a carboxyl carbon carries two oxygens, a
ring carbon bearing a phenol carries one. Typing them alike gives the carboxyl proton the
phenol's charge.

Both oxygens of a deprotonated carboxylate take the carboxylate type, because they are
equivalent by resonance — RDKit's `C(=O)[O-]` Kekulé form puts the formal charge on only one of
them, and typing the other from its double bond alone would leave one oxygen on the neutral acid
type.

```gherkin
Feature: Separate carboxyl oxygens from phenolic oxygens

  @rq-8669a35d
  Scenario: A carboxyl hydroxyl does not type as a phenol
    Given a structure carrying a carboxyl group
    When atom types are assigned
    Then the carboxyl hydroxyl oxygen and its hydrogen take the carboxylic-acid types

  @rq-b37f44d2
  Scenario: Both carboxylate oxygens take the anionic type
    Given a deprotonated carboxylate
    When atom types are assigned
    Then both of its oxygens type as carboxylate rather than one as a neutral carbonyl
```

## Cross-references <!-- rq-ef429ed9 -->

- The functional groups typed here are placed by `heteroatom-assignment.md`.
- The full reasoning behind the three verification depths is recorded in
  `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md`; do not re-derive it.
- The convention for holding a deferred gap open with `xfail(strict)` is in
  `docs/solutions/conventions/xfail-strict-as-deferred-gap-tripwire.md`.
