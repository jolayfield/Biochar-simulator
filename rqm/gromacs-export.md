# Feature: GROMACS Topology and Structure Export <!-- rq-14ed3cdc -->

`src/biochar/export/gromacs_export.py` turns a typed, charged molecule into the three files
GROMACS reads: a `.gro` of coordinates, an `.itp` defining the molecule, and a `.top` that
includes a real forcefield and names the system.

The export path inherits the failure characteristic that `opls-typing.md` describes and extends
it one file further. The `.itp` carries **no parameters** — no masses beyond the atom line, no
Lennard-Jones, no force constants. Every physical quantity is resolved by `grompp` from the
`opls_XXX` names and the bonded terms this module emits. A term that is wrong, or a term that is
simply absent, produces a topology that `grompp` accepts and a simulation that runs and reports
plausible numbers while integrating the wrong physics.

Absence is the harder half. A wrong name eventually collides with something; an omitted
interaction never announces itself. The requirements here therefore pin what the topology
**must contain** at least as tightly as what it must get right.

## Feature API <!-- rq-0e187d3e -->

- `GROFileWriter.write(filepath: str, mol: Chem.Mol, coords: np.ndarray, molecule_name: str = "BC", box_vectors: Optional[np.ndarray] = None, title: Optional[str] = None) -> None` <!-- rq-312047ff -->
  - Writes fixed-width `.gro`: title, atom count, one line per atom, box vectors last.
  - Converts coordinates from ångström to nanometres. Written to three decimal places, so
    positions carry 0.001 nm of quantisation.
  - Truncates the residue name to five characters, the `.gro` format limit, and warns when it
    does.
  - Emits three box values for an orthogonal cell and nine in GROMACS triclinic order otherwise.
  - Wraps residue and atom numbers modulo 100000, the width of their fields.

- `ITPFileWriter.write(filepath: str, mol: Chem.Mol, atom_types: Dict[int, str], charges: Dict[int, float], molecule_name: str = "BIOCHAR", include_dihedrals: bool = True) -> None` <!-- rq-8340860d -->
  - Writes `[ moleculetype ]`, `[ atoms ]`, `[ bonds ]`, `[ angles ]`, and proper `[ dihedrals ]`.
  - Emits no `[ atomtypes ]`. The including `.top` resolves types from the forcefield, and a
    local definition would shadow it.
  - Emits bonded terms as connectivity plus a function number, with no parameters, so every
    constant comes from `ffbonded.itp`.
  - States the molecule's running and total charge as comments after `[ atoms ]`.

- `TOPFileWriter.write(filepath: str, mol: Chem.Mol, atom_types: Dict[int, str], charges: Dict[int, float], molecule_name: str = "BIOCHAR", forcefield_path: str = "oplsaa.ff/forcefield.itp", include_dihedrals: bool = True) -> None` <!-- rq-330a61eb -->
  - Writes a self-contained topology: the forcefield `#include`, the molecule definition,
    `[ system ]`, and `[ molecules ]`.
  - Names a non-zero total charge and the step that neutralises it.

- `GromacsExporter(output_directory: str = ".")` <!-- rq-e17f705d -->
  - `export(mol, coords, atom_types, charges, molecule_name="BC", basename="structure", include_periodic_box=False, box_size=None) -> Tuple[Path, Path, Path]`
    returns the `.gro`, `.top`, and `.itp` paths in that order.
  - Sizes a box around the molecule when one is not supplied, and translates coordinates into
    it. Exported positions are therefore related to the input by a scale and a rigid shift, not
    equal to them.

## Coordinates Convert to Nanometres <!-- rq-b0c8554c -->

RDKit works in ångström and the `.gro` format is defined in nanometres. Coordinates are divided
by ten on the way out.

The consequence of getting this wrong is a box ten times too large in every dimension, a density
a thousand times too low, and a trajectory that runs to completion looking like a dilute gas.
Nothing in the toolchain flags it, because both magnitudes are dimensionally plausible.

A test for this must compare a quantity that is invariant under the export's rigid shift.
Absolute positions are not, since box placement translates the molecule; interatomic distances
are. Distances also make the check unambiguous — a factor of ten cannot hide inside the
0.001 nm write precision.

```gherkin
Feature: Write coordinates in the units the format defines

  @rq-6024903d
  Scenario: Every interatomic distance scales by exactly one tenth
    Given a molecule whose coordinates are in Angstrom
    When the structure file is written
    Then every interatomic distance in the file is one tenth of the Angstrom distance
    And the agreement is within the file's three-decimal write precision
```

## Atom Names Identify One Atom Each <!-- rq-c817effd -->

An atom's name appears in both the `.gro` and the `.itp`, and `grompp` compares them. The two
files agree on every atom, and within a molecule no two atoms share a name.

Names are built from the element symbol and the atom index, so their length grows with the
molecule. The `.gro` atom-name field is five characters wide. A molecule large enough to need a
sixth character is exactly where two independent things break at once: the name in the `.gro`
stops matching the name in the `.itp`, and distinct atoms collapse onto the same `.gro` name —
the ten-thousandth carbon and the thousandth carbon both become `C1000`.

`grompp` reports the mismatch as a warning rather than an error. That warning is one of the
suppressed classes described in `md-setup.md`, so it reaches nobody.

```gherkin
Feature: Keep atom names consistent and unique across the exported files

  @rq-07b27c79
  Scenario: Atom names agree between the structure and the topology at any size
    Given a molecule large enough that its atom names exceed five characters
    When the structure and topology are written
    Then the atom name for every atom is the same in both files

  @rq-b2e50e9d
  Scenario: No two atoms share a name in the structure file
    Given a molecule of more than ten thousand atoms
    When the structure file is written
    Then no atom name appears twice
```

## The Topology Carries Every Non-Bonded Interaction It Excludes <!-- rq-25d4feed -->

`[ moleculetype ]` declares `nrexcl = 3`, which is what OPLS-AA expects: non-bonded interactions
between atoms separated by three bonds or fewer are excluded from the plain non-bonded loop.

OPLS-AA does not discard those 1–4 interactions. It reintroduces them at reduced strength — the
forcefield's `[ defaults ]` sets `gen-pairs = yes` with `fudgeLJ = 0.5` and `fudgeQQ = 0.5`, and
those scaled parameters are applied to the pairs listed in `[ pairs ]`. `gen-pairs` generates the
pair *parameters*; it does not generate the pair *list*. A topology that declares `nrexcl = 3`
and omits `[ pairs ]` therefore excludes every 1–4 interaction and restores none of them, which
is a different force field from the one it claims to use.

```gherkin
Feature: Exclude and restore 1-4 interactions the way OPLS-AA defines them

  @rq-a22b8e8d
  Scenario: The exclusion count matches what OPLS-AA expects
    Given a molecule exported to a topology
    When the moleculetype line is read
    Then the exclusion count is 3

  @rq-fb164878
  Scenario: Every excluded 1-4 pair is listed for rescaled treatment
    Given a molecule whose bonded graph contains at least one 1-4 pair
    When the topology is written
    Then it contains a pairs section
    And every 1-4 pair in the molecule appears in it
```

## Aromatic Rings Carry the Terms That Keep Them Planar <!-- rq-6659180d -->

An aromatic ring is planar because its sp2 centres are held planar, and under OPLS-AA that is
the job of improper dihedrals. Proper dihedrals constrain rotation about a bond; they do not stop
a substituted ring carbon from pyramidalising out of the ring plane.

A topology whose only torsional terms are proper dihedrals is free to buckle at every sp2 centre
during dynamics. The structure this package generates is planar when written and has no term
requiring it to stay that way. The stock forcefield's own residue templates carry an
`[ impropers ]` block for exactly this reason.

```gherkin
Feature: Constrain aromatic centres out of plane

  @rq-2f3c2379
  Scenario: A ring carbon carries an improper term
    Given a molecule containing an aromatic ring
    When the topology is written
    Then it contains an improper dihedral for each aromatic ring carbon
```

## A Charged Molecule Says So <!-- rq-190f0328 -->

Deprotonation at a requested pH leaves a molecule with a non-zero formal charge, and a charged
system under particle-mesh Ewald needs a neutralising counter-ion before it will simulate
correctly.

The topology states the total charge it carries and names the step that balances it, rather than
leaving the reader to sum the `[ atoms ]` column. The partial charges sum to the formal charge:
the molecule's charge is a property of its protonation state, not an artefact of rounding the
charge assignment.

```gherkin
Feature: Report a non-zero system charge rather than leaving it implicit

  @rq-78e7477d
  Scenario: The topology states a non-zero total charge and how to neutralise it
    Given a molecule carrying a non-zero formal charge
    When the topology is written
    Then it states the total charge
    And it names the neutralisation step

  @rq-f256712c
  Scenario: Partial charges sum to the formal charge
    Given a molecule carrying a non-zero formal charge
    When the topology is written
    Then the sum of the partial charges equals the formal charge
```

## The Residue Name Fits the Format <!-- rq-86ce8bbe -->

The `.gro` residue-name field is five characters. A longer `molecule_name` is truncated to its
first five characters and the caller is told, rather than the file being written with a shifted
column layout that GROMACS would misparse.

```gherkin
Feature: Keep the residue name within the format's field

  @rq-5340daab
  Scenario: A residue name longer than five characters is truncated
    Given a molecule_name of more than five characters
    When the structure file is written
    Then the residue name in the file is its first five characters
    And the caller is warned that truncation happened
```

## Cross-references <!-- rq-66a857ff -->

- Atom typing, the `opls_XXX` mapping, and the three depths of forcefield verification —
  including that every emitted bond, angle, and dihedral resolves in `ffbonded.itp` — are
  specified in `opls-typing.md`.
- The `grompp` warning classes this module's output can trip, and their suppression, are
  specified in `md-setup.md`.
- Multi-sheet export and the relationship between the `.itp` count and the `[ molecules ]`
  count are specified in `surface-stacking.md`.
