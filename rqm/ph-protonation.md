# Feature: pH-Dependent Protonation <!-- rq-67aceac9 -->

`pipeline/protonation.py` decides, per ionizable site, whether that site is protonated or ionized
at a given pH, and applies the decision to the molecule. It runs after every heteroatom exists, so
there are sites to titrate, and before hydrogens are filled in, which owns acidic-H placement and
would otherwise fight the decision.

The states are **discrete and assigned once, at build time**, sampled from Henderson-Hasselbalch.
This is not constant-pH MD: there is no lambda dynamics and no state change during a run. That
single fact drives the requirements here. A generated structure is not *the* structure at that pH,
it is **one sample from the ensemble** at that pH, and every promise below exists to keep that
visible — in the census the generator reports, in the warning a caller gets when the sample is
close to a coin flip, and in the files the structure is written to.

The second organising fact is that the charge has a sign, and the sign depends on which side of
the pair the proton sits. An acidic group is neutral when protonated and becomes an anion by losing
H⁺, so its ionized fraction **rises** with pH. A basic group is neutral when deprotonated and
becomes a cation by gaining H⁺, so its ionized fraction **falls** with pH. Collapsing those two
cases is the easiest way to get this feature silently wrong: it inverts the charge on every
nitrogen while still producing plausible output that runs.

## Feature API <!-- rq-2fcb3d7a -->

- `fraction_ionized(pKa: float, pH: float, kind: str) -> float` <!-- rq-f1bddaf4 -->
  - Probability in `[0, 1]` that one site is in the ionized form. `pKa` always belongs to the
    protonated member of the pair; `kind` is `"acidic"` or `"basic"` and decides the direction.

- `ProtonationAssigner(seed: Optional[int] = None)` <!-- rq-9e66167d -->
  - `assign(mol, pH, result=None) -> Tuple[Chem.Mol, CompositionResult]` draws every site
    independently against its Henderson-Hasselbalch probability and applies the outcome.
  - Uses an instance-local RNG, so a caller's own random stream is untouched.

- `PROTONATION_STATES` <!-- rq-1713a8a2 -->
  - The titratable groups: pKa, kind, the OPLS type of the neutral form, of the ionized form, and
    of the exchangeable hydrogen. A group absent from this table is never titrated.

- `PH_MIN`, `PH_MAX` <!-- rq-2bb117a3 -->
  - The aqueous range a request must fall in.

- `CompositionResult.net_charge`, `.ionized_counts`, `.titratable_counts` <!-- rq-b5909c8c -->
  - The census: the formal charge the molecule ended up carrying, what actually ionized, and what
    could have.

## The Direction of Ionization Follows the Kind of Group <!-- rq-33ac898e -->

`fraction_ionized` is the whole model, and its one subtlety is the sign. For both kinds the
fraction deprotonated is `1 / (1 + 10 ** (pKa - pH))`; what differs is whether *deprotonated* means
*ionized*.

For an acid it does, so raising the pH makes the structure more anionic. For a base it is the
reverse: the ionized form is the protonated one, so raising the pH makes the structure *less*
cationic. A molecule carrying both is cationic at low pH, anionic at high pH, and passes through
its isoelectric point in between.

```gherkin
Feature: Titrate acids and bases in opposite directions

  @rq-952d1150
  Scenario: An acidic group ionizes as pH rises
    Given an acidic group
    When the pH is raised past its pKa
    Then the fraction ionized increases

  @rq-9bf47e81
  Scenario: A basic group deionizes as pH rises
    Given a basic group
    When the pH is raised past its pKa
    Then the fraction ionized decreases

  @rq-5c5a5c84
  Scenario: A molecule carrying both is cationic at low pH and anionic at high
    Given a molecule with both an acidic and a basic group
    When it is titrated at a low pH and again at a high pH
    Then its net charge is positive at the low pH and negative at the high pH
```

## A Structure Is One Sample, Not the Ensemble Mean <!-- rq-072134ee -->

Each site is drawn independently against its own probability rather than flipped on a `pH > pKa`
threshold. A threshold would ionize every carboxyl on every flake at the same pH, producing a step
function where the chemistry is a curve, and no site-to-site variation within one sheet.

The cost of sampling is that a single structure is one draw. Away from a pKa this barely matters —
at pH 7 a carboxyl is ionized with probability 0.998 — but within about a unit of the pKa each site
is close to a coin flip, and one structure is then a poor stand-in for the ensemble. Averaging
across seeds recovers the Henderson-Hasselbalch fraction; a single structure does not.

So the sampling is made visible rather than left to be inferred. The census reports what was
actually placed, the seed reproduces a sample exactly, and a request that lands in a transition
band says so.

```gherkin
Feature: Keep the sampling visible

  @rq-e9369729
  Scenario: The same seed and pH reproduce the same sample
    Given a molecule titrated at a pH with a given seed
    When it is titrated again with the same seed and pH
    Then the same sites are ionized

  @rq-7286d3d9
  Scenario: Titration does not disturb the caller's random stream
    Given a caller that has seeded the global random module
    When a molecule is titrated
    Then the caller's next random draw is unchanged

  @rq-7e778b79
  Scenario: A pH inside a group's transition band is reported
    Given a molecule carrying a group whose pKa is within a unit of the requested pH
    When it is titrated at that pH
    Then a warning names the group and says one structure is a sample of the ensemble

  @rq-9ee7e1df
  Scenario: A pH far from every pKa is not reported
    Given a molecule whose every titratable group has a pKa more than a unit away
    When it is titrated at that pH
    Then no transition-band warning is issued
```

## The Census Counts What Happened, Not What Was Decided <!-- rq-191c700f -->

A decision that cannot be carried out contributes no charge — a site with no acidic hydrogen to
give up stays neutral. Counting the decision rather than the outcome would put `ionized_counts` and
`net_charge` into disagreement, and the census is the only place a caller can see what the sampling
produced.

`titratable_counts` records what could have ionized. The difference between the two is the
information a caller needs to tell an unlikely draw from a site that was never able to titrate.

```gherkin
Feature: Report the protonation the molecule actually carries

  @rq-633ccda5
  Scenario: The census never claims an ionization the molecule does not carry
    Given a titrated molecule
    When its ionized counts and net formal charge are compared
    Then the counts do not imply more charge than the molecule carries
```

## Only Groups in the Table Are Titrated <!-- rq-99f30fb0 -->

Site detection asks the atom typer what each atom is and looks the answer up in
`PROTONATION_STATES`, so "what is a carboxyl" is answered in one place rather than re-derived here.
An ether oxygen, a thioether sulfur and a graphitic nitrogen have no entry, carry no acidic proton,
and are never touched at any pH.

That lookup inverts the table, mapping neutral OPLS type back to group. The inversion is only sound
while those types are unique: two groups sharing one would let the dictionary silently shadow one
of them, leaving it permanently untitratable at every pH with nothing to indicate it. The table is
checked when the assigner is constructed rather than trusted.

```gherkin
Feature: Titrate exactly the groups the table names

  @rq-4897529c
  Scenario: A group with no entry in the table is never titrated
    Given a molecule whose only heteroatom is an ether oxygen
    When it is titrated at an extreme pH
    Then it carries no charge

  @rq-09a84056
  Scenario: Two groups sharing a neutral type are refused
    Given a protonation table in which two groups share a neutral OPLS type
    When an assigner is constructed over it
    Then the table is refused rather than one group silently shadowing the other
```

## A Request Outside the Aqueous Range Is Refused <!-- rq-e07d2059 -->

Henderson-Hasselbalch saturates well inside pH 0–14: at three units past a pKa the fraction is
already within a thousandth of its limit, so nothing about the chemistry changes beyond the bounds.
A value outside them is far more likely a unit error than a request, and is refused rather than
silently clamped.

```gherkin
Feature: Refuse a pH that is not a pH

  @rq-af2d7d20
  Scenario: A pH outside the aqueous range is refused
    Given a pH below zero or above fourteen
    When a molecule is titrated at it
    Then the request is refused
```

## Applying a Decision Leaves the Rest of the Molecule Intact <!-- rq-bba3fc5f -->

Ionization edits the molecule — removing an acidic hydrogen, adding one to a base — and every edit
is followed by a sanitisation pass. Sanitisation re-perceives aromaticity, and a bridging ether
C–O closing a ring that satisfies Hückel is marked `AROMATIC` when it is not. Left uncorrected the
bond reaches typing as an aromatic C–O and takes a type that does not exist for it.

Index order matters for the same reason. Added atoms are appended and disturb nothing; removed
atoms shift every later index, so removals are deferred and applied in descending order rather than
inline, which would invalidate the index of every site not yet visited.

```gherkin
Feature: Leave untitrated chemistry unchanged

  @rq-7e9935a1
  Scenario: An ether bond is not left aromatic by the protonation pass
    Given a molecule with both an ether bridge and a titratable group
    When it is titrated
    Then the ether C-O bonds are still single bonds
```

## The Files Record the pH They Were Titrated At <!-- rq-479efbdd -->

Two structures generated at pH 3 and pH 11 differ in their protonation and in nothing else a reader
can see. They have the same skeleton, the same functional groups, and file headers that differ only
by timestamp.

The pH is not recoverable from the files either. A charge can be summed from the topology, but a
charge of zero is equally consistent with pH 2 and with no pH having been requested at all, and a
non-zero charge does not say which pH produced it or which of the ensemble's samples this is. The
structure is one draw from a distribution, and the draw is only interpretable alongside the pH and
the seed that produced it.

So the structure files state it: the pH, the net charge, and the seed the sample came from, written
where a reader opening the file will see them.

```gherkin
Feature: Say what pH a structure was built at

  @rq-7ed1460a
  Scenario: A titrated structure names its pH in the coordinate file
    Given a structure generated at a requested pH
    When its coordinate file is read
    Then the title states the pH it was titrated at and the charge it carries

  @rq-4189fd49
  Scenario: A titrated structure names its pH and seed in the topology
    Given a structure generated at a requested pH
    When its include topology is read
    Then the header states the pH and the seed the protonation was sampled with

  @rq-e0932c52
  Scenario: An untitrated structure claims no pH
    Given a structure generated with no pH requested
    When its coordinate file and topology are read
    Then neither claims a pH
```

## Cross-references <!-- rq-f72df2c8 -->

- The partial charges an ionized structure carries, the total-charge statement its topology makes,
  and the neutralising step that statement names are specified in `gromacs-export.md`.
- The refusal of a pH request combined with the ML charge backend is specified in
  `generation-config.md`.
- Per-sheet titration of a stacked surface, and the seed range it consumes, are specified in
  `surface-stacking.md`.
- The OPLS types the neutral and ionized forms take are specified in `opls-typing.md`.
