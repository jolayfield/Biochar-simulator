# Feature: Valence Validation <!-- rq-a6d2743c -->

`pipeline/valence.py` answers one question — how many bonds may this atom have, and how many does
it have — and the rest of the pipeline is built on the answer. `OxygenAssigner` picks its sites by
asking which carbons have free valence, `HydrogenAssigner` fills what is still needed, and
`StructureValidator` reports what is left over. A wrong answer here does not produce a valence
error; it produces a structure built around one.

The module is also a documented public surface (`docs/guides/VALENCE_SYSTEM.md`), so its helpers
are used by callers who are not this pipeline and who have no other check behind them.

Two things make the question harder than the table of standard valences suggests.

**Aromatic bonds have no integer order.** Summing bond orders gives 1.5 per aromatic bond, which is
right for a ring member that accepts into the π system — an aromatic carbon, a pyridinic nitrogen —
and wrong for one that donates a lone pair into it. A furan oxygen has two σ bonds and a valence of
two, but its bond orders sum to three.

**A formal charge changes the range, in both directions.** An anion trades a bond for a lone pair,
so it has one fewer bond *and* no room for another. Getting only the first half right leaves a
deprotonated phenolate looking like it needs a hydrogen, and the hydrogen assigner then puts one
back.

## Feature API <!-- rq-5e6bf91f -->

- `get_valence_range(atomic_num: int, formal_charge: int = 0) -> Tuple[int, int]` <!-- rq-9a268ae7 -->
  - The (minimum, maximum) bond count for an atom of this element and charge.

- `ValenceInfo` <!-- rq-de2d882d -->
  - One atom's answer: the range, the bonds it currently has, how many more it may take, how many
    it still needs, and whether it is currently valid.

- `ValenceValidator.get_valence_info(mol, atom_idx) -> ValenceInfo` <!-- rq-a802935e -->
  - The count for one atom. `available_valence` is what site-selection asks; `needed_valence` is
    what hydrogen saturation asks.

- `ValenceValidator.validate_molecule(mol) -> Tuple[bool, List[str]]` <!-- rq-dfa9c5ea -->
  - Every atom checked, with one message per violation naming the atom, the count and the bound.

- `SafeBondAdder.can_add_bond(mol, i, j, bond_type=1) -> Tuple[bool, str]` <!-- rq-8156f09f -->
  - Whether a bond of this order fits in what both atoms have left, and whether they are already
    bonded.

- `SafeBondAdder.add_bond_safe(emol, mol, i, j, bond_type=SINGLE) -> Tuple[bool, str]` <!-- rq-a05b950d -->
  - Adds the bond if it fits. Judged against the molecule as the editor currently has it, so a
    sequence of additions to one atom is bounded by that atom's valence rather than each addition
    being weighed against the state before any of them.

- `SafeBondAdder.add_atom_safe(emol, atomic_num, formal_charge=0) -> int` <!-- rq-48295a76 -->
  - Appends an atom and returns its index.

- `ValenceReport.get_summary(mol) -> Dict` <!-- rq-6b827a77 -->
  - Atom counts, per-element counts, and whether everything is valid.

## An Aromatic Ring Member Is Counted by Its Bonds, Not Its Bond Orders <!-- rq-3fcb12da -->

Aromatic bond orders are 1.5, and a heteroatom that donates its lone pair into the ring holds two
of them while having a valence of two. Summing orders reports three, so a furan oxygen, a thiophene
sulfur and a pyrrolic nitrogen are all reported as exceeding their maximum — three of the ring
chemistries this package exists to build.

The consequence is not a miscounted number in a report. `validate_molecule` feeds
`StructureValidator`, so under strict validation a structure containing a furan bridge is refused
for a violation that is not there.

An aromatic atom's valence is the count RDKit's own aromaticity model derives, which distinguishes
the donor from the acceptor. Everything else is counted by bond order as before.

```gherkin
Feature: Count an aromatic heteroatom correctly

  @rq-c70dd417
  Scenario: A furan-type ring oxygen is valid
    Given a molecule with an oxygen as an aromatic ring member
    When its valence is validated
    Then no violation is reported for that oxygen

  @rq-302b177e
  Scenario: A pyrrolic ring nitrogen carrying a hydrogen is valid
    Given a molecule with a pyrrole-type nitrogen
    When its valence is validated
    Then no violation is reported for that nitrogen

  @rq-e9202195
  Scenario: A thiophene-type ring sulfur is valid
    Given a molecule with a sulfur as an aromatic ring member
    When its valence is validated
    Then no violation is reported for that sulfur

  @rq-c179a6e7
  Scenario: A genuinely over-bonded atom is still reported
    Given a neutral nitrogen carrying four bonds
    When its valence is validated
    Then a violation is reported naming the atom and the bound it exceeds
```

## An Anion Has One Fewer Bond and No Room for Another <!-- rq-9cac4dad -->

A negatively charged heteroatom trades a bond for a lone pair. Both halves of that matter, and they
fail in opposite directions.

If the *minimum* is not lowered, a phenolate oxygen with one bond looks like it still needs one, and
hydrogen saturation puts a hydrogen back — silently undoing the deprotonation that the requested pH
asked for. If the *maximum* is not lowered, a complete anion looks like it has room, and a
bond-adding caller attaches to it.

```gherkin
Feature: Give an anion the range it actually has

  @rq-eb28d0a1
  Scenario: A deprotonated oxygen needs no further bond
    Given an oxygen carrying one bond and a negative formal charge
    When its valence is read
    Then it needs no further bonds

  @rq-a467c647
  Scenario: A deprotonated oxygen has no free valence
    Given an oxygen carrying one bond and a negative formal charge
    When its valence is read
    Then it has no valence available
```

## A Substituted Ring Carbon Is Not a Pyrrolic Nitrogen Site <!-- rq-53ec84de -->

Nitrogen doping replaces a ring carbon, and a pyrrolic nitrogen is defined by what it carries: two
σ bonds to its ring neighbours and one N–H, three bonds in total.

That leaves no room for anything the carbon was already carrying. Substituting a five-ring carbon
that already holds a functional-group oxygen produces a nitrogen bonded to two ring neighbours, that
oxygen, and the N–H added afterwards — four bonds on a neutral nitrogen, which is not a molecule.
Heteroatom placement runs before nitrogen substitution, so decorated carbons are exactly what the
substituter is choosing from.

The candidate must therefore be an unsubstituted five-ring carbon, in the same way that a graphitic
site must be an interior junction with no hydrogen — and only one per ring, since two nitrogens in
one pentagon is a pyrazole or an imidazole, where only one of them carries the hydrogen. A request
that cannot find enough qualifying carbons places fewer and says so, rather than taking a carbon
that does not qualify.

The requirement is about the nitrogen in the finished structure, not only about the site it was
placed on, and the remaining way to miss it is the ring rather than the carbon. A skeleton is not
aromatic everywhere: growth and aliphatic decoration leave pockets aromaticity perception refuses,
and a pentagon in one of them is a cyclopentadiene whose bonds are the kekulé single and double the
perception fell back to. Every free carbon in such a ring holds one of those double bonds, so a
substitution that edits only the atom gives a nitrogen with a C=N, and the N–H makes four again.

Refusing those rings would refuse the chemistry. A cyclopenta-fused pentagon is not aromatic; the
same ring with an N–H in it is a pyrrole, which is, and the nitrogen's own lone pair is what
aromatises it. The substitution therefore puts the ring into the state the new nitrogen implies —
the pentagon's bonds become aromatic, the double bond the site was holding included — and the
carbon that gave up that double bond is left an aromatic carbon with a free valence for its
hydrogen.

That is checked rather than assumed, because a pentagon can also hold an sp3 carbon from aliphatic
decoration, which no amount of flagging makes aromatic. The substitution is made on a copy and the
whole ring is read back; a ring that comes out with an atom over its maximum is discarded untouched
and the search moves on. Over the maximum and only that: hydrogen saturation has not run at this
point, so every aromatic edge carbon is legitimately *under* its minimum by the hydrogen it is
about to be given.

```gherkin
Feature: Dope only carbons that can become the nitrogen requested

  @rq-ee235774
  Scenario: A doped structure's nitrogen carries no more bonds than nitrogen can
    Given a request for pyrrolic nitrogen on a structure carrying oxygen groups
    When the structure is generated
    Then every nitrogen in it has a valid valence

  @rq-070a2c21
  Scenario: Pyrrolic nitrogen is still placed where the sites qualify
    Given a structure with unsubstituted five-ring carbons
    When pyrrolic nitrogen is requested
    Then it is placed

  @rq-ca3487ca
  Scenario: A non-aromatic pentagon becomes a pyrrole rather than being skipped
    Given a five-ring whose bonds are kekulé single and double
    When a pyrrolic nitrogen is placed in it
    Then the ring's bonds are aromatic
    And the nitrogen carries two ring bonds and one hydrogen

  @rq-171f79aa
  Scenario: A ring that cannot carry the nitrogen is left as it was
    Given a five-ring that would put an atom over its maximum valence
    When a pyrrolic nitrogen is offered that ring
    Then the ring keeps the bonds and atoms it had
```

## A Bond Is Judged Against the Molecule Being Edited <!-- rq-30fef9a1 -->

`add_bond_safe` takes the editable molecule it writes to. What it must judge against is the state of
that editor, not a snapshot from before the editing began.

Judging against a stale molecule fails in two ways, and the workflow the guide documents — add an
atom, then bond it — hits both. An atom appended to the editor does not exist in the snapshot, so
looking it up raises rather than answering. And a sequence of bonds to one atom is each weighed
against the same original free valence, so four bonds are added to a carbon that had room for two,
which is the exact outcome the helper exists to prevent.

The bond order must also be the order, not the enumeration value that carries it. RDKit numbers its
aromatic bond type 12, so an aromatic bond asked for by name demands twelve free valences and is
refused on every atom in every molecule.

```gherkin
Feature: Judge a bond against the molecule as it now stands

  @rq-e0d6b415
  Scenario: A newly added atom can be bonded
    Given an atom appended to an editable molecule
    When a bond to it is added safely
    Then the bond is added rather than the lookup failing

  @rq-998f84fb
  Scenario: A sequence of bonds to one atom is bounded by its valence
    Given a carbon with room for two more bonds
    When four bonds are added to it one at a time
    Then only the two that fit are added

  @rq-95aaba31
  Scenario: An aromatic bond is weighed as one bond
    Given two aromatic atoms each with a free valence
    When an aromatic bond between them is considered
    Then it is judged against one unit of valence, not the bond type's numeric value
```

## Every Violation Is Named <!-- rq-6d766a3c -->

A validation report lists one message per violating atom, naming the atom, its symbol, the count it
has and the bound it broke. A caller acting on the report needs to know which atom, and a count
without the bound it violated does not say which direction the fault lies in.

```gherkin
Feature: Report each violation individually

  @rq-08130b9d
  Scenario: A molecule with no violations validates
    Given a correctly built molecule
    When it is validated
    Then it reports valid with no messages

  @rq-b8611a4d
  Scenario: Each violating atom is named
    Given a molecule with an atom below its minimum valence
    When it is validated
    Then a message names that atom and the minimum it falls short of
```

## Cross-references <!-- rq-0c391221 -->

- The site selection this module's `available_valence` drives, and the functional-group fallbacks
  it produces, are specified in `heteroatom-assignment.md`.
- The deprotonation whose formal charges this module must respect is specified in
  `ph-protonation.md`.
- The structure-level validation report these violations are folded into is specified in
  `geometry-embedding.md`.
