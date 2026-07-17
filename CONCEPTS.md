# Concepts

Shared domain vocabulary for this project — entities, named processes, and status concepts
with project-specific meaning. Seeded with core domain vocabulary, then accretes as
ce-compound and ce-compound-refresh process learnings; direct edits are fine. Glossary
only, not a spec or catch-all.

## Atom typing

Three distinct names refer to the same atom at different stages of the pipeline. Conflating
them is the root of a recurring class of defect, so the distinction is load-bearing.

### Internal atom type
The generator's own short name for an atom's chemical role — aromatic carbon, ether oxygen,
thioether sulfur, and so on. It is this project's vocabulary, not the force field's, and it
is what the typing stage assigns by inspecting each atom's environment.

An internal type is chosen from chemistry alone and carries no guarantee that the force
field can parameterise it. Every internal type must be translated to an OPLS type name
before export.

### OPLS type name
The force field's identifier for a parameterised atom, written into the exported topology.
This is the entire contract between this project and the force field: the topology asserts
names and delegates all physics.

A name that exists but names the wrong chemistry produces a topology that runs and
simulates the wrong molecule — nothing rejects it. Existence is therefore not evidence of
correctness.

### Bonded type
The coarser key the force field indexes its bond, angle, and dihedral tables by — distinct
from the OPLS type name, and many OPLS type names collapse onto one bonded type.

This is the concept most often missed. An OPLS type name can be the right element with the
right mass while the *combination* of bonded types it participates in has no entry in the
force field's bonded tables, which fails only when a simulation is actually started.

## Structure

### Functional group
An oxygen-, nitrogen-, or sulfur-bearing decoration attached to the aromatic skeleton's
edge to hit a target composition.

Some requested groups cannot be placed on a pure aromatic edge, which has no free valence
for them; these silently degrade to a simpler group rather than failing. A requested
composition is therefore a target, not a promise, and the realised structure must be
measured rather than assumed.

### Ring nitrogen substitution
Replacing a skeleton carbon with nitrogen, as distinct from attaching nitrogen as a
Functional group. The position determines the species: nitrogen at a six-ring edge, in a
five-ring bearing hydrogen, or in the interior are three different chemistries with three
different parameterisations and three different Internal atom types.

## Flagged ambiguities

- "Atom type" had been used for both the Internal atom type and the OPLS type name — these
  are distinct, and a third concept, Bonded type, sits between them. Every layer of the
  export depends on which one is meant.
