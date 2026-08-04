---
title: A swallowed SanitizeMol failure leaves the molecule without ring information
date: 2026-08-03
category: bugs
module: geometry_3d
problem_type: bug
component: embedding-preparation
severity: high
applies_when:
  - "RuntimeError: Pre-condition Violation / RingInfo not initialized"
  - "Wrapping Chem.SanitizeMol in try/except and continuing"
  - "A crash reports a missing precondition rather than the thing that went wrong"
  - "Adding a GetRingInfo() call to a pipeline stage"
related_components:
  - opls_typing
  - heteroatom_assignment
tags:
  - rdkit
  - sanitization
  - ring-perception
  - error-handling
  - action-at-a-distance
---

# A swallowed SanitizeMol failure leaves the molecule without ring information

## Context

`geometry_3d._kekulize_or_dearomatize` prepares a molecule for embedding: kekulise it, or —
when the graph will not kekulise, which is routine for large biochar PAHs — strip aromaticity
and sanitise what is left. The sanitisation is allowed to fail. A graph carrying one bad
valence is still a graph the force field can be handed, so the exception is swallowed and the
geometry goes ahead.

## Symptom

```
RuntimeError: Pre-condition Violation
	RingInfo not initialized
	Violation occurred on line 28 in file .../GraphMol/RingInfo.cpp
```

raised from `opls_typing.AtomTyper._is_aromatic_carbon`, three pipeline stages after the
molecule was prepared, on some seeds and not others. Atom typing reports it because typing is
the first stage that asks every atom whether it is in a ring.

## Cause

`Chem.SanitizeMol` clears the molecule's **computed properties — ring information among
them — before it does any of its work**, and its valence check (`SANITIZE_PROPERTIES`) runs
before ring perception (`SANITIZE_SYMMRINGS`). A sanitisation that raises on the valence check
therefore returns a molecule that knows *less* about itself than the one that went in:

```python
def ri(m):
    try:    return m.GetRingInfo().NumRings()
    except RuntimeError: return "UNINIT"

# seed 5, at each pipeline stage
skel=ok(19)  O=ok(19)  N=ok(19)  prot=ok(19)  H=ok(19)  geom=UNINIT
```

Nothing near the failure notices, because nothing near the failure asks about rings. The
molecule embeds, coordinates come out, and the generation dies later reporting a missing
precondition instead of the valence error that actually caused it.

On the seeds where this fired, the valence error was a pyrrolic nitrogen at four bonds — see
`substitution-that-edits-an-atom-but-not-its-ring.md`. Two defects, one visible symptom, and
the visible one named neither the module nor the atom at fault.

## Fix

`_kekulize_or_dearomatize` perceives rings before returning, on both branches:

```python
def _with_rings_perceived(mol):
    try:
        mol.GetRingInfo().NumRings()
    except RuntimeError:
        Chem.FastFindRings(mol)
    return mol
```

`FastFindRings` answers ring membership without sanitising, so it cannot disturb the bond types
the heteroatom stages worked to set. It is idempotent — a molecule that already knows its rings
is returned untouched.

`AtomTyper.assign_atom_types` carries the same guard (`rq-c6ab7cbe`) and keeps it. That one is
a backstop for any molecule arriving unperceived from anywhere; this one removes the source.

## Prevention

- `except Exception: pass` around `SanitizeMol` is not "carry on unchanged". Sanitisation
  mutates before it validates, so the swallowed path leaves state the happy path never
  produces. Restore the invariant explicitly in the handler rather than assuming the input
  survived.
- When a crash names a missing precondition, the stage that raised it is rarely the stage that
  broke it. Walking the invariant backwards stage by stage found the real site in one pass;
  reading the traceback did not.
- Guarding at the point of use is a legitimate backstop and is not a fix. It converts a crash
  into silence, which is worse if the underlying state is wrong for other reasons too.

## See also

- `rqm/geometry-embedding.md` — `rq-07eab6b8` and its scenarios `rq-a157df17`, `rq-ac433b55`.
- `rqm/opls-typing.md` — `rq-c6ab7cbe`, the backstop in atom typing.
