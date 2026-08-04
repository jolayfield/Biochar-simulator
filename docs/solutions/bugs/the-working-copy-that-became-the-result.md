---
title: A working copy made for the force field became the structure that was returned
date: 2026-08-03
category: bugs
module: geometry_3d
problem_type: bug
component: embedding
severity: high
applies_when:
  - "Every ring carbon in a structure reports Valence 3 below minimum 4"
  - "A finished structure has no aromatic bonds but its atoms are flagged aromatic"
  - "Rewriting a molecule to satisfy an RDKit API and continuing with the rewrite"
  - "A validator, report, or check is handed the caller's own molecule"
related_components:
  - validation
  - opls_typing
  - valence
tags:
  - kekulization
  - aromaticity
  - working-copy
  - mutation
  - rdkit
  - silent-corruption
---

# A working copy made for the force field became the structure that was returned

## Context

RDKit's embedders and force fields need integer bond orders. A large fused biochar sheet
frequently has no kekulé structure at all, so `geometry_3d._kekulize_or_dearomatize` falls back
to stripping **every** aromatic bond to `SINGLE` and clearing **every** aromatic flag. That is a
sound scaffold for computing coordinates. It is not a molecule anyone asked for.

## Symptom

`ValidationEngine` reporting `Atom N (C): Valence 3 below minimum 4` for *every* ring carbon at
once — 57 or 61 errors in a single structure, not a handful:

```
Validation FAILED with 61 error(s): Atom 0 (C): Valence 3 below minimum 4;
Atom 1 (C): Valence 3 below minimum 4; Atom 2 (C): ... [58 more]
```

12 of 40 seeds at `target_num_carbons=60, defect_fraction=0.3, O_C_ratio=0.15`. The tell is the
count: a real valence bug touches one atom, and this touched all of them. The other tell is the
log line above it, `Can't kekulize mol`.

## Cause — two independent places, one symptom

**One.** `generate_3d_coordinates` set `best_mol = mol_copy` and returned that copy. Coordinates
were the result; the bond orders came along with them. Everything downstream then read a sheet
whose chemistry was gone: valence validation counted aromatic carbons as three-bonded, atom
typing fell through to its ring-and-degree proxy for "is this an aromatic carbon" instead of
reading the flag, and `_enforce_aromatic_planarity` found no aromatic ring to be planar about.

**Two.** `ChemicalFeasibilityValidator.validate` asked its connectivity question by calling
`Chem.SanitizeMol(mol)` on the caller's molecule. Sanitisation is not a read — it kekulises, and
kekulisation of a sheet with no kekulé structure rewrites aromatic bonds to single *on its way to
raising*. The exception was caught and turned into a warning; the damage stayed. Worse, the
validator's own valence check runs before the sanitisation, so a single call reported the
structure sound and handed back one that was not:

```
   H=    {'A': 80, 'S': 45}/ar62     <- 80 aromatic bonds entering geometry
   geom= {'A': 80, 'S': 45}/ar62     <- fix one: they survive geometry
   valid={'S': 125}/ar62             <- fix two needed: validation flattened them
```

Fixing only the first left a molecule in a state neither half produces on its own — atoms still
flagged aromatic, every bond single — which passes the valence check (it reads the flag) while
being wrong in the topology.

## Fix

**Embedding returns the molecule it was given.** `_adopt_conformer(mol, best_mol)` copies the
computed coordinates onto the caller's molecule; the rewritten copy stays inside the function.

**Refinement builds its own working copy.** `validate_and_relax` returns coordinates and nothing
else, so nothing of the caller's molecule needs to survive into it. It now calls
`_kekulize_or_dearomatize` itself. This is not optional once embedding stops handing its copy
back: MMFF and UFF both refuse aromatic bond orders, both refusals are swallowed, and the pass
would have returned the coordinates untouched and refined-looking.

**Validation sanitises a copy.** `Chem.SanitizeMol(Chem.Mol(mol))` answers the same question
about the same graph without answering it destructively.

Across 40 seeds: structures with valence errors 12 → 0, and every structure keeps its aromatic
bonds through to export.

## Prevention

- A representation built to satisfy an API is a working copy. Name it as one and keep it inside
  the function that made it — the moment it is returned, everything downstream inherits a
  chemistry nobody chose.
- A read that mutates is worse than a write, because nothing downstream expects to re-perceive.
  `SanitizeMol`, `Kekulize`, and `SetAromaticity` all rewrite in place; asking them a question
  about a molecule you do not own means asking a copy.
- Check the *count* in a validation failure. Errors on every atom of a class are a state
  problem, not a chemistry problem, and the atom named first is rarely the one to look at.
- When two independent mutations produce one symptom, fixing either one leaves a state that is
  neither the bug nor correct. Trace the invariant through every stage before believing a fix —
  the stage-by-stage bond-type dump above is what showed the second one existed.

## See also

- `rqm/geometry-embedding.md` — `rq-2a552a1a` (embedding returns what it was given) and
  `rq-0066e1ac` (validation does not rewrite what it validates).
- `docs/solutions/bugs/failed-sanitisation-forgets-the-rings.md` — the other way a swallowed
  sanitisation failure damaged the molecule it was handed.
