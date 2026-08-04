---
title: A substitution that edits the atom but not its ring inherits the ring's bond orders
date: 2026-08-03
category: bugs
module: heteroatom_assignment
problem_type: bug
component: nitrogen-doping
severity: high
applies_when:
  - "A neutral nitrogen is reported at four bonds after pyrrolic doping"
  - "Adding or tightening a site rule for ring substitution"
  - "A site screen fixes the valence but the placement rate collapses"
  - "Reasoning about which parts of a biochar skeleton are aromatic"
related_components:
  - valence
  - geometry_3d
  - opls_typing
tags:
  - nitrogen-doping
  - pyrrolic
  - aromaticity
  - kekulization
  - valence
  - site-selection
---

# A substitution that edits the atom but not its ring inherits the ring's bond orders

## Context

`NitrogenSubstitutor._substitute_pyrrolic` turns a five-ring carbon into a nitrogen and
attaches an explicit N–H. Three σ bonds is the whole budget for a neutral nitrogen: two ring
neighbours and the hydrogen.

Two site rules had already been added to protect that budget — the carbon must be undecorated
with exactly two carbon neighbours, and only one nitrogen per ring. Both are about the *atom*.
The remaining failure was about the *ring*, and no atom-level rule could see it.

## Symptom

`ValenceValidator.validate_molecule` reported `Atom NN (N): Valence 4 exceeds maximum 3` on
roughly 40% of seeds, at:

```python
GeneratorConfig(target_num_carbons=60, seed=0, strict=False, defect_fraction=0.3,
                H_C_tolerance=1.0, O_C_tolerance=1.0, num_pyrrolic=2, O_C_ratio=0.15)
```

The nitrogen came out of `_substitute_pyrrolic` already wrong. The prevailing theory — that
the pentagon lost its aromaticity somewhere downstream, in H/C shaping or one of the
`_safe_sanitize` passes — was wrong, and a stage-by-stage trace of the nitrogen from
substitution to export is what disproved it: the bonds were `SINGLE`/`DOUBLE` from the moment
they were written and never changed afterwards.

## Cause

A biochar skeleton is not aromatic everywhere. On the failing seeds the skeleton carried a
non-aromatic pocket — 13 of 65 atoms, with four `DOUBLE` bonds — and both pentagons sat
inside it. Aromaticity perception refuses a cyclopenta-fused pentagon (it is cross-conjugated,
as in acenaphthylene), so its bonds are the kekulé single and double the perception fell back
to:

```
pentagon: C52[D,S]  C14*[A,A,D]  C13*[A,A,S]  C54[S,D]  C53[S,D]      (* = aromatic)
```

Every carbon with two neighbours — every carbon the site rule would offer — holds one of the
ring's double bonds. Substituting any of them gives a nitrogen carrying a C=N, and the N–H
makes four. No choice of site avoids it, which is exactly why an atom-level screen could not.

`_candidate_carbons` was documented as returning "aromatic ring carbons" and only ever checked
ring membership, which is how the gap survived review.

## The fix that did not work

Screening the site on aromaticity — `atom.GetIsAromatic()` and all bonds `AROMATIC` — makes
every nitrogen valid and drops placement from 14 to 6 across the seven seeds the contract
pins, failing `rq-070a2c21`, the counter-requirement that exists to catch exactly this. The
rule was not too strict by accident: on those skeletons *no* qualifying site exists.

## Fix

Change the ring instead of rejecting it. A cyclopenta-fused pentagon is not aromatic; the same
ring with an N–H in it is a **pyrrole**, which is, and the nitrogen's own lone pair is what
aromatises it. So the substitution puts the ring into the state the new nitrogen implies —
the pentagon's bonds become `AROMATIC`, the double bond the site was holding included — and
the carbon that gave up that double bond is left an aromatic carbon with a free valence for
the hydrogen `HydrogenAssigner` will give it.

Whether that lands is checked rather than assumed: the substitution happens on a copy, the
whole pentagon is read back, and a ring that comes out with an atom over its maximum is
discarded untouched. A pentagon holding an sp3 carbon from aliphatic decoration is the case
that check catches — no amount of flagging makes that carbon aromatic.

Result across 40 seeds: nitrogen violations 17 → 0; placement across the pinned seeds 12 of a
possible 14.

## The trap in the read-back

The obvious check — `ValenceValidator.get_valence_info(...).is_valid` on each ring atom —
rejects every ring in the molecule. Nitrogen substitution runs *before* hydrogen saturation,
so an aromatic edge carbon is legitimately **under** its minimum by exactly the hydrogen it is
about to be given. Only the upper bound is meaningful at this point in the pipeline:

```python
v.current_bonds > v.max_valence
```

This cost one full debugging cycle: the fix was correct, every nitrogen came out at 3/3, and
the placement rate still read 4 of 14 because the *carbons* were being rejected.

## Prevention

- A site rule protects a budget the substitution spends. When the budget is spent on bonds the
  site did not choose — ring bond orders, fused-ring context — the rule cannot be the whole
  answer, and tightening it trades one requirement for its counter-requirement.
- Trace the defect to the stage that writes the bad state before theorising about which later
  stage corrupted it. Here the theory named three stages, and the state was wrong before any
  of them ran.
- A helper whose docstring claims a check it does not make is a defect with a delayed fuse.
  `_candidate_carbons` said "aromatic ring carbons" for as long as every skeleton happened to
  be aromatic.

## See also

- `rqm/valence-validation.md` — `rq-53ec84de`, and the scenarios `rq-ee235774`, `rq-070a2c21`,
  `rq-ca3487ca`, `rq-171f79aa`.
- `docs/solutions/bugs/failed-sanitisation-forgets-the-rings.md` — the crash this defect used
  to cause three stages downstream.
