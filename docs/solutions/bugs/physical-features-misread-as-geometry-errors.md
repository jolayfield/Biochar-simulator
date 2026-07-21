---
title: A clash threshold that overlaps a real physical feature rejects correct structures
date: 2026-07-20
category: bugs
module: geometry_3d
problem_type: bug
component: validation
severity: high
applies_when:
  - "Every point in a sweep degrades to the fallback path"
  - "Adding or tuning a geometry threshold in GeometryValidator or ClashResolver"
  - "Strict validation fails on all max_retries seeds for one region of parameter space"
  - "Deciding whether a reported steric clash is real"
related_components:
  - validation
  - biochar_generator
  - sweep
tags:
  - steric-clash
  - hydrogen-bond
  - vdw-radii
  - false-positive
  - threshold
  - masked-bug
  - fallback
---

# A clash threshold that overlaps a real physical feature rejects correct structures

## Context

`GeometryValidator._check_steric_clashes` scored every atom pair against one rule:

```python
min_distance = 0.75 * (r_vdw_i + r_vdw_j)
```

For an O/H pair that evaluates to `0.75 * (1.52 + 1.20)` = **2.04 Å**. An intramolecular
O–H···O hydrogen bond sits at **1.6–2.2 Å**. The threshold did not sit near the physical
range — it sat *inside* it, so a correct hydrogen bond and a genuine atom overlap were
indistinguishable by construction.

Adjacent phenolic –OH groups on a char sheet form exactly that contact, and MMFF94 forms it
*correctly* during refinement. The validator then rejected the structure the force field had
just improved.

## Symptom

Every 400 °C softwood point in a sweep came back `status="fallback"`. On a cluster this read
as an environment problem; it was not, and it reproduced identically on a laptop.

The tell was that composition was always **on target** — H/C 0.644–0.656 against 0.650, O/C
0.197–0.203 against 0.205, both well inside the 10% tolerance. Only the geometry verdict was
wrong.

The failure scaled with oxygen loading, which is what identified the cause:

| T (°C) | O/C | Formula | False clashes (3 seeds) |
|---|---|---|---|
| 400 | 0.205 | C61H40O12 | 2–5 every seed |
| 500 | 0.120 | C50H22O6 | 0–1 |
| 600 | 0.071 | C50H22O4 | 0 |

## Why seed-retry could not save it

`sweep.build_point` retries up to `max_retries` seeds before falling back. That strategy
assumes the failure is *stochastic*. This one was **deterministic in O/C**: twelve oxygens on
a ~60-carbon sheet make adjacent –OH pairs unavoidable, so the probability of a seed with zero
hydrogen bonds was effectively zero. Eight retries bought nothing but wall-clock.

> When every seed fails identically, suspect a systematic constraint, not bad luck. Retrying is
> only a remedy for variance.

## Fix

Hydrogen-bonded donor–acceptor pairs get a **reduced floor**, not an exemption:

- `_get_hbond_pairs` finds polar H (bonded to N/O) near an N/O acceptor, gated on a D–H···A
  angle so an acceptor jammed into the *side* of a D–H bond is still a clash.
- `_clash_floor` returns `HBOND_MIN_H_ACCEPTOR_DISTANCE` (1.5 Å) for those pairs, the generic
  0.75 × vdW sum for everything else.
- Both `GeometryValidator` and `ClashResolver` use it. Sharing matters: otherwise the resolver
  pushes apart the very H-bonds the next MMFF pass re-forms.

A real overlap between a donor H and an acceptor is still reported, now labelled `type: H-bond`.

This mirrors the pre-existing 1-2/1-3 exclusion in `_get_excluded_pairs`, added because flat-sheet
1-3 C–C contacts at 2.46 Å fell under the 2.55 Å C–C threshold. **Same bug class, second instance.**

## The bug this was hiding

Removing the false clashes surfaced a second defect immediately: seed 0 now failed with an
unrelaxed 1.16 Å aromatic C–C bond.

The force-field refinement in `_generate_geometry` sat behind `if steric_clashes:`. Because every
high-oxygen structure carried at least one false clash, **that pass always ran** — so the coupling
between "was a clash detected" and "does geometry get relaxed" was invisible. Fixing the clash rule
would have silently disabled refinement.

Two masking mechanisms were at work, and both are worth recognising elsewhere:

1. **Gating** — cleanup work hidden behind a condition that happened to always be true.
2. **Truncation** — `StructureValidator` reports only `geom_errors[:3]`, and clashes filled all
   three slots, so the bond-length error never appeared in the report at all.

Refinement is now unconditional on the non-hex-lattice path; clash *resolution* still runs only
when there is a clash to resolve.

## Verifying a change like this

Do not compare against intuition — compare against the old code path directly. Monkeypatching the
new helper to a no-op reproduces the previous behaviour exactly, without stashing files:

```python
import biochar.geometry_3d as G
G._get_hbond_pairs = lambda mol, coords: set()
```

That produced the decisive evidence: rule **off** → 5 clashes, 0 bad bonds; rule **on** → 0 clashes,
1 bad bond. Without it, the new bond error looked pre-existing rather than self-inflicted.

## Related

`_validate_bond_lengths` had a quieter version of the same flaw: it compared against the sum of
*single-bond* covalent radii, so it reported "expected 1.52" for an aromatic C–C that should be
1.40. A wide 0.8–1.5 tolerance band absorbed the error, which is why it went unnoticed. Expected
lengths are now scaled by bond order (`BOND_ORDER_LENGTH_FACTORS`), with tightened tolerance
factors so the corrected — and lower — expectation does not raise the absolute floor.

## Rule of thumb

Before trusting a geometry threshold, ask what real physical features fall in its range. A cutoff
derived from one interaction (vdW overlap) will misjudge any other interaction that occupies the
same distances — hydrogen bonds, peri contacts, ring 1-3 distances. Encode the exception, keep a
reduced bound rather than exempting outright, and make the detector and any resolver share it.
