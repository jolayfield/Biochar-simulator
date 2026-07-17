---
title: Verify OPLS type mappings against a real oplsaa.ff, at all three depths
date: 2026-07-16
category: conventions
module: constants
problem_type: convention
component: testing_framework
severity: high
applies_when:
  - "Adding or changing an entry in GROMACS_OPLS_TYPE_MAP"
  - "Adding a functional group that introduces a bond or angle never emitted before"
  - "Adding any force-field constant to biochar/constants.py"
  - "Reviewing a claim that force-field data has been verified"
related_components:
  - opls_typing
  - gromacs_export
  - constants
tags:
  - opls-aa
  - force-field
  - gromacs
  - ffbonded
  - silent-failure
  - verification
---

# Verify OPLS type mappings against a real oplsaa.ff, at all three depths

## Context

An exported `.top` carries **no force-field parameters of its own**. `TOPFileWriter.write`
emits `[ bonds ]` and `[ angles ]` with a `funct` code and no `c0`/`c1`:

```python
f.write(f"{i:5d} {j:5d} 1\n")                    # bond: indices + funct, no params
f.write(f"{n1+1:5d} {j+1:5d} {n2+1:5d} 1\n")     # angle: same
```

So 100% of the physics comes from the `#include "oplsaa.ff/forcefield.itp"`, resolved by
the `opls_XXX` names written into `[ atoms ]` via `GROMACS_OPLS_TYPE_MAP`. Those names are
the entire contract.

That contract has now been checked at three successive depths, and **each depth was
believed sufficient until the next defect proved it wasn't**:

| Depth | Check | Defeated by |
|-------|-------|-------------|
| 1 | The mapped type **exists** | Every wrong mapping named a real type |
| 2 | The mapped type has the right **element/mass** | `CA-S-CA` — right element, right mass, still broken |
| 3 | Every emitted **bond/angle combination resolves** in `ffbonded.itp` | (holds so far) |

**Depth 1** (fixed in c26bacd): `SS` (thioether sulfur) mapped to `opls_209`, which is a
*carbon* — bonded type `CT`, mass 12.011. Pyridinic and pyrrolic N mapped to carbon types
too. `grompp` accepts any name that exists, so nothing errored; the simulation ran the
wrong chemistry.

**Depth 2**: the fix added element/mass consistency checks
(`TestTypesAreTheRightElement` in `tests/test_opls_type_map.py`). This is the current
ceiling.

**Depth 3** (open): `SS -> opls_222` is genuinely sulfur (mass 32.06, bonded type `S`), so
depth 2 passes cleanly — yet a thioether bridges two aromatic carbons and emits a
`CA-S-CA` angle, and `ffbonded.itp` has **no such angletype**. Only `CA S CT` (thioanisole)
and `CA S CM` exist, both `104.200` deg / `518.816` kJ/mol/rad². `grompp` rejects the
topology with "No default Angle types". Reproduced end-to-end with
`functional_groups={"thioether": 2}`. Of everything the generator emits, this is the only
hole — all bonds and every other angle resolve.

## Guidance

When a generator maps internal types onto an external force field and relies on that force
field to supply the parameters, verification must reach all three depths. They catch
**non-overlapping** failure classes; passing depth 2 says nothing about depth 3.

The depth-3 lookup needs an indirection that is easy to get wrong: `ffbonded.itp` is keyed
by **bonded type** (column 2 of `ffnonbonded.itp`), *not* by the `opls_XXX` name. Several
opls names collapse onto one bonded type — `opls_145` and `opls_521` are both `CA`. The
chain is:

```
internal type -> GROMACS_OPLS_TYPE_MAP -> opls_XXX
              -> ffnonbonded.itp column 2 -> bonded type
              -> ffbonded.itp [bondtypes] / [angletypes]
```

Corollary: **never hand-copy force-field values into the codebase.** Three tables that did
(`OPLS_LJ_PARAMS`, `OPLS_BOND_PARAMS`, `OPLS_ANGLE_PARAMS`) drifted into an AMBER/OPLS
chimera with no single provenance and were deleted. When the force field genuinely lacks a
parameter, supply it from a local `.itp` with a provenance comment — do not reinstate a
table.

## Why This Matters

The two failure modes are asymmetric, and the dangerous one is silent:

- **Wrong-but-existing type name** — runs to completion, emits a `.top` that looks
  perfectly normal, and simulates the wrong molecule. Nothing flags it: not `grompp`, not
  the suite, not reading the output. It surfaces only as physically implausible results,
  if anyone notices. This shipped, undetected, for an unknown period.
- **Missing bonded combination** — `grompp` refuses outright. Loud, but only at simulation
  time, which the test suite never reaches. A contributor can add a functional group, watch
  every test pass, and hand a user a topology that cannot start.

Both trace to the same design: the `.top` asserts names and delegates all physics. The
names are unverified data pretending to be code.

## When to Apply

- Adding or changing any `GROMACS_OPLS_TYPE_MAP` entry.
- Adding a functional group or heteroatom that can introduce a bond/angle combination never
  emitted before — the new combination has never been checked by anything.
- Any time "verified against the force field" is claimed: ask *at what depth*. The phrase
  has now been true-but-insufficient twice.
- Generally: when a data-driven mapping feeds a downstream lookup table, the *combinations*
  are a separate verification surface from the *entries*.

## Examples

### The depth-3 test

Slots into `tests/test_opls_type_map.py`, reusing its existing `_find_oplsaa()` /
`requires_oplsaa` skip machinery. Verified against the stock force field — it reports
`CA-S-CA` and nothing else (no false positives):

```python
from itertools import combinations


def _parse_bonded_types(ff):
    """opls_XXX -> bonded type (column 2 of ffnonbonded.itp)."""
    out = {}
    for line in (ff / "ffnonbonded.itp").read_text().splitlines():
        parts = line.split(";", 1)[0].split()
        if len(parts) >= 2 and parts[0].startswith("opls_"):
            out[parts[0]] = parts[1]
    return out


def _parse_section(ff, header, n):
    """[ bondtypes ] is `i j func b0 kb`; [ angletypes ] is `i j k func th0 cth`."""
    out, inside = set(), False
    for raw in (ff / "ffbonded.itp").read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = line.startswith(header)
            continue
        if not inside or not line:
            continue
        p = line.split()
        if len(p) >= n:
            # bonds are order-free; angles fix the middle atom
            out.add(tuple(sorted(p[:2])) if n == 2
                    else (min(p[0], p[2]), p[1], max(p[0], p[2])))
    return out


@requires_oplsaa
class TestBondedResolution:
    """Depth 3: element and mass can both be right while the combination the
    molecule actually emits is still missing from ffbonded.itp."""

    def test_every_emitted_bond_and_angle_resolves(self):
        bonded = _parse_bonded_types(OPLSAA)
        bondtypes = _parse_section(OPLSAA, "[ bondtypes", 2)
        angletypes = _parse_section(OPLSAA, "[ angletypes", 3)

        def to_bonded(internal):
            return bonded[GROMACS_OPLS_TYPE_MAP[internal]]

        gen = BiocharGenerator(GeneratorConfig(
            target_num_carbons=20,
            functional_groups={"thioether": 2},
            strict=False,          # composition tolerances are not what we assert here
        ))
        mol, _, _ = gen.generate()
        types = gen.atom_types          # instance attr, set by generate()

        missing = []
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if tuple(sorted([to_bonded(types[i]), to_bonded(types[j])])) not in bondtypes:
                missing.append(("bond", types[i], types[j]))

        for atom in mol.GetAtoms():
            nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
            j = to_bonded(types[atom.GetIdx()])
            for x, y in combinations(nbrs, 2):
                i, k = to_bonded(types[x]), to_bonded(types[y])
                if (min(i, k), j, max(i, k)) not in angletypes:
                    missing.append(("angle", f"{i}-{j}-{k}"))

        assert not missing, f"combinations absent from ffbonded.itp: {missing}"
```

**Two caveats that must ship with this test, not be discovered after:**

1. **It fails on day one.** `CA-S-CA` is a real, currently-shipping gap, so this cannot
   land as a bare assertion — it lands *with* the fix (a local `.itp` supplying the angle;
   closest stock value is the `CA-S-CT` thioanisole angle, `104.200` / `518.816`). The
   day-one failure is the test working, not the test being wrong.
2. **It has no teeth in CI.** CI never installs GROMACS, so every `requires_oplsaa` test
   skips on every run. Until CI gets a force field — vendor a fixture, or set
   `BIOCHAR_OPLSAA_FF` — this only protects people who run it locally with one installed.

### Two traps that look like the answer

Both are tempting, confident, and wrong. Recorded because the next reader will reach for
them:

- **"The HA/HC epsilon gap is a clean factor of 2 — it's a kcal/kJ conversion bug."**
  It isn't. Every epsilon in the deleted table decodes as a round number in kcal/mol, so
  the conversion was fine. `HA`'s `(0.2600, 0.0630)` is *exactly* AMBER parm94's HA;
  `NA` is exactly AMBER's N; `CA`/`CT`/`C` are genuine OPLS. The table was a chimera, not a
  mis-conversion — there was no single error to undo, which is precisely why it was deleted
  rather than corrected.
- **"A wrong force constant means the table's convention is wrong."**
  Also no. The bond table's convention was *sound* — OPLS literature `kr` = half the
  GROMACS `kb` — and `CA-HA` (367), `CT-CT` (268), `CA-OH` (450) matched the stock file
  exactly. Individual values were still wrong (`CA-S` 340 vs 250; `CA-CA` 770 vs 469).
  Internal self-consistency across most rows is not evidence any row is right.

## Related

- PR [#20](https://github.com/jolayfield/Biochar-simulator/pull/20) — depth-1 fix
  (wrong-element mappings); added `tests/test_opls_type_map.py`.
- PR [#21](https://github.com/jolayfield/Biochar-simulator/pull/21) — deleted the three
  hand-copied tables and the unreachable `[ atomtypes ]` fallback.
- Issue #3 — S-doping (thiol/thioether); origin of the types behind the `CA-S-CA` gap.
- `biochar/constants.py` — `GROMACS_OPLS_TYPE_MAP` and the `SS -> opls_222` note, which
  records the `CA-S-CA` gap and its closest stock value.
- `docs/api/constants.rst` — states that LJ/bonded parameters live in `oplsaa.ff`, not in
  this module.
- `docs/user_guide/gromacs_workflow.rst` — force-field notes; does not yet mention the
  thioether caveat.

**Open reconciliation risk:** branch `feat/ph-protonation` (`b21a8bf`) independently redid
these mappings and **regresses `SS` back to `opls_209`** — the carbon that depth 1 fixed.
If it merges unreconciled, the wrong-element bug returns. The depth-3 test would catch it.

**Note on a divergent test file:** notes from the `feat/ph-protonation` line refer to
`tests/test_constants_ff.py` for the element/mass check; the merged lineage names it
`tests/test_opls_type_map.py`. Same purpose, different branches. Verify which exists before
relying on either. *(auto memory [claude])*
