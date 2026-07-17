# fix: supply the NC-CA-OH angle for hydroxyl adjacent to pyridinic N

**Status:** active
**Type:** fix
**Date:** 2026-07-17

## Summary

Supply the one bonded-angle parameter stock `oplsaa.ff` lacks for a hydroxyl on a ring
carbon next to a pyridinic N, using the phenol angle as the nearest analog, then retire the
`xfail(strict)` that has been holding the gap. This closes the last known depth-3 gap.

---

## Problem Frame

`tests/test_opls_type_map.py::TestBondedResolution::test_ring_nitrogen_emits_only_resolvable_terms[num_pyridinic]`
is `xfail(strict)` on main. A pyridinic-N-doped sheet that also carries a phenolic hydroxyl
(which the default O/C ratio adds) puts an `-OH` on a ring carbon adjacent to the pyridinic
N. That emits an `NC-CA-OH` angle — bonded types `NC` (pyridine N), `CA` (aromatic C), `OH`
(hydroxyl O). Stock `ffbonded.itp` has no such angletype, so `grompp` would reject the
topology with "No default Angle types".

This is a genuine force-field gap, not a typing bug: 3-hydroxypyridine is real chemistry
OPLS simply omits. It is therefore the legitimate case for `SUPPLEMENTARY_ANGLE_PARAMS` —
the same mechanism that already supplies `CA-S-CA` for thioethers. The hard part was never
the code; it was choosing a defensible value with provenance. That is now settled.

The convention governing this table is
`docs/solutions/conventions/verify-opls-types-against-real-forcefield.md`: the supplement may
hold **only** combinations stock `oplsaa.ff` does not define, every entry needs a provenance
comment, and `TestSupplementDoesNotShadowForcefield` enforces the no-duplication half.

---

## Key Technical Decisions

**KTD1 — Value derives from the phenol angle: `CA-CA-OH` = `120.000` deg / `585.760`
kJ/mol/rad^2.** Verified present in `ffbonded.itp` `[ angletypes ]`. It is the closest stock
analog: it holds the hydroxyl-on-aromatic-carbon geometry (`C-O` bond, ring-`C`-bearing-`OH`)
exactly, and differs from the target only in the far ring atom (neutral `CA` vs. pyridinic
`NC`). The value goes into the table verbatim — no averaging or adjustment — so its
provenance is a single traceable stock entry.

Rejected alternative: a pyridine-ring-derived angle (e.g. deriving from `NC-CA-HA`, 116.000).
That keeps the nitrogen's ring environment but changes the *other* leg of the angle (the
`H`/`O` swap), so it is not obviously closer, and it would be a derived rather than a
directly-transcribed value — weaker provenance. Deferred; the phenol value is the honest
nearest analog and is marked, like the whole supplement, as not QM-validated.

**KTD2 — Scope is exactly the combination the check caught: `NPY-CA-OH`.** Do not
pre-emptively supplement sibling combinations (ether or carbonyl adjacent to pyridinic N, or
any functional group adjacent to `NGR`/`NPR`). The depth-3 check varies one doping mode at a
time against the default composition, so it currently cannot see those siblings; adding
speculative entries would put values into the table that no test exercises — exactly the
untested-drift the convention warns against. If a sibling gap is real it will surface as its
own `xfail` when the check is later widened, and get its own provenanced value then.

**KTD3 — Supplement entry and xfail removal land in one atomic change.** Adding the
supplement while the `xfail(strict)` marker remains would make the test XPASS and fail the
suite; removing the marker before adding the supplement leaves a hard failure. They are two
files but one logical change and must be committed together.

---

## Implementation Units

### U1. Supply the NC-CA-OH angle and retire its xfail

**Goal:** `grompp` can resolve every angle a pyridinic-N + phenolic-OH sheet emits, and the
depth-3 check proves it rather than expecting failure.

**Requirements:** Closes the `num_pyridinic` depth-3 gap (KTD1, KTD2, KTD3).

**Dependencies:** none.

**Files:**
- `biochar/constants.py` — add one entry to `SUPPLEMENTARY_ANGLE_PARAMS`.
- `tests/test_opls_type_map.py` — remove the `_PYRIDINIC_XFAIL` constant and the
  `pytest.param("num_pyridinic", marks=...)` wrapper, so `num_pyridinic` becomes a plain
  parametrize entry alongside `num_pyrrolic` and `num_graphitic`.

**Approach:**
- The new entry is keyed by internal atom type, outer two sorted, matching the existing
  `("CA", "SS", "CA")` entry and the keying `_angle_suffix` performs: central atom is the
  ring carbon `CA`; the outer pair is the pyridinic N internal type and the hydroxyl `OH`.
  Value is the phenol angle from KTD1 in GROMACS units `(theta0_deg, k_kJ/mol/rad^2)`.
- The entry carries a provenance comment in the same shape as the `CA-S-CA` one: what the
  force field lacks, which stock angle the value is transcribed from (`CA-CA-OH`, phenol),
  and that it is a nearest analog, not QM-validated.
- Confirm the exact internal type on the pyridinic ring carbon's nitrogen neighbour before
  writing the key. The check parametrizes `num_pyridinic`, which routes through
  `NitrogenSubstitutor` to the `NPY` internal type; verify by reading the emitted
  `atom_types` for a generated pyridinic sheet rather than assuming, since the key must match
  what `_angle_suffix` looks up at export time.

**Patterns to follow:** the existing `("CA", "SS", "CA"): (104.200, 518.816)` entry in
`biochar/constants.py` and its provenance comment — same table, same comment shape, same
"approximate, not QM-validated" honesty.

**Test scenarios:**
- `test_ring_nitrogen_emits_only_resolvable_terms[num_pyridinic]` now **passes** (no longer
  xfail): a generated pyridinic-N sheet at default composition emits no angle absent from
  both `ffbonded.itp` and `SUPPLEMENTARY_ANGLE_PARAMS`.
- `TestSupplementDoesNotShadowForcefield::test_no_supplementary_angle_duplicates_a_stock_angletype`
  still passes: the new `NPY-CA-OH` entry resolves to bonded `NC-CA-OH`, which is confirmed
  absent from stock `[ angletypes ]`, so it shadows nothing.
- Regression guard (informal, verify during work): temporarily reverting just the supplement
  entry makes `[num_pyridinic]` fail with `NC-CA-OH` in the message — i.e. the entry is what
  closes the gap, not some unrelated change. Do not leave this in the suite.
- Export smoke: a pyridinic-N sheet carrying a phenolic OH writes an `[ angles ]` line for
  the `NPY-CA-OH` triple with the two explicit parameters appended, in both `.top` and
  `.itp` (mirrors how `CA-S-CA` was verified). `funct`-only lines remain for every resolvable
  angle.

**Verification:** `BIOCHAR_OPLSAA_FF=<stock ff> pytest tests/test_opls_type_map.py` runs with
**zero xfails remaining** (was 1). Full suite green. The depth-3 guard now runs in CI (it has
`gromacs-data`), so this also holds on CI, not only locally.

### U2. Reconcile the convention doc

**Goal:** the convention doc reflects that `NC-CA-OH` is supplied, not open — so it does not
describe a gap that no longer exists.

**Requirements:** documentation consistency with U1.

**Dependencies:** U1 (lands in the same commit or immediately after).

**Files:**
- `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md` — the
  "What the check found immediately" section currently lists `num_pyridinic` / `NC-CA-OH` as
  a live gap "with a defensible value, which needs its own change." Update it to past tense:
  the gap was real, it was supplied from the phenol analog, and it now stands as a worked
  example of the supplement's intended use rather than an open item.

**Approach:** small prose edit. Keep the entry as a teaching example — it is a good
illustration of "genuine force-field gap, so it belongs in the supplement," which contrasts
with the carboxyl case ("typing bug, so it does NOT"). Do not delete the example; convert it
from open-gap to resolved-example. Preserve the CI-skip caveat elsewhere in the doc (still
partly true — `test_constants_ff.py` runs in CI, the forcefield-backed tests now do too via
`gromacs-data`).

**Test scenarios:** `Test expectation: none -- documentation-only edit.`

**Verification:** `python3` frontmatter validator passes on the doc; no remaining prose calls
`NC-CA-OH` an open gap.

---

## Scope Boundaries

**In scope:** the single `NPY-CA-OH` angle, its xfail retirement, and the doc reconciliation.

### Deferred to Follow-Up Work
- **Widening the depth-3 check to combinations** (nitrogen-doping × functional-group
  crossed, rather than one variable at a time). This is the real way to find sibling gaps
  like ether-next-to-pyridinic-N or any functional group next to `NGR`/`NPR`. It is a
  larger, separate change with its own value-provenance work per gap found, and was
  explicitly held out of this fix. Worth its own plan.
- **QM validation** of the supplement's values (`CA-S-CA` and now `NPY-CA-OH`). Both are
  nearest-analog transcriptions, chemically reasonable but not publication-grade.

---

## Risks & Dependencies

- **Wrong key, silent miss.** If the supplement key does not match what `_angle_suffix`
  looks up (wrong internal type, wrong outer-pair ordering), the entry is inert and
  `[num_pyridinic]` stays failing after the xfail is removed — a loud failure, caught
  immediately by U1's own test run. Mitigated by reading the generated `atom_types` before
  writing the key (U1 Approach).
- **CI now has teeth.** Because `gromacs-data` is installed in CI, an unresolvable angle no
  longer only fails locally — it fails the PR. This is a benefit, but it means the branch's
  own CI is the real gate: confirm green before merge, not just the local run.
