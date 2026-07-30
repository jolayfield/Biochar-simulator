# docs: expand the rqm/ requirements set to the modules the 2026-07-30 review identified

**Status:** Phase 0 complete; Phases 1-5 not started
**Type:** docs / traceability
**Date:** 2026-07-30

> **Progress.** Phase 0 landed all five items, including the optional one. Vendoring upstream's
> `test_rqm.sh` immediately earned itself: 27 of its 28 tests pass against our `rqm.sh`, and the
> one failure is real. `stamp --fix-duplicates` is broken on macOS — `rqm.sh:271` is the file's
> only `sed -i`, and BSD `sed -i` requires an argument GNU `sed -i` does not. All nine `stamp`
> tests pass, so the ordinary `rq stamp` path every later phase depends on is unaffected. Not
> patched locally: the vendored files stay byte-for-byte so they can be re-synced, and the bug
> belongs upstream. Documented in `tools/rqm/README.md`.
>
> **D5 resolved: `rr-plan` is vendored.** `SKILL.md` sits at `.claude/skills/rr-plan/SKILL.md`
> byte-for-byte, with every local path redirect and house convention in
> `.riprap/user/skills/rr-plan/local.md` — the extension point upstream reads last and lets take
> precedence. This mirrors the rule already established for `rqm.sh`: vendor verbatim, configure
> beside it.
>
> **One consequence is now live.** CI hard-fails on any scenario without a defending test, so
> from Phase 1 onward a requirements file cannot merge alone — scenario and test land in the same
> PR, with `xfail(strict=True)` when the behaviour does not exist yet. That constraint is written
> into `local.md` so it survives this document.
**Source review:** `docs/reviews/2026-07-30-001-riprap-coverage-and-unadopted-capabilities.md`
**Commit reviewed:** `2ad77a161bf5997e058cf97cf2680632ab99a4b8` (`2ad77a1`)

## Summary

The 2026-07-30 review found that `rqm/` specifies three modules of roughly twenty, and that the
migration's remaining value is coverage rather than the deferred container phase. This plan
sequences that coverage work.

It is built by applying **riprap's `rr-plan` procedure** — fetched from upstream
`MolSSI/riprap@main`, `.riprap/managed/skills/rr-plan/SKILL.md` — to each review finding. The
useful output of that procedure is not prose but a *classification*: for every finding, whether
it becomes a new requirements file, an in-place edit to an existing one, a non-requirements
fix, or nothing at all because it is already specified. Section 3 is that classification and is
the substance of this plan.

## 1. What `rr-plan` is, and what this plan therefore is

Worth stating plainly, because the name misleads. **`rr-plan` is not a project-planning skill.**
Its first line is *"Do not start implementation. Focus only on planning and documentation. Do
not write to anything except requirements file(s)."* It authors and edits `rqm/*.md` — feature
description, optional Feature API, Gherkin scenarios — and then stamps IDs. It produces
specifications, not schedules.

So `rr-plan` cannot itself produce this document. What it produces is each *phase* below. This
document is a plan in the repo's own `docs/plans/` sense — a sequence of PRs — whose individual
units are right-sized `rr-plan` invocations.

Four pieces of `rr-plan` discipline are adopted here and are binding on every phase:

**Right-sizing.** The skill's default action is to *modify an existing requirements file*, not
create one. A new file is justified only when the feature is "large and clearly distinct in
scope from anything described under `rqm/`". A one-bullet change gets a one-bullet diff. This is
why section 3 has an in-place column at all, and why Phase 1 exists.

**Current-state framing.** Requirements describe the current *desired* state and must be
Markovian — readable with no knowledge of what came before. No "this adds", no "previously X,
now Y", no comparison to prior versions. This has a sharp consequence for this plan: several
findings in the review are **defects**. A requirement must never be written as "fix the
`wet.top` ordering bug". It states what is true when the system is right — *the topology the
solvation step consumes exists before that step runs* — and the defect then shows up as a
failing test. Specification first, red test, then fix.

**No placeholder IDs.** Draft with no `rq-` annotations at all; never `@rq-PENDING` or similar.
After writing or editing any requirements file, run in order:

```bash
bash tools/rqm/rq stamp rqm/<file>.md && bash tools/rqm/rq index
```

(Upstream's skill names `.riprap/managed/skills/rr-plan/rqm.sh`; this repo's equivalent is the
`tools/rqm/rq` wrapper. See Phase 0.)

**One scenario, one test.** Scenarios are written so that a single unit test can be constructed
for each. That is what makes back-annotation mechanical rather than interpretive.

## 2. Decisions to make

`rr-plan` requires unresolved material decisions to be surfaced before requirements are written,
and recorded when deferred. These are recorded, not resolved — each changes what gets written.

**D1 — Do we adopt Feature API sections?** `rr-plan` requires them for new files describing
interfaces reachable from elsewhere in the code. This repo's four existing `rqm/` documents have
none, and `rqm/registry.json` holds zero API entries, so half the tooling's entity model is
unexercised. *Recommendation:* adopt for new files covering `export/` and `workflows/`, whose
writer and builder classes are consumed across subpackage boundaries; do not retrofit the three
existing documents.

**D2 — What is the charge-sum tolerance?** The suite asserts charge-sums-to-formal-charge at
fifteen-plus sites using 1e-4, 1e-5, and 1e-6 with no stated rationale. A requirement must pick
one. *This is a numerical-methods judgement, not a documentation choice* — it needs an answer
before Phase 2 writes the `qtot` scenarios.

**D3 — Which annealing peak temperature is correct?** `export/md_setup.py:204-205,244-245`
hardcodes 1000 K for every structure; `workflows/condensation.py:57-61` scales 1000/2000/3000 K
by HTT. Both claim the Wood et al. protocol. **A requirements document cannot be written over
this until someone says which is right** — writing the current `md_setup` behaviour down would
specify a probable bug as intended. Blocks Phase 3.

**D4 — Are improper dihedrals in scope?** The review found none are emitted. Specifying that
aromatic rings carry the impropers that keep them planar is a genuine force-field change with
simulation consequences, not a documentation exercise. *Recommendation:* write the scenario in
Phase 2 and let it fail, rather than quietly specifying the absence as intended.

**D5 — Install `rr-plan`, or follow it by hand?** See Phase 0.

## 3. Decision-tree classification

`rr-plan`'s tree applied to every review finding. This is the core output.

| Review finding | Branch | Destination |
|---|---|---|
| Phenolic fallback for carbonyl/quinone/lactone | **already covered — stop** | `rq-a9492b86` + 4 scenarios, back-annotated |
| Unconditional force-field refinement (solutions doc) | in-place | `rqm/geometry-embedding.md` |
| `geom_errors[:3]` truncation hides defects | in-place | `rqm/geometry-embedding.md` |
| Depth-3 verification omits dihedrals | in-place | `rqm/opls-typing.md` |
| `nrexcl`, impropers, `cgnr`, dihedral dedup | new file | `rqm/gromacs-export.md` |
| Å→nm for atom coordinates, not just boxes | new file | `rqm/gromacs-export.md` |
| `.gro`/`.top` atom-name divergence >10k atoms | new file | `rqm/gromacs-export.md` |
| 5-character residue limit | new file | `rqm/gromacs-export.md` |
| `qtot` / genion charge contract (needs D2) | new file | `rqm/gromacs-export.md` |
| `wet.top` ordering; `-maxwarn 2`; stale genion `.tpr` | new file | `rqm/md-setup.md` |
| Peak-temperature disagreement (blocked on D3) | new file | `rqm/md-setup.md` |
| Config precedence, aromaticity clamp, strict semantics | new file | `rqm/generation-config.md` |
| Sheet spacing, SVD flattening, amorphous fallback title | new file | `rqm/surface-stacking.md` |
| Seed retry is a variance remedy only | new file | `rqm/parameter-sweep.md` (deferred) |
| `get_valid_range` is property-agnostic | new file | `rqm/temperature-model.md` (deferred) |
| CI branch filter; `rq check` warn-only; README staleness | **not a requirement** | Phase 0 |
| `check-secrets.sh` gaps; `importorskip`; `test_pah_quality` | **not a requirement** | Phase 0 / spin-off |
| Changelog and `README.md:665` staleness | **not a requirement** | spin-off |

Three findings are deliberately **not** specified. The `xfail(strict=True)` convention is
process, not product behaviour, and `docs/solutions/` is its right home. The dead
`ATOMIC_MASSES` table and three unused `get_*` helpers in `constants.py` should be deleted
rather than specified — writing a requirement for dead code with a already-divergent sulfur mass
would preserve it. `FUNCTIONAL_GROUPS`' 106 lines of unread `atoms`/`connectivity` payload needs
a decision — wire it up or reduce it to the name registry it actually is — **before** anything
specifies it, or the requirement will specify fiction.

## 4. Phases

One PR per phase, branched from `main`. Per `AGENTS.md`, **do not stack** — a PR based on
another branch runs zero checks while reporting `MERGEABLE / CLEAN`.

### Phase 0 — Make the process runnable; retire the stale config

No requirements written. Everything here has its rationale already recorded elsewhere.

1. Drop `branches: [main]` from the `pull_request` trigger,
   `.github/workflows/ci.yml:6-7`. Rationale in
   `docs/solutions/workflow-issues/stacked-pr-ci-and-base-branch-traps.md`.
2. Promote `rq check` from `|| true` (`ci.yml:103`) to a hard failure, **scoped to
   `type: scenario`**. Unscoped it fails immediately: 30 of the current warnings are `section`
   and `file` headings that no test can reference. The plan's own exit criterion — zero
   unreferenced — is met for scenarios and will never be met for headings.
3. Correct `tools/rqm/README.md:15-18`, which still describes the pre-Phase-3 flat layout and
   claims the wrapper sets `SRC_DIR=biochar`. It sets `SRC_DIR=src`.
4. **Resolve D5.** Upstream `rr-plan` hardcodes `.riprap/managed/skills/rr-plan/rqm.sh`, which
   does not exist here. Riprap's own extension point is
   `.riprap/user/skills/rr-plan/local.md`, which the skill reads last and which takes
   precedence. Adopting the skill therefore means vendoring `SKILL.md` plus a `local.md` that
   redirects every tooling path to `tools/rqm/rq`. *Recommendation:* vendor it — the redirect
   is a handful of lines, and the alternative is that this plan's conventions live only in this
   document and rot.
5. Optionally vendor upstream's `tests/test_rqm.sh` (584 lines, one test per scenario, testing
   the 739-line script this repo already depends on). Cheap, and it removes the review's
   sharpest irony.

### Phase 1 — In-place edits (three scenarios, three small diffs)

The `rr-plan` warm-up: real gaps, trivial diffs, no new files. Each closes a lesson that today
exists only in `docs/solutions/`.

- `rqm/geometry-embedding.md` — force-field refinement runs unconditionally on the
  non-hex-lattice path. It sat behind `if steric_clashes:` and only ever ran because false
  clashes were always present; fixing the clash rule would have silently disabled it. The code
  comment survives at `src/biochar/pipeline/biochar_generator.py:876`.
- `rqm/geometry-embedding.md` — a geometry report names every error, not the first three.
  `geom_errors[:3]` truncation hid the bond-length defect entirely.
- `rqm/opls-typing.md` — depth-3 verification resolves emitted **dihedrals** as well as bonds
  and angles. The solution doc names this gap itself; it is the same failure class that produced
  eight fix commits on angles.

### Phase 2 — `rqm/gromacs-export.md` (new file)

Highest value. Pairs with `rqm/opls-typing.md`, which already owns the argument that only the
`opls_XXX` name reaches GROMACS; this extends it one file downstream, to the writers that emit
those names. Feature API section per D1, covering `GROFileWriter`, `ITPFileWriter`,
`TOPFileWriter`, `GromacsExporter`.

Scenario sketch — each is one unit test:

```
Scenario: The exclusion count matches what OPLS-AA expects
Scenario: An aromatic ring carries the improper terms that keep it planar
Scenario: Two distinct dihedrals over the same four atoms are both emitted
Scenario: Atom names agree between the .gro and the .top for any system size
Scenario: Atom coordinates are written in nanometres
Scenario: A residue name longer than five characters is truncated to five
Scenario: A charged molecule states its charge and names the step that balances it
```

Six existing tests already defend the 5-character limit and can carry IDs with no new test work
(`tests/test_gromacs_export.py:109`, `tests/test_sweep.py:60`, `tests/test_cli.py:111`,
`tests/test_ph_integration.py:364`, `tests/test_generator.py:852`,
`tests/test_surface_builder.py:37`).

Expected to fail on first run, by design: the improper scenario (D4), and the coordinate-unit
scenario — `tests/test_gromacs_export.py:74` asserts `abs(x_nm) < 5.0`, which benzene satisfies
in ångström too, so a real assertion has to replace it.

### Phase 3 — `rqm/md-setup.md` (new file) — **blocked on D3**

About 360 lines of literal `.mdp` physics currently defended by zero value assertions. Writing
this surfaces the `wet.top` ordering defect and the `-maxwarn 2` suppression as scenarios rather
than guesses.

```
Scenario: The topology the solvation step consumes exists before that step runs
Scenario: Every ion species is placed against the coordinates the previous species produced
Scenario: A suppressed grompp warning is named, not counted
Scenario: A structure that only passed as a fallback is marked as such in its run directory
Scenario: An include file the topology needs is resolved from the manifest, not a filename guess
Scenario: The pressure-coupling axis convention preserves the vacuum gap
```

The last one is why D3 blocks: `md_setup` and `condensation` currently use opposite
semi-isotropic conventions and disagree on peak temperature, and a requirement cannot be written
over a contradiction.

### Phase 4 — `rqm/generation-config.md` (new file)

The orchestrator's config resolution — the densest source of user surprise, and the place a
caller's request quietly becomes something else.

```
Scenario: An explicitly requested ratio wins over the temperature-derived one
Scenario: A temperature-derived aromatic fraction below the buildable floor is raised to the floor
Scenario: A partially placed functional group does not satisfy strict validation
Scenario: An unattainable H/C is refused at construction, not at generation
Scenario: pH combined with the ML charge backend is refused rather than silently neutralised
```

### Phase 5 — `rqm/surface-stacking.md` (new file)

`tests/test_ph_integration.py:341-422` is strong enough to carry IDs immediately. The
highest-value scenario is the probable defect: on amorphous→slit fallback, `config.pore_type` is
never updated and the `.gro` is written titled `"Amorphous surface"` while being a slit stack.

```
Scenario: A surface that fell back to slit geometry is not labelled amorphous
Scenario: Sheet centroids are separated by the pore diameter plus the graphene spacing
Scenario: A flattened sheet is planar within embedding noise
Scenario: The .top molecule count matches the number of .itp files written
```

### Deferred

`rqm/carbon-skeleton.md`, `rqm/parameter-sweep.md`, `rqm/temperature-model.md`. All three are
justified — carbon skeleton has five rewrites in history and is first on the project's own
not-yet-specified list — but they are not blocking, and the plan is already five PRs long.

## 5. Out of scope

Riprap's container and launcher layers, unchanged from the original migration plan. `rr-quiz`
and `rr-architecture`: this repo already has a hand-written `rqm/ARCHITECTURE.md` that is better
than a generated one would be. `rr-implement` is deliberately excluded here — it belongs to the
implementation cycle that follows a requirement, not to writing one, and adopting it should be
its own decision after at least one phase above has been through the full write-fail-fix loop.

Fixing the defects this plan specifies is also out of scope *for the requirements PRs*. Each
phase lands a specification and its tests; a failing test is the intended state at merge, marked
with `xfail(strict=True)` per
`docs/solutions/conventions/xfail-strict-as-deferred-gap-tripwire.md`, so that the fix retires
the exemption by XPASS.

## 6. Exit criteria

- `rq check` reports zero unreferenced scenarios and runs as a hard failure in CI.
- `rqm/ARCHITECTURE.md:126-141` — the requirements-set table and the not-yet-specified list —
  reflects the documents that exist.
- Every scenario written by this plan has exactly one test, and every `xfail(strict=True)` left
  behind names the defect it defers.
- `docs/reviews/2026-07-30-001-...` gains a status line pointing at this plan, so a later reader
  can tell which findings were acted on.

## 7. Correction to the source review

Section 1.6 of the review originally listed the carbonyl/quinone/lactone → phenolic fallback as
an unspecified cheap win. It is already fully specified and back-annotated. The review has been
corrected in place with the error recorded rather than deleted. Anyone working from a copy of
the review taken before 2026-07-30 should drop that item.
