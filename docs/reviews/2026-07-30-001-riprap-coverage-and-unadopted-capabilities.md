# review: what in this codebase would still benefit from the riprap framework

**Type:** review — point-in-time analysis, not a plan
**Date:** 2026-07-30
**Acted on by:** `docs/plans/2026-07-30-001-docs-rqm-coverage-expansion-plan.md`, which classifies
every finding below via riprap's `rr-plan` decision tree
**Commit:** `2ad77a161bf5997e058cf97cf2680632ab99a4b8` (`2ad77a1`, main, clean tree)
**Scope:** the whole repository, not the phases in
`docs/plans/2026-07-28-001-refactor-riprap-scaffold-migration-plan.md`

> **What this document is.** A dated, commit-pinned assessment produced by a coding agent, in
> the spirit of riprap's own `review/` convention. It is **not** a plan and nothing in it is
> scheduled. It records what was true at the commit above so that a later reader can tell
> whether a finding still holds. Findings are split into *verified here* and *reported by
> survey* — the second group is well-cited but was not independently re-run, and should be
> confirmed before anyone acts on it.

## Summary

The migration plan's phase list tracked **structural** adoption — vendoring, `src/` layout,
subpackage split — and that work is genuinely complete. It never tracked **coverage**. Treating
"only Phase 5 (container) remains" as the state of the migration is therefore misleading, and
this document exists to replace that framing.

At this commit `rqm/` specifies three modules out of roughly twenty: `pipeline/geometry_3d.py`,
`pipeline/heteroatom_assignment.py`, and `pipeline/opls_typing.py`. The frequently-quoted
"41 of 41 scenarios covered" is a coverage statistic *within those three documents*. Across the
test suite, 18 of 27 files carry zero `rq-` references.

The encouraging half of that finding: by git archaeology, the three specified modules are
exactly the three with the worst historical fix density (`opls_typing` 46% fix commits,
`constants` 45%, `geometry_3d` three threshold retunes inside 48 hours). The spine demonstrably
goes where the pain is. It has simply been pointed at a fifth of the code.

Separately, several riprap capabilities were never evaluated — not deferred with a reason, but
absent from the plan altogether. Those are listed in section 3.

---

## 1. Coverage: which modules should get a requirements document next

Ranked by fix history × current exposure × absence of specification. This ordering deliberately
differs from the not-yet-specified list at `rqm/ARCHITECTURE.md:134-136`, which ranks `valence`
second — `valence` is the best-tested module in the set — and omits the orchestrator entirely.

### 1.1 `export/gromacs_export.py` — first

The `.top` file delegates one hundred percent of the physics to the names it asserts. That is
the same argument `rqm/opls-typing.md` already makes and wins, extended one file downstream.

*Verified here:*

- `nrexcl` is written as a bare literal `3` at `src/biochar/export/gromacs_export.py:197` and
  again at `:338` (`f"{molecule_name:20s} 3\n\n"`). No named constant, no rationale comment.
- No improper dihedrals are emitted anywhere in the module (`grep -nE 'nrexcl|improper'` over
  the file returns nothing). For planar aromatic sheets under OPLS-AA, impropers are what keep
  sp2 centres from pyramidalizing during MD.
- `grep -rnE 'nrexcl|improper|cgnr|qtot' tests/` returns nothing — none of these fields is
  asserted anywhere in the suite.

*Reported by survey, not re-verified:* dihedrals are deduplicated by `tuple(sorted(...))` at
`:279` and `:412`, which discards atom ordering and collapses two genuinely distinct dihedrals
over the same four atoms — routine in fused rings — into one. `cgnr` is `1` for every atom
(`:213`, `:354`). Atom names diverge between `.gro` (truncated, `:140`) and `.top` (not
truncated, `:211`, `:352`) above 10,000 atoms.

*Historical precedent.* The worst silent-wrongness this project has shipped is in this file:
the `.gro` coordinate unit error (Å written where nm was required, `CHANGELOG.md:167`), which
meant every simulation before the fix ran a 10×-oversized box. `tests/test_gromacs_export.py:74`
still cannot detect a recurrence — it asserts `abs(x_nm) < 5.0`, and benzene at ~1.4 Å satisfies
that in either unit. The Å→nm conversion is asserted for *box sizes* at `:362` and never for
atom coordinates.

*Also relevant:* `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md` names
its own gap — dihedrals are emitted with `funct` only and are never checked against the
forcefield, the same failure class as the angles that produced eight fix commits. That gap has
no requirement.

### 1.2 `export/md_setup.py` — second, and writing it will surface bugs

Roughly 360 lines of the module are literal GROMACS `.mdp` physics. `tests/test_md_setup.py`
(22 tests) makes no assertion about any `.mdp` value — no test opens `dry_em.mdp`,
`anneal_npt.mdp`, or any of the seven templates.

*Verified here:*

- The local driver script calls `gmx solvate ... -p "$SIM/wet.top"` at
  `src/biochar/export/md_setup.py:716`, **before** anything creates `wet.top`; line `:717`
  writes `wet.top.base`. Under `set -euo pipefail` this path cannot pass the solvation stage.
  The SLURM path does it correctly — `cp "$in_top" ./wet.top` at `:852`.
- Every `grompp` invocation passes `-maxwarn 2` (`:722`, `:740`, `:745`, `:747`, `:858`, and
  others). That suppresses exactly the two science-fatal warnings: net charge with PME, and
  `.gro`/`.top` atom-name mismatch — the latter being the divergence described in 1.1.

*Reported by survey, not re-verified:* two tests assert the dead `wet.top.base` line
(`tests/test_md_setup.py:170`, `:177`). Sequential `genion` calls share one stale `.tpr`
(`:722-730`), so each rebuilds from pre-ion coordinates. The annealing peak temperature is
hardcoded at 1000 K for every structure (`:204-205`, `:244-245`) while `workflows/condensation.py:57-61`
scales it 1000/2000/3000 K by HTT — two implementations of the same Wood et al. protocol that
disagree. The two modules also use **opposite** semi-isotropic compressibility conventions
(`md_setup.py:459-465` vs `condensation.py:670`); swapping either produces a run that completes
and reports a plausible density while being physically wrong.

### 1.3 `pipeline/biochar_generator.py` config resolution — third

The densest source of user surprise in the codebase: where a caller's request quietly becomes
something else. All *reported by survey*:

- Three-level composition precedence — explicit beats `(temperature, feedstock)`-derived beats
  hard defaults `0.5 / 0.1 / 90.0` (`:264-307`). Undiscoverable at the call site.
- Aromaticity is silently clamped **up** to `MIN_BUILDABLE_AROMATICITY = 70.0`
  (`:292-299`, constant at `constants.py:667`), making it a floor rather than a prediction and
  discarding the data-derived value.
- `strict` raises only when `placed_counts == 0` (`:745-754`). Requesting ten carboxyls and
  receiving one passes strict validation.
- `charge_method="ml"` combined with `pH` is rejected with a five-line rationale (`:250-262`),
  but the guard lives here rather than in `charges/ml_charges.py` — calling the refiner directly
  on an ion silently neutralises it.

### 1.4 `workflows/surface_builder.py` — fourth

Zero requirement coverage; `tests/test_surface_builder.py` carries no `rq-` refs and appears
nowhere in `rqm/registry.json`. The module concentrates the material `rqm/` exists to hold: a
physically-motivated constant driving geometry (`pore_diameter + 3.4 Å`), an identity/deep-copy
optimisation with a non-obvious correctness condition, and a GROMACS external contract split
across two files.

*Reported by survey, not re-verified — treat as a probable defect:* when amorphous packing fails
and falls back to a slit stack (`:360-372`), `config.pore_type` is never updated, and the `.gro`
title is selected from it at `:431-439`. A fallen-back system is therefore written to disk
titled `"Amorphous surface"` while actually being a slit stack. Neither the returned
`SheetResult`s nor `box_vectors` record the fallback; the only trace is a Python warning the
caller may have suppressed. `tests/test_surface_builder.py:527-538` exercises the path and never
reads the title.

Two pieces of this module already read as requirement prose and are the model for the rest: the
`_sheets_identical` docstring (`:458-471`) and the pH-seed comment (`:488-490`). The pH tests at
`tests/test_ph_integration.py:341-422` are strong enough to carry `rq-` IDs today with no new
test work.

### 1.5 Remaining, in order

`pipeline/carbon_skeleton.py` — five rewrites in history, two H/C fixes on the same day
(`887739a`, `408d064`), and `target_aromaticity` is accepted and explicitly ignored (`:832-833`,
"Unused (kept for backward compatibility)"). Listed first in the project's own not-yet-specified
list. Then `workflows/sweep.py` (what `status="fallback"` promises; retry-vs-fallback semantics,
named in the clash-threshold solution doc as the thing that masked a systematic failure), then
`models/temperature_model.py` (see 4.3), then a single **silent-degradation** document covering
`amorphous_fallback`, `on_error="skip"`, and the warning-only config guards — the standing
principle at `rqm/ARCHITECTURE.md:122` ("Requested composition is a target, not a promise")
currently has no scenario.

### 1.6 The cheapest traceability win

The 5-character GROMACS residue-name limit. Six existing tests across four files already defend
it — `tests/test_gromacs_export.py:109`, `tests/test_sweep.py:60`, `tests/test_cli.py:111`,
`tests/test_ph_integration.py:364`, `tests/test_generator.py:852`,
`tests/test_surface_builder.py:37` — and it is already stated in `CLAUDE.md` and at
`rqm/ARCHITECTURE.md:117`. Cost: one scenario and six comment lines.

One further high-redundancy behaviour worth a single named requirement:
charge-sum-equals-formal-charge, defended at fifteen-plus sites using **four different
tolerances** (1e-4, 1e-5, 1e-6, and `test_generator.py:1118`'s 1e-5) with no stated rationale
for the spread. Writing that requirement forces the tolerance question nobody has answered.

> **Correction, same day.** An earlier revision of this section also listed the
> carbonyl/quinone/lactone → phenolic fallback as an unspecified cheap win. That is wrong. It is
> already fully specified at `rqm/heteroatom-assignment.md:27-62` under `rq-a9492b86`, with four
> scenarios (`rq-d3ef932d`, `rq-7d3e9432`, `rq-2f5ab3f8`, `rq-fbb6a9c0`) all back-annotated to
> `tests/test_heteroatom_assignment_gaps.py`. The claim came from a survey pass that searched
> `CLAUDE.md` and the tests but not `rqm/`, and it survived into the first draft because the
> contradicting pass was not reconciled against it. Recorded rather than silently deleted,
> because it is a concrete instance of this document's own stated limitation: *reported by
> survey* claims were not independently re-executed, and at least one of them was false.

---

## 2. Solution docs whose lesson has no requirement

`rqm/` and `docs/solutions/` are complements — one states the invariant, the other records how
it was learned. Three lessons currently exist only on the retrospective side.

From `docs/solutions/bugs/physical-features-misread-as-geometry-errors.md`:

1. **Force-field refinement must be unconditional on the non-hex-lattice path.** Refinement sat
   behind `if steric_clashes:` and only ever ran because false clashes were always present;
   fixing the clash rule would have silently disabled it. The code comment survives at
   `src/biochar/pipeline/biochar_generator.py:876`; no scenario in `rqm/geometry-embedding.md`
   mentions it. One scenario would guard a defect that would otherwise re-appear invisibly.
2. **Seed retry is a remedy for variance only.** `sweep.build_point`'s `max_retries` → fallback
   masked a failure that was deterministic in O/C.
3. **`StructureValidator` truncates to `geom_errors[:3]`,** which hid the bond-length defect
   entirely. No requirement covers error-report completeness.

From `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md`: the dihedral gap
described in 1.1, plus the absence of an end-to-end `grompp` smoke test — depth 3 proves a table
entry exists, not that GROMACS starts.

`docs/solutions/conventions/xfail-strict-as-deferred-gap-tripwire.md` has no requirement either,
but that is correct — it is a process convention about test scaffolding, not product behaviour.
It is also, on inspection, the single healthiest convention in the repo: both `xfail(strict=True)`
sites (`tests/test_constants_ff.py:112`, `tests/test_opls_type_map.py:352`) are currently empty
dicts, meaning every exemption has been retired by XPASS exactly as designed.

---

## 3. Riprap capabilities never evaluated

Distinct from the container layer, which the plan defers with a stated reason. These were never
named in the plan at all.

**Three of riprap's four agent skills.** Riprap ships `rr-architecture` (writes
`rqm/ARCHITECTURE.md`, intended to run first), `rr-plan`, `rr-implement` (code from requirements,
with the traceability discipline built into the loop), and `rr-quiz` (comprehension quiz on your
own code). The plan names `rr-plan` only, and only as the directory `rqm.sh` was copied out of.
Given the `compound-engineering` plugin already covers plan/work/review, `rr-implement` is the
interesting one — it is the mechanism that would keep requirements moving with code
automatically rather than by the `AGENTS.md` rule alone.

**Upstream's `tests/test_rqm.sh` — 584 lines, one test per Gherkin scenario, testing the very
script this repo vendored.** Not copied. `tools/rqm/ids.md:155` points readers at it; that path
does not exist here. The result is an unverified 739-line bash dependency at the centre of a
repository whose thesis is traceability.

**The `## Feature API` entity type.** Riprap's `bse.md` demonstrates it; none of this repo's four
`rqm/*.md` files has such a section (`grep` returns zero hits), and `rqm/registry.json` holds 4
`file`, 29 `section`, and 41 `scenario` entries with no API entries. Half the tool is unexercised.
Not necessarily wrong — but it is a choice nobody has consciously made.

**The `ownership-exceptions` manifest and "Riprap-managed" marker convention.** A
Copier-independent way to mark vendored files as not locally owned. `tools/rqm/README.md` does
this in prose only.

**The `review/` convention itself** — dated, commit-pinned analyses. This document is the first
instance.

**Windows support and the CI secret scan.** Riprap ships `.ps1`/`.bat` throughout and CI-tests on
`windows-latest`; this repo's tooling is bash/macOS/Linux and CI is ubuntu-only. Riprap's CI
template also runs `check-secrets.sh --repository` as a job step — see 4.4.

**No drift to fix.** `tools/rqm/rqm.sh` and `tools/rqm/ids.md` were confirmed byte-for-byte
identical to upstream `main` (MD5 `0489ba6d1515d32883eeef27b208e31c` on both sides). Upstream's
last change to `rqm.sh` predates the vendoring. Nothing needs re-syncing.

---

## 4. Adopted parts that have gone stale

All four verified at this commit.

### 4.1 The `rq check` CI gate is still warn-only, and its exit criterion is met

`.github/workflows/ci.yml:103` runs `bash tools/rqm/rq check || true`, and the surrounding
comment still says 16 scenarios have no defending test. That was true at Phase 1; it is now
zero. The migration plan's own exit criterion — promote `check` to a hard failure once the
unreferenced count reaches zero — is satisfied.

One caveat for whoever does it: `rq check` currently emits 30 warnings at this commit, but every
one is a `type: section` or `type: file` heading (`rq-12b97d18` "## Public surface",
`rq-18307b81` "## Layering", and so on) that no test could reasonably reference. All 41
`type: scenario` entries are referenced. A hard failure must be scoped to scenarios or the job
fails immediately on structural headings.

### 4.2 `tools/rqm/README.md` describes the pre-Phase-3 layout

Lines 15-18 say the wrapper sets `SRC_DIR=biochar` "for this repo's flat package layout" and
that calling `rqm.sh` directly "would scan a nonexistent `src/`". `tools/rqm/rq:19` sets
`SRC_DIR=src`, and has since Phase 3. The wrapper's own header comment is correct and current;
only the README was missed. The stated failure mode is now inverted — the real reason to go
through the wrapper is the `cd "$REPO_ROOT"` at `:16`, which makes the tool work from a
subdirectory.

### 4.3 The CI base-branch filter the solution doc told us to remove

`.github/workflows/ci.yml:6-7` still reads `pull_request: branches: [main]`. That is the exact
configuration `docs/solutions/workflow-issues/stacked-pr-ci-and-base-branch-traps.md` documents
as the cause of stacked PRs running zero checks while reporting `MERGEABLE / CLEAN`. The doc's
own words: with the filter in place, every future stacked PR re-earns the trap and the only
defence is remembering to look. The compensating control that shipped instead is the last line
of `AGENTS.md` ("Do not stack PRs") — a memory-dependent rule, which is precisely the mechanism
the sibling convention doc argues always rots.

This is a one-line deletion whose rationale is already written down. It has the best
effort-to-value ratio in this review.

### 4.4 `check-secrets.sh` is a rewrite, not riprap's, and has real gaps

`tools/hooks/check-secrets.sh:32` matches `AKIA`, `gh[pousr]_`, `xox`, `AIza`, and JWTs. It does
**not** match `sk-ant-`. In a repository developed with Claude Code, an Anthropic API key pasted
into a file passes the gate. There is also no `--repository` mode, so nothing in CI can invoke
it — the scan exists only in a hook that `git commit --no-verify` bypasses.

Riprap's version is complementary rather than superior: it has path-based rules (rejecting
`.env`, `.claude/.credentials.json`, `.codex/auth.json` by filename regardless of content) and
`--repository`, but misses AWS, Slack, Google, and JWT patterns. Neither is a superset of the
other.

---

## 5. Findings worth acting on independent of traceability

Each verified at this commit unless noted.

**`tests/test_pah_quality.py` contributes zero collected tests.** Its only class is
`PAHQualityTester` (`:76`) and `pyproject.toml:88` sets `python_classes = ["Test*"]`, so none of
its seven `test_*` methods is collected. `CLAUDE.md` documents this as intentional ("not pytest
— run directly"), and that is a legitimate choice — but `grep` for `pah_quality` across
`.github/workflows/` and `tools/hooks/pre-commit` returns nothing, so it runs in neither CI nor
the commit gate. *Reported by survey:* the file contains no `assert` statements at all; its
thresholds are computed and printed, and several bodies swallow exceptions and continue. Net
effect: the quality bar for stage 1 of the pipeline is enforced nowhere automatic.

**Two `importorskip` calls guard a first-party module.** `tests/test_condensation.py:206` and
`:342` call `pytest.importorskip("biochar.pipeline.biochar_generator")`. If that module ever
fails to import — syntax error, rename, broken transitive dependency — both turn
green-as-skipped rather than red, and they are the only coverage of the `generate_and_condense`
path. Other files import the same module unconditionally
(`tests/test_temperature_model.py:16`, `tests/test_ml_charges.py:17`), so the guard buys nothing.
It should be deleted.

**`TemperatureModel.get_valid_range` is property-agnostic by signature.** It takes only
`feedstock` (`src/biochar/models/temperature_model.py:403-404`) and returns the first hit from
`_OVERRIDE_PROPS = ("H_C_ratio", "O_C_ratio")` — in practice always H/C's 100-1000 °C range. The
generator's out-of-range temperature warning therefore consults the wrong property's support:
`temperature=1000` warns about nothing even though `pH` (200-900) and `ec_dsm` (220-900) are
100 °C beyond their data. *Reported by survey:* the one test that exercises this uses
`temperature=50`, below every range, so the mismatch is invisible
(`tests/test_generator.py:1039-1044`).

**Both changelogs are stale.** `CHANGELOG.md:5`'s last entry is 0.4.0 (2026-07-11) while roughly
twenty PRs of behaviour change have landed since, including pH protonation, non-zero net charge,
four corrected wrong-element type mappings, and three geometry-threshold changes.
`docs/changelog.rst:4` stops at 0.2.0 while `pyproject.toml:7` declares 0.4.0. `CHANGELOG.md`'s
footer still reads "Current Version: 0.1.4 / Last Updated: May 31, 2026", and four of the seven
unchecked "Future Enhancements" boxes have shipped.

**`README.md:665` claims amorphous porous surfaces are not yet implemented.** They shipped in
0.1.3 (`8df0dee`) and gained `amorphous_fallback` in 0.3.0 (`069b86c`).

**`.claude/settings.json` still carries the `PreToolUse` commit hook.** `CLAUDE.md` correctly
describes it as not the gate, but it remains live configuration invoking the full test suite
under a 600s timeout the suite cannot meet.

---

## 6. Suggested order, if this is picked up

Nothing here is scheduled. Presented as a recommendation only.

1. **Retire the two stale-config items** — drop `branches: [main]` from the `pull_request`
   trigger (4.3), and promote `rq check` to a hard failure scoped to `type: scenario` (4.1).
   Both are one-line changes with the reasoning already written down.
2. **`rqm/gromacs-export.md`** (1.1) — pairs with the document that already owns "only the
   `opls_XXX` name reaches GROMACS", and covers the worst historical failure class in the repo.
   Should absorb the dihedral gap the OPLS solution doc names in section 2.
3. **`rqm/md-setup.md`** (1.2) — writing it surfaces the `wet.top` ordering bug and the
   `-maxwarn 2` suppression as scenarios rather than guesses.
4. **The cheap wins** (1.6) — the 5-character residue limit, the phenolic fallback, and the
   charge-sum tolerance question.
5. Then the config-resolution and surface-stacking documents (1.3, 1.4).

Independently of all of the above, section 5's four items are small and self-contained.

## Method and limitations

Produced by five parallel agent surveys over the working tree at `2ad77a1` — module audit,
riprap upstream inventory, test-suite audit, git/PR archaeology, and a follow-up on `constants.py`
and `surface_builder.py`. Claims marked *verified here* were re-run directly against the commit
above. Claims marked *reported by survey* carry citations but were not independently re-executed;
confirm before acting.

This is a snapshot. Line numbers will drift. Where a finding matters, the file and symbol names
are given alongside the line so it stays locatable.
