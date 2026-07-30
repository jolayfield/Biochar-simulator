# refactor: map Biochar-simulator onto the MolSSI riprap scaffold

**Status:** phases 0-3 complete; phase 4's one functional item (the surface_builder cycle) done,
its subpackage-split remainder optional; phase 5 deferred
**Type:** refactor
**Date:** 2026-07-28

> **Progress.** Phase 0 and 1 merged in PR #34, phase 2 in PR #35, phase 3 in PR #37. PR #36 (not
> in the original plan) fixed the embedder-retry defect that profiling turned up while doing the
> phase 2 follow-up work. PR #38 (also not in the original plan) closed the 16-scenario `rqm/`
> coverage gap Phase 1 left open — `bash tools/rqm/rq check` now reports 41 of 41 scenarios
> covered. PR #40 broke the `biochar_generator` <-> `surface_builder` import cycle, the one item
> in Phase 4 with a real functional payoff. Remaining: the rest of Phase 4 (subpackage split across
> five 1,000+-line modules, optional and `rqm/`-guided) and phase 5 (container, deferred). Neither
> is blocking.

## Summary

Riprap (https://github.com/MolSSI/riprap, MIT, Taylor Barnes) is a Copier template that
scaffolds a scientific software project around four separable ideas: a requirements
traceability system (`rqm/`), a hardened Python project skeleton, a containerized agent
sandbox, and a managed/user file-ownership split that lets `copier update` push upstream
improvements into your repo.

This repo already has better *content* than a fresh riprap project — 730 tests, real CI, a
Sphinx docs tree, and a `docs/solutions/` learning corpus. What it lacks is riprap's
**spine**: a machine-checkable link from "this is the behavior we promised" to "this is the
test that defends it." That single missing piece is the highest-value thing riprap has to
offer here, and it can be adopted without moving a line of code.

**Recommendation: vendor, don't `copier copy`.** Take the `rqm/` traceability system whole,
take the project-skeleton fixes selectively, defer the container, and skip the
managed/user ownership split. Reasoning in "Why not full Copier adoption" below.

---

## What riprap actually prescribes

Four layers, in descending order of value to this repo:

**1. Requirements traceability (`rqm/`) — adopt fully.**
Feature-level markdown docs under `rqm/`, each containing prose requirements plus Gherkin
scenarios. Every entity — file, `##`/`###` section, Feature-API bullet, and `Scenario:`
block — carries a stable opaque ID of the form `rq-XXXXXXXX` (8 lowercase hex, assigned
once, never re-derived). Headings and API bullets embed it as a trailing HTML comment;
Gherkin scenarios use an `@rq-XXXXXXXX` tag. Source and test files reference an ID in a
comment line. A ~750-line standalone bash script, `rqm.sh` (needs only `grep`/`find`/`sed`/
`jq`/`od`), maintains `rqm/registry.json`, which maps each ID to its type, file, title,
declaration line, and the list of files referencing it. Four subcommands: `stamp` (assign
missing IDs), `index` (rebuild registry, fail on duplicates), `check` (stale refs are
errors, unreferenced requirements are warnings), `clean` (drop dead entries).

**2. Python project skeleton — adopt selectively.**
`src/` layout, `py.typed`, `devtools/conda-envs/test_env.yaml`, `.codecov.yml`,
`.gitattributes`, `.github/CONTRIBUTING.md`, `.github/PULL_REQUEST_TEMPLATE.md`, a CodeQL
workflow, and a real git `pre-commit` hook plus a `check-secrets` hook.

**3. Agent-neutral configuration — adopt cheaply.**
Riprap generates `CLAUDE.md` *and* `AGENTS.md`, and ships each skill in four parallel
locations (`.claude/skills/`, `.agents/skills/`, `.opencode/skills/`, and the canonical
`.riprap/managed/skills/`) so the project is not locked to one agent vendor.

**4. Containerized agent sandbox — defer.**
Rootless Podman, `rr.sh`/`rr.bat` launchers, per-agent credential isolation in named
volumes, Apptainer `.sif` export for clusters. This is the bulk of `.riprap/managed/` and
is orthogonal to maintainability.

---

## Structural mapping

### Where this repo already satisfies riprap

| Riprap element | Biochar-simulator | Verdict |
|---|---|---|
| `rqm/ARCHITECTURE.md` | `CLAUDE.md` "Architecture" + `CONCEPTS.md` | Content exists, wrong home |
| Sphinx `docs/` + ReadTheDocs | `docs/` + `.readthedocs.yaml` | Already better than the template |
| `tests/` + CI matrix | 25 files, ~730 tests, 4-version matrix | Already better |
| Packaging metadata | `pyproject.toml` with entry points, extras | Good |
| Agent permission defaults | `.claude/settings.json` allow/ask/deny | Present, Claude-only |
| (no riprap equivalent) | `docs/solutions/` learning corpus | Keep — riprap has nothing like it |

### Where the gap is real

| Riprap element | Biochar-simulator | Gap |
|---|---|---|
| `rqm/*.md` + `registry.json` | absent | **No link from behavior spec to test** |
| `src/{{package}}/` layout | flat `biochar/` | Tests can pass against the working tree, not the installed package |
| `py.typed` + type hints | absent | No typing contract for downstream users |
| Lint config + CI lint job | `.ruff_cache/` exists, no `[tool.ruff]`, no CI step | Ruff was run once by hand and never wired up |
| `devtools/conda-envs/test_env.yaml` | absent | The **recommended** install path (conda-forge) is never tested |
| Real git `pre-commit` hook | `.claude/settings.json` `PreToolUse` hook | **Gate only fires when Claude commits** |
| `check-secrets` hook | absent | — |
| `.codecov.yml` | absent, but CI uploads to Codecov | Upload with no config |
| `.gitattributes` | absent | — |
| `CONTRIBUTING.md`, PR template, CodeQL | absent | — |
| `AGENTS.md` | absent | Agent config is Claude-only |
| Root kept minimal | 8 root `.md` files | Doc sprawl |

---

## The three findings worth acting on first

These came out of the audit and stand on their own, independent of whether you adopt
riprap at all.

**F1 — The pre-commit test gate is bypassable.** `CLAUDE.md` states that a `PreToolUse`
hook "runs the full suite before any `git commit` and blocks the commit if it fails." That
is true only for commits *Claude* makes. A human running `git commit` in a terminal never
touches the hook. Riprap puts this in `.riprap/managed/hooks/pre-commit`, a real git hook,
which catches both. The current arrangement gives a false sense of a gate.

**F2 — CI never resolves the declared dependencies, and never tests the recommended
install path.** `.github/workflows/ci.yml` installs with a hand-written
`pip install rdkit networkx numpy scipy scikit-learn` followed by
`pip install -e . --no-deps`. The version floors declared in `pyproject.toml`
(`numpy>=1.24.0`, `rdkit>=2023.9.1`, `scipy>=1.10.0`, `networkx>=3.1`) are therefore never
exercised — a bad floor would ship undetected. Separately, `CLAUDE.md` and the README
recommend `conda install -c conda-forge biochar`, and no CI leg installs that way.

**F3 — `biochar_generator` and `surface_builder` are mutually dependent.**
`surface_builder.py:24` imports `BiocharGenerator` at module scope;
`biochar_generator.py:1363` imports `SurfaceBuilder` *inside a function* to break the
resulting cycle. The deferred import works, but it marks a layering violation:
`generate_surface` is a surface-level entry point living in the single-molecule module.

Everything else in the tree is in decent shape. Of 121 tracked files, the simulation
output directories (`sim/`, `sim_200*/`, `output/`, `sweep_out/`, `build/`, `data/`) are
correctly untracked and gitignored — the working directory looks messy but the repository
is not.

---

## Why not full Copier adoption

Running `copier copy gh:MolSSI/riprap .` over this repo is technically possible —
`copier.yml` marks 50+ paths skip-if-exists, so `README.md`, `LICENSE`, and the CI
workflows would survive. It is still the wrong trade:

- **The bulk of what Copier would manage is the container layer**, which is the layer being
  deferred. You would take on `copier update` as an ongoing obligation to receive updates
  to `.riprap/managed/container/` and `.riprap/managed/launch/` — files that would sit
  unused.
- **`.claude/settings.json` would collide.** Riprap owns that file
  (`settings.json.jinja`); this repo has a custom `PreToolUse` hook and a deny-list in it.
  Riprap's escape hatch is `.local.json` overrides, which means restructuring working
  config to fit the template's ownership model.
- **The `src/` layout divergence would fight every update.** Until the layout move in Phase
  3 lands, template-generated paths would not match the tree.
- **The piece with real value is one MIT-licensed bash script and a documented ID
  convention.** Vendoring `rqm.sh` and `ids.md` into `tools/rqm/`, with attribution, costs
  nothing and carries no update obligation.

Revisit this if the container layer later becomes attractive — the MOPAC and GROMACS
external-binary dependencies are a genuine reproducibility argument for it.

---

## Migration plan

### Phase 0 — Clear the decks (cheap; the hard part is already done)

Phase 3 rewrites every path in the repo, so any unmerged branch would be left nearly
unmergeable. **Verified 2026-07-28: there is no unmerged work.** Every remote branch —
`feat/heptagon-curvature`, `feat/qm-cm1a-charges`, `worktree-tidy-enchanting-valley` — is
`ahead=0` relative to `main` and is a stale ref only. The pH-protonation work landed via
PR #25 (`reconcile/ph-protonation`, `c11787f`), which is also where the `SS`-as-carbon
regression was corrected and the carboxyl/phenolic typing split (`1c85e33`) was kept; its
plan doc is marked `status: completed`. This is the cleanest possible moment for a layout
move — the window is open, and it closes again as soon as long-lived branches appear.

1. **Delete `test_valence_comprehensive.py` and `test_valence_update.py`** from the repo
   root. `CLAUDE.md` already documents them as "older scratch tests, not part of the
   suite"; they are outside `testpaths` and run by nobody.
2. **Prune the three merged remote branches** so the branch list stops implying there is
   outstanding work.
3. **Do Phase 3 before opening any long-lived branch**, or accept re-doing this check.

*Exit criterion:* stale refs pruned, scratch tests gone, `git status` clean.

---

### Phase 1 — Traceability spine (highest value, zero code movement)

This is the phase that justifies the whole exercise. The repo has ~730 tests and no way to
answer "which promised behavior does this test defend?" or "if I change this threshold,
what breaks?" Back-annotating existing tests is far cheaper than writing them was.

1. **Vendor the tooling.** Copy `rqm.sh` and `ids.md` from
   `.riprap/managed/skills/rr-plan/` into `tools/rqm/`, preserving the MIT notice and
   crediting MolSSI/Taylor Barnes. Note that `rqm.sh` defaults `SRC_DIR=src/` — set
   `SRC_DIR=biochar` until Phase 3 lands. Add `jq` to the dev environment.
2. **Write `rqm/ARCHITECTURE.md`.** Lift the pipeline description from `CLAUDE.md` and the
   vocabulary from `CONCEPTS.md`. `CLAUDE.md` then shrinks to agent instructions and points
   at `rqm/ARCHITECTURE.md`, matching riprap's convention.
3. **Write feature docs, most-churned modules first.** Do not attempt all thirteen at once.
   Start where rework has actually concentrated:
   - `rqm/geometry-embedding.md` — clash thresholds, the H-bond floor, the hex-lattice
     skip, the clash-severity tolerance. This area has been retuned repeatedly.
   - `rqm/heteroatom-assignment.md` — O/C placement, aliphatic hydroxyls, and the
     **silent** carbonyl/quinone/lactone fallback to phenolic. That fallback is documented
     in `CLAUDE.md` prose and defended by no named scenario; it is exactly what a Gherkin
     scenario should pin down.
   - `rqm/opls-typing.md` — the nearest-analog force-field parameters (`CA-S-CA`,
     `NPY-CA-OH`, `SS`, `NGR`) that are deliberately not QM-validated. Recording that as an
     explicit, ID-bearing requirement is strictly better than the current arrangement,
     where the decision lives in a memory note.

   Then, as each area is next touched: `carbon-skeleton`, `valence-validation`
   (from `VALENCE_SYSTEM.md`), `gromacs-export` (from `GROMACS_WORKFLOW.md`),
   `parameter-sweep` (from `BATCH_GENERATION_GUIDE.md`), `qm-charge-backend` (from
   `docs/qm-charge-backend.md`), `surface-stacking`, `md-setup`, `condensation-annealing`,
   `ph-protonation`, `temperature-feedstock-model`.
4. **Stamp and index.** `rqm.sh stamp` assigns IDs; `rqm.sh index` builds
   `rqm/registry.json`.
5. **Back-annotate tests.** Add `# rq-xxxxxxxx` comments to the existing tests that already
   defend each scenario. This is the step that pays: `rqm.sh check` then reports every
   requirement with no test behind it, which is a real coverage map rather than a line-count
   one.
6. **Wire `rqm.sh check` into CI as warn-only.** Promote to a hard failure once the
   unreferenced-requirement count reaches zero for the documented modules.

*Exit criterion:* `rqm.sh check` runs green in CI; the three priority modules have zero
unreferenced scenarios.

**Relationship to `docs/solutions/`:** these are complements, not duplicates. Keep both.
`docs/solutions/` is retrospective — "we hit this bug, here is what we learned."
`rqm/` is prospective — "this is the behavior we promise, here is the test that proves it."
A recurring pattern in this repo is a fix landing in `docs/solutions/` with no durable
statement of the invariant it established; `rqm/` is where that invariant goes.

---

### Phase 2 — CI and hygiene (independent of Phase 1, can run in parallel)

1. **Install a real git `pre-commit` hook** (fixes F1). Port riprap's
   `.riprap/managed/hooks/pre-commit` and `check-secrets.sh`. Keep the Claude `PreToolUse`
   hook — belt and braces — but the git hook becomes the authority. Update `CLAUDE.md`,
   which currently overstates the existing gate's coverage.
2. **Fix the CI dependency install** (fixes F2). Replace the hand-written pip line with a
   plain `pip install -e ".[dev,ml]"` so declared floors are actually resolved. Add
   `devtools/conda-envs/test_env.yaml` and one CI leg that installs via conda-forge, so the
   recommended path is tested.
3. **Add ruff.** A `[tool.ruff]` block in `pyproject.toml` and a CI lint job. `.ruff_cache/`
   shows ruff has been run by hand; wire it up or delete the cache.
4. **Add the community and security files:** `.gitattributes`, `.codecov.yml` (CI already
   uploads coverage with no config), `.github/CONTRIBUTING.md`,
   `.github/PULL_REQUEST_TEMPLATE.md`, and a CodeQL workflow.
5. **Add `AGENTS.md`.** Riprap's own is short — coding-standards pointers plus a reference
   to `rqm/ARCHITECTURE.md`. Makes the project legible to Codex and OpenCode, not just
   Claude.

*Exit criterion:* a human `git commit` is gated by the suite; CI resolves declared
dependencies and exercises the conda path; lint runs on every PR.

---

### Phase 3 — Layout and doc consolidation

1. **Move `biochar/` to `src/biochar/`.** Set `[tool.setuptools.packages.find]
   where = ["src"]`. Test imports (`from biochar.x import y`) are unchanged under an
   editable install, so this is lower-risk than it looks — but it rewrites every path in
   the repo, which is why Phase 0 comes first. Payoff: tests can no longer accidentally
   import the working tree instead of the installed package, which is the specific failure
   mode that lets a packaging bug reach users with CI green.
2. **Add `py.typed`** and begin annotating the public API surface declared in
   `__init__.py.__all__`.
3. **Consolidate root markdown into `docs/`.** `BATCH_GENERATION_GUIDE.md`,
   `BEST_PRACTICES.md`, `GROMACS_WORKFLOW.md`, `VALENCE_SYSTEM.md`, and
   `RELEASE_SUMMARY.md` move under `docs/`; their normative content graduates to `rqm/`.
   Root keeps `README.md`, `CHANGELOG.md`, `LICENSE`, `CLAUDE.md`, `AGENTS.md`, and
   `CONCEPTS.md`.
4. **Retire `plans/issue-status-*.md`.** Six dated status snapshots are a symptom of having
   no durable requirements registry. Once `rqm/registry.json` exists, `rqm.sh check` answers
   what those snapshots were approximating. Keep the design plans
   (`wood-condensation-pipeline.md`, `temperature_model_plan.md`,
   `split-pfas-workflow.md`); fold them into `docs/plans/`.
5. **Set `SRC_DIR` back to the default** in the vendored `rqm.sh` invocation.

*Exit criterion:* `pip install .` from a clean checkout, then the full suite, passes against
the installed package.

---

### Phase 4 — Subpackage decomposition (optional, `rqm/`-guided)

Only worth doing after Phase 1, because the `rqm/` feature docs are what tell you where the
seams actually are. Five modules exceed 1,000 lines (`biochar_generator.py` at 1,400,
`geometry_3d.py` at 1,211, `carbon_skeleton.py` at 1,173, `md_setup.py` at 1,131,
`heteroatom_assignment.py` at 1,112) against riprap's "favor composable designs" guidance.

The import graph is already an almost-clean DAG, which makes this tractable:

- **Layer 0** (no intra-package imports): `constants`, `valence`, `qm_charges`
- **Layer 1** (→ `constants`): `carbon_skeleton`, `geometry_3d`, `opls_typing`,
  `temperature_model`
- **Layer 2**: `heteroatom_assignment`, `gromacs_export`, `ml_charges`, `protonation`,
  `validation`
- **Layer 3**: `biochar_generator` (imports 12 modules)
- **Layer 4**: `surface_builder`, `condensation`, `sweep`, `md_setup`
- **Layer 5**: the four CLI modules

Proposed subpackages, following those layers: `biochar/pipeline/` (the five generation
steps), `biochar/charges/` (`qm_charges`, `ml_charges`), `biochar/export/`
(`gromacs_export`, `md_setup`), `biochar/workflows/` (`sweep`, `condensation`,
`surface_builder`), `biochar/models/` (`temperature_model`), `biochar/cli/`.

**Break the `biochar_generator` ↔ `surface_builder` cycle** (fixes F3) by moving
`generate_surface` out of `biochar_generator.py` into the surface module, retiring the
deferred import at `biochar_generator.py:1363`. `__init__.py` re-exports it either way, so
the public API is unchanged.

> **Done.** Merged 2026-07-30 as PR #40 (`1065060`) — the one item in Phase 4 with a real
> functional payoff, split out from the line-count tidying that makes up the rest of the
> phase. The rest of Phase 4 (subpackage split across five 1,000+-line modules) remains
> optional and not started.

Sequencing note: `__init__.py.__all__` is the public API contract and `pyproject.toml`
`[project.scripts]` declares the entry points. Both must keep working unchanged through
this phase — subpackage moves stay invisible to users.

---

### Phase 5 — Container (deferred; revisit, do not schedule)

Riprap's rootless-Podman sandbox has a genuine argument here that it lacks for a pure-Python
project: this package shells out to an external MOPAC binary, needs GROMACS force-field data
at test time, and recommends a conda install. A container would pin all three. But it is
orthogonal to maintainability, it is the heaviest piece to adopt, and Phases 1–3 deliver
more per unit effort. Revisit after they land — and if the answer is yes, that is the point
to reconsider full Copier adoption, since the container is what `copier update` would
actually be maintaining.

---

## Effort and value

| Phase | Value | Effort | Risk |
|---|---|---|---|
| 0 — Clear the decks | Prerequisite | Very low | None |
| 1 — `rqm/` traceability | **Highest** | Medium, incremental | None (no code moves) |
| 2 — CI and hygiene | High | Low | Low |
| 3 — Layout | Medium | Medium | Medium (touches every path) |
| 4 — Subpackages | Medium | High | Medium |
| 5 — Container | Low, for now | High | — |

Phases 1 and 2 are independent and can proceed in parallel. Phase 1 alone delivers most of
the long-term maintainability gain, because the durable problem in this repo is not code
organization — it is that hard-won invariants (clash thresholds, force-field analog choices,
silent fallbacks) live in prose, memory notes, and retrospective solution docs rather than in
anything a test run can check.

---

## Attribution

Riprap is MIT-licensed, Copyright (c) 2026 Taylor Barnes, MolSSI.
Any vendored `rqm.sh` / `ids.md` must retain that notice.
