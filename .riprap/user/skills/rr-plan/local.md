# rr-plan — local overrides for Biochar-simulator

`SKILL.md` is vendored byte-for-byte from MolSSI/riprap so it can be re-synced by
re-downloading. Everything this repository needs to change about it lives here. The skill reads
this file last, and these instructions take precedence where they conflict.

## Tooling paths

Upstream assumes a Copier-managed `.riprap/managed/` tree. This project vendored only the
traceability tooling, so every path the skill names is different here.

| SKILL.md says | Use instead |
|---|---|
| `.riprap/managed/skills/rr-plan/rqm.sh stamp <path>` | `bash tools/rqm/rq stamp <path>` |
| `.riprap/managed/skills/rr-plan/rqm.sh index` | `bash tools/rqm/rq index` |
| `.riprap/managed/skills/rr-plan/rqm.sh show rq-XXXXXXXX` | `bash tools/rqm/rq show rq-XXXXXXXX` |
| `.riprap/managed/skills/rr-plan/bse.md` | `tools/rqm/bse.md` |

Always go through the `rq` wrapper, never `rqm.sh` directly — it `cd`s to the repository root,
without which a run from a subdirectory scans nothing and reports every requirement as
unreferenced.

## Orientation

Read `AGENTS.md` first, then `rqm/ARCHITECTURE.md`. `CLAUDE.md` exists but points at `AGENTS.md`
for everything except Claude-specific configuration.

The language is **Python**. `SKILL.md`'s Feature API example and `tools/rqm/bse.md` are both
Rust; the structure they demonstrate applies unchanged, the syntax does not.

Before writing, check `docs/solutions/` for a documented past problem in the same area.
`rqm/` and `docs/solutions/` are complements: one states the invariant, the other records how it
was learned. A requirement that contradicts a solution doc means one of them is wrong — resolve
that before drafting.

## Scenarios land with their tests

CI fails on any `type: scenario` entry with no defending test. A requirements file therefore
**cannot merge on its own** — the PR that adds a scenario must also add the test that references
it. This is stricter than upstream and it is deliberate: an unreferenced scenario is the exact
gap this system exists to surface, and letting it sit in `main` unmonitored defeats the purpose.

When the scenario describes behaviour the code does not yet have — which is normal, since
specifying a defect correctly means writing what *should* be true — the test still lands in the
same PR, marked `pytest.mark.xfail(strict=True)` with a comment naming the defect it defers. See
`docs/solutions/conventions/xfail-strict-as-deferred-gap-tripwire.md`. `strict=True` matters: the
fix retires the exemption by XPASS rather than leaving a stale marker behind.

Write scenarios so that exactly one test can be constructed per scenario. A scenario needing
three tests is really three scenarios.

## Feature API sections

Whether new files carry a `## Feature API` section is an open decision (D1 in
`docs/plans/2026-07-30-001-docs-rqm-coverage-expansion-plan.md`). No existing `rqm/` document has
one, and `rqm/registry.json` currently holds zero API entries. Until that is settled, add one
only for a document covering interfaces consumed across subpackage boundaries — the `export/`
writers and the `workflows/` builders — and not for `pipeline/` internals.

## Domain constraints worth knowing before drafting

These are easy to get wrong and expensive to get wrong:

- Requested composition is a **target, not a promise**. Measure the realised structure. A
  scenario asserting that a request is honoured exactly is usually wrong.
- `molecule_name` and residue names are **5 characters or fewer** — a GROMACS `.gro` limit.
- Only the `opls_XXX` name reaches GROMACS. A name identifying the wrong element silently
  simulates the wrong chemistry and nothing complains.
- Clash warnings on hex-lattice structures are artefacts of a geometrically exact lattice, not
  defects.
- A green test run is not automatically evidence: force-field-backed tests skip silently without
  a discoverable `oplsaa.ff`. Check the skip count — 3 is normal, ~21 means those tests verified
  nothing.

## Pull requests

Branch from `main`. **Do not stack** — historically a PR based on another branch ran zero checks
while reporting `MERGEABLE / CLEAN`. The `pull_request` branch filter that caused it is gone, but
the one-PR-per-change habit stands.
