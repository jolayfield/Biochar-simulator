# rqm — requirements traceability tooling

`rqm.sh`, `ids.md`, `bse.md`, and `tests/test_rqm.sh` are vendored verbatim from
[MolSSI/riprap](https://github.com/MolSSI/riprap), MIT licensed,
Copyright (c) 2026 Taylor Barnes. The licence text is in `LICENSE-riprap`.
The companion `rr-plan` skill is vendored to `.claude/skills/rr-plan/SKILL.md`,
with this repository's overrides in `.riprap/user/skills/rr-plan/local.md`.

They are vendored rather than adopted via `copier` because the traceability
system is the only part of the riprap template this project uses; taking a
`copier update` dependency would also make us track riprap's container and
launcher layers, which we do not use. See
`docs/plans/2026-07-28-001-refactor-riprap-scaffold-migration-plan.md`.

## Use

Go through the `rq` wrapper rather than calling `rqm.sh` directly. Since the
src-layout move its `RQM_DIR` / `SRC_DIR` / `TESTS_DIR` values match `rqm.sh`'s
own defaults, so what the wrapper actually buys you is the `cd` to the
repository root — without it, running from a subdirectory scans nothing and
silently reports every requirement as unreferenced.

```bash
tools/rqm/rq stamp [file] # assign IDs to entities that lack one
tools/rqm/rq index        # rebuild rqm/registry.json; fails on duplicate IDs
tools/rqm/rq check        # stale refs are errors; unreferenced requirements are warnings
tools/rqm/rq clean        # drop registry entries for deleted files / missing IDs
tools/rqm/rq show ID      # print one requirement and everything referencing it
```

Requires `jq` and bash >= 4.0. macOS `/bin/bash` is 3.2 — `brew install bash`.

## Workflow

1. Write or edit a feature document under `rqm/`. Draft with **no** `rq-`
   annotations — not even placeholders like `@rq-PENDING`, which stop `stamp`
   from assigning a real ID and then have to be cleaned up by hand.
2. `tools/rqm/rq stamp rqm/<file>.md` to give new headings, API bullets, and
   scenarios their IDs.
3. Add `# rq-xxxxxxxx` comments to the tests that defend each scenario.
4. `tools/rqm/rq index` then `tools/rqm/rq check`.

`check` reporting a **scenario** as unreferenced means no test claims to defend
it. That is the signal this system exists to produce — treat it as a coverage
gap, not as noise, and note that CI fails on it. `check` also warns about
unreferenced `section` and `file` entries; those are document headings that no
test can reasonably cite, and CI ignores them.

## Verifying the tooling

`tests/test_rqm.sh` is upstream's own suite, one test per Gherkin scenario:

```bash
bash tools/rqm/tests/test_rqm.sh
```

**One test fails on macOS, and the failure is real.**
`t_fix_duplicates_resolves_copy` fails with `sed: -I or -i may not be used with
stdin`. The cause is the single `sed -i` in the file, at `rqm.sh:271` inside
`_fix_duplicates()`: GNU `sed -i` takes no argument, BSD/macOS `sed -i` requires
one. So **`stamp --fix-duplicates` does not work on macOS.**

The blast radius is contained. That is the only `sed -i` in `rqm.sh`, and all
nine `stamp` tests pass — the ordinary `rq stamp` path this project depends on
is unaffected, on macOS as elsewhere. `--fix-duplicates` is only reachable after
`index` has already refused a duplicate ID, which is rare. When it happens on
macOS, edit one of the two colliding IDs by hand and re-run `index`. Installing
GNU sed ahead of BSD sed on `PATH` should also work, but has not been tested
here.

Not patched locally, on purpose — see below.

## Local modifications

None. `rqm.sh`, `ids.md`, `bse.md`, and `tests/test_rqm.sh` are byte-for-byte
upstream so they can be re-synced by re-downloading, and the macOS `sed`
defect above is deliberately left in place rather than patched: a local fix
would silently be reverted by the next re-sync, and the bug belongs upstream.
All local configuration lives in the `rq` wrapper and, for the skill, in
`.riprap/user/skills/rr-plan/local.md`.

Note that `ids.md` and the vendored `SKILL.md` both describe the upstream
install path (`.riprap/managed/skills/rr-plan/`) and Rust examples. The
conventions they document apply unchanged here; only the paths differ, and
`local.md` states the mapping.
