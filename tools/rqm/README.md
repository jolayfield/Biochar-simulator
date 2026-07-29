# rqm — requirements traceability tooling

`rqm.sh` and `ids.md` are vendored verbatim from
[MolSSI/riprap](https://github.com/MolSSI/riprap), MIT licensed,
Copyright (c) 2026 Taylor Barnes. The licence text is in `LICENSE-riprap`.

They are vendored rather than adopted via `copier` because the traceability
system is the only part of the riprap template this project uses; taking a
`copier update` dependency for two files would also make us track riprap's
container and launcher layers, which we do not use. See
`docs/plans/2026-07-28-001-refactor-riprap-scaffold-migration-plan.md`.

## Use

Always go through the `rq` wrapper — it sets `SRC_DIR=biochar` for this repo's
flat package layout. Calling `rqm.sh` directly would scan a nonexistent `src/`
and report every requirement as unreferenced.

```bash
tools/rqm/rq stamp    # assign IDs to entities that lack one
tools/rqm/rq index    # rebuild rqm/registry.json; fails on duplicate IDs
tools/rqm/rq check    # stale refs are errors; unreferenced requirements are warnings
tools/rqm/rq clean    # drop registry entries for deleted files / missing IDs
tools/rqm/rq show ID  # print one requirement and everything referencing it
```

Requires `jq` and bash >= 4.0. macOS `/bin/bash` is 3.2 — `brew install bash`.

## Workflow

1. Write or edit a feature document under `rqm/`.
2. `tools/rqm/rq stamp` to give new headings and scenarios their IDs.
3. Add `# rq-xxxxxxxx` comments to the tests that defend each scenario.
4. `tools/rqm/rq index` then `tools/rqm/rq check`.

`check` reporting a requirement as unreferenced means no test claims to defend
it. That is the signal this system exists to produce — treat it as a coverage
gap, not as noise.

## Local modifications

None. `rqm.sh` and `ids.md` are byte-for-byte upstream so they can be re-synced
by re-downloading. All local configuration lives in the `rq` wrapper. Note that
`ids.md` describes the upstream install path (`.riprap/managed/skills/rr-plan/`)
and Rust examples; the conventions it documents apply unchanged here.
