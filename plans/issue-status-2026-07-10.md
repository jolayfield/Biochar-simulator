# Biochar Simulator — Issue Status & Outstanding Work Plan

**Generated:** 2026-07-10 (automated)
**Version:** 0.3.0 (unchanged since last report)
**Open GitHub issues:** 0 (all #1–#5 closed)

---

## Summary

No open GitHub issues — consistent with the 2026-06-26 and 2026-07-07 reports.
The entire Issues A–I backlog proposed on 2026-07-07 (sweep driver, GROMACS MD
setup pipeline, PFAS-ligand integration layer) remains **completely unactioned**:
same files untracked/modified, nothing committed, no GitHub issues filed. On top
of that, a new chunk of work has landed since 07-07 and is also uncommitted: a
three-stage fix for the H/C-ratio shortfall in generated structures
(`biochar/carbon_skeleton.py`, `biochar/heteroatom_assignment.py`,
`biochar/biochar_generator.py`, `tests/test_hc_control.py`). This report carries
the prior backlog forward unchanged and adds one new issue (J) for the H/C fix.

---

## Re-verification of the 2026-07-07 backlog

Spot-checked the two cheapest, most "should be quick to just fix" items to
confirm the backlog is genuinely stale, not just unfiled:

- **Issue G** (stale README claim, `README.md:592`, "Amorphous porous surfaces
  not yet implemented") — **still present, unfixed.**
- **Issue H** (`sweep_out/` missing from `.gitignore`) — **still missing**;
  `.gitignore` only excludes `output/` / `examples/output/`.

`git diff --stat` confirms the same modified-file set as 07-07
(`.claude/settings.json`, `README.md`, `biochar/__init__.py`, `pyproject.toml`)
plus the same untracked new modules/tests (`biochar/sweep.py`, `sweep_cli.py`,
`md_setup.py`, `md_setup_cli.py`, `pfas_ligands.py`, `PFAS_HANDOFF.md`,
`examples/sweeps/`, `sweep_out/`, `tests/test_md_setup.py`, `tests/test_sweep.py`).
Nothing here has changed in three days — no commits landed against any of it.

---

## New since 2026-07-07: H/C-shortfall fix

**New, untracked/modified, and test-covered:**
- `biochar/carbon_skeleton.py` (+226/-? lines) — `_grow_graph` gained a
  `compactness` param; a two-pass compact-vs-elongated build
  (`_build_pure_hex_compact` / `_build_pure_hex_elongated`) picks whichever
  edge fraction lands closest to the requested H/C target.
- `biochar/heteroatom_assignment.py` (+67 lines) — `assign_hydrogens` now logs
  a warning and sets `CompositionResult.h_c_ceiling` /
  `.h_c_target_unreachable` when a target exceeds the structural ceiling
  (previously it silently capped). New `attach_aliphatic_carbons` adds pendant
  sp3 methyl groups (OPLS types `opls_135`/`opls_140`) on a smaller aromatic
  core so total C stays on target while H/C climbs to 0.6–0.8.
- `biochar/biochar_generator.py` (+88 lines) — new `GeneratorConfig.allow_aliphatic`
  flag (default `True`); set `False` or `aromaticity_percent >= 99` to force a
  pure-aromatic structure (old behavior).
- `tests/test_hc_control.py` (new) — covers all three stages; per project
  notes, targets 0.4–0.8 now land within 10% tolerance at 30/50/100 C. Known
  remaining limit: an H/C *floor* (~0.44 at 50 C for a maximally condensed
  small flake) — not a regression, just a documented boundary.

**Combined diff:** 3 files, +337/-44 lines, plus one new test file. All
uncommitted.

---

## Proposed GitHub Issues

### Issues A–I (carried forward from 2026-07-07, unchanged — still open)

- **A** — Commit the sweep / MD-setup / PFAS-ligand feature work
  (`chore`, 30–45 min) — `biochar/sweep.py`, `sweep_cli.py`, `md_setup.py`,
  `md_setup_cli.py`, `pfas_ligands.py`, both new test files, README additions,
  `__init__.py` exports, two new `pyproject.toml` entry points. Suggest 2–3
  commits split by feature (sweep; md_setup; PFAS glue).
- **B** — No unit tests for `biochar/pfas_ligands.py` (`test`, 2–3 hr) — add
  `tests/test_pfas_ligands.py` covering `render_ligpargen_molecules_txt`,
  `render_ligpargen_build_script`, `_split_biochar_top`,
  `merge_biochar_pfas_topology`, happy path + `PFASLigandError` cases.
- **C** — Extend `tests/test_md_setup.py` for the `cfg.ligands` code path
  (`test`, 1–2 hr) — `MDSetupConfig.ligands`/`.ligand_system_dir` and the
  ligand-insertion staging have no coverage yet.
- **D** — No local end-to-end smoke test for the PFAS pipeline (`test`, ~1 hr)
  — fabricate a fake `ligand_system_dir` and confirm `merged.top`,
  `run_pipeline.sh`, `solvate_ions.slurm`/`submit_chain.sh` come out
  well-formed without touching `gmx`/`ligpargen`/hpc.
- **E** — Expose `ligands`/`ligand_system_dir` on the CLI (`enhancement`,
  2–3 hr) — `md_setup_cli.py`/`sweep_cli.py` are Python-API-only right now.
- **F** — Add a top-level `examples/pfas_binding/` walkthrough (`docs`,
  1–2 hr) — tie generate → LigParGen build script → merge → MD pipeline
  together end-to-end (3-species × structure × ion-profile matrix from
  `PFAS_HANDOFF.md` §4).
- **G** — Stale README claim: "Amorphous porous surfaces not yet implemented"
  (`docs`, 10 min) — `README.md` ~line 592 still contradicts
  `SurfaceConfig(pore_type="amorphous")`/`amorphous_fallback`, implemented and
  documented since v0.3.0 (`069b86c`); issue #1 closed since 2026-06-01.
- **H** — `sweep_out/`-style run directories not covered by `.gitignore`
  (`chore`, 10 min) — add `sweep_out/` (or the sweep driver's documented
  default output convention) to `.gitignore`.
- **I** — Review/decide on the `.claude/settings.json` plugin change (`chore`,
  5 min decision only) — uncommitted diff enabling the `compound-engineering`
  plugin; tooling config, unrelated to A–H; flag separately so it isn't
  silently swept into or dropped from an unrelated commit.

### Issue J — Commit and document the H/C-shortfall fix (new this report)

**Label:** `bug` (fix) + `docs` | **Effort:** 45–60 min

The three-stage H/C fix described above is complete and passing its own test
suite (`tests/test_hc_control.py`), but:

1. **Not committed** — recommend a commit separate from Issue A's
   sweep/PFAS-glue commits, since this is a correctness fix to the core
   generation path rather than new opt-in feature code.
2. **Not documented** — README has no mention of `GeneratorConfig.allow_aliphatic`,
   the new `CompositionResult.h_c_ceiling`/`.h_c_target_unreachable` fields, or
   the known H/C-floor limitation (~0.44 at 50 C for maximally condensed small
   flakes).
3. **No CHANGELOG entry** — last version bump (`bf9a6ed`, v0.3.0) predates this
   fix; consider whether a v0.3.1 patch release is warranted given it corrects
   a composition-accuracy defect on the main (non-opt-in) generation path.

**Touchpoints:** `biochar/carbon_skeleton.py`, `biochar/heteroatom_assignment.py`,
`biochar/biochar_generator.py`, `tests/test_hc_control.py`, `README.md`,
`CHANGELOG.md`, `pyproject.toml` (if a release is warranted).

---

## Recommended execution order

1. **Issue J first** — a correctness fix to generation internals that the
   sweep/PFAS work (Issues B–F) builds on top of; committing and releasing it
   before the sweep/PFAS commits keeps `git bisect`/review meaningful.
2. **A, G, H, I** — quick, mostly mechanical, no interdependencies.
3. **B → C → D → E → F** in order — each builds on the previous landing.

---

## Future candidates (not urgent, carried over from 2026-06-26 / 2026-07-07)

- **LBCC bond-charge corrections** — `qm_charges.py` notes 1.14*CM1A-LBCC
  corrections aren't applied. Low priority; base 1.14*CM1A already matches
  LigParGen's default.
- **N-doped amorphous structures** — combining `pore_type="amorphous"` with
  pyridinic/graphitic N substitution is untested; may surface ring-parity
  edge cases in `carbon_skeleton.py` (now also touched by the H/C fix above —
  worth re-checking for interaction once Issue J lands).
- **hpc `scratch_root` misconfiguration** (from `PFAS_HANDOFF.md` §5) —
  `/scratch/<user>` should likely be `/home0/<user>` or
  `/extra0/<user>` in the Compute panel. Ops/config only; matters once
  someone is ready to actually execute the PFAS pipeline on hpc, not before.
- **`tests/test_ml_charges.py` sklearn gap** — 18 failures from a missing
  `scikit-learn` dependency in the `biochar-md` conda env; environment issue,
  not a code defect, so not proposed as a numbered issue above.
