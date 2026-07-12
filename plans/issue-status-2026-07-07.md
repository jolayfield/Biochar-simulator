# Biochar Simulator — Issue Status & Outstanding Work Plan

**Generated:** 2026-07-07 (automated)
**Version:** 0.3.0 (unchanged since last report)
**Open GitHub issues:** 0 (all #1–#5 closed)

---

## Summary

No open GitHub issues — consistent with the 2026-06-26 report. However, a
substantial chunk of new work has accumulated in the working tree since then
and is still entirely uncommitted: a parameter-sweep driver, a GROMACS MD
setup pipeline, and a PFAS-ligand integration layer (~2,200 lines across 5
new modules + 2 test files). None of it is on GitHub as an issue or a commit
yet. This report proposes 9 issues to track it to done, mirroring how the
2026-06-09 backlog (ISSUE-A..K) was tracked to the 2026-06-26 release.

---

## Changes Since 2026-06-26 Status Report

**New, untracked, and functional (has passing tests per `PFAS_HANDOFF.md`):**
- `biochar/sweep.py` / `biochar/sweep_cli.py` — factorial parameter-sweep
  driver (`biochar-sweep template|run`), documented in README, tests in
  `tests/test_sweep.py`.
- `biochar/md_setup.py` / `biochar/md_setup_cli.py` — GROMACS MD pipeline
  script generator (dry anneal → solvate → ion → wet-equilibrate), tests in
  `tests/test_md_setup.py`.

**New, untracked, and not yet test-covered:**
- `biochar/pfas_ligands.py` — PFOA⁻/PFOS⁻/PFBS⁻ LigParGen build-script
  rendering + biochar/ligand topology merge. No `tests/test_pfas_ligands.py`
  exists yet.

**New, untracked, supporting files:**
- `PFAS_HANDOFF.md` (design/session handoff notes — not code)
- `examples/sweeps/pfas_temperature_grid.yaml`, `examples/sweeps/oxygen_group_grid.yaml`
- `sweep_out/pfas_temperature_grid/` (a generated sweep run — build output)

**Modified but uncommitted:**
- `README.md` — new "Parameter sweeps" + MD-setup documentation sections.
- `biochar/__init__.py` — exports for `sweep`, `md_setup`, `pfas_ligands`.
- `pyproject.toml` — two new console-script entry points
  (`biochar-sweep`, `biochar-md-setup`).
- `.claude/settings.json` — plugin enablement (tooling, not package code).

**Known pre-existing, unrelated test failures** (per last session's full-suite
run, not re-verified this run since `rdkit`/`biochar` aren't installed in this
shell): 18 failures in `tests/test_ml_charges.py` — `ModuleNotFoundError: No
module named 'sklearn'`, an environment gap, not a code defect.

---

## Proposed GitHub Issues

### Issue A — Commit the sweep / MD-setup / PFAS-ligand feature work

**Label:** `chore` | **Effort:** 30–45 min

`biochar/sweep.py`, `sweep_cli.py`, `md_setup.py`, `md_setup_cli.py`,
`pfas_ligands.py`, both new test files, the README additions, the
`__init__.py` exports, and the two new `pyproject.toml` script entries are
all sitting as uncommitted/untracked working-tree changes. Nothing is on a
branch or PR. Suggest 2–3 commits split by feature (sweep; md_setup; PFAS
glue) rather than one giant commit, so `git bisect`/review stays meaningful.

**Touchpoints:** all files listed in "New, untracked" and "Modified but
uncommitted" above.

---

### Issue B — No unit tests for `biochar/pfas_ligands.py`

**Label:** `test` | **Effort:** 2–3 hr

Add `tests/test_pfas_ligands.py` covering `render_ligpargen_molecules_txt`,
`render_ligpargen_build_script`, `_split_biochar_top`, and
`merge_biochar_pfas_topology` — happy path plus the documented
`PFASLigandError` cases (missing `atomtypes.itp`, missing per-species
`.itp`/`.gro`, duplicate species, `n_copies < 1`, malformed biochar `.top`).

**Touchpoints:** `biochar/pfas_ligands.py`, new `tests/test_pfas_ligands.py`.

---

### Issue C — Extend `tests/test_md_setup.py` for the `cfg.ligands` code path

**Label:** `test` | **Effort:** 1–2 hr

`MDSetupConfig.ligands` / `.ligand_system_dir` and the ligand-insertion
staging in `_render_pipeline_script`, `_render_solvate_ions_slurm`, and
`_render_hpc_chain_script` currently have no test coverage — only the
pre-existing non-PFAS paths are exercised.

**Touchpoints:** `biochar/md_setup.py`, `tests/test_md_setup.py`.

---

### Issue D — No local end-to-end smoke test for the PFAS pipeline

**Label:** `test` | **Effort:** ~1 hr

Fabricate a minimal fake `ligand_system_dir` (`atomtypes.itp` +
`PFOA.itp`/`PFOA.gro`) and run
`setup_one_structure(..., config=MDSetupConfig(ligands=...))` against the
already-generated `sweep_out/pfas_temperature_grid/structures/000_T300_softwood/`
structure to confirm `merged.top`, `run_pipeline.sh`, and (with
`cluster="hpc"`) `solvate_ions.slurm`/`submit_chain.sh` all come out
well-formed — without touching `gmx`/`ligpargen`/hpc.

**Touchpoints:** new test/script under `tests/` or `examples/`.

---

### Issue E — Expose `ligands`/`ligand_system_dir` on the CLI

**Label:** `enhancement` | **Effort:** 2–3 hr

`md_setup_cli.py` and `sweep_cli.py` don't expose the new PFAS-ligand fields
— the feature is Python-API-only right now. Add the equivalent CLI flags
(species list + copy counts + ligand system dir) once Issues B–D land.

**Touchpoints:** `biochar/md_setup_cli.py`, `biochar/sweep_cli.py`.

---

### Issue F — Add a top-level `examples/pfas_binding/` walkthrough

**Label:** `docs` | **Effort:** 1–2 hr

No example currently ties the full workflow together end-to-end: generate a
biochar structure → render the LigParGen build script for PFOA/PFOS/PFBS →
merge topology → build the MD pipeline. Add an example directory
demonstrating the full 3-species × structure × ion-profile matrix described
in `PFAS_HANDOFF.md` §4.

**Touchpoints:** new `examples/pfas_binding/` directory.

---

### Issue G — Stale README claim: "Amorphous porous surfaces not yet implemented"

**Label:** `docs` | **Effort:** 10 min

`README.md`'s "Known Limitations" table (~line 562) still lists "Amorphous
porous surfaces not yet implemented ... reserved for a future release." This
is stale: `SurfaceConfig(pore_type="amorphous")` plus `amorphous_fallback`
have been implemented and documented since v0.3.0 (`069b86c`), and GitHub
issue #1 ("Amorphous porous packing (Phase 2)") has been closed since
2026-06-01. This row should be removed or rewritten to describe the current
(fully implemented) behavior.

**Touchpoints:** `README.md` "Known Limitations" section.

---

### Issue H — `sweep_out/`-style run directories aren't covered by `.gitignore`

**Label:** `chore` | **Effort:** 10 min

`.gitignore` already excludes `output/` and `examples/output/`, and excludes
GROMACS file extensions individually, but a sweep run's `manifest.csv`,
`manifest.json`, and any QC plots (e.g. `sweep_out/.../qc_oxygen_trend.png`)
aren't covered by any existing pattern, so the whole `sweep_out/` directory
shows up as untracked bulk noise in `git status`. Add `sweep_out/` (or
whatever the sweep driver's documented default output directory convention
is) to `.gitignore`.

**Touchpoints:** `.gitignore`.

---

### Issue I — Review/decide on the `.claude/settings.json` plugin change

**Label:** `chore` | **Effort:** 5 min (decision only)

`.claude/settings.json` has an uncommitted diff enabling the
`compound-engineering` plugin. This is tooling configuration rather than
package code, and unrelated to the sweep/md_setup/PFAS work above — flagging
it separately so it doesn't get silently swept into an unrelated commit (or
silently dropped) when Issue A is executed.

**Touchpoints:** `.claude/settings.json`.

---

## Future candidates (not urgent, carried over from 2026-06-26)

- **LBCC bond-charge corrections** — `qm_charges.py` notes 1.14*CM1A-LBCC
  corrections aren't applied. Low priority; base 1.14*CM1A already matches
  LigParGen's default.
- **N-doped amorphous structures** — combining `pore_type="amorphous"` with
  pyridinic/graphitic N substitution is untested; may surface ring-parity
  edge cases in `carbon_skeleton.py`.
- **hpc `scratch_root` misconfiguration** (from `PFAS_HANDOFF.md` §5) —
  `/scratch/<user>` should likely be `/home0/<user>` or
  `/extra0/<user>` in the Compute panel. Ops/config only; matters once
  someone is ready to actually execute the PFAS pipeline on hpc, not before.
- **`tests/test_ml_charges.py` sklearn gap** — 18 failures from a missing
  `scikit-learn` dependency in the `biochar-md` conda env; environment issue,
  not a code defect, so not proposed as a numbered issue above.
