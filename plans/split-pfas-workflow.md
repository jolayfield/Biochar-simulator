# Plan: Split biochar-simulator (generic) from the PFAS sorption workflow

**Goal:** keep `biochar-simulator` an independent tool that only generates and
equilibrates biochar surfaces, and move all PFAS-specific work into a separate
project that consumes `biochar-simulator` as a library/driver.

**Key simplifier:** all the affected work is still **uncommitted** (untracked
files + uncommitted edits), so there is no git history to surgically preserve.
The split is just deciding where each file lands and doing one refactor before
anything is committed.

---

## 1. Target end-state

Two independent repos with a one-way dependency (`biochar-pfas` → `biochar`):

```
┌─────────────────────────────┐         ┌──────────────────────────────────┐
│  biochar-simulator  (this)  │◄────────│  biochar-pfas  (new repo)        │
│  a library + CLIs           │ imports │  an application / driver         │
│                             │         │                                  │
│  • generate biochar surfaces│         │  • PFAS species (PFOA/PFOS/PFBS) │
│  • parameter sweeps         │         │  • LigParGen build glue          │
│  • MD equilibration pipeline│         │  • biochar+ligand topology merge │
│    (dry anneal→solvate→ion  │         │  • orchestration: sweep surfaces │
│    →wet equilibrate)        │         │    → equilibrate → insert PFAS   │
│  • GENERIC injection seam   │         │  • consumes the injection seam   │
└─────────────────────────────┘         └──────────────────────────────────┘
   knows nothing about PFAS               depends on `biochar>=X.Y`
```

**Boundary rule of thumb:** anything that names a chemical species, ligand,
LigParGen, or PFAS is application code → new repo. Anything that generates or
equilibrates a bare biochar surface is generic → stays here.

**Coupling today (the thing being broken):** `md_setup.py` imports
`biochar.pfas_ligands` (line 39) and threads ligand params through its
renderers. The dependency is one-directional (`md_setup` → `pfas_ligands`);
`pfas_ligands.py` itself imports only stdlib, so it moves out cleanly.

---

## 2. Where every uncommitted file goes

| File | Destination | Notes |
|---|---|---|
| `biochar/sweep.py`, `sweep_cli.py`, `tests/test_sweep.py` | **stays** (biochar) | Generic factorial generation; not PFAS-specific |
| `biochar/md_setup.py`, `md_setup_cli.py`, `tests/test_md_setup.py` | **stays**, but **refactored** | Base equilibration pipeline is generic; PFAS ligand hooks extracted (§3) |
| `IonProfile` / `ION_PROFILES` (in `md_setup.py`) | **stays** | Generic solvation water chemistry (rename `mn_calcareous_default` if you want it un-branded) |
| `examples/sweeps/oxygen_group_grid.yaml` | **stays** | Generic example |
| `examples/sweeps/pfas_temperature_grid.yaml` | **new repo** | PFAS-named config; keep a generic `temperature_grid.yaml` here instead |
| `biochar/pfas_ligands.py` | **new repo** | Self-contained (stdlib-only imports); moves cleanly |
| `PFAS_HANDOFF.md` | **new repo** | Application design notes |
| `sweep_out/` | **neither** | Build output → `.gitignore` here; regenerate in the new repo as needed |
| `plans/issue-status-2026-07-07.md` | **stays** (or delete) | Status report, not code |
| `.claude/settings.json` (+3) | decide separately | Tooling; don't sweep into a feature commit |
| `biochar/__init__.py` (+50), `pyproject.toml` (+2) | **split** | Generic exports/entry-points stay; PFAS exports move (§4) |
| `README.md` (+177) | **split** | "Parameter sweeps" + "GROMACS run setup" stay (minus PFAS ligand bits); PFAS workflow docs move |
| *(missing)* `tests/test_pfas_ligands.py` | **write, in new repo** | The untested merge logic gets its coverage there |

---

## 3. The one real refactor: a generic injection seam in `md_setup.py`

Today `md_setup.py` knows about ligands directly. Replace that PFAS knowledge
with a **domain-neutral extension point** so biochar-simulator stays generic but
the PFAS repo can still splice its stage into the middle of the pipeline
(ligands must go in *after* dry-anneal, *before* solvation).

**Chosen seam shape: A1 — declarative data object.** The caller declares *what*
to insert; biochar renders the `gmx insert-molecules` commands itself (inserting
molecules into a box is a generic MD operation, not PFAS-specific). biochar owns
its own path/shell conventions, so the caller can't break by guessing them; the
object is serializable (drops into a manifest) and trivially unit-testable.

**In biochar-simulator**, add two neutral dataclasses and thread a
`pre_solvation_stage` through instead of
`ligands`/`ligand_system_dir`/`ligand_gro_names`/`merged_top_name`:

```python
@dataclass
class MoleculeInsertion:
    gro_file: str            # a .gro filename staged into the run dir
    n_copies: int
    n_try: int = 500         # gmx insert-molecules -try

@dataclass
class PreSolvationStage:
    """Extra molecule insertions spliced in after dry-anneal, before solvation."""
    name: str                              # banner label in the generated script
    insertions: list[MoleculeInsertion]    # biochar renders the insert-molecules chain
    solvation_top: str                     # topology to use from solvation onward
    extra_files: list[str] = field(default_factory=list)  # files copied into the run dir
```

- `MDSetupConfig` drops `ligands` / `ligand_system_dir`; gains
  `pre_solvation_stage: Optional[PreSolvationStage] = None`.
- `_render_pipeline_script` / `_render_solvate_ions_slurm` /
  `_render_slurm_chain_script` take the neutral `pre_solvation_stage`, render the
  `gmx insert-molecules` chain from `insertions` (using biochar's own path
  exprs), copy `extra_files`, and switch the topology to `solvation_top` from
  solvation onward — exactly the current logic, minus the words "ligand"/"PFAS".
- **Delete** `from biochar.pfas_ligands import ...` (line 39),
  `_render_ligand_insertion_lines`, and the `merge_biochar_pfas_topology` call
  in `setup_one_structure`.

**In biochar-pfas**, `pfas_ligands.py` grows the piece that was in `md_setup`:
- `merge_biochar_pfas_topology(...)` (unchanged),
- a new `build_pre_solvation_stage(placements, ligand_system_dir, merged) ->
  PreSolvationStage` that maps PFAS placements to a list of `MoleculeInsertion`s
  and returns biochar's neutral `PreSolvationStage`.

Net effect: biochar's equilibration pipeline is fully generic and independently
testable; the PFAS repo owns 100% of the ligand knowledge and passes it in
through one typed seam.

---

## 4. The interface contract (what biochar-pfas imports)

```python
from biochar.sweep     import run_sweep, load_sweep_config           # generate surfaces
from biochar.md_setup  import (setup_md_from_manifest, MDSetupConfig,  # equilibrate
                               PreSolvationStage, IonProfile)
```

That's the entire public surface the PFAS repo needs. biochar-simulator imports
*nothing* from biochar-pfas.

---

## 5. Execution phases

**Phase 0 — prep (this repo, small):**
1. `echo "sweep_out/" >> .gitignore`; remove `sweep_out/` from the tree (regenerable).
2. Decide names: new repo (`biochar-pfas`? `pfas-sorption`?) and whether to
   rebrand `mn_calcareous_default`.

**Phase 1 — land the generic layers here (2 PRs):**
3. **Sweep PR** — `sweep.py` + `sweep_cli.py` + `tests/test_sweep.py` + generic
   example + the `__init__`/`pyproject`/README hunks for sweep only (per-hunk
   staging of the shared files, like the H/C README).
4. **MD-setup PR** — refactor `md_setup.py` to `PreSolvationStage` (§3), drop the
   PFAS import, keep `IonProfile`; add `md_setup_cli.py` + `tests/test_md_setup.py`
   (now PFAS-free) + its `__init__`/`pyproject`/README hunks. This PR is where
   the decoupling actually happens.

After Phase 1, biochar-simulator is a clean, independent, fully-tested
surface-generation + equilibration tool with a generic extension seam.
**Tag/release it** (e.g. v0.4.0) so the PFAS repo can pin a version.

**Phase 2 — scaffold biochar-pfas (new repo):**
5. `git init` new repo; `pyproject.toml` with `dependencies = ["biochar>=0.4.0"]`,
   a `pfas-*` console entry point.
6. Move in `pfas_ligands.py` (+ the extracted `render_ligand_insertion_stage`),
   `PFAS_HANDOFF.md`, `pfas_temperature_grid.yaml`.
7. **Write `tests/test_pfas_ligands.py`** — the merge/`PFASLigandError` cases
   (the current gap).
8. Add the orchestration driver: sweep →
   `setup_md_from_manifest(..., MDSetupConfig(pre_solvation_stage=render_ligand_insertion_stage(...)))`.

**Phase 3 — verify end-to-end:**
9. In biochar-pfas, run the PFOA/PFOS/PFBS × structure × ion-profile matrix
   through to generated (not executed) run dirs; confirm `merged.top` +
   `run_pipeline.sh` come out well-formed against a fabricated fake ligand
   system (no `gmx`/LigParGen needed — matches the "setup only, never run"
   constraint from `PFAS_HANDOFF.md`).

---

## 6. Decisions (confirmed)

1. **New repo name** — `biochar-pfas`. ✅
2. **Ion profiles** — **keep** `IonProfile`/`ION_PROFILES` (incl. the MN
   groundwater presets) in biochar-simulator as generic solvation config. ✅
3. **Seam shape** — **A1**, the declarative `PreSolvationStage` +
   `MoleculeInsertion` data object (§3). biochar renders the insert-molecules
   chain; no caller code runs inside biochar. ✅
4. **New repo location** — `~/Claude_Cowork/biochar-pfas`. ✅

## 7. Progress

- [x] **Phase 0** — `sweep_out/` added to `.gitignore` (files left on disk;
      regenerable). Names/decisions locked above.
- [x] **Phase 1 — Sweep PR** → PR #9 (merged). Generic sweep driver + CLI + 26
      tests + `temperature_grid.yaml` / `oxygen_group_grid.yaml` examples.
- [x] **Phase 1 — MD-setup PR** → PR #10 (`feat/md-setup-pipeline`, off `main`).
      Refactored `md_setup.py` to the generic `PreSolvationStage` +
      `MoleculeInsertion` seam; removed the `biochar.pfas_ligands` import and all
      ligand/PFAS coupling; kept `IonProfile`/`ION_PROFILES`. 21 tests (4 new for
      the seam). md_setup exports + `biochar-md-setup` entry.
- [x] **Phase 2 — scaffold `~/Claude_Cowork/biochar-pfas`** (git-init'd, initial
      commit `287ce4b`). `biochar_pfas/` package: `pfas_ligands.py` (+ new
      `build_pre_solvation_stage` → biochar `PreSolvationStage`),
      `orchestrate.py` (`setup_pfas_md`), `cli.py` (`biochar-pfas build-inputs`).
      `tests/test_pfas_ligands.py` — 16 tests (merge happy path + every
      `PFASLigandError` case + the adapter), all passing. Moved `PFAS_HANDOFF.md`
      + `pfas_temperature_grid.yaml`. `pyproject.toml` depends on `biochar>=0.4.0`.
- [x] **PR #10 merged** → biochar-simulator `main` (`50c848a`) is fully
      independent and PFAS-free.
- [x] **Phase 3 — end-to-end verify**: biochar-pfas `tests/test_orchestrate.py`
      (commit `d263e4c`) generates a real biochar structure + fake ligand system,
      runs `setup_pfas_md`, and asserts the run dir is well-formed (insertion
      before solvation, merged topology used onward). Setup-only; no gmx/MD.
      Full biochar-pfas suite: 17 passing.
- [x] **Finalize the split (biochar-simulator working tree)**: removed
      `biochar/pfas_ligands.py`, `PFAS_HANDOFF.md`,
      `examples/sweeps/pfas_temperature_grid.yaml`; cleaned the pfas import/exports
      out of the working-tree `__init__.py`. `main` itself was already PFAS-free.
- [ ] Optional local tidy: switch the local checkout off the merged
      `fix/hc-ratio-control` branch onto `main` and drop the redundant untracked
      copies of files now committed on `main` (blocked by the destructive-op
      permission gate this session — safe for the user to do directly).
- [ ] Tag biochar-simulator `v0.4.0` so biochar-pfas's `biochar>=0.4.0`
      dependency resolves against a real release (local editable installs work
      regardless for now).
