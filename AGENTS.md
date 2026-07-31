# AGENTS.md

Guidance for coding agents working in this repository. Vendor-neutral: Claude
Code reads `CLAUDE.md`, which points here for everything except its own
tool-specific configuration.

## Orient here first

- **`rqm/ARCHITECTURE.md`** — system design and the index of the requirements
  set. Read this before changing anything.
- **`rqm/*.md`** — feature requirements with Gherkin scenarios, each carrying a
  stable `rq-XXXXXXXX` ID. These are **prospective**: the behaviour we promise.
- **`docs/solutions/`** — documented solutions to past problems. **Retrospective**:
  what we learned, and why.
- **`CONCEPTS.md`** — domain vocabulary. The internal atom type / OPLS type name
  / bonded type distinction is load-bearing; read it before touching export.
- **`docs/reviews/`** — dated, commit-pinned assessments. **Point-in-time**: what
  was true at one commit. Nothing in a review is scheduled, and line numbers in
  one drift; re-verify a finding before acting on it.

`rqm/` and `docs/solutions/` are complements, not duplicates. One states the
invariant, the other records how it was learned.

## Working rules

**Requirements move with the code.** Changing behaviour in a specified area
means updating its `rqm/` document in the same change, then `bash tools/rqm/rq
index`. A requirement that no longer matches the code is worse than no
requirement.

**Run the fast tier constantly, the full suite at boundaries.**

```bash
pytest -m "not slow" -n auto   # ~24s
pytest -n auto                 # full suite
```

Always `-n auto`. 15 of 842 tests are ~98% of the runtime; they are marked
`slow` so the tier you run constantly stays usable.

**A green run is not automatically evidence.** Two silent-skip traps here:

- Forcefield-backed tests skip without a discoverable `oplsaa.ff`. **Check the
  skip count: 1 is normal (MOPAC), 21 means the forcefield was missing and those
  20 tests verified nothing.** The conda env already ships GROMACS, but nothing
  points the tests at it, so the default is the bad case:

  ```bash
  export GMXDATA="$CONDA_PREFIX/share/gromacs"   # note: parent of top/, not top/ itself
  ```

  `BIOCHAR_OPLSAA_FF` also works but points straight at the `oplsaa.ff`
  directory, and `tests/test_constants_ff.py` does not honour it — that file
  reads only `GMXDATA` and `gmx` on `PATH`. Prefer `GMXDATA`, which un-skips
  both files.
- `python3` without rdkit silently is not the project environment.

**Measure the structure, do not assume it.** A requested composition is a
target, not a promise: `carbonyl`, `quinone`, and `lactone` are silently placed
as phenolic, and oxygen can spill onto sp3 carbons to reach an O/C target.

**Prefer removing an unreachable branch to leaving it.** The `N`/`NT` atom types
sat unreachable-but-defined until a kekulization failure made them reachable
again, and they had no GROMACS mapping. Dead code in a lookup path is a latent
bug, not neutral.

## Force-field changes

Read `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md`
first. Verification has three depths, each of which has caught a bug the
previous one missed: the name exists; its element and mass are right; every
emitted bond and angle resolves in `ffbonded.itp`.

Only the `opls_XXX` name reaches GROMACS. A name identifying the wrong element
simulates the wrong chemistry, produces plausible output, and nothing complains.

Several parameters (`CA-S-CA`, `NPY-CA-OH`, `SS`, `NGR`, the ionized types) are
nearest-analog choices with provenance comments, deliberately not QM-validated.
They are recorded in `rqm/opls-typing.md`. **Do not "restore" them** on the
assumption they drifted, and do not pick up QM validation unprompted.

## Style

- Match the surrounding code's idiom, naming, and comment density.
- Comments should say *why*, especially for a threshold that was chosen rather
  than derived. This codebase's constants carry their rationale inline; keep
  that up.
- Prefer composable functions over inheritance. Avoid adding a dependency for
  something small.
- Do not assume happy paths.

## Constraints that are easy to violate locally

- `molecule_name` / residue name: **5 characters or fewer** (GROMACS `.gro`).
- `_fix_heteroatom_bond_types` must run after **any** `SanitizeMol` pass on a
  molecule with ether oxygens.
- Clash warnings on hex-lattice structures are artefacts of an exact lattice,
  not defects.

## Pull requests

Branch from `main`. **Do not stack PRs** — a PR based on another branch runs
zero checks while reporting `MERGEABLE / CLEAN`.
