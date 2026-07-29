# Contributing

## Setup

```bash
conda install -c conda-forge rdkit          # rdkit is why conda is recommended
pip install -e ".[dev,ml,viz]"
bash tools/hooks/install.sh                 # installs the pre-commit gate
```

`tools/hooks/install.sh` sets `core.hooksPath`, so hook changes reach everyone
without a reinstall.

## Tests

The suite has two tiers. **15 of 842 tests account for roughly 98% of the
runtime** — they are end-to-end regressions that generate real structures across
several seeds. They are marked `slow` so that the tier you run constantly stays
fast.

```bash
pytest -m "not slow" -n auto   # ~24s  -- what the pre-commit hook runs
pytest -n auto                 # full suite -- what CI runs
pytest -m slow                 # just the regressions
```

Always pass `-n auto`. Serially the fast tier is 5m39s instead of 24s.

Two things to know:

- **Tests need rdkit**, so a bare `python3` will not work. Use the environment
  you installed into; set `BIOCHAR_PYTHON` if the hook cannot find it.
- **Forcefield-backed tests skip silently** without a discoverable `oplsaa.ff`.
  A green run that skipped them has verified nothing. Set `GMXDATA` or
  `BIOCHAR_OPLSAA_FF` and check the skip count.

If you add a test that takes more than a few seconds, mark it `slow`. Mark the
individual case, not the whole class — a class-level marker exiles its fast
siblings from the tier that runs before every commit.

## Requirements

Behaviour is specified in `rqm/`, with stable `rq-XXXXXXXX` IDs linking each
requirement to the tests defending it. Start at `rqm/ARCHITECTURE.md`.

**When you change behaviour in a specified area, update the requirement in the
same change.** Then:

```bash
bash tools/rqm/rq stamp    # give new headings/scenarios their IDs
bash tools/rqm/rq index    # rebuild the registry; fails on duplicate IDs
bash tools/rqm/rq check    # stale refs are errors, unreferenced are warnings
```

Reference a requirement from the test defending it with a comment:

```python
# rq-3789a326
def test_amino_nitrogen_types_as_aniline_n(self): ...
```

`rq check` reporting a requirement as unreferenced means no test claims to
defend it. That is the signal the system exists to produce — a coverage gap, not
noise. Some are currently open and being closed incrementally.

`rqm/` and `docs/solutions/` are complements: `rqm/` states the invariant,
`docs/solutions/` records how it was learned. A fix usually earns an entry in
both.

## Force-field changes

Read `docs/solutions/conventions/verify-opls-types-against-real-forcefield.md`
before touching atom types or the GROMACS type map. Verification has three
depths and each has caught a bug the previous one missed:

1. The mapped `opls_XXX` name exists in `atomtypes.atp`.
2. Its element and mass match the internal type's intent.
3. Every bond and angle the topology emits resolves in `ffbonded.itp`.

Only the `opls_XXX` name reaches GROMACS. A name identifying the wrong element
simulates the wrong chemistry and nothing downstream complains.

## Pull requests

CI runs lint, a conda-forge install check, requirements traceability, and the
full suite on Python 3.9–3.12.

**Do not stack PRs.** A PR based on another branch runs *zero* checks while
reporting `MERGEABLE / CLEAN` — that status is meaningless there. Branch from
`main`.
