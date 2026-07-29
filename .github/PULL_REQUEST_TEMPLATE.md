## What and why

<!-- What changes, and what problem it solves. If it fixes a bug, state the
     failure it produced -- a reader should be able to tell whether the fix is
     right without reconstructing the bug themselves. -->

## Requirements

<!-- Did behaviour change in an area specified under rqm/? If so, the
     requirement should be updated in this same PR, with `rq index` re-run.
     If this closes a previously-unreferenced scenario, name its rq- ID. -->

- [ ] `rqm/` updated, or this changes no specified behaviour
- [ ] `bash tools/rqm/rq index` re-run if `rqm/` changed

## Verification

<!-- What you actually ran, and what it reported. "Tests pass" is not evidence;
     a count is. Note the skip count for forcefield-backed tests -- a green run
     that skipped them verified nothing. -->

- [ ] `pytest -n auto` full suite, with `GMXDATA` or `BIOCHAR_OPLSAA_FF` set
- [ ] `ruff check .` clean

## Notes for the reviewer

<!-- Anything non-obvious: a threshold chosen rather than derived, a nearest-
     analog force-field value, a deliberate silent fallback, a limitation
     accepted for now. These are the things that get lost. -->
