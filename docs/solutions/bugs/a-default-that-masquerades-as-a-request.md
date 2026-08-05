---
title: A default that masquerades as a request silently discards a loaded config
date: 2026-08-04
category: bugs
module: cli
problem_type: silent_data_loss
component: argument_parsing
severity: high
applies_when:
  - "Merging a config file with command-line overrides"
  - "Adding a flag with a non-None default to biochar-gen"
  - "Reviewing anything that reads argparse values to decide precedence"
tags:
  - cli
  - argparse
  - configuration
  - reproducibility
  - silent-failure
---

# A default that masquerades as a request silently discards a loaded config

## Context

`biochar-gen --save-config` writes the resolved `GeneratorConfig` to JSON, and `--load-config`
reads one back. The pair exists so a run can be reproduced later.

`main()` merged the two sources with a set named `_ALWAYS_OVERRIDE`: twelve fields — including
`target_num_carbons`, `molecule_name`, `seed`, `temperature`, `feedstock` and `charge_method` —
were taken from the parsed arguments unconditionally, and the loaded config supplied only the
rest.

## The problem

A parsed argparse namespace cannot tell a flag the user typed from one left at its default. Both
arrive as plain values. So "the command line always wins" was implemented as "the parser's
defaults always win", and a reload with no other flags overwrote the file with them:

```
saved:    carbons=10  name=TST  seed=1234  temperature=600  feedstock=softwood
reloaded: carbons=50  name=BC   seed=None  temperature=None feedstock=None
```

Nothing warned. The config file was read, parsed, and then largely thrown away.

Two things made it worse than a plain override bug:

- **`seed` is the reproducibility mechanism**, and it is precisely the field a user does not
  re-type on the reload, because carrying it is what the file is for.
- **`H_C_ratio`, `O_C_ratio` and `aromaticity_percent` survived**, because they took a different
  code path that only overrode when the flag was not `None`. So the reload kept the numbers
  derived from 600 °C while dropping the 600 °C that explained them — a config that still builds a
  plausible structure and no longer records what it rests on.

The three-way split was itself the tell: `pH` had already been special-cased with a comment
noting that `None` means both "not supplied" and "the neutral default". That comment describes the
whole bug, applied to one field.

## The fix

Ask argparse which options were actually present, by re-parsing the same argv against a parser
whose defaults are all `SUPPRESS`:

```python
def _supplied_dests(argv) -> set:
    probe = _build_parser()
    for action in probe._actions:
        action.default = argparse.SUPPRESS
    return set(vars(probe.parse_args(argv)))
```

With that, one rule replaces all three cases: the loaded file is the base, a typed flag overrides
it, and an untyped flag only fills a field the file does not carry.

```python
for dest, field in _CONFIG_FIELDS.items():
    value = getattr(args, dest)
    if dest in supplied:
        cfg_dict[field] = value
    elif field not in cfg_dict and value is not None:
        cfg_dict[field] = value
```

Note what this gets right that a "is the value non-default?" heuristic does not: typing
`--carbons 50` when 50 *is* the default still overrides a loaded 14. The user made a request;
that it coincides with the default is irrelevant.

Two smaller points came with it. `GeneratorConfig.from_dict` replaced `GeneratorConfig(**cfg_dict)`,
because `from_dict` converts a JSON list back to the `np.ndarray` that `box_size` expects. And the
six functional-group flags are resolved as one set — any of them replaces a loaded group
dictionary rather than merging into it, since the counts compete for the same edge sites and a set
assembled from two sources is a composition nobody chose.

## How to avoid it

- **A default is a fallback, never an override.** Any time you merge a config file with command-line
  arguments, you need the supplied-vs-default distinction. There is no way to recover it from the
  parsed namespace afterwards.
- **A field-by-field special case is a smell.** `pH`'s `if args.pH is not None` was a correct fix to
  one instance of a general bug. When you find yourself writing the second one, the rule is missing.
- **A round-trip deserves a round-trip test.** `test_save_config_creates_json` and
  `test_load_config_used_as_base` both passed throughout: one checked that saving writes a file, the
  other that loading does not crash. Neither compared a reloaded config against the saved one, which
  is the only assertion that could have caught this.

`tests/test_cli_arguments_contract.py` now guards the rule, and
`rqm/cli-arguments.md` states it. It also guards the mapping itself: a new flag that reaches neither
`_CONFIG_FIELDS`, the group set, nor the declared local options fails the test, because an option
that parses and then does nothing is this layer's most likely defect.
