"""Command-line interface for the biochar parameter-sweep driver.

Run a declarative grid of biochar structures from a YAML/JSON config:

    biochar-sweep run sweep.yaml
    biochar-sweep run sweep.yaml --output-dir runs/grid1 --quiet
    biochar-sweep template > sweep.yaml      # write a starter config

The config describes the sweep axes; the driver expands their Cartesian
product, builds each structure with seed-retry + strict/fallback handling,
and writes a manifest (manifest.csv + manifest.json) recording the achieved
composition and validation status of every structure.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

_TEMPLATE = """\
# Biochar structure-sweep config
# Expands the Cartesian product of every list under `axes`, applies `fixed`
# to all points, and writes GROMACS files + a manifest under output_directory.
name: pfas_grid
output_directory: sweep_out

# Applied identically to every grid point (any GeneratorConfig field).
fixed:
  target_num_carbons: 100
  charge_method: opls          # opls | gasteiger | qm (qm needs MOPAC)

# Each key is a GeneratorConfig field; the sweep takes the product of the lists.
# temperature + feedstock derive H/C, O/C, and aromaticity from the built-in
# temperature model (see biochar.temperature_model.VALID_FEEDSTOCKS).
axes:
  temperature: [300, 400, 500, 600, 700]
  feedstock: [softwood, hardwood]

seed: 0                  # base random seed; retries use seed, seed+1, ...
max_retries: 8           # strict-mode seeds to try before falling back
on_validation_fail: fallback   # fallback | skip | strict
"""


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="biochar-sweep",
        description="Run a declarative parameter sweep of biochar structures.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="Run a sweep from a YAML/JSON config file.")
    p_run.add_argument("config", help="Path to a YAML or JSON sweep config.")
    p_run.add_argument(
        "--output-dir", default=None, metavar="DIR",
        help="Override output_directory from the config.",
    )
    p_run.add_argument(
        "--quiet", action="store_true",
        help="Suppress per-structure generator chatter (still prints progress).",
    )
    p_run.add_argument("--verbose", action="store_true", help="INFO-level logging.")
    p_run.add_argument("--debug", action="store_true", help="DEBUG-level logging.")

    sub.add_parser("template", help="Print a starter sweep config to stdout.")
    return parser


def _progress(i: int, n: int, res) -> None:
    tag = {
        "strict_pass": "ok   ",
        "fallback": "fallbk",
        "failed": "FAIL ",
    }.get(res.status, res.status)
    extra = ""
    if res.status != "failed":
        extra = f"  {res.molecular_formula}  H/C={res.H_C_ratio}  O/C={res.O_C_ratio}"
    elif res.error:
        extra = f"  {res.error[:60]}"
    print(f"[{i:>3}/{n}] {tag}  {res.label:<28}{extra}", file=sys.stderr)


def main(argv=None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "template":
        sys.stdout.write(_TEMPLATE)
        return 0

    if getattr(args, "debug", False):
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s %(name)s: %(message)s")
    elif getattr(args, "verbose", False):
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    from biochar.sweep import run_sweep, SweepError

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        print(f"Config not found: {cfg_path}", file=sys.stderr)
        return 1

    try:
        summary = run_sweep(
            str(cfg_path),
            output_directory=args.output_dir,
            progress_callback=_progress,
            quiet=args.quiet,
        )
    except SweepError as exc:
        print(f"Sweep config error: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:  # pragma: no cover
        print(f"Sweep failed: {exc}", file=sys.stderr)
        return 1

    print()
    print(f"Sweep '{summary['name']}' complete.")
    print(f"  grid points : {summary['n_points']}")
    print(f"  built       : {summary['n_built']}  "
          f"(strict {summary['n_strict_pass']}, fallback {summary['n_fallback']})")
    print(f"  failed      : {summary['n_failed']}")
    print(f"  output      : {summary['output_directory']}")
    print(f"  manifest    : {summary['manifest_csv']}")
    print(f"                {summary['manifest_json']}")
    return 0 if summary["n_failed"] == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
