"""
Command-line interface for the Biochar Simulator.

Usage:
    biochar-gen --carbons 80 --hc-ratio 0.4 --oc-ratio 0.1 --name BC600 --seed 42
    biochar-gen --carbons 50 --phenolic 3 --carboxyl 1 --amino 2 --output-dir ./output
"""

import argparse
import json
import logging
import sys
from pathlib import Path


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="biochar-gen",
        description="Generate biochar molecular structures for GROMACS simulations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Structure size
    parser.add_argument(
        "--carbons", type=int, default=50, metavar="N",
        help="Target number of carbon atoms",
    )

    # Composition ratios
    parser.add_argument("--hc-ratio", type=float, default=0.5, help="Target H/C ratio")
    parser.add_argument("--oc-ratio", type=float, default=0.1, help="Target O/C ratio")
    parser.add_argument(
        "--aromaticity", type=float, default=90.0,
        help="Target aromaticity percent (0–100)",
    )

    # Structural defects
    parser.add_argument(
        "--defects", type=float, default=0.0, metavar="FRAC",
        help="Pentagon defect fraction [0, 1) — 0 = pure hexagonal PAH",
    )

    # Explicit functional groups (overrides --oc-ratio when any are specified)
    fg_group = parser.add_argument_group(
        "functional groups",
        "Explicit group counts. When any are given, --oc-ratio is ignored.",
    )
    fg_group.add_argument("--phenolic", type=int, default=None, metavar="N",
                          help="Number of phenolic (Ar–OH) groups")
    fg_group.add_argument("--carboxyl", type=int, default=None, metavar="N",
                          help="Number of carboxyl (Ar–COOH) groups")
    fg_group.add_argument("--ether", type=int, default=None, metavar="N",
                          help="Number of ether (Ar–O–Ar) bridges")
    fg_group.add_argument("--amino", type=int, default=None, metavar="N",
                          help="Number of amino (Ar–NH2) groups")

    # Identity and output
    parser.add_argument(
        "--name", type=str, default="BC", metavar="NAME",
        help="Residue name written to GROMACS files (max 5 chars)",
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--output-dir", type=str, default=".", metavar="DIR",
        help="Output directory for GROMACS files",
    )
    parser.add_argument(
        "--basename", type=str, default="biochar",
        help="Base filename stem for output files",
    )

    # Config I/O
    parser.add_argument(
        "--save-config", type=str, default=None, metavar="FILE",
        help="Save the resolved GeneratorConfig to a JSON file",
    )
    parser.add_argument(
        "--load-config", type=str, default=None, metavar="FILE",
        help="Load a GeneratorConfig from a JSON file (command-line args override)",
    )

    # Verbosity
    parser.add_argument(
        "-v", "--verbose", action="store_true",
        help="Enable INFO-level logging from the biochar pipeline",
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable DEBUG-level logging",
    )

    return parser


def main(argv=None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Configure logging so pipeline messages are visible when requested
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format="%(levelname)s %(name)s: %(message)s")
    elif args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    from biochar import generate_biochar, GeneratorConfig

    # Base config — start from loaded file if given, then apply CLI overrides
    if args.load_config:
        try:
            with open(args.load_config) as fh:
                base = json.load(fh)
        except (OSError, json.JSONDecodeError) as exc:
            print(f"Error loading config: {exc}", file=sys.stderr)
            return 1
    else:
        base = {}

    # Collect explicit functional groups
    fg_map = {
        "phenolic": args.phenolic,
        "carboxyl": args.carboxyl,
        "ether": args.ether,
        "amino": args.amino,
    }
    functional_groups = {k: v for k, v in fg_map.items() if v is not None} or None

    # Build config dict, CLI values win over loaded config
    cfg_dict = {
        "target_num_carbons": args.carbons,
        "H_C_ratio": args.hc_ratio,
        "O_C_ratio": args.oc_ratio,
        "aromaticity_percent": args.aromaticity,
        "defect_fraction": args.defects,
        "molecule_name": args.name,
        "seed": args.seed,
        "functional_groups": functional_groups,
        **{k: v for k, v in base.items() if k not in {
            "target_num_carbons", "H_C_ratio", "O_C_ratio",
            "aromaticity_percent", "defect_fraction", "molecule_name",
            "seed", "functional_groups",
        }},
    }

    try:
        config = GeneratorConfig(**cfg_dict)
    except (ValueError, TypeError) as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 1

    if args.save_config:
        try:
            with open(args.save_config, "w") as fh:
                json.dump(config.to_dict(), fh, indent=2)
            print(f"Config saved to {args.save_config}")
        except OSError as exc:
            print(f"Could not save config: {exc}", file=sys.stderr)
            return 1

    try:
        mol, coords, gro_path, top_path, itp_path = generate_biochar(
            target_num_carbons=config.target_num_carbons,
            H_C_ratio=config.H_C_ratio,
            O_C_ratio=config.O_C_ratio,
            aromaticity_percent=config.aromaticity_percent,
            functional_groups=config.functional_groups,
            defect_fraction=config.defect_fraction,
            molecule_name=config.molecule_name,
            seed=config.seed,
            output_directory=args.output_dir,
            basename=args.basename,
        )
    except Exception as exc:
        print(f"Generation failed: {exc}", file=sys.stderr)
        return 1

    print(f"Structure:  {gro_path}")
    print(f"Topology:   {top_path}")
    print(f"Include:    {itp_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
