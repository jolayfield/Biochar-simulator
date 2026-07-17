"""
Command-line interface for the Biochar Simulator.

Usage:
    biochar-gen --temperature 600 --feedstock softwood --carbons 80 --name BC600 --seed 42
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

    # Data-driven composition from pyrolysis conditions (mutually consistent with
    # explicit ratio flags — explicit flags win over temperature-derived values)
    parser.add_argument(
        "--temperature", type=float, default=None, metavar="TEMP_C",
        help=(
            "Pyrolysis temperature in °C (100–1000). Derives H/C ratio, O/C ratio, "
            "and aromaticity from the UC Davis Biochar Database model. Any explicit "
            "--hc-ratio / --oc-ratio / --aromaticity flag overrides the derived value."
        ),
    )
    parser.add_argument(
        "--feedstock", type=str, default=None,
        choices=["corn_stover", "grass", "hardwood", "manure", "softwood", "wood"],
        help=(
            "Feedstock type for temperature-model lookup (requires --temperature). "
            "Low-data feedstocks fall back to the pooled curve with a warning."
        ),
    )

    # Composition ratios — default None so --temperature can fill them
    parser.add_argument(
        "--hc-ratio", type=float, default=None,
        help="Target H/C ratio (default: 0.5, or derived from --temperature)",
    )
    parser.add_argument(
        "--oc-ratio", type=float, default=None,
        help="Target O/C ratio (default: 0.1, or derived from --temperature)",
    )
    parser.add_argument(
        "--aromaticity", type=float, default=None,
        help="Target aromaticity percent (0–100; default: 90.0, or derived from --temperature)",
    )

    parser.add_argument(
        "--pH", type=float, default=None, metavar="PH",
        help=(
            "Environmental pH [0–14]. Omit (default) to leave every group "
            "neutral. When set, each titratable site (carboxyl, phenolic, "
            "thiol, aniline N, pyridinic N) is independently ionized with its "
            "Henderson-Hasselbalch probability, and the topology carries a real "
            "net charge for `gmx genion -neutral` to balance"
        ),
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
    fg_group.add_argument("--thiol", type=int, default=None, metavar="N",
                          help="Number of thiol (Ar–SH) groups")
    fg_group.add_argument("--thioether", type=int, default=None, metavar="N",
                          help="Number of thioether (Ar–S–Ar) bridges")

    # Ring-substituting nitrogen (replaces ring carbons; not pendant groups)
    n_group = parser.add_argument_group(
        "ring nitrogen doping",
        "Replace ring carbons with nitrogen (pyridinic/pyrrolic/graphitic).",
    )
    n_group.add_argument("--pyridinic", type=int, default=0, metavar="N",
                         help="Number of pyridinic N (edge 6-ring, no H)")
    n_group.add_argument("--pyrrolic", type=int, default=0, metavar="N",
                         help="Number of pyrrolic N (5-ring N-H; needs --defects > 0)")
    n_group.add_argument("--graphitic", type=int, default=0, metavar="N",
                         help="Number of graphitic/quaternary N (interior 6-ring, no H)")

    # Partial charge method
    parser.add_argument(
        "--charge-method", type=str, default="opls", choices=["opls", "ml", "qm"],
        dest="charge_method",
        help="Partial charge assignment: 'opls' (default, static table), 'ml' "
             "(environment-aware GPR; requires scikit-learn), or 'qm' "
             "(LigParGen-style 1.14*CM1A via AM1; requires a MOPAC binary)",
    )

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
        "thiol": args.thiol,
        "thioether": args.thioether,
    }
    functional_groups = {k: v for k, v in fg_map.items() if v is not None} or None

    # Build config dict: start from loaded config, CLI args override.
    # Composition fields (H_C_ratio etc.) are only set from CLI when explicitly
    # provided (not None); otherwise GeneratorConfig derives them from
    # --temperature / --feedstock or falls back to historical defaults.
    _ALWAYS_OVERRIDE = {
        "target_num_carbons", "defect_fraction", "molecule_name", "seed",
        "functional_groups", "num_pyridinic", "num_pyrrolic", "num_graphitic",
        "charge_method", "temperature", "feedstock",
    }
    cfg_dict = {k: v for k, v in base.items() if k not in _ALWAYS_OVERRIDE}
    cfg_dict.update({
        "target_num_carbons": args.carbons,
        "defect_fraction": args.defects,
        "molecule_name": args.name,
        "seed": args.seed,
        "functional_groups": functional_groups,
        "num_pyridinic": args.pyridinic,
        "num_pyrrolic": args.pyrrolic,
        "num_graphitic": args.graphitic,
        "charge_method": args.charge_method,
        "temperature": args.temperature,
        "feedstock": args.feedstock,
    })
    # Composition: CLI wins only when explicitly supplied
    if args.hc_ratio is not None:
        cfg_dict["H_C_ratio"] = args.hc_ratio
    if args.oc_ratio is not None:
        cfg_dict["O_C_ratio"] = args.oc_ratio
    if args.aromaticity is not None:
        cfg_dict["aromaticity_percent"] = args.aromaticity
    # pH: None is both "not supplied" and the neutral default, so a loaded
    # config's pH survives unless --pH is given explicitly.
    if args.pH is not None:
        cfg_dict["pH"] = args.pH

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
            pH=config.pH,
            defect_fraction=config.defect_fraction,
            num_pyridinic=config.num_pyridinic,
            num_pyrrolic=config.num_pyrrolic,
            num_graphitic=config.num_graphitic,
            charge_method=config.charge_method,
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
