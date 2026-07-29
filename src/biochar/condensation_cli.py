"""`biochar-condense` — set up a Wood et al. 2024 condensation run.

Two entry points, both writing a ready-to-run (setup-only) directory:

    biochar-condense generate   --output-dir DIR --copies N --htt T [--carbons ...]
    biochar-condense from-files --output-dir DIR --gro mol.gro --itp mol.itp --copies N --htt T

Each writes: the packed-topology + four Wood `.mdp`s + `run_condensation.sh`
(pack N copies of ONE molecule, then EM -> NVT -> NPT anneal -> NPT final),
plus (unless disabled) the surface `run_surface.sh` and the `analyze.sh`
validation scripts. No `gmx` is invoked.
"""

from __future__ import annotations

import argparse

from .condensation import (
    add_surface_and_validation,
    generate_and_condense,
    setup_condensation,
)


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--output-dir", required=True, help="run directory to create")
    p.add_argument("--copies", type=int, required=True, help="number of molecule copies to pack")
    p.add_argument("--htt", type=float, required=True,
                   help="pyrolysis HTT (°C); sets the anneal peak temperature/timestep")
    p.add_argument("--box", type=float, default=None,
                   help="cubic packing-box side (nm); default = a loose auto-estimate")
    p.add_argument("--repeats", type=int, default=3, help="independent annealing repeats (default 3)")
    p.add_argument("--gap", type=float, default=10.0, help="surface vacuum gap in nm (default 10)")
    p.add_argument("--which-repeat", type=int, default=1,
                   help="repeat whose condensed bulk becomes the surface (default 1)")
    p.add_argument("--no-surface", action="store_true", help="skip the surface-creation setup")
    p.add_argument("--no-validation", action="store_true", help="skip the validation-analysis setup")
    p.add_argument("--gmx", default="gmx", help="gmx binary name in the generated scripts")
    p.add_argument("--ntomp", type=int, default=8, help="OpenMP threads in the generated scripts")


def _finish(args) -> int:
    if not args.no_surface or not args.no_validation:
        add_surface_and_validation(
            args.output_dir,
            which_repeat=args.which_repeat,
            gap_nm=args.gap,
            surface=not args.no_surface,
            validation=not args.no_validation,
            gmx_bin=args.gmx,
            ntomp=args.ntomp,
        )
    print(f"Condensation setup written to {args.output_dir}")
    print("Review, then run:  ./run_condensation.sh"
          + ("  →  ./run_surface.sh" if not args.no_surface else "")
          + ("  →  ./analyze.sh" if not args.no_validation else ""))
    return 0


def _cmd_generate(args) -> int:
    from .biochar_generator import GeneratorConfig
    from .constants import WOOD_PENTAGON_FRACTION, WOOD_HEPTAGON_FRACTION

    if args.wood_curvature:
        defect_fraction = WOOD_PENTAGON_FRACTION
        heptagon_fraction = WOOD_HEPTAGON_FRACTION
    else:
        defect_fraction = args.defect_fraction
        heptagon_fraction = args.heptagon_fraction

    cfg = GeneratorConfig(
        target_num_carbons=args.carbons,
        H_C_ratio=args.hc,
        O_C_ratio=args.oc,
        defect_fraction=defect_fraction,
        heptagon_fraction=heptagon_fraction,
        molecule_name=args.name,
        seed=args.seed,
        strict=False,
    )
    generate_and_condense(
        args.output_dir, args.copies, generator_config=cfg, htt_c=args.htt,
        box_nm=args.box, n_repeats=args.repeats, gmx_bin=args.gmx, ntomp=args.ntomp,
    )
    return _finish(args)


def _cmd_from_files(args) -> int:
    setup_condensation(
        args.output_dir, args.gro, args.itp, args.copies, htt_c=args.htt,
        box_nm=args.box, n_repeats=args.repeats, gmx_bin=args.gmx, ntomp=args.ntomp,
    )
    return _finish(args)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="biochar-condense", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    gen = sub.add_parser("generate", help="generate one molecule fresh, then condense N copies")
    gen.add_argument("--carbons", type=int, default=100, help="target carbon count for the molecule")
    gen.add_argument("--hc", type=float, default=0.5, help="target H/C ratio")
    gen.add_argument("--oc", type=float, default=0.1, help="target O/C ratio")
    gen.add_argument("--name", default="BC", help="residue/moleculetype name (<=5 chars)")
    gen.add_argument("--seed", type=int, default=None, help="RNG seed for the molecule")
    gen.add_argument("--defect-fraction", type=float, default=0.0,
                     help="per-ring pentagon probability (positive curvature); default 0")
    gen.add_argument("--heptagon-fraction", type=float, default=0.0,
                     help="per-ring heptagon probability (negative curvature); default 0")
    gen.add_argument("--wood-curvature", action="store_true",
                     help="use Wood et al. 2024 curvature ratios "
                          "(pentagon 2/13, heptagon 1/13); overrides the two fractions")
    _add_common(gen)
    gen.set_defaults(func=_cmd_generate)

    ff = sub.add_parser("from-files", help="condense N copies of an existing .gro/.itp molecule")
    ff.add_argument("--gro", required=True, help="single-molecule .gro")
    ff.add_argument("--itp", required=True, help="single-molecule .itp")
    _add_common(ff)
    ff.set_defaults(func=_cmd_from_files)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
