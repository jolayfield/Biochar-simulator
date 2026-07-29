"""``biochar-md-setup`` — generate GROMACS run directories from a sweep manifest."""

from __future__ import annotations

import argparse
import sys

from biochar.md_setup import (
    ION_PROFILES,
    MDSetupConfig,
    MDSetupError,
    setup_md_from_manifest,
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="biochar-md-setup",
        description="Turn a biochar.sweep manifest into ready-to-submit GROMACS run directories.",
    )
    parser.add_argument("manifest_csv", help="Path to a sweep manifest.csv")
    parser.add_argument("--output-root", required=True, help="Directory to write one subfolder per structure")
    parser.add_argument(
        "--ion-profile",
        default="mn_calcareous_default",
        choices=sorted(ION_PROFILES),
        help="Background electrolyte for the wet stage (default: mn_calcareous_default)",
    )
    parser.add_argument("--solvent-pad-nm", type=float, default=1.2, help="editconf -d padding, nm")
    parser.add_argument("--water-model", default="spce", help="spce | tip3p | tip4p")
    parser.add_argument("--ntomp", type=int, default=8, help="OpenMP threads per mdrun call")
    parser.add_argument("--gmx-bin", default="gmx", help="gmx executable name/path on the run host")
    args = parser.parse_args(argv)

    cfg = MDSetupConfig(
        solvent_pad_nm=args.solvent_pad_nm,
        water_model=args.water_model,
        ion_profile=args.ion_profile,
        ntomp=args.ntomp,
        gmx_bin=args.gmx_bin,
    )

    try:
        results = setup_md_from_manifest(args.manifest_csv, args.output_root, config=cfg)
    except MDSetupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    n_written = sum(1 for r in results if r["run_dir"])
    n_skipped = len(results) - n_written
    for r in results:
        if r["run_dir"]:
            print(f"  [ok]      {r['label']:30s} -> {r['run_dir']}", file=sys.stderr)
        else:
            print(f"  [skipped] {r['label']:30s} -> {r['skipped_reason']}", file=sys.stderr)
    print(f"Wrote {n_written} run director{'y' if n_written == 1 else 'ies'}, skipped {n_skipped}.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
