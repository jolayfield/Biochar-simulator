"""
Biochar Structure Generator

Main API for generating biochar building blocks for GROMACS simulations.
"""

from typing import Optional, Tuple, List, Dict
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem
import numpy as np

from .carbon_skeleton import PAHAssembler, SkeletonValidator
from .heteroatom_assignment import (
    OxygenAssigner,
    HydrogenAssigner,
    HeteroatomValidator,
    CompositionInfo,
    _fix_heteroatom_bond_types,
)
from .geometry_3d import CoordinateGenerator, GeometryValidator, ClashResolver
from .opls_typing import AtomTyper, ChargeAssigner
from .gromacs_export import GromacsExporter
from .validation import ValidationEngine


@dataclass
class GeneratorConfig:
    """Configuration for biochar generator."""

    # Size parameters
    target_num_carbons: int = 50
    size_tolerance: float = 0.10

    # Composition parameters
    H_C_ratio: float = 0.5
    H_C_tolerance: float = 0.10
    O_C_ratio: float = 0.1
    O_C_tolerance: float = 0.10

    # Structural parameters
    aromaticity_percent: float = 90.0
    aromaticity_tolerance: float = 5.0

    # Functional groups — dict mapping group name → exact count to place.
    # e.g. {"phenolic": 3, "carboxyl": 1}
    # Valid keys: phenolic, hydroxyl, carboxyl, ether,
    #             carbonyl*, quinone*, lactone*   (* fall back to phenolic)
    # If None, O_C_ratio controls total oxygen using phenolic groups.
    functional_groups: Optional[Dict[str, int]] = None

    # System setup
    periodic_box: bool = False
    box_size: Optional[np.ndarray] = None

    # Naming and identification
    # Residue name (max 5 chars for GROMACS .gro format)
    # Suggested: BC400, BC600, BC800 (temperature), BCH05, BCO10 (composition)
    molecule_name: str = "BC"

    # Random seed for reproducibility
    seed: Optional[int] = None

    # Ring defects: probability [0, 1) that each ring addition is a pentagon.
    # 0.0 = pure hexagonal PAH (default).  0.1 ≈ 10% pentagons.
    defect_fraction: float = 0.0

    # Maximum graph distance (in bonds along the carbon skeleton) between the
    # two ring carbons that an ether oxygen may bridge.  The ring formed by the
    # C-O-C bridge has (max_ether_span + 2) members.  Minimum enforced = 3
    # (5-membered, furan-like) to avoid strained smaller rings.
    #   3 → 5-membered ring (furan/benzofuran-like) ← default, always flat
    #   4 → 6-membered ring (pyran/chromene-like)   — may cause minor strain
    #   5 → 7-membered ring                         — risk of sheet folding
    # Values > 5 risk cross-sheet bridges that fold the sheet into a nanotube.
    max_ether_span: int = 3

    def __post_init__(self):
        # functional_groups defaults to None → O_C_ratio-driven phenolic placement
        # (no default list here; OxygenAssigner handles the None case)

        # Validate molecule name length
        if len(self.molecule_name) > 5:
            raise ValueError(f"molecule_name must be ≤5 characters (GROMACS .gro format), got '{self.molecule_name}' ({len(self.molecule_name)} chars)")


class BiocharGenerator:
    """Main biochar structure generator."""

    def __init__(self, config: Optional[GeneratorConfig] = None):
        """
        Initialize generator.

        Args:
            config: GeneratorConfig object
        """
        self.config = config or GeneratorConfig()
        self.mol = None
        self.coords = None
        self.composition = None
        self.atom_types = None
        self.charges = None
        self.validation_report = None

    def generate(self) -> Tuple[Chem.Mol, np.ndarray, CompositionInfo]:
        """
        Generate biochar structure.

        Returns:
            (molecule, coordinates, composition_info)
        """
        # Step 1: Generate carbon skeleton
        print(f"Generating carbon skeleton with {self.config.target_num_carbons} carbons...")
        skeleton = self._generate_carbon_skeleton()

        # Step 2: Assign heteroatoms (O, then H)
        print("Assigning oxygen atoms...")
        mol, o_composition = self._assign_oxygens(skeleton.mol)

        print("Assigning hydrogen atoms...")
        mol, composition = self._assign_hydrogens(mol)

        # Carry functional-group bookkeeping from the oxygen step into the
        # final composition (HydrogenAssigner resets it to {}).
        composition.functional_groups = o_composition.functional_groups

        # Step 3: Generate 3D coordinates
        print("Generating 3D coordinates...")
        mol, coords = self._generate_geometry(mol)

        # Force field computation (inside geometry) internally re-sanitizes the
        # mol with RDKit's default aromaticity model, which can mark ether C-O
        # bonds as AROMATIC.  Restore correct single-bond types before typing.
        mol = _fix_heteroatom_bond_types(mol)

        # Step 4: Assign OPLS types and charges
        print("Assigning OPLS-AA atom types and charges...")
        self._assign_opls_properties(mol)

        # Step 5: Validate
        print("Validating structure...")
        self._validate(mol, composition, coords)

        # validation.py calls Chem.SanitizeMol() which can re-mark ether C-O
        # bonds as AROMATIC using RDKit's default model.  Fix again.
        mol = _fix_heteroatom_bond_types(mol)

        # Store results
        self.mol = mol
        self.coords = coords
        self.composition = composition

        return mol, coords, composition

    def export_gromacs(
        self,
        output_directory: str = ".",
        basename: str = "biochar",
    ) -> Tuple[Path, Path, Path]:
        """
        Export to GROMACS files.

        Args:
            output_directory: Output directory path
            basename: Base filename

        Returns:
            (gro_path, top_path, itp_path)
        """
        if self.mol is None or self.coords is None:
            raise RuntimeError("Must call generate() before export_gromacs()")

        exporter = GromacsExporter(output_directory)
        gro_path, top_path, itp_path = exporter.export(
            self.mol,
            self.coords,
            self.atom_types,
            self.charges,
            molecule_name=self.config.molecule_name,
            basename=basename,
            include_periodic_box=self.config.periodic_box,
            box_size=self.config.box_size,
        )

        print(f"\nGROMACS files written:")
        print(f"  Structure:  {gro_path}")
        print(f"  Topology:   {top_path}")
        print(f"  Include:    {itp_path}")

        return gro_path, top_path, itp_path

    def print_summary(self):
        """Print summary of generated structure."""
        if self.composition is None:
            print("No structure generated yet. Call generate() first.")
            return

        print("\n" + "=" * 60)
        print("BIOCHAR STRUCTURE SUMMARY")
        print("=" * 60)
        print(f"\nComposition:")
        print(f"  Carbons:     {self.composition.num_carbons}")
        print(f"  Hydrogens:   {self.composition.num_hydrogens}")
        print(f"  Oxygens:     {self.composition.num_oxygens}")
        print(f"\nRatios:")
        print(f"  H/C ratio:   {self.composition.H_C_ratio:.3f} (target: {self.config.H_C_ratio:.3f})")
        print(f"  O/C ratio:   {self.composition.O_C_ratio:.3f} (target: {self.config.O_C_ratio:.3f})")
        print(f"\nFunctional Groups:")
        if self.composition.functional_groups:
            for group_name, count in self.composition.functional_groups.items():
                print(f"  {group_name}: {count}")
        else:
            print("  None")

        if self.validation_report:
            print(f"\nValidation:")
            print(f"  Status: {'VALID' if self.validation_report[0] else 'INVALID'}")
            if self.validation_report[1]:
                print(f"  Errors: {len(self.validation_report[1])}")
                for error in self.validation_report[1][:3]:
                    print(f"    - {error}")
            if self.validation_report[2]:
                print(f"  Warnings: {len(self.validation_report[2])}")
                for warning in self.validation_report[2][:3]:
                    print(f"    - {warning}")

        print("=" * 60 + "\n")

    # Private methods

    def _generate_carbon_skeleton(self):
        """Generate aromatic carbon skeleton."""
        assembler = PAHAssembler(seed=self.config.seed)
        skeleton = assembler.generate(
            self.config.target_num_carbons,
            self.config.aromaticity_percent,
            defect_fraction=self.config.defect_fraction,
        )

        # Validate skeleton
        valid, errors = SkeletonValidator.validate(skeleton)
        if not valid:
            print(f"Warning: Skeleton validation issues: {errors}")

        return skeleton

    def _assign_oxygens(self, mol: Chem.Mol) -> Tuple[Chem.Mol, CompositionInfo]:
        """Assign oxygen atoms. Returns (mol, composition) including functional-group counts."""
        assigner = OxygenAssigner(seed=self.config.seed,
                                  max_ether_span=self.config.max_ether_span)
        mol, comp = assigner.assign_oxygens(
            mol,
            self.config.O_C_ratio,
            O_C_tolerance=self.config.O_C_tolerance,
            functional_group_preference=self.config.functional_groups,
        )
        return mol, comp

    def _assign_hydrogens(self, mol: Chem.Mol) -> Tuple[Chem.Mol, CompositionInfo]:
        """Assign hydrogen atoms."""
        assigner = HydrogenAssigner(seed=self.config.seed)
        mol, composition = assigner.assign_hydrogens(
            mol,
            self.config.H_C_ratio,
            H_C_tolerance=self.config.H_C_tolerance,
        )

        # Validate
        valid, errors = HeteroatomValidator.validate_ratios(
            composition,
            self.config.H_C_ratio,
            self.config.O_C_ratio,
            self.config.H_C_tolerance,
            self.config.O_C_tolerance,
        )
        if not valid:
            print(f"Warning: Composition validation issues: {errors}")

        return mol, composition

    def _generate_geometry(self, mol: Chem.Mol) -> Tuple[Chem.Mol, np.ndarray]:
        """Generate 3D coordinates with clash resolution."""
        generator = CoordinateGenerator(seed=self.config.seed)
        mol, coords = generator.generate_3d_coordinates(
            mol,
            force_aromatic_planarity=True,
        )

        # Validate and detect clashes
        valid, errors = GeometryValidator.validate_geometry(mol, coords)
        steric_clashes = [e for e in errors if "Steric clash" in e]

        # Resolve clashes if found.
        #
        # When hex-lattice positions are used (generator.used_hex_lattice is
        # True), skip ALL clash resolution and FF passes.  The hex lattice
        # gives perfect 1.42 Å CC bonds; the "clashes" reported are physical
        # features of fused PAH geometry:
        #   * Peri C-C at 2.44 Å (bay regions) — real PAH geometry, not a
        #     clash; below the 0.75×vdW threshold (2.55 Å) but chemically fine.
        #   * Peri H-C: also within expected PAH ranges.
        # Displacing ring carbons via the resolver would shatter the ring
        # geometry.  GROMACS energy minimisation (gmx grompp + em) will relax
        # H/O positions further before any production MD run.
        if steric_clashes and not generator.used_hex_lattice:
            # First pass: iterative clash resolution
            coords = ClashResolver.resolve_clashes(
                mol, coords, max_iterations=15, displacement_step=0.15, use_vdw_radii=True
            )

            # Second pass: force-field refinement.
            coords, _ = generator.validate_and_relax(mol, coords, max_iterations=200)

            # Third pass: final clash resolution if needed
            valid_check, errors_check = GeometryValidator.validate_geometry(mol, coords)
            steric_clashes_after_ff = [e for e in errors_check if "Steric clash" in e]

            if steric_clashes_after_ff:
                coords = ClashResolver.resolve_clashes(
                    mol, coords, max_iterations=10, displacement_step=0.08, use_vdw_radii=True
                )
                coords, _ = generator.validate_and_relax(mol, coords, max_iterations=200)

            # Final validation
            valid, errors = GeometryValidator.validate_geometry(mol, coords)
            steric_clashes_after = [e for e in errors if "Steric clash" in e]
            clashes_resolved = len(steric_clashes) - len(steric_clashes_after)
            print(f"  Clash resolution: {clashes_resolved}/{len(steric_clashes)} clashes resolved")

        if not valid:
            print(f"Warning: Geometry validation issues:")
            for error in errors[:3]:
                print(f"  - {error}")

        # Measure planarity
        planarity, assessment = GeometryValidator.measure_ring_planarity(mol, coords)
        print(f"  Ring planarity: {assessment} (deviation: {planarity:.3f} Å)")

        return mol, coords

    def _assign_opls_properties(self, mol: Chem.Mol):
        """Assign OPLS-AA atom types and charges."""
        typer = AtomTyper()
        self.atom_types = typer.assign_atom_types(mol)

        charger = ChargeAssigner()
        self.charges = charger.assign_charges(mol, self.atom_types)

        print(f"  Atom types assigned: {len(set(self.atom_types.values()))} unique types")
        print(f"  Total charge: {sum(self.charges.values()):.3f} e")

    def _validate(
        self, mol: Chem.Mol, composition: CompositionInfo, coords: np.ndarray
    ):
        """Validate complete structure."""
        is_valid, errors, warnings, metrics = ValidationEngine.validate_complete(
            mol,
            composition,
            coords,
            self.config.H_C_ratio,
            self.config.O_C_ratio,
            self.config.aromaticity_percent,
            self.config.H_C_tolerance,
            self.config.O_C_tolerance,
        )

        self.validation_report = (is_valid, errors, warnings, metrics)

        if not is_valid:
            print(f"  Validation FAILED with {len(errors)} error(s):")
            for error in errors:
                print(f"    - {error}")
        else:
            print(f"  Validation PASSED")
            if warnings:
                print(f"  {len(warnings)} warning(s)")


def generate_biochar(
    target_num_carbons: int = 50,
    H_C_ratio: float = 0.5,
    O_C_ratio: float = 0.1,
    aromaticity_percent: float = 90.0,
    functional_groups: Optional[Dict[str, int]] = None,
    defect_fraction: float = 0.0,
    max_ether_span: int = 5,
    output_directory: str = ".",
    basename: str = "biochar",
    molecule_name: str = "BC",
    seed: Optional[int] = None,
) -> Tuple[Chem.Mol, np.ndarray, Path, Path, Path]:
    """
    Convenience function to generate and export biochar in one call.

    Args:
        target_num_carbons: Target number of carbon atoms.
        H_C_ratio: Target hydrogen-to-carbon ratio.
        O_C_ratio: Target oxygen-to-carbon ratio.  Used to determine total
            oxygen when *functional_groups* is None.
        aromaticity_percent: Target aromaticity percentage.
        functional_groups: Dict mapping functional group name → exact count,
            e.g. ``{"phenolic": 3, "carboxyl": 1}``.

            Supported groups:

            * ``phenolic``  — aromatic C–OH             (1 O per group)
            * ``hydroxyl``  — same as phenolic for pure aromatic PAH (1 O)
            * ``carboxyl``  — aromatic C–C(=O)(OH)      (2 O per group)
            * ``ether``     — aromatic C–O–C bridge     (1 O per group)
            * ``carbonyl``  — not supported; substituted with phenolic
            * ``quinone``   — not supported; substituted with phenolic
            * ``lactone``   — not supported; substituted with phenolic

            If ``None`` (default), the total oxygen count is derived from
            *O_C_ratio* and placed as phenolic (–OH) groups.
        defect_fraction: Probability [0, 1) that each ring added during
            skeleton growth is a 5-membered (pentagon) ring rather than a
            hexagon.  0.0 (default) = pure hexagonal PAH.  Values ~0.1–0.2
            introduce realistic topological disorder.
        max_ether_span: Maximum number of C–C bonds between the two ring
            carbons bridged by each ether oxygen.  Controls the ring size of
            the C–O–C bridge (ring size = max_ether_span + 2).  Default 3
            (5-membered furan/benzofuran-like ring — always stays flat).
            Use 4 for pyran/chromene-like (6-membered) or 5 for 7-membered;
            larger values risk cross-sheet bridges that fold the molecule.
        output_directory: Output directory for GROMACS files.
        basename: Base filename for output files.
        molecule_name: Residue name (max 5 chars). Suggested: BC400, BC600,
            BCH05, BCO10.
        seed: Random seed for reproducibility.

    Returns:
        (molecule, coordinates, gro_path, top_path, itp_path)
    """
    config = GeneratorConfig(
        target_num_carbons=target_num_carbons,
        H_C_ratio=H_C_ratio,
        O_C_ratio=O_C_ratio,
        aromaticity_percent=aromaticity_percent,
        functional_groups=functional_groups,
        defect_fraction=defect_fraction,
        max_ether_span=max_ether_span,
        molecule_name=molecule_name,
        seed=seed,
    )

    generator = BiocharGenerator(config)
    mol, coords, composition = generator.generate()
    generator.print_summary()

    gro_path, top_path, itp_path = generator.export_gromacs(
        output_directory=output_directory,
        basename=basename,
    )

    return mol, coords, gro_path, top_path, itp_path


def generate_biochar_series(
    configurations: List[Dict],
    output_directory: str = ".",
    create_combined_top: bool = True,
    verbose: bool = True,
) -> Dict[str, Tuple[Path, Path, Path]]:
    """
    Generate multiple biochar structures for mixed simulations.

    This function is ideal for creating temperature series, composition series,
    or mixed biochar systems for GROMACS simulations.

    Args:
        configurations: List of configuration dictionaries. Each dict should contain:
            - 'molecule_name' (str, required): Residue name (max 5 chars, e.g., 'BC400')
            - 'target_num_carbons' (int, optional): Default 50
            - 'H_C_ratio' (float, optional): Default 0.5
            - 'O_C_ratio' (float, optional): Default 0.1
            - 'aromaticity_percent' (float, optional): Default 90.0
            - 'seed' (int, optional): For reproducibility
            - 'functional_groups' (dict, optional): e.g. {"phenolic": 3, "carboxyl": 1}

        output_directory: Output directory for all files
        create_combined_top: If True, generate a combined topology for all structures
        verbose: If True, print progress information

    Returns:
        Dictionary mapping molecule_name -> (gro_path, top_path, itp_path)

    Example:
        >>> configs = [
        ...     {'molecule_name': 'BC400', 'H_C_ratio': 0.65, 'O_C_ratio': 0.20},
        ...     {'molecule_name': 'BC600', 'H_C_ratio': 0.55, 'O_C_ratio': 0.12},
        ...     {'molecule_name': 'BC800', 'H_C_ratio': 0.40, 'O_C_ratio': 0.05},
        ... ]
        >>> results = generate_biochar_series(configs, output_directory='output')
    """
    # Create output directory
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    results = {}
    itp_files = {}
    molecule_names = []

    if verbose:
        print("\n" + "=" * 70)
        print(f"BATCH BIOCHAR GENERATION - {len(configurations)} structures")
        print("=" * 70 + "\n")

    for i, config in enumerate(configurations, 1):
        # Extract configuration
        molecule_name = config.get("molecule_name")
        if not molecule_name:
            raise ValueError(f"Configuration {i} missing required 'molecule_name'")

        if len(molecule_name) > 5:
            raise ValueError(
                f"molecule_name '{molecule_name}' exceeds 5 character limit "
                f"(GROMACS .gro format requirement)"
            )

        target_carbons = config.get("target_num_carbons", 50)
        h_c = config.get("H_C_ratio", 0.5)
        o_c = config.get("O_C_ratio", 0.1)
        arom = config.get("aromaticity_percent", 90.0)
        seed = config.get("seed", None)
        fg = config.get("functional_groups", None)

        if verbose:
            print(f"[{i}/{len(configurations)}] Generating {molecule_name}...")
            print(f"  Parameters: {target_carbons}C, H/C={h_c:.2f}, O/C={o_c:.2f}")

        try:
            mol, coords, gro_path, top_path, itp_path = generate_biochar(
                target_num_carbons=target_carbons,
                H_C_ratio=h_c,
                O_C_ratio=o_c,
                aromaticity_percent=arom,
                functional_groups=fg,
                output_directory=str(output_path),
                basename=molecule_name.lower(),
                molecule_name=molecule_name,
                seed=seed,
            )

            results[molecule_name] = (gro_path, top_path, itp_path)
            itp_files[molecule_name] = itp_path
            molecule_names.append(molecule_name)

            if verbose:
                print(f"  ✓ Successfully generated {gro_path.name}\n")

        except Exception as e:
            if verbose:
                print(f"  ✗ Failed: {e}\n")
            raise

    # Generate combined topology file if requested
    if create_combined_top and len(results) > 1:
        combined_top_path = output_path / "combined.top"
        _create_combined_topology(
            itp_files, molecule_names, combined_top_path, verbose
        )

    if verbose:
        print("=" * 70)
        print(f"BATCH GENERATION COMPLETE - {len(results)} structures generated")
        print("=" * 70 + "\n")

        if create_combined_top and len(results) > 1:
            print(f"Combined topology: {combined_top_path}")
            print("\nTo use in GROMACS simulations:")
            print(f"  gmx grompp -f md.mdp -c combined.gro -p combined.top -o topol.tpr")
            print()

    return results


def _create_combined_topology(
    itp_files: Dict[str, Path],
    molecule_names: List[str],
    output_path: Path,
    verbose: bool = True,
) -> Path:
    """
    Create a combined topology file for mixed biochar simulations.

    Args:
        itp_files: Dictionary mapping molecule_name -> itp_path
        molecule_names: List of molecule names in order
        output_path: Path to write combined .top file
        verbose: If True, print status information

    Returns:
        Path to combined topology file
    """
    if verbose:
        print(f"\nCreating combined topology file...")

    with open(output_path, "w") as f:
        # Header
        f.write("; Combined topology for mixed biochar simulation\n")
        f.write("; Auto-generated by Biochar Generator\n")
        f.write(f"; Contains {len(molecule_names)} biochar structures\n\n")

        # Include forcefield
        f.write('#include "oplsaa.ff/forcefield.itp"\n\n')

        # Include all molecule topologies
        f.write("; Molecule definitions\n")
        for mol_name in molecule_names:
            itp_file = itp_files[mol_name]
            f.write(f'#include "{itp_file.name}"\n')

        f.write("\n")

        # System section
        f.write("[ system ]\n")
        f.write("Mixed Biochar System\n\n")

        # Molecules section
        f.write("[ molecules ]\n")
        f.write("; Name            #molecules\n")
        for mol_name in molecule_names:
            f.write(f"{mol_name:20s} 1\n")

    if verbose:
        print(f"✓ Combined topology written to {output_path.name}")

    return output_path


def generate_surface(
    target_num_carbons: int = 50,
    H_C_ratio: float = 0.3,
    O_C_ratio: float = 0.05,
    functional_groups: Optional[Dict[str, int]] = None,
    defect_fraction: float = 0.0,
    pore_diameter: float = 10.0,
    num_sheets: int = 2,
    sheet_overrides: Optional[List[Dict]] = None,
    output_directory: str = ".",
    basename: str = "surface",
    system_name: str = "SLIT",
    seed: Optional[int] = None,
) -> Tuple[list, Path, Path, list]:
    """
    Generate a slit-pore surface system and export to GROMACS files.

    Creates *num_sheets* parallel graphene-like sheets separated by
    *pore_diameter* Ångströms, applies functional groups to each sheet,
    and writes GROMACS-ready ``.gro`` / ``.top`` / ``.itp`` files.

    Args:
        target_num_carbons: Number of carbon atoms per sheet.
        H_C_ratio: Target H/C ratio for each sheet.
        O_C_ratio: Target O/C ratio for each sheet (used when
            *functional_groups* is ``None``).
        functional_groups: Functional groups applied to every sheet,
            e.g. ``{'phenolic': 2, 'ether': 1}``.  Overridden per-sheet
            if *sheet_overrides* is provided.
        pore_diameter: Gap between sheet inner van-der-Waals surfaces,
            in Ångströms.
        num_sheets: Number of parallel sheets (default 2 → one slit pore).
        sheet_overrides: List of per-sheet config dicts (length must equal
            *num_sheets*).  Accepted keys: ``target_num_carbons``,
            ``H_C_ratio``, ``O_C_ratio``, ``functional_groups``,
            ``aromaticity_percent``, ``seed``.  If ``None``, all sheets
            are chemically identical.
        output_directory: Directory for output files.
        basename: Base filename for ``.gro``/``.top``/``.itp`` files.
        system_name: Name written to the ``[ system ]`` section in .top.
        seed: Random seed for reproducibility.

    Returns:
        ``(sheets, gro_path, top_path, itp_paths)``

        * *sheets* — ``List[SheetResult]`` with mol, coords, composition.
        * *gro_path*, *top_path* — :class:`pathlib.Path` objects.
        * *itp_paths* — list of :class:`pathlib.Path` (one per unique
          sheet type).

    Examples::

        # Simple slit pore — two identical sheets, 10 Å pore
        sheets, gro, top, itps = generate_surface(
            target_num_carbons=40,
            functional_groups={'phenolic': 2, 'ether': 1},
            pore_diameter=10.0,
        )

        # Asymmetric pore — different chemistry on each wall
        sheets, gro, top, itps = generate_surface(
            pore_diameter=8.0,
            sheet_overrides=[
                {'functional_groups': {'phenolic': 3}, 'target_num_carbons': 40},
                {'functional_groups': {'carboxyl': 2}, 'target_num_carbons': 50},
            ],
        )
    """
    # Import here to avoid circular imports at module level
    from .surface_builder import SurfaceBuilder, SurfaceConfig

    surface_config = SurfaceConfig(
        target_num_carbons=target_num_carbons,
        H_C_ratio=H_C_ratio,
        O_C_ratio=O_C_ratio,
        functional_groups=functional_groups,
        aromaticity_percent=95.0,
        defect_fraction=defect_fraction,
        pore_type="slit",
        num_sheets=num_sheets,
        pore_diameter=pore_diameter,
        sheet_overrides=sheet_overrides,
        system_name=system_name,
        seed=seed,
    )

    builder = SurfaceBuilder(surface_config)
    sheets, box_vectors = builder.build()

    gro_path, top_path, itp_paths = builder.export_gromacs(
        output_directory=output_directory,
        basename=basename,
    )

    print(f"\nSurface GROMACS files written:")
    print(f"  Structure:  {gro_path}")
    print(f"  Topology:   {top_path}")
    for itp in itp_paths:
        print(f"  Include:    {itp}")
    print(f"\nTo run with GROMACS:")
    print(f"  gmx grompp -f md.mdp -c {gro_path.name} -p {top_path.name} -o topol.tpr")

    return sheets, gro_path, top_path, itp_paths
