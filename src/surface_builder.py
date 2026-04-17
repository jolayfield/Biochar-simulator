"""
Surface Builder for Porous Biochar Systems

Generates slit-pore (Phase 1) and amorphous (Phase 2, deferred) surface
systems consisting of multiple parallel PAH sheets. Each sheet is an
independent biochar molecule produced by the existing BiocharGenerator
pipeline. The builder positions the sheets in 3D, applies periodic
boundary conditions, and exports GROMACS-ready files.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem

from .biochar_generator import BiocharGenerator, GeneratorConfig
from .constants import CARBON_VDW_DIAMETER
from .gromacs_export import ITPFileWriter
from .heteroatom_assignment import _fix_heteroatom_bond_types
from .opls_typing import AtomTyper, ChargeAssigner


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SheetResult:
    """Result of generating a single sheet."""

    mol: Chem.Mol
    coords: np.ndarray  # Angstroms
    composition: object  # CompositionInfo — avoid circular import
    atom_types: Dict[int, str]
    charges: Dict[int, float]
    molecule_name: str  # e.g. "SHT" (identical) or "SHT1" (distinct)


@dataclass
class SurfaceConfig:
    """Configuration for porous surface generation."""

    # --- Sheet chemistry (passed through to BiocharGenerator) ---------------
    target_num_carbons: int = 50
    H_C_ratio: float = 0.3
    O_C_ratio: float = 0.05
    functional_groups: Optional[Dict[str, int]] = None
    aromaticity_percent: float = 95.0

    # --- Pore geometry ------------------------------------------------------
    pore_type: str = "slit"  # "slit" (Phase 1); "amorphous" reserved
    num_sheets: int = 2
    pore_diameter: float = 10.0  # Angstroms — gap between inner vdW surfaces

    # --- Per-sheet overrides ------------------------------------------------
    # If provided, must be a list of dicts (one per sheet) with keys matching
    # GeneratorConfig fields.  If None, all sheets are chemically identical.
    sheet_overrides: Optional[List[Dict]] = None

    # --- Box padding (nm) ---------------------------------------------------
    box_padding_xy: float = 1.0
    box_padding_z: float = 1.0

    # --- Naming -------------------------------------------------------------
    system_name: str = "SLIT"
    sheet_base_name: str = "SHT"  # max 3 chars (+ digit ≤ 5 GROMACS limit)

    # --- Reproducibility ----------------------------------------------------
    seed: Optional[int] = None

    # --- Ring defects -------------------------------------------------------
    # Probability [0, 1) that each ring addition is a 5-membered pentagon.
    # 0.0 = pure hexagonal (default).
    defect_fraction: float = 0.0

    def __post_init__(self):
        self._validate()

    def _validate(self):
        if self.pore_type not in ("slit",):
            raise ValueError(
                f"pore_type must be 'slit' (got '{self.pore_type}'). "
                "Amorphous packing is not yet implemented."
            )
        if self.num_sheets < 2:
            raise ValueError(f"num_sheets must be >= 2 (got {self.num_sheets})")
        if self.pore_diameter <= 0:
            raise ValueError(f"pore_diameter must be > 0 (got {self.pore_diameter})")
        if len(self.sheet_base_name) > 3:
            raise ValueError(
                f"sheet_base_name must be <= 3 chars for GROMACS compatibility "
                f"(got '{self.sheet_base_name}', {len(self.sheet_base_name)} chars)"
            )
        if self.sheet_overrides is not None:
            if len(self.sheet_overrides) != self.num_sheets:
                raise ValueError(
                    f"sheet_overrides length ({len(self.sheet_overrides)}) "
                    f"must equal num_sheets ({self.num_sheets})"
                )


# ---------------------------------------------------------------------------
# Rotation helper
# ---------------------------------------------------------------------------

def _rotation_matrix_from_vectors(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Compute rotation matrix *R* such that ``R @ a = b`` (both unit vectors).

    Uses Rodrigues' rotation formula.  Handles degenerate cases where
    *a* and *b* are parallel or anti-parallel.
    """
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)

    if np.linalg.norm(v) < 1e-8:
        if c > 0:
            return np.eye(3)
        # 180-degree rotation: pick any perpendicular axis
        perp = np.array([1.0, 0.0, 0.0]) if abs(a[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
        perp = perp - np.dot(perp, a) * a
        perp = perp / np.linalg.norm(perp)
        return 2.0 * np.outer(perp, perp) - np.eye(3)

    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0],
    ])
    return np.eye(3) + vx + vx @ vx / (1.0 + c)


# ---------------------------------------------------------------------------
# Surface builder
# ---------------------------------------------------------------------------

class SurfaceBuilder:
    """Build porous surface systems from parallel biochar sheets."""

    def __init__(self, config: SurfaceConfig):
        self.config = config
        self.sheets: List[SheetResult] = []
        self._box_vectors: Optional[np.ndarray] = None  # nm

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def build(self) -> Tuple[List[SheetResult], np.ndarray]:
        """
        Generate all sheets and compute their positioned coordinates.

        Returns:
            ``(list_of_SheetResult, box_vectors_nm)``

        Each :class:`SheetResult`.coords is updated **in place** to world
        coordinates (translated along *z* so sheets are stacked with the
        requested pore gap).
        """
        self.sheets = []
        for i in range(self.config.num_sheets):
            print(f"\n{'='*60}")
            print(f"Generating sheet {i + 1}/{self.config.num_sheets}...")
            print(f"{'='*60}")
            sheet = self._generate_single_sheet(i)
            self.sheets.append(sheet)

        self._position_sheets()
        self._box_vectors = self._compute_box_vectors()
        self._centre_in_box(self._box_vectors)
        return self.sheets, self._box_vectors

    def export_gromacs(
        self,
        output_directory: str = ".",
        basename: str = "surface",
    ) -> Tuple[Path, Path, List[Path]]:
        """
        Export the positioned surface system to GROMACS files.

        Returns:
            ``(gro_path, top_path, list_of_itp_paths)``
        """
        if not self.sheets:
            raise RuntimeError("Must call build() before export_gromacs()")

        # Import here to avoid circular imports at module level
        from .gromacs_export import MultiSheetGROWriter, SurfaceTopologyWriter

        output_dir = Path(output_directory)
        output_dir.mkdir(parents=True, exist_ok=True)

        sheets_identical = self.config.sheet_overrides is None

        # --- Write .itp files -----------------------------------------------
        itp_paths: List[Path] = []
        if sheets_identical:
            itp_path = output_dir / f"{basename}_sheet.itp"
            ITPFileWriter.write(
                str(itp_path),
                self.sheets[0].mol,
                self.sheets[0].atom_types,
                self.sheets[0].charges,
                molecule_name=self.sheets[0].molecule_name,
            )
            itp_paths.append(itp_path)
        else:
            for sheet in self.sheets:
                itp_path = output_dir / f"{basename}_{sheet.molecule_name.lower()}.itp"
                ITPFileWriter.write(
                    str(itp_path),
                    sheet.mol,
                    sheet.atom_types,
                    sheet.charges,
                    molecule_name=sheet.molecule_name,
                )
                itp_paths.append(itp_path)

        # --- Write combined .gro --------------------------------------------
        gro_path = output_dir / f"{basename}.gro"
        MultiSheetGROWriter.write(
            str(gro_path),
            sheets=self.sheets,
            box_vectors=self._box_vectors,
            title=(
                f"Slit pore surface, pore_diameter="
                f"{self.config.pore_diameter:.1f} A, "
                f"{self.config.num_sheets} sheets"
            ),
        )

        # --- Write .top -----------------------------------------------------
        top_path = output_dir / f"{basename}.top"
        SurfaceTopologyWriter.write(
            str(top_path),
            sheets=self.sheets,
            itp_paths=itp_paths,
            sheets_identical=sheets_identical,
            system_name=self.config.system_name,
        )

        return gro_path, top_path, itp_paths

    # ------------------------------------------------------------------
    # Sheet generation
    # ------------------------------------------------------------------

    def _generate_single_sheet(self, sheet_index: int) -> SheetResult:
        """Generate a single biochar sheet (or deep-copy if identical)."""
        # Identical-sheet optimisation: copy first sheet for indices > 0
        if self.config.sheet_overrides is None and sheet_index > 0:
            base = self.sheets[0]
            return SheetResult(
                mol=Chem.Mol(base.mol),
                coords=base.coords.copy(),
                composition=base.composition,
                atom_types=dict(base.atom_types),
                charges=dict(base.charges),
                molecule_name=base.molecule_name,
            )

        # Build GeneratorConfig for this sheet
        gen_kwargs = dict(
            target_num_carbons=self.config.target_num_carbons,
            H_C_ratio=self.config.H_C_ratio,
            O_C_ratio=self.config.O_C_ratio,
            aromaticity_percent=self.config.aromaticity_percent,
            functional_groups=self.config.functional_groups,
            defect_fraction=self.config.defect_fraction,
            seed=self.config.seed,
        )

        # Apply per-sheet overrides if present
        if self.config.sheet_overrides is not None:
            overrides = self.config.sheet_overrides[sheet_index]
            gen_kwargs.update(overrides)

        # Molecule name: same for identical sheets, indexed for distinct
        if self.config.sheet_overrides is None:
            mol_name = self.config.sheet_base_name
        else:
            mol_name = f"{self.config.sheet_base_name}{sheet_index + 1}"

        gen_kwargs["molecule_name"] = mol_name

        config = GeneratorConfig(**gen_kwargs)
        generator = BiocharGenerator(config)
        mol, coords, composition = generator.generate()

        # Ensure clean bond types after pipeline
        mol = _fix_heteroatom_bond_types(mol)

        # OPLS typing (using public classes directly)
        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)
        charger = ChargeAssigner()
        charges = charger.assign_charges(mol, atom_types)

        # Flatten sheet to xy plane
        coords = self._flatten_to_xy(coords)

        return SheetResult(
            mol=mol,
            coords=coords,
            composition=composition,
            atom_types=atom_types,
            charges=charges,
            molecule_name=mol_name,
        )

    # ------------------------------------------------------------------
    # Geometry helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _flatten_to_xy(coords: np.ndarray) -> np.ndarray:
        """
        Rotate coordinates so the best-fit plane coincides with z=0,
        then translate centroid to origin.
        """
        centroid = coords.mean(axis=0)
        centred = coords - centroid

        # SVD: smallest singular vector = plane normal
        _, _, vt = np.linalg.svd(centred)
        normal = vt[2]

        # Rotation to align normal with z-axis
        z_axis = np.array([0.0, 0.0, 1.0])
        R = _rotation_matrix_from_vectors(normal, z_axis)
        rotated = (R @ centred.T).T

        # Zero residual z drift
        rotated[:, 2] -= rotated[:, 2].mean()
        return rotated

    def _position_sheets(self):
        """Translate each sheet along z to create slit pore gaps."""
        for i, sheet in enumerate(self.sheets):
            z_offset = i * (self.config.pore_diameter + CARBON_VDW_DIAMETER)
            sheet.coords[:, 2] += z_offset

    def _compute_box_vectors(self) -> np.ndarray:
        """
        Compute orthogonal periodic box vectors in nm.

        x, y: span of all sheets + 2 * box_padding_xy
        z: span of all sheets + 2 * box_padding_z
        """
        all_coords = np.vstack([s.coords for s in self.sheets])  # Angstroms
        all_coords_nm = all_coords * 0.1

        min_xyz = all_coords_nm.min(axis=0)
        max_xyz = all_coords_nm.max(axis=0)
        extent = max_xyz - min_xyz

        box = np.array([
            extent[0] + 2 * self.config.box_padding_xy,
            extent[1] + 2 * self.config.box_padding_xy,
            extent[2] + 2 * self.config.box_padding_z,
        ])
        return box

    def _centre_in_box(self, box_nm: np.ndarray):
        """
        Translate all sheet coordinates so the system is centred in the box.

        Coordinates remain in Angstroms internally; *box_nm* is in nm.
        """
        all_coords = np.vstack([s.coords for s in self.sheets])
        min_xyz_A = all_coords.min(axis=0)
        max_xyz_A = all_coords.max(axis=0)
        current_centre = (min_xyz_A + max_xyz_A) / 2.0
        box_centre_A = (box_nm * 10.0) / 2.0  # nm -> Angstroms
        shift = box_centre_A - current_centre
        for sheet in self.sheets:
            sheet.coords += shift
