"""
Surface Builder for Porous Biochar Systems

Generates slit-pore and amorphous porous surface systems consisting of multiple
parallel PAH sheets. Each sheet is an
independent biochar molecule produced by the existing BiocharGenerator
pipeline. The builder positions the sheets in 3D, applies periodic
boundary conditions, and exports GROMACS-ready files.
"""

from __future__ import annotations

import logging
import copy
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem

logger = logging.getLogger(__name__)

from ..pipeline.biochar_generator import BiocharGenerator, GeneratorConfig
from ..constants import CARBON_VDW_DIAMETER
from ..export.gromacs_export import ITPFileWriter
from ..pipeline.heteroatom_assignment import _fix_heteroatom_bond_types
from ..pipeline.opls_typing import AtomTyper, ChargeAssigner


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class SheetResult:
    """
    Result of generating and positioning a single sheet.

    Attributes:
        mol: RDKit molecule with 3-D conformer and OPLS-AA types.
        coords: Atomic coordinates in Ångströms, shape ``(N, 3)``.
            Updated in-place by :meth:`SurfaceBuilder.build` to world
            coordinates (translated along *z* after flattening).
        composition: Atom counts and ratios
            (:class:`~heteroatom_assignment.CompositionInfo`).
        atom_types: Mapping of atom index → OPLS-AA type string.
        charges: Mapping of atom index → partial charge (electrons).
        molecule_name: Residue name used in ``.gro`` and ``.itp`` files.
            Identical for all sheets when *sheet_overrides* is ``None``;
            indexed (e.g. ``SHT1``, ``SHT2``) for distinct sheets.
    """

    mol: Chem.Mol
    coords: np.ndarray  # Angstroms
    composition: object  # CompositionInfo — avoid circular import
    atom_types: Dict[int, str]
    charges: Dict[int, float]
    molecule_name: str  # e.g. "SHT" (identical) or "SHT1" (distinct)


@dataclass
class SurfaceConfig:
    """
    Configuration for :class:`SurfaceBuilder`.

    Attributes:
        target_num_carbons: Carbon atoms per sheet (passed to
            :class:`~biochar_generator.BiocharGenerator`).
        H_C_ratio: Target H/C ratio for each sheet.
        O_C_ratio: Target O/C ratio (used when *functional_groups* is
            ``None``).
        functional_groups: Functional groups applied to every sheet, e.g.
            ``{"phenolic": 2, "ether": 1}``.  Overridden per-sheet via
            *sheet_overrides*.
        aromaticity_percent: Target aromatic carbon fraction (percentage).
        pore_type: Pore geometry.  ``"slit"`` stacks parallel sheets along
            *z*; ``"amorphous"`` performs random rigid-body placement with
            steric rejection (each sheet independently rotated/translated).
        num_sheets: Number of parallel sheets.  Must be ≥ 2.
        pore_diameter: Gap between the inner van-der-Waals surfaces of
            adjacent sheets, in Ångströms.  The centre-to-centre sheet
            separation is ``pore_diameter + 3.4 Å`` (slit pore only).
        max_attempts: Maximum random placement attempts per sheet for
            ``pore_type="amorphous"`` before raising ``RuntimeError`` or
            triggering *amorphous_fallback*.
        amorphous_fallback: What to do when amorphous packing fails within
            *max_attempts*.  ``None`` (default) re-raises ``RuntimeError``.
            ``"slit"`` emits a warning and falls back to slit-pore geometry.
        min_separation: Minimum allowed inter-sheet atom-atom distance
            (Ångströms) used as the steric-rejection criterion for
            ``pore_type="amorphous"``.
        sheet_overrides: Per-sheet configuration overrides.  If provided,
            must be a list of dicts (one per sheet) with keys matching
            :class:`~biochar_generator.GeneratorConfig` fields.  If ``None``,
            all sheets are chemically identical (one ``.itp`` with
            ``count = num_sheets``).
        box_padding_xy: Padding added to each side of the bounding box in
            *x* and *y*, in nm.
        box_padding_z: Padding added to each side of the bounding box in
            *z*, in nm.
        system_name: Name written to the ``[ system ]`` section of the
            ``.top`` file.
        sheet_base_name: Residue name prefix for sheet molecules.  Must be
            ≤ 3 characters so that the full name (prefix + digit) fits within
            GROMACS' 5-character residue-name limit.
        seed: Random seed passed to each :class:`~biochar_generator.BiocharGenerator`.
        defect_fraction: Pentagon insertion probability for each sheet (see
            :attr:`~biochar_generator.GeneratorConfig.defect_fraction`).
    """

    # --- Sheet chemistry (passed through to BiocharGenerator) ---------------
    target_num_carbons: int = 50
    H_C_ratio: float = 0.3
    O_C_ratio: float = 0.05
    functional_groups: Optional[Dict[str, int]] = None
    aromaticity_percent: float = 95.0

    # --- Pore geometry ------------------------------------------------------
    pore_type: str = "slit"  # "slit" (parallel stack) or "amorphous"
    num_sheets: int = 2
    pore_diameter: float = 10.0  # Angstroms — gap between inner vdW surfaces

    # --- Amorphous packing --------------------------------------------------
    # Maximum random placement attempts per sheet before giving up.
    max_attempts: int = 500
    # Graceful degradation: None = raise RuntimeError; "slit" = warn + fall back.
    amorphous_fallback: Optional[str] = None
    # Minimum allowed inter-sheet atom-atom distance (Angstroms).
    min_separation: float = 3.0

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

    # --- Environmental pH ---------------------------------------------------
    # None (default) leaves every sheet neutral.  When set, each sheet's
    # titratable sites are independently ionized (see GeneratorConfig.pH) and
    # the stack carries a real net charge for `gmx genion -neutral` to balance.
    #
    # Setting pH disables the identical-sheet optimisation: protonation is a
    # per-site random draw, so two sheets built at the same pH are independent
    # samples rather than copies.  Each sheet is generated from its own seed and
    # gets its own .itp accordingly.
    pH: Optional[float] = None

    # --- Ring defects -------------------------------------------------------
    # Probability [0, 1) that each ring addition is a 5-membered pentagon.
    # 0.0 = pure hexagonal (default).
    defect_fraction: float = 0.0

    # --- Strict validation --------------------------------------------------
    # Propagated to each sheet's GeneratorConfig.  Set False for lenient mode.
    strict: bool = True

    def __post_init__(self):
        self._validate()

    def _validate(self):
        if self.pore_type not in ("slit", "amorphous"):
            raise ValueError(
                f"pore_type must be 'slit' or 'amorphous' "
                f"(got '{self.pore_type}')."
            )
        if self.amorphous_fallback not in (None, "slit"):
            raise ValueError(
                f"amorphous_fallback must be None or 'slit' "
                f"(got '{self.amorphous_fallback}')."
            )
        if self.num_sheets < 2:
            raise ValueError(f"num_sheets must be >= 2 (got {self.num_sheets})")
        if self.pore_diameter <= 0:
            raise ValueError(f"pore_diameter must be > 0 (got {self.pore_diameter})")
        if self.max_attempts < 1:
            raise ValueError(f"max_attempts must be >= 1 (got {self.max_attempts})")
        if self.min_separation <= 0:
            raise ValueError(
                f"min_separation must be > 0 (got {self.min_separation})"
            )
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
        if self.box_padding_xy <= 0:
            raise ValueError(
                f"box_padding_xy must be > 0 nm (got {self.box_padding_xy})"
            )
        if self.box_padding_z <= 0:
            raise ValueError(
                f"box_padding_z must be > 0 nm (got {self.box_padding_z})"
            )
        if self.box_padding_xy > 10:
            logger.warning(
                "box_padding_xy=%.1f nm is unusually large (> 10 nm). "
                "Possible unit error? Expected value in nm, not Å.",
                self.box_padding_xy,
            )
        if self.box_padding_z > 10:
            logger.warning(
                "box_padding_z=%.1f nm is unusually large (> 10 nm). "
                "Possible unit error? Expected value in nm, not Å.",
                self.box_padding_z,
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


def _random_rotation_matrix(rng: np.random.Generator) -> np.ndarray:
    """
    Draw a uniformly-distributed random rotation matrix in SO(3).

    Samples a uniform random unit quaternion (Shoemake's method) and
    converts it to a 3x3 rotation matrix.  Self-contained (numpy only).
    """
    u1, u2, u3 = rng.random(3)
    q1 = np.sqrt(1.0 - u1) * np.sin(2.0 * np.pi * u2)
    q2 = np.sqrt(1.0 - u1) * np.cos(2.0 * np.pi * u2)
    q3 = np.sqrt(u1) * np.sin(2.0 * np.pi * u3)
    q4 = np.sqrt(u1) * np.cos(2.0 * np.pi * u3)
    # Quaternion (x, y, z, w) = (q1, q2, q3, q4) -> rotation matrix
    x, y, z, w = q1, q2, q3, q4
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


# ---------------------------------------------------------------------------
# Surface builder
# ---------------------------------------------------------------------------

class SurfaceBuilder:
    """
    Build slit-pore surface systems from parallel biochar sheets.

    The builder orchestrates a three-phase pipeline:

    1. **Sheet generation** — creates each sheet via :class:`BiocharGenerator`,
       assigns OPLS-AA types and charges, then flattens the sheet to the *xy*
       plane using SVD.
    2. **Positioning** — translates sheet *i* to
       ``z = i × (pore_diameter + 3.4 Å)`` so consecutive sheets are separated
       by the requested pore gap.
    3. **Box computation** — wraps all sheets in an orthogonal periodic box
       padded by :attr:`SurfaceConfig.box_padding_xy` (nm) in *x*/*y* and
       :attr:`SurfaceConfig.box_padding_z` (nm) in *z*, then centres the
       system inside the box.

    When :attr:`SurfaceConfig.sheet_overrides` is ``None`` **and**
    :attr:`SurfaceConfig.pH` is ``None`` all sheets are chemically identical:
    only the first sheet is generated; the rest are deep-copied.  A single
    ``.itp`` is produced and the ``.top`` lists it with ``count = num_sheets``.

    Setting ``pH`` breaks that identity — protonation is a per-site random draw,
    so sheets built at the same pH are independent samples, not copies.  Each
    sheet is then generated from its own seed and written to its own ``.itp``.

    Use :func:`generate_surface` (below) for a one-call convenience wrapper.

    Examples::

        config = SurfaceConfig(
            target_num_carbons=60,
            functional_groups={"phenolic": 2},
            pore_diameter=12.0,
            num_sheets=2,
            seed=42,
        )
        builder = SurfaceBuilder(config)
        sheets, box_nm = builder.build()
        gro, top, itps = builder.export_gromacs("output", "slit_pore")
    """

    def __init__(self, config: SurfaceConfig):
        self.config = config
        self.sheets: List[SheetResult] = []
        self._box_vectors: Optional[np.ndarray] = None  # nm
        # config.pore_type is the request; this is what was actually built. They
        # differ when amorphous packing fails and amorphous_fallback substitutes
        # a slit stack, and everything describing the surface -- the .gro title
        # above all -- must read the second, not the first.
        self._realised_pore_type: str = config.pore_type
        self._packing_fell_back: bool = False

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

        if self.config.pore_type == "amorphous":
            # Random rigid-body packing with steric rejection.  The box is
            # sized to give enough free volume for all sheets, then each
            # sheet is randomly rotated and translated inside it.
            self._box_vectors = self._compute_amorphous_box_vectors()
            try:
                self._pack_amorphous(self._box_vectors)
            except RuntimeError as exc:
                if self.config.amorphous_fallback == "slit":
                    warnings.warn(
                        f"Amorphous packing failed ({exc}). "
                        f"Falling back to slit-pore geometry "
                        f"(amorphous_fallback='slit').",
                        stacklevel=2,
                    )
                    self._position_sheets()
                    self._box_vectors = self._compute_box_vectors()
                    self._centre_in_box(self._box_vectors)
                    self._realised_pore_type = "slit"
                    self._packing_fell_back = True
                else:
                    raise
        else:
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
        from ..export.gromacs_export import MultiSheetGROWriter, SurfaceTopologyWriter

        output_dir = Path(output_directory)
        output_dir.mkdir(parents=True, exist_ok=True)

        sheets_identical = self._sheets_identical

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
            title=self._structure_title(),
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

    def _structure_title(self) -> str:
        """The .gro title line: the geometry built, and why if it was substituted.

        This is the only human-readable record of what the geometry is, and it is
        read weeks later by someone deciding whether a trajectory answers their
        question.
        """
        if self.realised_pore_type == "amorphous":
            return f"Amorphous surface, {self.config.num_sheets} sheets"
        suffix = (
            " (amorphous packing failed; slit fallback)"
            if self._packing_fell_back
            else ""
        )
        return (
            f"Slit pore surface, pore_diameter={self.config.pore_diameter:.1f} A, "
            f"{self.config.num_sheets} sheets{suffix}"
        )

    @property
    def realised_pore_type(self) -> str:
        """The pore geometry this surface actually has.

        Equal to ``config.pore_type`` unless amorphous packing failed and
        ``amorphous_fallback`` substituted another geometry. Reading the request
        instead of this is how a slit stack ends up labelled amorphous.
        """
        return self._realised_pore_type

    @property
    def packing_fell_back(self) -> bool:
        """True when the geometry built is not the geometry requested."""
        return self._packing_fell_back

    @property
    def _sheets_identical(self) -> bool:
        """
        True when every sheet is chemically the same molecule.

        Two things break identity: explicit per-sheet overrides, and pH.
        Protonation is an independent per-site draw, so sheets built at the same
        pH are samples rather than copies -- deep-copying one across the stack
        would report N times the same ionization pattern instead of an ensemble,
        understating variance and biasing net charge.

        Sheet generation, molecule naming, and .itp export all depend on this,
        so it is answered once here rather than re-derived at each site.
        """
        return self.config.sheet_overrides is None and self.config.pH is None

    def _generate_single_sheet(self, sheet_index: int) -> SheetResult:
        """Generate a single biochar sheet (or deep-copy if identical)."""
        # Identical-sheet optimisation: copy first sheet for indices > 0
        if self._sheets_identical and sheet_index > 0:
            base = self.sheets[0]
            return SheetResult(
                mol=Chem.Mol(base.mol),
                coords=base.coords.copy(),
                # deepcopy, not the dataclass itself and not a shallow copy:
                # CompositionResult carries functional_groups, placed_counts and
                # requested_counts as dicts, so anything less leaves copied
                # sheets sharing mutable state and one annotation changes them all.
                composition=copy.deepcopy(base.composition),
                atom_types=dict(base.atom_types),
                charges=dict(base.charges),
                molecule_name=base.molecule_name,
            )

        # Each sheet must be an independent protonation draw, so a pH stack
        # offsets the seed per sheet. Reusing one seed would make every sheet
        # titrate identically -- the copy problem again, just via the RNG.
        seed = self.config.seed
        if seed is not None and self.config.pH is not None:
            seed += sheet_index

        # Build GeneratorConfig for this sheet
        gen_kwargs = dict(
            target_num_carbons=self.config.target_num_carbons,
            H_C_ratio=self.config.H_C_ratio,
            O_C_ratio=self.config.O_C_ratio,
            aromaticity_percent=self.config.aromaticity_percent,
            functional_groups=self.config.functional_groups,
            pH=self.config.pH,
            defect_fraction=self.config.defect_fraction,
            seed=seed,
            strict=self.config.strict,
        )

        # Apply per-sheet overrides if present
        if self.config.sheet_overrides is not None:
            overrides = self.config.sheet_overrides[sheet_index]
            gen_kwargs.update(overrides)

        # Molecule name: same for identical sheets, indexed for distinct
        if self._sheets_identical:
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

    # ------------------------------------------------------------------
    # Amorphous packing
    # ------------------------------------------------------------------

    def _compute_amorphous_box_vectors(self) -> np.ndarray:
        """
        Compute a cubic periodic box (nm) for amorphous packing.

        The box must hold *num_sheets* arbitrarily-oriented sheets without
        forcing overlap, so it is sized from the largest sheet diameter and
        scaled by the cube root of the sheet count plus padding.
        """
        # Largest pairwise extent of any single sheet (its longest diagonal
        # bounds the sheet under any rotation).
        max_diag_A = 0.0
        for sheet in self.sheets:
            extent = sheet.coords.max(axis=0) - sheet.coords.min(axis=0)
            diag = float(np.linalg.norm(extent))
            max_diag_A = max(max_diag_A, diag)

        max_diag_nm = max_diag_A * 0.1
        # Linear scale so volume grows with sheet count; +1 sheet of margin.
        scale = (self.config.num_sheets + 1) ** (1.0 / 3.0)
        side = max_diag_nm * scale + 2 * self.config.box_padding_xy
        return np.array([side, side, side])

    def _pack_amorphous(self, box_nm: np.ndarray):
        """
        Place each sheet with a random rigid-body transform, rejecting
        placements that bring any two sheets closer than
        :attr:`SurfaceConfig.min_separation`.

        Each accepted sheet's :attr:`SheetResult.coords` is updated in place
        to its world coordinates (Ångströms).  Raises ``RuntimeError`` if a
        sheet cannot be placed within :attr:`SurfaceConfig.max_attempts`.
        """
        rng = np.random.default_rng(self.config.seed)
        box_A = box_nm * 10.0  # nm -> Angstroms
        min_sep = self.config.min_separation

        placed_coords: List[np.ndarray] = []

        for i, sheet in enumerate(self.sheets):
            # Centre the sheet at the origin so rotation is about its centroid.
            base = sheet.coords - sheet.coords.mean(axis=0)
            radius = float(np.linalg.norm(base, axis=1).max())

            placed = False
            best_sep = 0.0  # Track closest non-clashing separation achieved
            for _ in range(self.config.max_attempts):
                R = _random_rotation_matrix(rng)
                rotated = (R @ base.T).T

                # Keep the sheet fully inside the box: translate the rotated
                # centroid (origin) to a point at least `radius` from each
                # wall.  If the sheet is larger than the box, fall back to a
                # uniform translation across the whole box.
                lo = np.minimum(radius, box_A / 2.0)
                hi = np.maximum(box_A - radius, box_A / 2.0)
                translation = lo + rng.random(3) * (hi - lo)
                candidate = rotated + translation

                sep = self._min_separation(candidate, placed_coords)
                if sep > best_sep:
                    best_sep = sep
                if sep >= min_sep:
                    sheet.coords = candidate
                    placed_coords.append(candidate)
                    placed = True
                    break

            if not placed:
                raise RuntimeError(
                    f"Failed to place sheet {i + 1}/{self.config.num_sheets} "
                    f"after {self.config.max_attempts} attempts. "
                    f"Best achieved separation: {best_sep:.2f} Å "
                    f"(needed {min_sep:.2f} Å). "
                    f"Suggestions: increase box_padding_xy, reduce num_sheets, "
                    f"lower min_separation, or raise max_attempts."
                )

    @staticmethod
    def _min_separation(
        candidate: np.ndarray,
        placed: List[np.ndarray],
    ) -> float:
        """
        Return the minimum atom-atom distance (Å) between *candidate* and all
        already-placed sheets.  Returns ``inf`` when *placed* is empty (no
        constraint yet).
        """
        if not placed:
            return float("inf")
        min_d = float("inf")
        for other in placed:
            # Cheap centroid/radius pre-screen.
            c_centre = candidate.mean(axis=0)
            o_centre = other.mean(axis=0)
            c_rad = float(np.linalg.norm(candidate - c_centre, axis=1).max())
            o_rad = float(np.linalg.norm(other - o_centre, axis=1).max())
            centre_dist = float(np.linalg.norm(c_centre - o_centre))
            if centre_dist > c_rad + o_rad + min_d:
                continue
            diff = candidate[:, None, :] - other[None, :, :]
            d = float(np.sqrt((diff * diff).sum(axis=2)).min())
            if d < min_d:
                min_d = d
        return min_d

    @staticmethod
    def _is_clash_free(
        candidate: np.ndarray,
        placed: List[np.ndarray],
        min_sep: float,
    ) -> bool:
        """
        Return ``True`` if *candidate* (N, 3) keeps every atom at least
        *min_sep* Å from all atoms of already-*placed* sheets.
        """
        return SurfaceBuilder._min_separation(candidate, placed) >= min_sep


def generate_surface(
    target_num_carbons: int = 50,
    H_C_ratio: float = 0.3,
    O_C_ratio: float = 0.05,
    functional_groups: Optional[Dict[str, int]] = None,
    pH: Optional[float] = None,
    defect_fraction: float = 0.0,
    pore_diameter: float = 10.0,
    num_sheets: int = 2,
    pore_type: str = "slit",
    max_attempts: int = 500,
    min_separation: float = 3.0,
    amorphous_fallback: Optional[str] = None,
    aromaticity_percent: float = 95.0,
    box_padding_xy: float = 1.0,
    box_padding_z: float = 1.0,
    sheet_overrides: Optional[List[Dict]] = None,
    output_directory: str = ".",
    basename: str = "surface",
    system_name: str = "SLIT",
    seed: Optional[int] = None,
    strict: bool = True,
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
        pore_type: ``"slit"`` (parallel stacked sheets) or ``"amorphous"``
            (random rigid-body packing with steric rejection).
        max_attempts: Max random placement attempts per sheet for
            ``pore_type="amorphous"`` before raising ``RuntimeError``.
        min_separation: Minimum inter-sheet atom-atom distance (Å) for
            ``pore_type="amorphous"``.
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
    surface_config = SurfaceConfig(
        target_num_carbons=target_num_carbons,
        H_C_ratio=H_C_ratio,
        O_C_ratio=O_C_ratio,
        functional_groups=functional_groups,
        pH=pH,
        aromaticity_percent=aromaticity_percent,
        defect_fraction=defect_fraction,
        pore_type=pore_type,
        num_sheets=num_sheets,
        pore_diameter=pore_diameter,
        max_attempts=max_attempts,
        min_separation=min_separation,
        amorphous_fallback=amorphous_fallback,
        box_padding_xy=box_padding_xy,
        box_padding_z=box_padding_z,
        sheet_overrides=sheet_overrides,
        system_name=system_name,
        seed=seed,
        strict=strict,
    )

    builder = SurfaceBuilder(surface_config)
    sheets, box_vectors = builder.build()

    gro_path, top_path, itp_paths = builder.export_gromacs(
        output_directory=output_directory,
        basename=basename,
    )

    print("\nSurface GROMACS files written:")
    print(f"  Structure:  {gro_path}")
    print(f"  Topology:   {top_path}")
    for itp in itp_paths:
        print(f"  Include:    {itp}")
    print("\nTo run with GROMACS:")
    print(f"  gmx grompp -f md.mdp -c {gro_path.name} -p {top_path.name} -o topol.tpr")

    return sheets, gro_path, top_path, itp_paths
