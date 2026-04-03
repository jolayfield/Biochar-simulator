"""
3D Geometry Generation and Validation

Generate 3D coordinates for biochar molecules and validate geometry.
"""

from typing import Tuple, List, Optional
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial.distance import pdist, squareform

from .constants import VDW_RADII, COVALENT_RADII


class CoordinateGenerator:
    """Generate valid 3D coordinates for molecules."""

    def __init__(self, seed: int = None):
        self.seed = seed
        if seed is not None:
            np.random.seed(seed)

    def generate_3d_coordinates(
        self,
        mol: Chem.Mol,
        force_aromatic_planarity: bool = True,
        num_conformers: int = 1,
        use_distance_geometry: bool = True,
        max_embedding_attempts: int = 3,
        max_ff_iterations: int = None,
    ) -> Tuple[Chem.Mol, np.ndarray]:
        """
        Generate 3D coordinates for molecule with enhanced relaxation.

        Args:
            mol: RDKit molecule
            force_aromatic_planarity: If True, enforce planarity on aromatic rings
            num_conformers: Number of conformers to generate
            use_distance_geometry: Use distance geometry for coordinate generation
            max_embedding_attempts: Number of independent embedding attempts (default 3)
            max_ff_iterations: Max force field iterations (default: 500 for large molecules, 300 for small)

        Returns:
            (mol_with_coords, coords_array)
        """
        # Add hydrogens if not already present
        if mol.GetNumHeavyAtoms() == mol.GetNumAtoms():
            mol = Chem.AddHs(mol)

        # Adaptive force field iterations based on molecule size
        if max_ff_iterations is None:
            num_heavy_atoms = mol.GetNumHeavyAtoms()
            if num_heavy_atoms > 300:
                max_ff_iterations = 500
            elif num_heavy_atoms > 100:
                max_ff_iterations = 400
            else:
                max_ff_iterations = 300

        # Try multiple independent embeddings and keep the best one
        best_mol = None
        best_energy = float('inf')
        best_converged = False

        for attempt in range(max_embedding_attempts):
            # Use different seed for each attempt to escape local minima
            seed = (self.seed + attempt) if self.seed is not None else attempt

            # Create a copy for this attempt
            mol_copy = Chem.RWMol(mol)

            # Try ETKDGv3 first (best quality), fall back to ETKDGv2, then basic
            result = -1
            for params_fn in [
                lambda: AllChem.ETKDGv3(),
                lambda: AllChem.ETKDGv2(),
                lambda: AllChem.EmbedParameters(),
            ]:
                params = params_fn()
                params.randomSeed = seed
                params.maxIterations = 200
                result = AllChem.EmbedMolecule(mol_copy, params)
                if result == 0:
                    break

            # If embedding failed, skip this attempt
            if result != 0:
                continue

            # Apply force field relaxation with enhanced iterations
            ff_converged = False
            ff_energy = float('inf')

            try:
                # Try MMFF94 with tighter convergence
                ff = AllChem.MMFFGetMoleculeForceField(mol_copy)
                if ff is not None:
                    minimize_result = ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-6)
                    ff_converged = minimize_result == 0
                    ff_energy = ff.CalcEnergy()
            except Exception:
                pass

            # Fallback to UFF if MMFF94 unavailable
            if ff_energy == float('inf'):
                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol_copy)
                    if ff is not None:
                        minimize_result = ff.Minimize(maxIts=max_ff_iterations, forceTol=1e-6)
                        ff_converged = minimize_result == 0
                        ff_energy = ff.CalcEnergy()
                except Exception:
                    pass

            # Keep track of best embedding (lowest energy)
            if ff_energy < best_energy:
                best_energy = ff_energy
                best_mol = mol_copy
                best_converged = ff_converged

        # Use best embedding, or fall back to 2D if all attempts failed
        if best_mol is not None:
            mol = best_mol
        else:
            # Fallback: use 2D coordinates
            AllChem.Compute2DCoords(mol)

        # Extract coordinates
        try:
            conf = mol.GetConformer(0)
            coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        except Exception:
            coords = np.zeros((mol.GetNumAtoms(), 3))

        # Apply planarity constraints to aromatic rings if requested
        if force_aromatic_planarity:
            coords = self._enforce_aromatic_planarity(mol, coords)

        return mol, coords

    def _enforce_aromatic_planarity(
        self, mol: Chem.Mol, coords: np.ndarray, tolerance: float = 0.5
    ) -> np.ndarray:
        """
        Enforce planarity on aromatic ring systems.

        Args:
            mol: RDKit molecule
            coords: Current 3D coordinates
            tolerance: Maximum deviation from plane (Angstrom)

        Returns:
            Adjusted coordinates
        """
        # Find aromatic rings
        ring_info = mol.GetRingInfo()
        aromatic_rings = [
            ring for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        adjusted_coords = coords.copy()

        for ring_atoms in aromatic_rings:
            if len(ring_atoms) < 3:
                continue

            # Get ring atom coordinates
            ring_coords = coords[list(ring_atoms)]

            # Fit plane to ring atoms
            plane_normal, plane_point = self._fit_plane(ring_coords)

            # Project ring atoms onto plane
            for i, atom_idx in enumerate(ring_atoms):
                point = adjusted_coords[atom_idx]
                projected = point - np.dot(point - plane_point, plane_normal) * plane_normal
                adjusted_coords[atom_idx] = projected

        return adjusted_coords

    @staticmethod
    def _fit_plane(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit a plane to 3D points using PCA.

        Returns:
            (plane_normal, plane_point)
        """
        center = points.mean(axis=0)
        centered_points = points - center

        # SVD
        u, s, vt = np.linalg.svd(centered_points)
        normal = vt[2]  # Normal is the smallest singular vector

        return normal, center

    def validate_and_relax(
        self, mol: Chem.Mol, coords: np.ndarray, max_iterations: int = 100
    ) -> Tuple[np.ndarray, bool]:
        """
        Validate and relax coordinates using force field.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            max_iterations: Maximum relaxation iterations

        Returns:
            (relaxed_coords, converged)
        """
        # Set coordinates in molecule
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, coord in enumerate(coords):
            conf.SetAtomPosition(i, tuple(coord))
        mol.RemoveAllConformers()
        mol.AddConformer(conf)

        # Try MMFF94 force field
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol)
            if ff is not None:
                result = ff.Minimize(maxIts=max_iterations)
                converged = result == 0
                conf = mol.GetConformer(0)
                relaxed_coords = np.array([
                    list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())
                ])
                return relaxed_coords, converged
        except Exception:
            pass

        # Fallback to UFF
        try:
            ff = AllChem.UFFGetMoleculeForceField(mol)
            if ff is not None:
                result = ff.Minimize(maxIts=max_iterations)
                converged = result == 0
                conf = mol.GetConformer(0)
                relaxed_coords = np.array([
                    list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())
                ])
                return relaxed_coords, converged
        except Exception:
            pass

        return coords, False


class ClashResolver:
    """Resolve steric clashes in molecular geometries."""

    @staticmethod
    def resolve_clashes(
        mol: Chem.Mol,
        coords: np.ndarray,
        max_iterations: int = 10,
        displacement_step: float = 0.1,
        use_vdw_radii: bool = True,
    ) -> np.ndarray:
        """
        Iteratively resolve steric clashes by displacing atoms away from each other.

        Algorithm:
        - Detect pairs of atoms that are too close (violating Van der Waals radii)
        - Displace the lighter atom (usually H) away from its neighbor
        - Repeat until no clashes remain or max iterations reached

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            max_iterations: Maximum number of clash resolution iterations
            displacement_step: Step size for atom displacement (Angstrom)
            use_vdw_radii: If True, use Van der Waals radii; else use fixed threshold

        Returns:
            Adjusted coordinates with reduced clashes
        """
        adjusted_coords = coords.copy()

        for iteration in range(max_iterations):
            clashes = ClashResolver._detect_clashes(mol, adjusted_coords, use_vdw_radii)

            if not clashes:
                # No more clashes
                break

            # Displace clashing atoms
            for atom_i, atom_j in clashes:
                # Prefer to move the lighter atom (usually H)
                mass_i = mol.GetAtomWithIdx(atom_i).GetMass()
                mass_j = mol.GetAtomWithIdx(atom_j).GetMass()

                if mass_i < mass_j:
                    moving_atom = atom_i
                else:
                    moving_atom = atom_j

                # Calculate displacement direction (away from the other atom)
                if moving_atom == atom_i:
                    other_atom = atom_j
                else:
                    other_atom = atom_i

                direction = adjusted_coords[moving_atom] - adjusted_coords[other_atom]
                distance = np.linalg.norm(direction)

                if distance > 1e-6:  # Avoid division by zero
                    direction_normalized = direction / distance
                    adjusted_coords[moving_atom] += displacement_step * direction_normalized

        return adjusted_coords

    @staticmethod
    def _detect_clashes(
        mol: Chem.Mol,
        coords: np.ndarray,
        use_vdw_radii: bool = True,
    ) -> List[Tuple[int, int]]:
        """
        Detect steric clashes between atoms.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            use_vdw_radii: If True, use Van der Waals radii; else use fixed 2.0 Å threshold

        Returns:
            List of (atom_i, atom_j) pairs with clashes
        """
        clashes = []
        distances = squareform(pdist(coords))

        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                # Skip bonded atoms
                if mol.GetBondBetweenAtoms(i, j) is not None:
                    continue

                if use_vdw_radii:
                    # Use Van der Waals radii with 0.2 Å buffer
                    atom_i = mol.GetAtomWithIdx(i)
                    atom_j = mol.GetAtomWithIdx(j)
                    symbol_i = atom_i.GetSymbol()
                    symbol_j = atom_j.GetSymbol()
                    r_vdw_i = VDW_RADII.get(symbol_i, 1.70)
                    r_vdw_j = VDW_RADII.get(symbol_j, 1.70)
                    min_distance = r_vdw_i + r_vdw_j + 0.2
                else:
                    # Use fixed threshold
                    min_distance = 2.0

                if distances[i, j] < min_distance:
                    clashes.append((i, j))

        return clashes


class GeometryValidator:
    """Validate 3D molecular geometry."""

    @staticmethod
    def validate_geometry(
        mol: Chem.Mol, coords: np.ndarray
    ) -> Tuple[bool, List[str]]:
        """
        Validate 3D geometry for issues.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            (is_valid, error_messages)
        """
        errors = []

        # Check for NaN coordinates
        if np.any(np.isnan(coords)):
            errors.append("NaN values found in coordinates")
            return False, errors

        # Check for steric clashes
        clash_errors = GeometryValidator._check_steric_clashes(mol, coords)
        errors.extend(clash_errors)

        # Check bond lengths
        bond_errors = GeometryValidator._validate_bond_lengths(mol, coords)
        errors.extend(bond_errors)

        # Check for unphysical geometries
        if np.max(coords) > 1000 or np.min(coords) < -1000:
            errors.append("Coordinates are unphysically large")

        return len(errors) == 0, errors

    @staticmethod
    def _check_steric_clashes(
        mol: Chem.Mol, coords: np.ndarray, clash_threshold: float = 2.0
    ) -> List[str]:
        """
        Check for steric clashes between atoms with enhanced reporting.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            clash_threshold: Minimum allowed distance (Angstrom)

        Returns:
            List of error messages with clash severity
        """
        errors = []
        clashes = []

        # Calculate pairwise distances
        distances = squareform(pdist(coords))

        for i in range(mol.GetNumAtoms()):
            for j in range(i + 1, mol.GetNumAtoms()):
                # Skip bonded atoms
                if mol.GetBondBetweenAtoms(i, j) is not None:
                    continue

                atom_i = mol.GetAtomWithIdx(i)
                atom_j = mol.GetAtomWithIdx(j)

                # Get Van der Waals radii
                symbol_i = atom_i.GetSymbol()
                symbol_j = atom_j.GetSymbol()
                r_vdw_i = VDW_RADII.get(symbol_i, 1.70)
                r_vdw_j = VDW_RADII.get(symbol_j, 1.70)
                min_distance = r_vdw_i + r_vdw_j + 0.2  # Use realistic Van der Waals sum

                distance = distances[i, j]
                if distance < min_distance:
                    severity = min_distance - distance  # How far below minimum
                    clash_type = "H-H" if symbol_i == "H" and symbol_j == "H" else \
                                 "H-C" if symbol_i in ["H", "C"] and symbol_j in ["H", "C"] else \
                                 "Other"

                    clashes.append((i, j, distance, min_distance, severity, clash_type))
                    errors.append(
                        f"Steric clash: atoms {i} and {j} "
                        f"(distance: {distance:.2f}, min: {min_distance:.2f}, "
                        f"severity: {severity:.2f} Å, type: {clash_type})"
                    )

        return errors

    @staticmethod
    def _validate_bond_lengths(
        mol: Chem.Mol, coords: np.ndarray
    ) -> List[str]:
        """
        Validate bond lengths are reasonable.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            List of error messages
        """
        errors = []

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            atom_i = mol.GetAtomWithIdx(i)
            atom_j = mol.GetAtomWithIdx(j)

            symbol_i = atom_i.GetSymbol()
            symbol_j = atom_j.GetSymbol()

            # Calculate distance
            distance = np.linalg.norm(coords[i] - coords[j])

            # Expected distance: sum of covalent radii
            r_cov_i = COVALENT_RADII.get(symbol_i, 0.76)
            r_cov_j = COVALENT_RADII.get(symbol_j, 0.76)
            expected_distance = r_cov_i + r_cov_j

            # Allow 20% deviation
            min_dist = expected_distance * 0.8
            max_dist = expected_distance * 1.5

            if distance < min_dist or distance > max_dist:
                errors.append(
                    f"Unusual bond length: atoms {i}-{j} "
                    f"(distance: {distance:.2f}, expected: {expected_distance:.2f})"
                )

        return errors

    @staticmethod
    def check_aromaticity(mol: Chem.Mol) -> Tuple[bool, int]:
        """
        Check if aromaticity is preserved.

        Returns:
            (is_valid, num_aromatic_atoms)
        """
        num_aromatic = sum(
            1 for atom in mol.GetAtoms() if atom.GetIsAromatic()
        )
        return num_aromatic > 0, num_aromatic

    @staticmethod
    def measure_ring_planarity(mol: Chem.Mol, coords: np.ndarray) -> Tuple[float, str]:
        """
        Measure planarity of aromatic rings.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates

        Returns:
            (avg_deviation, assessment)
        """
        ring_info = mol.GetRingInfo()
        aromatic_rings = [
            ring for ring in ring_info.AtomRings()
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        ]

        if not aromatic_rings:
            return 0.0, "No aromatic rings"

        deviations = []

        for ring_atoms in aromatic_rings:
            if len(ring_atoms) < 3:
                continue

            ring_coords = coords[list(ring_atoms)]
            normal, point = CoordinateGenerator._fit_plane(ring_coords)

            # Calculate deviations from plane
            for coord in ring_coords:
                deviation = abs(np.dot(coord - point, normal))
                deviations.append(deviation)

        avg_deviation = np.mean(deviations) if deviations else 0.0
        if avg_deviation < 0.1:
            assessment = "Good planarity"
        elif avg_deviation < 0.3:
            assessment = "Acceptable planarity"
        else:
            assessment = "Poor planarity"

        return avg_deviation, assessment
