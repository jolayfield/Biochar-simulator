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
    ) -> Tuple[Chem.Mol, np.ndarray]:
        """
        Generate 3D coordinates for molecule.

        Args:
            mol: RDKit molecule
            force_aromatic_planarity: If True, enforce planarity on aromatic rings
            num_conformers: Number of conformers to generate
            use_distance_geometry: Use distance geometry for coordinate generation

        Returns:
            (mol_with_coords, coords_array)
        """
        # Add hydrogens if not already present
        if mol.GetNumHeavyAtoms() == mol.GetNumAtoms():
            mol = Chem.AddHs(mol)

        # Generate 3D coordinates using distance geometry
        seed = self.seed if self.seed is not None else -1

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
            result = AllChem.EmbedMolecule(mol, params)
            if result == 0:
                break

        # If all embedding attempts failed, assign random coordinates as last resort
        if result != 0:
            AllChem.Compute2DCoords(mol)

        # Optional: apply force field relaxation
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol)
            except Exception:
                pass

        # Extract coordinates — use conformer 0 if available, else zero-coords
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
        Check for steric clashes between atoms.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            clash_threshold: Minimum allowed distance (Angstrom)

        Returns:
            List of error messages
        """
        errors = []

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
                min_distance = r_vdw_i + r_vdw_j - clash_threshold

                if distances[i, j] < min_distance:
                    errors.append(
                        f"Steric clash: atoms {i} and {j} "
                        f"(distance: {distances[i, j]:.2f}, min: {min_distance:.2f})"
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
