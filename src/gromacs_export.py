"""
GROMACS File Export

Write biochar structures to GROMACS format files (.gro, .top, .itp).
"""

from typing import Dict, List, Optional, Tuple
from pathlib import Path
from datetime import datetime

import numpy as np
from rdkit import Chem

from .opls_typing import OPLSPropertyTable


class GROFileWriter:
    """Write GROMACS structure (.gro) files."""

    @staticmethod
    def write(
        filepath: str,
        mol: Chem.Mol,
        coords: np.ndarray,
        molecule_name: str = "BC",
        box_vectors: Optional[np.ndarray] = None,
        title: Optional[str] = None,
    ):
        """
        Write .gro structure file.

        Args:
            filepath: Output file path
            mol: RDKit molecule
            coords: 3D coordinates (in Ångströms, will be converted to nanometers)
            molecule_name: Name of molecule in GRO file (max 5 chars for .gro format)
            box_vectors: Box vectors (3x3 matrix) for periodic systems (in nm)
            title: Optional title line

        Note:
            Residue names in GROMACS .gro format are limited to 5 characters.
            Suggested naming conventions for mixed biochar simulations:
            - Temperature-based: BC400, BC600, BC800
            - Composition: BCH05 (H/C=0.5), BCO10 (O/C=0.10)
            - Sequential: BC001, BC002, BC003

            Coordinate Units:
            - Input: RDKit coordinates (Ångströms, Å)
            - Output: GROMACS coordinates (nanometers, nm)
            - Conversion: 1 Å = 0.1 nm
        """
        if title is None:
            title = f"Biochar structure generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"

        num_atoms = mol.GetNumAtoms()

        # Ensure residue name is max 5 characters (GROMACS requirement)
        residue_name = molecule_name[:5]
        if len(molecule_name) > 5:
            print(f"Warning: Residue name '{molecule_name}' exceeds 5 character limit, truncating to '{residue_name}'")

        with open(filepath, "w") as f:
            # Title line
            f.write(f"{title}\n")

            # Number of atoms
            f.write(f"{num_atoms:5d}\n")

            # Atom lines
            # Note: RDKit outputs coordinates in Ångströms, but GROMACS expects nanometers
            # Convert: 1 Å = 0.1 nm
            for atom_idx, atom in enumerate(mol.GetAtoms()):
                residue_num = 1
                atom_name = f"{atom.GetSymbol()}{atom_idx:d}"[:5]  # Truncate to 5 chars max
                atom_num = atom_idx + 1

                # Convert from Ångströms (RDKit) to nanometers (GROMACS)
                x, y, z = coords[atom_idx] * 0.1

                # Format: residue number, residue name, atom name, atom number, x, y, z
                # Columns: 1-5 (resnum), 6-10 (resname), 11-15 (atomname), 16-20 (atomnum)
                #         21-28 (x), 29-36 (y), 37-44 (z)
                # GROMACS expects coordinates in nanometers (nm)
                f.write(
                    f"{residue_num:5d}{residue_name:5s}{atom_name:5s}{atom_num:5d}"
                    f"{x:8.5f}{y:8.5f}{z:8.5f}\n"
                )

            # Box vectors
            # Note: box_vectors should be in nanometers (nm)
            if box_vectors is None:
                # Default: no periodic box
                f.write("0.0 0.0 0.0\n")
            else:
                # Write box vectors (3 numbers: v1x v2y v3z [v1y v1z v2x v2z v3x v3y])
                # For orthogonal box: just write diagonal (in nm)
                if box_vectors.shape == (3,):
                    f.write(f"{box_vectors[0]:.5f} {box_vectors[1]:.5f} {box_vectors[2]:.5f}\n")
                else:
                    # Triclinic box (in nm)
                    f.write(
                        f"{box_vectors[0,0]:.5f} {box_vectors[1,1]:.5f} "
                        f"{box_vectors[2,2]:.5f} {box_vectors[0,1]:.5f} "
                        f"{box_vectors[0,2]:.5f} {box_vectors[1,0]:.5f} "
                        f"{box_vectors[1,2]:.5f} {box_vectors[2,0]:.5f} "
                        f"{box_vectors[2,1]:.5f}\n"
                    )


class TOPFileWriter:
    """Write GROMACS topology (.top) files."""

    @staticmethod
    def write(
        filepath: str,
        mol: Chem.Mol,
        atom_types: Dict[int, str],
        charges: Dict[int, float],
        molecule_name: str = "BIOCHAR",
        forcefield_path: str = "oplsaa.ff/forcefield.itp",
        include_dihedrals: bool = True,
    ):
        """
        Write .top topology file.

        Args:
            filepath: Output file path
            mol: RDKit molecule
            atom_types: Dictionary of {atom_idx: opls_type}
            charges: Dictionary of {atom_idx: charge}
            molecule_name: Name of molecule/residue
            forcefield_path: Path to forcefield file (relative or absolute)
            include_dihedrals: If True, include dihedral section
        """
        with open(filepath, "w") as f:
            # Header
            f.write("; Generated by Biochar Generator\n")
            f.write(f"; Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            # Include forcefield
            f.write(f"#include \"{forcefield_path}\"\n\n")

            # Moleculetype section
            f.write("[ moleculetype ]\n")
            f.write(f"{molecule_name:20s} 3\n\n")

            # Atoms section
            f.write("[ atoms ]\n")
            f.write("; nr type resnr residue atom cgnr charge mass\n")

            prop_table = OPLSPropertyTable(mol, atom_types, charges)

            for prop in prop_table.get_properties():
                atom_num = prop.atom_idx + 1
                opls_type = prop.opls_type
                residue_num = 1
                residue_name = molecule_name
                atom_name = f"{prop.symbol}{prop.atom_idx}"
                cgnr = 1
                charge = prop.charge
                mass = prop.mass

                f.write(
                    f"{atom_num:5d} {opls_type:5s} {residue_num:5d} {residue_name:5s} "
                    f"{atom_name:5s} {cgnr:5d} {charge:10.6f} {mass:10.4f}\n"
                )

            # Bonds section
            f.write("\n[ bonds ]\n")
            f.write("; ai aj funct c0 c1\n")

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx() + 1
                j = bond.GetEndAtomIdx() + 1
                # GROMACS harmonic bond: funct=1
                f.write(f"{i:5d} {j:5d} 1\n")

            # Angles section
            f.write("\n[ angles ]\n")
            f.write("; ai aj ak funct c0 c1\n")

            ring_info = mol.GetRingInfo()
            ring_info.NumRings()  # Force computation

            # Extract angles from rings and bonds
            angles_written = set()
            for atom in mol.GetAtoms():
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                if len(neighbors) >= 2:
                    for i, n1 in enumerate(neighbors):
                        for n2 in neighbors[i + 1 :]:
                            angle_tuple = tuple(sorted([n1, atom.GetIdx(), n2]))
                            if angle_tuple not in angles_written:
                                f.write(f"{n1 + 1:5d} {atom.GetIdx() + 1:5d} {n2 + 1:5d} 1\n")
                                angles_written.add(angle_tuple)

            # Dihedrals section (if requested)
            if include_dihedrals:
                f.write("\n[ dihedrals ]\n")
                f.write("; ai aj ak al funct\n")

                # Simple approach: write dihedrals for all 4-atom chains
                dihedrals_written = set()
                for atom in mol.GetAtoms():
                    for neighbor in atom.GetNeighbors():
                        for third in neighbor.GetNeighbors():
                            if third.GetIdx() == atom.GetIdx():
                                continue
                            for fourth in third.GetNeighbors():
                                if fourth.GetIdx() == neighbor.GetIdx():
                                    continue
                                dihedral = (
                                    atom.GetIdx(),
                                    neighbor.GetIdx(),
                                    third.GetIdx(),
                                    fourth.GetIdx(),
                                )
                                # Only write unique dihedrals
                                dih_sorted = tuple(sorted(dihedral))
                                if dih_sorted not in dihedrals_written:
                                    f.write(
                                        f"{dihedral[0] + 1:5d} {dihedral[1] + 1:5d} "
                                        f"{dihedral[2] + 1:5d} {dihedral[3] + 1:5d} 3\n"
                                    )
                                    dihedrals_written.add(dih_sorted)

            # System section
            f.write("\n[ system ]\n")
            f.write(f"{molecule_name}\n\n")

            # Molecules section
            f.write("[ molecules ]\n")
            f.write(f"{molecule_name} 1\n")


class ITPFileWriter:
    """Write GROMACS include topology (.itp) files."""

    @staticmethod
    def write(
        filepath: str,
        mol: Chem.Mol,
        atom_types: Dict[int, str],
        charges: Dict[int, float],
        molecule_name: str = "BIOCHAR",
        include_dihedrals: bool = True,
    ):
        """
        Write .itp molecule definition file.

        Args:
            filepath: Output file path
            mol: RDKit molecule
            atom_types: Dictionary of {atom_idx: opls_type}
            charges: Dictionary of {atom_idx: charge}
            molecule_name: Name of molecule
            include_dihedrals: If True, include dihedral section
        """
        with open(filepath, "w") as f:
            # Header
            f.write("; Molecular definition for biochar\n")
            f.write(f"; Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

            # Moleculetype section
            f.write("[ moleculetype ]\n")
            f.write(f"{molecule_name:20s} 3\n\n")

            # Atoms section
            f.write("[ atoms ]\n")
            f.write("; nr type resnr residue atom cgnr charge mass\n")

            prop_table = OPLSPropertyTable(mol, atom_types, charges)

            for prop in prop_table.get_properties():
                atom_num = prop.atom_idx + 1
                opls_type = prop.opls_type
                residue_num = 1
                residue_name = molecule_name
                atom_name = f"{prop.symbol}{prop.atom_idx}"
                cgnr = 1
                charge = prop.charge
                mass = prop.mass

                f.write(
                    f"{atom_num:5d} {opls_type:5s} {residue_num:5d} {residue_name:5s} "
                    f"{atom_name:5s} {cgnr:5d} {charge:10.6f} {mass:10.4f}\n"
                )

            # Bonds section
            f.write("\n[ bonds ]\n")
            f.write("; ai aj funct c0 c1\n")

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx() + 1
                j = bond.GetEndAtomIdx() + 1
                f.write(f"{i:5d} {j:5d} 1\n")

            # Angles section
            f.write("\n[ angles ]\n")
            f.write("; ai aj ak funct c0 c1\n")

            angles_written = set()
            for atom in mol.GetAtoms():
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                if len(neighbors) >= 2:
                    for i, n1 in enumerate(neighbors):
                        for n2 in neighbors[i + 1 :]:
                            angle_tuple = tuple(sorted([n1, atom.GetIdx(), n2]))
                            if angle_tuple not in angles_written:
                                f.write(f"{n1 + 1:5d} {atom.GetIdx() + 1:5d} {n2 + 1:5d} 1\n")
                                angles_written.add(angle_tuple)

            if include_dihedrals:
                f.write("\n[ dihedrals ]\n")
                f.write("; ai aj ak al funct\n")

                dihedrals_written = set()
                for atom in mol.GetAtoms():
                    for neighbor in atom.GetNeighbors():
                        for third in neighbor.GetNeighbors():
                            if third.GetIdx() == atom.GetIdx():
                                continue
                            for fourth in third.GetNeighbors():
                                if fourth.GetIdx() == neighbor.GetIdx():
                                    continue
                                dihedral = (
                                    atom.GetIdx(),
                                    neighbor.GetIdx(),
                                    third.GetIdx(),
                                    fourth.GetIdx(),
                                )
                                dih_sorted = tuple(sorted(dihedral))
                                if dih_sorted not in dihedrals_written:
                                    f.write(
                                        f"{dihedral[0] + 1:5d} {dihedral[1] + 1:5d} "
                                        f"{dihedral[2] + 1:5d} {dihedral[3] + 1:5d} 3\n"
                                    )
                                    dihedrals_written.add(dih_sorted)


class GromacsExporter:
    """High-level interface for exporting to GROMACS files."""

    def __init__(self, output_directory: str = "."):
        self.output_dir = Path(output_directory)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def export(
        self,
        mol: Chem.Mol,
        coords: np.ndarray,
        atom_types: Dict[int, str],
        charges: Dict[int, float],
        molecule_name: str = "BC",
        basename: str = "structure",
        include_periodic_box: bool = False,
        box_size: Optional[np.ndarray] = None,
    ) -> Tuple[Path, Path, Path]:
        """
        Export molecule to GROMACS files.

        Args:
            mol: RDKit molecule
            coords: 3D coordinates
            atom_types: Atom type mapping
            charges: Charge mapping
            molecule_name: Molecule name
            basename: Base filename (without extension)
            include_periodic_box: If True, create periodic box
            box_size: Box dimensions [Lx, Ly, Lz]

        Returns:
            (gro_path, top_path, itp_path)
        """
        # Determine box vectors
        box_vectors = None
        if include_periodic_box:
            if box_size is None:
                # Auto-determine box size from coordinates
                box_size = self._calculate_box_size(coords)
            box_vectors = np.diag(box_size)

        # Write files
        gro_path = self.output_dir / f"{basename}.gro"
        GROFileWriter.write(
            str(gro_path),
            mol,
            coords,
            molecule_name=molecule_name,
            box_vectors=box_vectors,
        )

        top_path = self.output_dir / f"{basename}.top"
        TOPFileWriter.write(
            str(top_path),
            mol,
            atom_types,
            charges,
            molecule_name=molecule_name,
        )

        itp_path = self.output_dir / f"{basename}.itp"
        ITPFileWriter.write(
            str(itp_path),
            mol,
            atom_types,
            charges,
            molecule_name=molecule_name,
        )

        return gro_path, top_path, itp_path

    @staticmethod
    def _calculate_box_size(coords: np.ndarray, padding: float = 1.0) -> np.ndarray:
        """
        Calculate appropriate box size based on coordinates.

        Args:
            coords: 3D coordinates (expected in Ångströms from RDKit)
            padding: Padding around molecule (nm)

        Returns:
            Box size [Lx, Ly, Lz] in nanometers (nm)
        """
        # Coordinates are in Ångströms, convert to nanometers
        coords_nm = coords * 0.1
        min_coords = coords_nm.min(axis=0)
        max_coords = coords_nm.max(axis=0)
        extent = max_coords - min_coords
        box_size = extent + 2 * padding
        return box_size
