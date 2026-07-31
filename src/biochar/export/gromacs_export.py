"""
GROMACS File Export

Write biochar structures to GROMACS format files (.gro, .top, .itp).
"""

from typing import Dict, Optional, Tuple
from pathlib import Path
from datetime import datetime

import numpy as np
from rdkit import Chem

from ..pipeline.opls_typing import OPLSPropertyTable
from ..constants import GROMACS_OPLS_TYPE_MAP, SUPPLEMENTARY_ANGLE_PARAMS


def _angle_suffix(atom_types: Dict[int, str], i: int, j: int, k: int) -> str:
    """Explicit parameters for an angle stock oplsaa.ff cannot resolve, else "".

    Angles normally carry funct only and are resolved from the #included
    forcefield by type name. A handful of combinations have no entry there at all
    (see SUPPLEMENTARY_ANGLE_PARAMS); for those, grompp would fail with "No
    default Angle types" unless the parameters are written inline.
    """
    outer = sorted([atom_types.get(i, ""), atom_types.get(k, "")])
    key = (outer[0], atom_types.get(j, ""), outer[1])
    params = SUPPLEMENTARY_ANGLE_PARAMS.get(key)
    if params is None:
        return ""
    theta0, force_constant = params
    return f" {theta0:.3f} {force_constant:.3f}"


_GRO_ATOM_NAME_WIDTH = 5
_BASE36 = "0123456789abcdefghijklmnopqrstuvwxyz"

# OPLS-AA keeps sp2 centres planar with an improper torsion, supplied as a
# #define in ffbonded.itp rather than as a [ dihedraltypes ] entry. The .rtp
# [ bondedtypes ] row declares improper function type 1 -- the periodic *proper*
# form -- because, as ffbonded.itp puts it, these are "implemented as proper
# dihedrals [1+cos(2*x+180)] for the moment, to keep things compatible". The
# macro name is written where parameters would otherwise go; the preprocessor
# substitutes it when the forcefield is #included.
_AROMATIC_IMPROPER = "improper_Z_CA_X_Y"
_IMPROPER_FUNCT = 1


def _atom_name(symbol: str, idx: int) -> str:
    """A per-atom name that fits the .gro field and still identifies one atom.

    The .gro atom-name column is five characters. Truncating a longer name is not
    an option: ``C10000`` and ``C1000`` would collapse onto the same name, so one
    name would identify several atoms *and* stop matching the .itp, which does not
    truncate. Both files call this, so they agree by construction.

    Decimal is used while it fits, which keeps the familiar C0/C1/H12 names for
    everything below the width limit. Past that the index is written in base 36,
    buying 36^4 = 1.68M single-letter-element atoms in the same five columns.
    """
    budget = _GRO_ATOM_NAME_WIDTH - len(symbol)
    if budget < 1:
        raise ValueError(f"element symbol {symbol!r} leaves no room for an index")
    if len(str(idx)) <= budget:
        return f"{symbol}{idx}"

    encoded, n = "", idx
    while n:
        encoded = _BASE36[n % 36] + encoded
        n //= 36
    encoded = encoded or "0"
    if len(encoded) > budget:
        raise ValueError(
            f"atom index {idx} does not fit {budget} character(s) after {symbol!r}; "
            f"the .gro atom-name field cannot identify this many atoms uniquely"
        )
    return f"{symbol}{encoded}"


def _pairs_1_4(mol: Chem.Mol):
    """Atom index pairs exactly three bonds apart, 0-based and ordered.

    nrexcl = 3 removes these from the plain non-bonded loop. OPLS-AA does not
    discard them -- it restores them at half strength, and [ defaults ]
    gen-pairs=yes generates the *parameters* for the pairs a topology lists. It
    does not generate the list. Omitting [ pairs ] therefore silently drops every
    1-4 interaction instead of rescaling it.

    Walked over the bond graph rather than read off Chem.GetDistanceMatrix: that
    matrix is N x N, so a 10,000-atom sheet would allocate 100M entries and scan
    50M pairs to find a list that is linear in the number of torsions.

    Candidates come from torsion paths i-j-k-l, then any that are also 1-2 or 1-3
    bonded are dropped. In a fused ring the two ends of a torsion can be closer
    than three bonds by another route, and such a pair is already excluded by
    nrexcl -- listing it would double-count.
    """
    adjacency = [
        {n.GetIdx() for n in atom.GetNeighbors()} for atom in mol.GetAtoms()
    ]
    # 1-2 and 1-3 partners of each atom, by the shortest route available.
    near = []
    for i, direct in enumerate(adjacency):
        two_away = set().union(*(adjacency[n] for n in direct)) if direct else set()
        near.append((direct | two_away) - {i})

    pairs = set()
    for bond in mol.GetBonds():
        j, k = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        for i in adjacency[j] - {k}:
            for l_ in adjacency[k] - {j}:
                if i == l_ or l_ in near[i]:
                    continue
                pairs.add((min(i, l_), max(i, l_)))
    return sorted(pairs)


def _aromatic_impropers(mol: Chem.Mol):
    """Improper torsions holding aromatic ring carbons in their ring plane.

    Proper dihedrals restrain rotation about a bond; nothing in them stops a
    substituted ring carbon pyramidalising out of plane. Without these a sheet is
    planar only because it was written that way.

    Atom order follows the convention in oplsaa.ff's own residue templates
    (``CG CE2 CD2 HD2`` for a phenylalanine ring CH): the central atom sits
    third, its three neighbours occupy the other positions, and the substituent
    -- the neighbour outside the ring where there is one -- goes last.
    """
    out = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C" or not atom.GetIsAromatic() or not atom.IsInRing():
            continue
        neighbours = list(atom.GetNeighbors())
        if len(neighbours) != 3:
            continue
        substituent = next(
            (n for n in neighbours if not (n.GetIsAromatic() and n.IsInRing())),
            neighbours[-1],
        )
        others = [n.GetIdx() for n in neighbours if n.GetIdx() != substituent.GetIdx()]
        out.append((others[0], others[1], atom.GetIdx(), substituent.GetIdx()))
    return out


def _write_pairs(f, mol: Chem.Mol) -> None:
    """Emit [ pairs ], the 1-4 list nrexcl excludes and OPLS-AA rescales."""
    f.write("\n[ pairs ]\n")
    f.write("; ai aj funct   ; 1-4 interactions, rescaled by fudgeLJ/fudgeQQ\n")
    for i, j in _pairs_1_4(mol):
        f.write(f"{i + 1:5d} {j + 1:5d} 1\n")


def _write_impropers(f, mol: Chem.Mol) -> None:
    """Emit the improper torsions that hold aromatic centres planar.

    A separate [ dihedrals ] block from the proper torsions, which is how
    pdb2gmx writes them too.
    """
    impropers = _aromatic_impropers(mol)
    if not impropers:
        return
    f.write("\n[ dihedrals ]\n")
    f.write("; ai aj ak al funct   ; impropers: keep aromatic centres planar\n")
    for i, j, k, l_ in impropers:
        f.write(
            f"{i + 1:5d} {j + 1:5d} {k + 1:5d} {l_ + 1:5d} "
            f"{_IMPROPER_FUNCT}    {_AROMATIC_IMPROPER}\n"
        )


def _write_qtot_footer(f, mol: Chem.Mol, qtot: float) -> None:
    """
    Record the molecule's net charge after an ``[ atoms ]`` block.

    A non-zero qtot is legitimate -- a deprotonated carboxyl really does carry
    -1 -- so this states the charge rather than warning about it, and names the
    step responsible for balancing it. Without the note, a reader hitting
    ``qtot -3.0000`` has no way to tell a real anion from a topology bug.
    """
    formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    f.write(f"; total charge: {qtot:+.4f} e (formal {formal:+d} e)\n")
    if formal != 0:
        f.write(
            "; NOTE: this molecule is intentionally charged. Neutralise the "
            "system with `gmx genion -neutral` at solvation.\n"
        )


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

        # Convert coordinates from Ångströms (RDKit) to nanometers (GROMACS)
        coords_nm = coords * 0.1

        # Compute a valid bounding box and shift so the molecule sits inside it.
        # A zero-size box (0.0 0.0 0.0) is invalid for VMD/PyMOL/GROMACS — they can't
        # draw the periodic cell or run minimization. When no box is supplied, build
        # an orthogonal box around the molecule with 1 nm padding and shift the
        # coordinates so all atoms are at positive (in-box) positions.
        if box_vectors is None:
            mins = coords_nm.min(axis=0)
            maxs = coords_nm.max(axis=0)
            extent = maxs - mins
            padding = 1.0  # nm on each side
            box = np.maximum(extent + 2.0 * padding, 1.0)
            shift = padding - mins  # places mins at +padding inside the box
            coords_nm = coords_nm + shift
            box_line = f"{box[0]:.5f} {box[1]:.5f} {box[2]:.5f}\n"
        else:
            if box_vectors.shape == (3,):
                box_line = (
                    f"{box_vectors[0]:.5f} {box_vectors[1]:.5f} "
                    f"{box_vectors[2]:.5f}\n"
                )
            else:
                box_line = (
                    f"{box_vectors[0,0]:.5f} {box_vectors[1,1]:.5f} "
                    f"{box_vectors[2,2]:.5f} {box_vectors[0,1]:.5f} "
                    f"{box_vectors[0,2]:.5f} {box_vectors[1,0]:.5f} "
                    f"{box_vectors[1,2]:.5f} {box_vectors[2,0]:.5f} "
                    f"{box_vectors[2,1]:.5f}\n"
                )

        with open(filepath, "w") as f:
            # Title line
            f.write(f"{title}\n")

            # Number of atoms
            f.write(f"{num_atoms:5d}\n")

            # Atom lines
            for atom_idx, atom in enumerate(mol.GetAtoms()):
                residue_num = 1
                atom_name = _atom_name(atom.GetSymbol(), atom_idx)
                atom_num = atom_idx + 1

                x, y, z = coords_nm[atom_idx]

                # Format: residue number, residue name, atom name, atom number, x, y, z
                # Columns: 1-5 (resnum), 6-10 (resname), 11-15 (atomname), 16-20 (atomnum)
                #         21-28 (x), 29-36 (y), 37-44 (z)
                # Note: GROMACS format requires residue_name LEFT-aligned, atom_name RIGHT-aligned
                # GROMACS expects coordinates in nanometers (nm) with 3 decimal places
                f.write(
                    f"{residue_num:5d}{residue_name:<5s}{atom_name:>5s}{atom_num:5d}"
                    f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
                )

            # Box vectors line (nm). Always emit a valid non-zero box.
            f.write(box_line)


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

            # Include forcefield. This must resolve to a real oplsaa.ff: every atom
            # type, LJ and bonded parameter is looked up there by the opls_XXX name
            # written below, and [ bonds ]/[ angles ] here carry no parameters of
            # their own. grompp fails loudly if it cannot find the include.
            f.write(f"#include \"{forcefield_path}\"\n\n")

            # Moleculetype section
            f.write("[ moleculetype ]\n")
            f.write(f"{molecule_name:20s} 3\n\n")

            # Atoms section
            f.write("[ atoms ]\n")
            f.write("; nr type resnr residue atom cgnr charge mass\n")

            prop_table = OPLSPropertyTable(mol, atom_types, charges)

            qtot = 0.0
            for prop in prop_table.get_properties():
                atom_num = prop.atom_idx + 1
                opls_type = GROMACS_OPLS_TYPE_MAP.get(prop.opls_type, prop.opls_type)
                residue_num = 1
                residue_name = molecule_name
                atom_name = _atom_name(prop.symbol, prop.atom_idx)
                cgnr = 1
                charge = prop.charge
                mass = prop.mass
                qtot += charge

                f.write(
                    f"{atom_num:5d} {opls_type:5s} {residue_num:5d} {residue_name:5s} "
                    f"{atom_name:5s} {cgnr:5d} {charge:10.6f} {mass:10.4f}"
                    f"   ; qtot {qtot:+.4f}\n"
                )

            _write_qtot_footer(f, mol, qtot)

            # Bonds section
            f.write("\n[ bonds ]\n")
            f.write("; ai aj funct c0 c1\n")

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx() + 1
                j = bond.GetEndAtomIdx() + 1
                # GROMACS harmonic bond: funct=1
                f.write(f"{i:5d} {j:5d} 1\n")

            _write_pairs(f, mol)

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
                                suffix = _angle_suffix(atom_types, n1, atom.GetIdx(), n2)
                                f.write(
                                    f"{n1 + 1:5d} {atom.GetIdx() + 1:5d} {n2 + 1:5d} 1{suffix}\n"
                                )
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

                _write_impropers(f, mol)

            # System section
            f.write("\n[ system ]\n")
            f.write(f"{molecule_name}\n\n")

            # Molecules section
            f.write("[ molecules ]\n")
            f.write(f"{molecule_name} 1\n")


class ITPFileWriter:
    """Write GROMACS include topology (.itp) files.

    Note: .itp files are typically included in a .top file that defines atom types.
    We don't include [ atomtypes ] here to avoid duplication.
    """

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

        Note:
            The .itp file includes atom type definitions for standalone use.
            When included in a .top file, the atom types will be shared.
        """
        with open(filepath, "w") as f:
            # Header
            f.write("; Molecular definition for biochar\n")
            f.write(f"; Created: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("; Include this file in a .top file that defines atom types\n\n")

            # Note: We don't include [ atomtypes ] here. The .itp file is meant to be
            # included in a .top file that already defines all atom types via #include.

            # Moleculetype section
            f.write("[ moleculetype ]\n")
            f.write(f"{molecule_name:20s} 3\n\n")

            # Atoms section
            f.write("[ atoms ]\n")
            f.write("; nr type resnr residue atom cgnr charge mass\n")

            prop_table = OPLSPropertyTable(mol, atom_types, charges)

            qtot = 0.0
            for prop in prop_table.get_properties():
                atom_num = prop.atom_idx + 1
                opls_type = GROMACS_OPLS_TYPE_MAP.get(prop.opls_type, prop.opls_type)
                residue_num = 1
                residue_name = molecule_name
                atom_name = _atom_name(prop.symbol, prop.atom_idx)
                cgnr = 1
                charge = prop.charge
                mass = prop.mass
                qtot += charge

                f.write(
                    f"{atom_num:5d} {opls_type:5s} {residue_num:5d} {residue_name:5s} "
                    f"{atom_name:5s} {cgnr:5d} {charge:10.6f} {mass:10.4f}"
                    f"   ; qtot {qtot:+.4f}\n"
                )

            _write_qtot_footer(f, mol, qtot)

            # Bonds section
            f.write("\n[ bonds ]\n")
            f.write("; ai aj funct c0 c1\n")

            for bond in mol.GetBonds():
                i = bond.GetBeginAtomIdx() + 1
                j = bond.GetEndAtomIdx() + 1
                f.write(f"{i:5d} {j:5d} 1\n")

            _write_pairs(f, mol)

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
                                suffix = _angle_suffix(atom_types, n1, atom.GetIdx(), n2)
                                f.write(
                                    f"{n1 + 1:5d} {atom.GetIdx() + 1:5d} {n2 + 1:5d} 1{suffix}\n"
                                )
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

                _write_impropers(f, mol)


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


class MultiSheetGROWriter:
    """Write a multi-residue .gro file for surface systems.

    Each sheet occupies a separate residue with an incrementing residue
    number.  Atom numbering is global and contiguous across all sheets.
    """

    @staticmethod
    def write(
        filepath: str,
        sheets: list,               # List[SheetResult] — avoid circular import
        box_vectors: np.ndarray,    # 1-D (3,) in nm, or 3×3 in nm
        title: Optional[str] = None,
    ):
        """
        Write a single .gro containing atoms from multiple sheets.

        Args:
            filepath: Output path.
            sheets: List of SheetResult objects.  Each must expose
                ``.mol`` (RDKit Mol), ``.coords`` (Å, shape N×3),
                and ``.molecule_name`` (str ≤ 5 chars).
            box_vectors: Periodic box vectors in nm (orthogonal 1-D or
                full 3×3 triclinic).
            title: Optional title line.
        """
        if title is None:
            title = (
                f"Multi-sheet surface — "
                f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
            )

        total_atoms = sum(s.mol.GetNumAtoms() for s in sheets)

        with open(filepath, "w") as f:
            f.write(f"{title}\n")
            f.write(f"{total_atoms:5d}\n")

            global_atom_num = 0
            for residue_num, sheet in enumerate(sheets, start=1):
                residue_name = sheet.molecule_name[:5]
                for atom_idx, atom in enumerate(sheet.mol.GetAtoms()):
                    global_atom_num += 1
                    atom_name = _atom_name(atom.GetSymbol(), atom_idx)
                    x_nm, y_nm, z_nm = sheet.coords[atom_idx] * 0.1  # Å → nm
                    # Format: residue_name LEFT-aligned, atom_name RIGHT-aligned
                    f.write(
                        f"{residue_num % 100000:5d}"
                        f"{residue_name:<5s}"
                        f"{atom_name:>5s}"
                        f"{global_atom_num % 100000:5d}"
                        f"{x_nm:8.3f}{y_nm:8.3f}{z_nm:8.3f}\n"
                    )

            # Box vectors line
            bv = np.asarray(box_vectors)
            if bv.ndim == 1 and bv.shape[0] == 3:
                f.write(
                    f"{bv[0]:.5f} {bv[1]:.5f} {bv[2]:.5f}\n"
                )
            else:
                # Full 3×3 triclinic in GROMACS 9-value order:
                # v1x v2y v3z v1y v1z v2x v2z v3x v3y
                f.write(
                    f"{bv[0,0]:.5f} {bv[1,1]:.5f} {bv[2,2]:.5f} "
                    f"{bv[0,1]:.5f} {bv[0,2]:.5f} {bv[1,0]:.5f} "
                    f"{bv[1,2]:.5f} {bv[2,0]:.5f} {bv[2,1]:.5f}\n"
                )


class SurfaceTopologyWriter:
    """Write .top files for multi-sheet surface systems."""

    @staticmethod
    def write(
        filepath: str,
        sheets: list,              # List[SheetResult]
        itp_paths: list,           # List[Path]
        sheets_identical: bool,
        system_name: str = "Slit Pore Surface",
        forcefield_path: str = "oplsaa.ff/forcefield.itp",
    ):
        """
        Write a combined .top for a surface system.

        When *sheets_identical* is ``True``:
          - One ``#include`` for the single .itp.
          - ``[ molecules ]`` lists the molecule name once with count = N.

        When *sheets_identical* is ``False``:
          - One ``#include`` per distinct .itp.
          - ``[ molecules ]`` lists each molecule name once with count = 1.
        """
        # Total system charge. Each sheet's own qtot lives in its .itp, but the
        # stack's charge is what genion has to balance, and reading it off the
        # .itp means multiplying by the [molecules] count by hand.
        per_sheet = [
            sum(a.GetFormalCharge() for a in sheet.mol.GetAtoms())
            for sheet in sheets
        ]
        total_charge = sum(per_sheet)

        with open(filepath, "w") as f:
            f.write("; Surface topology — GROMACS simulation\n")
            f.write("; Auto-generated by Biochar Surface Builder\n")
            f.write(
                f"; Total system charge: {total_charge:+d} e "
                f"({len(sheets)} sheets)\n"
            )
            if total_charge != 0:
                f.write(
                    "; NOTE: this system is intentionally charged. Neutralise "
                    "it with `gmx genion -neutral` at solvation.\n"
                )
            f.write("\n")
            f.write(f'#include "{forcefield_path}"\n\n')

            if sheets_identical:
                f.write(f'#include "{Path(itp_paths[0]).name}"\n\n')
                f.write("[ system ]\n")
                f.write(f"{system_name}\n\n")
                f.write("[ molecules ]\n")
                f.write("; Name                #molecules\n")
                f.write(
                    f"{sheets[0].molecule_name:<20s} {len(sheets)}\n"
                )
            else:
                for itp_path in itp_paths:
                    f.write(f'#include "{Path(itp_path).name}"\n')
                f.write("\n[ system ]\n")
                f.write(f"{system_name}\n\n")
                f.write("[ molecules ]\n")
                f.write("; Name                #molecules\n")
                for sheet in sheets:
                    f.write(f"{sheet.molecule_name:<20s} 1\n")
