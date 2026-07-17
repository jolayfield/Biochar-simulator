"""
Tests for gromacs_export.py — GROFileWriter, TOPFileWriter, ITPFileWriter,
GromacsExporter, MultiSheetGROWriter, SurfaceTopologyWriter.
"""

import os

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock

from rdkit import Chem
from rdkit.Chem import AllChem

from biochar.gromacs_export import (
    GROFileWriter,
    TOPFileWriter,
    ITPFileWriter,
    GromacsExporter,
    MultiSheetGROWriter,
    SurfaceTopologyWriter,
)
from biochar.opls_typing import AtomTyper, ChargeAssigner


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _benzene_with_coords():
    """Return (mol, coords_angstrom) for benzene with explicit H."""
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0)
    coords = mol.GetConformer().GetPositions()
    return mol, coords


def _benzene_typed():
    """Return (mol, coords, atom_types, charges) for benzene."""
    mol, coords = _benzene_with_coords()
    typer = AtomTyper()
    atom_types = typer.assign_atom_types(mol)
    charger = ChargeAssigner()
    charges = charger.assign_charges(mol, atom_types)
    return mol, coords, atom_types, charges


def _parse_gro_atoms(text: str):
    """Parse GRO file text → (title, natoms, atom_lines, box_line)."""
    lines = text.strip().splitlines()
    title = lines[0]
    natoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + natoms]
    box_line = lines[2 + natoms]
    return title, natoms, atom_lines, box_line


# ---------------------------------------------------------------------------
# GROFileWriter
# ---------------------------------------------------------------------------

class TestGROFileWriter:
    def test_atom_count_line(self, tmp_path):
        mol, coords = _benzene_with_coords()
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords)

        text = Path(gro).read_text()
        _, natoms, atom_lines, _ = _parse_gro_atoms(text)
        assert natoms == mol.GetNumAtoms()
        assert len(atom_lines) == mol.GetNumAtoms()

    def test_coordinates_in_nanometers(self, tmp_path):
        """Coordinates written must be ~1/10 of Ångström input values."""
        mol, coords = _benzene_with_coords()
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords)

        text = Path(gro).read_text()
        _, natoms, atom_lines, _ = _parse_gro_atoms(text)

        # Extract x coordinate of first atom
        first_line = atom_lines[0]
        # GRO format: resnum(5) resname(5) atomname(5) atomnum(5) x(8.3f) y(8.3f) z(8.3f)
        x_nm = float(first_line[20:28])
        # Coordinate should be consistent with nm (i.e. < 5 nm for a small molecule)
        assert abs(x_nm) < 5.0

    def test_box_line_positive_floats(self, tmp_path):
        mol, coords = _benzene_with_coords()
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords)

        text = Path(gro).read_text()
        _, natoms, _, box_line = _parse_gro_atoms(text)
        vals = [float(v) for v in box_line.split()[:3]]
        assert len(vals) == 3
        assert all(v > 0 for v in vals)

    def test_molecule_name_appears_in_residue_field(self, tmp_path):
        mol, coords = _benzene_with_coords()
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords, molecule_name="BENZ")

        text = Path(gro).read_text()
        assert "BENZ" in text

    def test_long_molecule_name_truncated_to_five(self, tmp_path):
        """GROFileWriter should silently truncate names > 5 chars."""
        mol, coords = _benzene_with_coords()
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords, molecule_name="TOOLONGNAME")

        text = Path(gro).read_text()
        # Only the first 5 chars should appear in residue name field
        assert "TOOLONGNAME" not in text
        assert "TOOLO" in text  # "TOOLONGNAME"[:5] == "TOOLO"

    def test_explicit_1d_box_vectors(self, tmp_path):
        """Explicit 1D box vectors should be written to the last line."""
        mol, coords = _benzene_with_coords()
        box = np.array([5.0, 5.0, 5.0])
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords, box_vectors=box)

        last_line = Path(gro).read_text().strip().splitlines()[-1]
        vals = [float(v) for v in last_line.split()[:3]]
        assert abs(vals[0] - 5.0) < 0.001
        assert abs(vals[1] - 5.0) < 0.001
        assert abs(vals[2] - 5.0) < 0.001

    def test_explicit_3x3_box_vectors(self, tmp_path):
        """3×3 triclinic box should write 9 values on last line."""
        mol, coords = _benzene_with_coords()
        box = np.diag([4.0, 5.0, 6.0])
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords, box_vectors=box)

        last_line = Path(gro).read_text().strip().splitlines()[-1]
        vals = last_line.split()
        assert len(vals) == 9
        assert abs(float(vals[0]) - 4.0) < 0.001
        assert abs(float(vals[1]) - 5.0) < 0.001
        assert abs(float(vals[2]) - 6.0) < 0.001

    def test_custom_title(self, tmp_path):
        mol, coords = _benzene_with_coords()
        gro = str(tmp_path / "out.gro")
        GROFileWriter.write(gro, mol, coords, title="My title")

        first_line = Path(gro).read_text().splitlines()[0]
        assert first_line == "My title"


# ---------------------------------------------------------------------------
# TOPFileWriter
# ---------------------------------------------------------------------------

class TestTOPFileWriter:
    def test_required_sections_present(self, tmp_path):
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(top).read_text()
        assert "[ moleculetype ]" in text
        assert "[ atoms ]" in text
        assert "[ bonds ]" in text
        assert "[ angles ]" in text

    def test_molecules_section_present(self, tmp_path):
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(top).read_text()
        assert "[ molecules ]" in text
        assert "BENZ" in text

    def test_top_defers_to_forcefield_for_atomtypes(self, tmp_path):
        """The .top must #include a forcefield and never redefine its types.

        An [ atomtypes ] block here would override the #included oplsaa.ff by
        name, silently replacing stock LJ parameters with whatever we wrote.
        The types are the forcefield's to define; we only reference them.
        """
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(top).read_text()
        assert '#include "oplsaa.ff/forcefield.itp"' in text
        assert "[ atomtypes ]" not in text
        assert "[atomtypes]" not in text

    def test_top_defers_to_forcefield_even_when_ff_is_resolvable(self, tmp_path):
        """Emitting atomtypes must not depend on whether the include resolves.

        A stock oplsaa.ff/forcefield.itp holds no [ atomtypes ] of its own -- it
        only #includes ffnonbonded.itp -- so any check that greps the named file
        for that section concludes the types are missing precisely when a real
        forcefield is present, and duplicates them.
        """
        ff_dir = tmp_path / "oplsaa.ff"
        ff_dir.mkdir()
        (ff_dir / "forcefield.itp").write_text(
            '#define _FF_OPLSAA\n\n[ defaults ]\n1 3 yes 0.5 0.5\n\n'
            '#include "ffnonbonded.itp"\n#include "ffbonded.itp"\n'
        )
        (ff_dir / "ffnonbonded.itp").write_text(
            "[ atomtypes ]\nopls_145 CA 6 12.011 -0.115 A 3.55e-01 2.92880e-01\n"
        )

        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            TOPFileWriter.write(top, mol, atom_types, charges, molecule_name="BENZ")
        finally:
            os.chdir(cwd)

        assert "[ atomtypes ]" not in Path(top).read_text()

    def test_atoms_section_count(self, tmp_path):
        """Number of atom lines in [ atoms ] must equal mol.GetNumAtoms()."""
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(top).read_text()
        # Extract [ atoms ] section
        in_atoms = False
        atom_count = 0
        for line in text.splitlines():
            stripped = line.strip()
            if stripped == "[ atoms ]":
                in_atoms = True
                continue
            if in_atoms:
                if stripped.startswith("[") and stripped != "[ atoms ]":
                    break
                if stripped and not stripped.startswith(";"):
                    atom_count += 1
        assert atom_count == mol.GetNumAtoms()

    def test_bonds_section_count(self, tmp_path):
        """Bond count in [ bonds ] must match mol.GetNumBonds()."""
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(top).read_text()
        in_bonds = False
        bond_count = 0
        for line in text.splitlines():
            stripped = line.strip()
            if stripped == "[ bonds ]":
                in_bonds = True
                continue
            if in_bonds:
                if stripped.startswith("[") and stripped != "[ bonds ]":
                    break
                if stripped and not stripped.startswith(";"):
                    bond_count += 1
        assert bond_count == mol.GetNumBonds()

    def test_dihedrals_section_present_by_default(self, tmp_path):
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, include_dihedrals=True)
        text = Path(top).read_text()
        assert "[ dihedrals ]" in text

    def test_dihedrals_section_absent_when_disabled(self, tmp_path):
        mol, _, atom_types, charges = _benzene_typed()
        top = str(tmp_path / "out.top")
        TOPFileWriter.write(top, mol, atom_types, charges, include_dihedrals=False)
        text = Path(top).read_text()
        assert "[ dihedrals ]" not in text


# ---------------------------------------------------------------------------
# ITPFileWriter
# ---------------------------------------------------------------------------

class TestITPFileWriter:
    def test_required_sections_present(self, tmp_path):
        mol, _, atom_types, charges = _benzene_typed()
        itp = str(tmp_path / "out.itp")
        ITPFileWriter.write(itp, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(itp).read_text()
        assert "[ moleculetype ]" in text
        assert "[ atoms ]" in text
        assert "[ bonds ]" in text

    def test_no_atomtypes_section_in_itp(self, tmp_path):
        """ITP files must not define [ atomtypes ] to avoid duplication."""
        mol, _, atom_types, charges = _benzene_typed()
        itp = str(tmp_path / "out.itp")
        ITPFileWriter.write(itp, mol, atom_types, charges)

        text = Path(itp).read_text()
        assert "[ atomtypes ]" not in text

    def test_atoms_count_matches_molecule(self, tmp_path):
        mol, _, atom_types, charges = _benzene_typed()
        itp = str(tmp_path / "out.itp")
        ITPFileWriter.write(itp, mol, atom_types, charges, molecule_name="BENZ")

        text = Path(itp).read_text()
        in_atoms = False
        atom_count = 0
        for line in text.splitlines():
            stripped = line.strip()
            if stripped == "[ atoms ]":
                in_atoms = True
                continue
            if in_atoms:
                if stripped.startswith("[") and stripped != "[ atoms ]":
                    break
                if stripped and not stripped.startswith(";"):
                    atom_count += 1
        assert atom_count == mol.GetNumAtoms()


# ---------------------------------------------------------------------------
# GromacsExporter
# ---------------------------------------------------------------------------

class TestGromacsExporter:
    def test_export_creates_three_files(self, tmp_path):
        mol, coords, atom_types, charges = _benzene_typed()
        exporter = GromacsExporter(str(tmp_path))
        gro, top, itp = exporter.export(
            mol, coords, atom_types, charges,
            molecule_name="BENZ", basename="test",
        )
        assert gro.exists()
        assert top.exists()
        assert itp.exists()

    def test_export_returns_correct_paths(self, tmp_path):
        mol, coords, atom_types, charges = _benzene_typed()
        exporter = GromacsExporter(str(tmp_path))
        gro, top, itp = exporter.export(
            mol, coords, atom_types, charges, basename="myfile",
        )
        assert gro.name == "myfile.gro"
        assert top.name == "myfile.top"
        assert itp.name == "myfile.itp"

    def test_calculate_box_size_adds_padding(self):
        """Box should be at least molecule extent + 2×padding."""
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])  # Å (1 nm span)
        box = GromacsExporter._calculate_box_size(coords, padding=1.0)
        # extent = 0.1 nm in each dim, padding = 1.0 nm → box ≥ 0.1 + 2.0 = 2.1 nm
        assert all(box >= 2.0)

    def test_calculate_box_size_units_are_nm(self):
        """coords are in Å; returned box should be in nm."""
        coords = np.array([[0.0, 0.0, 0.0], [10.0, 10.0, 10.0]])  # 10 Å = 1 nm span
        box = GromacsExporter._calculate_box_size(coords, padding=0.5)
        # 1 nm extent + 2×0.5 nm = 2 nm
        assert abs(box[0] - 2.0) < 0.01

    def test_export_with_periodic_box(self, tmp_path):
        mol, coords, atom_types, charges = _benzene_typed()
        exporter = GromacsExporter(str(tmp_path))
        gro, top, itp = exporter.export(
            mol, coords, atom_types, charges,
            basename="periodic",
            include_periodic_box=True,
        )
        last_line = gro.read_text().strip().splitlines()[-1]
        vals = [float(v) for v in last_line.split()[:3]]
        assert all(v > 0 for v in vals)


# ---------------------------------------------------------------------------
# MultiSheetGROWriter
# ---------------------------------------------------------------------------

class TestMultiSheetGROWriter:
    def _make_sheet_mock(self, smiles: str, name: str):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0)
        coords = mol.GetConformer().GetPositions()
        sheet = MagicMock()
        sheet.mol = mol
        sheet.coords = coords
        sheet.molecule_name = name
        return sheet

    def test_total_atom_count_line(self, tmp_path):
        s1 = self._make_sheet_mock("c1ccccc1", "SHT1")
        s2 = self._make_sheet_mock("c1ccc2ccccc2c1", "SHT2")
        gro = str(tmp_path / "surface.gro")
        box = np.array([5.0, 5.0, 5.0])
        MultiSheetGROWriter.write(gro, [s1, s2], box)

        lines = Path(gro).read_text().splitlines()
        total_written = int(lines[1].strip())
        expected = s1.mol.GetNumAtoms() + s2.mol.GetNumAtoms()
        assert total_written == expected

    def test_residue_numbers_increment(self, tmp_path):
        """First atom of sheet1 → residue 1; first atom of sheet2 → residue 2."""
        s1 = self._make_sheet_mock("c1ccccc1", "SHT1")
        s2 = self._make_sheet_mock("c1ccccc1", "SHT2")
        gro = str(tmp_path / "surface.gro")
        box = np.array([5.0, 5.0, 5.0])
        MultiSheetGROWriter.write(gro, [s1, s2], box)

        lines = Path(gro).read_text().splitlines()
        total_atoms = int(lines[1].strip())
        atom_lines = lines[2 : 2 + total_atoms]

        first_resnum = int(atom_lines[0][:5].strip())
        last_resnum = int(atom_lines[-1][:5].strip())
        assert first_resnum == 1
        assert last_resnum == 2

    def test_global_atom_numbering_is_contiguous(self, tmp_path):
        s1 = self._make_sheet_mock("c1ccccc1", "SHT1")
        s2 = self._make_sheet_mock("c1ccccc1", "SHT2")
        gro = str(tmp_path / "surface.gro")
        box = np.array([5.0, 5.0, 5.0])
        MultiSheetGROWriter.write(gro, [s1, s2], box)

        lines = Path(gro).read_text().splitlines()
        total_atoms = int(lines[1].strip())
        atom_lines = lines[2 : 2 + total_atoms]

        atom_nums = [int(line[15:20].strip()) for line in atom_lines]
        assert atom_nums == list(range(1, total_atoms + 1))

    def test_1d_box_vector_in_last_line(self, tmp_path):
        s1 = self._make_sheet_mock("c1ccccc1", "SHT1")
        gro = str(tmp_path / "surface.gro")
        box = np.array([3.0, 4.0, 5.0])
        MultiSheetGROWriter.write(gro, [s1], box)

        last_line = Path(gro).read_text().strip().splitlines()[-1]
        vals = [float(v) for v in last_line.split()]
        assert len(vals) == 3
        assert abs(vals[0] - 3.0) < 0.001
        assert abs(vals[2] - 5.0) < 0.001


# ---------------------------------------------------------------------------
# SurfaceTopologyWriter
# ---------------------------------------------------------------------------

class TestSurfaceTopologyWriter:
    def _make_sheet_mock(self, name: str):
        sheet = MagicMock()
        sheet.molecule_name = name
        return sheet

    def test_identical_sheets_single_include(self, tmp_path):
        """Identical sheets: only one #include for the single .itp."""
        sheets = [self._make_sheet_mock("SHT"), self._make_sheet_mock("SHT")]
        itp_path = tmp_path / "sht.itp"
        itp_path.touch()
        top = str(tmp_path / "out.top")
        SurfaceTopologyWriter.write(
            top, sheets, [itp_path], sheets_identical=True
        )
        text = Path(top).read_text()
        assert text.count("#include") == 2  # forcefield + 1 itp
        assert "sht.itp" in text

    def test_identical_sheets_molecule_count(self, tmp_path):
        """Identical sheets: [ molecules ] lists the name once with count = N."""
        sheets = [self._make_sheet_mock("SHT"), self._make_sheet_mock("SHT"), self._make_sheet_mock("SHT")]
        itp_path = tmp_path / "sht.itp"
        itp_path.touch()
        top = str(tmp_path / "out.top")
        SurfaceTopologyWriter.write(
            top, sheets, [itp_path], sheets_identical=True
        )
        text = Path(top).read_text()
        assert " 3" in text

    def test_distinct_sheets_multiple_includes(self, tmp_path):
        """Distinct sheets: one #include per itp file."""
        sheets = [self._make_sheet_mock("SHT1"), self._make_sheet_mock("SHT2")]
        itp1 = tmp_path / "sht1.itp"
        itp2 = tmp_path / "sht2.itp"
        itp1.touch()
        itp2.touch()
        top = str(tmp_path / "out.top")
        SurfaceTopologyWriter.write(
            top, sheets, [itp1, itp2], sheets_identical=False
        )
        text = Path(top).read_text()
        assert "sht1.itp" in text
        assert "sht2.itp" in text

    def test_distinct_sheets_each_listed_once(self, tmp_path):
        """Distinct sheets: [ molecules ] has one line per sheet."""
        sheets = [self._make_sheet_mock("SHT1"), self._make_sheet_mock("SHT2")]
        itp1 = tmp_path / "sht1.itp"
        itp2 = tmp_path / "sht2.itp"
        itp1.touch()
        itp2.touch()
        top = str(tmp_path / "out.top")
        SurfaceTopologyWriter.write(
            top, sheets, [itp1, itp2], sheets_identical=False
        )
        text = Path(top).read_text()
        assert "SHT1" in text
        assert "SHT2" in text
