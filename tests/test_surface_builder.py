"""
Tests for src/surface_builder.py — porous slit-pore surface generation.
"""

import pytest
from pathlib import Path

import numpy as np

from biochar_simulator.surface_builder import SurfaceBuilder, SurfaceConfig, SheetResult
from biochar_simulator.biochar_generator import generate_surface


# ---------------------------------------------------------------------------
# SurfaceConfig validation
# ---------------------------------------------------------------------------

class TestSurfaceConfigValidation:
    def test_default_config_valid(self):
        """Default config should construct without errors."""
        SurfaceConfig()

    def test_invalid_pore_type(self):
        with pytest.raises(ValueError, match="pore_type"):
            SurfaceConfig(pore_type="amorphous")

    def test_invalid_num_sheets_one(self):
        with pytest.raises(ValueError, match="num_sheets"):
            SurfaceConfig(num_sheets=1)

    def test_invalid_num_sheets_zero(self):
        with pytest.raises(ValueError, match="num_sheets"):
            SurfaceConfig(num_sheets=0)

    def test_invalid_pore_diameter_negative(self):
        with pytest.raises(ValueError, match="pore_diameter"):
            SurfaceConfig(pore_diameter=-1.0)

    def test_invalid_pore_diameter_zero(self):
        with pytest.raises(ValueError, match="pore_diameter"):
            SurfaceConfig(pore_diameter=0.0)

    def test_sheet_base_name_too_long(self):
        with pytest.raises(ValueError, match="sheet_base_name"):
            SurfaceConfig(sheet_base_name="SHEET")  # 5 chars > 3

    def test_sheet_base_name_max_allowed(self):
        """Exactly 3 chars should be fine."""
        SurfaceConfig(sheet_base_name="SHT")

    def test_sheet_overrides_length_mismatch(self):
        with pytest.raises(ValueError, match="sheet_overrides"):
            SurfaceConfig(num_sheets=2, sheet_overrides=[{}])  # only 1 override for 2 sheets

    def test_sheet_overrides_correct_length(self):
        """Overrides list matching num_sheets should be valid."""
        SurfaceConfig(num_sheets=2, sheet_overrides=[{}, {}])


# ---------------------------------------------------------------------------
# Sheet generation
# ---------------------------------------------------------------------------

class TestSheetGeneration:
    def test_single_sheet_flatness(self):
        """After flatten_to_xy, z range of coords should be < 2 Å."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        sheet = builder._generate_single_sheet(0)

        z_range = sheet.coords[:, 2].max() - sheet.coords[:, 2].min()
        assert z_range < 2.0, f"Sheet not flat: z range = {z_range:.2f} Å"

    def test_sheet_has_atoms(self):
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        sheet = builder._generate_single_sheet(0)
        assert sheet.mol.GetNumAtoms() > 0

    def test_sheet_has_opls_types(self):
        """atom_types dict should be populated."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        sheet = builder._generate_single_sheet(0)
        assert len(sheet.atom_types) == sheet.mol.GetNumAtoms()

    def test_sheet_with_functional_groups(self):
        """Sheet with phenolic groups should have oxygens."""
        config = SurfaceConfig(
            target_num_carbons=24,
            functional_groups={"phenolic": 2},
            num_sheets=2,
            seed=42,
        )
        builder = SurfaceBuilder(config)
        sheet = builder._generate_single_sheet(0)
        num_O = sum(1 for a in sheet.mol.GetAtoms() if a.GetAtomicNum() == 8)
        assert num_O >= 2, f"Expected >= 2 oxygens, got {num_O}"

    def test_identical_sheets_copy(self):
        """Second sheet (identical) should be a deep copy, not the same object."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.sheets = [builder._generate_single_sheet(0)]
        sheet1 = builder._generate_single_sheet(1)

        # Modifying copy's coords should not affect original
        original_z = builder.sheets[0].coords[0, 2]
        sheet1.coords[0, 2] += 999.0
        assert builder.sheets[0].coords[0, 2] == original_z


# ---------------------------------------------------------------------------
# Sheet positioning
# ---------------------------------------------------------------------------

class TestSheetPositioning:
    def test_two_sheet_z_gap(self):
        """Centroid separation for 2 sheets should be pore_diameter + 3.4 Å."""
        config = SurfaceConfig(
            target_num_carbons=16,
            pore_diameter=10.0,
            num_sheets=2,
            seed=42,
        )
        builder = SurfaceBuilder(config)
        sheets, _ = builder.build()

        z0 = sheets[0].coords[:, 2].mean()
        z1 = sheets[1].coords[:, 2].mean()
        expected = 10.0 + 3.4
        assert abs((z1 - z0) - expected) < 0.5, (
            f"Sheet z separation {z1 - z0:.2f} Å ≠ expected {expected:.2f} Å"
        )

    def test_three_sheet_equal_gaps(self):
        """Three sheets should have equal inter-centroid gaps."""
        config = SurfaceConfig(
            target_num_carbons=16,
            pore_diameter=8.0,
            num_sheets=3,
            seed=42,
        )
        builder = SurfaceBuilder(config)
        sheets, _ = builder.build()

        z_means = [s.coords[:, 2].mean() for s in sheets]
        gap_01 = z_means[1] - z_means[0]
        gap_12 = z_means[2] - z_means[1]
        assert abs(gap_01 - gap_12) < 0.5, (
            f"Unequal gaps: {gap_01:.2f} vs {gap_12:.2f} Å"
        )

    def test_box_vectors_positive(self):
        """All three box dimensions should be positive."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        _, box = builder.build()
        assert all(v > 0 for v in box), f"Non-positive box vector: {box}"

    def test_coords_inside_box(self):
        """After centring, all atom coords should lie within the box (in nm)."""
        config = SurfaceConfig(
            target_num_carbons=16,
            num_sheets=2,
            box_padding_xy=1.0,
            box_padding_z=1.0,
            seed=42,
        )
        builder = SurfaceBuilder(config)
        sheets, box_nm = builder.build()

        all_coords_nm = np.vstack([s.coords for s in sheets]) * 0.1
        assert all_coords_nm.min() >= -0.01, "Atoms outside box (negative side)"
        assert (all_coords_nm[:, 0].max() <= box_nm[0] + 0.01), "Atoms outside box (x)"
        assert (all_coords_nm[:, 1].max() <= box_nm[1] + 0.01), "Atoms outside box (y)"
        assert (all_coords_nm[:, 2].max() <= box_nm[2] + 0.01), "Atoms outside box (z)"


# ---------------------------------------------------------------------------
# GROMACS file export — identical sheets
# ---------------------------------------------------------------------------

class TestGROMACSExportIdentical:
    def test_gro_exists(self, tmp_path):
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        gro, top, itps = builder.export_gromacs(str(tmp_path), "test")
        assert gro.exists()

    def test_gro_atom_count(self, tmp_path):
        """Second line of .gro must equal total atom count."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        gro, _, _ = builder.export_gromacs(str(tmp_path), "test")

        lines = gro.read_text().splitlines()
        total_atoms = int(lines[1].strip())
        expected = sum(s.mol.GetNumAtoms() for s in builder.sheets)
        assert total_atoms == expected

    def test_gro_residue_numbering(self, tmp_path):
        """First atom should be residue 1, last atom should be residue 2."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        gro, _, _ = builder.export_gromacs(str(tmp_path), "test")

        lines = gro.read_text().splitlines()
        total_atoms = int(lines[1].strip())
        atom_lines = lines[2:2 + total_atoms]

        first_resnum = int(atom_lines[0][:5].strip())
        last_resnum = int(atom_lines[-1][:5].strip())
        assert first_resnum == 1
        assert last_resnum == 2

    def test_gro_box_vectors_positive(self, tmp_path):
        """Last line of .gro must contain three positive floats."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        gro, _, _ = builder.export_gromacs(str(tmp_path), "test")

        last_line = gro.read_text().splitlines()[-1]
        vals = [float(x) for x in last_line.split()[:3]]
        assert all(v > 0 for v in vals), f"Box vectors not positive: {vals}"

    def test_top_exists(self, tmp_path):
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        _, top, _ = builder.export_gromacs(str(tmp_path), "test")
        assert top.exists()

    def test_single_itp_for_identical_sheets(self, tmp_path):
        """Identical sheets should produce exactly one .itp file."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=2, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        _, _, itps = builder.export_gromacs(str(tmp_path), "test")
        assert len(itps) == 1
        assert itps[0].exists()

    def test_top_molecule_count(self, tmp_path):
        """Identical sheets: .top [ molecules ] should list count = num_sheets."""
        config = SurfaceConfig(target_num_carbons=16, num_sheets=3, seed=42)
        builder = SurfaceBuilder(config)
        builder.build()
        _, top, _ = builder.export_gromacs(str(tmp_path), "test")

        top_text = top.read_text()
        # Should have " 3" in the molecules section
        assert " 3" in top_text, f"Expected count 3 in .top:\n{top_text}"


# ---------------------------------------------------------------------------
# GROMACS file export — distinct sheets
# ---------------------------------------------------------------------------

class TestGROMACSExportDistinct:
    def test_one_itp_per_sheet(self, tmp_path):
        """Distinct sheets should produce one .itp per sheet."""
        config = SurfaceConfig(
            target_num_carbons=16,
            num_sheets=2,
            sheet_overrides=[
                {"functional_groups": {"phenolic": 1}},
                {"functional_groups": {"ether": 1}},
            ],
            seed=42,
        )
        builder = SurfaceBuilder(config)
        builder.build()
        _, _, itps = builder.export_gromacs(str(tmp_path), "test")
        assert len(itps) == 2
        for itp in itps:
            assert itp.exists()

    def test_top_lists_each_sheet(self, tmp_path):
        """Distinct sheets: .top should have one entry per sheet type."""
        config = SurfaceConfig(
            target_num_carbons=16,
            num_sheets=2,
            sheet_overrides=[
                {"functional_groups": {"phenolic": 1}},
                {"functional_groups": {"ether": 1}},
            ],
            seed=42,
        )
        builder = SurfaceBuilder(config)
        builder.build()
        _, top, _ = builder.export_gromacs(str(tmp_path), "test")

        top_text = top.read_text()
        assert "SHT1" in top_text
        assert "SHT2" in top_text

    def test_distinct_sheet_names_in_gro(self, tmp_path):
        """Distinct sheets: residue names should be SHT1 and SHT2."""
        config = SurfaceConfig(
            target_num_carbons=16,
            num_sheets=2,
            sheet_overrides=[
                {"functional_groups": {"phenolic": 1}},
                {"functional_groups": {"phenolic": 1}},
            ],
            seed=42,
        )
        builder = SurfaceBuilder(config)
        builder.build()
        gro, _, _ = builder.export_gromacs(str(tmp_path), "test")

        gro_text = gro.read_text()
        assert "SHT1" in gro_text
        assert "SHT2" in gro_text


# ---------------------------------------------------------------------------
# Convenience function
# ---------------------------------------------------------------------------

class TestGenerateSurface:
    def test_basic_call(self, tmp_path):
        """generate_surface() should return files that all exist."""
        sheets, gro, top, itps = generate_surface(
            target_num_carbons=16,
            pore_diameter=10.0,
            output_directory=str(tmp_path),
            basename="surf",
            seed=42,
        )
        assert len(sheets) == 2
        assert gro.exists()
        assert top.exists()
        assert all(p.exists() for p in itps)

    def test_returns_sheet_results(self, tmp_path):
        """Each element of sheets list should be a SheetResult."""
        sheets, *_ = generate_surface(
            target_num_carbons=16,
            output_directory=str(tmp_path),
            seed=42,
        )
        for s in sheets:
            assert isinstance(s, SheetResult)

    def test_three_sheets(self, tmp_path):
        """Can generate more than 2 sheets."""
        sheets, gro, top, itps = generate_surface(
            target_num_carbons=16,
            num_sheets=3,
            pore_diameter=8.0,
            output_directory=str(tmp_path),
            seed=42,
        )
        assert len(sheets) == 3

        lines = gro.read_text().splitlines()
        total_atoms = int(lines[1].strip())
        expected = sum(s.mol.GetNumAtoms() for s in sheets)
        assert total_atoms == expected

    def test_asymmetric_pore(self, tmp_path):
        """Asymmetric pore (distinct overrides) should produce 2 .itp files."""
        sheets, _, _, itps = generate_surface(
            pore_diameter=8.0,
            sheet_overrides=[
                {"functional_groups": {"phenolic": 2}, "target_num_carbons": 16},
                {"functional_groups": {"carboxyl": 1}, "target_num_carbons": 16},
            ],
            output_directory=str(tmp_path),
            seed=42,
        )
        assert len(itps) == 2
