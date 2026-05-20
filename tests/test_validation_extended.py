"""
Extended tests for validation.py — CompositionValidator, ChemicalFeasibilityValidator,
StructureValidator, and edge cases in ValidationEngine.

The existing tests/test_generator.py::TestValidation provides only a smoke test of
ValidationEngine.validate_complete(). These tests cover the individual validators
and their error/warning paths in detail.
"""

import pytest
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem

from biochar.heteroatom_assignment import CompositionInfo
from biochar.validation import (
    CompositionValidator,
    ChemicalFeasibilityValidator,
    StructureValidator,
    ValidationEngine,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_composition(num_C=20, num_H=10, num_O=0, H_C=0.5, O_C=0.0, **kwargs) -> CompositionInfo:
    return CompositionInfo(
        num_carbons=num_C,
        num_hydrogens=num_H,
        num_oxygens=num_O,
        H_C_ratio=H_C,
        O_C_ratio=O_C,
        functional_groups=kwargs.get("functional_groups", {}),
    )


def _benzene_3d():
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0)
    coords = mol.GetConformer().GetPositions()
    return mol, coords


# ---------------------------------------------------------------------------
# CompositionValidator
# ---------------------------------------------------------------------------

class TestCompositionValidator:
    def test_valid_composition_passes(self):
        comp = _make_composition(num_C=20, num_H=10, H_C=0.5, O_C=0.0)
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.0)
        assert report.is_valid
        assert report.errors == []

    def test_h_c_ratio_too_high_produces_error(self):
        # Target 0.5, actual 0.9 → >10% error
        comp = _make_composition(num_C=20, num_H=18, H_C=0.9, O_C=0.0)
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.0)
        assert not report.is_valid
        assert any("H/C" in e for e in report.errors)

    def test_o_c_ratio_too_high_produces_error(self):
        # Target 0.1, actual 0.5 → >10% error
        comp = _make_composition(num_C=20, num_H=10, num_O=10, H_C=0.5, O_C=0.5)
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.1)
        assert not report.is_valid
        assert any("O/C" in e for e in report.errors)

    def test_zero_target_o_c_but_oxygens_present_is_warning(self):
        comp = _make_composition(num_C=20, num_H=10, num_O=2, H_C=0.5, O_C=0.1)
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.0)
        # Should produce a warning, not necessarily an error
        assert any("oxygen" in w.lower() or "O/C" in w for w in report.warnings)

    def test_too_few_carbons_produces_error(self):
        comp = _make_composition(num_C=2, num_H=2, H_C=1.0, O_C=0.0)
        report = CompositionValidator.validate(comp, target_H_C=1.0, target_O_C=0.0)
        assert any("carbon" in e.lower() for e in report.errors)

    def test_no_hydrogens_produces_error(self):
        comp = _make_composition(num_C=20, num_H=0, H_C=0.0, O_C=0.0)
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.0)
        assert any("hydrogen" in e.lower() for e in report.errors)

    def test_metrics_dict_contains_ratio_fields(self):
        comp = _make_composition(num_C=20, num_H=10, H_C=0.5, O_C=0.0)
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.0)
        assert "actual_H_C" in report.metrics
        assert "target_H_C" in report.metrics
        assert "H_C_error_pct" in report.metrics

    def test_report_is_valid_field_is_bool(self):
        comp = _make_composition()
        report = CompositionValidator.validate(comp, target_H_C=0.5, target_O_C=0.0)
        assert isinstance(report.is_valid, bool)


# ---------------------------------------------------------------------------
# ChemicalFeasibilityValidator
# ---------------------------------------------------------------------------

class TestChemicalFeasibilityValidator:
    def test_valid_benzene_passes(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        report = ChemicalFeasibilityValidator.validate(mol)
        assert report.is_valid
        assert report.errors == []

    def test_none_mol_fails(self):
        report = ChemicalFeasibilityValidator.validate(None)
        assert not report.is_valid
        assert any("None" in e or "none" in e.lower() for e in report.errors)

    def test_empty_mol_fails(self):
        mol = Chem.RWMol().GetMol()
        report = ChemicalFeasibilityValidator.validate(mol)
        assert not report.is_valid

    def test_aromaticity_percent_in_metrics(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        report = ChemicalFeasibilityValidator.validate(mol)
        assert "aromaticity_percent" in report.metrics

    def test_aromatic_fraction_benzene(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        report = ChemicalFeasibilityValidator.validate(mol)
        # All 6 atoms are aromatic
        assert report.metrics["aromaticity_percent"] == pytest.approx(100.0)

    def test_ethanol_is_feasible(self):
        mol = Chem.MolFromSmiles("CCO")
        mol = Chem.AddHs(mol)
        report = ChemicalFeasibilityValidator.validate(mol)
        assert report.is_valid


# ---------------------------------------------------------------------------
# StructureValidator
# ---------------------------------------------------------------------------

class TestStructureValidator:
    def test_valid_benzene_no_coords(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        report = StructureValidator.validate(mol)
        assert report.is_valid

    def test_valid_benzene_with_coords(self):
        mol, coords = _benzene_3d()
        report = StructureValidator.validate(mol, coords)
        assert report.is_valid

    def test_none_mol_fails(self):
        report = StructureValidator.validate(None)
        assert not report.is_valid

    def test_no_carbon_fails(self):
        mol = Chem.MolFromSmiles("O")
        mol = Chem.AddHs(mol)
        report = StructureValidator.validate(mol)
        assert not report.is_valid
        assert any("carbon" in e.lower() for e in report.errors)

    def test_no_hydrogen_fails(self):
        # Pure carbon with no H — just a carbon radical, not realistic
        mol = Chem.RWMol()
        mol.AddAtom(Chem.Atom(6))
        mol = mol.GetMol()
        report = StructureValidator.validate(mol)
        assert not report.is_valid
        assert any("hydrogen" in e.lower() for e in report.errors)

    def test_disconnected_molecule_fails(self):
        # Two separate molecules in one RDKit Mol → disconnected graph
        mol = Chem.MolFromSmiles("C.C")
        mol = Chem.AddHs(mol)
        report = StructureValidator.validate(mol)
        assert not report.is_valid
        assert any("connect" in e.lower() for e in report.errors)

    def test_metrics_contains_num_atoms(self):
        mol, coords = _benzene_3d()
        report = StructureValidator.validate(mol, coords)
        assert "num_atoms" in report.metrics
        assert report.metrics["num_atoms"] == mol.GetNumAtoms()

    def test_is_connected_single_atom(self):
        """Single atom should be considered connected."""
        assert StructureValidator._is_connected(Chem.MolFromSmiles("[CH4]"))

    def test_is_connected_ethane(self):
        mol = Chem.MolFromSmiles("CC")
        assert StructureValidator._is_connected(mol)

    def test_is_not_connected_two_fragments(self):
        mol = Chem.MolFromSmiles("C.C")
        assert not StructureValidator._is_connected(mol)


# ---------------------------------------------------------------------------
# ValidationEngine (edge cases beyond the smoke test in test_generator.py)
# ---------------------------------------------------------------------------

class TestValidationEngine:
    def test_returns_four_tuple(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        comp = _make_composition(num_C=6, num_H=6, H_C=1.0, O_C=0.0)
        result = ValidationEngine.validate_complete(mol, comp)
        assert len(result) == 4

    def test_all_valid_benzene(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        comp = _make_composition(num_C=6, num_H=6, H_C=1.0, O_C=0.0)
        is_valid, errors, warnings, metrics = ValidationEngine.validate_complete(
            mol, comp, target_H_C=1.0, target_O_C=0.0
        )
        assert isinstance(is_valid, bool)
        assert isinstance(errors, list)
        assert isinstance(warnings, list)
        assert isinstance(metrics, dict)

    def test_invalid_composition_propagates_error(self):
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        # H/C=5.0 vs target 0.5 → composition error
        comp = _make_composition(num_C=6, num_H=30, H_C=5.0, O_C=0.0)
        is_valid, errors, _, _ = ValidationEngine.validate_complete(
            mol, comp, target_H_C=0.5, target_O_C=0.0
        )
        assert not is_valid
        assert len(errors) > 0

    def test_none_mol_propagates_error(self):
        comp = _make_composition()
        is_valid, errors, _, _ = ValidationEngine.validate_complete(None, comp)
        assert not is_valid


# ---------------------------------------------------------------------------
# Constants sanity
# ---------------------------------------------------------------------------

class TestConstantsSanity:
    def test_functional_groups_has_required_entries(self):
        from biochar.constants import FUNCTIONAL_GROUPS
        for name in ("phenolic", "carboxyl", "ether", "amino"):
            assert name in FUNCTIONAL_GROUPS, f"Missing functional group: {name}"

    def test_each_functional_group_has_o_per_group(self):
        from biochar.constants import FUNCTIONAL_GROUPS
        for name, data in FUNCTIONAL_GROUPS.items():
            assert "O_per_group" in data, f"{name} missing O_per_group"

    def test_gromacs_type_map_covers_core_types(self):
        from biochar.constants import GROMACS_OPLS_TYPE_MAP
        for t in ("CA", "HA", "CT", "HC", "OH", "OS"):
            assert t in GROMACS_OPLS_TYPE_MAP, f"Missing GROMACS mapping for {t}"

    def test_opls_lj_params_sigma_non_negative(self):
        from biochar.constants import OPLS_LJ_PARAMS
        for atom_type, (sigma, epsilon) in OPLS_LJ_PARAMS.items():
            assert sigma >= 0, f"{atom_type}: negative sigma {sigma}"
            assert epsilon >= 0, f"{atom_type}: negative epsilon {epsilon}"

    def test_opls_atom_types_entries_have_three_fields(self):
        from biochar.constants import OPLS_ATOM_TYPES
        for key, val in OPLS_ATOM_TYPES.items():
            assert len(val) == 3, f"OPLS_ATOM_TYPES['{key}'] should have 3 fields"


# ---------------------------------------------------------------------------
# OPLSTyping — specific type assignments
# ---------------------------------------------------------------------------

class TestOPLSTypingSpecific:
    def test_aromatic_carbons_typed_as_CA(self):
        from biochar.opls_typing import AtomTyper
        mol = Chem.MolFromSmiles("c1ccccc1")
        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and atom.GetAtomicNum() == 6:
                assert atom_types[atom.GetIdx()] == "CA", (
                    f"Aromatic C at {atom.GetIdx()} should be CA, got {atom_types[atom.GetIdx()]}"
                )

    def test_aromatic_hydrogens_typed_as_HA(self):
        from biochar.opls_typing import AtomTyper
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                neighbor = list(atom.GetNeighbors())[0]
                if neighbor.GetIsAromatic():
                    assert atom_types[atom.GetIdx()] == "HA", (
                        f"H on aromatic C should be HA, got {atom_types[atom.GetIdx()]}"
                    )

    def test_property_table_length_matches_atom_count(self):
        from biochar.opls_typing import AtomTyper, ChargeAssigner, OPLSPropertyTable
        mol = Chem.MolFromSmiles("c1ccccc1O")
        mol = Chem.AddHs(mol)
        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)
        charger = ChargeAssigner()
        charges = charger.assign_charges(mol, atom_types)
        table = OPLSPropertyTable(mol, atom_types, charges)
        props = table.get_properties()
        assert len(props) == mol.GetNumAtoms()

    def test_charge_sum_is_approximately_zero(self):
        """Total charge of a neutral molecule should be near zero."""
        from biochar.opls_typing import AtomTyper, ChargeAssigner, OPLSPropertyTable
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)
        charger = ChargeAssigner()
        charges = charger.assign_charges(mol, atom_types)
        table = OPLSPropertyTable(mol, atom_types, charges)
        total_q = table.get_total_charge()
        assert abs(total_q) < 0.3  # Approximate neutrality after OPLS equilibration


# ---------------------------------------------------------------------------
# GeometryValidator — uncovered public methods
# ---------------------------------------------------------------------------

class TestGeometryValidatorPublicMethods:
    def test_check_aromaticity_benzene(self):
        from biochar.geometry_3d import GeometryValidator
        mol = Chem.MolFromSmiles("c1ccccc1")
        is_aromatic, num_aromatic = GeometryValidator.check_aromaticity(mol)
        assert is_aromatic
        assert num_aromatic == 6

    def test_check_aromaticity_ethane_not_aromatic(self):
        from biochar.geometry_3d import GeometryValidator
        mol = Chem.MolFromSmiles("CC")
        is_aromatic, num_aromatic = GeometryValidator.check_aromaticity(mol)
        assert not is_aromatic
        assert num_aromatic == 0

    def test_measure_ring_planarity_benzene(self):
        from biochar.geometry_3d import GeometryValidator
        mol, coords = _benzene_3d()
        planarity, assessment = GeometryValidator.measure_ring_planarity(mol, coords)
        assert isinstance(planarity, float)
        assert isinstance(assessment, str)
        # Benzene is a planar ring; deviation should be very small
        assert planarity < 0.5

    def test_measure_ring_planarity_no_rings(self):
        from biochar.geometry_3d import GeometryValidator
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=0)
        coords = mol.GetConformer().GetPositions()
        planarity, _ = GeometryValidator.measure_ring_planarity(mol, coords)
        # No rings: should return 0 or some default
        assert planarity >= 0
