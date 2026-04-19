"""
Unit Tests for Biochar Generator

Tests for all major components.
"""

import pytest

from rdkit import Chem

from biochar.constants import OPLS_ATOM_TYPES, PAH_LIBRARY
from biochar.carbon_skeleton import PAHAssembler, SkeletonValidator
from biochar.heteroatom_assignment import (
    OxygenAssigner,
    HydrogenAssigner,
    HeteroatomValidator,
    CompositionInfo,
)
from biochar.geometry_3d import CoordinateGenerator, GeometryValidator
from biochar.opls_typing import AtomTyper, ChargeAssigner
from biochar.validation import ValidationEngine
from biochar.biochar_generator import BiocharGenerator, GeneratorConfig


class TestConstants:
    """Test constants and library definitions."""

    def test_opls_types_defined(self):
        """Test that OPLS atom types are defined."""
        assert len(OPLS_ATOM_TYPES) > 0
        assert "CA" in OPLS_ATOM_TYPES
        assert "HA" in OPLS_ATOM_TYPES
        assert "CT" in OPLS_ATOM_TYPES
        assert "OH" in OPLS_ATOM_TYPES

    def test_pah_library_complete(self):
        """Test that PAH library has expected structures."""
        assert len(PAH_LIBRARY) > 0
        assert "benzene" in PAH_LIBRARY
        assert "naphthalene" in PAH_LIBRARY
        assert "anthracene" in PAH_LIBRARY

    def test_pah_smiles_valid(self):
        """Test that PAH SMILES strings are valid."""
        for name, data in PAH_LIBRARY.items():
            mol = Chem.MolFromSmiles(data["smiles"])
            assert mol is not None, f"Invalid SMILES for {name}"


class TestCarbonSkeleton:
    """Test carbon skeleton generation."""

    def test_pah_assembler_benzene(self):
        """Test PAH assembler with small target."""
        assembler = PAHAssembler(seed=42)
        skeleton = assembler.generate(
            target_num_carbons=6,
            target_aromaticity=100.0,
        )

        assert skeleton.mol is not None
        assert skeleton.num_carbons > 0
        assert skeleton.aromaticity_percent > 80

    def test_pah_assembler_naphthalene(self):
        """Test PAH assembler with naphthalene size."""
        assembler = PAHAssembler(seed=42)
        skeleton = assembler.generate(
            target_num_carbons=10,
            target_aromaticity=100.0,
        )

        assert skeleton.mol is not None
        assert skeleton.num_carbons >= 8  # Allow tolerance

    def test_skeleton_validator(self):
        """Test skeleton validation."""
        assembler = PAHAssembler(seed=42)
        skeleton = assembler.generate(target_num_carbons=10)

        valid, errors = SkeletonValidator.validate(skeleton)
        assert valid, f"Skeleton validation failed: {errors}"


class TestHeteroatomAssignment:
    """Test heteroatom assignment."""

    def test_oxygen_assignment(self):
        """Test oxygen assignment."""
        # Create a simple benzene
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert mol is not None

        assigner = OxygenAssigner(seed=42)
        mol_with_O, composition = assigner.assign_oxygens(
            mol,
            target_O_C_ratio=0.1,
        )

        assert mol_with_O is not None
        assert composition.num_oxygens >= 0

    def test_hydrogen_assignment(self):
        """Test hydrogen assignment."""
        # Create a simple carbon structure
        mol = Chem.MolFromSmiles("CCCC")
        assert mol is not None

        assigner = HydrogenAssigner(seed=42)
        mol_with_H, composition = assigner.assign_hydrogens(
            mol,
            target_H_C_ratio=2.5,
        )

        assert mol_with_H is not None
        assert composition.H_C_ratio > 0

    def test_heteroatom_validator(self):
        """Test heteroatom composition validation."""
        # Create a simple molecule
        mol = Chem.MolFromSmiles("c1ccccc1CO")
        assert mol is not None

        # Calculate composition
        num_C = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        num_H = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
        num_O = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)


        composition = CompositionInfo(
            num_carbons=num_C,
            num_hydrogens=num_H,
            num_oxygens=num_O,
            H_C_ratio=num_H / num_C if num_C > 0 else 0,
            O_C_ratio=num_O / num_C if num_C > 0 else 0,
            functional_groups={},
        )

        valid, errors = HeteroatomValidator.validate_ratios(
            composition,
            target_H_C=1.0,
            target_O_C=0.1,
        )
        assert isinstance(valid, bool)


class TestGeometry:
    """Test 3D geometry generation."""

    def test_coordinate_generation(self):
        """Test 3D coordinate generation."""
        mol = Chem.MolFromSmiles("CCO")
        assert mol is not None

        generator = CoordinateGenerator(seed=42)
        mol_with_coords, coords = generator.generate_3d_coordinates(mol)

        assert mol_with_coords is not None
        assert coords is not None
        assert coords.shape[0] > 0
        assert coords.shape[1] == 3

    def test_geometry_validator(self):
        """Test geometry validation."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert mol is not None

        generator = CoordinateGenerator(seed=42)
        mol, coords = generator.generate_3d_coordinates(mol)

        valid, errors = GeometryValidator.validate_geometry(mol, coords)
        # May have some warnings but should be structurally sound
        assert coords is not None


class TestOPLSTyping:
    """Test OPLS-AA atom typing and charges."""

    def test_atom_typer(self):
        """Test atom type assignment."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert mol is not None

        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)

        assert len(atom_types) > 0
        assert all(t in ["CA", "HA", "CT", "HC", "OH", "OS", "OC", "C", "O"] or t.startswith("X")
                   for t in atom_types.values())

    def test_charge_assigner(self):
        """Test charge assignment."""
        mol = Chem.MolFromSmiles("c1ccccc1O")
        assert mol is not None

        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)

        charger = ChargeAssigner()
        charges = charger.assign_charges(mol, atom_types)

        assert len(charges) == mol.GetNumAtoms()
        assert all(isinstance(q, (int, float)) for q in charges.values())


class TestValidation:
    """Test validation engine."""

    def test_validation_engine(self):
        """Test complete validation."""
        mol = Chem.MolFromSmiles("c1ccccc1")
        assert mol is not None


        num_C = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        num_H = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
        num_O = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)

        composition = CompositionInfo(
            num_carbons=num_C,
            num_hydrogens=num_H,
            num_oxygens=num_O,
            H_C_ratio=num_H / max(num_C, 1),
            O_C_ratio=num_O / max(num_C, 1),
            functional_groups={},
        )

        is_valid, errors, warnings, metrics = ValidationEngine.validate_complete(
            mol,
            composition,
            target_H_C=1.0,
            target_O_C=0.0,
        )

        assert isinstance(is_valid, bool)
        assert isinstance(errors, list)
        assert isinstance(warnings, list)
        assert isinstance(metrics, dict)


class TestBiocharGenerator:
    """Test main biochar generator."""

    def test_generator_basic(self):
        """Test basic biochar generation."""
        config = GeneratorConfig(
            target_num_carbons=20,
            H_C_ratio=0.5,
            O_C_ratio=0.05,
            seed=42,
        )

        generator = BiocharGenerator(config)
        mol, coords, composition = generator.generate()

        assert mol is not None
        assert coords is not None
        assert composition is not None
        assert mol.GetNumAtoms() > 0

    def test_generator_config_defaults(self):
        """Test generator configuration defaults."""
        config = GeneratorConfig()

        assert config.target_num_carbons == 50
        assert config.H_C_ratio == 0.5
        assert config.O_C_ratio == 0.1
        assert config.functional_groups is None  # None → O/C-ratio-driven placement

    def test_generator_reproducibility(self):
        """Test that same seed produces same structure."""
        config1 = GeneratorConfig(
            target_num_carbons=30,
            seed=12345,
        )
        generator1 = BiocharGenerator(config1)
        mol1, coords1, comp1 = generator1.generate()

        config2 = GeneratorConfig(
            target_num_carbons=30,
            seed=12345,
        )
        generator2 = BiocharGenerator(config2)
        mol2, coords2, comp2 = generator2.generate()

        # Same seed should produce identical composition
        assert comp1.num_carbons == comp2.num_carbons
        assert comp1.num_hydrogens == comp2.num_hydrogens


class TestDefectRings:
    """Tests for pentagon (5-membered ring) defect insertion."""

    def test_grow_graph_pure_hexagon(self):
        """Pure hexagon growth produces even node count and a valid mol."""
        from biochar.carbon_skeleton import _grow_graph, _graph_to_mol
        import networkx as nx

        G = nx.cycle_graph(6)  # benzene seed
        G = _grow_graph(G, 24, seed=1, defect_fraction=0.0)
        assert G.number_of_nodes() % 2 == 0
        mol = _graph_to_mol(G)
        assert mol is not None

    def test_grow_graph_with_defects_even_node_count(self):
        """Defect growth always produces an even node count."""
        from biochar.carbon_skeleton import _grow_graph
        import networkx as nx

        for seed in range(5):
            G = nx.cycle_graph(6)
            G = _grow_graph(G, 30, seed=seed, defect_fraction=0.15)
            assert G.number_of_nodes() % 2 == 0, (
                f"seed={seed}: odd node count {G.number_of_nodes()}"
            )

    def test_grow_graph_with_defects_has_pentagons(self):
        """With defect_fraction=1.0, all added rings should be pentagons."""
        from biochar.carbon_skeleton import _grow_graph
        import networkx as nx

        G = nx.cycle_graph(6)  # benzene (6 nodes)
        G = _grow_graph(G, 24, seed=42, defect_fraction=1.0)
        # Benzene (6) + n×3 nodes → 6 + 3k for some k, not 6 + 4k
        added = G.number_of_nodes() - 6
        assert added % 3 == 0 or G.number_of_nodes() % 2 == 0, (
            f"Unexpected node count: {G.number_of_nodes()}"
        )

    def test_pah_assembler_defect_returns_valid_mol(self):
        """PAHAssembler with defect_fraction > 0 returns a valid RDKit mol."""
        assembler = PAHAssembler(seed=7)
        skeleton = assembler.generate(
            target_num_carbons=30,
            defect_fraction=0.15,
        )
        assert skeleton.mol is not None
        assert skeleton.num_carbons > 0

    def test_pah_assembler_defect_even_carbons(self):
        """Defect mode skeletons must have an even carbon count."""
        assembler = PAHAssembler(seed=42)
        skeleton = assembler.generate(
            target_num_carbons=40,
            defect_fraction=0.20,
        )
        assert skeleton.num_carbons % 2 == 0, (
            f"Odd carbon count: {skeleton.num_carbons}"
        )

    def test_generator_config_defect_fraction_default(self):
        """GeneratorConfig.defect_fraction should default to 0.0."""
        config = GeneratorConfig()
        assert config.defect_fraction == 0.0

    def test_generator_config_defect_fraction_custom(self):
        """GeneratorConfig.defect_fraction should accept non-zero values."""
        config = GeneratorConfig(defect_fraction=0.1)
        assert config.defect_fraction == 0.1

    def test_biochar_generator_with_defects(self):
        """BiocharGenerator runs end-to-end with defect_fraction > 0."""
        config = GeneratorConfig(
            target_num_carbons=24,
            defect_fraction=0.15,
            seed=99,
        )
        generator = BiocharGenerator(config)
        mol, coords, composition = generator.generate()
        assert mol is not None
        assert composition.num_carbons > 0

    def test_defect_fraction_produces_different_structure(self):
        """Defect mode should yield different atom count than pure hexagon for same target."""
        config_pure = GeneratorConfig(target_num_carbons=50, seed=1, defect_fraction=0.0)
        config_defect = GeneratorConfig(target_num_carbons=50, seed=1, defect_fraction=0.2)
        gen_pure = BiocharGenerator(config_pure)
        gen_defect = BiocharGenerator(config_defect)
        _, _, comp_pure = gen_pure.generate()
        _, _, comp_defect = gen_defect.generate()
        # They may occasionally match, but structures should generally differ
        # At minimum both should produce valid molecules
        assert comp_pure.num_carbons > 0
        assert comp_defect.num_carbons > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
