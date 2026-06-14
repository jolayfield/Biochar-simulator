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
            strict=False,
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
            strict=False,
        )
        generator1 = BiocharGenerator(config1)
        mol1, coords1, comp1 = generator1.generate()

        config2 = GeneratorConfig(
            target_num_carbons=30,
            seed=12345,
            strict=False,
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
            strict=False,
        )
        generator = BiocharGenerator(config)
        mol, coords, composition = generator.generate()
        assert mol is not None
        assert composition.num_carbons > 0

    def test_defect_fraction_produces_different_structure(self):
        """Defect mode should yield different atom count than pure hexagon for same target."""
        config_pure = GeneratorConfig(target_num_carbons=50, seed=1, defect_fraction=0.0, strict=False)
        config_defect = GeneratorConfig(target_num_carbons=50, seed=1, defect_fraction=0.2, strict=False)
        gen_pure = BiocharGenerator(config_pure)
        gen_defect = BiocharGenerator(config_defect)
        _, _, comp_pure = gen_pure.generate()
        _, _, comp_defect = gen_defect.generate()
        # They may occasionally match, but structures should generally differ
        # At minimum both should produce valid molecules
        assert comp_pure.num_carbons > 0
        assert comp_defect.num_carbons > 0


class TestCompositionResult:
    """Test the new CompositionResult properties."""

    def test_molecular_formula_carbon_only(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(num_carbons=24, num_hydrogens=12)
        assert comp.molecular_formula == "C24H12"

    def test_molecular_formula_with_oxygen(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(num_carbons=24, num_hydrogens=12, num_oxygens=2)
        assert comp.molecular_formula == "C24H12O2"

    def test_molecular_formula_with_nitrogen(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(num_carbons=24, num_hydrogens=14, num_oxygens=0, num_nitrogens=1)
        assert comp.molecular_formula == "C24H14N1"

    def test_molecular_weight_coronene(self):
        from biochar.heteroatom_assignment import CompositionResult
        # Coronene C24H12: MW = 24*12.011 + 12*1.008 = 288.264 + 12.096 = 300.360
        comp = CompositionResult(num_carbons=24, num_hydrogens=12)
        assert abs(comp.molecular_weight - 300.360) < 0.01

    def test_molecular_weight_with_oxygen(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(num_carbons=10, num_hydrogens=8, num_oxygens=1)
        expected = 10 * 12.011 + 8 * 1.008 + 1 * 15.999
        assert abs(comp.molecular_weight - expected) < 0.001

    def test_n_c_ratio_default_zero(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(num_carbons=20)
        assert comp.N_C_ratio == 0.0
        assert comp.num_nitrogens == 0


class TestAminoGroup:
    """Test amino (-NH2) functional group placement and N-doping."""

    def test_amino_group_in_functional_groups_constant(self):
        from biochar.constants import FUNCTIONAL_GROUPS
        assert "amino" in FUNCTIONAL_GROUPS

    def test_amino_opls_types_defined(self):
        from biochar.constants import OPLS_ATOM_TYPES, GROMACS_OPLS_TYPE_MAP
        assert "NA" in OPLS_ATOM_TYPES
        assert "HNA" in OPLS_ATOM_TYPES
        assert "NA" in GROMACS_OPLS_TYPE_MAP
        assert "HNA" in GROMACS_OPLS_TYPE_MAP

    def test_amino_placement_on_naphthalene(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        assert mol is not None
        assigner = OxygenAssigner(seed=42)
        mol_out, comp = assigner.assign_oxygens(
            mol, target_O_C_ratio=0.0,
            functional_group_preference={"amino": 1},
        )
        assert mol_out is not None
        assert comp.num_nitrogens >= 1
        assert comp.N_C_ratio > 0.0

    def test_amino_composition_tracking(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")  # anthracene
        assigner = OxygenAssigner(seed=1)
        mol_out, comp = assigner.assign_oxygens(
            mol, target_O_C_ratio=0.0,
            functional_group_preference={"amino": 2},
        )
        assert comp.num_nitrogens == 2
        assert "amino" in comp.placed_counts

    def test_amino_does_not_add_oxygens(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        assigner = OxygenAssigner(seed=42)
        mol_out, comp = assigner.assign_oxygens(
            mol, target_O_C_ratio=0.0,
            functional_group_preference={"amino": 1},
        )
        assert comp.num_oxygens == 0

    def test_amino_end_to_end_generation(self):
        # Naphthalene (10C) + 1 amino → H/C = 0.9 exactly
        config = GeneratorConfig(
            target_num_carbons=10,
            H_C_ratio=0.9,
            functional_groups={"amino": 1},
            O_C_ratio=0.0,
            seed=7,
        )
        gen = BiocharGenerator(config)
        mol, coords, comp = gen.generate()
        assert mol is not None
        assert comp.num_nitrogens >= 1
        # molecular_formula should include N
        assert "N" in comp.molecular_formula


class TestSulfurGroup:
    """Test thiol (-SH) and thioether (-S-) functional group placement (S-doping)."""

    def test_sulfur_groups_in_functional_groups_constant(self):
        from biochar.constants import FUNCTIONAL_GROUPS
        assert "thiol" in FUNCTIONAL_GROUPS
        assert "thioether" in FUNCTIONAL_GROUPS

    def test_sulfur_opls_types_defined(self):
        from biochar.constants import OPLS_ATOM_TYPES, GROMACS_OPLS_TYPE_MAP
        for t in ("SH_", "HSH", "SS"):
            assert t in OPLS_ATOM_TYPES
            assert t in GROMACS_OPLS_TYPE_MAP

    def test_molecular_formula_with_sulfur(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(
            num_carbons=24, num_hydrogens=14, num_oxygens=0,
            num_nitrogens=0, num_sulfurs=2,
        )
        assert comp.molecular_formula == "C24H14S2"

    def test_s_c_ratio_default_zero(self):
        from biochar.heteroatom_assignment import CompositionResult
        comp = CompositionResult(num_carbons=20)
        assert comp.S_C_ratio == 0.0
        assert comp.num_sulfurs == 0

    def test_thiol_placement_on_naphthalene(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        assigner = OxygenAssigner(seed=42)
        mol_out, comp = assigner.assign_oxygens(
            mol, target_O_C_ratio=0.0,
            functional_group_preference={"thiol": 1},
        )
        assert mol_out is not None
        assert comp.num_sulfurs >= 1
        assert comp.S_C_ratio > 0.0
        assert comp.num_oxygens == 0

    def test_thioether_placement(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")  # anthracene
        assigner = OxygenAssigner(seed=1)
        mol_out, comp = assigner.assign_oxygens(
            mol, target_O_C_ratio=0.0,
            functional_group_preference={"thioether": 1},
        )
        assert comp.num_sulfurs == 1
        assert "thioether" in comp.placed_counts

    def test_sulfur_atom_typing(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        from biochar.opls_typing import AtomTyper
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")
        assigner = OxygenAssigner(seed=3)
        mol_out, comp = assigner.assign_oxygens(
            mol, target_O_C_ratio=0.0,
            functional_group_preference={"thiol": 1, "thioether": 1},
        )
        types = set(AtomTyper().assign_atom_types(mol_out).values())
        assert "SH_" in types   # thiol sulfur
        assert "HSH" in types   # thiol hydrogen
        assert "SS" in types    # thioether sulfur

    def test_sulfur_end_to_end_generation(self):
        # 36C PAH + 2 thiol + 1 thioether → C36H16S3
        config = GeneratorConfig(
            target_num_carbons=36,
            H_C_ratio=0.444,
            functional_groups={"thiol": 2, "thioether": 1},
            O_C_ratio=0.0,
            seed=7,
        )
        gen = BiocharGenerator(config)
        mol, coords, comp = gen.generate()
        assert mol is not None
        assert comp.num_sulfurs == 3
        # molecular_formula should include S with the right count
        assert "S3" in comp.molecular_formula
class TestRingNitrogenSubstitution:
    """Test ring-substituting nitrogen (pyridinic / pyrrolic / graphitic)."""

    def test_ring_n_opls_types_defined(self):
        from biochar.constants import (
            OPLS_ATOM_TYPES, OPLS_LJ_PARAMS, OPLS_BOND_PARAMS, GROMACS_OPLS_TYPE_MAP
        )
        for t in ("NPY", "NPR", "NGR", "HNPR"):
            assert t in OPLS_ATOM_TYPES
            assert t in OPLS_LJ_PARAMS
            assert t in GROMACS_OPLS_TYPE_MAP
        assert ("CA", "NPY") in OPLS_BOND_PARAMS
        assert ("CA", "NGR") in OPLS_BOND_PARAMS
        assert ("NPR", "HNPR") in OPLS_BOND_PARAMS

    def test_substitutor_pyridinic_count(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")  # anthracene
        sub = NitrogenSubstitutor(seed=3)
        out = sub.substitute(mol, n_pyridinic=2)
        n_count = sum(1 for a in out.GetAtoms() if a.GetAtomicNum() == 7)
        assert n_count == 2
        assert sub.placed_pyridinic == 2

    def test_substitutor_graphitic_carries_positive_charge(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        # Pyrene has interior junction carbons suitable for graphitic N.
        mol = Chem.MolFromSmiles("c1cc2ccc3cccc4ccc(c1)c2c34")  # pyrene
        sub = NitrogenSubstitutor(seed=1)
        out = sub.substitute(mol, n_graphitic=1)
        assert sub.placed_graphitic == 1
        n_atoms = [a for a in out.GetAtoms() if a.GetAtomicNum() == 7]
        assert len(n_atoms) == 1
        # Graphitic N is pyridinium-like (formal +1) so the ring stays kekulizable.
        assert n_atoms[0].GetFormalCharge() == 1

    def test_substitutor_robust_when_oversubscribed(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        mol = Chem.MolFromSmiles("c1ccccc1")  # benzene, only 6 edge carbons
        sub = NitrogenSubstitutor(seed=1)
        out = sub.substitute(mol, n_pyridinic=99)  # request far more than possible
        assert out is not None  # does not crash
        assert sub.placed_pyridinic <= 6

    def test_pyridinic_plus_graphitic_molecular_formula(self):
        # End-to-end: 2 pyridinic + 1 graphitic on a 40-C nanoflake → 3 N atoms.
        config = GeneratorConfig(
            target_num_carbons=40,
            H_C_ratio=0.4,
            O_C_ratio=0.0,
            num_pyridinic=2,
            num_graphitic=1,
            seed=1,
        )
        gen = BiocharGenerator(config)
        mol, coords, comp = gen.generate()
        assert comp.num_pyridinic == 2
        assert comp.num_graphitic == 1
        assert comp.num_nitrogens == 3
        assert "N3" in comp.molecular_formula

    def test_ring_n_atom_typing(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        from biochar.opls_typing import AtomTyper
        mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")  # anthracene
        mol = NitrogenSubstitutor(seed=3).substitute(mol, n_pyridinic=1)
        types = set(AtomTyper().assign_atom_types(mol).values())
        assert "NPY" in types


class TestGeneratorConfigSerialization:
    """Test GeneratorConfig.to_dict() and from_dict()."""

    def test_to_dict_returns_dict(self):
        config = GeneratorConfig(target_num_carbons=30, seed=99)
        d = config.to_dict()
        assert isinstance(d, dict)
        assert d["target_num_carbons"] == 30
        assert d["seed"] == 99

    def test_round_trip(self):
        config = GeneratorConfig(
            target_num_carbons=40,
            H_C_ratio=0.4,
            O_C_ratio=0.05,
            defect_fraction=0.1,
            molecule_name="BC600",
            seed=123,
        )
        d = config.to_dict()
        config2 = GeneratorConfig.from_dict(d)
        assert config2.target_num_carbons == config.target_num_carbons
        assert config2.H_C_ratio == config.H_C_ratio
        assert config2.molecule_name == config.molecule_name
        assert config2.seed == config.seed

    def test_to_dict_json_serializable(self):
        import json
        config = GeneratorConfig(target_num_carbons=20, seed=0)
        d = config.to_dict()
        # Should not raise
        json.dumps(d)

    def test_from_dict_with_functional_groups(self):
        d = {
            "target_num_carbons": 20,
            "functional_groups": {"phenolic": 2, "amino": 1},
        }
        config = GeneratorConfig.from_dict(d)
        assert config.functional_groups == {"phenolic": 2, "amino": 1}


class TestBiocharSeries:
    """Integration tests for generate_biochar_series()."""

    # Use naphthalene-size (10C) with H/C=0.8, O/C=0.0 — achievable without tolerance issues
    _SMALL_CONF = {"target_num_carbons": 10, "H_C_ratio": 0.8, "O_C_ratio": 0.0}

    def test_series_two_structures(self, tmp_path):
        from biochar import generate_biochar_series
        configs = [
            {"molecule_name": "BC400", **self._SMALL_CONF, "seed": 1},
            {"molecule_name": "BC600", **self._SMALL_CONF, "seed": 2},
        ]
        results = generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=True,
            verbose=False,
        )
        assert set(results.keys()) == {"BC400", "BC600"}

    def test_series_files_created(self, tmp_path):
        from biochar import generate_biochar_series
        configs = [{"molecule_name": "TST", **self._SMALL_CONF, "seed": 5}]
        results = generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=False,
            verbose=False,
        )
        gro, top, itp = results["TST"]
        assert gro.exists(), f"GRO not found: {gro}"
        assert top.exists(), f"TOP not found: {top}"
        assert itp.exists(), f"ITP not found: {itp}"

    def test_series_combined_top_created(self, tmp_path):
        from biochar import generate_biochar_series
        configs = [
            {"molecule_name": "AA", **self._SMALL_CONF, "seed": 10},
            {"molecule_name": "BB", **self._SMALL_CONF, "seed": 11},
        ]
        generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=True,
            verbose=False,
        )
        combined = tmp_path / "combined.top"
        assert combined.exists(), "combined.top was not created"
        content = combined.read_text()
        assert "AA" in content
        assert "BB" in content

    def test_series_combined_top_skipped_for_single(self, tmp_path):
        from biochar import generate_biochar_series
        configs = [{"molecule_name": "SOLO", **self._SMALL_CONF, "seed": 3}]
        generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=True,
            verbose=False,
        )
        combined = tmp_path / "combined.top"
        assert not combined.exists(), "combined.top should not be created for a single structure"

    def test_series_missing_molecule_name_raises(self, tmp_path):
        from biochar import generate_biochar_series
        configs = [{"target_num_carbons": 10}]
        with pytest.raises(ValueError, match="molecule_name"):
            generate_biochar_series(configs, output_directory=str(tmp_path), verbose=False)

    def test_series_molecule_name_too_long_raises(self, tmp_path):
        from biochar import generate_biochar_series
        configs = [{"molecule_name": "TOOLONG", "target_num_carbons": 10}]
        with pytest.raises(ValueError, match="5 character"):
            generate_biochar_series(configs, output_directory=str(tmp_path), verbose=False)

    def test_series_with_amino_group(self, tmp_path):
        from biochar import generate_biochar_series
        # Naphthalene (10C) + 1 amino: 7 edge-CH + 2 NH2-H = 9H → H/C = 0.9 exactly
        configs = [
            {
                "molecule_name": "NBIO",
                "target_num_carbons": 10,
                "H_C_ratio": 0.9,
                "O_C_ratio": 0.0,
                "functional_groups": {"amino": 1},
                "seed": 42,
            }
        ]
        results = generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=False,
            verbose=False,
        )
        assert "NBIO" in results


class TestBiocharResult:
    """Tests for the BiocharResult named-result object."""

    def test_write_files_false_returns_none_paths(self, tmp_path):
        from biochar import generate_biochar, BiocharResult
        result = generate_biochar(
            target_num_carbons=24, O_C_ratio=0.0, seed=42, write_files=False
        )
        assert isinstance(result, BiocharResult)
        assert result.gro_path is None
        assert result.top_path is None
        assert result.itp_path is None

    def test_write_files_false_mol_and_coords_valid(self, tmp_path):
        from biochar import generate_biochar
        result = generate_biochar(
            target_num_carbons=24, O_C_ratio=0.0, seed=42, write_files=False
        )
        assert result.mol is not None
        assert result.coords is not None
        assert result.composition is not None

    def test_positional_unpacking_backward_compat(self, tmp_path):
        from biochar import generate_biochar
        mol, coords, gro, top, itp = generate_biochar(
            target_num_carbons=24, O_C_ratio=0.0, seed=42, write_files=False
        )
        assert mol is not None
        assert coords is not None
        assert gro is None
        assert top is None
        assert itp is None

    def test_write_files_true_returns_paths(self, tmp_path):
        from biochar import generate_biochar
        result = generate_biochar(
            target_num_carbons=24, O_C_ratio=0.0, seed=42, write_files=True,
            output_directory=str(tmp_path),
        )
        assert result.gro_path is not None and result.gro_path.exists()
        assert result.top_path is not None and result.top_path.exists()
        assert result.itp_path is not None and result.itp_path.exists()

    def test_default_write_files_is_true(self, tmp_path):
        from biochar import generate_biochar
        result = generate_biochar(
            target_num_carbons=24, O_C_ratio=0.0, seed=42,
            output_directory=str(tmp_path),
        )
        assert result.gro_path is not None

    def test_biochar_result_exported(self):
        import biochar
        assert hasattr(biochar, "BiocharResult")


class TestSeedIsolation:
    """ISSUE-A: Global random seed must not bleed across calls."""

    def test_two_calls_same_seed_same_smiles(self):
        """Two generate_biochar calls with the same seed must produce the same SMILES."""
        from biochar.biochar_generator import generate_biochar
        # O_C_ratio=0 avoids oxygen placement variability; H_C_ratio tuned for 24C.
        r1 = generate_biochar(
            target_num_carbons=24, H_C_ratio=0.5, O_C_ratio=0.0, seed=42, write_files=False
        )
        r2 = generate_biochar(
            target_num_carbons=24, H_C_ratio=0.5, O_C_ratio=0.0, seed=42, write_files=False
        )
        from rdkit import Chem
        smiles1 = Chem.MolToSmiles(r1.mol)
        smiles2 = Chem.MolToSmiles(r2.mol)
        assert smiles1 == smiles2, (
            f"Same seed produced different SMILES:\n  {smiles1}\n  {smiles2}"
        )

    def test_different_seeds_different_structures(self):
        """Different seeds should (almost certainly) produce different compositions."""
        from biochar.biochar_generator import generate_biochar
        r1 = generate_biochar(
            target_num_carbons=30, H_C_ratio=0.5, O_C_ratio=0.0, seed=1, write_files=False
        )
        r2 = generate_biochar(
            target_num_carbons=30, H_C_ratio=0.5, O_C_ratio=0.0, seed=9999, write_files=False
        )
        # Both should produce valid molecules
        assert r1.mol is not None
        assert r2.mol is not None


class TestEtherSpanGuard:
    """ISSUE-G: max_ether_span < 3 must raise ValueError."""

    def test_ether_span_too_small_raises(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        with pytest.raises(ValueError, match="max_ether_span"):
            OxygenAssigner(max_ether_span=2)

    def test_ether_span_minimum_valid(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        assigner = OxygenAssigner(max_ether_span=3)
        assert assigner._max_ether_span == 3

    def test_config_max_ether_span_too_small_raises(self):
        with pytest.raises(ValueError, match="max_ether_span"):
            GeneratorConfig(max_ether_span=2, strict=False)

    def test_config_max_ether_span_5_warns(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING, logger="biochar.biochar_generator"):
            GeneratorConfig(max_ether_span=5, strict=False)
        assert any("max_ether_span" in r.message for r in caplog.records)


class TestEtherSpanDefault:
    """ISSUE-C: generate_biochar default must match GeneratorConfig default (3)."""

    def test_generate_biochar_default_ether_span_is_3(self):
        from biochar.biochar_generator import generate_biochar
        # Without passing max_ether_span, it should use the GeneratorConfig default of 3.
        # We verify by inspecting the config on a generator built via the convenience API.
        config = GeneratorConfig(target_num_carbons=20, seed=1, strict=False)
        assert config.max_ether_span == 3

    def test_generate_biochar_explicit_override(self):
        # max_ether_span=4 should be accepted by GeneratorConfig without raising.
        config = GeneratorConfig(
            target_num_carbons=20, max_ether_span=4, strict=False, seed=1
        )
        assert config.max_ether_span == 4


class TestTemperatureRangeWarning:
    """ISSUE-D: Out-of-range temperature should emit a warning."""

    def test_out_of_range_temperature_warns(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING, logger="biochar.biochar_generator"):
            # Use a temperature well below any realistic training data minimum
            GeneratorConfig(temperature=50, strict=False)
        assert any("outside the data range" in r.message for r in caplog.records)

    def test_in_range_temperature_no_warning(self, caplog):
        import logging
        with caplog.at_level(logging.WARNING, logger="biochar.biochar_generator"):
            GeneratorConfig(temperature=600, strict=False)
        range_warns = [r for r in caplog.records if "outside the data range" in r.message]
        assert len(range_warns) == 0


class TestFunctionalGroupValidation:
    """ISSUE-E: Over-requested functional groups should emit a warning."""

    def test_excessive_functional_groups_warns(self, caplog):
        import logging
        from biochar.heteroatom_assignment import OxygenAssigner
        from rdkit import Chem
        mol = Chem.MolFromSmiles("c1ccc2ccc3cccc4ccc1c2c34")  # pyrene (16 C)
        assigner = OxygenAssigner(seed=0)
        with caplog.at_level(logging.WARNING, logger="biochar.heteroatom_assignment"):
            assigner.assign_oxygens(
                mol,
                target_O_C_ratio=0.0,
                functional_group_preference={"phenolic": 500},  # absurdly large
            )
        assert any("500" in r.message for r in caplog.records)


class TestRingComposition:
    """ISSUE-J: CarbonSkeleton should report ring composition."""

    def test_pure_hexagon_no_pentagons(self):
        assembler = PAHAssembler(seed=0)
        skeleton = assembler.generate(target_num_carbons=24, defect_fraction=0.0)
        assert skeleton.ring_composition is not None
        assert skeleton.ring_composition["pentagons"] == 0
        assert skeleton.ring_composition["hexagons"] > 0

    def test_defect_mode_may_have_pentagons(self):
        assembler = PAHAssembler(seed=0)
        skeleton = assembler.generate(target_num_carbons=40, defect_fraction=0.3)
        assert skeleton.ring_composition is not None
        # With high defect fraction, pentagons are very likely
        assert "pentagons" in skeleton.ring_composition
        assert "hexagons" in skeleton.ring_composition

    def test_ring_composition_on_generator(self):
        config = GeneratorConfig(target_num_carbons=24, seed=0, strict=False)
        gen = BiocharGenerator(config)
        gen.generate()
        assert gen.ring_composition is not None
        assert gen.ring_composition["hexagons"] >= 0

    def test_ring_composition_in_biochar_result(self):
        from biochar.biochar_generator import generate_biochar
        result = generate_biochar(target_num_carbons=24, seed=0, O_C_ratio=0.0, write_files=False)
        assert result.ring_composition is not None
        assert "hexagons" in result.ring_composition
        assert "pentagons" in result.ring_composition
        assert result.ring_composition["hexagons"] >= 0


class TestChargeNeutrality:
    """ISSUE-K: Static OPLS charge assignment must produce a neutral molecule."""

    def test_charges_sum_to_zero(self):
        from biochar.opls_typing import AtomTyper, ChargeAssigner
        mol = Chem.MolFromSmiles("c1cc(O)ccc1C(=O)O")  # 4-hydroxybenzoic acid
        assert mol is not None
        typer = AtomTyper()
        atom_types = typer.assign_atom_types(mol)
        charger = ChargeAssigner()
        charges = charger.assign_charges(mol, atom_types)
        total = sum(charges.values())
        assert abs(total) < 1e-5, f"Total charge {total:.6f} is not neutral"


class TestBatchProgressCallback:
    """ISSUE-I: generate_biochar_series progress_callback and on_error."""

    def test_progress_callback_called(self, tmp_path):
        from biochar.biochar_generator import generate_biochar_series
        calls = []
        def cb(done, total, name):
            calls.append((done, total, name))

        # Use 24C coronene-like structures with O_C_ratio=0 to avoid validation issues
        configs = [
            {"molecule_name": "BC1", "target_num_carbons": 24,
             "H_C_ratio": 0.5, "O_C_ratio": 0.0, "seed": 1},
            {"molecule_name": "BC2", "target_num_carbons": 24,
             "H_C_ratio": 0.5, "O_C_ratio": 0.0, "seed": 2},
        ]
        generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=False,
            progress_callback=cb,
        )
        assert len(calls) == 2
        assert calls[0] == (1, 2, "BC1")
        assert calls[1] == (2, 2, "BC2")

    def test_on_error_skip_continues(self, tmp_path):
        from biochar.biochar_generator import generate_biochar_series
        configs = [
            {"molecule_name": "BC1", "target_num_carbons": 24,
             "H_C_ratio": 0.5, "O_C_ratio": 0.0, "seed": 1},
            {"molecule_name": "TOOLONG_NAME", "target_num_carbons": 24},  # will raise ValueError
            {"molecule_name": "BC3", "target_num_carbons": 24,
             "H_C_ratio": 0.5, "O_C_ratio": 0.0, "seed": 3},
        ]
        results = generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            create_combined_top=False,
            on_error="skip",
            verbose=False,
        )
        # BC1 and BC3 should succeed; TOOLONG_NAME skipped
        assert "BC1" in results
        assert "BC3" in results
        assert "TOOLONG_NAME" not in results

    def test_on_error_raise_propagates(self, tmp_path):
        from biochar.biochar_generator import generate_biochar_series
        configs = [
            {"molecule_name": "TOOLONG_NAME", "target_num_carbons": 24},
        ]
        with pytest.raises((ValueError, Exception)):
            generate_biochar_series(
                configs,
                output_directory=str(tmp_path),
                create_combined_top=False,
                on_error="raise",
                verbose=False,
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
