"""
Tests for ML-based partial charge refinement (issue #4).
"""

import pytest
import numpy as np
from rdkit import Chem

from biochar.ml_charges import MLChargeRefinement, _generate_training_data
from biochar.opls_typing import AtomTyper, ChargeAssigner
from biochar.biochar_generator import BiocharGenerator, GeneratorConfig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _benzene_with_types():
    mol = Chem.MolFromSmiles("c1ccccc1")
    mol = Chem.AddHs(mol)
    typer = AtomTyper()
    atom_types = typer.assign_atom_types(mol)
    return mol, atom_types


def _naphthalene_with_types():
    mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
    mol = Chem.AddHs(mol)
    typer = AtomTyper()
    atom_types = typer.assign_atom_types(mol)
    return mol, atom_types


# ---------------------------------------------------------------------------
# MLChargeRefinement — basic interface
# ---------------------------------------------------------------------------

class TestMLChargeRefinementInterface:
    def test_instantiation_uses_bundled_model(self):
        refiner = MLChargeRefinement()
        assert refiner._model is not None

    def test_custom_model_path_fallback(self, tmp_path):
        """Missing model path falls back to on-the-fly GPR without error."""
        refiner = MLChargeRefinement(model_path=tmp_path / "nonexistent.pkl")
        assert refiner._model is not None

    def test_refine_returns_dict(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, types)
        assert isinstance(charges, dict)
        assert set(charges.keys()) == set(range(mol.GetNumAtoms()))

    def test_refine_all_atoms_present(self):
        mol, types = _naphthalene_with_types()
        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, types)
        assert len(charges) == mol.GetNumAtoms()

    def test_refine_charge_neutrality(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, types)
        assert abs(sum(charges.values())) < 1e-6

    def test_refine_charge_neutrality_naphthalene(self):
        mol, types = _naphthalene_with_types()
        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, types)
        assert abs(sum(charges.values())) < 1e-6

    def test_refine_all_floats(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, types)
        for v in charges.values():
            assert isinstance(v, float)

    def test_refine_no_extreme_charges(self):
        """Predicted charges should be physically plausible (< 2 e)."""
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        charges = refiner.refine(mol, types)
        for v in charges.values():
            assert abs(v) < 2.0


# ---------------------------------------------------------------------------
# Featuriser
# ---------------------------------------------------------------------------

class TestFeaturiser:
    def test_featurize_shape(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        X = refiner._featurize(mol, types)
        assert X.shape == (mol.GetNumAtoms(), 8)

    def test_featurize_dtype_float(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        X = refiner._featurize(mol, types)
        assert X.dtype == float

    def test_featurize_aromatic_carbons(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        X = refiner._featurize(mol, types)
        # Benzene has 6 aromatic C (atomic_num=6) and 6 aromatic H
        c_rows = X[X[:, 0] == 6]
        assert len(c_rows) == 6
        assert all(c_rows[:, 1] == 1)  # is_aromatic
        assert all(c_rows[:, 2] == 1)  # in_ring

    def test_featurize_ring_size_benzene(self):
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement()
        X = refiner._featurize(mol, types)
        c_rows = X[X[:, 0] == 6]
        assert all(c_rows[:, 3] == 6)  # smallest_ring_size = 6

    def test_opls_group_carbons(self):
        for t in ("CA", "CT", "C"):
            assert MLChargeRefinement._opls_group(t) == 0

    def test_opls_group_hydrogens(self):
        for t in ("HA", "HO", "HC", "HNA", "HSH", "HNPR"):
            assert MLChargeRefinement._opls_group(t) == 1

    def test_opls_group_oxygens(self):
        for t in ("OH", "OS", "OC", "O", "OH2"):
            assert MLChargeRefinement._opls_group(t) == 2

    def test_opls_group_nitrogens(self):
        for t in ("NA", "NPY", "NPR", "NGR", "N", "NT"):
            assert MLChargeRefinement._opls_group(t) == 3

    def test_opls_group_sulfurs(self):
        for t in ("SH_", "SS"):
            assert MLChargeRefinement._opls_group(t) == 4

    def test_opls_group_unknown(self):
        assert MLChargeRefinement._opls_group("XYZ") == 5


# ---------------------------------------------------------------------------
# Training data
# ---------------------------------------------------------------------------

class TestTrainingData:
    def test_generate_returns_arrays(self):
        X, y = _generate_training_data()
        assert isinstance(X, np.ndarray)
        assert isinstance(y, np.ndarray)

    def test_generate_shapes_consistent(self):
        X, y = _generate_training_data()
        assert X.shape[0] == y.shape[0]
        assert X.shape[1] == 8

    def test_generate_non_empty(self):
        X, y = _generate_training_data()
        assert X.shape[0] > 0

    def test_generate_charges_plausible(self):
        _, y = _generate_training_data()
        assert np.all(np.abs(y) < 1.5)


# ---------------------------------------------------------------------------
# train_and_save
# ---------------------------------------------------------------------------

class TestTrainAndSave:
    def test_train_and_save_creates_file(self, tmp_path):
        X, y = _generate_training_data()
        out = tmp_path / "test_model.pkl"
        MLChargeRefinement.train_and_save(X, y, output_path=out)
        assert out.exists()

    def test_train_and_save_returns_refiner(self, tmp_path):
        X, y = _generate_training_data()
        out = tmp_path / "test_model.pkl"
        instance = MLChargeRefinement.train_and_save(X, y, output_path=out)
        assert isinstance(instance, MLChargeRefinement)

    def test_custom_model_loads_and_predicts(self, tmp_path):
        X, y = _generate_training_data()
        out = tmp_path / "custom.pkl"
        MLChargeRefinement.train_and_save(X, y, output_path=out)
        mol, types = _benzene_with_types()
        refiner = MLChargeRefinement(model_path=out)
        charges = refiner.refine(mol, types)
        assert abs(sum(charges.values())) < 1e-6


# ---------------------------------------------------------------------------
# Integration: GeneratorConfig.charge_method
# ---------------------------------------------------------------------------

class TestGeneratorConfigChargeMethod:
    def test_default_is_opls(self):
        config = GeneratorConfig()
        assert config.charge_method == "opls"

    def test_accepts_ml(self):
        config = GeneratorConfig(charge_method="ml")
        assert config.charge_method == "ml"

    def test_rejects_invalid(self):
        with pytest.raises(ValueError, match="charge_method"):
            GeneratorConfig(charge_method="dft")

    def test_to_dict_includes_charge_method(self):
        config = GeneratorConfig(charge_method="ml")
        d = config.to_dict()
        assert d["charge_method"] == "ml"

    def test_from_dict_roundtrip(self):
        config = GeneratorConfig(charge_method="ml")
        d = config.to_dict()
        config2 = GeneratorConfig.from_dict(d)
        assert config2.charge_method == "ml"


# ---------------------------------------------------------------------------
# Integration: end-to-end generate() with charge_method="ml"
# ---------------------------------------------------------------------------

class TestMLChargesEndToEnd:
    def test_generate_with_opls_unchanged(self):
        """OPLS path still works; charges are non-empty dict."""
        config = GeneratorConfig(
            target_num_carbons=6, H_C_ratio=1.0, O_C_ratio=0.0,
            charge_method="opls", seed=42, strict=False,
        )
        gen = BiocharGenerator(config)
        mol, coords, comp = gen.generate()
        assert isinstance(gen.charges, dict)
        assert len(gen.charges) == mol.GetNumAtoms()

    def test_generate_with_ml_returns_charges(self):
        config = GeneratorConfig(
            target_num_carbons=6, H_C_ratio=1.0, O_C_ratio=0.0,
            charge_method="ml", seed=42, strict=False,
        )
        gen = BiocharGenerator(config)
        mol, coords, comp = gen.generate()
        assert isinstance(gen.charges, dict)
        assert len(gen.charges) == mol.GetNumAtoms()

    def test_generate_with_ml_charge_neutrality(self):
        config = GeneratorConfig(
            target_num_carbons=6, H_C_ratio=1.0, O_C_ratio=0.0,
            charge_method="ml", seed=42, strict=False,
        )
        gen = BiocharGenerator(config)
        gen.generate()
        assert abs(sum(gen.charges.values())) < 1e-5

    def test_ml_and_opls_atom_count_identical(self):
        """Both paths produce charges for exactly the same number of atoms."""
        base = dict(target_num_carbons=6, H_C_ratio=1.0, O_C_ratio=0.0,
                    seed=42, strict=False)
        gen_opls = BiocharGenerator(GeneratorConfig(charge_method="opls", **base))
        gen_opls.generate()
        gen_ml = BiocharGenerator(GeneratorConfig(charge_method="ml", **base))
        gen_ml.generate()
        assert set(gen_opls.charges.keys()) == set(gen_ml.charges.keys())
