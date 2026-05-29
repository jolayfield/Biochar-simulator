"""
Tests for previously uncovered code paths identified by coverage analysis.

Covers: GeneratorConfig roundtrip with box_size, BiocharGenerator pre-generate
errors, print_summary variants, generate_biochar_series verbose mode,
SkeletonValidator rejection cases, NitrogenSubstitutor placement truncation,
valence edge cases (formal charge, unknown element, existing bond rejection),
ValenceReport helpers, OPLS atom typer fallback types, ChargeAssigner private
methods, OPLSPropertyTable unknown-type fallback, constants helper defaults,
_fix_heteroatom_bond_types, _safe_sanitize, and OxygenAssigner spec normalisation.
"""

import pytest
import numpy as np

from rdkit import Chem


# ---------------------------------------------------------------------------
# GeneratorConfig.to_dict / from_dict with box_size  (lines 164-165, 172-173)
# ---------------------------------------------------------------------------

class TestGeneratorConfigRoundtrip:
    def test_to_dict_converts_box_size_to_list(self):
        from biochar.biochar_generator import GeneratorConfig
        cfg = GeneratorConfig(box_size=np.array([5.0, 5.0, 5.0]))
        d = cfg.to_dict()
        assert isinstance(d["box_size"], list)
        assert d["box_size"] == [5.0, 5.0, 5.0]

    def test_from_dict_converts_box_size_list_to_array(self):
        from biochar.biochar_generator import GeneratorConfig
        d = GeneratorConfig().to_dict()
        d["box_size"] = [4.0, 4.0, 4.0]
        cfg = GeneratorConfig.from_dict(d)
        assert isinstance(cfg.box_size, np.ndarray)
        np.testing.assert_array_equal(cfg.box_size, [4.0, 4.0, 4.0])

    def test_to_dict_without_box_size_is_none(self):
        from biochar.biochar_generator import GeneratorConfig
        d = GeneratorConfig().to_dict()
        assert d["box_size"] is None


# ---------------------------------------------------------------------------
# BiocharGenerator pre-generate errors and no-op print_summary (298-299, 319-321)
# ---------------------------------------------------------------------------

class TestBiocharGeneratorPreGenerate:
    def test_export_gromacs_before_generate_raises(self):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig
        gen = BiocharGenerator(GeneratorConfig(seed=1))
        with pytest.raises(RuntimeError, match="generate()"):
            gen.export_gromacs()

    def test_print_summary_before_generate_prints_message(self, capsys):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig
        gen = BiocharGenerator(GeneratorConfig(seed=1))
        gen.print_summary()
        out = capsys.readouterr().out
        assert "No structure generated yet" in out


# ---------------------------------------------------------------------------
# print_summary variants — nitrogen/sulfur counts, ring-N, validation report
# (lines 330-351, 363-369)
# ---------------------------------------------------------------------------

class TestPrintSummaryVariants:
    @pytest.fixture(scope="class")
    def gen_with_pyridinic(self):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig
        # Naphthalene (C10H8): H/C = 0.8 — always passes strict validation.
        # We patch pyridinic/nitrogen counts on the composition object so that
        # print_summary exercises the ring-N and nitrogen print branches.
        cfg = GeneratorConfig(
            target_num_carbons=10,
            H_C_ratio=0.8,
            O_C_ratio=0.0,
            seed=42,
        )
        gen = BiocharGenerator(cfg)
        gen.generate()
        gen.composition.num_nitrogens = 1
        gen.composition.num_pyridinic = 1
        return gen

    def test_print_summary_shows_nitrogen_count(self, gen_with_pyridinic, capsys):
        gen_with_pyridinic.print_summary()
        out = capsys.readouterr().out
        assert "Nitrogens" in out

    def test_print_summary_shows_ring_nitrogen_section(self, gen_with_pyridinic, capsys):
        gen_with_pyridinic.print_summary()
        out = capsys.readouterr().out
        assert "Pyridinic" in out

    def test_print_summary_shows_validation_section(self, gen_with_pyridinic, capsys):
        gen_with_pyridinic.print_summary()
        out = capsys.readouterr().out
        assert "Validation" in out

    def test_print_summary_shows_sulfur(self, capsys):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig
        # Perylene-sized molecule; H/C=0.6 is achievable for 20-carbon PAH.
        cfg = GeneratorConfig(
            target_num_carbons=20,
            H_C_ratio=0.6,
            O_C_ratio=0.0,
            functional_groups={"thiol": 1},
            seed=7,
        )
        gen = BiocharGenerator(cfg)
        gen.generate()
        gen.print_summary()
        out = capsys.readouterr().out
        assert "Sulfurs" in out

    def test_print_summary_shows_errors_and_warnings_when_present(self, capsys):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig
        # Naphthalene (C10H8) has H/C=0.8; use 0.8 to pass validation.
        cfg = GeneratorConfig(target_num_carbons=10, H_C_ratio=0.8, seed=42)
        gen = BiocharGenerator(cfg)
        gen.generate()
        gen.validation_report = (False, ["synthetic error"], ["synthetic warning"], {})
        gen.print_summary()
        out = capsys.readouterr().out
        assert "synthetic error" in out
        assert "synthetic warning" in out


# ---------------------------------------------------------------------------
# generate_biochar_series verbose paths (lines 709-712, 733-735, 754-756,
# 769-778) and _create_combined_topology verbose (801-802, 831-832)
# ---------------------------------------------------------------------------

class TestGenerateBiocharSeriesVerbose:
    def test_verbose_header_and_footer_printed(self, tmp_path, capsys):
        from biochar.biochar_generator import generate_biochar_series
        # Naphthalene (C10H8): H/C=0.8; use that to satisfy strict validation.
        configs = [
            {
                "molecule_name": "BC001",
                "target_num_carbons": 10,
                "H_C_ratio": 0.8,
                "O_C_ratio": 0.0,
                "seed": 1,
            },
        ]
        generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            verbose=True,
            create_combined_top=False,
        )
        out = capsys.readouterr().out
        assert "BATCH BIOCHAR GENERATION" in out
        assert "BC001" in out
        assert "BATCH GENERATION COMPLETE" in out

    def test_combined_top_verbose_mentions_combined_topology(self, tmp_path, capsys):
        from biochar.biochar_generator import generate_biochar_series
        configs = [
            {
                "molecule_name": "BC001",
                "target_num_carbons": 10,
                "H_C_ratio": 0.8,
                "O_C_ratio": 0.0,
                "seed": 1,
            },
            {
                "molecule_name": "BC002",
                "target_num_carbons": 10,
                "H_C_ratio": 0.8,
                "O_C_ratio": 0.0,
                "seed": 2,
            },
        ]
        generate_biochar_series(
            configs,
            output_directory=str(tmp_path),
            verbose=True,
            create_combined_top=True,
        )
        out = capsys.readouterr().out
        assert "Combined topology" in out
        assert "combined.top" in out


# ---------------------------------------------------------------------------
# SkeletonValidator rejection cases (carbon_skeleton.py lines 908-930)
# ---------------------------------------------------------------------------

class TestSkeletonValidatorRejections:
    def test_validate_none_mol_fails(self):
        from biochar.carbon_skeleton import CarbonSkeleton, SkeletonValidator
        skeleton = CarbonSkeleton(
            mol=None, smiles="", num_carbons=0,
            num_aromatic_carbons=0, aromaticity_percent=0.0,
        )
        valid, errors = SkeletonValidator.validate(skeleton)
        assert not valid
        assert any("None" in e for e in errors)

    def test_validate_zero_carbons_fails(self):
        from biochar.carbon_skeleton import CarbonSkeleton, SkeletonValidator
        mol = Chem.MolFromSmiles("C")
        skeleton = CarbonSkeleton(
            mol=mol, smiles="C", num_carbons=0,
            num_aromatic_carbons=0, aromaticity_percent=0.0,
        )
        valid, errors = SkeletonValidator.validate(skeleton)
        assert not valid
        assert any("carbon" in e.lower() for e in errors)

    def test_validate_aromatic_carbons_exceed_total_fails(self):
        from biochar.carbon_skeleton import CarbonSkeleton, SkeletonValidator
        mol = Chem.MolFromSmiles("c1ccccc1")
        skeleton = CarbonSkeleton(
            mol=mol, smiles="c1ccccc1", num_carbons=6,
            num_aromatic_carbons=7,  # impossible: > num_carbons
            aromaticity_percent=95.0,
        )
        valid, errors = SkeletonValidator.validate(skeleton)
        assert not valid
        assert any("aromatic" in e.lower() for e in errors)

    def test_validate_aromaticity_out_of_range_fails(self):
        from biochar.carbon_skeleton import CarbonSkeleton, SkeletonValidator
        mol = Chem.MolFromSmiles("c1ccccc1")
        skeleton = CarbonSkeleton(
            mol=mol, smiles="c1ccccc1", num_carbons=6,
            num_aromatic_carbons=6,
            aromaticity_percent=110.0,  # > 100 is invalid
        )
        valid, errors = SkeletonValidator.validate(skeleton)
        assert not valid
        assert any("aromaticity" in e.lower() or "110" in e for e in errors)

    def test_validate_disconnected_graph_fails(self):
        from biochar.carbon_skeleton import CarbonSkeleton, SkeletonValidator
        # Two disconnected benzene rings in one mol object
        mol = Chem.MolFromSmiles("c1ccccc1.c1ccccc1")
        skeleton = CarbonSkeleton(
            mol=mol, smiles="c1ccccc1.c1ccccc1", num_carbons=12,
            num_aromatic_carbons=12, aromaticity_percent=100.0,
        )
        valid, errors = SkeletonValidator.validate(skeleton)
        assert not valid
        assert any("connected" in e.lower() for e in errors)


# ---------------------------------------------------------------------------
# NitrogenSubstitutor placement truncation (heteroatom_assignment 578-601)
# ---------------------------------------------------------------------------

class TestNDopingPlacementTruncation:
    """Substitutor places as many N as possible and records actual count."""

    def test_pyridinic_truncation_fewer_sites_than_requested(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        mol = Chem.MolFromSmiles("c1ccccc1")  # 6 edge carbons max
        sub = NitrogenSubstitutor(seed=1)
        result = sub.substitute(mol, n_pyridinic=100)
        assert result is not None
        assert sub.placed_pyridinic < 100

    def test_graphitic_zero_placed_when_no_interior_carbons(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        # Benzene has no interior carbons (all are edge) → 0 graphitic N possible
        mol = Chem.MolFromSmiles("c1ccccc1")
        sub = NitrogenSubstitutor(seed=1)
        result = sub.substitute(mol, n_graphitic=5)
        assert result is not None
        assert sub.placed_graphitic == 0

    def test_pyrrolic_zero_placed_when_no_five_membered_rings(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        # Naphthalene has only 6-membered rings — no 5-ring sites
        mol = Chem.MolFromSmiles("c1ccc2ccccc2c1")
        sub = NitrogenSubstitutor(seed=1)
        result = sub.substitute(mol, n_pyrrolic=10)
        assert result is not None
        assert sub.placed_pyrrolic == 0

    def test_placed_counts_always_set_even_on_truncation(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor
        mol = Chem.MolFromSmiles("c1ccc2ccc3ccccc3c2c1")  # anthracene — no 5-rings
        sub = NitrogenSubstitutor(seed=1)
        sub.substitute(mol, n_pyridinic=50, n_pyrrolic=50, n_graphitic=50)
        assert sub.placed_pyridinic <= 50
        assert sub.placed_pyrrolic == 0     # no 5-rings
        assert sub.placed_graphitic >= 0    # may place some on interior atoms


# ---------------------------------------------------------------------------
# get_valence_range with formal charge and unknown element (valence.py 63-78)
# ---------------------------------------------------------------------------

class TestGetValenceRange:
    def test_unknown_element_returns_default_1_4(self):
        from biochar.valence import get_valence_range
        min_v, max_v = get_valence_range(100)  # Fermium — not in STANDARD_VALENCES
        assert (min_v, max_v) == (1, 4)

    def test_nitrogen_positive_charge_extends_max_to_4(self):
        from biochar.valence import get_valence_range
        min_v, max_v = get_valence_range(7, formal_charge=1)
        assert max_v == 4  # ammonium-like: extended to 4

    def test_nitrogen_negative_charge_lowers_min_valence(self):
        from biochar.valence import get_valence_range
        min_v, max_v = get_valence_range(7, formal_charge=-1)
        standard_min = 3
        assert min_v < standard_min

    def test_sulfur_positive_charge_extends_max(self):
        from biochar.valence import get_valence_range
        min_v, max_v = get_valence_range(16, formal_charge=1)
        assert max_v > 2  # extended beyond standard S max of 2


# ---------------------------------------------------------------------------
# SafeBondAdder.can_add_bond rejects existing bonds (valence.py 273-274)
# ---------------------------------------------------------------------------

class TestSafeBondAdderExistingBond:
    def test_existing_bond_rejected(self):
        from biochar.valence import SafeBondAdder
        mol = Chem.MolFromSmiles("CC")
        can_add, reason = SafeBondAdder.can_add_bond(mol, 0, 1, bond_type=1)
        assert not can_add
        assert "already exists" in reason

    def test_nonexistent_bond_can_be_added(self):
        from biochar.valence import SafeBondAdder
        mol = Chem.MolFromSmiles("CCC")  # C0-C1-C2; no bond between C0 and C2
        can_add, _ = SafeBondAdder.can_add_bond(mol, 0, 2, bond_type=1)
        assert can_add


# ---------------------------------------------------------------------------
# ValenceReport helpers (valence.py lines 368-404)
# ---------------------------------------------------------------------------

class TestValenceReport:
    def test_get_summary_returns_expected_fields(self):
        from biochar.valence import ValenceReport
        mol = Chem.MolFromSmiles("CC")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert "total_atoms" in summary
        assert "valid_atoms" in summary
        assert "all_valid" in summary
        assert summary["element_counts"]["C"] == 2

    def test_get_summary_all_valid_for_valid_molecule(self):
        from biochar.valence import ValenceReport
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        summary = ValenceReport.get_summary(mol)
        assert summary["all_valid"]

    def test_print_summary_produces_output(self, capsys):
        from biochar.valence import ValenceReport
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        ValenceReport.print_summary(mol)
        out = capsys.readouterr().out
        assert "Valence Summary" in out


# ---------------------------------------------------------------------------
# ValenceValidator.print_valence_report (valence.py lines 215-244)
# ---------------------------------------------------------------------------

class TestValenceValidatorPrintReport:
    def test_print_valence_report_produces_output(self, capsys):
        from biochar.valence import ValenceValidator
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        ValenceValidator.print_valence_report(mol)
        out = capsys.readouterr().out
        assert "VALENCE VALIDATION REPORT" in out
        assert "OK" in out


# ---------------------------------------------------------------------------
# AtomTyper fallback type paths (opls_typing.py lines 83, 104-106, 164, 174)
# ---------------------------------------------------------------------------

class TestAtomTyperFallbacks:
    def test_aliphatic_carbon_typed_as_ct(self):
        from biochar.opls_typing import AtomTyper
        mol = Chem.MolFromSmiles("CC")  # ethane — all sp3
        mol = Chem.AddHs(mol)
        types = AtomTyper().assign_atom_types(mol)
        assert "CT" in types.values()

    def test_h_on_aliphatic_c_typed_as_hc(self):
        from biochar.opls_typing import AtomTyper
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        types = AtomTyper().assign_atom_types(mol)
        h_types = {t for i, t in types.items()
                   if mol.GetAtomWithIdx(i).GetAtomicNum() == 1}
        assert h_types == {"HC"}

    def test_n_with_three_nonaromatic_bonds_typed_as_n(self):
        from biochar.opls_typing import AtomTyper
        # Trimethylamine: N has 3 bonds to aliphatic C, no ring membership
        mol = Chem.MolFromSmiles("N(C)(C)C")
        mol = Chem.AddHs(mol)
        types = AtomTyper().assign_atom_types(mol)
        n_types = {t for i, t in types.items()
                   if mol.GetAtomWithIdx(i).GetAtomicNum() == 7}
        assert "N" in n_types

    def test_n_with_two_bonds_typed_as_nt(self):
        from biochar.opls_typing import AtomTyper
        # [NH2-]: N has formal -1 charge, 2 H bonds → not matching any special case
        mol = Chem.MolFromSmiles("[NH2-]")
        mol = Chem.AddHs(mol)
        types = AtomTyper().assign_atom_types(mol)
        n_types = {t for i, t in types.items()
                   if mol.GetAtomWithIdx(i).GetAtomicNum() == 7}
        assert "NT" in n_types

    def test_unknown_element_returns_x_n_fallback(self):
        from biochar.opls_typing import AtomTyper
        # Xenon (54) passes through all element-specific branches
        rw = Chem.RWMol()
        rw.AddAtom(Chem.Atom(54))
        mol = rw.GetMol()
        t = AtomTyper()._determine_atom_type(mol, mol.GetAtomWithIdx(0))
        assert t == "X54"


# ---------------------------------------------------------------------------
# ChargeAssigner._estimate_charge and _equilibrate_charges (lines 238-305)
# ---------------------------------------------------------------------------

class TestChargeAssignerPrivateMethods:
    def test_estimate_charge_carbon_returns_typical_value(self):
        from biochar.opls_typing import ChargeAssigner
        mol = Chem.MolFromSmiles("CC")
        charge = ChargeAssigner()._estimate_charge(mol, 0, "XUNKNOWN")
        assert charge == pytest.approx(-0.15)

    def test_estimate_charge_oxygen(self):
        from biochar.opls_typing import ChargeAssigner
        mol = Chem.MolFromSmiles("O")
        charge = ChargeAssigner()._estimate_charge(mol, 0, "XUNKNOWN")
        assert charge == pytest.approx(-0.50)

    def test_estimate_charge_nitrogen(self):
        from biochar.opls_typing import ChargeAssigner
        mol = Chem.MolFromSmiles("N")
        charge = ChargeAssigner()._estimate_charge(mol, 0, "XUNKNOWN")
        assert charge == pytest.approx(-0.30)

    def test_estimate_charge_sulfur(self):
        from biochar.opls_typing import ChargeAssigner
        mol = Chem.MolFromSmiles("S")
        charge = ChargeAssigner()._estimate_charge(mol, 0, "XUNKNOWN")
        assert charge == pytest.approx(-0.20)

    def test_estimate_charge_hydrogen(self):
        from biochar.opls_typing import ChargeAssigner
        mol = Chem.MolFromSmiles("[H][H]")
        charge = ChargeAssigner()._estimate_charge(mol, 0, "XUNKNOWN")
        assert charge == pytest.approx(0.10)

    def test_estimate_charge_unknown_element_returns_zero(self):
        from biochar.opls_typing import ChargeAssigner
        rw = Chem.RWMol()
        rw.AddAtom(Chem.Atom(54))  # Xe
        mol = rw.GetMol()
        charge = ChargeAssigner()._estimate_charge(mol, 0, "XUNKNOWN")
        assert charge == pytest.approx(0.0)

    def test_equilibrate_already_neutral_unchanged(self):
        from biochar.opls_typing import ChargeAssigner
        mol = Chem.MolFromSmiles("CC")
        charges = {0: 0.1, 1: -0.1}  # sums to 0
        result = ChargeAssigner()._equilibrate_charges(mol, charges)
        assert result == charges

    def test_equilibrate_imbalanced_charge_corrected_to_zero(self):
        from biochar.opls_typing import ChargeAssigner
        # Water molecule: give O a non-zero charge with no H's to spread to
        mol = Chem.MolFromSmiles("O")
        charges = {0: 0.4}
        result = ChargeAssigner()._equilibrate_charges(mol, charges)
        assert abs(sum(result.values())) < 1e-6


# ---------------------------------------------------------------------------
# OPLSPropertyTable unknown-type mass fallback (opls_typing.py lines 330-332)
# ---------------------------------------------------------------------------

class TestOPLSPropertyTableUnknownType:
    def test_unknown_type_uses_rdkit_atom_mass(self):
        from biochar.opls_typing import OPLSPropertyTable
        mol = Chem.MolFromSmiles("C")
        mol = Chem.AddHs(mol)
        n = mol.GetNumAtoms()
        atom_types = {i: "XUNKNOWN" for i in range(n)}
        charges = {i: 0.0 for i in range(n)}
        table = OPLSPropertyTable(mol, atom_types, charges)
        for prop in table.get_properties():
            rdkit_mass = mol.GetAtomWithIdx(prop.atom_idx).GetMass()
            assert prop.mass == pytest.approx(rdkit_mass, abs=0.1)


# ---------------------------------------------------------------------------
# constants.py helper defaults for unknown inputs (lines 464-478)
# ---------------------------------------------------------------------------

class TestConstantsHelperFallbacks:
    def test_get_atom_mass_unknown_type_falls_back_to_carbon(self):
        from biochar.constants import get_atom_mass
        assert get_atom_mass("XYZUNKNOWN") == pytest.approx(12.01)

    def test_get_vdw_radius_unknown_element_returns_default(self):
        from biochar.constants import get_vdw_radius
        assert get_vdw_radius("Xx") == pytest.approx(1.70)

    def test_get_covalent_radius_unknown_element_returns_default(self):
        from biochar.constants import get_covalent_radius
        assert get_covalent_radius("Xx") == pytest.approx(0.76)


# ---------------------------------------------------------------------------
# _fix_heteroatom_bond_types (heteroatom_assignment.py lines 111-157)
# ---------------------------------------------------------------------------

class TestFixHeteroatomBondTypes:
    def test_spurious_aromatic_co_bond_reset_to_single(self):
        from biochar.heteroatom_assignment import _fix_heteroatom_bond_types
        # Build benzene with a pendant O, then manually mark C-O bond as AROMATIC.
        rw = Chem.RWMol(Chem.MolFromSmiles("c1ccccc1"))
        o_idx = rw.AddAtom(Chem.Atom(8))
        c_idx = 0
        rw.AddBond(c_idx, o_idx, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        # Now set the C-O bond type to AROMATIC in a fresh RWMol
        rw2 = Chem.RWMol(mol)
        rw2.GetBondBetweenAtoms(c_idx, o_idx).SetBondType(Chem.BondType.AROMATIC)
        rw2.GetAtomWithIdx(o_idx).SetIsAromatic(True)
        mol_aromatic_co = rw2.GetMol()
        assert mol_aromatic_co.GetBondBetweenAtoms(c_idx, o_idx).GetBondTypeAsDouble() == 1.5

        fixed = _fix_heteroatom_bond_types(mol_aromatic_co)
        assert fixed.GetBondBetweenAtoms(c_idx, o_idx).GetBondTypeAsDouble() != 1.5

    def test_aromatic_ring_nitrogen_bonds_left_intact(self):
        from biochar.heteroatom_assignment import _fix_heteroatom_bond_types
        mol = Chem.MolFromSmiles("c1ccncc1")  # pyridine
        Chem.SanitizeMol(mol)
        fixed = _fix_heteroatom_bond_types(mol)
        n_idx = next(a.GetIdx() for a in fixed.GetAtoms() if a.GetAtomicNum() == 7)
        for bond in fixed.GetAtomWithIdx(n_idx).GetBonds():
            assert bond.GetBondTypeAsDouble() == 1.5  # still AROMATIC

    def test_already_correct_molecule_returned_unchanged(self):
        from biochar.heteroatom_assignment import _fix_heteroatom_bond_types
        mol = Chem.MolFromSmiles("c1ccccc1")  # no heteroatoms
        fixed = _fix_heteroatom_bond_types(mol)
        assert fixed is mol  # no changes → same object returned


# ---------------------------------------------------------------------------
# _safe_sanitize (heteroatom_assignment.py lines 99-108)
# ---------------------------------------------------------------------------

class TestSafeSanitize:
    def test_does_not_raise_on_simple_mol(self):
        from biochar.heteroatom_assignment import _safe_sanitize
        mol = Chem.MolFromSmiles("c1ccccc1")
        _safe_sanitize(mol)  # must not raise

    def test_does_not_raise_on_rw_mol(self):
        from biochar.heteroatom_assignment import _safe_sanitize
        mol = Chem.RWMol(Chem.MolFromSmiles("c1ccc2ccccc2c1"))
        _safe_sanitize(mol)  # must not raise

    def test_does_not_raise_on_mol_with_heteroatoms(self):
        from biochar.heteroatom_assignment import _safe_sanitize
        mol = Chem.MolFromSmiles("c1ccncc1")  # pyridine
        _safe_sanitize(mol)  # must not raise


# ---------------------------------------------------------------------------
# OxygenAssigner._validate_and_normalise_spec (lines 293-323)
# ---------------------------------------------------------------------------

class TestValidateAndNormaliseSpec:
    def test_unknown_group_keys_removed(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        result = OxygenAssigner(seed=1)._validate_and_normalise_spec(
            {"phenolic": 2, "bogusgroup": 3}
        )
        assert "bogusgroup" not in result
        assert result.get("phenolic") == 2

    def test_fallback_groups_mapped_to_phenolic(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        # "carbonyl" and "quinone" both map to "phenolic" (see _FALLBACK_GROUPS)
        result = OxygenAssigner(seed=1)._validate_and_normalise_spec(
            {"carbonyl": 1, "quinone": 2}
        )
        assert "carbonyl" not in result
        assert "quinone" not in result
        assert result.get("phenolic", 0) == 3

    def test_all_unknown_keys_returns_empty_dict(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        result = OxygenAssigner(seed=1)._validate_and_normalise_spec(
            {"fake1": 1, "fake2": 2}
        )
        assert result == {}

    def test_lactone_mapped_to_phenolic(self):
        from biochar.heteroatom_assignment import OxygenAssigner
        result = OxygenAssigner(seed=1)._validate_and_normalise_spec(
            {"lactone": 3, "phenolic": 1}
        )
        assert result.get("phenolic", 0) == 4  # 3 lactone + 1 phenolic


# ---------------------------------------------------------------------------
# validation.py gaps (lines 152, 163-164, 171, 221, 263)
# ---------------------------------------------------------------------------

class TestValidationModuleGaps:
    def test_chemical_feasibility_unusual_atom_generates_warning(self):
        from biochar.validation import ChemicalFeasibilityValidator
        # Gold (79): atomic_num > 18 and not in [16, 17, 35]
        rw = Chem.RWMol()
        c_idx = rw.AddAtom(Chem.Atom(6))
        h_idx = rw.AddAtom(Chem.Atom(1))
        au_idx = rw.AddAtom(Chem.Atom(79))
        rw.AddBond(c_idx, h_idx, Chem.BondType.SINGLE)
        rw.AddBond(c_idx, au_idx, Chem.BondType.SINGLE)
        mol = rw.GetMol()
        report = ChemicalFeasibilityValidator.validate(mol)
        assert any("Unusual atom" in w or "Au" in w for w in report.warnings)

    def test_chemical_feasibility_high_total_charge_generates_warning(self):
        from biochar.validation import ChemicalFeasibilityValidator
        # Three ammonium ions: total formal charge = 3 > 2
        mol = Chem.MolFromSmiles("[NH4+].[NH4+].[NH4+]")
        report = ChemicalFeasibilityValidator.validate(mol)
        assert any("charge" in w.lower() for w in report.warnings)

    def test_structure_validator_is_connected_empty_mol_returns_true(self):
        from biochar.validation import StructureValidator
        # _is_connected with 0-atom mol returns True (line 263)
        mol = Chem.RWMol().GetMol()
        assert StructureValidator._is_connected(mol)

    def test_structure_validator_large_molecule_no_warning(self):
        from biochar.validation import StructureValidator
        # num_atoms > 5000 triggers "very large" warning (line 221); skip —
        # small mol should produce no such warning
        mol = Chem.MolFromSmiles("c1ccccc1")
        mol = Chem.AddHs(mol)
        report = StructureValidator.validate(mol)
        assert not any("very large" in w for w in report.warnings)
