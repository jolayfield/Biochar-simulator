"""
End-to-end pH wiring: GeneratorConfig -> pipeline -> CLI -> sweep.

The most important test here is that omitting pH reproduces the pre-pH
generator exactly. Everything else is additive; that one is a promise.
"""

import pytest
from rdkit import Chem

from biochar import BiocharGenerator, GeneratorConfig


def gen(**kwargs):
    """
    Generate a structure for pH assertions.

    strict=False with wide tolerances on purpose: these tests are about the pH
    wiring, not about hitting composition targets. Explicit functional_groups
    override O_C_ratio for *placement* but are still checked against it during
    validation, and the flat hex-lattice path reports expected steric-clash
    artefacts -- neither says anything about protonation. Composition accuracy
    is covered by test_generator.py.
    """
    base = dict(
        target_num_carbons=36,
        H_C_ratio=0.45,
        O_C_ratio=0.15,
        seed=11,
        molecule_name="BCT",
        strict=False,
        H_C_tolerance=1.0,
        O_C_tolerance=1.0,
    )
    base.update(kwargs)
    g = BiocharGenerator(GeneratorConfig(**base))
    mol, coords, comp = g.generate()
    return g, mol, comp


def net_formal(mol) -> int:
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


class TestDefaultPathUnchanged:
    """R1: omitting pH must reproduce today's behaviour exactly."""

    def test_ph_defaults_to_none(self):
        assert GeneratorConfig().pH is None

    def test_no_ph_produces_a_neutral_structure(self):
        _, mol, comp = gen(functional_groups={"carboxyl": 3, "phenolic": 2})
        assert net_formal(mol) == 0
        assert comp.net_charge == 0

    def test_no_ph_leaves_the_protonation_census_empty(self):
        _, _, comp = gen(functional_groups={"carboxyl": 3})
        assert comp.ionized_counts == {}
        assert comp.titratable_counts == {}

    def test_no_ph_is_bit_identical_across_runs(self):
        _, a, _ = gen(functional_groups={"carboxyl": 2, "phenolic": 2})
        _, b, _ = gen(functional_groups={"carboxyl": 2, "phenolic": 2})
        assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)


class TestPhDrivesNetCharge:
    def test_neutral_ph_makes_a_carboxyl_structure_anionic(self):
        _, mol, comp = gen(functional_groups={"carboxyl": 3}, pH=7.0)
        assert net_formal(mol) < 0
        assert comp.net_charge < 0

    def test_low_ph_leaves_carboxyls_protonated(self):
        _, mol, comp = gen(functional_groups={"carboxyl": 3}, pH=1.0)
        assert comp.net_charge >= 0

    def test_high_ph_deprotonates_phenolics_too(self):
        _, _, comp = gen(functional_groups={"phenolic": 4}, pH=13.0)
        assert comp.net_charge < 0

    def test_neutral_ph_leaves_phenolics_alone(self):
        _, _, comp = gen(functional_groups={"phenolic": 4}, pH=7.0)
        assert comp.net_charge == 0

    def test_titration_series_is_monotonic(self):
        charges = [
            gen(functional_groups={"carboxyl": 3}, pH=p)[2].net_charge
            for p in (1.0, 3.0, 5.0, 7.0)
        ]
        assert charges == sorted(charges, reverse=True), charges

    def test_census_is_reported(self):
        _, _, comp = gen(functional_groups={"carboxyl": 3}, pH=7.0)
        assert comp.titratable_counts.get("carboxyl", 0) > 0
        assert comp.ionized_counts.get("carboxyl", 0) > 0


class TestPhReachesTheChargeModel:
    def test_partial_charges_sum_to_the_net_formal_charge(self):
        g, mol, comp = gen(functional_groups={"carboxyl": 3}, pH=7.0)
        assert comp.net_charge != 0, "fixture should be charged"
        total = sum(g.charges.values())
        assert total == pytest.approx(float(comp.net_charge), abs=1e-4)

    def test_no_atom_is_left_unrecognised(self):
        g, _, _ = gen(functional_groups={"carboxyl": 2, "phenolic": 2}, pH=7.0)
        assert not [t for t in g.atom_types.values() if t.startswith("X")]


class TestConfigValidation:
    @pytest.mark.parametrize("pH", [-0.5, 14.5, 100.0])
    def test_out_of_range_ph_is_rejected(self, pH):
        with pytest.raises(ValueError, match="pH"):
            GeneratorConfig(pH=pH)

    @pytest.mark.parametrize("pH", [0.0, 7.0, 14.0])
    def test_in_range_ph_is_accepted(self, pH):
        assert GeneratorConfig(pH=pH).pH == pH


# Composition a 24-carbon skeleton bearing 3 carboxyls actually reaches
# (measured, not derived): 27 C, 6 O -> O/C 0.222, and 12 H -> H/C 0.444.
#
# Deprotonation lowers the reachable H/C: each ionized carboxyl gives up its
# -OH hydrogen, so the same skeleton at pH 7 carries 9 H, not 12 -> H/C 0.333.
# Asking for more H than the structure can hold trips the ceiling and fails
# strict validation, which neither generate_biochar nor the CLI lets a caller
# switch off -- hence two constants rather than one.
REACHABLE_OC = 0.222
REACHABLE_HC_NEUTRAL = 0.444     # 12 H: 9 edge + 3 carboxyl -OH
REACHABLE_HC_IONIZED = 0.333     # 9 H: the 3 carboxyl -OH are gone


class TestGenerateBiocharFunction:
    def test_ph_reaches_the_convenience_function(self, tmp_path):
        from biochar import generate_biochar

        res = generate_biochar(
            target_num_carbons=24,
            functional_groups={"carboxyl": 3},
            H_C_ratio=REACHABLE_HC_IONIZED,
            O_C_ratio=REACHABLE_OC,
            pH=7.0,
            seed=5,
            molecule_name="BCT",
            output_directory=str(tmp_path),
            write_files=False,
        )
        assert res.composition.net_charge < 0


class TestCLI:
    def test_ph_flag_runs_and_charges_the_structure(self, tmp_path):
        from biochar.cli import main

        rc = main([
            "--carbons", "24", "--carboxyl", "3", "--pH", "7",
            "--hc-ratio", str(REACHABLE_HC_IONIZED),
            "--oc-ratio", str(REACHABLE_OC),
            "--name", "BCT", "--output", str(tmp_path), "--seed", "5",
        ])
        assert rc == 0
        itp = list(tmp_path.glob("*.itp"))
        assert itp, "no .itp written"

    def test_out_of_range_ph_exits_cleanly_without_a_traceback(self, tmp_path, capsys):
        from biochar.cli import main

        rc = main([
            "--carbons", "24", "--pH", "99",
            "--name", "BCT", "--output", str(tmp_path),
        ])
        assert rc == 1
        err = capsys.readouterr().err
        assert "pH" in err
        assert "Traceback" not in err

    def test_omitting_ph_still_works(self, tmp_path):
        from biochar.cli import main

        rc = main([
            "--carbons", "24", "--carboxyl", "3",
            "--hc-ratio", str(REACHABLE_HC_NEUTRAL),
            "--oc-ratio", str(REACHABLE_OC),
            "--name", "BCT", "--output", str(tmp_path), "--seed", "5",
        ])
        assert rc == 0


class TestSweepAxis:
    def test_ph_is_a_valid_sweep_axis(self):
        from biochar.sweep import expand_grid

        points = expand_grid({"pH": [3.0, 7.0, 11.0]})
        assert len(points) == 3

    def test_ph_axis_produces_readable_labels(self):
        from biochar.sweep import expand_grid

        labels = [p.label for p in expand_grid({"pH": [3.0, 7.0]})]
        assert all("pH" in lbl for lbl in labels), labels

    def test_ph_sweep_produces_a_titration_series(self, tmp_path):
        from biochar.sweep import run_sweep

        summary = run_sweep(
            {
                "name": "ph_titration",
                "output_directory": str(tmp_path),
                "axes": {"pH": [2.0, 7.0]},
                "fixed": {
                    "target_num_carbons": 36,
                    "functional_groups": {"carboxyl": 3},
                    "molecule_name": "BCT",
                    "seed": 5,
                    "strict": False,
                    "H_C_tolerance": 1.0,
                    "O_C_tolerance": 1.0,
                },
            }
        )
        assert summary["n_built"] == 2
        charges = [r.net_charge for r in summary["results"]]
        assert charges[0] > charges[1], (
            f"pH 2 should be less anionic than pH 7: {charges}"
        )

    def test_net_charge_is_recorded_in_the_manifest(self, tmp_path):
        """Downstream md_setup consumers must see the charge budget."""
        from biochar.sweep import run_sweep

        summary = run_sweep(
            {
                "name": "ph_manifest",
                "output_directory": str(tmp_path),
                "axes": {"pH": [7.0]},
                "fixed": {
                    "target_num_carbons": 36,
                    "functional_groups": {"carboxyl": 3},
                    "molecule_name": "BCT",
                    "seed": 5,
                    "strict": False,
                    "H_C_tolerance": 1.0,
                    "O_C_tolerance": 1.0,
                },
            }
        )
        row = summary["results"][0].to_row()
        assert "net_charge" in row
        assert row["net_charge"] < 0


class TestChargeMethodInteraction:
    """
    Not every charge backend can represent an ion. This must fail loudly rather
    than silently return the neutral structure the caller asked not to have.
    """

    def test_ml_charges_with_ph_is_rejected(self):
        with pytest.raises(ValueError, match="ml"):
            GeneratorConfig(pH=7.0, charge_method="ml")

    def test_ml_charges_without_ph_is_still_allowed(self):
        assert GeneratorConfig(charge_method="ml").charge_method == "ml"

    def test_opls_charges_with_ph_is_allowed(self):
        assert GeneratorConfig(pH=7.0, charge_method="opls").pH == 7.0

    def test_qm_charges_with_ph_is_allowed(self):
        """qm reads the formal charge and passes it to MOPAC as CHARGE=."""
        assert GeneratorConfig(pH=7.0, charge_method="qm").pH == 7.0


class TestItpExport:
    def test_itp_records_a_nonzero_qtot_for_an_ionized_structure(self, tmp_path):
        g, _, comp = gen(functional_groups={"carboxyl": 3}, pH=7.0)
        assert comp.net_charge < 0, "fixture should be charged"
        g.export_gromacs(output_directory=str(tmp_path), basename="ion")

        itp = (tmp_path / "ion.itp").read_text()
        assert "qtot" in itp
        assert "total charge:" in itp
        assert f"formal {comp.net_charge:+d} e" in itp

    def test_itp_names_the_step_that_neutralises_a_charged_molecule(self, tmp_path):
        g, _, comp = gen(functional_groups={"carboxyl": 3}, pH=7.0)
        g.export_gromacs(output_directory=str(tmp_path), basename="ion")
        itp = (tmp_path / "ion.itp").read_text()
        assert "genion" in itp

    def test_neutral_structure_reports_zero_and_no_warning(self, tmp_path):
        g, _, comp = gen(functional_groups={"carboxyl": 3})
        assert comp.net_charge == 0
        g.export_gromacs(output_directory=str(tmp_path), basename="neu")
        itp = (tmp_path / "neu.itp").read_text()
        assert "formal +0 e" in itp
        assert "genion" not in itp, "neutral molecule should not mention genion"

    def test_itp_qtot_matches_the_formal_charge(self, tmp_path):
        import re

        g, _, comp = gen(functional_groups={"carboxyl": 3}, pH=7.0)
        g.export_gromacs(output_directory=str(tmp_path), basename="ion")
        itp = (tmp_path / "ion.itp").read_text()

        qtots = [float(m) for m in re.findall(r"; qtot ([-+][\d.]+)", itp)]
        assert qtots, "no running qtot found"
        assert qtots[-1] == pytest.approx(float(comp.net_charge), abs=1e-3)
