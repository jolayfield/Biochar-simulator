"""
Tests for biochar/md_setup.py — GROMACS run-directory generation from a
sweep manifest (dry-anneal + solvate + ion + wet-equilibrate pipeline).

No `gmx` binary is invoked anywhere here: these tests only check that the
correct files, paths, and shell-script commands are produced.
"""

import csv
import shutil
import stat
from pathlib import Path

import pytest

from biochar.md_setup import (
    ION_PROFILES,
    IonProfile,
    MDSetupConfig,
    MDSetupError,
    MoleculeInsertion,
    PreSolvationStage,
    get_ion_profile,
    setup_md_from_manifest,
    setup_one_structure,
)
from biochar.sweep import load_sweep_config, run_sweep


# --------------------------------------------------------------------------- #
# Ion profiles
# --------------------------------------------------------------------------- #
class TestIonProfiles:
    def test_builtin_profiles_present(self):
        assert "mn_calcareous_default" in ION_PROFILES
        assert "pure_water" in ION_PROFILES

    def test_get_ion_profile_by_name(self):
        profile = get_ion_profile("mn_calcareous_default")
        assert isinstance(profile, IonProfile)
        assert profile.ca_mM > 0

    def test_get_ion_profile_passthrough(self):
        custom = IonProfile(name="custom", na_mM=10.0)
        assert get_ion_profile(custom) is custom

    def test_unknown_profile_raises(self):
        with pytest.raises(MDSetupError):
            get_ion_profile("not_a_real_profile")

    def test_pure_water_has_no_ions(self):
        profile = ION_PROFILES["pure_water"]
        assert profile.ca_mM == profile.mg_mM == profile.na_mM == profile.k_mM == 0.0


# --------------------------------------------------------------------------- #
# setup_one_structure
# --------------------------------------------------------------------------- #
@pytest.fixture()
def built_structure(tmp_path):
    """Build one real biochar structure+GROMACS export to feed into md_setup."""
    from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

    gen = BiocharGenerator(GeneratorConfig(
        target_num_carbons=30, H_C_ratio=0.5, O_C_ratio=0.1, strict=False, seed=7,
    ))
    result = gen.generate()
    gro, top, itp = gen.export_gromacs(str(tmp_path / "src"), basename="mini")
    return Path(gro), Path(top)


class TestSetupOneStructure:
    def test_writes_all_expected_files(self, built_structure, tmp_path):
        gro, top = built_structure
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini")
        written = {p.name for p in out.iterdir()}
        for expected in (
            "mini.gro", "mini.top", "mini.itp",
            "dry_em.mdp", "anneal_nvt.mdp", "anneal_npt.mdp", "final_npt.mdp",
            "wet_em.mdp", "wet_nvt.mdp", "wet_npt.mdp",
            "run_pipeline.sh",
        ):
            assert expected in written, f"missing {expected}"

    def test_pipeline_script_is_executable(self, built_structure, tmp_path):
        gro, top = built_structure
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini")
        script = out / "run_pipeline.sh"
        mode = script.stat().st_mode
        assert mode & stat.S_IXUSR

    def test_pipeline_script_stages_in_order(self, built_structure, tmp_path):
        gro, top = built_structure
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini")
        text = (out / "run_pipeline.sh").read_text()
        for marker in (
            "dry_em.mdp", "anneal_nvt.mdp", "anneal_npt.mdp", "final_npt.mdp",
            "editconf", "solvate", "genion", "wet_em.mdp", "wet_nvt.mdp", "wet_npt.mdp",
        ):
            assert marker in text
        # dry stages must precede solvation, which must precede wet stages
        # (search for the actual gmx invocation, not the header comment
        # that also mentions "solvate" in prose)
        solvate_cmd_idx = text.index('"$GMX" solvate')
        assert text.index("final_npt.mdp") < solvate_cmd_idx
        assert solvate_cmd_idx < text.rindex("wet_em.mdp")

    def test_missing_gro_raises(self, tmp_path):
        with pytest.raises(MDSetupError):
            setup_one_structure(tmp_path / "nope.gro", tmp_path / "nope.top", tmp_path / "run")

    def test_ion_valences_correct(self, built_structure, tmp_path):
        gro, top = built_structure
        cfg = MDSetupConfig(ion_profile="mn_calcareous_default")
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini", config=cfg)
        text = (out / "run_pipeline.sh").read_text()
        assert "-pname CA -pq 2" in text
        assert "-pname MG -pq 2" in text
        assert "-pname NA -pq 1" in text
        assert "-pname K -pq 1" in text

    def test_pure_water_profile_has_default_genion(self, built_structure, tmp_path):
        gro, top = built_structure
        cfg = MDSetupConfig(ion_profile="pure_water")
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini", config=cfg)
        text = (out / "run_pipeline.sh").read_text()
        assert "genion" in text
        assert "-neutral" in text

    def test_custom_solvent_pad_and_water_model(self, built_structure, tmp_path):
        gro, top = built_structure
        cfg = MDSetupConfig(solvent_pad_nm=2.0, water_model="tip3p")
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini", config=cfg)
        text = (out / "run_pipeline.sh").read_text()
        assert "-d 2.0" in text
        assert "tip3p.gro" in text


# --------------------------------------------------------------------------- #
# PreSolvationStage — the generic pre-solvation molecule-insertion seam
# --------------------------------------------------------------------------- #
class TestPreSolvationStage:
    def _stage(self, tmp_path):
        # Fabricate the files a downstream workflow would hand in.
        (tmp_path / "MOL.gro").write_text("dummy ligand gro\n")
        (tmp_path / "merged.top").write_text("dummy merged topology\n")
        return PreSolvationStage(
            name="Insert test molecules",
            insertions=[MoleculeInsertion("MOL.gro", n_copies=4, n_try=250)],
            solvation_top="merged.top",
            extra_files=[str(tmp_path / "MOL.gro"), str(tmp_path / "merged.top")],
        )

    def test_extra_files_copied_and_stage_rendered(self, built_structure, tmp_path):
        gro, top = built_structure
        stage = self._stage(tmp_path)
        cfg = MDSetupConfig(pre_solvation_stage=stage)
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini", config=cfg)

        # extra_files staged into the run dir
        assert (out / "MOL.gro").exists()
        assert (out / "merged.top").exists()

        text = (out / "run_pipeline.sh").read_text()
        # insertion stage rendered with the declared counts, before solvation
        assert "Insert test molecules" in text
        assert 'insert-molecules -f "$SIM/final_npt/npt.gro" -ci "$SIM/MOL.gro" -nmol 4 -try 250' in text
        ins_idx = text.index("insert-molecules")
        assert ins_idx < text.index('"$GMX" solvate')
        # solvation onward uses the stage's topology, not the bare biochar one
        assert 'cp "$SIM/merged.top" "$SIM/wet.top.base"' in text

    def test_no_stage_uses_bare_topology(self, built_structure, tmp_path):
        gro, top = built_structure
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini")
        text = (out / "run_pipeline.sh").read_text()
        assert "insert-molecules" not in text
        assert f'cp "$SIM/{top.name}" "$SIM/wet.top.base"' in text

    def test_missing_extra_file_raises(self, built_structure, tmp_path):
        gro, top = built_structure
        stage = PreSolvationStage(
            name="x", insertions=[MoleculeInsertion("MOL.gro", 1)],
            solvation_top="merged.top", extra_files=[str(tmp_path / "does_not_exist.gro")],
        )
        with pytest.raises(MDSetupError):
            setup_one_structure(gro, top, tmp_path / "run",
                                config=MDSetupConfig(pre_solvation_stage=stage))

    def test_mu3c_chain_threads_stage(self, built_structure, tmp_path):
        gro, top = built_structure
        stage = self._stage(tmp_path)
        cfg = MDSetupConfig(pre_solvation_stage=stage, cluster="mu3c")
        out = setup_one_structure(gro, top, tmp_path / "run", label="mini", config=cfg)
        assert (out / "solvate_ions.slurm").exists()
        assert (out / "submit_chain.sh").exists()
        slurm = (out / "solvate_ions.slurm").read_text()
        assert "insert-molecules" in slurm and "MOL.gro" in slurm
        chain = (out / "submit_chain.sh").read_text()
        assert "merged.top" in chain  # solvate_ions job handed the stage topology


# --------------------------------------------------------------------------- #
# setup_md_from_manifest — end to end from a real sweep
# --------------------------------------------------------------------------- #
@pytest.fixture()
def small_sweep_manifest(tmp_path):
    cfg = {
        "name": "mini_grid",
        "output_directory": str(tmp_path / "sweep"),
        "axes": {"temperature": [400, 600]},
        "fixed": {"target_num_carbons": 30, "feedstock": "softwood"},
        "seed": 3,
        "max_retries": 2,
        "on_validation_fail": "fallback",
    }
    summary = run_sweep(cfg, quiet=True)
    return Path(summary["manifest_csv"])


class TestSetupMdFromManifest:
    def test_writes_one_dir_per_row(self, small_sweep_manifest, tmp_path):
        results = setup_md_from_manifest(small_sweep_manifest, tmp_path / "md_runs")
        assert len(results) == 2
        assert all(r["run_dir"] is not None for r in results)
        for r in results:
            assert Path(r["run_dir"]).exists()
            assert (Path(r["run_dir"]) / "run_pipeline.sh").exists()

    def test_missing_manifest_raises(self, tmp_path):
        with pytest.raises(MDSetupError):
            setup_md_from_manifest(tmp_path / "no_manifest.csv", tmp_path / "out")

    def test_skips_rows_without_status_files(self, small_sweep_manifest, tmp_path):
        # Rewrite one row's status to "skipped" and confirm it's reported, not written.
        rows = list(csv.DictReader(small_sweep_manifest.open()))
        rows[0]["status"] = "skipped"
        fieldnames = rows[0].keys()
        edited = tmp_path / "edited_manifest.csv"
        with edited.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
        # gro/top paths in the manifest are relative to the sweep's own cwd,
        # not to `edited`'s directory -- copy inputs alongside so resolution
        # via candidate_bases still finds them for the row we didn't touch.
        results = setup_md_from_manifest(edited, tmp_path / "md_runs2")
        assert results[0]["run_dir"] is None
        assert "status=skipped" in results[0]["skipped_reason"]


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
class TestMdSetupCli:
    def test_cli_runs_end_to_end(self, small_sweep_manifest, tmp_path):
        from biochar.md_setup_cli import main

        out_root = tmp_path / "cli_runs"
        rc = main([str(small_sweep_manifest), "--output-root", str(out_root)])
        assert rc == 0
        assert out_root.exists()
        assert len(list(out_root.iterdir())) == 2

    def test_cli_rejects_bad_ion_profile(self, small_sweep_manifest, tmp_path):
        from biochar.md_setup_cli import main

        with pytest.raises(SystemExit):
            main([str(small_sweep_manifest), "--output-root", str(tmp_path / "x"),
                  "--ion-profile", "not_real"])
