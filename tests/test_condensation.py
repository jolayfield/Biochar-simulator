"""
Tests for biochar.condensation — the Wood et al. 2024 annealing setup (Phase 1).

Locks the exact protocol numbers (HTT-scaled peak temperature/timestep and the
NVT/anneal/final durations) against Wood Tables 6 & 7, and checks the run-script
+ file-writing behaviour. Setup-only; no gmx invoked.
"""

import stat

import pytest

from biochar.condensation import (
    AnnealSpec,
    CondensationError,
    anneal_spec_for_htt,
    estimate_box_nm,
    generate_and_condense,
    moleculetype_name,
    render_condensation_script,
    render_condensation_top,
    render_mdp_set,
    render_surface_script,
    setup_condensation,
    setup_surface,
    write_condensation_setup,
)


def _valid_gro() -> str:
    """A 3-atom molecule spanning 0.30 nm in x (largest extent)."""
    coords = [(0.0, 0.0, 0.0), (0.15, 0.0, 0.0), (0.30, 0.15, 0.10)]
    lines = ["test molecule", f"{len(coords):>5}"]
    for i, (x, y, z) in enumerate(coords, 1):
        lines.append(f"{1:>5}{'MOL':<5}{('C' + str(i)):>5}{i:>5}{x:8.3f}{y:8.3f}{z:8.3f}")
    lines.append("   1.000   1.000   1.000")
    return "\n".join(lines) + "\n"


_ITP = """[ moleculetype ]
; name  nrexcl
BCX     3

[ atoms ]
;  nr  type resnr res atom cgnr charge mass
    1  opls_145  1  BCX  C1  1  -0.115  12.011
"""


def _val(mdp: str, key: str) -> str:
    for ln in mdp.splitlines():
        s = ln.strip()
        if s.startswith(key) and "=" in s:
            return s.split("=", 1)[1].strip()
    raise KeyError(key)


class TestHTTScaling:
    @pytest.mark.parametrize("htt,peak,dt", [(400, 1000.0, 1.0),
                                             (600, 2000.0, 0.5),
                                             (800, 3000.0, 0.5)])
    def test_wood_anchors_exact(self, htt, peak, dt):
        s = anneal_spec_for_htt(htt)
        assert (s.peak_T_K, s.timestep_fs) == (peak, dt)

    def test_interpolates_between_anchors(self):
        assert anneal_spec_for_htt(500).peak_T_K == 1500.0  # midpoint 400-600
        assert anneal_spec_for_htt(500).timestep_fs == 0.5

    def test_clamps_outside_anchors(self):
        assert anneal_spec_for_htt(200).peak_T_K == 1000.0
        assert anneal_spec_for_htt(1200).peak_T_K == 3000.0

    def test_invalid_spec_rejected(self):
        with pytest.raises(CondensationError):
            AnnealSpec(peak_T_K=100.0, timestep_fs=1.0)
        with pytest.raises(CondensationError):
            AnnealSpec(peak_T_K=1000.0, timestep_fs=0.0)


class TestMdpNumbers:
    @pytest.mark.parametrize("htt,nvt_steps,anneal_steps,dt", [
        (400, "10000000", "25000000", "0.0010"),   # 1 fs
        (600, "20000000", "50000000", "0.0005"),   # 0.5 fs
        (800, "20000000", "50000000", "0.0005"),
    ])
    def test_durations(self, htt, nvt_steps, anneal_steps, dt):
        m = render_mdp_set(anneal_spec_for_htt(htt))
        assert _val(m["nvt.mdp"], "nsteps") == nvt_steps
        assert _val(m["nvt.mdp"], "dt") == dt
        assert _val(m["npt_anneal.mdp"], "nsteps") == anneal_steps

    def test_anneal_schedule(self):
        m = render_mdp_set(anneal_spec_for_htt(600))
        assert _val(m["npt_anneal.mdp"], "annealing-time") == "0 10000 20000 25000"
        assert _val(m["npt_anneal.mdp"], "annealing-temp") == "2000 2000 300 300"
        assert _val(m["npt_anneal.mdp"], "ref_p") == "100.0"      # Berendsen 100 bar
        assert _val(m["npt_anneal.mdp"], "pcoupltype") == "isotropic"

    def test_final_always_2fs_and_1bar(self):
        for htt in (400, 800):
            m = render_mdp_set(anneal_spec_for_htt(htt))
            assert _val(m["npt_final.mdp"], "dt") == "0.002"
            assert _val(m["npt_final.mdp"], "nsteps") == "5000000"   # 10 ns @ 2 fs
            assert _val(m["npt_final.mdp"], "ref_p") == "1.0"
            assert _val(m["npt_final.mdp"], "ref_t") == "300"

    def test_cutoffs_and_barostat(self):
        m = render_mdp_set(anneal_spec_for_htt(400))
        assert _val(m["em.mdp"], "rvdw") == "1.4"
        assert _val(m["nvt.mdp"], "rvdw") == "1.0"          # 1 nm during hot stages
        assert _val(m["npt_anneal.mdp"], "rvdw") == "1.0"
        assert _val(m["npt_final.mdp"], "rvdw") == "1.4"    # 1.4 nm for final
        assert _val(m["npt_anneal.mdp"], "tau_t") == "0.1"  # v-rescale 0.1 ps


class TestRunScript:
    def test_repeats_and_stage_order(self):
        sc = render_condensation_script("packed.gro", "system.top",
                                        anneal_spec_for_htt(400), n_repeats=3)
        assert "seq 1 3" in sc
        assert sc.index("em.mdp") < sc.index("nvt.mdp") < sc.index("npt_anneal.mdp") < sc.index("npt_final.mdp")
        assert "packed.gro" in sc and "system.top" in sc

    def test_zero_repeats_rejected(self):
        with pytest.raises(CondensationError):
            render_condensation_script("a.gro", "b.top", anneal_spec_for_htt(400), n_repeats=0)


class TestWriteSetup:
    def test_writes_all_files_and_executable_script(self, tmp_path):
        out = write_condensation_setup(tmp_path / "run", "packed.gro", "system.top", htt_c=600)
        for f in ("em.mdp", "nvt.mdp", "npt_anneal.mdp", "npt_final.mdp", "run_condensation.sh"):
            assert (out / f).exists()
        assert (out / "run_condensation.sh").stat().st_mode & stat.S_IXUSR

    def test_requires_htt_or_spec(self, tmp_path):
        with pytest.raises(CondensationError):
            write_condensation_setup(tmp_path / "run", "a.gro", "b.top")


# --------------------------------------------------------------------------- #
# Packing helpers (single molecule -> bulk of N copies)
# --------------------------------------------------------------------------- #
class TestPackingHelpers:
    def test_moleculetype_name(self):
        assert moleculetype_name(_ITP) == "BCX"

    def test_moleculetype_name_missing(self):
        with pytest.raises(CondensationError):
            moleculetype_name("[ atoms ]\n1 opls_145\n")

    def test_estimate_box_scales(self):
        # extent 0.30 -> cell max(0.30,0.5)*1.6 = 0.8; 8 copies -> ceil(2) grid -> 1.6 nm
        assert estimate_box_nm(_valid_gro(), 8) == 1.6
        # more copies -> larger box
        assert estimate_box_nm(_valid_gro(), 64) > estimate_box_nm(_valid_gro(), 8)

    def test_top_has_moltype_and_count(self):
        top = render_condensation_top("mol.itp", "BCX", 20)
        assert '#include "mol.itp"' in top
        assert '#include "oplsaa.ff/forcefield.itp"' in top
        assert "BCX" in top and " 20" in top
        assert "water" not in top.lower() and "sol" not in top.lower()  # dry


class TestSetupCondensation:
    @pytest.fixture()
    def molecule(self, tmp_path):
        (tmp_path / "mol.gro").write_text(_valid_gro())
        (tmp_path / "mol.itp").write_text(_ITP)
        return tmp_path / "mol.gro", tmp_path / "mol.itp"

    def test_setup_writes_and_packs(self, molecule, tmp_path):
        gro, itp = molecule
        out = setup_condensation(tmp_path / "run", gro, itp, n_copies=16, htt_c=600)
        for f in ("mol.gro", "mol.itp", "system.top", "em.mdp", "nvt.mdp",
                  "npt_anneal.mdp", "npt_final.mdp", "run_condensation.sh"):
            assert (out / f).exists()
        # dry .top with the right count
        assert "BCX" in (out / "system.top").read_text()
        assert " 16" in (out / "system.top").read_text()
        # per-repeat packing, seeded, before EM
        sc = (out / "run_condensation.sh").read_text()
        assert "insert-molecules" in sc and "-nmol 16" in sc and '-seed "$rep"' in sc
        assert sc.index("insert-molecules") < sc.index("em.mdp")
        # HTT-scaled anneal temperature (600 C -> 2000 K)
        assert "ref_t           = 2000" in (out / "nvt.mdp").read_text()

    def test_explicit_box_used(self, molecule, tmp_path):
        gro, itp = molecule
        out = setup_condensation(tmp_path / "run", gro, itp, n_copies=8, htt_c=400, box_nm=12.0)
        assert "-box 12 12 12" in (out / "run_condensation.sh").read_text()

    def test_missing_molecule_raises(self, tmp_path):
        with pytest.raises(CondensationError):
            setup_condensation(tmp_path / "run", tmp_path / "nope.gro",
                               tmp_path / "nope.itp", n_copies=4, htt_c=400)


class TestGenerateAndCondense:
    def test_end_to_end(self, tmp_path):
        gen = pytest.importorskip("biochar.biochar_generator")
        cfg = gen.GeneratorConfig(target_num_carbons=30, H_C_ratio=0.5, O_C_ratio=0.1,
                                  molecule_name="BC30", strict=False, seed=1)
        out = generate_and_condense(tmp_path / "run", n_copies=12,
                                    generator_config=cfg, htt_c=800)
        assert (out / "system.top").exists()
        assert (out / "molecule.gro").exists() and (out / "molecule.itp").exists()
        sc = (out / "run_condensation.sh").read_text()
        assert "-nmol 12" in sc
        assert "ref_t           = 3000" in (out / "nvt.mdp").read_text()  # 800 C -> 3000 K


# --------------------------------------------------------------------------- #
# Phase 4 — surface creation (z-expand + semi-isotropic NPT)
# --------------------------------------------------------------------------- #
class TestSurface:
    @pytest.fixture()
    def bulk(self, tmp_path):
        (tmp_path / "bulk.gro").write_text(
            "bulk\n    1\n    1MOL     C1    1   0.000   0.000   0.000\n   3.000   3.000   2.500\n"
        )
        (tmp_path / "system.top").write_text(
            '#include "oplsaa.ff/forcefield.itp"\n#include "mol.itp"\n\n[ molecules ]\nBCX 50\n'
        )
        (tmp_path / "mol.itp").write_text("[ moleculetype ]\nBCX 3\n")
        return tmp_path

    def test_writes_files_and_copies_inputs(self, bulk, tmp_path):
        out = setup_surface(tmp_path / "surf", bulk / "bulk.gro", bulk / "system.top",
                            itp=bulk / "mol.itp")
        for f in ("bulk.gro", "system.top", "mol.itp", "em.mdp", "surf_npt.mdp", "run_surface.sh"):
            assert (out / f).exists()
        assert (out / "run_surface.sh").stat().st_mode & stat.S_IXUSR

    def test_semiisotropic_z_frozen(self, bulk, tmp_path):
        out = setup_surface(tmp_path / "surf", bulk / "bulk.gro", bulk / "system.top")
        m = (out / "surf_npt.mdp").read_text()
        assert "pcoupltype      = semiisotropic" in m
        assert "compressibility = 4.5e-5 0" in m   # xy free, z frozen (keep gap)
        assert "ref_p           = 1.0 1.0" in m
        assert "nsteps          = 5000000" in m and "dt              = 0.002" in m

    def test_script_expands_z_then_em_then_npt(self, bulk, tmp_path):
        out = setup_surface(tmp_path / "surf", bulk / "bulk.gro", bulk / "system.top", gap_nm=8.0)
        sc = (out / "run_surface.sh").read_text()
        assert "editconf" in sc and "NEWZ" in sc and "GAP=8" in sc
        assert sc.index("editconf") < sc.index("em.mdp") < sc.index("surf_npt.mdp")

    def test_bad_gap_and_missing_input(self, bulk, tmp_path):
        with pytest.raises(CondensationError):
            render_surface_script("bulk.gro", "system.top", gap_nm=0)
        with pytest.raises(CondensationError):
            setup_surface(tmp_path / "surf", tmp_path / "nope.gro", bulk / "system.top")
