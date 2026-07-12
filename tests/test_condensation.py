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
    render_condensation_script,
    render_mdp_set,
    write_condensation_setup,
)


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
