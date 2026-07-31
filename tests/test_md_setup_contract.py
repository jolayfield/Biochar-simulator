"""
Tests for the rqm/md-setup.md scenarios -- the contract of the run directories
this module writes, as distinct from tests/test_md_setup.py's coverage of which
files appear and which stages are present.

Everything here reads rendered artifacts: the .mdp text and the driver script.
That is deliberate. What md_setup emits is never read by Python again -- the
first consumer is grompp, hours later and usually on another machine -- so the
text itself is the contract, and no GROMACS binary is needed to check it.

Every scenario here passes. The eight that once carried xfail(strict=True) were
retired by XPASS when their fixes landed.
"""

import csv

import pytest

from biochar.export.md_setup import MDSetupConfig, setup_md_from_manifest, setup_one_structure
from biochar.workflows.condensation import anneal_spec_for_htt

GRO = """minimal test structure
    3
    1BCX     C0    1   0.000   0.000   0.000
    1BCX     C1    2   0.150   0.000   0.000
    1BCX     C2    3   0.300   0.000   0.000
   2.00000   2.00000   2.00000
"""

TOP = """#include "oplsaa.ff/forcefield.itp"
#include "{itp_name}"

[ system ]
BCX

[ molecules ]
BCX 1
"""

ITP = """[ moleculetype ]
BCX                  3
"""


def _structure(tmp_path, stem="structure", itp_stem=None):
    """Write a minimal .gro/.top/.itp trio; return their paths.

    itp_stem defaults to matching the .top's stem. Passing a different one
    reproduces the surface layout, where the topology is <base>.top and the
    include is <base>_sheet.itp.
    """
    itp_stem = itp_stem or stem
    gro = tmp_path / f"{stem}.gro"
    top = tmp_path / f"{stem}.top"
    itp = tmp_path / f"{itp_stem}.itp"
    gro.write_text(GRO)
    top.write_text(TOP.format(itp_name=itp.name))
    itp.write_text(ITP)
    return gro, top, itp


def _manifest(tmp_path, rows):
    """Write a manifest.csv with the columns md_setup reads."""
    path = tmp_path / "manifest.csv"
    fields = [
        "index", "label", "status", "gro_path", "top_path", "itp_path",
        "axis_temperature",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({f: row.get(f, "") for f in fields})
    return path


def _run_dir(tmp_path, **cfg):
    gro, top, _ = _structure(tmp_path)
    out = tmp_path / "run"
    setup_one_structure(gro, top, out, label="t", config=MDSetupConfig(**cfg))
    return out


def _script(run_dir):
    return (run_dir / "run_pipeline.sh").read_text()


def _commands(run_dir):
    """Script lines that are commands, with comments and blanks dropped.

    Matching bare substrings against the whole script is too loose: a comment
    mentioning `solvate` or `genion.tpr` reads as the command itself.
    """
    return [
        ln for ln in _script(run_dir).splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]


def _mdp(run_dir, name):
    return (run_dir / name).read_text()


def _mdp_value(text, key):
    """The right-hand side of a `key = value` line, ignoring comments."""
    for raw in text.splitlines():
        line = raw.split(";", 1)[0]
        if "=" in line and line.split("=", 1)[0].strip() == key:
            return line.split("=", 1)[1].strip()
    return None


class TestStageOrdering:
    # rq-62e0732a
    def test_solvation_topology_exists_before_solvate_runs(self, tmp_path):
        lines = _commands(_run_dir(tmp_path))

        solvate_at = next(
            (i for i, ln in enumerate(lines) if "solvate -cp" in ln), None
        )
        assert solvate_at is not None, "no solvate stage in the rendered script"

        target = "wet.top"
        created_earlier = any(
            target in ln
            and any(verb in ln for verb in ("cp ", "cat ", "> ", "touch "))
            and "wet.top.base" not in ln.split(">")[-1]
            for ln in lines[:solvate_at]
        )
        assert created_earlier, (
            f"solvate updates {target!r} but nothing before line {solvate_at} "
            f"creates it; the script runs under `set -euo pipefail`"
        )


class TestIonPlacement:
    # rq-39ce3988
    def test_each_species_is_placed_against_the_previous_result(self, tmp_path):
        # mn_calcareous_default names more than one cation.
        lines = _commands(_run_dir(tmp_path, ion_profile="mn_calcareous_default"))

        # The invocations, not the grompp lines that build genion.tpr.
        genion_at = [i for i, ln in enumerate(lines) if "genion -s" in ln]
        assert len(genion_at) > 1, "need a profile with several cations for this check"

        for first, second in zip(genion_at, genion_at[1:]):
            between = lines[first + 1 : second]
            assert any("grompp" in ln for ln in between), (
                f"genion at line {second} reuses the run input built before the "
                f"genion at line {first}, so it cannot see the ions that one placed"
            )


class TestWarningSuppression:
    DRY_STAGES = ("dry_em.mdp", "anneal_nvt.mdp", "anneal_npt.mdp", "final_npt.mdp")

    # rq-dd80065b
    def test_dry_stages_suppress_nothing(self, tmp_path):
        lines = _commands(_run_dir(tmp_path))
        offenders = [
            ln.strip()
            for ln in lines
            if "grompp" in ln
            and "-maxwarn" in ln
            and any(stage in ln for stage in self.DRY_STAGES)
        ]
        assert not offenders, (
            "dry stages pass -maxwarn with no warning to excuse: "
            + "; ".join(o[:80] for o in offenders)
        )

    # rq-af56cae9
    def test_a_suppressed_warning_is_named_where_allowed(self, tmp_path):
        lines = _script(_run_dir(tmp_path)).splitlines()
        undocumented = []  # comments are the evidence here, so keep them
        for i, ln in enumerate(lines):
            if "grompp" not in ln or "-maxwarn" not in ln:
                continue
            context = " ".join(lines[max(0, i - 3) : i + 1]).lower()
            if not any(w in context for w in ("charge", "atom name", "warning")):
                undocumented.append(ln.strip())
        assert not undocumented, (
            "-maxwarn used with no nearby comment naming the warning: "
            + "; ".join(u[:80] for u in undocumented)
        )


class TestAnnealingSchedule:
    """The Wood scaling has one implementation, condensation.anneal_spec_for_htt.

    These read the expected values from it rather than restating 3000 K, so the
    test cannot drift away from the scaling it is checking.
    """

    # rq-5b7a5ab8
    def test_peak_temperature_follows_the_pyrolysis_temperature(self, tmp_path):
        gro, top, itp = _structure(tmp_path)
        manifest = _manifest(tmp_path, [{
            "index": "0", "label": "hot", "status": "strict_pass",
            "gro_path": str(gro), "top_path": str(top), "itp_path": str(itp),
            "axis_temperature": "800",
        }])
        setup_md_from_manifest(manifest, tmp_path / "runs")

        expected = anneal_spec_for_htt(800.0).peak_T_K
        temps = _mdp_value(_mdp(tmp_path / "runs" / "hot", "anneal_npt.mdp"), "annealing-temp")
        assert temps, "no annealing-temp line"
        peak = max(float(t) for t in temps.split())
        assert peak == pytest.approx(expected), (
            f"annealing peak is {peak} K; the Wood scaling for 800 C gives {expected} K"
        )

    # rq-997c1640
    def test_timestep_drops_with_the_schedule(self, tmp_path):
        gro, top, itp = _structure(tmp_path)
        manifest = _manifest(tmp_path, [{
            "index": "0", "label": "hot", "status": "strict_pass",
            "gro_path": str(gro), "top_path": str(top), "itp_path": str(itp),
            "axis_temperature": "800",
        }])
        setup_md_from_manifest(manifest, tmp_path / "runs")

        spec = anneal_spec_for_htt(800.0)
        assert spec.peak_T_K > 1500, "fixture should sit above the timestep threshold"
        dt = _mdp_value(_mdp(tmp_path / "runs" / "hot", "anneal_npt.mdp"), "dt")
        assert float(dt) == pytest.approx(spec.timestep_fs / 1000.0), (
            f"dt is {dt} ps; the schedule for this peak wants "
            f"{spec.timestep_fs / 1000.0} ps"
        )


class TestWetPressureCoupling:
    # rq-51c8f769
    def test_wet_stage_holds_xy_and_relaxes_z(self, tmp_path):
        text = _mdp(_run_dir(tmp_path), "wet_npt.mdp")

        assert _mdp_value(text, "pcoupltype") == "semiisotropic"
        xy, z = _mdp_value(text, "compressibility").split()
        # GROMACS reads the pair as "xy z": a zero holds that axis.
        assert float(xy) == 0.0, (
            "xy is compressible -- the barostat can squeeze the carbon lattice"
        )
        assert float(z) > 0.0, "z is held -- the solvent cannot reach its own density"


class TestRunDirectoryProvenance:
    # rq-9e4f2d8e
    def test_fallback_status_is_recorded_in_the_run_directory(self, tmp_path):
        gro, top, itp = _structure(tmp_path)
        manifest = _manifest(tmp_path, [{
            "index": "0", "label": "fb", "status": "fallback",
            "gro_path": str(gro), "top_path": str(top), "itp_path": str(itp),
        }])
        setup_md_from_manifest(manifest, tmp_path / "runs")

        run_dir = tmp_path / "runs" / "fb"
        recorded = any(
            "fallback" in p.read_text()
            for p in run_dir.iterdir()
            if p.is_file() and p.suffix in {".sh", ".json", ".txt", ".md", ""}
        )
        assert recorded, (
            f"nothing in {run_dir.name} records that this structure only "
            f"passed as a fallback"
        )

    # rq-da9561c4
    def test_a_row_with_no_files_gets_no_run_directory(self, tmp_path):
        manifest = _manifest(tmp_path, [
            {"index": "0", "label": "dead", "status": "failed"},
        ])
        results = setup_md_from_manifest(manifest, tmp_path / "runs")

        assert len(results) == 1
        assert results[0]["run_dir"] is None
        assert results[0]["skipped_reason"]
        assert not (tmp_path / "runs" / "dead").exists()


class TestIncludeFiles:
    # rq-2fdbb62b
    def test_include_is_resolved_from_the_manifest(self, tmp_path):
        # Surface layout: topology surf.top includes surf_sheet.itp.
        gro, top, itp = _structure(tmp_path, stem="surf", itp_stem="surf_sheet")
        assert itp.stem != top.stem

        manifest = _manifest(tmp_path, [{
            "index": "0", "label": "surface", "status": "strict_pass",
            "gro_path": str(gro), "top_path": str(top), "itp_path": str(itp),
        }])
        setup_md_from_manifest(manifest, tmp_path / "runs")

        copied = tmp_path / "runs" / "surface" / itp.name
        assert copied.exists(), (
            f"{itp.name} is #included by the topology but was not copied into the "
            f"run directory; grompp will not find the moleculetype"
        )
