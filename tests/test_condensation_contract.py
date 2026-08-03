"""
Tests for the rqm/condensation-annealing.md scenarios that no existing test
defends.

The published-protocol half of that document -- the anchors, the schedule, the
durations, the stage order -- is back-annotated onto tests/test_condensation.py,
which already checks it against Wood's tables. What is left here is the part
about a script that nobody watches: whether the packing step notices it placed
fewer molecules than the topology declares, whether a grompp warning it absorbs
was one anybody foresaw, and whether the directory says what it is a simulation
of.

Four scenarios carry xfail(strict=True); each names the defect it defers.
"""

import inspect
import json
import re
import subprocess
from pathlib import Path

import pytest

from biochar.workflows import condensation as cond
from biochar.workflows.condensation import (
    anneal_spec_for_htt,
    render_condensation_script,
    setup_condensation,
)

_ANCHOR_HTT = 400.0


def _script(**kw):
    spec = anneal_spec_for_htt(kw.pop("htt_c", 600.0))
    return render_condensation_script("packed.gro", "system.top", spec, **kw)


def _grompp_stages(script):
    """Each grompp invocation, paired with the comment line above it."""
    stages, comment = [], ""
    for line in script.splitlines():
        stripped = line.strip()
        if stripped.startswith("#"):
            comment = stripped
        elif "grompp" in stripped:
            stages.append((comment, stripped))
    return stages


# --------------------------------------------------------------------------- #
# Timestep rule
# --------------------------------------------------------------------------- #
class TestTheTimestepRule:
    # rq-0a4c0f30
    def test_the_timestep_halves_above_the_lowest_anchor(self):
        assert anneal_spec_for_htt(_ANCHOR_HTT).timestep_fs == 1.0
        assert anneal_spec_for_htt(_ANCHOR_HTT + 1).timestep_fs == 0.5

    # rq-fa9763e9
    @pytest.mark.xfail(
        strict=True,
        reason="the recorded rationale says the timestep drops once the peak "
               "exceeds ~1500 K, but the code compares the heat-treatment "
               "temperature with the 400 C anchor -- at 401 C the peak is about "
               "1005 K and the timestep has already halved",
    )
    def test_the_recorded_rationale_is_the_rule_the_code_follows(self):
        source = inspect.getsource(cond)
        explanation = source.split("@dataclass(frozen=True)", 1)[0]

        claims_a_peak_threshold = re.search(r"1500\s*K", explanation)
        tests_a_peak_threshold = re.search(r"peak_T\s*[<>]", source)

        assert not (claims_a_peak_threshold and not tests_a_peak_threshold), (
            "the explanation describes a peak-temperature threshold the code "
            "does not test; it keys off the lowest anchor instead"
        )
        # And the rule it does follow is stated where it is implemented.
        rule = inspect.getsource(anneal_spec_for_htt)
        assert "anchor" in rule.lower(), rule


# --------------------------------------------------------------------------- #
# Repeats
# --------------------------------------------------------------------------- #
class TestRepeats:
    # rq-516fd1c6
    def test_repeats_are_seeded_differently(self):
        script = _script(
            n_repeats=3,
            pack=cond.PackSpec(molecule_gro="m.gro", n_copies=8, box_nm=5.0),
        )
        insert = [ln for ln in script.splitlines() if "insert-molecules" in ln]
        assert insert, "the script does not pack"
        assert all("-seed" in ln for ln in insert), insert
        # One loop body, seeded from the loop variable rather than a constant.
        assert all(re.search(r'-seed\s+"?\$', ln) for ln in insert), (
            f"the insertion seed does not vary with the repeat: {insert}"
        )


# --------------------------------------------------------------------------- #
# Packing
# --------------------------------------------------------------------------- #
class TestThePackedBoxMatchesTheTopology:
    # rq-9a027a37
    @pytest.mark.xfail(
        strict=True,
        reason="gmx insert-molecules places what it can and exits successfully, "
               "so a box too tight silently yields fewer molecules than the .top "
               "declares; the first thing to notice is grompp, several commands "
               "later, reporting a coordinate count that does not match",
    )
    def test_the_script_verifies_the_packed_count(self):
        script = _script(
            pack=cond.PackSpec(molecule_gro="m.gro", n_copies=64, box_nm=6.0)
        )
        lines = script.splitlines()
        pack_at = next(i for i, ln in enumerate(lines) if "insert-molecules" in ln)
        following = "\n".join(lines[pack_at : pack_at + 12])

        assert "64" in following, (
            "nothing after the packing stage refers to the number of molecules "
            f"the topology declares:\n{following}"
        )
        assert re.search(r"\bexit\b|\breturn\b|\bdie\b", following), (
            f"the packing stage does not stop on a shortfall:\n{following}"
        )


# --------------------------------------------------------------------------- #
# grompp warnings
# --------------------------------------------------------------------------- #
class TestASuppressedWarningIsNamed:
    # rq-30534b7a
    @pytest.mark.xfail(
        strict=True,
        reason="every one of the four stages passes -maxwarn 2 though none of "
               "them expects a warning, so the allowance excuses nothing "
               "foreseen and absorbs anything that turns up",
    )
    def test_a_stage_with_no_expected_warning_suppresses_nothing(self):
        offenders = [
            (comment, stage)
            for comment, stage in _grompp_stages(_script())
            if "-maxwarn" in stage and "warn" not in comment.lower()
        ]
        assert not offenders, (
            "stages pass -maxwarn with no comment naming what it absorbs: "
            + "; ".join(stage for _, stage in offenders)
        )

    # rq-8410a38d
    @pytest.mark.xfail(
        strict=True,
        reason="the same defect from the other side: no stage names a warning, "
               "so no count can equal the number named",
    )
    def test_a_suppressed_warning_is_named_where_it_is_allowed(self):
        for comment, stage in _grompp_stages(_script()):
            match = re.search(r"-maxwarn\s+(\d+)", stage)
            if not match:
                continue
            named = re.findall(r"warning", comment, re.I)
            assert named, f"-maxwarn with no comment naming it: {stage}"
            assert int(match.group(1)) == len(named), (
                f"count {match.group(1)} does not match the {len(named)} "
                f"warning(s) named in: {comment}"
            )

    def test_the_stage_list_is_not_empty(self):
        """Guards both tests above: they pass vacuously over no stages."""
        assert len(_grompp_stages(_script())) >= 4


# --------------------------------------------------------------------------- #
# Provenance and setup-only
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def run_dir(tmp_path_factory):
    src = tmp_path_factory.mktemp("mol")
    (src / "m.gro").write_text(
        "one molecule\n    2\n"
        "    1BC      C0    1   0.000   0.000   0.000\n"
        "    1BC      C1    2   0.140   0.000   0.000\n"
        "   1.00000   1.00000   1.00000\n"
    )
    (src / "m.itp").write_text("[ moleculetype ]\nBC   3\n\n[ atoms ]\n")

    out = tmp_path_factory.mktemp("run")
    return setup_condensation(
        out, src / "m.gro", src / "m.itp", n_copies=12, htt_c=600.0
    )


class TestARunDirectoryRecordsItsSetup:
    # rq-75419e1a
    @pytest.mark.xfail(
        strict=True,
        reason="nothing records the heat-treatment temperature that produced "
               "the schedule; the .mdp carries the peak but not where it came "
               "from, and the .top a copy count but not the box it was packed "
               "into",
    )
    def test_the_run_directory_records_its_provenance(self, run_dir):
        record = Path(run_dir) / "condensation_provenance.json"
        assert record.exists(), (
            f"no provenance record; directory holds "
            f"{sorted(p.name for p in Path(run_dir).iterdir())}"
        )

        data = json.loads(record.read_text())
        assert data["htt_c"] == 600.0
        assert data["peak_T_K"] == anneal_spec_for_htt(600.0).peak_T_K
        assert data["timestep_fs"] == anneal_spec_for_htt(600.0).timestep_fs
        assert data["n_copies"] == 12
        assert data["biochar_version"]


class TestSetupOnly:
    # rq-ec21c31b
    def test_setting_up_invokes_no_external_program(self, tmp_path, monkeypatch):
        def refuse(*args, **kwargs):
            raise AssertionError(f"setup started an external process: {args}")

        for name in ("run", "check_call", "check_output", "Popen", "call"):
            monkeypatch.setattr(subprocess, name, refuse)

        src = tmp_path / "src"
        src.mkdir()
        (src / "m.gro").write_text(
            "one molecule\n    1\n"
            "    1BC      C0    1   0.000   0.000   0.000\n"
            "   1.00000   1.00000   1.00000\n"
        )
        (src / "m.itp").write_text("[ moleculetype ]\nBC   3\n")

        setup_condensation(
            tmp_path / "out", src / "m.gro", src / "m.itp",
            n_copies=4, htt_c=400.0,
        )
