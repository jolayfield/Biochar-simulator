"""
Tests for the rqm/parameter-sweep.md scenarios -- what the manifest claims
about the structures beside it.

Distinct from tests/test_sweep.py, which checks that grid expansion, building
and manifest writing work. These check that what is written down is true: that
one row means one structure, that the name in the row is the name in the
topology, that a status means one thing, and that a retry budget is spent on
failures a retry could actually change.

Most of these stub `_build_once` rather than generating structures. The subject
is the driver's bookkeeping around a build, not the build, and a stub makes the
failing seeds exact instead of hoped-for -- the alternative is a request chosen
because it usually fails, which is a slower test that proves less.

Every scenario here passes. The eight that carried xfail(strict=True) were
retired by XPASS when their fixes landed.
"""

import json
import re
from pathlib import Path

import pytest

import biochar
from biochar.pipeline.biochar_generator import ValidationError
from biochar.workflows import sweep as sweep_mod
from biochar.workflows.sweep import (
    SweepError,
    build_point,
    expand_grid,
    run_sweep,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- #
# Stubs
# --------------------------------------------------------------------------- #
class _FakeGen:
    """Just enough of BiocharGenerator for _ok_result to read."""

    def __init__(self, errors=()):
        self.validation_report = (not errors, list(errors), [], {})
        self.ring_composition = {"hexagon": 4}


class _FakeComp:
    molecular_formula = "C16H10"
    molecular_weight = 202.25
    num_carbons = 16
    H_C_ratio = 0.625
    O_C_ratio = 0.0
    functional_groups = {}
    net_charge = 0
    ionized_counts = {}


SHORTFALL = "functional groups were not fully placed (carboxyl 12/500)"


def _install(monkeypatch, behaviour):
    """Replace _build_once with `behaviour(seed, strict) -> paths tuple`."""

    def fake(config_kwargs, seed, strict, out_dir, basename,
             write_files=True, quiet=True):
        return behaviour(seed, strict)

    monkeypatch.setattr(sweep_mod, "_build_once", fake)


def _built(errors=()):
    return _FakeGen(errors), _FakeComp(), ("a.gro", "a.top", "a.itp")


def _strict_never_passes(seed, strict):
    """Fails strict at every seed with an unchanging report; builds loosely."""
    if strict:
        raise ValidationError(f"Strict mode: {SHORTFALL}")
    return _built([SHORTFALL])


def _point(**fixed):
    return expand_grid(axes={"temperature": [500]}, fixed=fixed or None)[0]


# --------------------------------------------------------------------------- #
# Naming
# --------------------------------------------------------------------------- #
class TestEveryPointIsDistinguishable:
    # rq-2a033e95
    def test_a_grid_whose_names_collide_is_refused(self):
        # BC0100 and BC1000 both truncate to BC100.
        axes = {"target_num_carbons": list(range(100, 1101))}
        with pytest.raises(SweepError) as exc:
            expand_grid(axes=axes)
        assert "BC100" in str(exc.value), (
            f"the refusal should name the colliding molecule name: {exc.value}"
        )

    def test_a_grid_below_the_collision_point_expands(self):
        """The check must not refuse grids that are actually fine."""
        pts = expand_grid(axes={"target_num_carbons": list(range(100, 1100))})
        names = [p.molecule_name for p in pts]
        assert len(set(names)) == len(names) == 1000


class TestTheRecordedNameIsTheWrittenName:
    # rq-4156245d
    def test_a_molecule_name_fixed_across_the_sweep_is_refused(self):
        with pytest.raises(SweepError) as exc:
            expand_grid(axes={"temperature": [400, 500]},
                        fixed={"molecule_name": "FIXED"})
        assert "molecule_name" in str(exc.value)

    # rq-a118d6d2
    @pytest.mark.parametrize("template", ["BC{i:03d}", "M{i}", "X{temperature}"])
    def test_the_reported_name_is_the_configured_name(self, template):
        pts = expand_grid(
            axes={"temperature": [400, 500, 600]},
            fixed={"target_num_carbons": 60},
            name_template=template,
        )
        mismatched = [
            (p.index, p.molecule_name, p.config_kwargs["molecule_name"])
            for p in pts
            if p.molecule_name != p.config_kwargs["molecule_name"]
        ]
        assert not mismatched, (
            f"points report a name their generator is not configured with: "
            f"{mismatched}"
        )


# --------------------------------------------------------------------------- #
# Seed retry
# --------------------------------------------------------------------------- #
class TestSeedRetryIsAVarianceRemedy:
    # rq-bedc92d0
    def test_a_validation_failure_is_retried_at_a_new_seed(self, monkeypatch, tmp_path):
        def behaviour(seed, strict):
            if strict and seed < 2:
                raise ValidationError(f"seed {seed} missed the O/C target")
            return _built()

        _install(monkeypatch, behaviour)
        res = build_point(_point(), tmp_path, base_seed=0, max_retries=8)

        assert res.status == "strict_pass"
        assert res.seed_used == 2, "the passing seed is not the one recorded"
        assert res.n_attempts == 3

    # rq-d600f3fc
    def test_a_failure_that_is_not_a_validation_failure_is_not_retried(
        self, monkeypatch, tmp_path
    ):
        def behaviour(seed, strict):
            raise RuntimeError("the skeleton could not be grown")

        _install(monkeypatch, behaviour)
        res = build_point(_point(), tmp_path, base_seed=0, max_retries=8)

        assert res.status == "failed"
        assert res.n_attempts == 1, (
            f"a deterministic failure was retried {res.n_attempts} times"
        )
        assert "skeleton could not be grown" in (res.error or "")

    # rq-a47e48e0
    def test_a_validation_failure_that_repeats_unchanged_stops_the_loop(
        self, monkeypatch, tmp_path
    ):
        _install(monkeypatch, _strict_never_passes)
        res = build_point(_point(), tmp_path, base_seed=0, max_retries=8,
                          on_validation_fail="skip")

        assert res.n_attempts <= 2, (
            f"the same report was produced {res.n_attempts} times; two identical "
            f"reports are enough to know a third seed will not differ"
        )


# --------------------------------------------------------------------------- #
# Failure modes
# --------------------------------------------------------------------------- #
class TestEachFailureModeIsADifferentAnswer:
    # rq-0cf2e5d4
    def test_the_default_mode_keeps_the_rejected_structure(self, monkeypatch, tmp_path):
        _install(monkeypatch, _strict_never_passes)
        res = build_point(_point(), tmp_path, base_seed=0, max_retries=2,
                          on_validation_fail="fallback")

        assert res.status == "fallback"
        assert res.gro_path, "a fallback point wrote no structure"
        assert any(SHORTFALL in e for e in res.validation_errors), (
            f"the fallback row does not carry why it was a fallback: "
            f"{res.validation_errors}"
        )

    # rq-086427f6
    def test_skip_records_the_point_as_skipped(self, monkeypatch, tmp_path):
        _install(monkeypatch, _strict_never_passes)
        res = build_point(_point(), tmp_path, base_seed=0, max_retries=2,
                          on_validation_fail="skip")

        assert res.status == "skipped"
        assert res.gro_path is None
        assert any(SHORTFALL in e for e in res.validation_errors)

    # rq-626e9fb8
    def test_strict_refuses_the_sweep(self, monkeypatch, tmp_path):
        _install(monkeypatch, _strict_never_passes)
        point = _point()

        with pytest.raises(SweepError) as exc:
            build_point(point, tmp_path, base_seed=0, max_retries=2,
                        on_validation_fail="strict")

        message = str(exc.value)
        assert point.label in message, f"the refusal does not name the point: {message}"
        assert SHORTFALL in message, f"the refusal does not say what was wrong: {message}"


# --------------------------------------------------------------------------- #
# Manifest
# --------------------------------------------------------------------------- #
class TestTheManifestRecordsWhatProducedIt:
    # rq-43924feb
    def test_the_manifest_names_the_package_version(self, monkeypatch, tmp_path):
        _install(monkeypatch, _strict_never_passes)
        summary = run_sweep({
            "name": "versioned",
            "output_directory": str(tmp_path / "out"),
            "axes": {"temperature": [400]},
            "max_retries": 1,
        })

        payload = json.loads(Path(summary["manifest_json"]).read_text())
        recorded = payload["meta"].get("biochar_version")
        assert recorded == biochar.__version__, (
            f"manifest meta records {recorded!r}, package is "
            f"{biochar.__version__!r}; meta keys were {sorted(payload['meta'])}"
        )


# --------------------------------------------------------------------------- #
# Status vocabulary
# --------------------------------------------------------------------------- #
def _readme_statuses():
    text = (REPO_ROOT / "README.md").read_text()
    section = re.search(
        r"### Validation and the `status` column(.*?)\n---\n", text, re.S
    )
    assert section, "the README's sweep status section moved; update this test"
    body = section.group(1).split("|---|---|", 1)
    assert len(body) == 2, "the README's status table is no longer a table"
    return set(re.findall(r"^\|\s*`(\w+)`\s*\|", body[1], re.M))


class TestAStatusMeansOneThing:
    @staticmethod
    def _observed(monkeypatch, tmp_path):
        """Every status build_point can assign, driven out of it directly."""
        seen = set()

        _install(monkeypatch, lambda seed, strict: _built())
        seen.add(build_point(_point(), tmp_path, max_retries=1).status)

        _install(monkeypatch, _strict_never_passes)
        for mode in ("fallback", "skip"):
            seen.add(build_point(_point(), tmp_path, max_retries=1,
                                 on_validation_fail=mode).status)

        def crash(seed, strict):
            raise RuntimeError("boom")

        _install(monkeypatch, crash)
        seen.add(build_point(_point(), tmp_path, max_retries=1).status)
        return seen

    # rq-2da595b7
    def test_every_status_a_build_assigns_is_declared(self, monkeypatch, tmp_path):
        declared = set(sweep_mod.POINT_STATUSES)
        observed = self._observed(monkeypatch, tmp_path)
        assert observed == declared, (
            f"statuses observed {sorted(observed)} but declared {sorted(declared)}"
        )

    # rq-7f22837b
    def test_the_readme_table_names_the_declared_statuses(self):
        assert _readme_statuses() == set(sweep_mod.POINT_STATUSES)
