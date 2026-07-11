"""
Tests for biochar/sweep.py — parameter-sweep driver and manifest.

Covers grid expansion, per-point build with seed-retry/fallback, config
validation, manifest writing, and the sweep-CLI entry point.
"""

import csv
import json
from pathlib import Path

import pytest

from biochar.sweep import (
    SweepError,
    GridPoint,
    PointResult,
    expand_grid,
    build_point,
    run_sweep,
    load_sweep_config,
)


# --------------------------------------------------------------------------- #
# Grid expansion
# --------------------------------------------------------------------------- #
class TestExpandGrid:
    def test_cartesian_product_size(self):
        pts = expand_grid(axes={"temperature": [300, 500], "feedstock": ["softwood", "hardwood"]})
        assert len(pts) == 4
        assert all(isinstance(p, GridPoint) for p in pts)

    def test_single_axis(self):
        pts = expand_grid(axes={"temperature": [300, 400, 500]})
        assert len(pts) == 3
        assert [p.axis_values["temperature"] for p in pts] == [300, 400, 500]

    def test_fixed_applied_to_all(self):
        pts = expand_grid(
            axes={"temperature": [300, 500]},
            fixed={"target_num_carbons": 80},
        )
        assert all(p.config_kwargs["target_num_carbons"] == 80 for p in pts)

    def test_axis_overrides_fixed(self):
        pts = expand_grid(
            axes={"target_num_carbons": [60, 100]},
            fixed={"target_num_carbons": 80},
        )
        assert sorted(p.config_kwargs["target_num_carbons"] for p in pts) == [60, 100]

    def test_labels_unique_and_safe(self):
        pts = expand_grid(axes={"temperature": [300, 500], "feedstock": ["softwood", "hardwood"]})
        labels = [p.label for p in pts]
        assert len(set(labels)) == len(labels)
        assert all("/" not in lbl and " " not in lbl for lbl in labels)
        assert pts[0].label.startswith("T300")

    def test_molecule_name_truncated_to_5(self):
        pts = expand_grid(axes={"temperature": [300]}, name_template="BIOCHAR{i}")
        assert len(pts[0].molecule_name) <= 5

    def test_functional_group_axis_labeled(self):
        pts = expand_grid(axes={"functional_groups": [{"phenolic": 4}, {"carboxyl": 8}]})
        assert len(pts) == 2
        assert all(p.label for p in pts)

    def test_empty_axes_rejected(self):
        with pytest.raises(SweepError):
            expand_grid(axes={})

    def test_unknown_field_rejected(self):
        with pytest.raises(SweepError):
            expand_grid(axes={"not_a_real_field": [1, 2]})

    def test_unknown_fixed_field_rejected(self):
        with pytest.raises(SweepError):
            expand_grid(axes={"temperature": [300]}, fixed={"bogus": 1})

    def test_non_list_axis_rejected(self):
        with pytest.raises(SweepError):
            expand_grid(axes={"temperature": 300})


# --------------------------------------------------------------------------- #
# Per-point build
# --------------------------------------------------------------------------- #
class TestBuildPoint:
    def test_builds_and_writes_files(self, tmp_path):
        pt = expand_grid(
            axes={"temperature": [500]},
            fixed={"target_num_carbons": 60, "feedstock": "softwood"},
        )[0]
        res = build_point(pt, tmp_path, base_seed=0, max_retries=3)
        assert isinstance(res, PointResult)
        assert res.status in ("strict_pass", "fallback")
        # files exist
        assert res.gro_path and Path(res.gro_path).exists()
        assert res.top_path and Path(res.top_path).exists()
        assert res.itp_path and Path(res.itp_path).exists()
        # composition captured
        assert res.molecular_formula and res.molecular_formula.startswith("C60")
        assert res.num_carbons == 60
        assert res.H_C_ratio is not None and res.O_C_ratio is not None

    def test_skip_mode_writes_no_files_on_fail(self, tmp_path):
        # tiny/dense structures reliably fail strict validation here
        pt = expand_grid(
            axes={"temperature": [500]},
            fixed={"target_num_carbons": 100, "feedstock": "softwood"},
        )[0]
        res = build_point(pt, tmp_path, base_seed=0, max_retries=2,
                          on_validation_fail="skip")
        if res.status == "failed":
            assert res.gro_path is None
        else:
            # if it happened to pass strict, that's also acceptable
            assert res.status == "strict_pass"

    def test_invalid_fail_mode_rejected(self, tmp_path):
        pt = expand_grid(axes={"temperature": [500]})[0]
        with pytest.raises(SweepError):
            build_point(pt, tmp_path, on_validation_fail="nonsense")

    def test_crash_recorded_not_raised(self, tmp_path):
        # an infeasible request (more oxygens than the scaffold can hold) should
        # be caught and recorded, not propagate
        pt = expand_grid(
            axes={"functional_groups": [{"carboxyl": 500}]},
            fixed={"target_num_carbons": 20},
        )[0]
        res = build_point(pt, tmp_path, base_seed=0, max_retries=2)
        assert res.status in ("failed", "fallback", "strict_pass")
        # must not raise; result object always returned
        assert isinstance(res, PointResult)


# --------------------------------------------------------------------------- #
# Full sweep + manifest
# --------------------------------------------------------------------------- #
class TestRunSweep:
    def test_end_to_end_manifest(self, tmp_path):
        cfg = {
            "name": "unit_grid",
            "output_directory": str(tmp_path / "out"),
            "fixed": {"target_num_carbons": 60},
            "axes": {"temperature": [400, 600], "feedstock": ["softwood"]},
            "seed": 0,
            "max_retries": 3,
        }
        summary = run_sweep(cfg)
        assert summary["n_points"] == 2
        assert summary["n_built"] + summary["n_failed"] == 2

        # manifest files exist and parse
        csv_path = Path(summary["manifest_csv"])
        json_path = Path(summary["manifest_json"])
        assert csv_path.exists() and json_path.exists()

        with csv_path.open() as fh:
            rows = list(csv.DictReader(fh))
        assert len(rows) == 2
        assert "molecular_formula" in rows[0]
        assert "axis_temperature" in rows[0]

        payload = json.loads(json_path.read_text())
        assert payload["meta"]["name"] == "unit_grid"
        assert len(payload["results"]) == 2

    def test_progress_callback_invoked(self, tmp_path):
        seen = []
        cfg = {
            "name": "cb",
            "output_directory": str(tmp_path / "out"),
            "fixed": {"target_num_carbons": 60, "feedstock": "softwood"},
            "axes": {"temperature": [500]},
            "max_retries": 2,
        }
        run_sweep(cfg, progress_callback=lambda i, n, r: seen.append((i, n)))
        assert seen == [(1, 1)]

    def test_output_dir_override(self, tmp_path):
        cfg = {
            "name": "ov",
            "output_directory": str(tmp_path / "ignored"),
            "fixed": {"target_num_carbons": 60, "feedstock": "softwood"},
            "axes": {"temperature": [500]},
            "max_retries": 2,
        }
        override = str(tmp_path / "used")
        summary = run_sweep(cfg, output_directory=override)
        assert summary["output_directory"] == override
        assert (Path(override) / "manifest.csv").exists()

    def test_missing_axes_rejected(self):
        with pytest.raises(SweepError):
            run_sweep({"name": "x", "fixed": {"target_num_carbons": 60}})

    def test_unknown_toplevel_key_rejected(self):
        with pytest.raises(SweepError):
            run_sweep({"axes": {"temperature": [500]}, "typo_key": 1})


# --------------------------------------------------------------------------- #
# Config loading
# --------------------------------------------------------------------------- #
class TestLoadConfig:
    def test_load_json(self, tmp_path):
        p = tmp_path / "cfg.json"
        p.write_text(json.dumps({"axes": {"temperature": [500]}}))
        cfg = load_sweep_config(p)
        assert cfg["axes"]["temperature"] == [500]

    def test_load_yaml(self, tmp_path):
        pytest.importorskip("yaml")
        p = tmp_path / "cfg.yaml"
        p.write_text("axes:\n  temperature: [400, 500]\n")
        cfg = load_sweep_config(p)
        assert cfg["axes"]["temperature"] == [400, 500]

    def test_unknown_extension_rejected(self, tmp_path):
        p = tmp_path / "cfg.txt"
        p.write_text("nope")
        with pytest.raises(SweepError):
            load_sweep_config(p)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
class TestSweepCLI:
    def test_template_subcommand(self, capsys):
        from biochar.sweep_cli import main
        rc = main(["template"])
        out = capsys.readouterr().out
        assert rc == 0
        assert "axes:" in out and "temperature" in out

    def test_run_subcommand(self, tmp_path, capsys):
        from biochar.sweep_cli import main
        cfg = tmp_path / "cfg.json"
        cfg.write_text(json.dumps({
            "name": "cli_grid",
            "output_directory": str(tmp_path / "out"),
            "fixed": {"target_num_carbons": 60, "feedstock": "softwood"},
            "axes": {"temperature": [500]},
            "max_retries": 2,
        }))
        rc = main(["run", str(cfg), "--quiet"])
        assert rc in (0, 2)  # 0 all built, 2 if any failed — both are clean exits
        assert (tmp_path / "out" / "manifest.csv").exists()

    def test_run_missing_config(self, capsys):
        from biochar.sweep_cli import main
        rc = main(["run", "/no/such/file.yaml"])
        assert rc == 1
