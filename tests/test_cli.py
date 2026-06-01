"""
Tests for biochar/cli.py — argument parser and main() entry point.
"""

import json
import pytest
from pathlib import Path

from biochar.cli import _build_parser, main


class TestParser:
    """Verify _build_parser() produces correct defaults and accepts all flags."""

    def test_defaults(self):
        parser = _build_parser()
        args = parser.parse_args([])
        assert args.carbons == 50
        assert args.hc_ratio is None        # None → GeneratorConfig applies default
        assert args.oc_ratio is None
        assert args.aromaticity is None
        assert args.temperature is None
        assert args.feedstock is None
        assert args.defects == 0.0
        assert args.seed is None
        assert args.output_dir == "."
        assert args.basename == "biochar"

    def test_temperature_flag(self):
        args = _build_parser().parse_args(["--temperature", "600"])
        assert args.temperature == 600.0

    def test_feedstock_flag(self):
        args = _build_parser().parse_args(["--temperature", "600", "--feedstock", "softwood"])
        assert args.feedstock == "softwood"

    def test_feedstock_invalid_rejected(self):
        """An unknown feedstock must be rejected at parse time."""
        import pytest
        with pytest.raises(SystemExit):
            _build_parser().parse_args(["--feedstock", "unicorn"])

    def test_carbons_flag(self):
        args = _build_parser().parse_args(["--carbons", "80"])
        assert args.carbons == 80

    def test_hc_oc_flags(self):
        args = _build_parser().parse_args(["--hc-ratio", "0.4", "--oc-ratio", "0.0"])
        assert args.hc_ratio == 0.4
        assert args.oc_ratio == 0.0

    def test_defects_flag(self):
        args = _build_parser().parse_args(["--defects", "0.2"])
        assert args.defects == 0.2

    def test_seed_flag(self):
        args = _build_parser().parse_args(["--seed", "99"])
        assert args.seed == 99

    def test_name_flag(self):
        args = _build_parser().parse_args(["--name", "BC400"])
        assert args.name == "BC400"

    def test_functional_group_flags(self):
        args = _build_parser().parse_args([
            "--phenolic", "2", "--carboxyl", "1", "--ether", "1", "--amino", "1"
        ])
        assert args.phenolic == 2
        assert args.carboxyl == 1
        assert args.ether == 1
        assert args.amino == 1

    def test_functional_group_defaults_are_none(self):
        """Without flags, all FG counts must be None (not 0)."""
        args = _build_parser().parse_args([])
        assert args.phenolic is None
        assert args.carboxyl is None
        assert args.ether is None
        assert args.amino is None

    def test_verbose_and_debug_flags(self):
        args = _build_parser().parse_args(["--verbose"])
        assert args.verbose is True
        args2 = _build_parser().parse_args(["--debug"])
        assert args2.debug is True


class TestMain:
    """Integration tests for main() — run with small target for speed."""

    _FAST_ARGS = [
        "--carbons", "10",
        "--hc-ratio", "0.8",
        "--oc-ratio", "0.0",
        "--seed", "42",
    ]

    def test_success_returns_zero(self, tmp_path):
        rc = main(self._FAST_ARGS + ["--output-dir", str(tmp_path), "--name", "TST"])
        assert rc == 0

    def test_output_files_created(self, tmp_path):
        main(self._FAST_ARGS + [
            "--output-dir", str(tmp_path),
            "--name", "TST",
            "--basename", "mybc",
        ])
        assert (tmp_path / "mybc.gro").exists()
        assert (tmp_path / "mybc.top").exists()
        assert (tmp_path / "mybc.itp").exists()

    def test_invalid_molecule_name_returns_one(self, tmp_path):
        """molecule_name > 5 chars must cause a config error (return 1)."""
        rc = main(self._FAST_ARGS + [
            "--output-dir", str(tmp_path),
            "--name", "TOOLONG",
        ])
        assert rc == 1

    def test_save_config_creates_json(self, tmp_path):
        cfg_file = str(tmp_path / "cfg.json")
        main(self._FAST_ARGS + [
            "--output-dir", str(tmp_path),
            "--name", "TST",
            "--save-config", cfg_file,
        ])
        assert Path(cfg_file).exists()
        data = json.loads(Path(cfg_file).read_text())
        assert "target_num_carbons" in data

    def test_save_config_json_values_match_flags(self, tmp_path):
        cfg_file = str(tmp_path / "cfg.json")
        main([
            "--carbons", "10",
            "--hc-ratio", "0.7",
            "--oc-ratio", "0.0",
            "--seed", "7",
            "--name", "TST",
            "--output-dir", str(tmp_path),
            "--save-config", cfg_file,
        ])
        data = json.loads(Path(cfg_file).read_text())
        assert data["target_num_carbons"] == 10
        assert abs(data["H_C_ratio"] - 0.7) < 1e-6
        assert data["seed"] == 7

    def test_load_config_used_as_base(self, tmp_path):
        """--load-config supplies non-CLI fields (e.g. strict=False) without crashing."""
        cfg_file = str(tmp_path / "base.json")
        Path(cfg_file).write_text(json.dumps({"strict": False}))
        rc = main([
            "--load-config", cfg_file,
            "--carbons", "10",
            "--hc-ratio", "0.8",
            "--oc-ratio", "0.0",
            "--seed", "42",
            "--name", "TST",
            "--output-dir", str(tmp_path),
        ])
        assert rc == 0

    def test_load_config_missing_file_returns_one(self, tmp_path):
        rc = main([
            "--load-config", str(tmp_path / "does_not_exist.json"),
            "--name", "TST",
            "--output-dir", str(tmp_path),
        ])
        assert rc == 1

    def test_load_config_invalid_json_returns_one(self, tmp_path):
        bad = tmp_path / "bad.json"
        bad.write_text("{not valid json}")
        rc = main([
            "--load-config", str(bad),
            "--name", "TST",
            "--output-dir", str(tmp_path),
        ])
        assert rc == 1

    def test_functional_group_flags_propagate(self, tmp_path):
        """--phenolic and --amino flags must result in a successful generation."""
        rc = main([
            "--carbons", "10",
            "--hc-ratio", "0.9",
            "--oc-ratio", "0.0",
            "--seed", "1",
            "--name", "TST",
            "--amino", "1",
            "--output-dir", str(tmp_path),
        ])
        assert rc == 0

    def test_output_dir_is_created(self, tmp_path):
        """main() should create the output directory if it doesn't exist."""
        new_dir = tmp_path / "subdir" / "deep"
        rc = main(self._FAST_ARGS + [
            "--output-dir", str(new_dir),
            "--name", "TST",
        ])
        assert rc == 0
        assert new_dir.exists()


class TestTemperatureAndFeedstock:
    """CLI integration for --temperature and --feedstock."""

    def test_temperature_only_succeeds(self, tmp_path):
        # 300 °C → H/C ~0.80 (pooled), achievable by C10 naphthalene.
        # --oc-ratio 0.0 overrides the temperature-derived O/C (high at 300 °C,
        # hard to place on a 10-carbon skeleton).
        rc = main([
            "--temperature", "300",
            "--oc-ratio", "0.0",
            "--carbons", "10",
            "--seed", "42",
            "--name", "TST",
            "--output-dir", str(tmp_path),
        ])
        assert rc == 0

    def test_temperature_with_feedstock_succeeds(self, tmp_path):
        rc = main([
            "--temperature", "300",
            "--feedstock", "softwood",
            "--oc-ratio", "0.0",
            "--carbons", "10",
            "--seed", "42",
            "--name", "TST",
            "--output-dir", str(tmp_path),
        ])
        assert rc == 0

    def test_temperature_derives_composition_in_saved_config(self, tmp_path):
        """With only --temperature, saved config should have resolved H/C and O/C."""
        cfg_file = str(tmp_path / "cfg.json")
        main([
            "--temperature", "600",
            "--feedstock", "softwood",
            "--carbons", "10",
            "--seed", "1",
            "--name", "TST",
            "--output-dir", str(tmp_path),
            "--save-config", cfg_file,
        ])
        data = json.loads(Path(cfg_file).read_text())
        assert data.get("temperature") == 600.0
        assert data.get("feedstock") == "softwood"
        # GeneratorConfig fills these from the model; they must be numeric and plausible
        assert 0.0 < data["H_C_ratio"] < 2.0
        assert 0.0 <= data["O_C_ratio"] < 1.0

    def test_explicit_hc_overrides_temperature(self, tmp_path):
        """When --hc-ratio is given alongside --temperature, it must win."""
        cfg_file = str(tmp_path / "cfg.json")
        main([
            "--temperature", "600",
            "--hc-ratio", "0.99",
            "--carbons", "10",
            "--seed", "1",
            "--name", "TST",
            "--output-dir", str(tmp_path),
            "--save-config", cfg_file,
        ])
        data = json.loads(Path(cfg_file).read_text())
        assert abs(data["H_C_ratio"] - 0.99) < 1e-6

    def test_no_composition_args_uses_defaults(self, tmp_path):
        """Without temperature or ratio flags, historical defaults (0.5/0.1/90) apply."""
        cfg_file = str(tmp_path / "cfg.json")
        main([
            "--carbons", "10",
            "--seed", "1",
            "--name", "TST",
            "--output-dir", str(tmp_path),
            "--save-config", cfg_file,
        ])
        data = json.loads(Path(cfg_file).read_text())
        assert abs(data["H_C_ratio"] - 0.5) < 1e-6
        assert abs(data["O_C_ratio"] - 0.1) < 1e-6
        assert abs(data["aromaticity_percent"] - 90.0) < 1e-6
