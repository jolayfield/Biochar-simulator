"""Contract tests for rqm/cli-arguments.md — the four console entry points.

Covers the decisions made at the CLI layer and nowhere else: precedence between
a loaded config and the command line, exit status, and cross-argument
consistency. Structure-building cases use tiny carbon counts so the file stays
in the fast tier.
"""

import json
from pathlib import Path

import pytest

from biochar.cli.cli import _build_parser, main as gen_main
from biochar.cli.condensation_cli import main as condense_main
from biochar.cli.md_setup_cli import main as md_setup_main
from biochar.cli.sweep_cli import _TEMPLATE, main as sweep_main

# Small enough to build in a second; H/C 0.8 and O/C 0.0 are reachable at 10 C.
_FAST = ["--carbons", "10", "--hc-ratio", "0.8", "--oc-ratio", "0.0", "--seed", "42"]


def _save_config(tmp_path, extra_args, cfg_name="cfg.json"):
    """Run biochar-gen end-to-end with --save-config; return the config it wrote."""
    cfg = tmp_path / cfg_name
    gen_main(extra_args + ["--output-dir", str(tmp_path), "--save-config", str(cfg)])
    return cfg, json.loads(cfg.read_text())


def _resolve(argv, base=None):
    """Resolve argv (plus an already-loaded config) to a GeneratorConfig.

    This is the step that holds the precedence rule, so testing it directly
    keeps these cases off the generator: what is under test is which value wins,
    not whether the structure builds.
    """
    from biochar import GeneratorConfig
    from biochar.cli.cli import _resolve_config_dict, _supplied_dests

    args = _build_parser().parse_args(argv)
    resolved = _resolve_config_dict(args, _supplied_dests(argv), base or {})
    return GeneratorConfig.from_dict(resolved)


class TestEntryPointAPI:
    """rq-b31f2ffe — every main() is callable with argv and returns a status."""

    def test_mains_accept_argv_and_return_status(self, tmp_path):
        # rq-fa8c5dbc, rq-46f0dae6, rq-a2c4ed0e, rq-c46e42dc
        rc = gen_main(_FAST + ["--name", "TST", "--output-dir", str(tmp_path)])
        assert rc == 0
        assert isinstance(rc, int)

    def test_parser_is_inspectable_without_running(self):
        # rq-fa8c5dbc — _build_parser() exposes the accepted surface on its own.
        dests = {a.dest for a in _build_parser()._actions}
        assert {"carbons", "seed", "temperature", "load_config"} <= dests

    def test_every_option_is_forwarded_or_declared_local(self):
        # rq-fa8c5dbc — an option missing from the mapping parses and then does
        # nothing, which is the accepted-and-ignored failure this layer is most
        # prone to. Adding a flag must mean placing it somewhere below.
        from biochar.cli.cli import _CONFIG_FIELDS, _FG_DESTS

        # Consumed by main() itself rather than forwarded to GeneratorConfig.
        local = {"help", "output_dir", "basename", "save_config", "load_config",
                 "verbose", "debug"}
        dests = {a.dest for a in _build_parser()._actions}
        unaccounted = dests - set(_CONFIG_FIELDS) - set(_FG_DESTS) - local
        assert not unaccounted, f"options that reach nothing: {sorted(unaccounted)}"


class TestSavedConfigRoundTrips:
    """rq-a71b7f66 — a saved config reloads as the config that was saved."""

    def test_reloaded_config_matches_saved(self, tmp_path):
        # rq-b4e2b5cc
        cfg, saved = _save_config(tmp_path, [
            "--carbons", "10", "--hc-ratio", "0.8", "--oc-ratio", "0.0",
            "--seed", "1234", "--name", "TST",
            "--temperature", "600", "--feedstock", "softwood",
        ])
        _, reloaded = _save_config(
            tmp_path, ["--load-config", str(cfg)], cfg_name="cfg2.json"
        )
        assert reloaded == saved

    def test_seed_survives_the_reload(self):
        # rq-b4e2b5cc — the seed is the whole reproducibility mechanism, and the
        # field a user is least likely to re-type because the file carries it.
        assert _resolve([], base={"seed": 1234}).seed == 1234

    def test_derived_composition_keeps_its_conditions(self):
        # rq-e2c50b04 — the numbers must not outlive the temperature they came from.
        saved = _resolve(["--temperature", "600", "--feedstock", "softwood"])
        reloaded = _resolve([], base={"temperature": saved.temperature,
                                      "feedstock": saved.feedstock,
                                      "H_C_ratio": saved.H_C_ratio})
        assert reloaded.temperature == 600.0
        assert reloaded.feedstock == "softwood"
        assert reloaded.H_C_ratio == pytest.approx(saved.H_C_ratio)


class TestDefaultIsNotARequest:
    """rq-e6faa000 — only flags the user typed override a loaded config."""

    def test_unsupplied_flag_leaves_loaded_value_alone(self):
        # rq-e6099cd4
        resolved = _resolve([], base={"target_num_carbons": 14, "molecule_name": "LDD"})
        assert resolved.target_num_carbons == 14
        assert resolved.molecule_name == "LDD"

    def test_supplied_flag_overrides_loaded_value(self):
        # rq-b8fb0aeb
        resolved = _resolve(["--carbons", "10"], base={"target_num_carbons": 14})
        assert resolved.target_num_carbons == 10

    def test_flag_supplied_at_its_default_still_overrides(self):
        # rq-1a4b6dae — typing the default value is a request like any other,
        # and is exactly the case a "is it non-default?" heuristic gets wrong.
        default_carbons = _build_parser().parse_args([]).carbons
        resolved = _resolve(["--carbons", str(default_carbons)],
                            base={"target_num_carbons": 14})
        assert resolved.target_num_carbons == default_carbons

    def test_supplied_dests_reports_only_typed_options(self):
        # rq-1a4b6dae — the mechanism the precedence rule rests on.
        from biochar.cli.cli import _supplied_dests

        assert _supplied_dests(["--carbons", "50"]) == {"carbons"}
        assert _supplied_dests([]) == set()
        assert "seed" in _supplied_dests(["--seed", "1"])

    def test_defaults_apply_without_a_config_file(self):
        # rq-b81e79b5 — a default is a fallback, and still reaches the config.
        parsed = _build_parser().parse_args([])
        resolved = _resolve([])
        assert resolved.target_num_carbons == parsed.carbons
        assert resolved.molecule_name == parsed.name
        assert resolved.charge_method == parsed.charge_method
        assert resolved.defect_fraction == parsed.defects
        assert resolved.num_pyridinic == parsed.pyridinic


class TestFunctionalGroupsAreASet:
    """rq-fd1e5f16 — the six group flags describe one dictionary."""

    def test_any_group_flag_replaces_the_loaded_set(self):
        # rq-92f2fb37 — not merged: the counts compete for the same edge sites,
        # so a set assembled from two sources is a request nobody made.
        resolved = _resolve(["--phenolic", "2"],
                            base={"functional_groups": {"carboxyl": 1}})
        assert resolved.functional_groups == {"phenolic": 2}

    def test_no_group_flag_leaves_the_loaded_set_intact(self):
        # rq-c72f3fd1
        resolved = _resolve([], base={"functional_groups": {"carboxyl": 1}})
        assert resolved.functional_groups == {"carboxyl": 1}


def _write_manifest(path, rows):
    header = "label,status,gro_path,top_path\n"
    path.write_text(header + "".join(",".join(r) + "\n" for r in rows))


class TestEmptyResultIsNotSuccess:
    """rq-1eae7bef — a command that wrote nothing does not exit 0."""

    def test_manifest_yielding_nothing_exits_non_zero(self, tmp_path):
        # rq-2b6f2c93
        manifest = tmp_path / "manifest.csv"
        _write_manifest(manifest, [["ptA", "failed", "", ""], ["ptB", "failed", "", ""]])
        out_root = tmp_path / "runs"
        rc = md_setup_main([str(manifest), "--output-root", str(out_root)])
        assert rc != 0
        assert not any(out_root.glob("*/")) if out_root.exists() else True

    def test_manifest_with_one_usable_row_succeeds(self, tmp_path):
        # rq-59fb4726 — partial skips stay a success.
        built = tmp_path / "built"
        built.mkdir()
        rc = gen_main(_FAST + ["--name", "TST", "--output-dir", str(built),
                               "--basename", "ptA"])
        assert rc == 0
        manifest = tmp_path / "manifest.csv"
        _write_manifest(manifest, [
            ["ptA", "strict_pass", str(built / "ptA.gro"), str(built / "ptA.top")],
            ["ptB", "failed", "", ""],
        ])
        out_root = tmp_path / "runs"
        rc = md_setup_main([str(manifest), "--output-root", str(out_root)])
        assert rc == 0
        assert (out_root / "ptA").is_dir()


class TestSurfaceRepeatIsInRange:
    """rq-a48afb56 — no setup pointing at a repeat that will not exist."""

    _COMMON = ["--copies", "4", "--htt", "600", "--repeats", "3"]
    _MOL = ["--carbons", "10", "--hc", "0.8", "--oc", "0.0", "--seed", "1", "--name", "TST"]

    def test_repeat_above_the_count_is_refused_before_generating(self, tmp_path):
        # rq-c53b3aeb
        out = tmp_path / "run"
        rc = condense_main(["generate", "--output-dir", str(out)]
                           + self._COMMON + self._MOL + ["--which-repeat", "7"])
        assert rc == 1
        assert not out.exists(), "refused request must not generate or write anything"

    def test_repeat_below_one_is_refused(self, tmp_path):
        # rq-a4f5a7d7
        out = tmp_path / "run"
        rc = condense_main(["generate", "--output-dir", str(out)]
                           + self._COMMON + self._MOL + ["--which-repeat", "0"])
        assert rc == 1
        assert not out.exists()

    def test_repeat_in_range_is_accepted(self, tmp_path):
        # rq-2ed54e2c
        out = tmp_path / "run"
        rc = condense_main(["generate", "--output-dir", str(out)]
                           + self._COMMON + self._MOL + ["--which-repeat", "3"])
        assert rc == 0
        surface = out / "run_surface.sh"
        assert surface.exists()
        assert "rep_3/final" in surface.read_text()


class TestTemplateNamesRealOptions:
    """rq-eb1c8214 — the starter template only offers what the config accepts."""

    def test_every_charge_method_named_in_the_template_is_accepted(self):
        # rq-cba5c290
        import re

        from biochar import GeneratorConfig

        named = set()
        for line in _TEMPLATE.splitlines():
            if "charge_method" in line and "#" in line:
                named.update(re.findall(r"\b(opls|ml|qm|gasteiger|am1bcc|mmff)\b",
                                        line.split("#", 1)[1]))
        assert named, "the template should document the charge backends it offers"
        for method in sorted(named):
            GeneratorConfig(target_num_carbons=10, charge_method=method,
                            molecule_name="TST")

    def test_template_expands_to_a_grid(self, tmp_path):
        # rq-32e0e0e6
        from biochar.workflows.sweep import expand_grid

        yaml = pytest.importorskip("yaml")
        cfg = yaml.safe_load(_TEMPLATE)
        points = expand_grid(cfg["axes"])
        assert len(points) == 10  # 5 temperatures x 2 feedstocks


class TestDiagnosticsStayOffStdout:
    """rq-ac71d20f — results on stdout, diagnostics on stderr."""

    def test_missing_config_is_reported_and_fails(self, tmp_path, capsys):
        # rq-b697d2e6
        rc = gen_main(["--load-config", str(tmp_path / "nope.json"),
                       "--name", "TST", "--output-dir", str(tmp_path)])
        assert rc == 1
        captured = capsys.readouterr()
        assert "Error loading config" in captured.err
        assert "Error loading config" not in captured.out

    def test_template_writes_only_the_template_to_stdout(self, capsys):
        # rq-c8d1f9e0 — `biochar-sweep template > sweep.yaml` depends on this:
        # anything else on stdout lands in the file the user then runs.
        yaml = pytest.importorskip("yaml")
        rc = sweep_main(["template"])
        assert rc == 0
        captured = capsys.readouterr()
        cfg = yaml.safe_load(captured.out)
        assert cfg["axes"] and cfg["output_directory"]


class TestCrossReferences:
    """rq-fb2e2ae7 — the neighbouring documents these flags resolve into."""

    def test_referenced_requirement_documents_exist(self):
        # rq-fb2e2ae7, rq-b9751811
        rqm = Path(__file__).resolve().parents[1] / "rqm"
        for name in ("generation-config", "temperature-model", "parameter-sweep",
                     "md-setup", "condensation-annealing"):
            assert (rqm / f"{name}.md").exists()
