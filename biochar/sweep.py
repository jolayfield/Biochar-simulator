"""Declarative parameter-sweep driver for biochar structure generation.

The single-molecule API (:func:`biochar.generate_biochar`) and the list API
(:func:`biochar.generate_biochar_series`) both require the caller to enumerate
every configuration by hand: one config dict per structure, a unique <=5-char
residue name for each, and no record of which structures actually passed
geometry/composition validation.

This module turns that hand-run loop into a declarative, reproducible pipeline.
A sweep is described by a small config (dict or YAML/JSON file) that names a set
of *axes*; the driver expands their Cartesian product, builds every grid point
with automatic seed-retry and a strict -> non-strict fallback, writes the
GROMACS files into an organised output tree, and emits a **manifest**
(``manifest.csv`` + ``manifest.json``) recording, for every structure: the axis
values, the achieved composition (formula, H/C, O/C, functional groups),
validation status, the seed that was used, and the output file paths.

Example
-------
>>> from biochar.sweep import run_sweep
>>> summary = run_sweep({
...     "name": "pfas_grid",
...     "output_directory": "sweep_out",
...     "fixed": {"target_num_carbons": 100},
...     "axes": {
...         "temperature": [400, 500, 600, 700],
...         "feedstock": ["softwood", "hardwood"],
...     },
... })
>>> summary["n_points"], summary["n_built"]
(8, 8)

Each grid point becomes one row in the manifest and one subdirectory under
``sweep_out/structures/``.
"""

from __future__ import annotations

import csv
import itertools
import json
import logging
import time
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

from .biochar_generator import (
    BiocharGenerator,
    GeneratorConfig,
    ValidationError,
)

logger = logging.getLogger(__name__)

__all__ = [
    "SweepError",
    "GridPoint",
    "PointResult",
    "expand_grid",
    "build_point",
    "run_sweep",
    "load_sweep_config",
]

# Fields of GeneratorConfig that may appear in `fixed`, `axes`, or a grid point.
# Anything outside this set in a sweep config is a typo and is rejected early.
_GENCONFIG_FIELDS = frozenset(GeneratorConfig.__dataclass_fields__.keys())

# Driver-level keys recognised at the top level of a sweep config.
_TOPLEVEL_KEYS = frozenset({
    "name", "output_directory", "axes", "fixed", "seed", "max_retries",
    "on_validation_fail", "name_template", "write_files",
})

_VALID_FAIL_MODES = ("fallback", "skip", "strict")


class SweepError(Exception):
    """Raised when a sweep configuration is malformed."""


# --------------------------------------------------------------------------- #
# Grid expansion
# --------------------------------------------------------------------------- #
@dataclass
class GridPoint:
    """One point in the parameter sweep.

    Attributes
    ----------
    index:
        0-based position in the expanded grid.
    label:
        Short, filesystem-safe identifier derived from the varying axis values
        (e.g. ``"T400_softwood"``). Used for the output subdirectory.
    molecule_name:
        <=5-char residue/molecule name written into the GROMACS topology.
    axis_values:
        Mapping of axis name -> value for *this* point (the varying part only).
    config_kwargs:
        Full keyword set passed to :class:`GeneratorConfig` (fixed + axis).
    """

    index: int
    label: str
    molecule_name: str
    axis_values: Dict[str, Any]
    config_kwargs: Dict[str, Any]


def _sanitize(value: Any) -> str:
    """Make an axis value safe for a filesystem label."""
    s = str(value)
    for bad in " /\\\t\n":
        s = s.replace(bad, "_")
    # collapse functional-group dicts to a compact tag
    return s.replace("{", "").replace("}", "").replace("'", "").replace(":", "").replace(",", "_")


def _make_label(axis_values: Dict[str, Any]) -> str:
    """Build a short, readable label from the varying axis values."""
    if not axis_values:
        return "point"
    parts = []
    for k, v in axis_values.items():
        tag = _sanitize(v)
        # Compact common axes: temperature -> T400, feedstock -> first token
        if k == "temperature":
            parts.append(f"T{tag}")
        elif k == "target_num_carbons":
            parts.append(f"C{tag}")
        else:
            parts.append(tag)
    return "_".join(parts)[:48]


def _validate_keys(name: str, mapping: Dict[str, Any]) -> None:
    unknown = set(mapping) - _GENCONFIG_FIELDS
    if unknown:
        raise SweepError(
            f"{name} contains keys that are not GeneratorConfig fields: "
            f"{sorted(unknown)}. Valid fields: {sorted(_GENCONFIG_FIELDS)}"
        )


def expand_grid(
    axes: Dict[str, Sequence[Any]],
    fixed: Optional[Dict[str, Any]] = None,
    name_template: str = "BC{i:03d}",
) -> List[GridPoint]:
    """Expand a set of sweep axes into a list of :class:`GridPoint`.

    Parameters
    ----------
    axes:
        Mapping of GeneratorConfig field name -> list of values to sweep.
        The Cartesian product of all axes is taken.
    fixed:
        Parameters applied identically to every grid point.
    name_template:
        ``str.format`` template for the <=5-char molecule/residue name. The
        running index is available as ``{i}``; axis values are available by
        name (e.g. ``"BC{temperature}"``). The result is truncated to 5 chars.

    Returns
    -------
    list of GridPoint
    """
    fixed = dict(fixed or {})
    _validate_keys("fixed", fixed)
    if not axes:
        raise SweepError("`axes` is empty: a sweep needs at least one axis.")
    _validate_keys("axes", axes)

    for k, vals in axes.items():
        if not isinstance(vals, (list, tuple)) or len(vals) == 0:
            raise SweepError(f"axis '{k}' must be a non-empty list, got {vals!r}")

    axis_names = list(axes.keys())
    value_lists = [list(axes[k]) for k in axis_names]

    points: List[GridPoint] = []
    for i, combo in enumerate(itertools.product(*value_lists)):
        axis_values = dict(zip(axis_names, combo))
        # axis values override fixed (axes win on conflict)
        cfg = dict(fixed)
        cfg.update(axis_values)

        try:
            mol_name = name_template.format(i=i, **axis_values)
        except (KeyError, IndexError) as exc:
            raise SweepError(
                f"name_template {name_template!r} references a field not in "
                f"axes/index: {exc}"
            ) from exc
        mol_name = mol_name[:5]

        label = _make_label(axis_values)
        cfg.setdefault("molecule_name", mol_name)
        points.append(
            GridPoint(
                index=i,
                label=label,
                molecule_name=mol_name,
                axis_values=axis_values,
                config_kwargs=cfg,
            )
        )
    return points


# --------------------------------------------------------------------------- #
# Per-point build with seed-retry + fallback
# --------------------------------------------------------------------------- #
@dataclass
class PointResult:
    """Outcome of building a single grid point."""

    index: int
    label: str
    molecule_name: str
    axis_values: Dict[str, Any]
    status: str                       # strict_pass | fallback | failed
    seed_used: Optional[int]
    n_attempts: int
    # achieved composition (None if the point could not be built at all)
    molecular_formula: Optional[str] = None
    molecular_weight: Optional[float] = None
    num_carbons: Optional[int] = None
    H_C_ratio: Optional[float] = None
    O_C_ratio: Optional[float] = None
    functional_groups: Optional[Dict[str, int]] = None
    ring_composition: Optional[Dict[str, int]] = None
    validation_errors: List[str] = field(default_factory=list)
    validation_warnings: List[str] = field(default_factory=list)
    gro_path: Optional[str] = None
    top_path: Optional[str] = None
    itp_path: Optional[str] = None
    error: Optional[str] = None       # exception text if build crashed
    elapsed_s: float = 0.0

    def to_row(self) -> Dict[str, Any]:
        """Flatten to a manifest row (JSON-friendly scalars)."""
        row = {
            "index": self.index,
            "label": self.label,
            "molecule_name": self.molecule_name,
            "status": self.status,
            "seed_used": self.seed_used,
            "n_attempts": self.n_attempts,
            "molecular_formula": self.molecular_formula,
            "molecular_weight": self.molecular_weight,
            "num_carbons": self.num_carbons,
            "H_C_ratio": self.H_C_ratio,
            "O_C_ratio": self.O_C_ratio,
            "functional_groups": json.dumps(self.functional_groups) if self.functional_groups else "",
            "ring_composition": json.dumps(self.ring_composition) if self.ring_composition else "",
            "n_validation_errors": len(self.validation_errors),
            "n_validation_warnings": len(self.validation_warnings),
            "gro_path": self.gro_path or "",
            "top_path": self.top_path or "",
            "itp_path": self.itp_path or "",
            "elapsed_s": round(self.elapsed_s, 3),
            "error": self.error or "",
        }
        # axis values as their own columns for easy pivoting
        for k, v in self.axis_values.items():
            row[f"axis_{k}"] = json.dumps(v) if isinstance(v, (dict, list)) else v
        return row


def _silence_rdkit() -> None:
    """Quiet RDKit's C++ logger (kekulization warnings bypass Python redirect)."""
    try:
        from rdkit import RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:  # pragma: no cover - rdkit always present in practice
        pass


def _build_once(
    config_kwargs: Dict[str, Any],
    seed: int,
    strict: bool,
    out_dir: Path,
    basename: str,
    write_files: bool,
    quiet: bool,
) -> Tuple[BiocharGenerator, Any, Tuple]:
    """Run one generation attempt. Returns (generator, composition, paths).

    Raises ValidationError in strict mode if the structure does not validate.
    """
    if quiet:
        _silence_rdkit()
    cfg = GeneratorConfig(**{**config_kwargs, "seed": seed, "strict": strict})
    gen = BiocharGenerator(cfg)
    sink = StringIO()
    if quiet:
        with redirect_stdout(sink), redirect_stderr(sink):
            mol, coords, comp = gen.generate()
    else:
        mol, coords, comp = gen.generate()

    paths: Tuple = (None, None, None)
    if write_files:
        if quiet:
            with redirect_stdout(sink), redirect_stderr(sink):
                paths = gen.export_gromacs(output_directory=str(out_dir), basename=basename)
        else:
            paths = gen.export_gromacs(output_directory=str(out_dir), basename=basename)
    return gen, comp, paths


def build_point(
    point: GridPoint,
    output_root: Path,
    base_seed: int = 0,
    max_retries: int = 8,
    on_validation_fail: str = "fallback",
    write_files: bool = True,
    quiet: bool = True,
) -> PointResult:
    """Build one grid point with seed-retry and strict -> fallback handling.

    Strategy
    --------
    1. Try ``strict=True`` with seeds ``base_seed, base_seed+1, ...`` up to
       ``max_retries`` attempts. The first seed that passes validation wins
       (status ``strict_pass``).
    2. If none pass and ``on_validation_fail == "fallback"`` (default), rebuild
       with ``strict=False`` at ``base_seed`` and keep the structure, recording
       the validation errors (status ``fallback``).
    3. If ``on_validation_fail == "skip"``, record the failure and write no
       files (status ``failed``).
    4. If ``on_validation_fail == "strict"``, the point is left unbuilt and
       reported as ``failed`` (the sweep continues; nothing is raised).

    A non-validation crash (e.g. an infeasible composition request) is caught
    and recorded as ``failed`` with the exception text; the sweep continues.
    """
    if on_validation_fail not in _VALID_FAIL_MODES:
        raise SweepError(
            f"on_validation_fail must be one of {_VALID_FAIL_MODES}, "
            f"got {on_validation_fail!r}"
        )

    point_dir = output_root / "structures" / f"{point.index:03d}_{point.label}"
    basename = point.label
    t0 = time.time()
    attempts = 0
    last_crash: Optional[str] = None

    # --- Phase 1: strict, seed-retry --------------------------------------- #
    for k in range(max(1, max_retries)):
        seed = base_seed + k
        attempts += 1
        try:
            gen, comp, paths = _build_once(
                point.config_kwargs, seed, True, point_dir, basename, write_files, quiet
            )
            # strict success
            return _ok_result(point, gen, comp, paths, "strict_pass", seed,
                              attempts, time.time() - t0)
        except ValidationError:
            continue
        except Exception as exc:  # infeasible request, etc. — do not seed-retry
            last_crash = f"{type(exc).__name__}: {exc}"
            logger.warning("Grid point %d crashed (non-validation): %s",
                           point.index, last_crash)
            break

    # --- Phase 2: handle non-passing point --------------------------------- #
    if last_crash is not None:
        return PointResult(
            index=point.index, label=point.label,
            molecule_name=point.molecule_name, axis_values=point.axis_values,
            status="failed", seed_used=None, n_attempts=attempts,
            error=last_crash, elapsed_s=time.time() - t0,
        )

    if on_validation_fail == "fallback":
        try:
            gen, comp, paths = _build_once(
                point.config_kwargs, base_seed, False, point_dir, basename,
                write_files, quiet
            )
            attempts += 1
            return _ok_result(point, gen, comp, paths, "fallback", base_seed,
                              attempts, time.time() - t0)
        except Exception as exc:
            return PointResult(
                index=point.index, label=point.label,
                molecule_name=point.molecule_name, axis_values=point.axis_values,
                status="failed", seed_used=None, n_attempts=attempts,
                error=f"fallback build crashed: {type(exc).__name__}: {exc}",
                elapsed_s=time.time() - t0,
            )

    # skip / strict: report failure without files; capture the last report
    errors: List[str] = []
    try:
        gen, comp, _ = _build_once(
            point.config_kwargs, base_seed, False, point_dir, basename,
            write_files=False, quiet=quiet,
        )
        rep = getattr(gen, "validation_report", None)
        if rep:
            errors = list(rep[1])
    except Exception:
        pass
    return PointResult(
        index=point.index, label=point.label,
        molecule_name=point.molecule_name, axis_values=point.axis_values,
        status="failed", seed_used=None, n_attempts=attempts,
        validation_errors=errors, elapsed_s=time.time() - t0,
    )


def _ok_result(point, gen, comp, paths, status, seed, attempts, elapsed) -> PointResult:
    rep = getattr(gen, "validation_report", None) or (True, [], [], {})
    _, errors, warnings, _metrics = rep
    return PointResult(
        index=point.index,
        label=point.label,
        molecule_name=point.molecule_name,
        axis_values=point.axis_values,
        status=status,
        seed_used=seed,
        n_attempts=attempts,
        molecular_formula=getattr(comp, "molecular_formula", None),
        molecular_weight=getattr(comp, "molecular_weight", None),
        num_carbons=getattr(comp, "num_carbons", None),
        H_C_ratio=getattr(comp, "H_C_ratio", None),
        O_C_ratio=getattr(comp, "O_C_ratio", None),
        functional_groups=dict(getattr(comp, "functional_groups", {}) or {}),
        ring_composition=dict(getattr(gen, "ring_composition", {}) or {}),
        validation_errors=list(errors),
        validation_warnings=list(warnings),
        gro_path=str(paths[0]) if paths[0] else None,
        top_path=str(paths[1]) if paths[1] else None,
        itp_path=str(paths[2]) if paths[2] else None,
        elapsed_s=elapsed,
    )


# --------------------------------------------------------------------------- #
# Config loading
# --------------------------------------------------------------------------- #
def load_sweep_config(path: Union[str, Path]) -> Dict[str, Any]:
    """Load a sweep config from a YAML (.yml/.yaml) or JSON (.json) file."""
    path = Path(path)
    text = path.read_text()
    if path.suffix.lower() in (".yml", ".yaml"):
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover
            raise SweepError(
                "PyYAML is required to read YAML sweep configs "
                "(`pip install pyyaml`), or use JSON."
            ) from exc
        return yaml.safe_load(text)
    if path.suffix.lower() == ".json":
        return json.loads(text)
    raise SweepError(f"Unrecognised config extension: {path.suffix} (use .yaml or .json)")


def _check_toplevel(config: Dict[str, Any]) -> None:
    unknown = set(config) - _TOPLEVEL_KEYS
    if unknown:
        raise SweepError(
            f"Unknown top-level keys in sweep config: {sorted(unknown)}. "
            f"Recognised: {sorted(_TOPLEVEL_KEYS)}"
        )
    if "axes" not in config:
        raise SweepError("sweep config must define `axes`.")


# --------------------------------------------------------------------------- #
# Manifest writers
# --------------------------------------------------------------------------- #
def _write_manifest(results: List[PointResult], output_root: Path,
                    meta: Dict[str, Any]) -> Tuple[Path, Path]:
    rows = [r.to_row() for r in results]
    # union of columns (axis columns vary by sweep)
    cols: List[str] = []
    for row in rows:
        for k in row:
            if k not in cols:
                cols.append(k)

    csv_path = output_root / "manifest.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in cols})

    json_path = output_root / "manifest.json"
    payload = {
        "meta": meta,
        "results": rows,
    }
    json_path.write_text(json.dumps(payload, indent=2, default=str))
    return csv_path, json_path


# --------------------------------------------------------------------------- #
# Top-level driver
# --------------------------------------------------------------------------- #
def run_sweep(
    config: Union[Dict[str, Any], str, Path],
    output_directory: Optional[str] = None,
    progress_callback: Optional[Callable[[int, int, PointResult], None]] = None,
    quiet: bool = True,
) -> Dict[str, Any]:
    """Run a full parameter sweep from a config dict or file.

    Parameters
    ----------
    config:
        Either a sweep-config dict, or a path to a YAML/JSON config file.
    output_directory:
        Overrides ``config["output_directory"]`` if given.
    progress_callback:
        Optional ``f(i, n, point_result)`` called after each grid point.
    quiet:
        Suppress the generator's per-structure stdout/stderr chatter.

    Returns
    -------
    dict
        Summary with keys ``name``, ``output_directory``, ``n_points``,
        ``n_built``, ``n_strict_pass``, ``n_fallback``, ``n_failed``,
        ``manifest_csv``, ``manifest_json``, and ``results`` (list of
        :class:`PointResult`).
    """
    if isinstance(config, (str, Path)):
        config = load_sweep_config(config)
    if not isinstance(config, dict):
        raise SweepError(f"config must be a dict or path, got {type(config)}")

    _check_toplevel(config)

    name = config.get("name", "sweep")
    out = output_directory or config.get("output_directory", name)
    output_root = Path(out)
    output_root.mkdir(parents=True, exist_ok=True)

    base_seed = int(config.get("seed", 0))
    max_retries = int(config.get("max_retries", 8))
    on_fail = config.get("on_validation_fail", "fallback")
    name_template = config.get("name_template", "BC{i:03d}")
    write_files = bool(config.get("write_files", True))

    points = expand_grid(
        axes=config["axes"],
        fixed=config.get("fixed"),
        name_template=name_template,
    )

    logger.info("Sweep '%s': %d grid points -> %s", name, len(points), output_root)

    results: List[PointResult] = []
    for j, point in enumerate(points):
        res = build_point(
            point, output_root,
            base_seed=base_seed,
            max_retries=max_retries,
            on_validation_fail=on_fail,
            write_files=write_files,
            quiet=quiet,
        )
        results.append(res)
        if progress_callback:
            progress_callback(j + 1, len(points), res)

    n_strict = sum(r.status == "strict_pass" for r in results)
    n_fallback = sum(r.status == "fallback" for r in results)
    n_failed = sum(r.status == "failed" for r in results)
    n_built = n_strict + n_fallback

    meta = {
        "name": name,
        "base_seed": base_seed,
        "max_retries": max_retries,
        "on_validation_fail": on_fail,
        "axes": {k: list(v) for k, v in config["axes"].items()},
        "fixed": config.get("fixed", {}),
        "n_points": len(points),
        "n_built": n_built,
        "n_strict_pass": n_strict,
        "n_fallback": n_fallback,
        "n_failed": n_failed,
    }
    csv_path, json_path = _write_manifest(results, output_root, meta)

    summary = dict(meta)
    summary.update({
        "output_directory": str(output_root),
        "manifest_csv": str(csv_path),
        "manifest_json": str(json_path),
        "results": results,
    })
    return summary
