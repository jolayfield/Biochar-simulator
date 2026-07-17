"""Tests that GROMACS_OPLS_TYPE_MAP agrees with a real oplsaa.ff.

When an exported .itp ``#include``s a real oplsaa.ff, only the ``opls_XXX`` name
reaches GROMACS -- mass, charge, LJ and bonded parameters are all looked up from
the forcefield by that name. A mapping that names an existing but chemically
wrong type therefore produces a topology that runs happily and simulates the
wrong molecule, so these tests check the mapped types against the installed
forcefield rather than against a hand-copied table.

The forcefield-backed tests run only when an oplsaa.ff can be located (set
``BIOCHAR_OPLSAA_FF`` to point at one directly, or have ``gmx`` on PATH /
``GMXDATA`` set). The regression tests below need no forcefield.
"""

from __future__ import annotations

import inspect
import os
import re
import shutil
from itertools import combinations
from pathlib import Path

import pytest

from biochar.constants import (
    FUNCTIONAL_GROUPS,
    GROMACS_OPLS_TYPE_MAP,
    OPLS_ATOM_TYPES,
    SUPPLEMENTARY_ANGLE_PARAMS,
)


def _find_oplsaa() -> Path | None:
    """Locate a stock oplsaa.ff directory, or return None if GROMACS is absent."""
    override = os.environ.get("BIOCHAR_OPLSAA_FF")
    if override:
        path = Path(override)
        return path if (path / "atomtypes.atp").is_file() else None

    candidates: list[Path] = []

    gmxdata = os.environ.get("GMXDATA")
    if gmxdata:
        candidates.append(Path(gmxdata) / "top" / "oplsaa.ff")

    gmxlib = os.environ.get("GMXLIB")
    if gmxlib:
        candidates.append(Path(gmxlib) / "oplsaa.ff")
        candidates.append(Path(gmxlib))

    for exe in ("gmx", "gmx_mpi"):
        found = shutil.which(exe)
        if not found:
            continue
        # <prefix>/bin/gmx -> <prefix>/share/gromacs/top/oplsaa.ff
        prefix = Path(found).resolve().parent.parent
        candidates.append(prefix / "share" / "gromacs" / "top" / "oplsaa.ff")

    for path in candidates:
        if (path / "atomtypes.atp").is_file():
            return path
    return None


OPLSAA = _find_oplsaa()
requires_oplsaa = pytest.mark.skipif(
    OPLSAA is None,
    reason="no oplsaa.ff found (install GROMACS or set BIOCHAR_OPLSAA_FF)",
)


def _parse_atomtypes(ff: Path) -> dict[str, float]:
    """Map opls_XXX -> mass from atomtypes.atp."""
    masses: dict[str, float] = {}
    for line in (ff / "atomtypes.atp").read_text().splitlines():
        match = re.match(r"\s*(opls_\w+)\s+([0-9.]+)", line)
        if match:
            masses[match.group(1)] = float(match.group(2))
    return masses


def _parse_atomtype_descriptions(ff: Path) -> dict[str, str]:
    """Map opls_XXX -> the trailing comment describing the type."""
    descriptions: dict[str, str] = {}
    for line in (ff / "atomtypes.atp").read_text().splitlines():
        match = re.match(r"\s*(opls_\w+)\s+[0-9.]+\s*;\s*(.*)", line)
        if match:
            descriptions[match.group(1)] = match.group(2).strip()
    return descriptions


def _parse_bonded_types(ff: Path) -> dict[str, str]:
    """Map opls_XXX -> bonded type (column 2 of ffnonbonded.itp).

    ffbonded.itp is keyed by this coarser name, not by the opls_XXX name --
    several opls types share one bonded type.
    """
    bonded: dict[str, str] = {}
    for line in (ff / "ffnonbonded.itp").read_text().splitlines():
        parts = line.split(";", 1)[0].split()
        if len(parts) >= 2 and parts[0].startswith("opls_"):
            bonded[parts[0]] = parts[1]
    return bonded


def _parse_ffbonded_section(ff: Path, header: str, arity: int) -> set:
    """Read [ bondtypes ] (`i j func b0 kb`) or [ angletypes ] (`i j k func th0 cth`).

    Bonds are order-free; angles fix the middle atom and are canonicalised on the
    outer two.
    """
    found, inside = set(), False
    for raw in (ff / "ffbonded.itp").read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            inside = line.startswith(header)
            continue
        if not inside or not line:
            continue
        parts = line.split()
        if len(parts) >= arity:
            if arity == 2:
                found.add(tuple(sorted(parts[:2])))
            else:
                i, j, k = parts[0], parts[1], parts[2]
                found.add((min(i, k), j, max(i, k)))
    return found


@requires_oplsaa
class TestTypesExist:
    def test_every_mapped_type_exists_in_forcefield(self):
        known = _parse_atomtypes(OPLSAA)
        missing = {
            internal: opls
            for internal, opls in GROMACS_OPLS_TYPE_MAP.items()
            if opls not in known
        }
        assert not missing, (
            f"mapped opls types absent from {OPLSAA}/atomtypes.atp: {missing}"
        )


@requires_oplsaa
class TestTypesAreTheRightElement:
    """Existence is not enough -- every wrong mapping found so far named a real
    type of the wrong element (e.g. a nitrogen typed as an aromatic carbon)."""

    def test_mapped_type_mass_matches_internal_type_mass(self):
        known = _parse_atomtypes(OPLSAA)
        descriptions = _parse_atomtype_descriptions(OPLSAA)

        mismatches = []
        for internal, opls in GROMACS_OPLS_TYPE_MAP.items():
            if opls not in known:
                continue  # covered by TestTypesExist
            expected = OPLS_ATOM_TYPES[internal][1]
            actual = known[opls]
            if abs(expected - actual) > 0.05:
                mismatches.append(
                    f"{internal} -> {opls}: internal mass {expected} but "
                    f"{opls} is {actual} ({descriptions.get(opls, '?')})"
                )

        assert not mismatches, "wrong-element OPLS mappings:\n" + "\n".join(mismatches)


class TestSulfurRegression:
    """opls_202 (sulfide/S=C S) and opls_209 (a *carbon*) were previously used for
    the thiol and thioether sulfurs. These pin the corrected types."""

    def test_thiol_sulfur_is_thiophenol_type(self):
        assert GROMACS_OPLS_TYPE_MAP["SH_"] == "opls_734"

    def test_thiol_hydrogen_is_thiol_h(self):
        assert GROMACS_OPLS_TYPE_MAP["HSH"] == "opls_204"

    def test_thioether_sulfur_is_aryl_sulfide_type(self):
        assert GROMACS_OPLS_TYPE_MAP["SS"] == "opls_222"

    @requires_oplsaa
    def test_sulfur_types_are_actually_sulfur(self):
        known = _parse_atomtypes(OPLSAA)
        for internal in ("SH_", "SS"):
            assert known[GROMACS_OPLS_TYPE_MAP[internal]] == pytest.approx(32.06)


class TestNitrogenRegression:
    """The N-doping and aniline types had the same defect: HNA/NPY/NPR/HNPR each
    named a real type of the wrong element."""

    def test_aniline_hydrogen_is_a_hydrogen_type(self):
        assert GROMACS_OPLS_TYPE_MAP["HNA"] == "opls_909"

    def test_pyridinic_n_is_pyridine_nitrogen(self):
        assert GROMACS_OPLS_TYPE_MAP["NPY"] == "opls_520"

    def test_pyrrolic_n_is_pyrrole_nitrogen(self):
        assert GROMACS_OPLS_TYPE_MAP["NPR"] == "opls_542"

    def test_pyrrolic_hydrogen_is_pyrrole_nh_hydrogen(self):
        assert GROMACS_OPLS_TYPE_MAP["HNPR"] == "opls_545"


# ---------------------------------------------------------------------------
# Depth 3 -- bonded resolution
# ---------------------------------------------------------------------------
#
# The checks above are depth 1 (the type exists) and depth 2 (it is the right
# element). Both passed while thioether topologies were still unusable: SS ->
# opls_222 is genuinely sulfur with the right mass, but a thioether emits a
# CA-S-CA angle and stock OPLS defines no such angletype, so grompp refused the
# topology with "No default Angle types".
#
# [ bonds ] and [ angles ] are emitted with funct only, so the forcefield is the
# sole source of bonded parameters -- except where SUPPLEMENTARY_ANGLE_PARAMS
# writes them inline. Every combination the generator can emit must resolve
# through one of those two routes.


def _config_kwargs(**extra):
    """Build a GeneratorConfig kwargs dict, tolerating signature drift."""
    from biochar.biochar_generator import GeneratorConfig

    params = inspect.signature(GeneratorConfig).parameters
    kwargs = {"target_num_carbons": 20, "strict": False, "seed": 1}
    kwargs.update(extra)
    unknown = set(kwargs) - set(params)
    assert not unknown, f"GeneratorConfig has no such parameter(s): {unknown}"
    return kwargs


def _unresolvable(mol, atom_types, bonded, bondtypes, angletypes):
    """Bonds/angles the molecule emits that nothing can parameterise."""

    def to_bonded(internal):
        return bonded[GROMACS_OPLS_TYPE_MAP[internal]]

    missing = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        pair = tuple(sorted([to_bonded(atom_types[i]), to_bonded(atom_types[j])]))
        if pair not in bondtypes:
            missing.append(f"bond {atom_types[i]}-{atom_types[j]} ({'-'.join(pair)})")

    for atom in mol.GetAtoms():
        j = atom.GetIdx()
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        for i, k in combinations(neighbours, 2):
            outer = sorted([atom_types[i], atom_types[k]])
            internal = (outer[0], atom_types[j], outer[1])
            bi, bk = to_bonded(atom_types[i]), to_bonded(atom_types[k])
            triple = (min(bi, bk), to_bonded(atom_types[j]), max(bi, bk))
            if triple in angletypes or internal in SUPPLEMENTARY_ANGLE_PARAMS:
                continue
            missing.append(f"angle {'-'.join(internal)} ({'-'.join(triple)})")
    return sorted(set(missing))


# Gaps this check found that are real but out of scope to fix here. Each is
# xfail(strict) so that fixing one fails the suite until its marker is removed --
# a silent pass would let the gap reopen unnoticed. carboxyl was the other one,
# and the pH work's typing fix retired it exactly that way.
_PYRIDINIC_XFAIL = (
    "A hydroxyl on a ring carbon adjacent to a pyridinic N emits NC-CA-OH, which "
    "stock OPLS does not define (3-hydroxypyridine is plausible chemistry OPLS "
    "just omits). Reachable with default settings, since the default O/C ratio "
    "adds phenolic OH. A genuine forcefield gap, so it belongs in "
    "SUPPLEMENTARY_ANGLE_PARAMS -- but only with a defensible value and "
    "provenance, which needs its own change."
)


def _group_params():
    """Functional-group cases, xfail-marked where a known gap makes them fail.

    carboxyl used to be marked here: AtomTyper never assigned the carboxylic-acid
    types, so a -COOH emitted O_2-C-OH instead of the real O_3-C-OH. The pH work
    fixed the typing, the marker started passing, and xfail(strict) failed the
    suite until it was removed.
    """
    known: dict[str, str] = {}
    params = []
    for group in sorted(FUNCTIONAL_GROUPS):
        reason = known.get(group)
        marks = [pytest.mark.xfail(strict=True, reason=reason)] if reason else []
        params.append(pytest.param(group, marks=marks))
    return params


@requires_oplsaa
class TestBondedResolution:
    """Depth 3: every bond and angle the generator emits must be parameterisable."""

    @pytest.mark.parametrize("group", _group_params())
    def test_functional_group_emits_only_resolvable_terms(self, group):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        gen = BiocharGenerator(
            GeneratorConfig(**_config_kwargs(functional_groups={group: 2}))
        )
        mol, _, _ = gen.generate()
        missing = _unresolvable(
            mol,
            gen.atom_types,
            _parse_bonded_types(OPLSAA),
            _parse_ffbonded_section(OPLSAA, "[ bondtypes", 2),
            _parse_ffbonded_section(OPLSAA, "[ angletypes", 3),
        )
        assert not missing, (
            f"'{group}' emits terms no forcefield entry and no "
            f"SUPPLEMENTARY_ANGLE_PARAMS covers: {missing}"
        )

    @pytest.mark.parametrize(
        "knob",
        [
            pytest.param(
                "num_pyridinic",
                marks=pytest.mark.xfail(strict=True, reason=_PYRIDINIC_XFAIL),
            ),
            "num_pyrrolic",
            "num_graphitic",
        ],
    )
    def test_ring_nitrogen_emits_only_resolvable_terms(self, knob):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        gen = BiocharGenerator(GeneratorConfig(**_config_kwargs(**{knob: 2})))
        mol, _, _ = gen.generate()
        missing = _unresolvable(
            mol,
            gen.atom_types,
            _parse_bonded_types(OPLSAA),
            _parse_ffbonded_section(OPLSAA, "[ bondtypes", 2),
            _parse_ffbonded_section(OPLSAA, "[ angletypes", 3),
        )
        assert not missing, f"'{knob}' emits unresolvable terms: {missing}"


@requires_oplsaa
class TestSupplementDoesNotShadowForcefield:
    """The supplement may only hold what oplsaa.ff lacks.

    A value that also exists in the forcefield is duplication, and duplication is
    what drifted the deleted OPLS_LJ_PARAMS/OPLS_BOND_PARAMS tables out of sync in
    the first place. A value that exists nowhere else cannot drift.
    """

    def test_no_supplementary_angle_duplicates_a_stock_angletype(self):
        bonded = _parse_bonded_types(OPLSAA)
        angletypes = _parse_ffbonded_section(OPLSAA, "[ angletypes", 3)

        shadowed = []
        for (i, j, k) in SUPPLEMENTARY_ANGLE_PARAMS:
            bi, bk = (bonded[GROMACS_OPLS_TYPE_MAP[i]], bonded[GROMACS_OPLS_TYPE_MAP[k]])
            triple = (min(bi, bk), bonded[GROMACS_OPLS_TYPE_MAP[j]], max(bi, bk))
            if triple in angletypes:
                shadowed.append(f"{i}-{j}-{k} ({'-'.join(triple)})")

        assert not shadowed, (
            "SUPPLEMENTARY_ANGLE_PARAMS shadows angletypes oplsaa.ff already "
            f"defines -- remove them and let the forcefield supply them: {shadowed}"
        )
