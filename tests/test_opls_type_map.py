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

import os
import re
import shutil
from pathlib import Path

import pytest

from biochar.constants import GROMACS_OPLS_TYPE_MAP, OPLS_ATOM_TYPES


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
