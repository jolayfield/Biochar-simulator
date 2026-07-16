"""
Force-field table integrity — GROMACS_OPLS_TYPE_MAP, OPLS_ATOM_TYPES,
OPLS_LJ_PARAMS, and the bonded tables.

Why this file exists
--------------------
Against stock ``oplsaa.ff`` the exporter does not emit its own ``[atomtypes]``
block, so GROMACS resolves mass, LJ parameters, and bonded terms from the
*mapped* ``opls_XXX`` name. ``grompp`` accepts any name that exists in the
forcefield, so mapping a nitrogen to a carbon type is silently accepted and
produces physically meaningless dynamics -- wrong mass, wrong radius, wrong
everything, with no error anywhere.

Asserting only that a mapped name *exists* does not catch that. Asserting the
mapped name refers to the *right element* does. Five real mappings were wrong in
exactly this way (a hydrogen carrying mass 14.007, sulfur typed as carbon) and
survived because nothing checked.
"""

import re
from pathlib import Path

import pytest

from biochar.constants import (
    GROMACS_OPLS_TYPE_MAP,
    OPLS_ANGLE_PARAMS,
    OPLS_ATOM_TYPES,
    OPLS_BOND_PARAMS,
    OPLS_LJ_PARAMS,
)


# Ground truth for every opls_XXX name this package maps to, transcribed from
# stock GROMACS share/top/oplsaa.ff (atomtypes.atp + ffnonbonded.itp).
#
#   opls_name: (atomic_number, mass, description from atomtypes.atp)
#
# This table is committed rather than parsed so the check runs everywhere, not
# only on machines with GROMACS installed. test_mapped_types_match_installed_
# forcefield below cross-checks it against a real install when one is present,
# which is what keeps this table honest.
STOCK_OPLS = {
    "opls_116": (8, 15.9994, "O SPC/E water"),
    "opls_135": (6, 12.0110, "alkane CH3"),
    "opls_140": (1, 1.0080, "alkane H"),
    "opls_145": (6, 12.0110, "benzene C"),
    "opls_146": (1, 1.0080, "benzene H"),
    "opls_154": (8, 15.9994, "alcohol OH"),
    "opls_155": (1, 1.0080, "alcohol HO"),
    "opls_202": (16, 32.0600, "all-atom S: sulfides, S=C"),
    "opls_204": (1, 1.0080, "all-atom H(S): thiols"),
    "opls_209": (6, 12.0110, "all-atom C: CH3, sulfides"),
    "opls_267": (6, 12.0110, "C: C in COOH"),
    "opls_268": (8, 15.9994, "O: OH in COOH"),
    "opls_269": (8, 15.9994, "O: =O in COOH"),
    "opls_270": (1, 1.0080, "H: H in COOH"),
    "opls_271": (6, 12.0110, "C in COO- carboxylate"),
    "opls_272": (8, 15.9994, "O in COO- carboxylate"),
    "opls_278": (8, 15.9994, "O: ketone C=O"),
    "opls_287": (7, 14.0067, "N (RNH3+)"),
    "opls_290": (1, 1.0080, "H (RNH3+)"),
    "opls_379": (7, 14.0067, "CytH+ N3 protonated cytosine"),
    "opls_417": (16, 32.0600, "S in CH3S- thiolate"),
    "opls_420": (8, 15.9994, "O in CH3O- alkoxide"),
    "opls_467": (8, 15.9994, "aryl ether O"),
    "opls_513": (1, 1.0080, "H on N in HIP"),
    "opls_520": (7, 14.0067, "N in pyridine"),
    "opls_542": (7, 14.0067, "N in pyrrole"),
    "opls_545": (1, 1.0080, "H1 in pyrrole (the N-H)"),
    "opls_900": (7, 14.0067, "N primary amines"),
    "opls_909": (1, 1.0080, "H(N) primary amines"),
}

# Element implied by each internal type name. Derived from what the type *is*
# chemically, independent of what it maps to -- that independence is the whole
# point, since a typo in the mapping must not be able to justify itself.
INTERNAL_ELEMENT = {
    "CA": 6, "HA": 1, "CT": 6, "HC": 1,
    "OH": 8, "HO": 1, "OS": 8, "OC": 8, "OW": 8,
    "C": 6, "O": 8, "OH2": 8, "HO2": 1,
    "N": 7, "NT": 7, "NA": 7, "HNA": 1,
    "SH_": 16, "HSH": 1, "SS": 16,
    "NPY": 7, "NPR": 7, "NGR": 7, "HNPR": 1,
    # ionized forms
    "CM": 6, "O2M": 8, "OM": 8, "SM": 16,
    "NPYP": 7, "HPYP": 1, "NAP": 7, "HNAP": 1,
}

SYMBOL = {1: "H", 6: "C", 7: "N", 8: "O", 16: "S"}

IONIZED_TYPES = ("CM", "O2M", "OM", "SM", "NPYP", "HPYP", "NAP", "HNAP")

# Mappings known to be wrong and deliberately not fixed here.
#
# SS (thioether S) points at opls_209, which is the *carbon* of a sulfide CH3
# group. The sulfur mappings are owned by a separate piece of work auditing
# SH_/SS against thiophenol's opls_734, so fixing SS here would collide with it.
# Marked xfail(strict) so that when that work lands, this passing forces the
# entry out of this dict rather than leaving a stale exemption behind.
KNOWN_BAD_MAPPINGS = {
    "SS": "thioether S maps to a carbon type; owned by the S-mapping audit",
}


def _mapping_params():
    """Parametrise over the mapping, xfailing entries known to be broken."""
    out = []
    for internal in sorted(GROMACS_OPLS_TYPE_MAP):
        marks = []
        if internal in KNOWN_BAD_MAPPINGS:
            marks.append(
                pytest.mark.xfail(strict=True, reason=KNOWN_BAD_MAPPINGS[internal])
            )
        out.append(pytest.param(internal, marks=marks))
    return out


class TestMappingElementConsistency:
    """The check that would have caught five silent, shipped mis-mappings."""

    @pytest.mark.parametrize("internal", _mapping_params())
    def test_mapped_type_has_the_right_element(self, internal):
        opls = GROMACS_OPLS_TYPE_MAP[internal]
        assert opls in STOCK_OPLS, (
            f"{internal} -> {opls}, which is not a known stock OPLS type. "
            f"If it is genuinely new, add it to STOCK_OPLS with its real "
            f"element and mass from oplsaa.ff."
        )
        expected = INTERNAL_ELEMENT[internal]
        actual, mass, desc = STOCK_OPLS[opls]
        assert actual == expected, (
            f"{internal} is {SYMBOL[expected]} but maps to {opls} "
            f"({desc}), which is {SYMBOL.get(actual, actual)} with mass {mass}. "
            f"grompp will accept this and silently simulate the wrong element."
        )

    @pytest.mark.parametrize("internal", _mapping_params())
    def test_mapped_type_mass_agrees_with_internal_table(self, internal):
        """OPLS_ATOM_TYPES' mass must match the mapped type's real mass."""
        opls = GROMACS_OPLS_TYPE_MAP[internal]
        _, stock_mass, desc = STOCK_OPLS[opls]
        _, internal_mass, _ = OPLS_ATOM_TYPES[internal]
        assert internal_mass == pytest.approx(stock_mass, abs=0.01), (
            f"{internal} carries mass {internal_mass} but maps to {opls} "
            f"({desc}) whose real mass is {stock_mass}"
        )

    def test_every_internal_type_has_a_declared_element(self):
        """Guard: a new type must be added to INTERNAL_ELEMENT, not skipped."""
        missing = set(GROMACS_OPLS_TYPE_MAP) - set(INTERNAL_ELEMENT)
        assert not missing, f"types with no declared element: {sorted(missing)}"

    def test_hydrogens_never_map_to_heavy_atoms(self):
        """The most egregious form of the bug, called out on its own."""
        for internal, opls in GROMACS_OPLS_TYPE_MAP.items():
            if INTERNAL_ELEMENT[internal] != 1:
                continue
            z, mass, desc = STOCK_OPLS[opls]
            assert z == 1 and mass < 2.0, (
                f"hydrogen {internal} maps to {opls} ({desc}) with mass {mass}"
            )


class TestIonizedTypesComplete:
    @pytest.mark.parametrize("t", IONIZED_TYPES)
    def test_ionized_type_is_fully_defined(self, t):
        assert t in OPLS_ATOM_TYPES, f"{t} missing from OPLS_ATOM_TYPES"
        assert t in OPLS_LJ_PARAMS, f"{t} missing from OPLS_LJ_PARAMS"
        assert t in GROMACS_OPLS_TYPE_MAP, f"{t} missing from GROMACS_OPLS_TYPE_MAP"

    def test_every_atom_type_has_lj_params(self):
        missing = set(OPLS_ATOM_TYPES) - set(OPLS_LJ_PARAMS)
        assert not missing, f"atom types with no LJ parameters: {sorted(missing)}"

    def test_bonded_tables_only_reference_known_types(self):
        for key in OPLS_BOND_PARAMS:
            for t in key:
                assert t in OPLS_ATOM_TYPES, f"bond {key} references unknown type {t}"
        for key in OPLS_ANGLE_PARAMS:
            for t in key:
                assert t in OPLS_ATOM_TYPES, f"angle {key} references unknown type {t}"

    def test_mapped_names_are_well_formed(self):
        for internal, opls in GROMACS_OPLS_TYPE_MAP.items():
            assert re.fullmatch(r"opls_\d+[A-Z]?", opls), (
                f"{internal} -> {opls!r} is not a valid opls name"
            )


class TestIonizedChargeDirection:
    """Ionized types must carry charge in the direction their sign implies."""

    def test_anionic_oxygens_are_more_negative_than_neutral(self):
        _, _, phenolate = OPLS_ATOM_TYPES["OM"]
        _, _, phenol = OPLS_ATOM_TYPES["OH"]
        assert phenolate < phenol, "phenolate O is not more negative than phenol O"

        _, _, carboxylate = OPLS_ATOM_TYPES["O2M"]
        _, _, carboxyl_oh = OPLS_ATOM_TYPES["OH2"]
        assert carboxylate < carboxyl_oh

    def test_anionic_sulfur_is_more_negative_than_neutral(self):
        _, _, thiolate = OPLS_ATOM_TYPES["SM"]
        _, _, thiol = OPLS_ATOM_TYPES["SH_"]
        assert thiolate < thiol

    def test_cationic_nitrogen_hydrogens_are_positive(self):
        """
        The N itself stays negative in OPLS (the + lives on the H atoms), so the
        cation's signature is its hydrogens, not its nitrogen.
        """
        for h_type in ("HPYP", "HNAP"):
            _, _, charge = OPLS_ATOM_TYPES[h_type]
            assert charge > 0, f"{h_type} on a cationic N is not positive"

    def test_cationic_nitrogen_h_is_more_positive_than_neutral_amine_h(self):
        _, _, anilinium_h = OPLS_ATOM_TYPES["HNAP"]
        _, _, aniline_h = OPLS_ATOM_TYPES["HNA"]
        assert anilinium_h > aniline_h


class TestAgainstInstalledForcefield:
    """
    Cross-check STOCK_OPLS against a real oplsaa.ff when one is installed.

    Skips cleanly when GROMACS is absent, mirroring how test_qm_charges.py
    handles a missing MOPAC binary. The committed STOCK_OPLS table above is the
    always-on check; this one keeps that table honest.
    """

    @staticmethod
    def _find_ffnonbonded():
        import os
        import shutil

        roots = []
        gmx = shutil.which("gmx") or shutil.which("gmx_mpi")
        if gmx:
            # <prefix>/bin/gmx -> <prefix>/share/gromacs/top
            roots.append(Path(gmx).resolve().parent.parent / "share/gromacs/top")
        if os.environ.get("GMXDATA"):
            roots.append(Path(os.environ["GMXDATA"]) / "top")
        for root in roots:
            candidate = root / "oplsaa.ff/ffnonbonded.itp"
            if candidate.is_file():
                return candidate
        return None

    def test_mapped_types_match_installed_forcefield(self):
        path = self._find_ffnonbonded()
        if path is None:
            pytest.skip("GROMACS oplsaa.ff not found; STOCK_OPLS is the fallback check")

        real = {}
        for line in path.read_text().splitlines():
            m = re.match(r"\s*(opls_\w+)\s+(\S+)\s+(\d+)\s+([\d.]+)", line)
            if m:
                real[m.group(1)] = (int(m.group(3)), float(m.group(4)))

        for internal, opls in sorted(GROMACS_OPLS_TYPE_MAP.items()):
            assert opls in real, f"{internal} -> {opls} absent from {path}"
            z, mass = real[opls]
            expected = INTERNAL_ELEMENT[internal]
            assert z == expected, (
                f"{internal} is {SYMBOL[expected]} but {opls} is "
                f"{SYMBOL.get(z, z)} in the installed forcefield at {path}"
            )
            # And the committed table must agree with the install.
            assert STOCK_OPLS[opls][0] == z, f"STOCK_OPLS is stale for {opls}"
            assert STOCK_OPLS[opls][1] == pytest.approx(mass, abs=0.01), (
                f"STOCK_OPLS mass is stale for {opls}"
            )
