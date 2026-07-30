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
    # rq-03a40f15
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

    # rq-03a40f15
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


class TestEveryAssignableTypeIsMapped:
    """GromacsExporter writes GROMACS_OPLS_TYPE_MAP.get(t, t), so an unmapped
    internal type is silently written out under its own name and grompp fails
    with "Atomtype N not found". Nothing in the export path defines these
    locally, so every type the pipeline can assign must have an entry."""

    # rq-70e5f033
    def test_no_opls_atom_type_lacks_a_gromacs_mapping(self):
        unmapped = sorted(set(OPLS_ATOM_TYPES) - set(GROMACS_OPLS_TYPE_MAP))
        assert not unmapped, (
            f"OPLS_ATOM_TYPES entries with no GROMACS_OPLS_TYPE_MAP entry: "
            f"{unmapped}. Either map them to a verified oplsaa.ff type or "
            f"remove them along with the code that assigns them."
        )

    # rq-3e169a31
    def test_map_has_no_entries_for_unknown_internal_types(self):
        extra = sorted(set(GROMACS_OPLS_TYPE_MAP) - set(OPLS_ATOM_TYPES))
        assert not extra, (
            f"GROMACS_OPLS_TYPE_MAP entries with no OPLS_ATOM_TYPES entry: {extra}"
        )


class TestNitrogenTypingIsSizeIndependent:
    """The tertiary/quaternary "N"/"NT" types were unmapped but reachable: RDKit
    kekulization fails on large fused sheets, so an aniline N whose neighbours
    all read is_aromatic=False fell past the NA branch to "N" (and its
    hydrogens to "HC"). The same Ar-NH2 therefore typed differently depending
    only on how big the sheet was."""

    # rq-4063a024
    def test_tertiary_and_quaternary_nitrogen_types_are_gone(self):
        for dead in ("N", "NT"):
            assert dead not in OPLS_ATOM_TYPES

    @pytest.mark.parametrize("num_carbons", [40, 150])
    # rq-3789a326 rq-c7715c26
    def test_amino_nitrogen_types_as_aniline_n_at_any_size(self, num_carbons):
        from biochar.biochar_generator import BiocharGenerator, GeneratorConfig

        generator = BiocharGenerator(
            GeneratorConfig(
                target_num_carbons=num_carbons,
                functional_groups={"amino": 1},
                molecule_name="NTEST",
                seed=1,
                strict=False,
            )
        )
        generator.generate()

        nitrogens = [
            atom.GetIdx()
            for atom in generator.mol.GetAtoms()
            if atom.GetAtomicNum() == 7
        ]
        assert nitrogens, "amino group placed no nitrogen"

        for idx in nitrogens:
            assert generator.atom_types[idx] == "NA"
            hydrogens = [
                n.GetIdx()
                for n in generator.mol.GetAtomWithIdx(idx).GetNeighbors()
                if n.GetAtomicNum() == 1
            ]
            assert hydrogens, "aniline N has no hydrogens"
            for h_idx in hydrogens:
                assert generator.atom_types[h_idx] == "HNA"


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
    # rq-03a40f15
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


# This check has found three gaps so far, all now closed. Each was held by an
# xfail(strict) marker until its fix landed -- at which point the marker XPASSed
# and failed the suite until removed, so a fix could never silently leave a stale
# exemption behind:
#   - carboxyl (O_2-C-OH): a typing bug, retired by the pH work's typing fix.
#   - num_pyridinic (NC-CA-OH): a genuine forcefield gap, retired by the
#     NPY-CA-OH entry in SUPPLEMENTARY_ANGLE_PARAMS (phenol-analog value).
# When the next gap appears, add its marker here the same way.


def _group_params():
    """Functional-group cases, xfail-marked where a known gap makes them fail.

    Currently none are marked. carboxyl used to be: AtomTyper never assigned the
    carboxylic-acid types, so a -COOH emitted O_2-C-OH instead of the real
    O_3-C-OH. The pH work fixed the typing, the marker started passing, and
    xfail(strict) failed the suite until it was removed.
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
    # rq-35a5407d
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
            "num_pyridinic",
            "num_pyrrolic",
            "num_graphitic",
        ],
    )
    # rq-35a5407d
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

    # rq-78f96107
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


class TestAnalogDerivedTypesAreStable:
    """NGR, SS, and the two SUPPLEMENTARY_ANGLE_PARAMS entries are nearest-
    analog choices, not QM-validated -- see rqm/opls-typing.md. NGR in
    particular was moved deliberately (2026-07-17, reviewed) away from the
    stock pyridine N it used to share with NPY; a value drifting back would be
    invisible to every other test here, since it would still resolve and still
    have the right element."""

    # rq-1548650e
    def test_analog_derived_values_match_their_recorded_provenance(self):
        assert GROMACS_OPLS_TYPE_MAP["SS"] == "opls_222"      # "S in thioanisoles"
        assert GROMACS_OPLS_TYPE_MAP["NGR"] == "opls_379"     # "CytH+ N3" -- do not restore opls_520

        assert SUPPLEMENTARY_ANGLE_PARAMS[("CA", "SS", "CA")] == (104.200, 518.816)
        assert SUPPLEMENTARY_ANGLE_PARAMS[("NPY", "CA", "OH")] == (120.000, 585.760)


class TestSupplementEntriesCarryProvenance:
    """SUPPLEMENTARY_ANGLE_PARAMS entries are hand-picked analog values, not
    machine-derivable ones, so the comment naming the analog IS the
    justification for why the entry is trustworthy at all. This inspects the
    source text directly since the provenance lives in a comment, not a field
    any runtime object carries."""

    # rq-8e4335b8
    def test_every_supplement_entry_has_a_preceding_comment(self):
        import inspect
        import re

        import biochar.constants as constants_module

        source = inspect.getsource(constants_module)
        block_match = re.search(
            r"^SUPPLEMENTARY_ANGLE_PARAMS.*?=\s*\{(.*?)\n\}",
            source, re.S | re.M,
        )
        assert block_match, "could not locate the SUPPLEMENTARY_ANGLE_PARAMS block"
        lines = block_match.group(1).split("\n")

        key_line = re.compile(r'^\s*\(\s*["\']')
        undocumented = []
        for i, line in enumerate(lines):
            if not key_line.match(line):
                continue
            j = i - 1
            has_comment = False
            while j >= 0 and lines[j].strip().startswith("#"):
                has_comment = True
                j -= 1
            if not has_comment:
                undocumented.append(line.strip())

        assert not undocumented, (
            f"SUPPLEMENTARY_ANGLE_PARAMS entries with no preceding provenance "
            f"comment: {undocumented}"
        )
        # Sanity check the parser actually found the real entries, so a
        # regex/formatting change can't silently make this pass vacuously.
        assert len(SUPPLEMENTARY_ANGLE_PARAMS) == sum(
            1 for line in lines if key_line.match(line)
        )


class TestForcefieldAbsentSkipIsVisible:
    """When no oplsaa.ff is discoverable, the forcefield-backed tests must
    skip -- loudly, by name -- rather than pass having verified nothing. This
    runs a real forcefield-backed test in a clean subprocess with every
    discovery path (BIOCHAR_OPLSAA_FF, GMXDATA, GMXLIB, `gmx` on PATH)
    removed, and inspects pytest's own report rather than trusting that the
    skipif mechanism (already covered by pytest itself) was wired up right."""

    # rq-b40b353d
    def test_a_forcefield_backed_test_skips_by_name_when_ff_is_absent(self):
        import subprocess
        import sys

        env = {
            k: v for k, v in os.environ.items()
            if k not in ("BIOCHAR_OPLSAA_FF", "GMXDATA", "GMXLIB")
        }
        # A PATH with no gmx/gmx_mpi on it, so shutil.which() genuinely fails
        # rather than happening to find a real GROMACS install on this host.
        env["PATH"] = "/usr/bin:/bin"

        result = subprocess.run(
            [
                sys.executable, "-m", "pytest",
                f"{__file__}::TestTypesExist::test_every_mapped_type_exists_in_forcefield",
                "-rs", "-q", "--no-header",
            ],
            cwd=str(Path(__file__).resolve().parent.parent),
            env=env, capture_output=True, text=True, timeout=60,
        )

        combined = result.stdout + result.stderr
        assert "1 skipped" in combined, (
            f"expected the test to skip, not pass/fail/error:\n{combined}"
        )
        assert "oplsaa.ff" in combined, (
            f"skip reason must name the missing forcefield:\n{combined}"
        )
        assert result.returncode == 0, (
            "a clean skip must not report as a suite failure"
        )
