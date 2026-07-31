"""
Tests for the rqm/gromacs-export.md scenarios -- the export contract itself,
as distinct from tests/test_gromacs_export.py's coverage of file structure and
section presence.

The .itp carries no parameters: every physical constant is resolved by grompp
from the opls_XXX names and the bonded terms written here. These tests therefore
check what the topology *contains*, not only that it parses -- an omitted
interaction produces a topology grompp accepts and a run that reports plausible
numbers while integrating different physics.

.gro is fixed-width. Columns, per the format definition:
    [0:5] residue number, [5:10] residue name, [10:15] atom name,
    [15:20] atom number, [20:28] x, [28:36] y, [36:44] z
"""

import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from biochar.export.gromacs_export import GROFileWriter, ITPFileWriter, GromacsExporter
from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig


def _export(**cfg):
    """Generate and export one molecule; return (mol, coords, charges, texts)."""
    kwargs = {"target_num_carbons": 20, "strict": False, "seed": 1}
    kwargs.update(cfg)
    gen = BiocharGenerator(GeneratorConfig(**kwargs))
    mol, coords, _ = gen.generate()
    out = Path(tempfile.mkdtemp())
    gro_p, top_p, itp_p = GromacsExporter(str(out)).export(
        mol, coords, gen.atom_types, gen.charges, molecule_name="BCX"
    )
    return mol, coords, gen.charges, {
        "gro": gro_p.read_text(),
        "top": top_p.read_text(),
        "itp": itp_p.read_text(),
    }


def _gro_atom_lines(gro_text):
    """The atom records: everything between the count line and the box line."""
    return gro_text.splitlines()[2:-1]


def _gro_positions(gro_text):
    return np.array(
        [
            [float(ln[20:28]), float(ln[28:36]), float(ln[36:44])]
            for ln in _gro_atom_lines(gro_text)
        ]
    )


def _gro_atom_names(gro_text):
    return [ln[10:15].strip() for ln in _gro_atom_lines(gro_text)]


def _itp_section(itp_text, header):
    """Lines of one section, comments and blanks dropped."""
    out, inside = [], False
    for raw in itp_text.splitlines():
        if raw.startswith("["):
            inside = raw.startswith(header)
            continue
        line = raw.strip()
        if inside and line and not line.startswith(";"):
            out.append(line)
    return out


def _itp_atom_names(itp_text):
    return [ln.split()[4] for ln in _itp_section(itp_text, "[ atoms ") if len(ln.split()) >= 5]


def _pairs_1_4(mol):
    """Atom index pairs exactly three bonds apart -- the 1-4 set."""
    dm = Chem.GetDistanceMatrix(mol)
    n = mol.GetNumAtoms()
    return {(i, j) for i in range(n) for j in range(i + 1, n) if dm[i][j] == 3}


def _big_carbon_mol(n):
    """n bare carbons. Enough for the writers, which read symbols and indices."""
    rw = Chem.RWMol()
    for _ in range(n):
        rw.AddAtom(Chem.Atom("C"))
    coords = np.zeros((n, 3))
    coords[:, 0] = np.arange(n) * 1.5
    return rw.GetMol(), coords


class TestCoordinateUnits:
    # rq-6024903d
    def test_every_interatomic_distance_scales_by_one_tenth(self):
        mol, coords, _, text = _export(O_C_ratio=0.15)
        written = _gro_positions(text["gro"])
        assert len(written) == mol.GetNumAtoms()

        def pdist(a):
            d = a[:, None, :] - a[None, :, :]
            return np.sqrt((d**2).sum(-1))

        iu = np.triu_indices(len(written), 1)
        d_nm, d_ang = pdist(written)[iu], pdist(coords)[iu]

        # Distances, not absolute positions: export translates the molecule into
        # its box, so positions are shifted but every distance is preserved.
        # atol covers the .gro's three-decimal write precision (0.001 nm).
        assert np.allclose(d_nm, 0.1 * d_ang, atol=0.003), (
            f"max deviation {np.abs(d_nm - 0.1 * d_ang).max():.5f} nm -- "
            f"coordinates are not in nanometres"
        )
        # A factor of ten cannot hide inside the write precision.
        assert d_nm.max() < 0.5 * d_ang.max(), "distances look like Angstrom, not nm"


class TestAtomNames:
    """Names are built as symbol + index, and the .gro field is five characters.

    Both scenarios need a molecule whose names exceed that width, which starts at
    the ten-thousandth atom of a single-letter element.
    """

    SIZE = 10050

    # rq-07b27c79
    def test_atom_names_agree_between_structure_and_topology(self):
        mol, coords = _big_carbon_mol(self.SIZE)
        out = Path(tempfile.mkdtemp())
        gro_p, itp_p = out / "big.gro", out / "big.itp"
        GROFileWriter.write(str(gro_p), mol, coords, molecule_name="BCX", title="t")
        ITPFileWriter.write(
            str(itp_p),
            mol,
            {i: "CA" for i in range(self.SIZE)},
            {i: 0.0 for i in range(self.SIZE)},
            molecule_name="BCX",
        )

        gro_names = _gro_atom_names(gro_p.read_text())
        itp_names = _itp_atom_names(itp_p.read_text())
        assert len(gro_names) == len(itp_names) == self.SIZE

        mismatched = [(a, b) for a, b in zip(gro_names, itp_names) if a != b]
        assert not mismatched, (
            f"{len(mismatched)} atom names differ between the .gro and the .itp, "
            f"first: {mismatched[:3]}"
        )

    # rq-b2e50e9d
    def test_no_two_atoms_share_a_name(self):
        mol, coords = _big_carbon_mol(self.SIZE)
        path = Path(tempfile.mkdtemp()) / "big.gro"
        GROFileWriter.write(str(path), mol, coords, molecule_name="BCX", title="t")

        names = _gro_atom_names(path.read_text())
        duplicated = {n for n, c in Counter(names).items() if c > 1}
        assert not duplicated, (
            f"{len(duplicated)} atom names identify more than one atom, "
            f"e.g. {sorted(duplicated)[:3]}"
        )


class TestExclusionsAndPairs:
    # rq-a22b8e8d
    def test_exclusion_count_is_three(self):
        _, _, _, text = _export(O_C_ratio=0.15)
        line = _itp_section(text["itp"], "[ moleculetype ")[0]
        assert line.split()[1] == "3", (
            f"moleculetype declares nrexcl {line.split()[1]}; OPLS-AA expects 3"
        )

    # rq-fb164878
    def test_every_1_4_pair_is_listed(self):
        mol, _, _, text = _export(O_C_ratio=0.15)
        expected = _pairs_1_4(mol)
        assert expected, "fixture should contain 1-4 pairs"

        assert "[ pairs ]" in text["itp"], "topology has no [ pairs ] section"
        listed = {
            tuple(sorted((int(p[0]) - 1, int(p[1]) - 1)))
            for p in (ln.split() for ln in _itp_section(text["itp"], "[ pairs "))
            if len(p) >= 2
        }
        assert listed == expected, (
            f"{len(expected - listed)} of {len(expected)} 1-4 pairs are missing"
        )


class TestImpropers:
    """OPLS-AA impropers are funct 1 carrying the improper_Z_CA_X_Y macro.

    Not funct 2 or 4, which is the natural guess. aminoacids.rtp's
    [ bondedtypes ] row declares improper dihedral type 1, and ffbonded.itp
    supplies the constants as a #define rather than a [ dihedraltypes ] entry --
    "implemented as proper dihedrals [1+cos(2*x+180)] ... to keep things
    compatible". A check for funct 2 or 4 would therefore never pass, however
    correct the topology was.
    """

    # rq-2f3c2379
    def test_aromatic_ring_carbons_carry_an_improper(self):
        mol, _, _, text = _export(O_C_ratio=0.15)
        substituted = [
            a.GetIdx()
            for a in mol.GetAtoms()
            if a.GetSymbol() == "C"
            and a.GetIsAromatic()
            and a.IsInRing()
            and a.GetDegree() == 3
        ]
        assert substituted, "fixture should contain three-connected aromatic carbons"

        impropers = [
            p
            for p in (ln.split() for ln in _itp_section(text["itp"], "[ dihedrals "))
            if len(p) >= 6 and p[4] == "1" and p[5].startswith("improper_")
        ]
        assert len(impropers) == len(substituted), (
            f"{len(substituted)} three-connected aromatic carbons but "
            f"{len(impropers)} improper terms -- rings can pyramidalise"
        )

        # The central atom sits third, per oplsaa.ff's own residue templates.
        centres = {int(p[2]) - 1 for p in impropers}
        assert centres == set(substituted), (
            "improper central atoms do not match the aromatic carbons"
        )


class TestChargeIsReported:
    # rq-78e7477d
    def test_topology_states_total_charge_and_neutralisation_step(self):
        mol, _, _, text = _export(O_C_ratio=0.2, pH=11.0)
        formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        assert formal != 0, "fixture should be charged; raise the pH if this fails"

        top = text["top"]
        assert "total charge" in top.lower(), "topology does not state its total charge"
        assert "genion" in top, "topology does not name the neutralisation step"

    # rq-f256712c
    def test_partial_charges_sum_to_the_formal_charge(self):
        mol, _, charges, _ = _export(O_C_ratio=0.2, pH=11.0)
        formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        assert formal != 0, "fixture should be charged; raise the pH if this fails"

        assert sum(charges.values()) == pytest.approx(float(formal), abs=1e-4), (
            f"partial charges sum to {sum(charges.values()):.4f} but the formal "
            f"charge is {formal}"
        )


class TestResidueNameFitsTheField:
    # rq-5340daab
    def test_long_residue_name_is_truncated_and_reported(self, capsys):
        mol, coords = _big_carbon_mol(3)
        path = Path(tempfile.mkdtemp()) / "n.gro"
        GROFileWriter.write(
            str(path), mol, coords, molecule_name="TOOLONGNAME", title="t"
        )

        names = {ln[5:10].strip() for ln in _gro_atom_lines(path.read_text())}
        assert names == {"TOOLO"}, f"residue names written: {names}"

        # The writer reports truncation on stdout rather than via warnings.warn.
        assert "TOOLONGNAME" in capsys.readouterr().out, (
            "truncation happened silently"
        )
