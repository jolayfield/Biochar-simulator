"""
Tests for the rqm/ph-protonation.md scenarios that no existing test defends.

The rest of that document is back-annotated onto tests/test_protonation.py and
tests/test_ph_integration.py, which already cover the Henderson-Hasselbalch
direction, the sampling, the census and the untouched groups. What is left here
is the part a reader of the *output* depends on: whether a structure says which
pH it was titrated at, and whether a caller is told when the sample they just
took was close to a coin flip.

Both matter for the same reason. A protonation state is one draw from an
ensemble, so a structure is only interpretable next to the pH and the seed that
produced it -- and a draw taken inside a transition band is the one most likely
to be mistaken for the ensemble mean.

Three scenarios carry xfail(strict=True); each names the defect it defers.
The other two are their complements, and pass both before and after.
"""

import warnings
from pathlib import Path

import pytest
from rdkit import Chem

from biochar.constants import PROTONATION_STATES
from biochar.pipeline.biochar_generator import BiocharGenerator, GeneratorConfig
from biochar.pipeline.protonation import ProtonationAssigner

# Benzoic acid: one carboxyl, pKa 4.20, and nothing else titratable.
BENZOIC = "O=C(O)c1ccccc1"


def _mol(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    Chem.SanitizeMol(mol)
    return mol


def _titrate(pH, smiles=BENZOIC, seed=1):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        ProtonationAssigner(seed=seed).assign(_mol(smiles), pH=pH)
    return [str(w.message) for w in caught]


class TestATransitionBandIsReported:
    # rq-7e778b79
    @pytest.mark.xfail(
        strict=True,
        reason="nothing tells a caller their pH sits on a pKa, though the "
               "module's own docstring says one structure is then a coin flip "
               "per site and needs replicates",
    )
    def test_a_ph_on_a_pka_is_reported(self):
        pKa = PROTONATION_STATES["carboxyl"].pKa
        messages = _titrate(pKa + 0.3)

        assert any("carboxyl" in m for m in messages), (
            f"a pH {pKa + 0.3} titration of a carboxyl (pKa {pKa}) warned "
            f"about nothing; messages were {messages}"
        )
        assert any("sample" in m.lower() for m in messages), (
            f"the warning does not say the structure is one sample: {messages}"
        )

    # rq-9ee7e1df
    def test_a_ph_far_from_every_pka_is_not_reported(self):
        """The complement: passes now because nothing warns, and must keep
        passing once the warning exists -- a band check that fires everywhere
        is no more use than one that fires nowhere."""
        # Every pKa in the table is below 9.5, and the only site here is a
        # carboxyl at 4.20.
        assert not _titrate(7.0), (
            "a pH 2.8 units from the only titratable pKa warned anyway"
        )

    def test_the_only_titratable_group_here_is_the_carboxyl(self):
        """Guards the two tests above: they mean nothing if the fixture drifts."""
        assigner = ProtonationAssigner(seed=1)
        sites = assigner._find_sites(_mol(BENZOIC))
        assert [group for _, group in sites] == ["carboxyl"], sites


# --------------------------------------------------------------------------- #
# File provenance
# --------------------------------------------------------------------------- #
def _export(tmp_path, name, **cfg):
    kwargs = dict(
        target_num_carbons=36,
        functional_groups={"carboxyl": 3},
        seed=5,
        strict=False,
        H_C_tolerance=1.0,
        O_C_tolerance=1.0,
        molecule_name="PHX",
    )
    kwargs.update(cfg)
    gen = BiocharGenerator(GeneratorConfig(**kwargs))
    gen.generate()
    gro, top, itp = gen.export_gromacs(
        output_directory=str(tmp_path / name), basename=name
    )
    return gen, Path(gro).read_text(), Path(itp).read_text()


@pytest.fixture(scope="module")
def titrated(tmp_path_factory):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return _export(tmp_path_factory.mktemp("ph"), "titrated", pH=11.0)


@pytest.fixture(scope="module")
def untitrated(tmp_path_factory):
    return _export(tmp_path_factory.mktemp("neutral"), "neutral")


class TestTheFilesRecordTheirPH:
    # rq-7ed1460a
    @pytest.mark.xfail(
        strict=True,
        reason="the .gro title is a bare timestamp, so a pH 3 and a pH 11 "
               "structure are indistinguishable to a reader opening either",
    )
    def test_the_coordinate_file_names_its_ph_and_charge(self, titrated):
        gen, gro, _ = titrated
        title = gro.splitlines()[0]

        assert gen.composition.net_charge != 0, (
            "fixture is void: raise the pH until the structure carries a charge"
        )
        assert "pH 11" in title, f"the title does not state the pH: {title!r}"
        assert str(gen.composition.net_charge) in title, (
            f"the title does not state the charge "
            f"{gen.composition.net_charge}: {title!r}"
        )

    # rq-4189fd49
    @pytest.mark.xfail(
        strict=True,
        reason="the .itp header records only a timestamp, so the sample a "
               "topology represents cannot be reproduced from it",
    )
    def test_the_topology_names_its_ph_and_seed(self, titrated):
        _, _, itp = titrated
        header = itp.split("[ moleculetype ]", 1)[0]

        assert "pH 11" in header, f"the header does not state the pH:\n{header}"
        assert "seed 5" in header.lower(), (
            f"the header does not state the seed the sample came from:\n{header}"
        )

    # rq-e0932c52
    def test_an_untitrated_structure_claims_no_ph(self, untitrated):
        """Pinned so the fix cannot answer the two above by stamping a pH on
        every structure whether one was requested or not."""
        gen, gro, itp = untitrated
        assert gen.config.pH is None, "fixture is void: it requested a pH"

        title = gro.splitlines()[0]
        header = itp.split("[ moleculetype ]", 1)[0]
        assert "pH" not in title, f"an untitrated structure claims a pH: {title!r}"
        assert "pH" not in header, f"an untitrated topology claims a pH:\n{header}"
