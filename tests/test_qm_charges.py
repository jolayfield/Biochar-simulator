"""Tests for the QM (1.14*CM1A) charge backend.

The pure-Python CM1A mapping and MOPAC-output parsers are tested with canned
data and need no external binary. An end-to-end test runs only when a ``mopac``
executable is on PATH (``conda install -c conda-forge mopac``).
"""

import shutil

import numpy as np
import pytest

from biochar.qm_charges import (
    QMChargeAssigner,
    QMChargeError,
    cm1a_from_am1,
    parse_bond_orders,
    parse_net_atomic_charges,
    scale_and_neutralize,
)

MOPAC = shutil.which("mopac")
requires_mopac = pytest.mark.skipif(MOPAC is None, reason="mopac not on PATH")


# Canned MOPAC AM1 output for methane (trimmed to the parsed sections).
METHANE_OUT = """\
              NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS

  ATOM NO.   TYPE          CHARGE      No. of ELECS.   s-Pop       p-Pop
    1          C          -0.265908        4.2659     1.23008     3.03583
    2          H           0.066477        0.9335     0.93352
    3          H           0.066477        0.9335     0.93352
    4          H           0.066477        0.9335     0.93352
    5          H           0.066477        0.9335     0.93352
 DIPOLE           X         Y         Z       TOTAL

            (VALENCIES)   BOND ORDERS

     1  C     (3.947)     2  H 0.987     3  H 0.987     4  H 0.987     5  H 0.987

     2  H     (0.996)     1  C 0.987

     3  H     (0.996)     1  C 0.987

     4  H     (0.996)     1  C 0.987

     5  H     (0.996)     1  C 0.987


          MULLIKEN POPULATION ANALYSIS
"""


class TestParsers:
    def test_net_atomic_charges(self):
        q = parse_net_atomic_charges(METHANE_OUT)
        assert len(q) == 5
        assert q[0] == pytest.approx(-0.265908)
        assert q[1] == pytest.approx(0.066477)

    def test_bond_orders(self):
        b = parse_bond_orders(METHANE_OUT)
        # four C-H bonds, no spurious H-H entries
        assert len(b) == 4
        for j in (1, 2, 3, 4):
            assert b[(0, j)] == pytest.approx(0.987)

    def test_missing_section_raises(self):
        with pytest.raises(QMChargeError):
            parse_net_atomic_charges("no charges here")
        with pytest.raises(QMChargeError):
            parse_bond_orders("no bonds here")


class TestCM1AMapping:
    def test_methane_is_unchanged(self):
        # Methane has no CM1A corrections (C and H-on-C parameters are all zero),
        # so CM1A == AM1 net charges.
        q = parse_net_atomic_charges(METHANE_OUT)
        b = parse_bond_orders(METHANE_OUT)
        cm = cm1a_from_am1(q, b, [6, 1, 1, 1, 1])
        assert cm == pytest.approx(q)

    def test_charge_is_conserved(self):
        # Arbitrary inputs: CM1A must conserve total charge exactly.
        qm = [-0.39, 0.20, 0.19]
        bonds = {(0, 1): 0.95, (0, 2): 0.90}
        cm = cm1a_from_am1(qm, bonds, [8, 1, 1])
        assert sum(cm) == pytest.approx(sum(qm), abs=1e-9)

    def test_water_matches_hand_trace(self):
        # AM1 net charges + bond order from a real MOPAC water run.
        qm = [-0.385400, 0.192700, 0.192700]
        bonds = {(0, 1): 0.963, (0, 2): 0.963}
        cm = cm1a_from_am1(qm, bonds, [8, 1, 1])
        # Hand-traced CM1A (see PR notes): O = -0.7083, H = +0.3542.
        assert cm[0] == pytest.approx(-0.7083, abs=1e-3)
        assert cm[1] == pytest.approx(0.3542, abs=1e-3)

    def test_scale_and_neutralize(self):
        cm = [-0.7083, 0.3542, 0.3542]
        scaled = scale_and_neutralize(cm, total_charge=0.0)
        assert sum(scaled) == pytest.approx(0.0, abs=1e-12)
        # 1.14 factor applied (within the tiny neutralization correction)
        assert scaled[0] == pytest.approx(1.14 * -0.7083, abs=1e-3)


class TestDriver:
    def test_missing_binary_raises(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("O"))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        coords = mol.GetConformer().GetPositions()
        assigner = QMChargeAssigner(mopac_cmd="definitely-not-a-real-mopac-binary")
        with pytest.raises(QMChargeError, match="not found"):
            assigner.assign(mol, coords)

    @requires_mopac
    def test_water_end_to_end(self):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("O"))
        AllChem.EmbedMolecule(mol, randomSeed=1)
        AllChem.MMFFOptimizeMolecule(mol)
        coords = mol.GetConformer().GetPositions()
        charges = QMChargeAssigner().assign(mol, coords)
        assert abs(sum(charges.values())) < 1e-6
        # O strongly negative, H positive; CM1A water is quite polar.
        o = next(charges[i] for i in range(mol.GetNumAtoms())
                 if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
        assert -0.85 < o < -0.65


class TestConfigWiring:
    def test_qm_accepted(self):
        from biochar.biochar_generator import GeneratorConfig

        cfg = GeneratorConfig(charge_method="qm", target_num_carbons=10)
        assert cfg.charge_method == "qm"

    def test_invalid_method_rejected(self):
        from biochar.biochar_generator import GeneratorConfig

        with pytest.raises(ValueError, match="charge_method"):
            GeneratorConfig(charge_method="bogus")
