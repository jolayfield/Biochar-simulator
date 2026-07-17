"""
Tests for pH-dependent protonation — the pKa table (U1) and the
ProtonationAssigner engine (U3).
"""

import pytest
from rdkit import Chem

from biochar.constants import (
    ACIDIC_GROUPS,
    BASIC_GROUPS,
    OPLS_ATOM_TYPES,
    PH_MAX,
    PH_MIN,
    PROTONATION_STATES,
    FUNCTIONAL_GROUPS,
    ProtonationState,
)


# Literature ranges each pKa must sit inside.  These are the *reported ranges*
# for biochar surface groups (Boehm / potentiometric titration), not the
# model-compound point values — a pKa drifting outside its range means the
# constant was mistyped, not that the chemistry is subtle.
LITERATURE_PKA_RANGES = {
    "carboxyl": (3.0, 6.0),
    "phenolic": (8.0, 11.0),
    "thiol": (6.0, 8.0),
    "amino": (4.0, 6.0),
    "pyridinic": (4.0, 6.0),
}


class TestProtonationStatesTable:
    def test_every_group_is_placeable_or_pyridinic(self):
        """
        A titratable group must be one the generator can actually build.
        'pyridinic' is the exception: NitrogenSubstitutor substitutes it into
        the ring rather than OxygenAssigner attaching it.
        """
        for group in PROTONATION_STATES:
            assert group in FUNCTIONAL_GROUPS or group == "pyridinic", (
                f"'{group}' is titratable but nothing can place it"
            )

    @pytest.mark.parametrize("group", sorted(LITERATURE_PKA_RANGES))
    def test_pka_within_literature_range(self, group):
        lo, hi = LITERATURE_PKA_RANGES[group]
        pKa = PROTONATION_STATES[group].pKa
        assert lo <= pKa <= hi, f"{group} pKa {pKa} outside literature range {lo}-{hi}"

    def test_table_covers_exactly_the_expected_groups(self):
        assert set(PROTONATION_STATES) == set(LITERATURE_PKA_RANGES)

    def test_graphitic_nitrogen_is_not_titratable(self):
        """Graphitic N is permanently +1 by construction — it must not titrate."""
        assert "graphitic" not in PROTONATION_STATES

    def test_kind_is_acidic_or_basic(self):
        for group, state in PROTONATION_STATES.items():
            assert state.kind in ("acidic", "basic"), f"{group} has kind {state.kind!r}"

    def test_acidic_and_basic_partition_the_table(self):
        assert ACIDIC_GROUPS | BASIC_GROUPS == set(PROTONATION_STATES)
        assert ACIDIC_GROUPS & BASIC_GROUPS == set()

    def test_oxygen_and_sulfur_groups_are_acidic(self):
        for group in ("carboxyl", "phenolic", "thiol"):
            assert PROTONATION_STATES[group].kind == "acidic"

    def test_nitrogen_groups_are_basic(self):
        for group in ("amino", "pyridinic"):
            assert PROTONATION_STATES[group].kind == "basic"

    def test_states_are_immutable(self):
        """Chemistry constants must not be mutable at runtime."""
        with pytest.raises(Exception):
            PROTONATION_STATES["carboxyl"].pKa = 1.0

    def test_neutral_and_ionized_types_differ(self):
        for group, state in PROTONATION_STATES.items():
            assert state.neutral_type != state.ionized_type, (
                f"{group} neutral and ionized map to the same type"
            )

    def test_ionized_types_are_unique_across_groups(self):
        ionized = [s.ionized_type for s in PROTONATION_STATES.values()]
        assert len(ionized) == len(set(ionized)), "two groups share an ionized type"

    def test_ph_bounds_span_the_aqueous_range(self):
        assert PH_MIN == 0.0
        assert PH_MAX == 14.0

    def test_neutral_types_exist_in_opls_table(self):
        """The neutral side is already implemented, so this must hold today."""
        for group, state in PROTONATION_STATES.items():
            assert state.neutral_type in OPLS_ATOM_TYPES, (
                f"{group} neutral type {state.neutral_type} missing from OPLS_ATOM_TYPES"
            )

    def test_ionized_types_exist_in_opls_table(self):
        """Every ionized type the table names must be a real OPLS type."""
        for group, state in PROTONATION_STATES.items():
            assert state.ionized_type in OPLS_ATOM_TYPES, (
                f"{group} ionized type {state.ionized_type} missing from OPLS_ATOM_TYPES"
            )

    def test_exchangeable_hydrogen_types_exist_in_opls_table(self):
        """The H removed (acidic) or added (basic) must be a real OPLS type."""
        for group, state in PROTONATION_STATES.items():
            assert state.h_type in OPLS_ATOM_TYPES, (
                f"{group} h_type {state.h_type} missing from OPLS_ATOM_TYPES"
            )


# ---------------------------------------------------------------------------
# U3 — the ProtonationAssigner engine
# ---------------------------------------------------------------------------

from biochar.protonation import ProtonationAssigner, fraction_ionized


class TestFractionIonized:
    """
    Henderson-Hasselbalch, with the sign convention that separates the two
    blocks. pKa always belongs to the protonated form; what "ionized" means
    flips between them.
    """

    def test_acidic_group_is_half_ionized_at_its_pka(self):
        assert fraction_ionized(4.2, 4.2, "acidic") == pytest.approx(0.5)

    def test_basic_group_is_half_ionized_at_its_pka(self):
        assert fraction_ionized(5.2, 5.2, "basic") == pytest.approx(0.5)

    def test_acidic_group_ionizes_as_ph_rises(self):
        low = fraction_ionized(4.2, 2.0, "acidic")
        high = fraction_ionized(4.2, 10.0, "acidic")
        assert low < 0.01 < 0.99 < high

    def test_basic_group_deionizes_as_ph_rises(self):
        """The inversion: a base is protonated at LOW pH, not high."""
        low = fraction_ionized(5.2, 2.0, "basic")
        high = fraction_ionized(5.2, 10.0, "basic")
        assert high < 0.01 < 0.99 < low

    def test_acidic_and_basic_are_complementary_at_equal_pka(self):
        for pH in (0.0, 3.0, 7.0, 11.0, 14.0):
            a = fraction_ionized(5.0, pH, "acidic")
            b = fraction_ionized(5.0, pH, "basic")
            assert a + b == pytest.approx(1.0)

    def test_carboxyl_is_essentially_ionized_at_neutral_ph(self):
        """The environmentally important case: biochar is anionic at pH 7."""
        assert fraction_ionized(4.2, 7.0, "acidic") > 0.99

    def test_phenolic_is_essentially_neutral_at_neutral_ph(self):
        assert fraction_ionized(9.5, 7.0, "acidic") < 0.05

    def test_pyridinic_is_mostly_neutral_at_neutral_ph(self):
        assert fraction_ionized(5.2, 7.0, "basic") < 0.05

    def test_fraction_is_always_a_probability(self):
        for pH in range(0, 15):
            for kind in ("acidic", "basic"):
                f = fraction_ionized(7.0, float(pH), kind)
                assert 0.0 <= f <= 1.0

    def test_unknown_kind_is_rejected(self):
        with pytest.raises(ValueError):
            fraction_ionized(4.2, 7.0, "amphoteric")


# --- fixtures ---------------------------------------------------------------

def _mol(smiles: str) -> Chem.Mol:
    return Chem.AddHs(Chem.MolFromSmiles(smiles))


def net_formal(mol) -> int:
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


# Benzene-1,3,5-tricarboxylic acid: three independent carboxyls on one ring.
TRIACID = "OC(=O)c1cc(C(=O)O)cc(C(=O)O)c1"
# A ring with many phenolic OH.
TRIPHENOL = "Oc1cc(O)cc(O)c1"


class TestProtonationAssignerAcidic:
    def test_carboxyls_deprotonate_at_high_ph(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=12.0)
        assert net_formal(out) == -3
        assert comp.ionized_counts["carboxyl"] == 3

    def test_carboxyls_stay_protonated_at_low_ph(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=1.0)
        assert net_formal(out) == 0
        assert comp.ionized_counts.get("carboxyl", 0) == 0

    def test_deprotonated_carboxyl_loses_its_acidic_hydrogen(self):
        out, _ = ProtonationAssigner(seed=1).assign(_mol("OC(=O)c1ccccc1"), pH=12.0)
        anionic = [a for a in out.GetAtoms()
                   if a.GetAtomicNum() == 8 and a.GetFormalCharge() == -1]
        assert len(anionic) == 1
        assert not [n for n in anionic[0].GetNeighbors() if n.GetAtomicNum() == 1]

    def test_deprotonated_carboxyl_keeps_its_carbonyl(self):
        """Only the -OH is removed; the C=O must survive untouched."""
        out, _ = ProtonationAssigner(seed=1).assign(_mol("OC(=O)c1ccccc1"), pH=12.0)
        assert sum(1 for a in out.GetAtoms() if a.GetAtomicNum() == 8) == 2

    def test_phenolics_deprotonate_at_high_ph(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol(TRIPHENOL), pH=13.5)
        assert net_formal(out) == -3
        assert comp.ionized_counts["phenolic"] == 3

    def test_phenolics_stay_protonated_at_neutral_ph(self):
        out, _ = ProtonationAssigner(seed=1).assign(_mol(TRIPHENOL), pH=7.0)
        assert net_formal(out) == 0

    def test_thiol_deprotonates_at_high_ph(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol("Sc1ccccc1"), pH=12.0)
        assert net_formal(out) == -1
        assert comp.ionized_counts["thiol"] == 1


class TestProtonationAssignerBasic:
    def test_aniline_protonates_at_low_ph(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol("Nc1ccccc1"), pH=1.0)
        assert net_formal(out) == +1
        assert comp.ionized_counts["amino"] == 1

    def test_protonated_aniline_gains_a_third_hydrogen(self):
        out, _ = ProtonationAssigner(seed=1).assign(_mol("Nc1ccccc1"), pH=1.0)
        n = next(a for a in out.GetAtoms() if a.GetAtomicNum() == 7)
        assert sum(1 for x in n.GetNeighbors() if x.GetAtomicNum() == 1) == 3

    def test_aniline_stays_neutral_at_neutral_ph(self):
        out, _ = ProtonationAssigner(seed=1).assign(_mol("Nc1ccccc1"), pH=7.0)
        assert net_formal(out) == 0

    def test_pyridinic_protonates_at_low_ph(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol("c1ccncc1"), pH=1.0)
        assert net_formal(out) == +1
        assert comp.ionized_counts["pyridinic"] == 1

    def test_protonated_pyridinic_gains_a_hydrogen(self):
        out, _ = ProtonationAssigner(seed=1).assign(_mol("c1ccncc1"), pH=1.0)
        n = next(a for a in out.GetAtoms() if a.GetAtomicNum() == 7)
        assert sum(1 for x in n.GetNeighbors() if x.GetAtomicNum() == 1) == 1


class TestSignInversion:
    """
    The single easiest way to get this feature silently wrong: treating the
    basic block like the acidic one inverts the charge on every nitrogen.
    """

    def test_raising_ph_increases_anions_but_decreases_cations(self):
        acid_low, _ = ProtonationAssigner(seed=2).assign(_mol(TRIACID), pH=1.0)
        acid_high, _ = ProtonationAssigner(seed=2).assign(_mol(TRIACID), pH=12.0)
        assert net_formal(acid_high) < net_formal(acid_low)

        base_low, _ = ProtonationAssigner(seed=2).assign(_mol("Nc1ccccc1"), pH=1.0)
        base_high, _ = ProtonationAssigner(seed=2).assign(_mol("Nc1ccccc1"), pH=12.0)
        assert net_formal(base_high) < net_formal(base_low)

    def test_titration_curve_is_monotonically_decreasing(self):
        charges = [
            net_formal(ProtonationAssigner(seed=7).assign(_mol(TRIACID), pH=p)[0])
            for p in (1.0, 3.0, 5.0, 7.0, 9.0, 12.0)
        ]
        assert charges == sorted(charges, reverse=True), charges
        assert charges[0] == 0 and charges[-1] == -3

    def test_amphoteric_molecule_is_cationic_low_and_anionic_high(self):
        smiles = "Nc1ccc(C(=O)O)cc1"  # 4-aminobenzoic acid
        low, _ = ProtonationAssigner(seed=3).assign(_mol(smiles), pH=1.0)
        high, _ = ProtonationAssigner(seed=3).assign(_mol(smiles), pH=12.0)
        assert net_formal(low) == +1
        assert net_formal(high) == -1


class TestDeterminism:
    def test_same_seed_and_ph_is_reproducible(self):
        a, _ = ProtonationAssigner(seed=42).assign(_mol(TRIACID), pH=4.2)
        b, _ = ProtonationAssigner(seed=42).assign(_mol(TRIACID), pH=4.2)
        assert net_formal(a) == net_formal(b)
        assert Chem.MolToSmiles(a) == Chem.MolToSmiles(b)

    def test_different_seeds_sample_differently_near_the_pka(self):
        """Proves it samples rather than thresholding."""
        results = {
            net_formal(ProtonationAssigner(seed=s).assign(_mol(TRIACID), pH=4.2)[0])
            for s in range(40)
        }
        assert len(results) > 1, "every seed gave the same answer at pH == pKa"

    def test_does_not_disturb_global_random_state(self):
        import random

        random.seed(123)
        expected = random.random()
        random.seed(123)
        ProtonationAssigner(seed=9).assign(_mol(TRIACID), pH=5.0)
        assert random.random() == expected


class TestStatistics:
    def test_about_half_ionized_at_the_pka(self):
        """Aggregate over seeds: the mean must track Henderson-Hasselbalch."""
        total, ionized = 0, 0
        for s in range(120):
            _, comp = ProtonationAssigner(seed=s).assign(_mol(TRIACID), pH=4.2)
            ionized += comp.ionized_counts.get("carboxyl", 0)
            total += 3
        assert 0.40 < ionized / total < 0.60, ionized / total


class TestUntouchedGroups:
    def test_ether_oxygen_is_never_titrated(self):
        out, comp = ProtonationAssigner(seed=1).assign(
            _mol("c1ccc2c(c1)Oc1ccccc1-2"), pH=14.0
        )
        assert net_formal(out) == 0
        assert comp.ionized_counts == {}

    def test_thioether_sulfur_is_never_titrated(self):
        out, _ = ProtonationAssigner(seed=1).assign(
            _mol("c1ccc2c(c1)Sc1ccccc1-2"), pH=14.0
        )
        assert net_formal(out) == 0

    def test_graphitic_nitrogen_is_never_titrated(self):
        from biochar.heteroatom_assignment import NitrogenSubstitutor

        mol = Chem.MolFromSmiles("c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61")
        sub = NitrogenSubstitutor(seed=1)
        doped = sub.substitute(mol, n_graphitic=1)
        assert sub.placed_graphitic == 1

        before = net_formal(doped)
        out, comp = ProtonationAssigner(seed=1).assign(doped, pH=1.0)
        assert net_formal(out) == before, "graphitic N was titrated"
        assert "pyridinic" not in comp.ionized_counts

    def test_molecule_with_no_titratable_sites_is_unchanged(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol("c1ccccc1"), pH=7.0)
        assert net_formal(out) == 0
        assert comp.net_charge == 0
        assert comp.ionized_counts == {}


class TestPhBounds:
    @pytest.mark.parametrize("pH", [-1.0, 15.0, 100.0])
    def test_ph_outside_the_aqueous_range_is_rejected(self, pH):
        with pytest.raises(ValueError, match="pH"):
            ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=pH)

    @pytest.mark.parametrize("pH", [PH_MIN, 7.0, PH_MAX])
    def test_bounds_themselves_are_accepted(self, pH):
        ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=pH)


class TestCompositionResult:
    def test_net_charge_matches_formal_charge(self):
        out, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=12.0)
        assert comp.net_charge == net_formal(out) == -3

    def test_ionized_counts_only_lists_groups_actually_ionized(self):
        _, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=1.0)
        assert all(v > 0 for v in comp.ionized_counts.values())

    def test_site_counts_report_what_was_available(self):
        _, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=7.0)
        assert comp.titratable_counts["carboxyl"] == 3


class TestOutputIsUsable:
    """The assigner's output must survive the rest of the pipeline."""

    def test_ionized_output_types_and_charges_correctly(self):
        from biochar.opls_typing import AtomTyper, ChargeAssigner

        out, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=12.0)
        types = AtomTyper().assign_atom_types(out)
        assert not [t for t in types.values() if t.startswith("X")]
        charges = ChargeAssigner().assign_charges(out, types)
        assert sum(charges.values()) == pytest.approx(comp.net_charge, abs=1e-6)

    def test_ionized_output_passes_valence_validation(self):
        from biochar.valence import ValenceValidator

        out, _ = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=12.0)
        ok, errors = ValenceValidator.validate_molecule(out)
        assert ok, errors

    def test_hydrogen_assignment_does_not_reprotonate_the_result(self):
        """End to end with U2: the anions must survive H saturation."""
        from biochar.heteroatom_assignment import HydrogenAssigner

        out, comp = ProtonationAssigner(seed=1).assign(_mol(TRIACID), pH=12.0)
        after, _ = HydrogenAssigner(seed=0).assign_hydrogens(out, target_H_C_ratio=1.0)
        assert net_formal(after) == comp.net_charge == -3
