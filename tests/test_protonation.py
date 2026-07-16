"""
Tests for pH-dependent protonation — the pKa table (U1) and the
ProtonationAssigner engine (U3).
"""

import pytest

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

    @pytest.mark.xfail(
        strict=True,
        reason="ionized OPLS types are defined in U4; U1 only declares the contract",
    )
    def test_ionized_types_exist_in_opls_table(self):
        """
        The ionized side is defined in U4. Until that lands this fails, which is
        the point: U1 declares the contract, U4 satisfies it.

        strict=True — when U4 lands, this passing becomes an error, forcing the
        marker off rather than leaving a stale xfail behind.
        """
        for group, state in PROTONATION_STATES.items():
            assert state.ionized_type in OPLS_ATOM_TYPES, (
                f"{group} ionized type {state.ionized_type} missing from OPLS_ATOM_TYPES"
            )

    @pytest.mark.xfail(
        strict=True,
        reason="cationic H types (HNAP/HPYP) are defined in U4",
    )
    def test_exchangeable_hydrogen_types_exist_in_opls_table(self):
        """Acidic H types exist today; basic (added) H types arrive in U4."""
        for group, state in PROTONATION_STATES.items():
            assert state.h_type in OPLS_ATOM_TYPES, (
                f"{group} h_type {state.h_type} missing from OPLS_ATOM_TYPES"
            )
