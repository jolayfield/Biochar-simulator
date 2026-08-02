"""
Tests for the rqm/carbon-skeleton.md scenarios -- what a caller asks the
assembler for, and how faithfully the request is answered.

Distinct from the size and quality checks elsewhere: these are about fidelity.
A skeleton is assembled from whole rings, so several requests can only be
approximated, and the point is which approximations are honest. A parameter
accepted and ignored, a count met by a different number, or a fallback returning
a fixed 16-carbon structure are all the same failure wearing different clothes.

Sizes are read from PAH_LIBRARY rather than hard-coded, so extending the library
moves the tests with it.
"""

import inspect

import pytest

from biochar.constants import PAH_LIBRARY
from biochar.pipeline.carbon_skeleton import PAHAssembler


def _library_sizes():
    return sorted({int(v["num_atoms"]) for v in PAH_LIBRARY.values()})


def _generate(target, **kw):
    return PAHAssembler(seed=1).generate(target_num_carbons=target, **kw)


def _carbons(skeleton):
    return sum(1 for a in skeleton.mol.GetAtoms() if a.GetAtomicNum() == 6)


def _ring_sizes(skeleton):
    return [len(r) for r in skeleton.mol.GetRingInfo().AtomRings()]


class TestLibraryTargets:
    # rq-4af71330
    @pytest.mark.parametrize("target", _library_sizes())
    def test_a_library_size_is_met_exactly(self, target):
        assert _carbons(_generate(target)) == target


class TestGrowthReachesFromAbove:
    # Targets chosen to sit between reachable ring counts.
    AWKWARD = [41, 43, 45, 47, 49, 53, 57, 61, 71, 83, 97]

    # rq-59e23c2e
    @pytest.mark.parametrize("target", AWKWARD)
    def test_a_grown_skeleton_is_never_smaller_than_requested(self, target):
        actual = _carbons(_generate(target))
        assert actual >= target, (
            f"requested {target} carbons and got {actual}; a skeleton short of "
            f"the request has fewer edge sites than the caller planned for"
        )

    # rq-5435c9d8
    def test_the_overshoot_stays_within_a_ring_or_two(self):
        deltas = {t: _carbons(_generate(t)) - t for t in self.AWKWARD}
        worst = max(deltas.values())
        assert worst <= 8, (
            f"overshoot of {worst} carbons is more than a ring or two: {deltas}"
        )


class TestDefectFractions:
    # rq-c6b675bb
    def test_zero_defect_fraction_gives_a_pure_hexagonal_skeleton(self):
        sizes = _ring_sizes(_generate(80, defect_fraction=0.0, heptagon_fraction=0.0))
        assert sizes, "skeleton has no rings"
        assert set(sizes) == {6}, f"expected only hexagons, got sizes {sorted(set(sizes))}"

    # rq-0c91f5d5
    def test_a_positive_defect_fraction_produces_pentagons_below_the_probability(self):
        requested = 0.3
        sizes = _ring_sizes(_generate(80, defect_fraction=requested))
        pentagons = sum(1 for s in sizes if s == 5)

        assert pentagons > 0, "a positive defect fraction produced no pentagons"
        realised = pentagons / len(sizes)
        # Only rings added during growth are subject to the draw; the seed is
        # all hexagons, so the realised share sits below the probability.
        assert realised <= requested, (
            f"realised pentagon share {realised:.3f} exceeds the requested "
            f"per-ring probability {requested}"
        )


class TestAromaticityIsNotARequest:
    # rq-b6eec273
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "generate() still accepts target_aromaticity, documented as 'Unused "
            "(kept for backward compatibility)'. It reads at the call site as a "
            "request being made. Retire this marker with the fix."
        ),
    )
    def test_the_assembler_takes_no_aromaticity_argument(self):
        params = inspect.signature(PAHAssembler.generate).parameters
        aromatic = [p for p in params if "aromatic" in p.lower()]
        assert not aromatic, (
            f"generate() accepts {aromatic}, which it cannot honour; aromaticity "
            f"is decided by the ring topology"
        )


class TestNoUnrelatedFallback:
    # rq-f953b9b0
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "When _build_from_seed returns None the assembler substitutes pyrene "
            "(16 carbons) regardless of the request. Retire this marker with the fix."
        ),
    )
    def test_a_skeleton_that_cannot_be_grown_raises(self, monkeypatch):
        monkeypatch.setattr(
            PAHAssembler, "_build_from_seed", lambda *a, **k: None
        )
        with pytest.raises(Exception) as exc:
            _generate(200)
        assert not isinstance(exc.value, AssertionError)


class TestLibraryLoadIsReported:
    # rq-4dc5ce35
    @pytest.mark.xfail(
        strict=True,
        reason=(
            "A library entry whose SMILES will not parse or sanitise is dropped "
            "into a bare `except Exception: pass`, so the working library shrinks "
            "with nothing said. Retire this marker with the fix."
        ),
    )
    def test_an_unparseable_library_entry_is_reported(self, monkeypatch, caplog):
        import logging

        broken = dict(PAH_LIBRARY)
        broken["not_a_molecule"] = {
            "smiles": "this is not smiles",
            "num_atoms": 99,
            "num_aromatic": 99,
            "molecular_formula": "C99",
            "references": "fixture",
        }
        monkeypatch.setattr(
            "biochar.pipeline.carbon_skeleton.PAH_LIBRARY", broken, raising=False
        )
        monkeypatch.setattr(PAHAssembler, "_SIZE_INDEX", None)

        with caplog.at_level(logging.WARNING, logger="biochar"):
            assembler = PAHAssembler(seed=1)

        assert "not_a_molecule" not in assembler.pahs, "fixture entry unexpectedly loaded"
        assert caplog.records, (
            "a library entry was dropped from the working set with nothing "
            "reported; PAH_LIBRARY claims every entry is validated"
        )
