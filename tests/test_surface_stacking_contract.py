"""
Tests for the rqm/surface-stacking.md scenarios -- what geometry a surface has,
and whether the exported files say the same thing about it.

Distinct from tests/test_surface_builder.py, which checks that the builder runs
and produces files of the right shape. These check the claims those files make:
the sheet separation a caller asked for, the geometry the .gro title reports, and
the agreement between the topology and the include files beside it.

Surfaces are expensive to build, so the fixtures are module-scoped and kept
small (16 carbons, 2-4 sheets).

Every scenario here passes. The three that carried xfail(strict=True) were
retired by XPASS when their fixes landed.
"""

import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest

from biochar.constants import CARBON_VDW_DIAMETER
from biochar.workflows.surface_builder import (
    SurfaceBuilder,
    SurfaceConfig,
    generate_surface,
)


def _slit(**extra):
    kwargs = dict(target_num_carbons=16, num_sheets=3, pore_diameter=10.0,
                  seed=1, strict=False)
    kwargs.update(extra)
    return SurfaceBuilder(SurfaceConfig(**kwargs))


def _unpackable(**extra):
    """An amorphous request that cannot converge: too much separation, few tries."""
    kwargs = dict(target_num_carbons=16, num_sheets=4, pore_type="amorphous",
                  min_separation=60.0, max_attempts=3, seed=1, strict=False)
    kwargs.update(extra)
    return SurfaceConfig(**kwargs)


@pytest.fixture(scope="module")
def slit_surface():
    builder = _slit()
    sheets, box = builder.build()
    return builder, sheets, box


class TestSheetSeparation:
    # rq-dbd8829a
    def test_adjacent_centroids_are_a_pore_diameter_plus_a_sheet_apart(self, slit_surface):
        builder, sheets, _ = slit_surface
        expected = builder.config.pore_diameter + CARBON_VDW_DIAMETER

        centroids = sorted(float(s.coords[:, 2].mean()) for s in sheets)
        gaps = np.diff(centroids)
        assert len(gaps) == len(sheets) - 1
        assert np.allclose(gaps, expected, atol=0.05), (
            f"centroid gaps {gaps} do not match pore_diameter + "
            f"{CARBON_VDW_DIAMETER} = {expected}; pore_diameter is the gap between "
            f"vdW surfaces, not between centres"
        )


class TestExportedFilesDescribeTheGeometry:
    @staticmethod
    def _build_with_fallback():
        builder = SurfaceBuilder(_unpackable(amorphous_fallback="slit"))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            builder.build()
        return builder, caught

    # rq-4e6abc16
    def test_the_fallback_is_reported(self):
        _, caught = self._build_with_fallback()
        messages = [str(w.message) for w in caught]
        assert any("slit" in m.lower() for m in messages), (
            f"packing failed and fell back with no warning naming the substituted "
            f"geometry; warnings were {messages}"
        )

    # rq-2c190d16
    def test_without_a_fallback_the_build_raises(self):
        with pytest.raises(RuntimeError):
            SurfaceBuilder(_unpackable(amorphous_fallback=None)).build()

    # rq-0a48ccd8
    def test_a_fallen_back_surface_is_not_labelled_amorphous(self):
        builder, _ = self._build_with_fallback()
        out = Path(tempfile.mkdtemp())
        builder.export_gromacs(str(out), basename="s")

        title = (out / "s.gro").read_text().splitlines()[0]
        # The claim the title leads with is the geometry it asserts; the word
        # "amorphous" may still appear afterwards as provenance.
        assert title.strip().lower().startswith("slit"), (
            f"a slit stack was written asserting a different geometry: {title!r}"
        )
        assert "fallback" in title.lower() or "failed" in title.lower(), (
            f"the title does not record that the requested geometry was "
            f"substituted: {title!r}"
        )

    # rq-6e562ff0
    def test_the_realised_geometry_is_readable_from_the_builder(self):
        builder, _ = self._build_with_fallback()

        assert builder.config.pore_type == "amorphous", "the request is unchanged"
        assert builder.realised_pore_type == "slit"
        assert builder.packing_fell_back is True


class TestCopiedSheetsShareNothing:
    # rq-1f914136
    def test_copied_sheets_do_not_share_coordinates(self, slit_surface):
        _, sheets, _ = slit_surface
        before = sheets[1].coords[0, 2]
        sheets[0].coords[0, 2] += 999.0
        try:
            assert sheets[1].coords[0, 2] == before, "sheets share a coordinate array"
        finally:
            sheets[0].coords[0, 2] -= 999.0

    # rq-9cdd2239
    def test_copied_sheets_do_not_share_a_composition_record(self, slit_surface):
        _, sheets, _ = slit_surface
        shared = [
            (i, j)
            for i in range(len(sheets))
            for j in range(i + 1, len(sheets))
            if sheets[i].composition is sheets[j].composition
        ]
        assert not shared, (
            f"sheet pairs {shared} share one composition object; mutating one "
            f"changes the others"
        )


class TestPerSheetProtonationSeeds:
    # rq-6ba78ac5
    def test_a_protonated_surface_is_reproducible_from_its_seed(self):
        def build():
            sheets, _ = _slit(num_sheets=2, pH=7.0, seed=42).build()
            return [s.coords.copy() for s in sheets]

        first, second = build(), build()
        assert len(first) == len(second) == 2
        for a, b in zip(first, second):
            np.testing.assert_allclose(a, b)

    # rq-682a21cd
    def test_the_surface_consumes_a_seed_range(self):
        """Seeds S..S+N-1, so ensembles at consecutive seeds overlap.

        Asserted on molecular identity, not coordinates: sheets are z-offset by
        their stack position after generation, so the same molecule in two
        surfaces never has equal coordinates. The precise consequence is
        positional -- sheet i+1 of the surface seeded S is the sheet seeded
        S+1, which is sheet i of the surface seeded S+1.
        """
        from rdkit import Chem

        def identities(seed):
            sheets, _ = _slit(num_sheets=3, pH=7.0, seed=seed).build()
            return [Chem.MolToSmiles(s.mol) for s in sheets]

        a, b = identities(100), identities(101)
        assert a[1:] == b[:-1], (
            f"sheet i+1 of seed 100 does not match sheet i of seed 101, so the "
            f"per-sheet seed is no longer S+index and the documented range is "
            f"stale.\n  seed 100: {[x[:24] for x in a]}\n  seed 101: {[x[:24] for x in b]}"
        )


class TestTopologyMatchesIncludes:
    @staticmethod
    def _export(builder):
        out = Path(tempfile.mkdtemp())
        gro, top, itps = builder.export_gromacs(str(out), basename="s")
        return Path(top).read_text(), [Path(p) for p in itps]

    @staticmethod
    def _molecules_section(top_text):
        rows, inside = [], False
        for raw in top_text.splitlines():
            if raw.strip().startswith("["):
                inside = raw.strip().startswith("[ molecules")
                continue
            line = raw.split(";", 1)[0].strip()
            if inside and line:
                parts = line.split()
                if len(parts) >= 2:
                    rows.append((parts[0], int(parts[1])))
        return rows

    # rq-17ffc86b
    def test_identical_sheets_are_written_once_and_counted(self, slit_surface):
        builder, sheets, _ = slit_surface
        top_text, itps = self._export(builder)

        assert len(itps) == 1, f"identical sheets produced {len(itps)} include files"
        rows = self._molecules_section(top_text)
        assert rows, "no [ molecules ] rows parsed"
        assert sum(count for _, count in rows) == len(sheets), (
            f"[ molecules ] totals {sum(c for _, c in rows)} for {len(sheets)} sheets"
        )

    # rq-884c2220
    def test_distinct_sheets_each_get_an_include_defining_their_moleculetype(self):
        # sheet_overrides defeats the identical-sheet optimisation.
        builder = _slit(
            num_sheets=2,
            sheet_overrides=[{"O_C_ratio": 0.05}, {"O_C_ratio": 0.20}],
        )
        builder.build()
        top_text, itps = self._export(builder)

        assert len(itps) == 2, f"distinct sheets produced {len(itps)} include files"

        defined = set()
        for path in itps:
            text = path.read_text()
            inside = False
            for raw in text.splitlines():
                if raw.startswith("["):
                    inside = raw.startswith("[ moleculetype")
                    continue
                line = raw.split(";", 1)[0].strip()
                if inside and line:
                    defined.add(line.split()[0])
                    inside = False

        named = {name for name, _ in self._molecules_section(top_text)}
        assert named <= defined, (
            f"the topology names moleculetypes no include defines: {named - defined}"
        )


class TestConvenienceWrapper:
    # rq-bfb7e848
    def test_the_fallback_is_reachable_from_the_convenience_function(self):
        out = Path(tempfile.mkdtemp())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            generate_surface(
                output_directory=str(out),
                target_num_carbons=16,
                num_sheets=4,
                pore_type="amorphous",
                min_separation=60.0,
                max_attempts=3,
                amorphous_fallback="slit",
                seed=1,
                strict=False,
            )
        assert list(out.glob("*.gro")), "no structure written"
