"""The version is declared in more than one place; these keep them in step.

`biochar.__version__` and the version pyproject.toml ships in the package
metadata drifted apart before -- docs/conf.py sat at 0.1.4 through three
releases, and conda-recipe/meta.yaml still does. A hand-copied version only ever
drifts in one direction, and nothing was checking.

conda-recipe/meta.yaml is deliberately not checked here: its version is paired
with a sha256 of the published sdist, so it can only move when a release is
actually uploaded to PyPI. Bumping it early would produce a recipe that fails to
build.
"""

from importlib.metadata import version as installed_version

import biochar


def test_dunder_version_matches_package_metadata():
    """__version__ agrees with what pyproject.toml declared at install time."""
    assert biochar.__version__ == installed_version("biochar"), (
        f"biochar.__version__ is {biochar.__version__} but the installed "
        f"distribution reports {installed_version('biochar')}.\n"
        f"Either pyproject.toml and src/biochar/__init__.py have drifted apart, "
        f"or a stale biochar.egg-info/ is shadowing the real install -- the "
        f"repository root is on sys.path under pytest, and a pre-src-layout "
        f"egg-info left there wins over site-packages. Check with:\n"
        f"  python -c \"from importlib.metadata import distributions; "
        f"[print(d.version, d._path) for d in distributions() "
        f"if (d.metadata['Name'] or '').lower() == 'biochar']\""
    )


def test_docs_version_is_derived_not_restated():
    """docs/conf.py reads the version rather than hard-coding one.

    Checked structurally rather than by value: importing conf.py would execute
    Sphinx configuration, and the failure this guards against is someone
    reintroducing a literal.
    """
    from pathlib import Path

    conf = Path(__file__).resolve().parents[1] / "docs" / "conf.py"
    if not conf.exists():  # docs are not shipped in every checkout
        return

    text = conf.read_text()
    assert "importlib.metadata" in text, (
        "docs/conf.py no longer derives the version from package metadata"
    )
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith(("release =", "version =")) and '"' in stripped:
            literal = stripped.split("=", 1)[1].strip().strip('"')
            assert not literal[:1].isdigit(), (
                f"docs/conf.py hard-codes a version again: {stripped!r}"
            )
