"""Sphinx configuration for Biochar Simulator documentation."""

import os
import sys

# Make the package importable without installation
sys.path.insert(0, os.path.abspath(".."))

# ---------------------------------------------------------------------------
# Project information
# ---------------------------------------------------------------------------

project = "Biochar Simulator"
copyright = "2026, jolayfield"
author = "jolayfield"
release = "0.1.1"
version = "0.1"

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",           # Pull docstrings automatically
    "sphinx.ext.napoleon",          # Google / NumPy style docstrings
    "sphinx.ext.viewcode",          # "View source" links
    "sphinx.ext.intersphinx",       # Cross-link to NumPy, RDKit, etc.
    "sphinx_autodoc_typehints",     # Type hints in signature + body
    "myst_parser",                  # Markdown support for .md files
]

# ---------------------------------------------------------------------------
# Autodoc settings
# ---------------------------------------------------------------------------

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"   # Show types in the description, not signature
autodoc_typehints_format = "short"  # Use short names (np.ndarray not numpy.ndarray)

# ---------------------------------------------------------------------------
# Napoleon (NumPy / Google docstring parser)
# ---------------------------------------------------------------------------

napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_param = False          # Suppress redundant "Parameters" header
napoleon_use_rtype = False          # Suppress "Return type" line (in description already)
napoleon_preprocess_types = True

# ---------------------------------------------------------------------------
# Intersphinx
# ---------------------------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "rdkit": ("https://www.rdkit.org/docs", None),
}

# ---------------------------------------------------------------------------
# MyST (Markdown) parser
# ---------------------------------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "includehidden": True,
    "titles_only": False,
}
html_static_path = ["_static"]
html_show_sourcelink = True
html_show_sphinx = False

# ---------------------------------------------------------------------------
# General
# ---------------------------------------------------------------------------

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
pygments_style = "friendly"
nitpicky = False                    # Don't warn on unresolved cross-refs
