# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import re
from pathlib import Path

# -- Version ------------------------------------------------------------------
# Read the version from setup.py rather than importing the installed package --
# the docs build does not install EEfinder itself (it would pull in BLAST,
# DIAMOND and bedtools, which are not pip-installable).
_setup_py = Path(__file__).parent.parent / "setup.py"
_match = re.search(
    r"""^\s*version\s*=\s*["']([^"']+)["']""",
    _setup_py.read_text(),
    re.MULTILINE,
)
if _match is None:  # pragma: no cover - defensive
    raise RuntimeError(f"Could not find a version= assignment in {_setup_py}")
release = _match.group(1)
version = release

# -- Project information ------------------------------------------------------
project = "EEfinder"
copyright = "2024, WallauBioinfo"
author = "Yago J. M. Dias, Filipe Z. Dezordi & Gabriel L. Wallau"

# -- General configuration ----------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "myst_parser",
    "sphinx_copybutton",
    "sphinxcontrib.mermaid",
]

# MyST Parser configuration for Markdown support
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "html_admonition",
    "html_image",
    "strikethrough",
    "substitution",
    "tasklist",
]

myst_heading_anchors = 3

# Source file suffixes
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Language settings
language = "en"

# -- Options for HTML output --------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Theme options
html_theme_options = {
    "logo_only": False,
    "prev_next_buttons_location": "bottom",
    "style_external_links": True,
    "style_nav_header_background": "#2980B9",
    # Toc options
    "collapse_navigation": False,
    "sticky_navigation": True,
    "navigation_depth": 4,
    "includehidden": True,
    "titles_only": False,
}

# Custom CSS
html_css_files = [
    "custom.css",
]
