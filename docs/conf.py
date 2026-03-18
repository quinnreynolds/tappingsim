"""Sphinx configuration for tappingsim documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath('..'))

# -- Project information -------------------------------------------------------

project = 'tappingsim'
copyright = '2024, Quinn Reynolds'
author = 'Quinn Reynolds'

from importlib.metadata import version as get_version
release = get_version('tappingsim')
version = '.'.join(release.split('.')[:2])

# -- General configuration -----------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'nbsphinx',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# -- Napoleon settings (NumPy-style docstrings) --------------------------------

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True

# -- Autodoc settings ----------------------------------------------------------

autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}
autodoc_member_order = 'bysource'

# -- Intersphinx ---------------------------------------------------------------

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
}

# -- Options for HTML output ---------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    'navigation_depth': 3,
    'titles_only': False,
}

# -- nbsphinx ------------------------------------------------------------------

nbsphinx_execute = 'never'  # don't re-run notebooks during docs build
