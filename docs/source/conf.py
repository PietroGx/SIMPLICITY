# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SIMPLICITY'
copyright = '2025 Pietro Gerletti, Jean-Baptiste Escudie' 
author = 'Pietro Gerletti, Jean-Baptiste Escudie'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',      # auto-generate docs from docstrings
    'sphinx.ext.napoleon',     # support for Google/NumPy-style docstrings
    'sphinx.ext.viewcode',     # show links to source code
    'sphinx.ext.todo',         # allow TODO directives (optional)
    'myst_parser'             # uncomment if using Markdown
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown'
}


templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
