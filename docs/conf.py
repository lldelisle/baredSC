# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#

import os
import sys

# from unittest.mock import Mock

import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

sys.path.insert(0, os.path.abspath('../'))

from baredSC._version import __version__

# import mock

# MOCK_MODULES = ['numpy',
#                  'matplotlib',
#                  'pandas',
#                  'itertools',
#                  'scipy',
# MOCK_MODULES = ['scaledAdaptiveMetropolis',
#                 'covarianceImportanceSampling',
#                 'logpriors',
#                 'acf',
#                 'lib.common',
#                 'lib.oned',
#                 'lib.twod']
# MOCK_MODULES = ['samsam']
# for mod_name in MOCK_MODULES:
#     sys.modules[mod_name] = Mock()

# autodoc_mock_imports = MOCK_MODULES

# -- Project information -----------------------------------------------------

project = 'baredSC'
copyright = '2021, Jean-Baptiste Delisle, Lucille Lopez-Delisle'
author = 'Jean-Baptiste Delisle, Lucille Lopez-Delisle'

# The full version, including alpha/beta/rc tags
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinxarg.ext',
              'sphinx_autorun',
              'IPython.sphinxext.ipython_console_highlighting',
              'IPython.sphinxext.ipython_directive'
              ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
