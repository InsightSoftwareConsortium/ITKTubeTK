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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'TubeTK'
copyright = '2021, Stephen R. Aylward'
author = 'Stephen R. Aylward'

# The full version, including alpha/beta/rc tags
release = '0.9.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
source_suffix = {".rst": "restructuredtext", ".txt": "restructuredtext", ".md": "markdown"}

extensions = [
    "recommonmark",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
    "sphinx_autodoc_typehints",
]

autoclass_content = "both"
add_module_names = True
autosectionlabel_prefix_document = True
napoleon_use_param = True
set_type_checking_flag = True

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
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_theme_options = {
    "collapse_navigation": True,
    "display_version": True,
    "sticky_navigation": True,  # Set to False to disable the sticky nav while scrolling.
    "logo_only": True,  # if we have a html_logo below, this shows /only/ the logo with no title text
    "style_nav_header_background": "#FBFBFB",
}
html_context = {
    "display_github": True,
    "github_user": "InsightSoftwareConsortium",
    "github_repo": "ITKTubeTK",
    "github_version": "master",
    "conf_py_path": "/docs/",
}
html_scaled_image_link = False
html_show_sourcelink = True
html_favicon = "images/favicon.ico"
html_logo = "images/TubeTK_Header.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = ["custom.css"]
