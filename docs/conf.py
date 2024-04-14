# Configuration file for the Sphinx documentation builder.
# Templated from both RTD tutorial and doped by SMTG-Bham.

# -- Project information

project = 'findTDE'
copyright = '2023, Alexander S. Hauck'
author = 'Alexander S. Hauck'

release = '0.1'
version = '1.0.3'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosectionlabel',
    'sphinx_book_theme',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_book_theme'

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = "findtde-logo.png"
html_title = "findTDE"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
# html_show_sphinx = True

html_theme_options = {
    "repository_url": "https://github.com/hauck-as/findTDE",
    "github_repo": "https://github.com/hauck-as/findTDE",
    "github_button": True,
    "github_user": "hauck-as", # Username
    "description": "Set of scripts to facilitate easy calculations of threshold displacement energies for materials in VASP/LAMMPS using ab initio/classical molecular dynamics.",
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_repository_button": True,
    "home_page_in_toc": True,
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "colab_url": "https://colab.research.google.com",
    },
}

# Adding “Edit Source” links on your Sphinx theme
html_context = {
    "display_github": True, # Integrate GitHub
    "github_user": "hauck-as", # Username
    "github_repo": "findTDE", # Repo name
    "github_version": "main", # Version
    "conf_py_path": "/docs/", # Path in the checkout to the docs root
}

# -- Options for EPUB output
epub_show_urls = 'footnote'