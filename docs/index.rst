Nuxsec Documentation
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents

   usage

Overview
--------

Nuxsec is a C++ analysis toolkit organised into top-level modules that each
provide a focused slice of functionality. Source code is split into ``include/``
headers and ``src/`` implementations, with executables living in ``apps/``.
Documentation is built with Sphinx, and this index provides a map of the
repository so newcomers can quickly locate relevant components.

Repository layout
-----------------

Modules live at the repository root and follow a consistent layout:

* ``ana/``: analysis-level utilities and orchestration.
* ``io/``: input/output helpers and data access.
* ``plot/``: plotting utilities and presentation helpers.
* ``sample/``: sample definitions and metadata management.
* ``pot/``: probability or potential-related utilities.
* ``rdf/``: ROOT RDataFrame integrations.
* ``sel/``: selection logic and filtering.
* ``syst/``: systematic variation handling.
* ``stat/``: statistical tooling.

Each module typically contains:

* ``include/`` for public headers (``.hh``).
* ``src/`` for implementations (``.cc``).
* ``macro/`` for ROOT-style macros and scripts.

Supporting folders include:

* ``apps/`` for small executables and command-line entry points.
* ``docs/`` for Sphinx configuration and documentation sources.
* ``scripts/`` for helper scripts and workflows.
* ``lib/`` for build outputs (shared libraries).

Conventions
-----------

The codebase follows consistent style rules:

* Headers begin with ``/* -- C++ -- */`` and use concise include guards.
* Library code lives in the ``nuxsec`` namespace.
* Class names use descriptive PascalCase, while functions prefer
  ``lowercase_with_underscores``.
* Member variables often use ``m_`` or ``p_`` prefixes, depending on ownership.
* British-English spelling is preferred in code and documentation.

Building documentation
----------------------

From the repository root, the documentation can be built with:

.. code-block:: console

   sphinx-build docs _build

The generated HTML will be located in ``_build/``.

Documentation deployment
------------------------

Documentation is automatically deployed to GitHub Pages from the default
branch. The GitHub Actions workflow builds the Sphinx site, adds the
``.nojekyll`` marker, and publishes the ``_build/`` directory to the
``gh-pages`` branch on every push to the default branch. That means
documentation changes must be merged to the default branch before the
public site is updated.
