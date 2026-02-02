Nuxsec Documentation
====================

.. toctree::
   :maxdepth: 2
   :caption: Contents

   quickstart
   definitions
   samples
   tutorial
   applications
   usage
   modules
   api

Overview
--------

Nuxsec is a C++ analysis toolkit organised into top-level modules that each
provide a focused slice of functionality. Source code is split into ``include/``
headers and ``src/`` implementations, with executables living in ``apps/``.
Documentation is built with Sphinx, and this index provides a map of the
repository so newcomers can quickly locate relevant components.

The documentation is written for analysts, developers, and reviewers who need
both a quick orientation and a precise map of the repository. It aims to
highlight module responsibilities, establish shared conventions, and provide a
repeatable workflow for extending the toolkit.

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

Architecture at a glance
------------------------

Nuxsec takes a layered approach to analysis development:

* **Data access and layout** live in ``io/`` and ``rdf/`` to provide a consistent
  interface to ROOT trees and ``RDataFrame`` workflows.
* **Selections and corrections** are captured in ``sel/`` and ``syst/`` so they
  can be reused across analyses and systematics studies.
* **Analysis orchestration** is grouped under ``ana/``, typically responsible
  for stitching together datasets, selections, and outputs.
* **Statistical tools** in ``stat/`` provide shared machinery for inference and
  limit-setting workflows.
* **Plots and presentation** live in ``plot/`` to keep visualisation logic
  separate from event processing.

This separation keeps low-level data access stable while allowing higher-level
analysis code to iterate rapidly.

Conventions
-----------

The codebase follows consistent style rules:

* Headers begin with ``/* -- C++ -- */`` and use concise include guards.
* Library code lives in the ``nuxsec`` namespace.
* Class names use descriptive PascalCase, while functions prefer
  ``lowercase_with_underscores``.
* Member variables often use ``m_`` or ``p_`` prefixes, depending on ownership.
* British-English spelling is preferred in code and documentation.

Navigating the documentation
----------------------------

The Sphinx site is designed to be scanned quickly:

* **Overview pages** explain the structure of the repository and how modules
  relate to each other.
* **Usage guides** provide step-by-step workflows for development tasks.
* **Reference material** (when present) captures module-specific classes,
  functions, and data conventions.

Use the left navigation menu to jump between sections. When adding new pages,
ensure they appear in the ``toctree`` above so they are discoverable.

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
