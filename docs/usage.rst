Nuxsec Usage Guide
=================

Getting started
---------------

Nuxsec is organised as a set of C++ modules with headers in ``include/`` and
implementation files in ``src/``. The repository also ships small executables
in ``apps/`` and helper scripts in ``scripts/``. When exploring the codebase,
start by identifying the module that matches your task (for example
``ana/`` for analysis-level tooling or ``sel/`` for selection logic).

Typical workflow
----------------

The common workflow when extending Nuxsec is:

1. Add or update headers in the module ``include/`` folder.
2. Implement the corresponding source code in ``src/``.
3. Update or add an executable in ``apps/`` if you need a command-line entry
   point.
4. Add or update documentation in ``docs/`` so that users can discover the
   feature.

Documentation checklist
-----------------------

Use this checklist when touching documentation:

* Update ``docs/index.rst`` with any new pages or navigation entries.
* Run ``sphinx-build docs _build`` to validate the site locally.
* Ensure your branch is merged to the default branch so the GitHub Pages
  deployment workflow publishes the new site.

Where to look for data definitions
----------------------------------

Some module-level cues to help orient a new contributor:

* ``io/`` and ``rdf/`` modules for input data layouts and ROOT
  ``RDataFrame`` usage.
* ``sel/`` for selection logic and filter definitions.
* ``syst/`` for systematic variation handling.
* ``stat/`` for statistical model helpers.

Contributing notes
------------------

Nuxsec uses British-English spelling in code and documentation, and applies
consistent naming for classes and functions. When adding new types, prefer
explicit PascalCase names (``Dataset``, ``Selection``, ``HistDef``), and use
``lowercase_with_underscores`` for methods.
