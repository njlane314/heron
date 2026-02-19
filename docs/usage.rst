heron Usage Guide
=================

Getting started
---------------

heron is organised as a set of C++ modules with headers in ``include/`` and
implementation files in ``src/``. The repository also ships small executables
in ``apps/`` and helper scripts in ``scripts/``. When exploring the codebase,
start by identifying the module that matches your task (for example
``ana/`` for analysis-level tooling or ``io/`` for data access).

The most effective way to learn the codebase is to trace a feature end-to-end:
begin with a command-line executable, locate the module that implements the
core logic, and then follow the associated headers to understand the data
structures involved. This approach reveals how datasets, selections, and
outputs are composed.

Prerequisites
-------------

Before extending heron, ensure you have:

* A working C++ toolchain that matches the project's build setup.
* ROOT installed and available in your environment.
* Familiarity with the module layout so you can place new logic in the
  appropriate directory.

If you are unsure where new code should live, start in ``ana/`` for analysis
orchestration and move outward to the supporting modules once responsibilities
are clear.

Typical workflow
----------------

The common workflow when extending heron is:

1. Add or update headers in the module ``include/`` folder.
2. Implement the corresponding source code in ``src/``.
3. Update or add an executable in ``apps/`` if you need a command-line entry
   point.
4. Add or update documentation in ``docs/`` so that users can discover the
   feature.

When refactoring, try to keep public interfaces stable. Prefer adding new
methods with clear, explicit parameters over changing existing signatures. This
minimises disruption for downstream analysis code.

Module orientation
------------------

Use these cues to decide where new functionality should live:

* **Data access**: ``io/`` handles file I/O, tree layouts, provenance capture,
  and normalisation helpers.
* **Analysis configuration and selections**: ``ana/`` provides ``RDataFrame``
  services, derived columns, and selection presets.
* **Plotting**: ``plot/`` keeps visualisation and plotting style logic separate
  from computation.
* **CLI drivers**: ``apps/`` hosts the command-line entry points.

Documentation checklist
-----------------------

Use this checklist when touching documentation:

* Update ``docs/index.rst`` with any new pages or navigation entries.
* Run ``sphinx-build docs _build`` to validate the site locally.
* Ensure your branch is merged to the default branch so the GitHub Pages
  deployment workflow publishes the new site.

To keep documentation consistent and polished:

* Use descriptive section headings and short introductory paragraphs.
* Prefer structured lists over dense paragraphs.
* Link back to relevant modules so readers can quickly navigate to the code.
* Keep spelling in British English, matching the codebase.

Where to look for data definitions
----------------------------------

Some module-level cues to help orient a new contributor:

* ``io/`` for input data layouts, provenance output, and SampleIO metadata.
* ``ana/`` for derived column definitions and selection logic.
* ``apps/config/event_columns.tsv`` for the event output schema reference.

If you need to locate column naming conventions or event layouts, search within
the ``io/`` module first, then trace any references into ``ana/``. The
``apps/config/event_columns.tsv`` file provides a reference list of tracked event
columns.

Contributing notes
------------------

heron uses British-English spelling in code and documentation, and applies
consistent naming for classes and functions. When adding new types, prefer
explicit PascalCase names (``Dataset``, ``Selection``, ``HistDef``), and use
``lowercase_with_underscores`` for methods.

When authoring documentation pages:

* Add a short overview section that states the purpose of the page.
* Provide a clear list of responsibilities or expected outcomes.
* Include references to relevant modules or executables to help readers pivot
  into the source code.
