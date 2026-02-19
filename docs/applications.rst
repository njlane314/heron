Applications and workflows
==========================

heron ships command-line applications under ``apps/``. The unified
``heron`` CLI dispatches to these tools, and each stage writes artefacts that
feed the next.

Unified CLI
-----------

The ``heron`` command exposes the main workflow stages:

* ``heron art``: scan art ROOT files and record provenance.
* ``heron sample``: aggregate SampleIO files with normalisation.
* ``heron event``: build event-level output with derived columns.
* ``heron macro``: run plotting macros.
* ``heron status``: print status for available binaries.
* ``heron paths``: print resolved workspace paths.
* ``heron env``: output environment exports for a workspace.

Standalone applications
-----------------------

The CLI maps to these executable entry points:

* ``heronArtFileIOdriver``: provenance scanning and metadata capture.
* ``heronSampleIOdriver``: sample aggregation and normalisation.
* ``heronEventIOdriver``: event tree production from SampleIO inputs.

These are useful when integrating into batch workflows or production scripts.

Workflow summary
----------------

A typical workflow is:

1. Scan art files to register provenance.
2. Build a SampleIO file and update the sample list.
3. Produce event output with column derivations.
4. Run macros to produce plots.

Each application writes its outputs under ``scratch/out/<set>`` by default, and
plotting artefacts are written under ``scratch/plot/<set>`` unless overridden by
environment variables. The ``<set>`` segment defaults to ``out`` and is
controlled by ``HERON_SET`` or ``heron --set``.
