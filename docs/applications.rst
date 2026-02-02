Applications and workflows
==========================

Nuxsec ships command-line applications under ``apps/``. The unified
``nuxsec`` CLI dispatches to these tools, and each stage writes artefacts that
feed the next.

Unified CLI
-----------

The ``nuxsec`` command exposes the main workflow stages:

* ``nuxsec art``: scan art ROOT files and record provenance.
* ``nuxsec sample``: aggregate SampleIO files with normalisation.
* ``nuxsec event``: build event-level output with derived columns.
* ``nuxsec macro``: run plotting macros.
* ``nuxsec status``: print status for available binaries.
* ``nuxsec paths``: print resolved workspace paths.
* ``nuxsec env``: output environment exports for a workspace.

Standalone applications
-----------------------

The CLI maps to these executable entry points:

* ``nuxsecArtFileIOdriver``: provenance scanning and metadata capture.
* ``nuxsecSampleIOdriver``: sample aggregation and normalisation.
* ``nuxsecEventIOdriver``: event tree production from SampleIO inputs.

These are useful when integrating into batch workflows or production scripts.

Workflow summary
----------------

A typical workflow is:

1. Scan art files to register provenance.
2. Build a SampleIO file and update the sample list.
3. Produce event output with column derivations.
4. Run macros to produce plots.

Each application writes its outputs under ``scratch/out`` by default, and
plotting artefacts are written under ``scratch/plot`` unless overridden by
environment variables.
