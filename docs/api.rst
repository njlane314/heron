API highlights
==============

This section lists the main libraries and classes that make up the Nuxsec API.
It is not exhaustive; instead it points to the primary entry points that are
commonly used by the applications and analysis workflows.

Analysis services (ana/)
------------------------

* ``AnalysisConfigService``: provides analysis configuration and processor
  selection logic.
* ``ColumnDerivationService``: defines derived columns for event trees.
* ``RDataFrameService``: loads samples into ROOT ``RDataFrame`` instances.
* ``SelectionService``: centralises selection configuration and predicates.

Input/output utilities (io/)
----------------------------

* ``ArtFileProvenanceIO``: reads and writes art provenance summaries.
* ``SampleIO``: stores sample metadata and resolves input ROOT files.
* ``NormalisationService``: builds sample normalisation and POT totals.
* ``EventListIO``: writes event-level output trees and metadata.
* ``RunDatabaseService``: queries run database information.
* ``SubRunInventoryService``: scans input files to collect run/subrun totals.

Plotting utilities (plot/)
--------------------------

* ``PlotEnv``: resolves plot output directories and formatting.
* ``PlotDescriptors``: defines standard labels, axis titles, and annotations.
* ``PlotChannels``: groups samples into stacked plot channels.
* ``Plotter``: orchestrates histogram styling and rendering.
* ``StackedHist``: helper for stacked histogram construction.

Command-line helpers (apps/)
----------------------------

The CLI layer relies on reusable helpers in ``apps/include`` such as ``AppUtils``
for path resolution, ``AppLog`` for structured logging, and module-specific
parsers (``ArtCLI``, ``SampleCLI``, and ``EventCLI``).
