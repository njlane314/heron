Module Technical Notes
======================

This page expands on the responsibilities of each module with technical
examples, equations, and snippets that reflect how analyses are typically
assembled. The intent is to give readers a compact, engineering-focused view of
how data, selections, and plots interact in the heron workflow.

io/ -- data access and provenance
---------------------------------

The ``io/`` module defines data access helpers and event layouts. It is where
ROOT trees, column names, and schema translations are centralised. The module
also handles provenance capture, run database queries, and sample
normalisation so later stages can rely on consistent metadata.

A representative event weight can be expressed as

.. math::

   w_{\mathrm{event}} = w_{\mathrm{gen}} \times \prod_{i=1}^{N_{\mathrm{corr}}} c_{i},

where the base generator weight is multiplied by corrections for pile-up,
trigger, and object efficiency scale factors. ``io/`` provides the consistent
column naming that keeps these terms discoverable across the codebase.

ana/ -- analysis configuration and RDataFrame services
------------------------------------------------------

The ``ana/`` module gathers analysis configuration, column derivations, and
selection logic. It wraps ROOT ``RDataFrame`` usage through ``RDataFrameService``
and pairs that with derived columns from ``ColumnDerivationService`` and
selection helpers in ``SelectionService``.

.. code-block:: c++

   ROOT::RDataFrame frame = RDataFrameService::load_sample(sample, "Events");
   auto enriched = ColumnDerivationService::instance().define(frame, processor);
   auto filtered = SelectionService::apply(
       enriched, Preset::Muon, selection_entry);

This module keeps analysis-level logic centralised so applications and macros
can use a consistent set of derived columns and selection labels.

plot/ -- plotting and presentation
----------------------------------

``plot/`` provides styling helpers and common plot construction patterns. It
keeps the presentation layer separate from analysis logic so that histogram
production and final rendering can evolve independently. The module typically
constructs standardised axis labels, colour palettes, and annotation helpers.

apps/ -- command-line drivers
-----------------------------

``apps/`` contains the executable entry points that orchestrate the analysis
pipeline. The unified ``heron`` CLI lives here and dispatches to the art,
sample, event, and macro drivers. Helper utilities such as logging and status
monitoring also live in this layer.

cols/ -- column schema references
---------------------------------

``cols/`` contains TSV references used by the event builder to define output
schemas. The ``cols/event_columns.tsv`` file enumerates event-level columns and
is consumed by the CLI when you override the compiled column list.

Putting it together
-------------------

A simplified flow that connects these modules is:

.. code-block:: c++

   ROOT::RDataFrame frame = RDataFrameService::load_sample(sample, "Events");
   auto enriched = ColumnDerivationService::instance().define(frame, processor);
   auto filtered = SelectionService::apply(
       enriched, Preset::Muon, selection_entry);

This pattern keeps the data definition and provenance in ``io/``, analysis
columns and selections in ``ana/``, and final presentation in ``plot/``.
