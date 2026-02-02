Definitions and conventions
===========================

This page defines the core analysis terms used throughout Nuxsec. These
definitions align with the production workflow that starts from FermiGrid art
outputs and ends in analysis-ready ROOT files.

Core objects
------------

* **Art provenance**: metadata extracted from production art ROOT files to
  capture run, subrun, and POT information.
* **Sample**: a logical grouping of production inputs that share a label,
  origin, and beam mode. Samples map to a single SampleIO ROOT file.
* **Event output**: an analysis-ready ROOT file containing an event tree with
  derived columns and weights.
* **Selection**: a predicate used to filter events for an analysis step.
* **Systematic variation**: a controlled shift of a quantity or weight used to
  evaluate modelling uncertainty.

Normalisation and weighting
---------------------------

Event weights are assembled from generator weights and corrections:

.. math::

   w_{\mathrm{event}} = w_{\mathrm{gen}} \times \prod_{i=1}^{N_{\mathrm{corr}}} c_{i}.

The integrated exposure is captured with POT sums, with sample normalisation
tracked explicitly. A simplified normalisation is

.. math::

   w_{\mathrm{norm}} = \frac{\mathrm{POT}_{\mathrm{target}}}{\mathrm{POT}_{\mathrm{sample}}}.

These terms appear in both SampleIO metadata and event-level columns.

Sample list entries
-------------------

Sample lists are tab-separated tables with a header row and four columns:

* ``sample_name``
* ``sample_origin``
* ``beam_mode``
* ``output_path``

Each entry points to a SampleIO ROOT file produced by the sample aggregation
step. Downstream applications read this list to build event outputs.

Selections and efficiencies
----------------------------

Selection efficiency is defined as

.. math::

   \varepsilon = \frac{N_{\mathrm{pass}}}{N_{\mathrm{total}}}.

Efficiencies are typically recorded per dataset, per selection, and per
systematic variation to maintain a traceable cutflow.
