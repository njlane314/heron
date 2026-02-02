Module Technical Notes
======================

This page expands on the responsibilities of each module with technical
examples, equations, and snippets that reflect how analyses are typically
assembled. The intent is to give readers a compact, engineering-focused view of
how data, selections, systematics, and statistical summaries interact in the
Nuxsec workflow.

io/ -- data access and layout
-----------------------------

The ``io/`` module defines data access helpers and event layouts. It is where
ROOT trees, column names, and schema translations are centralised. A typical
pattern is to map raw tree branches into a structured representation and
compute derived columns that downstream modules can share.

A representative event weight can be expressed as

.. math::

   w_{\mathrm{event}} = w_{\mathrm{gen}} \times \prod_{i=1}^{N_{\mathrm{corr}}} c_{i},

where the base generator weight is multiplied by corrections for pile-up,
trigger, and object efficiency scale factors. ``io/`` provides the consistent
column naming that keeps these terms discoverable across the codebase.

rdf/ -- RDataFrame integrations
-------------------------------

The ``rdf/`` module wraps ROOT ``RDataFrame`` usage for column definitions,
filters, and higher-level compositing. A typical workflow uses declarative
column definitions and explicit filters:

.. code-block:: c++

   ROOT::RDataFrame frame("Events", input_files);
   auto weighted = frame.Define("event_weight", "genWeight * puWeight")
                         .Filter("nJet >= 2", "dijet_preselection");

This module focuses on consistent naming, memory-safe definitions, and
reusable filter snippets that are shared with analysis code in ``ana/``.

sel/ -- selections and object definitions
-----------------------------------------

Selections are defined in ``sel/`` as explicit filters that can be combined.
The key metric is selection efficiency,

.. math::

   \varepsilon = \frac{N_{\mathrm{pass}}}{N_{\mathrm{total}}},

which is tracked per dataset and per selection stage. Typical filters are
structured as named predicates so they can be combined and reported.

.. code-block:: c++

   bool pass_dijet_selection(const Event& event) {
     if (event.n_jet < 2) {
       return false;
     }
     return event.jet_pt[0] > 30.0 && event.jet_pt[1] > 30.0;
   }

syst/ -- systematic variations
------------------------------

Systematic variations are captured in ``syst/`` by defining explicit
perturbations of nominal quantities. A common modelling pattern is a relative
shift,

.. math::

   x^{\pm} = x_{\mathrm{nominal}} \times (1 \pm \delta),

where :math:`\delta` is the fractional uncertainty. Systematics are applied
consistently to either event weights or kinematic quantities to propagate
uncertainty envelopes.

.. code-block:: c++

   double apply_scale(const double nominal, const double variation) {
     return nominal * (1.0 + variation);
   }

ana/ -- orchestration and workflows
-----------------------------------

The ``ana/`` module coordinates datasets, selections, and output definitions.
It typically receives pre-defined selection predicates, applies them to
``RDataFrame`` chains, and writes histograms or summary tables. An analysis
workflow often composes multiple selections into a cutflow.

.. math::

   N_{\mathrm{yield}} = \sum_{k \in \mathrm{events}} w_{k} \times I_{k},

with indicator :math:`I_{k}` set by the selection. This is the aggregation
stage where yields and efficiencies are recorded.

stat/ -- statistical tooling
----------------------------

The ``stat/`` module implements statistical summaries and inference helpers.
For discovery-style tests, a common approximation for the profile likelihood
significance is

.. math::

   Z = \sqrt{2\left[(s + b)\ln\left(1 + \frac{s}{b}\right) - s\right]},

where :math:`s` and :math:`b` are the signal and background yields. These tools
are used to standardise reporting and to keep analysis code concise.

plot/ -- plotting and presentation
----------------------------------

``plot/`` provides styling helpers and common plot construction patterns. It
keeps the presentation layer separate from analysis logic so that histogram
production and final rendering can evolve independently. The module typically
constructs standardised axis labels, colour palettes, and annotation helpers.

sample/ -- datasets and metadata
--------------------------------

The ``sample/`` module describes datasets, metadata, and data-taking periods.
It supplies consistent identifiers (for example dataset names, labels, or
luminosity blocks) so that analysis code can iterate over samples without
hard-coding file lists.

pot/ -- probabilities or potentials
-----------------------------------

The ``pot/`` module hosts probability- or potential-related utilities that are
used for likelihood construction or modelling of efficiencies. The exact
contents vary by analysis, but the intent is to keep mathematically oriented
utilities separate from selection logic and plotting.

Putting it together
-------------------

A simplified flow that connects these modules is:

.. code-block:: c++

   ROOT::RDataFrame frame("Events", input_files);
   auto enriched = frame.Define("event_weight", "genWeight * puWeight")
                        .Filter("nJet >= 2", "dijet_preselection");

   auto hist = enriched.Histo1D({"h_mjj", "m_{jj}", 50, 0.0, 2000.0},
                               "mjj", "event_weight");

This pattern keeps the data definition in ``io/``, selection logic in ``sel/``,
workflow wiring in ``ana/``, and plots in ``plot/``.
