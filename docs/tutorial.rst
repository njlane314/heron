Tutorial: end-to-end example
============================

This tutorial walks through an end-to-end workflow using the provided command
line applications and plotting macros. It assumes you have a FermiGrid
production file list ready.

Step 1: scan art provenance
---------------------------

.. code-block:: console

   nuxsec art nue_run1:data/run1_nue.list

This stage scans the art files and writes provenance metadata for later
normalisation.

Step 2: build a sample
----------------------

.. code-block:: console

   nuxsec sample nue_run1:data/run1_nue.list

The sample aggregation step writes a SampleIO ROOT file and updates the sample
list used by later stages.

Step 3: build event output
--------------------------

.. code-block:: console

   nuxsec event --list scratch/out/sample/samples.tsv \
     --output scratch/out/event/event_output.root

This stage applies analysis column derivations and writes an event tree that is
ready for selections and plotting.

Step 4: run a plotting macro
----------------------------

Macros live under ``plot/macro`` and can be run directly:

.. code-block:: console

   nuxsec macro plotFluxMinimal.C

A macro is a ROOT C++ script that expects the event output and workspace
configuration to be present. Adjust ``NUXSEC_PLOT_DIR`` or ``NUXSEC_PLOT_BASE``
if you want to separate outputs per production set.

Example macro snippet
---------------------

.. code-block:: c++

   void plotFluxMinimal() {
     TFile input("scratch/out/event/event_output.root", "READ");
     TTree *events = static_cast<TTree *>(input.Get("events"));
     events->Draw("nu_energy >> h_flux(50, 0.0, 5.0)");
   }

This snippet demonstrates how macros typically read the event tree and create a
histogram for later styling in the ``plot/`` module.
