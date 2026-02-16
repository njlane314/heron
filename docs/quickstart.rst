Quick start
===========

This quick start shows how to go from FermiGrid production outputs to
analysis-ready event trees and plots. It assumes you have art ROOT files from
production, a working heron build, and a configured workspace.

Prerequisites
-------------

* A local build of heron (see ``README.md`` for build steps).
* ROOT available in your environment.
* A file list of FermiGrid-produced art ROOT files for a sample.

Workflow overview
-----------------

heron uses a three-stage I/O pipeline for production outputs:

1. **Art provenance scan**: read the art ROOT files and capture run/subrun
   metadata for later POT normalisation.
2. **Sample aggregation**: combine provenance with beam database totals to
   produce a sample ROOT file and update the sample list.
3. **Event build**: load samples, define analysis columns, and write
   event-level output for analysis and plotting.

Quick start commands
--------------------

.. code-block:: console

   # 1) Register a production input with art provenance.
   heron art my_sample:data/filelist.txt

   # 2) Build a list of art provenance outputs from step (1).
   ls scratch/out/template/art/art_prov_my_sample*.root > scratch/out/template/lists/my_sample.txt

   # 3) Build a SampleIO output ROOT file and update samples.tsv.
   heron sample my_sample:scratch/out/template/lists/my_sample.txt

   # 4) Build the event-level output from the default samples list.
   heron event scratch/out/template/event/event_output.root

   # 5) Run a plotting macro from macros/plot/macro.
   heron macro plotFluxMinimal.C

   # 6) Run a standalone ROOT macro from macros/standalone/macro.
   heron macro plotOscPars.C

After step (4) you will have an event tree with analysis columns in the
specified output ROOT file. After step (5) plots are written under the plot
output directory (``scratch/plot/<set>`` by default). Step (6) shows a
standalone ROOT macro from ``macros/standalone/macro`` that does not require the
heron plotting libraries.

Workspace tips
--------------

Use ``--set`` or ``HERON_SET`` to switch workspace configurations. Plot outputs
are controlled by ``HERON_PLOT_BASE``, ``HERON_PLOT_DIR``, and
``HERON_PLOT_FORMAT``. Adjust these when comparing multiple productions or
analysis variants.
