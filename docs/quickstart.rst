Quick start
===========

This quick start shows how to go from FermiGrid production outputs to
analysis-ready event trees and plots. It assumes you have art ROOT files from
production, a working Nuxsec build, and a configured workspace.

Prerequisites
-------------

* A local build of Nuxsec (see ``README.md`` for build steps).
* ROOT available in your environment.
* A file list of FermiGrid-produced art ROOT files for a sample.

Workflow overview
-----------------

Nuxsec uses a three-stage I/O pipeline for production outputs:

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
   nuxsec art my_sample:data/filelist.txt

   # 2) Build a SampleIO output ROOT file and update samples.tsv.
   nuxsec sample my_sample:data/filelist.txt

   # 3) Build the event-level output from the samples list.
   nuxsec event --list scratch/out/sample/samples.tsv \
     --output scratch/out/event/event_output.root

   # 4) Run a plotting macro from plot/macro.
   nuxsec macro plotFluxMinimal.C

After step (3) you will have an event tree with analysis columns in the
specified output ROOT file. After step (4) plots are written under the plot
output directory (``scratch/plot`` by default).

Workspace tips
--------------

Use ``--set`` or ``NUXSEC_SET`` to switch workspace configurations. Plot outputs
are controlled by ``NUXSEC_PLOT_BASE``, ``NUXSEC_PLOT_DIR``, and
``NUXSEC_PLOT_FORMAT``. Adjust these when comparing multiple productions or
analysis variants.
