Samples and provenance
======================

heron models production inputs as samples built from art ROOT files. Each
sample tracks provenance, POT totals, and a normalisation factor to ensure
consistent scaling across analyses.

From FermiGrid production to SampleIO
-------------------------------------

FermiGrid production delivers art ROOT files. heron first scans them to record
run/subrun metadata, then aggregates the resulting provenance outputs into a
SampleIO file.

.. code-block:: console

   # Register art provenance for a production file list.
   heron art nue_run1:data/run1_nue.list

   # Build a list of art provenance outputs from the previous step.
   ls scratch/out/template/art/art_prov_nue_run1*.root > scratch/out/template/lists/nue_run1.txt

   # Build the SampleIO file and update samples.tsv.
   heron sample nue_run1:scratch/out/template/lists/nue_run1.txt

The sample step reads the beam database, sums POT, and writes a SampleIO ROOT
file that is referenced by ``samples.tsv``.

Sample list format
------------------

Sample lists live in ``scratch/out/<set>/sample/samples.tsv`` by default. They
are TSV files with a header row and per-sample lines such as:

.. code-block:: text

   # sample_name\tsample_origin\tbeam_mode\toutput_path
   nue_run1\tdata\tbeam\tscratch/out/template/sample/sample_root_nue_run1.root

Applications that build event outputs read this list to locate each sample
ROOT file and its metadata. The ``<set>`` segment defaults to ``out`` and
is controlled by ``HERON_SET`` or ``heron --set``.

Normalisation inputs
--------------------

SampleIO stores:

* The list of input files and their provenance.
* POT totals from the art provenance scan.
* Beam database totals (for example ``tortgt`` and ``tor101``).
* Derived normalisation factors used when constructing event weights.

This ensures consistent scaling when multiple samples are combined in an
analysis.
