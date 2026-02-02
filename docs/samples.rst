Samples and provenance
======================

Nuxsec models production inputs as samples built from art ROOT files. Each
sample tracks provenance, POT totals, and a normalisation factor to ensure
consistent scaling across analyses.

From FermiGrid production to SampleIO
-------------------------------------

FermiGrid production delivers art ROOT files. Nuxsec first scans them to record
run/subrun metadata, then aggregates them into a SampleIO file.

.. code-block:: console

   # Register art provenance for a production file list.
   nuxsec art nue_run1:data/run1_nue.list

   # Build the SampleIO file and update samples.tsv.
   nuxsec sample nue_run1:data/run1_nue.list

The sample step reads the beam database, sums POT, and writes a SampleIO ROOT
file that is referenced by ``samples.tsv``.

Sample list format
------------------

Sample lists live in ``scratch/out/sample/samples.tsv`` by default. They are
TSV files with a header row and per-sample lines such as:

.. code-block:: text

   # sample_name\tsample_origin\tbeam_mode\toutput_path
   nue_run1\tdata\tbeam\tscratch/out/sample/sample_root_nue_run1.root

Applications that build event outputs read this list to locate each sample
ROOT file and its metadata.

Normalisation inputs
--------------------

SampleIO stores:

* The list of input files and their provenance.
* POT totals from the art provenance scan.
* Beam database totals (for example ``tortgt`` and ``tor101``).
* Derived normalisation factors used when constructing event weights.

This ensures consistent scaling when multiple samples are combined in an
analysis.
