EventColumnProvider
===================

``EventColumnProvider`` loads event output column lists from the compiled
defaults or a user-provided TSV schema.

Location
--------

* Header: ``apps/include/EventColumnProvider.hh``
* Namespace: ``global``

Default schema
--------------

The default event output schema is loaded from
``apps/config/event_columns.tsv`` when no override TSV path is provided.
