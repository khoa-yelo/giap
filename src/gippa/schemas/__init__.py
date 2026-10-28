"""JSON Schemas shipped with GIPPA: run manifest and resolved config.

``OUTPUT_SCHEMA_VERSION`` (M11-14h; ADR-0060 §6) is the ``MAJOR.MINOR`` version of the run's
outputs (the tables and the manifest), written to ``run_manifest.json`` and carried by
``datapackage.json``. A new column, table or manifest key bumps the minor; a removed, renamed or
retyped one bumps the major. It is apart from ``schema_version``, the config format (Q-L13).
``1.0``: the first published output schema (spec §4, Q-L13, ADR-0060 §6). ``1.1``: the
manifest's ``annotation.translation_table`` (M11-16).
"""

from __future__ import annotations

__all__ = ["OUTPUT_SCHEMA_VERSION"]

#: keep equal to ``output_schema_version`` in ``datapackage.json``
OUTPUT_SCHEMA_VERSION = "1.1"
