`gippa compare` warns when the runs differ in GIPPA minor version, tool version or image digest,
database version or operon source (on stderr, in `logs/events.jsonl` and in `compare.html`).
Runs whose output schema majors differ stop with exit 1 unless `--allow-cross-major` is given.
`run_manifest.json` and `datapackage.json` now record `output_schema_version` (`1.0`), the version
of the output tables and manifest, apart from the config's `schema_version`.
