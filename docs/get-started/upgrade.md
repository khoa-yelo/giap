# Upgrade

Each GIAP release pins its own tool images (by digest, in the bundled
`workflow/conf/containers.config` and `workflow/conf/images.json`) and its own databases (version,
URL and checksum, in the bundled `src/giap/data/db_sources.yaml`). An upgrade is therefore a new
package plus whatever images and databases its pins name.

## 1. Install the new version

```bash
pip install -U giap==1.y.z
giap --version
```

Releases follow semantic versioning ([ADR-0060](../adr/ADR-0060-release-versioning-and-image-keys.md)).
Release candidates (`1.y.0rc1`, …) are uploaded to TestPyPI only, never to PyPI; `pip install -U`
from PyPI gets finals only.

## 2. Pull the images that changed

Run step 3 of [Install](install.md#the-five-steps) again with the same options:

```bash
giap setup --db-dir /path/to/db --images --profile apptainer   # or --profile docker
```

- An image whose pin is unchanged is kept: under Apptainer its `.img` file is already in the
  cache, and its `image:<tool>` record in `db_manifest.json` stays as it was.
- An image whose digest changed has a new file name, so it is pulled; the record is updated.
- `--check` verifies databases only. `giap setup --images --check` pulls nothing and does not
  compare images; to see which images a release changed, compare the bundled `images.json` and
  `containers.config` of the two versions, or read the release notes.
- **UniOP.** A release whose UniOP build inputs changed reports the old local build as stale in
  preflight (`local build <key_old> does not match this release (<key_new>)`), with the
  setup command to rerun. `setup --images` rebuilds it. A release with the same key reuses the
  build, and a rerun prints `image:uniop: ok (the <runtime> build is current)`.

## 3. Check the database pins

```bash
giap setup --db-dir /path/to/db --check
```

Each item prints `<name>: ok v<version>` or its problem. An installed Bakta database of another
major version than the pin, or a geNOMAD database of another version, is `mismatch` with
`installed <v_old>, this release pins <v_new>`, and the exit code is 2. Then run `giap setup`
without `--check` (same flags) to install the pinned version. An unchanged pin means the DB
directory is reused as it is.

## 4. Compare runs across versions

`giap compare` reads each run's `run_manifest.json` and prints one warning per difference to
stderr. The same lines go to `CMP/logs/events.jsonl` (step `compare`, level `WARNING`) and to
`CMP/report/compare.html`:

| Runs differ in | Line |
|---|---|
| package major or minor (a patch alone gives nothing) | `compare.cross_version: runs from 1.0.2 and 1.1.0` |
| `output_schema_version` major | `compare.output_schema: output schema 1.x and 2.x differ` |
| a tool's `version` or `container_digest` (a local UniOP build: `source_sha256`, `key`, `runtime`) | `compare.tool_versions: <tool> <version> sha256:<12 hex> vs …`, one line per tool |
| a database `version` | `compare.db_versions: <name> <v1> vs <v2>`, one line per database |
| `operon_source` | `compare.operon_source: run <label> <source>, …` |

A key missing from a run's manifest (older runs) is skipped. Runs whose `output_schema_version`
majors differ are refused with `compare: output schema 1.x and 2.x differ; rerun one with the same
major` (exit code 1); `--allow-cross-major` compares them anyway and keeps the warning.

## What each version number means

From [ADR-0060](../adr/ADR-0060-release-versioning-and-image-keys.md) §3. The deciding question is
whether your script or comparison breaks.

| Change | Bump |
|---|---|
| a column, table, manifest key, flag token or config key removed, renamed, retyped or changed in meaning; a CLI option removed; `output_schema_version` major; a DB pin whose old DB no longer works | **major** (2.0) |
| a tool version or image key change; a DB pin bump; a new column, table, flag, option or module; a changed default; any change to a default run's output bytes on the golden and TE7 runs; a calibration rerun that moves a gate number | **minor** (1.1) |
| a bug fix, crash fix, docs, performance with byte-identical goldens and TE7 | **patch** (1.0.1) |

`output_schema_version` (in `run_manifest.json` and `datapackage.json`) is `MAJOR.MINOR` and
starts at `1.0` with 1.0.0. A new column bumps its minor; its major equals the package major. The
config's `schema_version` is the config format and is unrelated.

So: a patch upgrade gives the same outputs; a minor one can change numbers, and `compare` says
which tools or databases moved; a major one can break scripts that read the tables.

## Roll back

```bash
pip install giap==1.x.z
giap setup --db-dir /path/to/db --images --profile apptainer
```

Every digest named by a released `containers.config` or `images.json` is kept in the registry
(ADR-0060 §2), so an older wheel's pins still resolve. Its databases are its own pins: run
`giap setup --check` after a rollback too, or keep one DB directory per version when a database
pin changed.
