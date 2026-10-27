# M11-09b: bare image references in `containers.config` (no `docker://` prefix)

Milestone M11 · **pass 2** · blocks: M11-09 (`docker-smoke`, draft PR #352) and the README
quickstart's `--profile docker` · depends on: M11-06 (merged; `update_digests.py`), M11-03 (merged;
readiness H4) · estimate ≈ 40 code lines + ≈ 15 lines of lead edits to existing tests + ≈ 150 new
test lines + ADR · not `stats` · roles: implementer (config, scripts, `src/gippa/setupdb/run.py`,
changelog fragment), lead (edits to existing tests, ADR acceptance), test-author (new tests
`tests/unit/test_m11_09b_*.py`, `tests/integration/test_m11_09b_container_resolution.py`), docs
(one sentence in `docs/get-started/install.md`) · ADR: [ADR-0059](../adr/ADR-0059-bare-image-refs.md)
(Proposed) · Khoa: no (Q5).

## Problem

Every pinned process fails under `-profile docker`. `workflow/conf/containers.config` writes each
image as `docker://<image>@sha256:<digest>`, and Nextflow hands that string to `docker run`
unchanged. The Docker CLI does not accept a `docker://` transport. The M11-09 implementer saw this
with a logging fake `docker`, which received
`docker run … docker://quay.io/biocontainers/miniprot@sha256:…`.

## Evidence (official Nextflow docs and source only, R6)

The pinned Nextflow is **26.10.0**: `scripts/bootstrap_farmshare.sh` (`nextflow=${NEXTFLOW_VERSION:-26.10.0}`),
ADR-0002, GIPPA_SPEC.md §3.4. `nextflow -version` in `gippa-wf` prints `version 26.10.0 build 14055`.
Sources were read at tag `v26.10.0` of https://github.com/nextflow-io/nextflow on 2026-10-25.

1. **Docs.** `docs/container/apptainer.mdx`, line 71: "You do not need to specify `docker://` to
   pull from a Docker registry. When Apptainer is enabled, Nextflow prepends it to your image name.
   The Docker engine does not work with containers specified as `docker://`."
2. **Docker engine.** `ContainerHandler.normalizeImageName` falls back to
   `normalizeDockerImageName` (`modules/nextflow/src/main/groovy/nextflow/container/ContainerHandler.groovy`,
   lines 99–100). With no `docker.registry` set, that function returns the name unchanged
   (lines 151–159). `DockerBuilder` then appends the image verbatim to the `docker run` line
   (`DockerBuilder.groovy`, lines 185–186). GIPPA sets no `docker.registry`
   (`workflow/nextflow.config`).
3. **Apptainer adds the prefix.** `normalizeApptainerImageName` (`ContainerHandler.groovy`,
   lines 229–256) treats its input in this order:
   - a leading `/` is a path;
   - `file://` is resolved against the launch directory;
   - any `scheme://` (`IMAGE_URL_PREFIX = ~/^[^\/:. ]+:\/\/(.*)/`, line 182) is returned as is;
   - an existing file at `<launch dir>/<name>` is used as that file;
   - anything else becomes `"docker://${normalizeDockerImageName(img)}"`.

   So a bare `quay.io/biocontainers/miniprot@sha256:<d>` becomes
   `docker://quay.io/biocontainers/miniprot@sha256:<d>`, the same string the config holds today.
   The `apptainer.ociAutoPull` branch (lines 79–80) receives the same `docker://` string either way.
4. **Cache file name.** The normalized `docker://…` URL then goes to
   `ApptainerCache(...).getCachePathFor` (lines 81–84, 108–111). `ApptainerCache` extends
   `SingularityCache`, and `SingularityCache.simpleName` (`SingularityCache.groovy`, lines 88–102)
   strips everything up to `://`, replaces `:` and `/` with `-`, and appends `.img`. The cache file
   is `<cacheDir>/<simpleName>` (`localCachePath`, lines 201–203). The bare reference and the
   `docker://` reference therefore name the **same** file. For example, both give
   `quay.io-biocontainers-miniprot@sha256-<d>.img`. The six images already in
   `$NXF_APPTAINER_CACHEDIR` (`$GIPPA_ROOT/containers`, listed 2026-10-25) are not pulled again.
5. **Floor version.** `normalizeApptainerImageName`, `IMAGE_URL_PREFIX` and `simpleName` are the
   same at tag `v25.04.8`. GIPPA_SPEC.md §3.4 allows 25.04 or newer.

The bare form is also what GIPPA_SPEC.md §6.3 already writes for the container field
(`ghcr.io/<owner>/gippa-<tool>@sha256:…`), and it is the nf-core convention.

## Definition

Definitions used below:

- **Bare digest reference**: a string that fully matches
  `^[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}$`. It has no scheme, no tag, no port and no uppercase.
  Every §3.2 registry fits this pattern.
- **Reader**: any code that reads a `container = '…'` value from `containers.config`.
- **Writer**: any code that writes such a value or prints a reference meant to be pasted there.

The rules:

1. **Config.** Every `withName: 'P' { container = '…' }` value in `workflow/conf/containers.config`
   is a bare digest reference. No non-comment line contains `docker://`. The header comment (line 2)
   says the same and cites ADR-0059, for example: "Digests only: every image is
   `<image>@sha256:<digest>` with no scheme; Nextflow adds `docker://` for Apptainer (ADR-0059)."
   The seven digests and the "Tags the digests were resolved from" block do not change.
2. **Writers emit the bare form.**
   - `scripts/release/update_digests.py` `pin_line()` →
     `  withName: '<P>' { container = '<ref>' }`. Its module docstring changes to match. Its input
     contract does not change: `digests.json` values are `ghcr.io/khoa-yelo/gippa-<tool>@sha256:…`,
     and a `docker://`-prefixed value is still rejected as the wrong image.
   - `scripts/resolve_digest.py` prints `<image>@sha256:<64 hex>`. Its docstring changes to match.
     It still accepts a `docker://` prefix on its input (`split_ref` strips it).
3. **Readers accept an optional prefix** (Q3). These readers match `container = '(?:docker://)?<ref>'`,
   so an older or user-supplied config with the prefix keeps working:
   - `src/gippa/setupdb/run.py` `ToolRunner.command`;
   - `scripts/slurm/record_tools.sbatch` `image()`;
   - `scripts/slurm/fetch_bakta_db.sbatch`;
   - `scripts/slurm/record_bakta_full.sbatch`;
   - `scripts/slurm/build_envs.sbatch`.

   The cache-name derivation does not change: drop a `docker://` prefix if present, map `/` and `:`
   to `-`, append `.img`. This is `simpleName`. `gippa.preflight.environment._image_file_name` and
   `tests/helpers/fake_env.image_file_name` already behave this way and only their docstrings change.
   `gippa.annotate.tools.read_containers_config` and `gippa.io.manifest_sources` (`sha256[:-]` in the
   trace `container` column) do not need to change.
4. **`build_envs.sbatch` pulls with the scheme Apptainer needs.** It selects every pinned
   reference (bare or prefixed, `@sha256:` required) and removes duplicates (BLAST and TBLASTN
   share one image). For each reference it either prints `present: <img>` or runs
   `apptainer pull <img> docker://<bare ref>`. The `apptainer pull` CLI itself needs the
   `docker://` URI.
5. **Readiness H4 is the regression guard** (Q1). `scripts/release/readiness.py` accepts a
   `container =` value only if it is a bare digest reference.
   - A value with a digest and a `docker://` prefix fails with the reason
     `docker:// prefix: <P>` (for a `withName` pin) or `docker:// prefix: <ref>` (for any other
     `container =`).
   - Any value without a digest, prefixed or not, keeps the reason `tag, not digest: …`.
   - The row's `expected` text becomes `<image>@sha256:<64 hex>, no docker:// prefix`.
6. **Spec amendment.** GIPPA_SPEC.md is read-only. ADR-0059 amends the reference *form* in §3
   (Containers row) and §6.x (the `apptainer` profile bullet). It also supersedes ADR-0008's
   "`docker://…@sha256:…`" wording for `containers.config`. The digests, tags and pinned set stay
   as ADR-0008 records them.
7. **Changelog.** Add `changelog.d/+bare-image-refs.bugfix.md`: "`-profile docker` runs the pinned
   images: `containers.config` drops the `docker://` prefix that `docker run` rejects."

### Map of every `docker://` occurrence

Found with `grep -rIn 'docker://'` over the repository, excluding `local/`, `.venv` and `.git`:
41 files on `origin/main` at `613d653c`.

| File | Class | Change |
|---|---|---|
| `workflow/conf/containers.config` | runtime | 7 pins and the line-2 comment drop the prefix (rule 1) |
| `src/gippa/setupdb/run.py:99` | runtime | regex takes an optional prefix (rule 3). **Breaks without the change**: `ToolUnavailable` for bakta, genomad and antismash downloads |
| `src/gippa/preflight/environment.py:156` | runtime | docstring only. `removeprefix` is already a no-op on bare refs |
| `scripts/slurm/build_envs.sbatch:22,25` | runtime (Slurm) | rule 4. **Breaks without the change**: its `grep` finds no refs, so nothing is pulled and the job still exits 0 |
| `scripts/slurm/record_tools.sbatch:42` | runtime (Slurm) | rule 3. **Breaks without the change**: the image path becomes `…/containers/.img` |
| `scripts/slurm/fetch_bakta_db.sbatch:29`, `record_bakta_full.sbatch:27` | runtime (Slurm) | rule 3. **Breaks without the change**: empty `DIGEST`, "missing …img" |
| `scripts/resolve_digest.py:5,76` | tooling | writer (rule 2). Line 32 (`removeprefix` on input) stays |
| `scripts/release/update_digests.py:10,84` | CI (`containers.yml` digest PR) | writer (rule 2) |
| `scripts/release/readiness.py:44,469` | CI / release | H4 (rule 5). **Breaks without the change**: H4 fails on all 7 bare pins |
| `GIPPA_SPEC.md:271,922` | spec | read-only; amended by ADR-0059 |
| `docs/adr/ADR-0008-tool-pins.md:62` | ADR | stays; ADR-0059 supersedes that wording. On acceptance the lead appends to its Status line: "; reference form superseded by ADR-0059" |
| `docs/get-started/install.md:97` | user docs (current) | one sentence (docs role): "named after the reference in `containers.config`, with `/` and `:` replaced by `-` and `.img` appended (a `docker://` prefix, if present, is dropped first)" |
| `docs/specs/M2-01-…`, `M2-06-…`, `M8-17-…`, `M11-03-…`, `M11-06-…` | old specs (historical) | leave. M8-17's future `FGSEA_IMAGE` pin must use the bare form when it lands (note for the lead) |
| `tests/data/**/PROVENANCE.md` (5 files) | test data (historical) | leave. They record the `apptainer pull docker://…` URI actually used, which is still correct for Apptainer |
| `tests/helpers/fake_env.py:33` | test helper | docstring only, optional (lead) |
| tests (unit and integration) | tests | see Test plan |

Files that do not contain `docker://` but are affected or were checked:

- `scripts/slurm/fetch_resources.sbatch` globs `antismash-standalone@sha256-*.img` (name unchanged).
- `record_tools.sbatch:45` (`digest()`) already accepts any prefix.
- PR #352 (`docker-smoke.yml`, `tests/unit/test_m11_09_*`) already accepts both forms:
  `sed -E "s|.*container = '(docker://)?…"`.

## Examples

`<d>` is the 64-hex digest already pinned for that process.

| # | Input | Expected output |
|---|---|---|
| 1 | `containers.config`, BAKTA line | `  withName: 'BAKTA' { container = 'oschwengers/bakta@sha256:1e24ae2f…07f18' }`. `grep -v '^\s*//' containers.config \| grep -c 'docker://'` prints `0` |
| 2 | `-profile docker`, a MINIPROT task (logging fake `docker` on PATH) | the logged `docker run` argv contains the word `quay.io/biocontainers/miniprot@sha256:<d>` and no word starting with `docker://` |
| 3 | `-profile apptainer`, `apptainer.cacheDir = $C` holding `$C/quay.io-biocontainers-miniprot@sha256-<d>.img`, logging fake `apptainer` | no `apptainer pull` call. The task's `.command.run` runs `apptainer exec … $C/quay.io-biocontainers-miniprot@sha256-<d>.img` |
| 4 | `_image_file_name("oschwengers/bakta@sha256:<d>")` and `_image_file_name("docker://oschwengers/bakta@sha256:<d>")` | both give `oschwengers-bakta@sha256-<d>.img` |
| 5 | `update_digests.py --digests {"OSTIR": "ghcr.io/khoa-yelo/gippa-ostir@sha256:<e>"}` on today's config | one new line `  withName: 'OSTIR' { container = 'ghcr.io/khoa-yelo/gippa-ostir@sha256:<e>' }`; the rest as M11-06 specifies |
| 6 | `resolve_digest.py quay.io/biocontainers/blast:2.16.0--h66d330f_5` (registry HEAD answers `Docker-Content-Digest: sha256:9dfc…959c`) | prints `quay.io/biocontainers/blast@sha256:9dfc…959c`. The same input with a `docker://` prefix gives the same output |
| 7 | readiness H4 on the new config | `pass`, value `7 digests; …` |
| 8 | H4 after BAKTA's value is changed to `docker://oschwengers/bakta@sha256:<d>` | `fail`, value contains `docker:// prefix: BAKTA` |
| 9 | H4 after BAKTA's value is changed to `oschwengers/bakta:v1.12.2` | `fail`, value contains `tag, not digest: BAKTA` (M11-03 example 13 still holds) |
| 10 | `build_envs.sbatch` pull loop with fake `apptainer`, empty `$NXF_APPTAINER_CACHEDIR` | 6 calls (BLAST and TBLASTN share one), each `pull $C/<name>.img docker://<bare ref>`. A second run with the files present prints 6 `present:` lines and makes no `pull` call |

## Edge cases

- **Legacy prefixed config** (an older checkout, or a user's `-c` override that still has `docker://`):
  - readers (rule 3) still find the image and derive the same `.img` name;
  - `update_digests.py` replaces a prefixed line for a pinned process in place with a bare line and
    leaves the other lines byte-identical (the M11-06 contract);
  - H4 flags any remaining prefixed lines.
- **Mixed config** (some pins bare, some prefixed): runs under Apptainer. Under Docker, only the
  prefixed processes fail. H4 fails and names each prefixed process.
- **Docker Hub short names** (`oschwengers/bakta`, `antismash/standalone`). With no
  `docker.registry` or `apptainer.registry` set, both engines resolve them on Docker Hub. A user
  config that sets `docker.registry` or `apptainer.registry` (as nf-core pipeline configs do,
  for example `quay.io`) now gets those two refs rewritten to `<registry>/oschwengers/bakta@…`,
  which do not exist there. The `docker://` form was immune to this for Apptainer. The default keeps
  the short names, because they are what `resolve_digest.py` prints and they keep today's cache
  names. See Q2.
- **A file named like the reference in the launch directory.** Under Apptainer, Nextflow uses
  `<launch dir>/<ref>` if that path exists (`ContainerHandler.groovy`, lines 247–251). A path such as
  `quay.io/biocontainers/blast@sha256:…` will not exist in practice. Documented only.
- **`apptainer.ociAutoPull = true`.** Nextflow returns `docker://<ref>` for both forms (lines
  79–80). No difference.
- **`-profile conda` / `standard` / `test` (stub).** No container engine is enabled, so the
  container string is not used. No change.
- **Duplicate image** (BLAST and TBLASTN). One pull and one cache file, as today.
- **Empty or absent pins.** A config with no `withName … container` lines: readers return nothing,
  `build_envs.sbatch` pulls nothing and exits with the env-build status, and H4 reports the
  existing "not pinned, not conda-only" reasons. An empty `digests.json` stays an error (M11-06).
- **The `//` in comments.** `test_base_config_retry._strip_comments` comments on "the // of
  docker:// URIs". With bare refs no `//` is left on pin lines, so the helper is unaffected.
- Genomic edge cases (circular/linear, minus strand, ties) do not apply: no coordinate, sequence or
  statistic changes.

## Acceptance criteria

Each criterion is checkable by a test.

1. `workflow/conf/containers.config`: every value returned by `read_containers_config` is a bare
   digest reference, and no non-comment line contains `docker://`. The process set, the seven
   digests and the tag comment block are unchanged from `origin/main` at `613d653c` (compare the
   digests as sets).
2. Example 4: the cache names derived from the new config by `_image_file_name` are exactly the
   six names derived from the old (prefixed) config, so nothing is re-pulled.
3. Example 2: under `-profile docker`, the argv of a pinned process contains no `docker://`.
4. Example 3: under `-profile apptainer` with a seeded cache, Nextflow makes no pull call and runs
   the seeded `.img`.
5. Examples 5 and 6: the writers emit bare references, and `resolve_digest.py` still accepts a
   prefixed input.
6. Examples 7–9: H4 passes on the repository config and fails with the stated reasons.
7. Rule 3: `ToolRunner.command("bakta", …)` returns
   `apptainer exec … <cache>/oschwengers-bakta@sha256-<d>.img …` for a bare config and for a
   prefixed one.
8. Example 10: `build_envs.sbatch` pulls each distinct pinned image once, with a `docker://` URI,
   into the `simpleName` file.
9. No file under `scripts/`, `src/` or `workflow/` (excluding comments) requires the literal
   `container = 'docker://` to find a pin.
10. The existing suite passes after the lead's edits below:
    `pytest -m "not slow and not calibration_full and not e2e"`, plus the two edited integration
    tests where Nextflow is available. PR #352's `docker-smoke` job is green once its branch merges main with
    this change. That is M11-09's AC 2, unblocked here.

## Test plan

### Existing tests the lead must edit

Test-authors and implementers never edit existing tests. These edits go into this PR. Without
them the suite is red, or weakened without anyone noticing.

| Test | Line | Why it changes | Edit |
|---|---|---|---|
| `tests/unit/test_nf_modules_lint.py` | 234 | fullmatch requires `docker://` | `re.fullmatch(r"[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}", …)` |
| `tests/unit/annotate/test_tool_env_files.py` | 14, 288–289 | `DIGEST` and `startswith("docker://…")` | `DIGEST = re.compile(r"^[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}$")`; `startswith("oschwengers/bakta@")`, `startswith("antismash/standalone@")` |
| `tests/unit/annotate/test_tool_registry.py` | 159 | `startswith("docker://oschwengers/bakta@sha256:")` | `startswith("oschwengers/bakta@sha256:")` |
| `tests/unit/test_setupdb_review.py` | 380 | `_image` regex requires the prefix (`assert m is not None` fails) | `container = '(?:docker://)?([^']+)'` |
| `tests/unit/test_m11_06_update_digests.py` | 57 | `_pin_line` expects the prefix | `f"  withName: '{process}' {{ container = '{ref}' }}"` |
| `tests/unit/test_m11_06_pins_profile.py` | 80 | `findall(r"docker://(ghcr…)")` would silently match nothing | `re.findall(r"container = '(?:docker://)?(ghcr\.io/[^']+)'", CONTAINERS)` |
| `tests/unit/test_m11_03_rows.py` | 463–464 | the `re.sub` pattern no longer matches, so the `assert "bakta:v1.12.2'" in text` fails | pattern `container = 'oschwengers/bakta@sha256:[0-9a-f]{64}'` → `container = 'oschwengers/bakta:v1.12.2'` |
| `tests/unit/test_m11_03_review_fixes.py` | 476 | under rule 5 `docker://ubuntu@<digest>` is now a fail | `f"process.container = 'ubuntu@{DIGEST}'"`. Lines 460, 461 and 475 stay (a prefixed tag still fails; a comment is still ignored) |
| `tests/integration/test_nf_conda_profile.py` | 112 | `startswith("docker://")` on the flattened config | `re.fullmatch(r"[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}", container.strip("'\""))` |
| `tests/integration/test_tblastn_live.py` | 35 | the regex requires the prefix, so the assert fails rather than skipping | `r"container = '(?:docker://)?quay\.io/biocontainers/blast@sha256:([0-9a-f]{64})'"` |

These tests pass unchanged and are left as they are:

- `test_m11_06_conda_off.py:130`: duplicate detection keys on `withName`; the prefixed fixture
  line is still a duplicate.
- `test_tool_env_files.py:266` (example 8): the parser is ref-agnostic.
- `test_base_config_retry.py:28`: only a comment mentions `docker://`.
- `test_recording_provenance.py:145`: a PROVENANCE fixture; `import_recording.py` reads only
  `sha256:<hex>`.
- `manifest/test_manifest_cli_nextflow.py:116`, `test_manifest_sources.py:129`,
  `test_manifest_review_236.py:50`: the trace `container` column. `sha256[:-]` matches either
  form. Under `-profile docker` the trace now holds the bare ref, which the test-author covers
  below.
- `fake_env.py`: docstring only.

### New tests (test-author)

`tests/unit/test_m11_09b_bare_refs.py` (fast):

- AC 1, AC 2, AC 9 as static checks. For AC 2, the "old" config is the current text with
  `docker://` restored on each pin, built in the test; no git access.
- Example 4 on `_image_file_name` and `fake_env.image_file_name`.
- Examples 5 and 6:
  - `update_digests.py` through a subprocess on a temporary copy;
  - `resolve_digest.resolve` with `urllib.request.urlopen` monkeypatched to return a
    `Docker-Content-Digest` header, imported by path with `-I` semantics;
  - a prefixed input to `split_ref`.
- Examples 7–9 through the `tests/helpers/release_readiness.py` fixture.
- AC 7: `ToolRunner.command` with a bare config and with a prefixed config, through `GIPPA_REPO`
  pointing at a temporary repository.
- One trace row with a bare `container` value gives the same `container_digest`.

`tests/unit/test_m11_09b_build_envs_pull.py`: example 10. Run `build_envs.sbatch` with bash
(not sbatch) in a sandbox like `annotate/test_build_envs_review.py`: a throw-away `GIPPA_REPO`,
stub `scripts/build_envs.sh` that exits 0, and a fake `apptainer` that logs its argv and touches
the target. Assert six distinct `pull <C>/<simpleName>.img docker://<bare ref>` lines, then
`present:` lines on the second run. Also run it against a prefixed config to check the
optional-prefix reading.

`tests/integration/test_m11_09b_container_resolution.py` (`needs_nextflow`; no Docker or Apptainer
needed): a one-process fixture pipeline in `tmp_path`, process `MINIPROT`, a `stub:` block that
touches its output, and `includeConfig` of the real `workflow/conf/containers.config`.

- (a) Example 2: `docker.enabled = true` and a logging fake `docker` first on PATH. The fake runs
  the remaining command after the image argument, or exits 0. Assert on the logged argv.
- (b) Example 3: `apptainer.enabled = true`, a seeded `apptainer.cacheDir`, and a logging fake
  `apptainer` (`exec <img> <cmd…>` runs `<cmd…>`). Assert on `.command.run` and that no `pull`
  call was logged.

Both runs use `-stub-run` and the local executor and stay in the session allocation (under a
minute). If Nextflow on PATH is not 26.10.x, the test records the version it ran under in its
output.

## Size

Code ≈ 40 lines:

- `containers.config`: 8
- `update_digests.py`: 2
- `resolve_digest.py`: 2
- `readiness.py`: ≈ 10
- `setupdb/run.py`: 1
- four sbatch scripts: ≈ 8
- docstrings: 3
- changelog fragment: 2

Lead test edits ≈ 15 lines. New tests ≈ 150 lines. ADR ≈ 70 prose lines. The PR is well under
~400 lines. Single PR, because the config change and the lead's test edits must land together.

## Questions for the lead (reversible defaults in use)

- **Q1. H4 strict or tolerant?** Strict means a `docker://` prefix on any `container =` fails H4,
  so the Docker breakage cannot come back unnoticed. Default: **strict**. Reversible: one regex
  in `readiness.py`.
- **Q2. Docker Hub short names, or `docker.io/…`?** Fully qualified
  `docker.io/oschwengers/bakta@…` and `docker.io/antismash/standalone@…` would be immune to a
  user's `docker.registry`/`apptainer.registry`. But they change the cache names to
  `docker.io-oschwengers-bakta@…img` and `docker.io-antismash-standalone@…img`. That re-pulls Bakta
  and the ~9 GB antiSMASH image on FarmShare, and needs edits to the hard-coded names in
  `fetch_bakta_db.sbatch`, `record_bakta_full.sbatch` and `fetch_resources.sbatch` and to
  `test_tool_registry.py`. Default: **short names** (no re-pull). The registry caveat goes in
  `docs/user-guide/platforms.md` under M11-09.
- **Q3. Readers tolerate an optional `docker://`, or accept bare only?** Default: **tolerate**
  (rule 3). User `-c` overrides and older configs keep working, and H4 guards the repository file.
- **Q4. `install.md`.** Edit its one sentence (default, docs role, in this PR), or leave it? It
  stays literally true either way.
- **Q5. Does ADR-0059 need Khoa?** Default: **no**. §6.3 already writes the bare form, and §3 and
  §6.x still describe how Apptainer runs the images. The change only fixes the profile the spec
  keeps "for users elsewhere". If you read it as an interface change, mark it Needs Khoa: yes and
  batch it with the other Khoa questions. It does not block: the fix is reversible.
- **Q6. ADR number.** ADR-0059 was the next free number on `origin/main` and in local worktrees
  on 2026-10-25. Renumber if a concurrent branch takes it first.
- **Q7. Historical files.** Leave the five `tests/data/**/PROVENANCE.md` files and the M2-01,
  M2-06, M8-17, M11-03 and M11-06 specs unchanged (default). They record what was done at the time.
