# M11-14: container-complete install with versioned artifacts

Milestone M11 · **pass 2, after the rebrand** · depends on: M11-01 (wheel bundles `workflow/`),
M11-04 (`release.yml`), M11-06 (`containers.yml`, `update_digests.py`), M11-09b (bare refs,
ADR-0059), M11-03 (readiness H4); the GIAP → GIPPA rename PR (backlog "Rebrand") before any
publish · estimate ≈ 1,150 code lines + ≈ 1,300 test lines in 8 code PRs (14a–14h) and one
docs PR (14i), plus one Khoa action (14j) · not `stats` (the operon fallback in 14d reuses the
existing distance model; no statistic changes) · roles: implementer (`scripts/release/`,
`.github/workflows/containers.yml`, `src/giap/preflight/`, `src/giap/setupdb/`, `src/giap/io/fgsea_r.py`,
`src/giap/compare/`, `workflow/`), test-author (`tests/unit/test_m11_14_*.py`), docs
(`docs/get-started/install.md`, `docs/get-started/upgrade.md`, `docs/release.md`) · versioning
policy: [ADR-0060](../adr/ADR-0060-release-versioning-and-image-keys.md) (proposed) · **Khoa: yes**
(the first GHCR publish, Q-K13–Q-K19).

Khoa's direction (2026-10-25): "`pip install` + preflight passes" with every tool in a container,
and "iterate on the final version fast. Make sure that whatever is built can be iterated to have
1.1, 1.2, 2.0, and so on."

## Spec sections implemented

- §3 Containers row ("referenced by digest"), §6.3 container field
  `ghcr.io/<owner>/<pkg>-<tool>@sha256:…`, §6.x profiles (`docker`, `apptainer`, "for users
  elsewhere"), G19.
- §7.0 Environment rows (`environment.tools`: "binaries, environments and DB versions of enabled
  modules found → else error (exit 2)").
- §9.2 Versioning row: "`schema_version` in the manifest; semantic versioning; breaking column
  changes only in major versions"; §12 Releases: "semantic versioning … containers pinned by
  digest"; §12 `containers.yml`; §14.6 "`containers.config` by digest"; §17 "prefer `ghcr.io`
  images built by `containers.yml`".
- §3.10 and AGENTS.md: agents never push tags and never publish; only Khoa publishes.
- R5/R11 (UniOP unlicensed: never redistributed), R6 (official sources only), R14 (cluster caps).
- Backlog: "Off-FarmShare install (release blocker, Khoa question)" (this issue answers it) and
  "Rebrand" (images are published only under the new name).

## The target 1.0 user path

```text
1. pip install gippa                       # Python package; workflow/ and containers.config inside
2. Docker or Apptainer, Nextflow ≥ 25.04.8 + Java ≥ 17 on PATH   (prerequisites, Q-K17)
3. gippa setup --db-dir DB [--images --profile apptainer]   # databases (+ pre-pulled images)
4. gippa preflight ... --profile docker|apptainer            # passes
5. gippa run ... --profile docker|apptainer
```

No repository checkout, no `build_envs.sbatch`, no `$GIAP_ROOT/micromamba`. Nextflow pulls every
tool image by digest (docker), or runs it from the cache that step 3 filled (apptainer). Below,
`<pkg>` is the PyPI project name read from `pyproject.toml` (`giap` today, `gippa` after the
rename); `<owner>` is the GitHub owner (`khoa-yelo`, Q-K14).

## Definition

### 1. Image coverage (every non-in-process tool)

| Tool | Image today | After M11-14 | PR |
|---|---|---|---|
| bakta, transtermhp, genomad, antismash, miniprot, blast (BLAST, TBLASTN) | upstream, by digest | unchanged | none |
| promotech, promoter_calculator, terminatornet, bactermfinder, ostir | Dockerfile, never built (M11-06) | `ghcr.io/<owner>/<pkg>-<tool>` by digest, content-addressed tag (item 2) | 14a, 14b, 14j |
| fgsea | conda-only, started by the `function` step (M8-17 decision 9, ADR-0038) | `ghcr.io/<owner>/<pkg>-fgsea` built from `workflow/containers/fgsea/Dockerfile` that installs `workflow/envs/locks/fgsea.explicit.txt`; `giap.io.fgsea_r` runs it through the profile's runtime when no prebuilt env is found | 14c |
| uniop | conda-only; unlicensed | **never published.** Default (Q-K13): the operon layer degrades to the distance model with a recorded reason when UniOP is not installed (14d); opt-in local build `<pkg> setup --tools uniop` (14e); Khoa asks the author for a licence in parallel | 14d, 14e |
| codon_bias | in-process | unchanged | none |
| cascadia, timsconvert (`PILOT_TOOLS`), DIA-NN | pilot only | out of scope: not on the user path | none |

**fgsea (14c).** `giap.io.fgsea_r` resolves the engine in this order: `GIAP_FGSEA_PREFIX`; the
prebuilt `giap-fgsea` env (marker equal to the env file's sha256, as today); the pinned image
from `containers.config` (a normal `withName: 'FGSEA' { container = '<ref>' }` line, read by
`read_containers_config`, although no Nextflow process of that name exists) run as
`docker run --rm -v <work>:<work> <ref> Rscript …` or `apptainer exec -B <work> <cache>/<file>.img
Rscript …`, the runtime given by `giap run` in the env var `GIAP_CONTAINER_RUNTIME`
(`docker`|`apptainer`|unset). No GSEA Nextflow process is added (M8-17 decision 9 stays: the
step starts R itself). R is GPL-2 and fgsea MIT; the image label and NOTICE say so.

**UniOP (14d, 14e).** Three options, not exclusive:

| Option | What happens | Cost | Default |
|---|---|---|---|
| (a) licence from the author | publish `<pkg>-uniop` like the others | Khoa emails hongsua/UniOP; unknown delay | asked in parallel (Q-K13) |
| (b) local build | `<pkg> setup --tools uniop [--runtime docker|apptainer|micromamba]` fetches the pinned commit's GitHub archive (official, R6) on the user's machine and builds a local image (`docker build`, or `apptainer build --fakeroot` from a local `docker-daemon`/OCI archive) or a micromamba env from the bundled `uniop.yml` + locks. Tagged `localhost/<pkg>-uniop:<key>`, never pushed; the user is told UniOP has no licence | 14e ≈ 160 lines; needs root/fakeroot or micromamba on the user's side | opt-in |
| (c) degrade | config `annotation.uniop.enabled: auto` (new; `auto` \| `true` \| `false`). Preflight resolves `auto` to `true` when UniOP is installed for the profile, else `false` with warning `environment.tools: uniop not installed (unlicensed, not distributed); operons use the distance model`. `giap run` passes `params.run_uniop`; the `UNIOP` process is skipped; `operon_pairs` already falls back to the distance model (`DISTANCE_FLAG`); `operon_meta.json` reason becomes `UniOP not installed` and the manifest's `layers_missing` gains `uniop` | 14d ≈ 140 lines | **on** (the 1.0 pip path) |

`enabled: true` with UniOP missing stays an error (exit 2), so a user who asks for UniOP gets it or
a clear failure. `--publication` does not change this (Q-K13 asks whether it should require UniOP).

### 2. Content-addressed, immutable image tags (14a, 14b; ADR-0060 §1)

- **Image key.** For tool `t`, the inputs are every file under `workflow/containers/t/`, plus
  `workflow/envs/t.yml`, plus `workflow/envs/locks/t.explicit.txt` and `t.pip.txt` when present.
  `key(t)` = the first 12 hex of sha256 over the UTF-8 text made of one line
  `<repo-relative POSIX path>\t<sha256 hex of the file bytes>\n` per input, sorted by path.
  A static test parses each Dockerfile's `COPY` sources and asserts each one is an input (so the
  key cannot miss a file the build reads).
- **Tag.** `<version>-<key>`, where `<version>` is the registry version with `git:<40 hex>` written
  `git-<first 8 hex>`, and every character outside `[A-Za-z0-9_.-]` replaced by `_`. Example:
  `ostir:1.1.3-3f9c0a1b2d4e`, `promotech:git-56251ad9-0b7e19c4aa21`. It fully matches the OCI tag
  rule `[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}`.
- **Locked builds.** Each Dockerfile installs from the explicit lock (and pip lock) when one exists,
  as `scripts/build_envs.sh` does, so one key gives one set of package builds. BacTermFinder has
  no lock today: 14b adds `locks/bactermfinder.{explicit,pip}.txt` (resolved on FarmShare by the
  existing lock script) or the PR states why it cannot (Q-L12).
- **Immutable.** `containers.yml` pushes only the `<version>-<key>` tag. It drops the `:<sha7>`
  and the bare `:<tool version>` tags of M11-06 (the latter is mutable). Before building, the plan
  job asks the registry (`docker buildx imagetools inspect`, anonymous, in Actions only) whether
  `<pkg>-t:<version>-<key>` exists. If it exists, the tool is not rebuilt and its existing digest
  goes into the digest PR. A tag is never pushed twice: if the push would overwrite an existing
  tag (race), the job fails `tag exists: <ref>`.
- **Only changed tools rebuild.** A GIAP release in which `ostir`'s inputs did not change keeps the
  same `ostir` tag and digest as the release before.
- **Lock file.** `workflow/conf/images.json` (new, bundled in the wheel) lists, per GHCR tool:
  `{"key", "tag", "ref" (bare digest ref), "processes", "built_from" (commit), "licence"}`.
  `update_digests.py` writes it in the same PR as `containers.config`; `scripts/release/image_key.py
  --check` recomputes every key and fails `stale image: <tool> (inputs changed since <key>)` when a
  Dockerfile, env file or lock changed without a new pin.
- **No deletion, no retagging (ADR-0060 §2).** Agents never delete or retag GHCR packages or
  versions. Retention: every digest named in a released `containers.config` (any `v*` tag) is kept
  forever. Digests that no release names may be pruned by Khoa only, not before 180 days after
  their build, and never the latest one per tool. Public GHCR storage is free, so the default is
  "prune nothing" (Q-K16).

### 3. Preflight and setup (14f, 14d)

`missing_tools(tools, profile, giap_root)` becomes `missing_tools(tools, profile, image_cache,
env_root)`:

| Tool kind | `docker` | `apptainer` | `conda` / `standard` / FarmShare profiles |
|---|---|---|---|
| pinned image (upstream or GHCR) | present (Nextflow pulls by digest) | present iff `<image_cache>/<simpleName>.img` exists | as today |
| pinned image **and** a prebuilt env | the image is used; the env is ignored | same | as today (env) |
| conda-only (none after M11-14 except a local UniOP env) | needs the prebuilt env **and** a `conda` profile; under `docker` alone → missing with `conda-only tool under a container profile: <tool>` | same | as today |
| uniop with `enabled: auto` | not an error; warning (item 1c) | same | as today |
| in-process | never missing | never missing | never missing |

- `image_cache` = `NXF_APPTAINER_CACHEDIR`, else `$GIAP_ROOT/containers`, else
  `<db_dir>/containers`, so a pip user with no `$GIAP_ROOT` has a defined place (Q-L11). `giap run`
  passes the same directory to Nextflow (`apptainer.cacheDir`).
- The "no prebuilt env needed" rule applies exactly to tools with a `withName` pin in the bundled
  `containers.config`. The current FarmShare behaviour (envs under `$GIAP_ROOT/micromamba`) is
  unchanged for the `farmshare*` profiles.
- **`<pkg> setup --images --profile apptainer|docker [--cache-dir DIR]`** (14f): for every pinned
  process needed by the enabled modules, `apptainer pull <cache>/<simpleName>.img docker://<ref>`
  (the file name Nextflow looks for, ADR-0059) or `docker pull <ref>`. It verifies the pulled
  digest equals the pin and records `{"name": "image:<tool>", "version": "<tag>", "checksum":
  "<digest>", "status": "ok"}` items in `db_manifest.json`. This is the only command that touches a
  registry. `ToolRunner` (the Bakta/antiSMASH downloaders) gains the `docker run` path, so
  `<pkg> setup` itself needs no checkout and no micromamba.
- **Database pins per release.** `src/giap/data/db_sources.yaml` (in the wheel) is the per-release
  pin of every database: version, URL, checksum. `<pkg> setup --check` reports an installed item
  whose version or checksum differs from this release's pin as `mismatch` with
  `installed <v_old>, this release pins <v_new>` (today's status set; no new status). An unchanged
  pin across releases means the DB directory is reused as is.

### 4. Versioning, release candidates and compatibility (14g, 14h; ADR-0060 §3–§6)

**What bumps what** (ADR-0060 §3; the decisive column is "a user's script or comparison breaks"):

| Change | Bump |
|---|---|
| a column, table, manifest key, flag token or config key is removed, renamed, retyped or changes meaning; a CLI option removed; `schema_version` major; a DB pin whose old DB no longer works | **major** (2.0) |
| a tool version or image key change; a DB pin bump; a new column, table, flag, option or module; a default changed; any change that alters a default run's output bytes on the golden/TE7 runs; a calibration re-run that moves a gate number | **minor** (1.1) |
| a bug fix, crash fix, docs, performance with byte-identical goldens and TE7, an image rebuild with the same key (none: keys are content-addressed, so that cannot happen) | **patch** (1.0.1) |

- `schema_version` (manifest and tables) is `MAJOR.MINOR`: a new column bumps the minor, a
  breaking change bumps the major, and the schema major equals the package major (2.x writes
  `schema_version: "2.<n>"`). Today `schema_version` is the **config's** field copied into the
  manifest (`src/giap/config/model.py:362`); 14h separates `output_schema_version` (tables +
  manifest) from the config's (Q-L13).
- A static test (14h) compares the column dictionaries in `src/giap/schemas/tables/` against the
  previous release tag's copy (read with `git show v<prev>:…`, no network) and fails when a column
  was removed or renamed without a major bump in `__version__`, or a column was added without a
  schema minor bump.

**Release candidates** (ADR-0060 §4). `1.1.0rc1`, `1.1.0rc2`, … are cut on main by the same
release-prep PR as M11-13 step 1. Agents dispatch only `dry_run=true`; Khoa dispatches
`dry_run=false target=testpypi` for an rc (rc versions are uploaded to TestPyPI only, never to
PyPI; Q-K15). The final `1.1.0` follows on the same commit as the last rc plus the version line.

| Gate (M11-03 rows) | patch | minor | major |
|---|---|---|---|
| unit/lint/types, package smoke (CI) | rerun | rerun | rerun |
| `image_key.py --check`, H4, static schema diff test | rerun | rerun | rerun |
| calibration (nightly subset + attacks; `calibration_full`) | reuse if R..C touches only `NON_STATS_PATHS` | rerun (subset); full if `stats`-labelled code, `calibration.yaml` or a tool feeding a statistic changed | full rerun |
| weekly e2e (K-12, Nissle, TE8) | rerun one (K-12 standard) | rerun | rerun |
| A2 pilot gate (TE1, mapping FDR) | reuse | rerun if `mapping/`, `novel/`, DB or tool pins changed; else reuse | rerun |
| benchmark k12_standard window | reuse | rerun | rerun |
| Khoa approval (§14.8 item 5) | yes | yes | yes |

`NON_STATS_PATHS` (new constant in `readiness.py`, 14g): `docs/`, `changelog.d/`, `README.md`,
`src/giap/cli.py`, `src/giap/preflight/`, `src/giap/setupdb/`, `src/giap/runner/`,
`src/giap/report/`, `src/giap/viz/`, `scripts/release/`, `scripts/docs/`, `.github/`. Calibration
evidence from commit C is reused for R when `git diff --name-only C R` lies entirely in it (today
`CALIB_PATHS` = all of `src`, so any `src` change reruns calibration). Readiness prints
`reused from <C7> (diff ⊆ NON_STATS_PATHS)` in the row.

**Compatibility check** (14h). The manifest already records `giap_version`, each tool's
`container_digest` and the databases. `<pkg> compare` adds `compare_versions` to
`CMP/report/compare.html` and the event log:

| Runs differ in | Effect |
|---|---|
| package patch only | nothing |
| package minor | warning `compare.cross_version: runs from 1.0.2 and 1.1.0` |
| any tool `container_digest` or version | warning `compare.tool_versions: ostir 1.1.3-3f9c… vs 1.1.4-77a0…` (one line per tool) |
| database version | warning `compare.db_versions: bakta_db 6.0 vs 6.1` |
| output schema major | error, exit 1: `compare: output schema 1.x and 2.x differ; rerun one with the same major` (override `--allow-cross-major`, which keeps the warning) |

### 5. The release flow (14a, 14i; ADR-0060 §5)

```text
feature PR edits workflow/containers/<t>/ or envs/<t>.yml (or locks)
   └─ PR CI: image_key.py --plan prints "will rebuild: <t> (<old key> → <new key>)"; no push
merge to main
   └─ containers.yml (push trigger): plan → keys whose tag is absent in GHCR
        ├─ vars.GHCR_PUBLISH != 'true'  → stop after plan ("publishing off"); nothing is pushed
        └─ build (environment `ghcr`, Khoa approves until Q-K18 says otherwise) → push <version>-<key>
             └─ digest PR ci/digests-<sha7>: containers.config + images.json (+ nextflow.config selector)
                  └─ FarmShare Apptainer smoke (M11-06 item 5) → job IDs in the PR → merge (waitmerge)
release-prep PR (version, changelog, CITATION) → evidence by the table above → readiness
   └─ release.yml dry_run=true (lead): verify adds `image_key.py --check` and an anonymous
      `imagetools inspect` of every GHCR pin (fails `not public: <ref>`)
        └─ Khoa: approval → testpypi (rc) → pypi (final) → tag v<v> by the workflow
```

- **Agents never push tags, never dispatch `release.yml` with `dry_run=false`, never set
  `GHCR_PUBLISH`, never approve the `ghcr` environment, never change package visibility, and never
  delete packages.**
- The **first** `containers.yml` publish is a Khoa action (14j): after the rename is merged, Khoa
  sets `GHCR_PUBLISH=true`, dispatches `tools=all`, approves the `ghcr` environment and makes the
  packages public (Q-K8). Until then, `containers.yml` runs its plan job only.
- The wheel ships `workflow/conf/containers.config` and `images.json` (M11-01 bundles `workflow/`),
  so each PyPI release pins exact digests forever.

### 6. Sequencing with the rebrand

- 14a–14h land in either order relative to the rename PR (code only, nothing published). The rename
  PR's script also renames the image prefix (`giap-` → `gippa-`) in `containers.yml`,
  `update_digests.py`, the Dockerfile labels and the `IMAGE_NS`-derived names.
- `containers.yml` derives the prefix from `pyproject.toml` `[project] name`, and a static test
  asserts every GHCR pin in `containers.config`/`images.json` starts with
  `ghcr.io/<owner>/<project name>-`. So nothing can be published as `giap-*` after the rename,
  and no `giap-*` pin can survive it.
- 14j (the first publish) is gated on "rename merged" and `GHCR_PUBLISH=true` (Khoa). No image is
  ever published under `giap-`.

### 7. Docs (14i)

- `docs/get-started/install.md` rewritten for the five-step path above: prerequisites table
  (Python, Nextflow, Java, Docker or Apptainer), `pip install`, `setup --images`, preflight, run;
  FarmShare from a checkout becomes a second section; the "Supported platform" warning is replaced
  by the tested-platforms table (M11-09). UniOP: what is lost, and how to opt in (14e), with the
  licence warning.
- `docs/get-started/upgrade.md` (new; nav after install): `pip install -U <pkg>==1.y`; which
  images change (`images.json` diff, shown by `<pkg> setup --images --check`); `setup --check` for
  DB pins; comparing runs across versions (the compare warnings); the bump table from ADR-0060;
  rollback (`pip install <pkg>==1.x`: its own pins still resolve, item 2).
- `docs/release.md`: step 5 rewritten for `GHCR_PUBLISH` and the `ghcr` environment; a new
  "Cutting 1.x.y" section with the gate table; the backlog "Off-FarmShare install" entry closed.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `image_key.py ostir` on a tree where `workflow/containers/ostir/Dockerfile`, `workflow/envs/ostir.yml`, `locks/ostir.explicit.txt`, `locks/ostir.pip.txt` exist | 12 hex key over 4 sorted `path\tsha256` lines; tag `1.1.3-<key>` |
| 2 | the same tree with one byte changed in `ostir.yml` | different key; `--check` against the old `images.json` exits 1 `stale image: ostir (inputs changed since <old key>)`; other tools unaffected |
| 3 | a comment added to `workflow/envs/bakta.yml` (upstream-image tool) | no GHCR key changes; `--plan` prints `will rebuild: none` |
| 4 | `promotech` registry version `git:56251ad9b883…` | tag `git-56251ad9-<key>` |
| 5 | `containers.yml` plan where `ghcr.io/khoa-yelo/gippa-ostir:1.1.3-<key>` exists (fake inspect in the test) | `ostir` not built; its existing digest goes into `digests.json` |
| 6 | `missing_tools([ostir], "docker", …)` with an `OSTIR` pin and no env | `[]` |
| 7 | `missing_tools([ostir], "apptainer", cache=fake)` with the pin and no `.img` in the fake cache | `["ostir"]`; message names `<pkg> setup --images --profile apptainer` |
| 8 | preflight, `annotation.uniop.enabled: auto`, profile `docker`, no UniOP env | exit 0; warning `uniop not installed …`; resolved `params.run_uniop = false`; after the run `layers_missing` contains `uniop`, `operon_meta.json` reason `UniOP not installed` |
| 9 | the same with `enabled: true` | exit 2, `environment.tools: uniop` missing |
| 10 | `compare` of a 1.0.0 run and a 1.1.0 run whose OSTIR digests differ | exit 0; warnings `compare.cross_version` and `compare.tool_versions: ostir …` |
| 11 | `compare` of `schema_version` 1.3 and 2.0 runs | exit 1 unless `--allow-cross-major` |
| 12 | readiness for 1.0.1 where `R..C` touches only `src/giap/cli.py` and `docs/` | calibration rows `pass (reused from <C7>)` |
| 13 | the same with `src/giap/stats/fdr.py` in the diff | calibration rows `pending` until rerun on R |
| 14 | `update_digests.py` with `{"UNIOP": "ghcr.io/…/gippa-uniop@sha256:…"}` | exit 1 (M11-06 example 5 still holds) |

## Edge cases

- **Two processes, one image** (`PROMOTECH`/`PROMOTECH_NOVEL`, `OSTIR`/`OSTIR_NOVEL`,
  `BLAST`/`TBLASTN`): one key, one digest, both lines; `images.json` lists both processes.
- **Key unchanged, upstream source changed**: impossible for the four git tools (the archive URL
  names the pinned commit, which is in the Dockerfile and so in the key). For pip/conda packages the
  locks are in the key. The base image is pinned by digest inside the Dockerfile.
- **Rebuilding the same key gives a different digest** (non-reproducible layers, timestamps): the
  tag already exists, so no rebuild happens; the first digest wins forever.
- **A tool fails to build** (Promotech's Python 3.6 env, §17): that tool keeps its previous pin (or
  stays conda-only on FarmShare); the digest PR says so; `--check` keeps failing until it builds,
  and readiness H4 is `fail` (no release with a stale pin).
- **Private package** (Q-K8 not done): release.yml verify fails `not public: <ref>`; the FarmShare
  smoke reports `private package`.
- **Offline compute nodes** (HPC): `setup --images` runs where the network is; preflight checks
  the cache, not the registry.
- **Docker on arm64 (Apple Silicon)**: images are `linux/amd64` only; Docker runs them emulated and
  slowly; the platforms page says so (M11-06 question, unchanged).
- **A user with an old `$GIAP_ROOT/micromamba` env and a new image**: under container profiles the
  image wins (table in item 3); under `farmshare*` the env wins, as today.
- **Empty input** (no enabled tool needs an image, e.g. every layer precomputed): `setup --images`
  pulls nothing and records nothing; preflight passes.
- **GHCR outage**: R6 (retry once later, then PROGRESS › Blocked).
- Genomic edge cases (circular/linear, minus strand, ties): not applicable.

## Acceptance criteria

All tests are offline: no test pulls from or queries a registry; `imagetools inspect`, `docker
pull` and `apptainer pull` are replaced by fakes on `PATH` or injected callables; preflight tests use
temporary fake `GIAP_ROOT`/cache/DB directories.

1. `image_key.py`: examples 1–4 on a temporary tree; keys are stable across runs and independent of
   file mtime and directory listing order; a Dockerfile `COPY` of a file outside the inputs fails
   the static test.
2. `containers.yml` static test: only the `<version>-<key>` tag is pushed; the plan job skips
   existing tags (example 5 through a fake inspect script); `build` needs `vars.GHCR_PUBLISH ==
   'true'` and the `ghcr` environment; the image prefix is derived from `pyproject.toml`; no `uniop`;
   every third-party action is SHA-pinned; nothing pushes to main.
3. `update_digests.py` writes `images.json` and `containers.config` consistently (idempotent;
   M11-06 examples 1–5 still pass); `image_key.py --check` passes on the repo after the digest PR.
4. A static test: every GHCR pin starts with `ghcr.io/<owner>/<project name>-`; every
   non-in-process tool in `TOOLS` other than `uniop` has a `withName` pin (after 14j merges; before
   that the test lists the expected-missing set explicitly and must shrink only).
5. Preflight: the item 3 table, one test per cell, on fake roots; examples 6–9.
6. `setup --images`: with fake `apptainer`/`docker` binaries, the file names equal Nextflow's
   `simpleName` (ADR-0059), a digest mismatch fails `pulled <d1>, pinned <d2>`, the manifest items
   are written; `setup --check` reports a DB pin change as `mismatch` with both versions.
7. fgsea: `fgsea_r` resolution order (prefix → env → image) with fakes; the `docker run` and
   `apptainer exec` argv are exact; the GSEA golden (`golden-r.yml`) is unchanged when run through
   the image on FarmShare (one Slurm job, ID in the PR).
8. UniOP: examples 8–9; the stub-run integration test with `run_uniop=false` produces
   `layers_missing: [uniop]` and a distance-model `operon_pairs.tsv`; `setup --tools uniop` with a
   fake runtime builds `localhost/<pkg>-uniop:<key>` and never calls `push`.
9. Compare: examples 10–11 on two fake run directories whose manifests differ only in the listed
   fields.
10. Readiness: examples 12–13; H4 fails on a stale key; `NON_STATS_PATHS` is a constant with a test
    that every listed path exists.
11. Schema diff test (14h): a removed column with an unchanged major fails; an added column with
    an unchanged schema minor fails; run against fixture schema trees, not the real tag.
12. Docs: `mkdocs build --strict -d <scratch>` passes with `upgrade.md` in the nav; the install
    page's commands match the CLI options (the generated CLI reference, M11-07).
13. After 14j (Khoa): on a FarmShare node with only the wheel installed in a fresh venv (no
    checkout on `PYTHONPATH`), `setup --images --profile apptainer`, `preflight` and a 100-kb
    `run` pass in one Slurm job; job ID in PROGRESS.

## Estimate and PR split (each ≤ ~400 lines of code + tests counted separately, ADR-0007)

| PR | Scope | Code | Tests | Depends on |
|---|---|---|---|---|
| 14a | `scripts/release/image_key.py` (key, tag, `--plan`, `--check`), `workflow/conf/images.json`, `update_digests.py` writes it | 170 | 200 | M11-06 |
| 14b | `containers.yml`: key tags only, skip-existing, `GHCR_PUBLISH` + `ghcr` environment, prefix from pyproject; Dockerfiles install from locks; BacTermFinder locks | 140 | 120 | 14a |
| 14c | fgsea Dockerfile; `fgsea_r` image path; `GIAP_CONTAINER_RUNTIME` from `giap run`; NOTICE line | 130 | 150 | 14a |
| 14d | `annotation.uniop.enabled` (config model + defaults.yaml + schema), preflight resolution, `params.run_uniop`, workflow skip, manifest/meta reason | 140 | 170 | none |
| 14e | `setup --tools uniop` local build (docker / apptainer fakeroot / micromamba), never pushed | 160 | 130 | 14d |
| 14f | `missing_tools` new signature and table; image-cache resolution; `setup --images`; `ToolRunner` docker path | 190 | 220 | 14a |
| 14g | readiness: `NON_STATS_PATHS` reuse, H4 key check; `release.yml` verify: `--check` + anonymous inspect | 110 | 130 | 14a |
| 14h | compare version warnings; `output_schema_version`; schema diff test | 110 | 180 | none |
| 14i | docs: install, upgrade, release.md, backlog entry closed; changelog fragments | 0 (≈ 350 prose) | 20 | 14a–14h |
| 14j | **Khoa**: after the rename, `GHCR_PUBLISH=true`, dispatch `tools=all`, approve, make public; then the digest PR and AC 13 | 0 | 0 | rename, 14b |

Total ≈ 1,150 code, ≈ 1,320 test lines, ≈ 350 prose. Waves: {14a, 14d, 14h} → {14b, 14c, 14f,
14g, 14e} → 14i → (rename) → 14j.

## Questions

Numbering continues M11-index (Q-K12, Q-L10). Each has a reversible default the PRs use.

| # | Question | Default in use | PR |
|---|---|---|---|
| Q-K13 | UniOP for pip users: (a) ask the author for a licence, (b) opt-in local build, (c) degrade to the distance model when absent? And should `--publication` require UniOP? | all three: (c) on by default (`enabled: auto`), (b) opt-in, Khoa sends (a); `--publication` does not require UniOP but the report states the operon source | 14d, 14e |
| Q-K14 | GHCR owner after the rename: `ghcr.io/khoa-yelo/gippa-<tool>` or an organisation? | `khoa-yelo` (same as the repo owner today) | 14b, 14j |
| Q-K15 | Release candidates: TestPyPI only, or also PyPI pre-releases (`pip install --pre`)? | TestPyPI only; PyPI gets finals only (M11 Q-K4 for 1.0 extended to every release) | 14g, 14i |
| Q-K16 | GHCR retention: prune unreleased digests after 180 days, or never prune? | never prune (public storage is free); the 180-day rule exists only if Khoa chooses to prune | ADR-0060 |
| Q-K17 | Nextflow and Java are not pip-installable. Accept them as prerequisites (step 2), or add `<pkg> setup --nextflow` that downloads the pinned launcher from the official release (still needs Java)? | prerequisites, documented with a micromamba one-liner; no download command in 1.0 | 14i |
| Q-K18 | Every `containers.yml` build needs Khoa's approval of the `ghcr` environment, or only the first one (Khoa then removes the reviewer)? | every build, until Khoa removes the reviewer (it is a publish) | 14b |
| Q-K19 | fgsea in an image (14c) vs. keeping M8-17 decision 9 (conda-only; GSEA `layers_missing` for pip users)? | image (14c); GSEA otherwise breaks "preflight passes" for `--abundance` runs | 14c |
| Q-L11 | Image cache without `$GIAP_ROOT` | `NXF_APPTAINER_CACHEDIR` → `$GIAP_ROOT/containers` → `<db_dir>/containers` | 14f |
| Q-L12 | BacTermFinder has no lock files; resolve one on FarmShare? | yes (14b), with the build log in the PR | 14b |
| Q-L13 | `schema_version` is today the config's field copied into the manifest; split out `output_schema_version`? | yes, starting at `1.0`; the config's `schema_version` keeps its meaning (config format) | 14h |
| Q-L14 | Keep the M11-06 `:<sha7>` tags for traceability? | no: `images.json` `built_from` and the OCI `revision` label record the commit | 14b |
| Q-L15 | `NON_STATS_PATHS` contents (calibration reuse) | the list in item 4; the statistician reviews it in 14g | 14g |

### Gaps found while writing (for the lead)

- Under `-profile docker` alone, a conda-only process has its `conda` directive ignored
  (`conda.enabled` is off), so it runs on the host without its tool. Preflight today reports it
  present when the env exists. Item 3 makes that a `missing` with a reason; until 14f merges, the
  docker profile silently breaks for any conda-only tool.
- The M11-06 Dockerfiles solve `workflow/envs/<tool>.yml` at build time instead of installing the
  locks, so an image built today and one built next month can differ under the same env file.
  14b fixes this before the first publish.
- `GIAP_STEP` processes run on the host (no container) with the host's `giap`. On a Slurm
  executor every node needs the same Python environment. The install page must say so (Q-L8 keeps
  "no core image").
- The distance-model fallback for operons has not been checked against the TE and calibration
  gates as a *default* path. If (c) becomes the 1.0 default for pip users, the statistician should
  confirm that the operon gate rows still hold or are reported as "distance model" (14d review).
