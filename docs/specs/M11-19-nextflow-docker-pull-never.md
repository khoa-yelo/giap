# M11-19: Nextflow's docker step images are never pulled at run time (N6b)

Milestone M11 · pass 2 · depends on: M11-14 (14f1 preflight image rule, 14f2 `setup --images`,
14p docker profile parity), #465 (N6: `--pull=never` in the UniOP probe and the R fgsea engine;
`_docker_image_exists`), M11-09 (`docker-smoke`) · estimate ≈ 80 code lines + ≈ 380 test lines
+ ≈ 40 prose lines, one PR (19a) · roles: test-author (new `tests/unit/test_m11_19_*.py`,
`tests/integration/test_m11_19_*.py`), lead (the existing tests that pin the old rule, see
"Existing tests"), implementer (`workflow/nextflow.config`, `src/gippa/preflight/environment.py`,
docs, changelog) · source: docs/backlog.md, 2026-10-28 (lead), "#453 review notes", N6 and its
follow-up N6b (named in #465's `missing_tool_reasons` docstring).

## Spec sections implemented

- §3.2 (container images "pulled by digest"), §6.3 (the Nextflow profiles), §7.0 Preflight
  ("every evaluated rule yields one check"; an environment error stops `gippa run` with exit 2
  before Nextflow starts), §12 ("containers pinned by digest"), §16 risk row "Docker Hub pull
  rate limit".
- **Policy (lead, 2026-10-28): only `gippa setup` may pull images.** `gippa run` and
  `gippa preflight` never pull: they use what `setup` put on the machine, or stop and say which
  `setup` command fetches it.

## The gap

#465 gave the two `docker run` calls GIPPA itself makes `--pull=never`, and under docker
preflight's fgsea row now checks the pinned `images.json` ref with `docker image inspect`. The
step images that **Nextflow** starts are not covered:

1. `workflow/nextflow.config` › `profiles.docker` sets
   `docker.runOptions = (params.db_dir ? "-v <db>:<db>:ro " : '') + '-u $(id -u):$(id -g)'`.
   With no `--pull` flag, Docker's default policy (`missing`) applies: when a pinned image is
   absent, `docker run` pulls it in the middle of the run (antiSMASH is ≈ 9 GB), from a task on
   whatever node runs it, and hits registry rate limits or a disk-full there instead of in
   `setup`.
2. `missing_tool_reasons` (`src/gippa/preflight/environment.py`, the `ref is not None` branch)
   checks a `containers.config` pin only under apptainer (`<image_cache>/<simpleName>.img`).
   Under docker a pinned tool is always "present" (install.md: "present (Nextflow pulls by
   digest)"), so preflight passes on a machine that has none of the images.
3. `gippa setup --images --profile docker` already pulls every pinned ref the configuration's
   modules need and checks the pulled digest (`src/gippa/setupdb/images.py` `_pull`), so the
   fetch step exists; only the run side still pulls.

## Changes

### 1. `--pull=never` in the docker profile (workflow)

`workflow/nextflow.config` › `profiles.docker`:

```groovy
docker.runOptions   = '--pull=never ' + (params.db_dir ? "-v ${params.db_dir}:${params.db_dir}:ro " : '') + '-u $(id -u):$(id -g)'
```

- `--pull=never` is the **first** token, exactly once, in the equals form (`--pull=never`, not
  `--pull never`: one token, so no option-value split can separate it). The expression stays on
  **one line** (`test_m11_14_profile_parity.py` reads `runOptions = (.+)` from one line).
- Nextflow's `DockerBuilder` puts `runOptions` (and the per-process `containerOptions`, the ps
  shim) before the image name of every `docker run` it writes into `.command.run`, so every
  containerised task of the docker profile carries it: the seven pinned processes and the local
  UniOP image of `local_tools.config` (whose docker block sets only `container`).
- Nextflow issues no `docker pull` of its own for the docker engine (it pulls only for
  Apptainer/Singularity/Charliecloud), so after this line no Nextflow code path pulls under
  docker.
- Tasks with no `container` (the conda-only tools under `docker,conda`, every `GIPPA_STEP`) run
  on the host and never call docker: unchanged.
- `apptainer`, `conda`, `standard`, `farmshare*`, `slurm`, `test`, `stub` blocks: unchanged.

### 2. Preflight: every pinned ref the run uses must be local under docker

In `missing_tool_reasons`, the `ref is not None` branch (a tool pinned in the bundled
`containers.config`) gains a docker arm. **Runtime docker** means
`container_runtime(profile) == "docker"` (`docker` in the profile list and no `apptainer`), with
or without `conda`: the conda-off selector of `container_common.config` keeps the pinned processes
in their images under `docker,conda` too, so the check applies there as well.

```python
DOCKER_MIN: tuple[int, int] = (20, 10)    # first Docker whose `docker run` has --pull
PODMAN_MIN: tuple[int, int] = (3, 0)      # Q3
SETUP_DOCKER = f"run: {_PKG} setup --images --profile docker"

def parse_docker_version(text: str) -> tuple[Literal["docker", "podman"], tuple[int, int]] | None:
    """Pure: the CLI family and (major, minor) from the first line of `docker --version`;
    None when absent or unrecognised."""
def docker_client_version() -> tuple[Literal["docker", "podman"], tuple[int, int]] | None:
    """`docker --version` (client only, never the daemon) through shutil.which("docker"),
    timeout 30 s; None when docker is not on PATH, the call fails or the output does not parse."""
def pull_never_problem(version: tuple[Literal["docker", "podman"], tuple[int, int]] | None) -> str:
    """Pure: "" when `docker run --pull=never` is supported or the version is unknown, else the
    too-old reason below."""
def docker_image_reason(ref: str) -> str:
    """"" when the image is local (`_docker_image_exists(ref)`), else the absent / no-docker /
    no-daemon reason below. Never pulls."""
```

`_docker_image_exists(tag)` stays the presence oracle, unchanged in signature and behaviour
(#465 and 14c tests monkeypatch it). `docker_image_reason` calls it first and only classifies a
`False`: `shutil.which("docker") is None` → no-docker; else one `docker image inspect <ref>` whose
stderr contains `Cannot connect to the Docker daemon`, `permission denied while trying to connect
to the Docker daemon` or `error during connect`, or which times out → no-daemon; else absent.

**`parse_docker_version`.** The first line of stdout must match
`^(Docker|podman) version (\d+)\.(\d+)` (case-sensitive; podman-docker's "Emulate Docker CLI
using podman" notice goes to stderr and is ignored). `Docker` → `"docker"`, `podman` →
`"podman"`. Anything else (`nerdctl version 1.7.6`, empty, garbage) → `None`. The Docker CLI
version is what matters: `--pull` on `docker run` is handled by the client (it decides whether to
pull after the daemon reports the image missing), so the server version is not probed.

**Order per tool, first match wins** (only under runtime docker, for a pinned tool the run needs):

| # | Condition | Reason (`<ref>` = the pinned ref as written in `containers.config`) |
|---|---|---|
| 1 | `pull_never_problem(v)` is not empty | `docker <X.Y> is older than <MIN>, the first version whose docker run has --pull=never, and runs never pull images; upgrade Docker (Docker >= 20.10, or podman >= 3.0 with podman-docker)` (`docker` reads `podman` for a podman CLI; `<MIN>` is `20.10` or `3.0`) |
| 2 | `_docker_image_exists(ref)` and `v` known and supported, or unknown | present: no reason |
| 3 | docker not on `PATH` | `image <ref> cannot be checked: docker is not on PATH (runs never pull images); install Docker >= 20.10, then run: gippa setup --images --profile docker` |
| 4 | daemon unreachable / inspect timed out | `image <ref> cannot be checked: the docker daemon is not reachable (<first stderr line, or "timed out">); start Docker, then run: gippa setup --images --profile docker if the image is missing` |
| 5 | otherwise | `image <ref> is not present locally and runs never pull it (docker image inspect <ref> failed); run: gippa setup --images --profile docker` |

`gippa` in the texts is `<pkg>` (`_PKG`), as in every existing reason. Rule 1 comes first so a
user with an old Docker is not sent to `setup` and back. A version that cannot be parsed never
blocks (Q2): such a run fails loudly at the first task (`unknown flag: --pull`, exit 125) if the
CLI really lacks the flag.

**Shared with the other docker-run tools.** The same helpers serve, under runtime docker:

- the fgsea row of #465 (a step-started tool whose `images.json` ref is docker-run with
  `--pull=never` by `gippa.io.fgsea_r`): its present/absent test is unchanged, its reason becomes
  `docker_image_reason(step_ref)` (rows 3–5; the absent text is #465's with the cause made
  precise), and rule 1 applies to it too;
- the local UniOP docker image (`uniop_local_problem`), whose task now runs with `--pull=never`
  through `runOptions`: rule 1 applies when that build is otherwise usable; its image ID check is
  unchanged.

**Probing cost.** `docker --version` runs at most once per `missing_tool_reasons` call, and only
when at least one tool needs a docker image (none for an all-conda `docker,conda` selection).
`docker image inspect` runs at most once per distinct ref (BLAST and TBLASTN share one), in
sorted tool order; after the first no-daemon result the remaining refs get reason 4 without
another call. No call ever pulls (`docker pull`, `docker run` without `--pull=never`, `docker
create`), so the `tests/helpers/container_guard.py` fake never refuses one.

**Which tools.** The tools are `required_tools(cfg, abundance)`, as today: a module that is off,
or a layer given precomputed (`annotation.precomputed.bakta_json` skips Bakta,
`precomputed.regions` skips geNOMAD and antiSMASH), needs no image and is not checked. This is the
same selection `gippa setup --images --profile docker` pulls for a configuration.

**Surfaces.** The reasons reach the `environment.tools` check (`level: error`, `exit_code: 2`,
message `missing: <tools>; <tool>: <reason>; …`, `value.missing`) through the existing
`_tools_check`; no new rule id, no schema change. `gippa run --profile docker` therefore exits 2
before Nextflow starts; `gippa preflight --profile docker` reports the same check.

### 3. Where it does not apply

| Profile list | `--pull=never` in task argv | Preflight image check under docker |
|---|---|---|
| `docker`, `standard,docker`, `docker,conda`, `slurm,docker` | yes | yes (rules 1–5) |
| `apptainer`, `farmshare,apptainer`, `slurm,apptainer`, `apptainer,docker` | no (apptainer runs the tasks) | no: the apptainer cache rule, unchanged |
| `conda`, `farmshare`, `farmshare_local`, `standard`, `slurm` | no container | no (pinned processes are not checked; unchanged) |
| `test,docker`, `test,docker,stub`, `stub,docker` | yes (the docker block applies) | no: `test`/`stub` need no tools (`_NO_TOOL_PROFILES`), unchanged |
| `test,stub` | no container | no |

`test,docker` (the `docker-smoke` workflow) therefore needs its images local:
`.github/workflows/docker-smoke.yml` already `docker pull`s each image before `nf-test`, so it
keeps passing and now also proves the tasks start without a pull. `gippa run` rejects `test` and
`stub` (`_PROFILES`), so no user run reaches that row.

Apptainer is out of scope: Nextflow pulls a missing image into `apptainer.cacheDir` and has no
"never pull" switch; preflight's cache rule is the only guard there (Q7).

### 4. podman-docker

With the `podman-docker` package, `docker` is a wrapper that runs `podman` with the same
arguments. `podman run --pull=never` and `podman image inspect` exist (Podman ≥ 3.0, Q3), and the
missing-image error of `podman run --pull=never` is `Error: <ref>: image not known` (exit 125,
already matched by `gippa.io.fgsea_r._NO_IMAGE`). Consequences:

- `docker --version` prints `podman version X.Y.Z`: rule 1 compares against `PODMAN_MIN`, never
  `DOCKER_MIN` (podman 4.9 is not "older than 20.10").
- Presence is decided by the same CLI and the same ref string the task uses, so podman's local
  short-name resolution (`oschwengers/bakta`, `antismash/standalone`) gives the same answer in
  preflight and in the task.
- Podman has no daemon: rule 4's markers are Docker's and do not occur; a podman failure falls to
  rule 5.
- Podman stays **untested** (no CI, no FarmShare podman); the docs say "best effort".

### 5. Docs and changelog (implementer)

- `docs/get-started/install.md`: the prerequisites row "Docker or Apptainer" gains "Docker 20.10
  or later (podman 3.0 or later through podman-docker: best effort, untested)"; the profile table's
  `docker` row adds "never pulls (`--pull=never`): images come from `gippa setup --images
  --profile docker`"; in "What `gippa run` needs", the pinned row's `docker` cell becomes "present
  iff `docker image inspect <ref>` succeeds (runs never pull); Docker < 20.10 is refused" and the
  `with conda` cell says "with `docker` in the list (no `apptainer`): as in the `docker` column".
  One sentence states the policy: only `gippa setup` pulls images.
- `docs/user-guide/platforms.md` › Linux + Docker and macOS: the Docker floor, "images must be
  pulled by `gippa setup --images --profile docker` first", the podman note, and that a
  `docker.runOptions` of your own Nextflow config replaces GIPPA's and must keep `--pull=never`
  (Q6).
- `docs/get-started/setup-databases.md`: one line that `setup` is the only step that pulls.
- `docs/specs/M11-14-container-install.md` item 3 table and the line "Nextflow pulls every tool
  image by digest (docker)": an "Amended by M11-19" note (the text stays, as for other amended
  specs).
- `missing_tool_reasons` docstring: drop "Nextflow still pulls … (N6b)"; say the docker rule.
- `changelog.d/<PR>.bugfix.md`: "Under `--profile docker`, Nextflow's tool containers no longer
  pull images at run time (`docker run --pull=never`), and preflight stops the run when a pinned
  image the selected modules need is not present locally, naming `gippa setup --images --profile
  docker`. Docker 20.10 or later is required."

## Examples (input → expected output)

`parse_docker_version` / `pull_never_problem`:

1. `Docker version 27.3.1, build ce12230` → `("docker", (27, 3))` → `""`.
2. `Docker version 20.10.24+dfsg1, build 297e128` → `("docker", (20, 10))` → `""` (the floor
   itself passes).
3. `Docker version 19.03.15, build 99e3ed8919` → `("docker", (19, 3))` → reason starting
   `docker 19.3 is older than 20.10` (`X.Y` is rendered from the parsed ints, so `19.03` reads `19.3`).
4. `Docker version 1.13.1, build 7d71120/1.13.1` (RHEL 7's docker) → `("docker", (1, 13))` →
   too old.
5. `podman version 4.9.3` → `("podman", (4, 9))` → `""`; `podman version 2.2.1` → reason starting
   `podman 2.2 is older than 3.0`.
6. `Docker version 24.0.7, build 24.0.7-0ubuntu4.1` → `("docker", (24, 0))` → `""`.
7. `nerdctl version 1.7.6`, `""`, `docker: command not found` → `None` → `""` (unknown never
   blocks).

`missing_tool_reasons` under `--profile docker` (fake `docker` first on `PATH`: `--version` prints
a configurable line, `image inspect <ref>` exits 0 iff `<ref>` is in a tags file, every call is
logged; default modules `core,mobile,bgc,function`, no precomputed layers):

8. Tags file holds all six distinct pinned refs, `Docker version 27.3.1` → no reason for bakta,
   transtermhp, genomad, antismash, miniprot, blast; the log has exactly one `--version` and six
   `image inspect` calls (blast's ref once), and no `pull`, `run` or `create`.
9. Tags file lacks the antiSMASH ref → `{"antismash": "image antismash/standalone@sha256:d3f3…
   is not present locally and runs never pull it (docker image inspect antismash/standalone@sha256:d3f3…
   failed); run: gippa setup --images --profile docker"}` (full digest); `run_preflight`'s
   `environment.tools` is `error`, `exit_code` 2, `value.missing` contains `antismash`.
10. Example 9 with `--modules core,function` → antiSMASH and geNOMAD are not required, so neither
    is inspected and neither is missing.
11. Example 9 with `annotation.precomputed.regions` set → antiSMASH and geNOMAD not inspected.
12. Example 8 with `Docker version 19.03.15` → every pinned tool has the too-old reason, and no
    `image inspect` is logged (rule 1 first).
13. No `docker` on `PATH` (empty bin dir only) → every pinned tool has reason 3 (contains
    `docker is not on PATH` and `gippa setup --images --profile docker`).
14. `image inspect` prints `Cannot connect to the Docker daemon at unix:///var/run/docker.sock.
    Is the docker daemon running?` and exits 1 → the first ref gets reason 4 quoting that line;
    the log has exactly one `image inspect`; every other pinned tool also has reason 4.
15. `--profile docker,conda`, tags file empty → the pinned tools are missing with reason 5 (the
    check applies under `docker,conda`); the conda-only tools follow their existing rows.
16. `--profile apptainer` with an empty cache → unchanged apptainer reasons; the docker log is
    empty. `--profile test,docker` → `{}`, the log is empty.
17. Fake `docker --version` prints `podman version 4.9.3`, the tags file holds all refs → `{}`.

Workflow (`nextflow config workflow -flat -profile <p>`, no daemon):

18. `-profile docker` → `docker.runOptions` = `--pull=never -u $(id -u):$(id -g)`.
19. `-profile docker` with `--db_dir /db` (or `-params-file` setting `db_dir`) →
    `--pull=never -v /db:/db:ro -u $(id -u):$(id -g)`.
20. `-profile docker,conda` → the same as 18; `-profile apptainer` → `apptainer.runOptions`
    contains no `--pull` and `docker.runOptions` is absent; `-profile test,stub` and
    `-profile conda` → no `docker.runOptions`.
21. `-stub-run` of `main.nf` under `-profile test,docker,stub` with the logging fake docker of
    `test_m11_14_docker_profile_stub.py` → every recorded `docker run` argv contains `--pull=never`
    exactly once, before the image ref (BAKTA, ANTISMASH and the other containerised processes).

## Edge cases

- **The image disappears after preflight** (`docker image prune` during the run): the task fails
  with docker's `No such image: <ref>` (exit 125) and `gippa run` reports the failed step as
  for any task failure; it is not pulled. A specific hint is deferred (Q8).
- **Same digest, other name.** `docker image inspect repo@sha256:<d>` succeeds only if the local
  image's `RepoDigests` holds that name and digest. An image pulled by tag whose digest equals the
  pin counts (Docker records the repo digest on pull); an image `docker load`ed from a tarball has
  no `RepoDigests` and does not count. `docker run --pull=never repo@sha256:<d>` answers the same
  way, so preflight and the task agree.
- **A wrong digest locally** (an older `oschwengers/bakta` pulled by tag) does not satisfy the pin:
  inspect by digest fails → reason 5; `setup --images` pulls the pin.
- **`docker.registry` set in the user's Nextflow config** (nf-core style): Nextflow prefixes the
  two Docker Hub short names, so the task runs a ref preflight did not inspect; already documented
  as unsupported in platforms.md › Image references (Q11).
- **A user `docker.runOptions`** (`$NXF_HOME/config`, launch-directory `nextflow.config`, a future
  `--nf-config`, M11-15) can replace GIPPA's string and drop `--pull=never`; GIPPA does not merge
  user Nextflow config. Documented, not enforced (Q6).
- **A repeated `--pull`** in `containerOptions` of a user config: the Docker CLI takes the last
  value; same remedy as above.
- **Remote daemons and contexts** (`DOCKER_HOST`, `docker context`, Docker Desktop): preflight
  runs `docker` with the environment `gippa run` passes to Nextflow, so it inspects the daemon the
  tasks will use. A Slurm executor with `docker` is refused earlier (M11-15 edge case (a)).
- **Rootless docker / `-u $(id -u):$(id -g)`**: unchanged (14p).
- **Old Docker without preflight** (`gippa run` with a probe that skips tools, or a hand-run
  Nextflow): `unknown flag: --pull` fails the first task, exit 125; nothing is pulled.
- **Empty selection.** No pinned tool required (every layer precomputed, `--modules` without the
  tools' modules) → no `docker --version`, no inspect, no reason.
- **Ties / duplicates.** BLAST and TBLASTN share a ref: one inspect; the reason is keyed by the
  registry tool `blast` (as today).
- **The step images of M11-15g** (ADR-0065, not merged) travel as params, not through
  `containers.config`; under docker they are docker-run by Nextflow and get `--pull=never` from
  item 1. Whichever of 15g and 19a merges second adds the step image refs to the docker check
  (Q10).
- **The test guard.** `tests/helpers/container_guard.py` refuses `docker run` without
  `--pull=never`; after item 1, a real-docker host running the stub suite through the guard no
  longer trips it on Nextflow's task argv.

## Acceptance criteria

1. Examples 1–7 hold for `parse_docker_version` and `pull_never_problem`; property: for any
   `major, minor ≥ 0`, `pull_never_problem(("docker", (major, minor)))` is `""` iff
   `(major, minor) >= (20, 10)`, and for podman iff `>= (3, 0)`; `pull_never_problem(None) == ""`.
2. Examples 8–17 hold for `missing_tool_reasons` and `missing_tools` (keys of the one equal the
   other). Every docker reason except rule 1 contains `<pkg> setup --images --profile docker`
   and the full ref.
3. Under every profile of the docker row of item 3's table, a pinned required tool whose ref is
   not local is missing; under the other rows the result equals today's (the f1 apptainer and
   conda rows pass unchanged).
4. No test or code path in this change runs `docker pull`, `docker run` or `docker create`: the
   fake docker's log holds only `--version` and `image inspect` lines, and `--version` at most
   once per call; `image inspect` at most once per distinct ref; after a no-daemon result, no
   further inspect.
5. #465's fgsea tests (`tests/unit/preflight/test_n6_fgsea_docker_image_preflight.py`,
   `tests/unit/test_m11_14c_runtime_and_preflight.py`) pass unchanged; the fgsea reason for an
   absent image is rule 5's text with the `images.json` ref.
6. `gippa preflight --profile docker` (CliRunner, fake docker, tags file without the Bakta ref)
   exits 2 and `preflight.md` shows `bakta: image oschwengers/bakta@sha256:… is not present
   locally …`; with all refs in the tags file it reports no missing pinned tool; `gippa run
   --profile docker` in the first case exits 2 and never starts the (fake) Nextflow.
7. `nextflow config -flat` gives examples 18–20 (integration, `needs_nextflow`); the
   docker-profile `runOptions` line in `workflow/nextflow.config` holds exactly one `--pull=never`,
   as its first token, and no other `--pull` (unit, text parse, no Nextflow).
8. Example 21: in the `-stub-run` under `test,docker,stub`, every `docker run` argv carries
   `--pull=never` once, before the image ref.
9. The apptainer profile block and `container_common.config` contain no `--pull`.
10. The docs and changelog of item 5 are present; `mkdocs build --strict -d <scratch>` passes;
    `scripts/ci/check_changelog.py` passes.
11. No network, no real docker/apptainer, no `$GIPPA_ROOT`/`$GIPPA_DB_DIR`; ruff, format and mypy
    clean.

No statistical claim, so no calibration test and no `stats` review.

## Existing tests that pin the old behaviour

Edited by the **lead** (test-author and implementer never edit existing tests), in the same PR,
before the implementation commit:

- `tests/unit/test_m11_14f1_preflight_images.py`:
  - `test_signature_takes_cache_env_root_and_local_tools_by_keyword` (bakta under `docker`
    `== []` with an empty tags file): put the Bakta ref in the fake's tags file;
  - `test_example_6_pinned_ostir_under_docker_needs_no_env` ("docker was queried for a pinned
    image"): becomes "present iff the ref is local; only `image inspect` is queried";
  - `test_row_pinned_image_docker_is_present_with_nothing_installed[docker|standard,docker]`:
    present iff every ref is in the tags file; otherwise missing with reason 5;
  - `test_row_pinned_image_and_env_the_env_is_ignored[docker]`: the image decides under docker
    too (inspect), still never the env;
  - `test_property_docker_depends_on_nothing_installed`: missing = conda-only tools ∪ pinned tools
    whose ref is not in the tags file (draw the tags set);
  - the module docstring (pinned under docker "present").
- Unchanged but re-run: `tests/unit/test_m11_14_profile_parity.py::test_docker_runs_as_the_calling_user`
  (one-line `runOptions`), `tests/integration/test_m11_14_docker_profile_stub.py`
  (`-u` present, no `-v`/`null` without `db_dir`), `tests/unit/test_m11_14_install_gaps.py::test_preflight_profile_matches_run_mode[docker]`,
  `tests/unit/test_m11_15a_cli_auto.py` (`--profile docker` runs), `tests/integration/test_m11_14u2_uniop_stub.py`
  (filters to UniOP). The implementer reports any other failure under `-m "not slow …"` to the
  lead instead of editing the test.

## Test plan

- `tests/unit/test_m11_19_docker_version.py`: examples 1–7, AC 1 (Hypothesis over
  `(major, minor)` and both families; a fake `docker --version` for `docker_client_version`,
  including no docker on `PATH` and a non-zero exit).
- `tests/unit/test_m11_19_pinned_images_docker.py`: examples 8–17, AC 2–5 (the f1 fake-tree
  fixtures: `workflow_main` patched, a logging fake docker with a tags file and a configurable
  `--version` line, put first on `PATH` after the autouse guard).
- `tests/unit/test_m11_19_preflight_cli.py`: AC 6 (CliRunner; fake Nextflow as in the 15a tests).
- `tests/unit/test_m11_19_runoptions_text.py`: AC 7 second half and AC 9 (parse
  `workflow/nextflow.config` and `container_common.config` text; no Nextflow).
- `tests/integration/test_m11_19_docker_pull_never.py` (`needs_nextflow`): examples 18–21, AC 7
  first half and AC 8 (reuse the `_config` and fake-docker `-stub-run` helpers of
  `test_m11_14_docker_profile_stub.py`; no daemon, no network). An nf-test is not needed: the
  flat config and the stub-run argv cover the profile, and `docker-smoke.yml` exercises real
  images on Actions (images pre-pulled).
- Targeted runs in the session; the full `-m "not slow and not calibration_full and not e2e"`
  suite on Slurm (`test_branch.sbatch`) before review.

## Size and PR split

| Slice | Content | Code | Tests | Prose |
|---|---|---|---|---|
| 19a (one PR) | item 1 (1 config line), item 2 (≈ 70 lines in `environment.py`: parse, version probe, reason classification, wiring into the pinned, fgsea and local-UniOP arms, docstring), item 5 docs, changelog fragment, the lead's edits of the pinned tests | ≈ 80 | ≈ 380 (≈ 300 new, ≈ 80 edited) | ≈ 40 |

Feature loop on branch `feat/m11-19-docker-pull-never`: test-author commits the new failing tests
→ lead edits the pinned existing tests → implementer makes them pass, docs and fragment →
review. Well under ≈ 400 code lines, so no split. A runtime hint for a task that fails with a
missing image (Q8) would be a separate later slice (19b), not planned.

## Questions for the lead

- **Q1** Report the Docker floor as a per-tool reason inside `environment.tools`, or as a new
  check `environment.docker`? **Proposed: per-tool reason** (no new rule id, no schema or writer
  change; the user sees it next to the tools it blocks).
- **Q2** An unparseable `docker --version` (nerdctl, a wrapper, an odd build string): block or
  pass? **Proposed: pass** (the run then fails at the first task with `unknown flag: --pull` if
  the flag is truly missing; blocking would lock out working CLIs such as nerdctl, which supports
  `--pull=never`).
- **Q3** Podman floor. `podman run --pull` predates 3.0, but the exact first version is not
  verified here and no podman is available on FarmShare. **Proposed: 3.0** (every distribution
  release since 2021 ships ≥ 3.0: Debian 11 3.0.1, Ubuntu 22.04 3.4.4, RHEL 8.4+), podman
  documented as best effort and untested.
- **Q4** A distinct "daemon not reachable" reason (rule 4) from Docker's stderr markers, or fold
  it into "not present locally"? **Proposed: distinct**; "not present, run setup" misleads a user
  whose daemon is stopped.
- **Q5** Keep `--pull=never` under `test,docker` too (the docker block has no test carve-out), so
  `docker-smoke` needs its images local? **Proposed: yes**; `docker-smoke.yml` already pre-pulls,
  and a pull inside a test is what the policy forbids.
- **Q6** A user's own `docker.runOptions` (Nextflow config outside GIPPA, or M11-15's future
  `--nf-config`) replaces GIPPA's string. Enforce (e.g. `gippa run` appends `--pull=never`
  through `containerOptions`), or document? **Proposed: document** in platforms.md; revisit when
  `--nf-config` lands.
- **Q7** Apptainer: Nextflow pulls a missing image into `apptainer.cacheDir` and offers no
  never-pull setting. Is preflight's cache rule enough? **Proposed: yes for 1.0**; add a backlog
  item (N6c) to look for a hard switch (e.g. a read-only cache dir or `NXF_OFFLINE`).
- **Q8** A task that fails because its image vanished after preflight: add a `gippa run` hint
  (task exit 125 and `No such image` / `image not known` in `.command.err` → "run `gippa setup
  --images --profile docker`")? **Proposed: defer** to the backlog (rare; the Docker error already
  names the image).
- **Q9** Record "only `gippa setup` pulls" as an ADR? **Proposed: no ADR**; it is stated in this
  spec, the install page and the `missing_tool_reasons` docstring. An ADR only if the lead wants
  it to bind later features (step images, singularity).
- **Q10** M11-15g step images (ADR-0065) under docker: whichever of 15g and 19a merges second
  adds the step refs to the docker presence check. **Proposed: yes**, noted on both PRs.
- **Q11** `docker.registry` set by the user: inspect the prefixed name instead? **Proposed: no**
  (already documented as unsupported; GIPPA does not read user Nextflow config).
- **Q12** `apptainer,docker` keeps resolving to apptainer (`container_runtime`), so no docker
  check there. **Proposed: unchanged.**
