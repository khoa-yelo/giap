# ADR-0059: Bare image references in `containers.config` (no `docker://` prefix)
- Status: proposed
- Date: 2026-10-25
- Spec sections: §3 (Containers row), §6.3 (container field), §6.x Profiles (`apptainer` bullet), §12, G19; supersedes in part ADR-0008 ("Image digests" wording)
- Needs Khoa: no (default of M11-09b Q5: the change makes the config match §6.3's own form and repairs the `docker` profile; Apptainer behaviour and cache names are unchanged)

## Context
GIAP_SPEC.md §3 says images are "referenced by digest and run as `docker://…@sha256:…` through
Apptainer". The `apptainer` profile bullet says "images `docker://…@sha256:…` cached in
`$NXF_APPTAINER_CACHEDIR`". ADR-0008 records the pins as "`docker://…@sha256:…`" in
`workflow/conf/containers.config`. That is how M2-01 wrote them.

The same spec keeps a `docker` profile "for users elsewhere", and the README quickstart and
M11-09's `docker-smoke` job use it. Under that profile Nextflow passes the container string to
`docker run` unchanged, and the Docker CLI rejects the `docker://` transport. Every pinned process
fails. The M11-09 implementer's logging fake `docker` received
`docker run … docker://quay.io/biocontainers/miniprot@sha256:…`.

§6.3 already writes the container field bare (`ghcr.io/<owner>/giap-<tool>@sha256:…`). The bare
form is also the nf-core convention.

## Decision
1. Each `container =` value in `workflow/conf/containers.config` is a bare digest reference. It
   fully matches `^[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}$`, with no scheme and no tag. The
   seven digests and the tags they were resolved from stay as ADR-0008 records them.
2. Code that writes a reference for that file emits the bare form: `scripts/release/update_digests.py`
   and `scripts/resolve_digest.py`. Code that reads one accepts an optional `docker://`: `giap setup`'s
   `ToolRunner` and the Slurm scripts. `build_envs.sbatch` adds `docker://` itself when it calls
   `apptainer pull`.
3. Readiness check H4 rejects a `docker://` prefix (`docker:// prefix: <process>`), so the
   breakage cannot return.
4. Reading of the spec: §3 and the `apptainer` bullet describe how Apptainer *runs* the images,
   and that is still `docker://…@sha256:…`, because Nextflow adds the prefix. They do not prescribe
   the text of `containers.config`. Where the two disagree, the config holds the bare form, as in
   §6.3. ADR-0008's "`docker://…@sha256:…`" is superseded on this point only.
5. Docker Hub images keep their short names (`oschwengers/bakta`, `antismash/standalone`), as
   `resolve_digest.py` prints them. A user config that sets `docker.registry` or
   `apptainer.registry` would rewrite these two names. This is documented, not guarded (M11-09b Q2).

## Evidence (before → after: metrics, tests, simulations)
Nextflow 26.10.0 is the version pinned by `scripts/bootstrap_farmshare.sh` and ADR-0002, and
`nextflow -version` in `giap-wf` prints `26.10.0 build 14055`. The sources below were read at tag
`v26.10.0` of https://github.com/nextflow-io/nextflow (official, R6) on 2026-10-25:

- `docs/container/apptainer.mdx`, line 71: "You do not need to specify `docker://` to pull from a
  Docker registry. When Apptainer is enabled, Nextflow prepends it to your image name. The Docker
  engine does not work with containers specified as `docker://`."
- `modules/nextflow/src/main/groovy/nextflow/container/ContainerHandler.groovy`:
  - lines 99–100 and 151–159: under Docker with no `docker.registry`, the name is passed through
    unchanged;
  - lines 229–256 (`normalizeApptainerImageName`): a name with no scheme, no path and no matching
    local file becomes `"docker://" + name`;
  - line 182: `IMAGE_URL_PREFIX`;
  - lines 81–84 and 108–111: the result goes to `ApptainerCache.getCachePathFor`.
- `DockerBuilder.groovy`, lines 185–186: the image is appended verbatim to `docker run`.
- `SingularityCache.groovy`, which `ApptainerCache` extends:
  - lines 88–102 (`simpleName`): drop everything up to `://`, map `:` and `/` to `-`, append `.img`;
  - lines 201–203: the cache file is `<cacheDir>/<simpleName>`.

  So `docker://X` and bare `X` name the same cache file.
- The same functions are unchanged at `v25.04.8`, the §3.4 floor.

Before the change:

- `-profile docker` fails on every pinned process.
- `-profile apptainer` uses `$GIAP_ROOT/containers/<simpleName>.img`. Six such files are
  present on FarmShare on 2026-10-25, for example `quay.io-biocontainers-miniprot@sha256-2eb5…f.img`.

After the change:

- `docker run … quay.io/biocontainers/miniprot@sha256:…` works.
- Apptainer resolves the same six file names, so nothing is pulled again.

The tests are M11-09b AC 1–10. They include an `-stub-run` integration test with fake `docker`
and `apptainer` binaries and a seeded cache.

## Consequences
- The `docker` profile works for the pinned tools. M11-09's `docker-smoke` job and the README
  quickstart are unblocked.
- No re-pull and no change under Apptainer on FarmShare. No digest or tool version changes.
- Ten existing tests assert the `docker://` form and are edited by the lead in the M11-09b PR.
  The historical records keep the old form, because it is the URI that Apptainer pulled:
  `tests/data/**/PROVENANCE.md` and the M2-01, M2-06, M8-17, M11-03 and M11-06 specs.
- Future pins use the bare form, including M8-17's deferred `FGSEA_IMAGE`. H4 enforces this.
- On acceptance, ADR-0008's Status line gains "; reference form superseded by ADR-0059".
