# M11-09: platforms page and a Docker-profile smoke on GitHub Actions

Milestone M11 · **pass 2** · depends on: M11-07 (site), M11-06 (soft: more images to smoke) ·
estimate ≈ 90 code lines + ≈ 60 test lines + ≈ 80 prose lines · not `stats` · roles: implementer
(`docker-smoke` job in `.github/workflows/ci.yml`, `workflow/tests/docker_smoke.nf.test`), docs
(`docs/user-guide/platforms.md`), test-author (`tests/unit/test_m11_09_*.py`) · Khoa: Q-K10 (Actions minutes).

## Spec sections implemented

- Backlog (Khoa, 2026-10-24): "running GIAP on macOS is untested. The Python package passes on
  macOS CI; the full pipeline would need `-profile docker`, which spec §3 keeps but doesn't test.
  Candidate for M11: a Docker-profile test run, and documenting the platforms."
- §3 container row: "The `docker` profile is kept for users elsewhere but is not tested here"
  ("here" = FarmShare; GitHub's ubuntu runners have Docker). G19. §12 Actions minutes (≈ 15 min per
  run; heavy suites stay on FarmShare).

## Definition

1. **`docs/user-guide/platforms.md`** (in the nav under User guide): a support matrix, each cell one of
   `tested` (with where: CI job or FarmShare suite), `smoke-tested`, `untested`:

   | Platform | Python package (`giap preflight`, `convert`, `report`, `compare`) | Full pipeline |
   |---|---|---|
   | Linux x86-64 + Apptainer (FarmShare) | tested (CI ubuntu; nightly) | tested (`-profile farmshare,apptainer`; nightly, weekly) |
   | Linux x86-64 + conda | tested (CI) | tested for the conda-only tools (nightly) |
   | Linux x86-64 + Docker | tested (CI) | smoke-tested (this issue: two tools) |
   | macOS (arm64 / x86-64) | tested (CI macOS 3.12) | untested; Docker Desktop runs amd64 images under emulation; Bakta full DB size and antiSMASH ~9 GB image are practical limits |
   | Windows | untested; WSL2 as Linux | untested |

   Plus: the conda-only tools (UniOP, fgsea, and any not yet on GHCR) need `conda` alongside
   `docker`; DIA-NN is the user's own.
2. **`docker-smoke` job** (ubuntu-latest, Java 21, Nextflow as in `nf-test-stub`): runs
   `nf-test test workflow/tests/docker_smoke.nf.test --profile test,docker`, which runs the real
   `MINIPROT` and `BLAST`/`TBLASTN` processes (images already pinned, small: < 200 MB each) on the
   100-kb genome's fixture inputs under `tests/data/` and asserts the outputs parse with the
   existing parsers (non-empty, expected header). Runs on PRs that touch `workflow/**` or
   `workflow/conf/containers.config`, and on `schedule` weekly on main (Q-K10). Time budget: job
   `timeout-minutes: 10`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | PR touching `workflow/conf/containers.config` | `docker-smoke` runs; pulls by digest; both processes succeed; ≤ 10 min |
| 2 | PR touching only `src/giap/stats/` | `docker-smoke` skipped (paths filter) |
| 3 | a digest in `containers.config` replaced by a non-existent one (test fixture branch) | the job fails at pull with the image name in the log |
| 4 | `docs/user-guide/platforms.md` | every `tested` cell names a CI job or a Slurm script that exists (test) |

## Edge cases

- **Docker Hub / quay.io rate limits** on shared runners: retry the pull once (`docker pull` with
  one retry), then fail with `rate-limited` in the summary — not a code failure, rerun.
- **Images not multi-arch**: irrelevant on ubuntu x86-64; documented for macOS arm64.
- **User inside the container** (images expecting root, `$HOME`): `docker.runOptions = '-u
  $(id -u):$(id -g)'` set in the `docker` profile if the smoke shows permission errors; recorded in
  the module's `meta.yml` (§17).
- Genomic edge cases: the fixture genome is the existing circular-tagged 100-kb genome; no new
  coordinates logic.

## Acceptance criteria

1. `docs/user-guide/platforms.md` exists with the matrix; the test of example 4 passes.
2. The `docker-smoke` job is green on its PR; its run time is recorded in the PR body.
3. A static test asserts the job's paths filter, timeout and that it is not in the `summary` job's
   required set (it is advisory until Khoa says otherwise).

## Estimate

Code ≈ 90 (job ≈ 35, nf-test ≈ 55). Tests ≈ 60. Prose ≈ 80.

## Questions for the lead (reversible defaults in use)

- Q-K10 (Khoa) Actions minutes for the smoke. Default: path-filtered PRs + weekly.
- Required check or advisory? Default: advisory (not in `summary.needs`), so registry outages do not
  block merges.
