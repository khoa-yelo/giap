# M11-06: `containers.yml` and GHCR images by digest for the conda-only tools

Milestone M11 · **pass 2** · depends on: M11-01 (soft) · estimate ≈ 210 code lines + ≈ 160 test
lines · not `stats` · roles: implementer (`.github/workflows/containers.yml`,
`scripts/release/update_digests.py`, `workflow/conf/containers.config`, `nextflow.config` profile
lines), test-author (`tests/unit/test_m11_06_*.py`, the Apptainer smoke under
`tests/integration/`) · Khoa: Q-K8 (make the packages public).

## Spec sections implemented

- §13 row M11 "containers by digest"; §12 `containers.yml`: "changes to `workflow/containers/**` or
  `workflow/envs/**`, manual | build per-tool images, push to `ghcr.io/<owner>/giap-<tool>`, open a
  PR updating the digests in `workflow/conf/containers.config`"; §12 Releases "containers pinned by
  digest"; §14.6 "`containers.config` by digest"; §17 "prefer `ghcr.io` images built by
  `containers.yml`".
- §6.3 container field `ghcr.io/<owner>/giap-<tool>@sha256:…`; G19 (Apptainer runs OCI images by
  digest).
- R5/R11: UniOP is unlicensed — its Dockerfile says "never pushed publicly"; it is **not** built by
  this workflow. ADR-0008 (pins), M8-17 decision 9 (fgsea conda-only).

## Definition

1. **Images built**: `promotech`, `promoter_calculator`, `terminatornet`, `bactermfinder`, `ostir`
   (5). Not `uniop` (licence) and not `fgsea` (deferred by M8-17 decision 9; no Dockerfile).
2. **`containers.yml`**: triggers `push` to main with paths `workflow/containers/**`,
   `workflow/envs/**`, and `workflow_dispatch` (input `tools`, default all five). A matrix job per
   tool: `docker/build-push-action` (SHA-pinned, Q-L5), context = repository root, `file =
   workflow/containers/<tool>/Dockerfile`, platform `linux/amd64`, tags
   `ghcr.io/khoa-yelo/giap-<tool>:<git-sha7>` and `:<tool version>` (from the env file pin),
   labels `org.opencontainers.image.source`, `.revision`, `.licenses` (the NOTICE licence). Permissions
   `packages: write` on that job only. Output: the pushed manifest digest.
   A final job collects the digests and runs `update_digests.py` and opens a PR
   (`ci: pin GHCR digests (<sha7>)`, branch `ci/digests-<sha7>`) with `peter-evans/create-pull-request`
   (SHA-pinned) or `gh pr create`; never pushes to main. The job passes
   `--nextflow-config workflow/nextflow.config`, so the newly pinned processes join the apptainer
   profile's `withName: '…' { conda = '' }` selector in the same PR (lead decision on test-author
   ambiguity 5: the pinned set equals the conda-off set).
3. **`scripts/release/update_digests.py --digests digests.json [--config workflow/conf/containers.config] [--nextflow-config workflow/nextflow.config] --date YYYY-MM-DD --source-commit SHA`**:
   for each `{process: "ghcr.io/khoa-yelo/giap-<tool>@sha256:<hex>"}` adds or replaces exactly the
   `withName: '<PROCESS>' { container = 'docker://…' }` line, keeps every other line and comment
   byte-identical, updates the "Conda-only for now" comment list, and appends the build date +
   source commit to the "Tags the digests were resolved from" comment. Process names: `PROMOTECH`,
   `PROMOTER_CALCULATOR`, `TERMINATORNET`, `BACTERMFINDER`, `OSTIR` (+ any second-round process that
   reuses an image, e.g. the novel-ORF Promotech/OSTIR runs: the implementer lists them from
   `workflow/modules/local/*/main.nf`). Idempotent.
4. **Profile**: under `-profile apptainer` the five processes get `conda = ''` like the upstream-image
   tools (the `withName` regex in `nextflow.config`), once pinned. UniOP and fgsea keep conda.
5. **Smoke on FarmShare** (Slurm, after the digest PR's branch exists): `apptainer exec
   docker://…@sha256:… <tool> --help` (or the module's version command) for each image, then the
   existing contract tests of the five tools under `-profile farmshare,apptainer` on the 100-kb
   genome via `scripts/slurm/nf_process_real.sbatch`. Results and job IDs go into the digest PR.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `update_digests.py` on today's `containers.config` with `{"OSTIR": "ghcr.io/khoa-yelo/giap-ostir@sha256:<64 hex>"}` | one new line `withName: 'OSTIR' { container = 'docker://ghcr.io/khoa-yelo/giap-ostir@sha256:…' }` inside `process { }`; `ostir` removed from the conda-only comment; all 7 existing lines unchanged |
| 2 | the same input run twice | second run: no diff |
| 3 | an existing OSTIR line with an older digest | replaced in place (same line number) |
| 4 | digest `sha256:abc` (not 64 hex) | exit 1, `bad digest for OSTIR` |
| 5 | `{"UNIOP": …}` | exit 1, `UNIOP must not be pinned to a public image (unlicensed; R5)` |
| 6 | a PR changing `workflow/envs/ostir.yml` merged to main | `containers.yml` builds `ostir` only (paths filter maps env file → tool); a digest PR follows |
| 7 | dispatch with `tools=uniop` | the matrix rejects it before building: `uniop is never built here` |

## Edge cases

- **Two processes, one image** (Promotech first and second round): both lines get the same digest.
- **Base image digest** (`mambaorg/micromamba@sha256:…`) is already pinned; the build is
  reproducible only up to conda solves, which the env files pin exactly; the image label records
  the env file's sha256.
- **Package still private** (Q-K8 not done): Apptainer pull fails with 401/403 on FarmShare; the
  smoke reports `private package` and the digest PR stays draft; PROGRESS › Questions for Khoa.
- **Promotech's frozen Python 3.6 env** fails to build (§17): fallback `python:3.6-slim` image per
  §17; the tool stays conda-only until it builds; the PR says so.
- **GHCR rate limits / outage**: R6 — retry once later, then Blocked.
- Genomic edge cases: not applicable.

## Acceptance criteria

1. `update_digests.py` tests cover examples 1–5 and the two-process edge case on copies of the real
   config (no network).
2. A static test of `containers.yml`: triggers and paths as in item 2; `uniop` absent from the
   matrix; `packages: write` only on the build job; every non-`actions/*` action SHA-pinned; no push
   to `main`.
3. After merge and one dispatch, the digest PR exists; the FarmShare smoke passed for all five
   (job IDs in the PR body); `readiness.py` H4 passes with UniOP and fgsea as the only conda-only
   tools.
4. NOTICE.md gains one sentence: the five GHCR images contain GPL-3.0 tools (Promotech, Promoter
   Calculator, OSTIR, BacTermFinder) fetched from upstream, with source links (GPL source
   availability).

## Estimate

Code ≈ 210 (workflow ≈ 110, `update_digests.py` ≈ 90, config lines ≈ 10). Tests ≈ 160.

## Questions for the lead (reversible defaults in use)

- Q-K8 (Khoa) public GHCR packages. Default: yes, five packages.
- Q-L8 no `giap` core image. Q-L5 SHA pins.
- `linux/arm64` images? Default: no (FarmShare and the tested path are amd64; M11-09 documents it).
