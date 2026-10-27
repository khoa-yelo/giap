# ADR-0065: GIPPA step images `<pkg>-step` and `<pkg>-step-r`, tagged `<version>`
- Status: proposed (Q-K39; lead default: accept)
- Date: 2026-10-27
- Spec sections: §0.2 R1, R5, R6; §3 (Containers row); §12; §13; M11-06 (Q-L8); M11-14 14c;
  M11-15 item 8 (8a, slice 15g; Q-K37, Q-K39); ADR-0038, ADR-0059, ADR-0060
- Needs Khoa: yes (reverses an earlier decision, amends ADR-0060 decision 1, publishes new
  images; Q-K39)

## Context
`GIPPA_STEP` tasks run `gippa` itself. Today they run on the host, so every compute node needs a
Python environment with GIPPA installed (on FarmShare the shared venv). M11-15 makes the pipeline
portable to clusters where the user has only Nextflow and a container runtime, so the steps need
an image too. M11-06 decided "no `gippa` core image" (Q-L8, reading R1 as "per-tool images
only"); M11-15 item 8 reverses that.

ADR-0060 decision 1 names every GHCR image `ghcr.io/<owner>/<pkg>-<tool>` with the single tag
`<version>-<key>`, where `<key>` hashes the image's build inputs, and pins images by digest from
inside the release (`containers.config`, `images.json`). A step image cannot follow that rule:
it contains the release wheel, so the wheel cannot carry the image's digest (or a key over
inputs that include the wheel itself).

## Decision
1. **A package image.** Two images run GIPPA's own steps, so compute nodes need no Python
   environment: `workflow/containers/step/Dockerfile` (`<pkg>-step`) for every `GIPPA_STEP`
   task, and `workflow/containers/step-r/Dockerfile` (`<pkg>-step-r`) for `FUNCTION_STEP` and
   `SENSITIVITY`, which start R fgsea (a container cannot start a sibling container). This
   reverses M11-06 Q-L8.
2. **Names and tags (amends ADR-0060 decision 1).** The names follow ADR-0060:
   `ghcr.io/<owner>/<pkg>-step` and `ghcr.io/<owner>/<pkg>-step-r`. The tag is the package
   version alone, `<version>` = `gippa.__version__`, not `<version>-<key>`. The image is a
   function of the release wheel and the committed lock, and a PyPI version is immutable, so the
   version is the content key. The wheel names the tag; it cannot pin the digest of an image that
   contains it. `gippa setup` pulls by tag, checks the `org.opencontainers.image.version` label,
   verifies the digest against the release asset `step-images-<version>.json` (else records it
   on first use with a warning) and records it in `db_manifest.json`; `gippa run` passes
   `<ref>:<version>@<digest>`, so every run is still pinned by digest (§12). ADR-0060's no-retag
   rule applies unchanged. Dev builds (`scripts/dev/build_step_image.sh`) use the same names
   with the tag `dev-<sha12 of the wheel>`. The step images are not keyed tool images:
   `scripts/release/image_key.py` and the readiness H4 check skip them.
3. **Contents (GIPPA stays MIT, Q-K34).** `<pkg>-step` holds GIPPA and its permissive Python
   dependencies only, installed with `pip install --require-hashes -r requirements.lock`
   (`uv export --frozen --no-dev --no-emit-project`, checked against `uv.lock` by a test), then
   the wheel with `--no-deps`. `scripts/release/check_image_licences.py` runs inside the built
   image and fails on any Python `dist-info` whose licence is not in
   `workflow/containers/step/licences.allow`. MPL-2.0 (certifi) is allowed: file-level
   copyleft, shipped unmodified. GPL, LGPL, AGPL, non-commercial and unlicensed Python packages
   are not. No bioinformatics tool and no conda in `<pkg>-step`.
4. **OS packages are aggregated, not linked.** The Debian `python:3.12-slim` base (pinned by
   digest), `procps` (GPL-2; Nextflow's task metrics need `ps`) and `ca-certificates` are
   separate system programs. In `<pkg>-step-r`, R (GPL-2) and fgsea (MIT) come unmodified from
   the fgsea image's explicit lock and run only as a separate process started by `subprocess`,
   never imported (R5), as in the published `gippa-fgsea` image (Q-K37). NOTICE.md documents
   both images and where the corresponding sources are.
5. **Marker file.** Both images carry `/opt/gippa/STEP_IMAGE`, whose content is the image's own
   fgsea prefix (`/opt/conda` in step-r, empty in step), so a host `GIPPA_FGSEA_PREFIX` exported
   by Nextflow's `env` scope cannot leak into the container (M11-15 item 8a; wired in 15h).

## Evidence (before → after: metrics, tests, simulations)
- `tests/unit/test_m11_15g_step_dockerfiles.py`: base by digest, procps, the hashed lock before
  the wheel, the marker, labels, no conda/R/tool in step, exactly the fgsea lock in step-r, not
  keyed by `image_key.py` or readiness H4.
- `tests/unit/test_m11_15g_requirements_lock.py`: the lock equals the runtime closure of
  `uv.lock`, hashes included, and re-exports identically.
- `tests/unit/test_m11_15g_image_licences.py`: one GPL dist fails and is named (AC 7); the
  locked dependency set passes the committed allowlist.

## Consequences
- Every release builds two more images from the release wheel before the PyPI upload and
  publishes their digests as a release asset (15n). A release whose step images fail to build
  or fail the licence check does not ship.
- A dependency bump re-exports `requirements.lock`; a new licence needs an allowlist entry with a
  reason, reviewed like any other change. An LGPL or GPL runtime dependency blocks the image.
- The wheel names a tag rather than a digest; the digest pin moves to setup and the manifest.
