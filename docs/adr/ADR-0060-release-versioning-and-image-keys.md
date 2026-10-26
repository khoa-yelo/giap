# ADR-0060: Release versioning and immutable, content-addressed tool images
- Status: proposed
- Date: 2026-10-25
- Spec sections: §3 (Containers row), §6.3, §9.2 Versioning row, §12 (`containers.yml`, `release.yml`, Releases), §14.6, §14.8 item 5, §3.10; R5, R6, R11; ADR-0008, ADR-0022, ADR-0059; M11-06, M11-13, M11-14
- Needs Khoa: yes (release policy, GHCR publishing and retention, rc channel; Q-K13–Q-K19 of M11-14)

## Context
Khoa (2026-10-25) wants a `pip install` path in which every tool runs from a container and
preflight passes. He also wants the final version to be easy to iterate: "Make sure that whatever is
built can be iterated to have 1.1, 1.2, 2.0, and so on."

The spec gives the principles: "semantic versioning; breaking column changes only in major
versions" (§9.2), and "containers pinned by digest" (§12). It does not say:

- what each bump covers;
- how image tags are named;
- whether images may change under a release;
- how quickly a patch can be cut.

M11-06 tags images `:<git-sha7>` and `:<tool version>`. The second tag is mutable: a rebuild with
the same tool version overwrites it. Every push to the paths rebuilds the tool, and
`containers.yml` has never run. The calibration evidence of M11-03 is reused only when nothing
under `src/` changed, so any patch reruns the full calibration.

## Decision
1. **Image identity.** A GHCR tool image is named `ghcr.io/<owner>/<pkg>-<tool>`, where `<pkg>`
   is the PyPI project name. Its only tag is `<version>-<key>`:
   - `<key>` is the first 12 hex of sha256 over the `path\tmode\tsha256` lines of the tool's
     build inputs, sorted by the code points of the path. The inputs are its
     `workflow/containers/<tool>/` files, `workflow/envs/<tool>.yml` and its lock files; fgsea also
     includes `workflow/bin/giap_fgsea.R`, which is baked into the image.
   - Dockerfiles install from the locks, and every fetched source archive is checksummed
     (`ADD --checksum`).
   - A tag is pushed once and never overwritten. A tool whose key already exists in the registry
     is not rebuilt, and its digest is reused.
   - Tag existence is checked **authenticated**, and only a registry "not found" counts as absent;
     any other error fails the run. The check is repeated immediately before the push, and the
     workflow runs in one global concurrency group. The anonymous check is used only to verify
     that a package is public.
   - `containers.config` pins bare digest references (ADR-0059). `workflow/conf/images.json`
     records key, tag, digest, processes, source commit and licence.
     An image that no Nextflow process uses (fgsea) is pinned in `images.json` only.
   - Both files ship in the wheel, so a release pins its exact images forever.
2. **No deletion, no retagging.** Agents never delete, retag or change the visibility of a package
   or a version. Every digest named by a released `containers.config` or `images.json` is kept forever. Pruning
   unreleased digests is Khoa's choice alone. The default is never; if Khoa prunes, the floor is 180
   days, keeping the latest digest per tool. UniOP (unlicensed, R5/R11) is never published. Each
   user's `<pkg> setup --uniop` (included in `setup --images`) fetches the pinned commit from the
   authors' GitHub with a sha256 check and builds it locally, after the user accepts the terms
   (M11-14 item 1b, 14u1/14u2). If the authors grant a licence, it is published like the others and
   the local build becomes a fallback.
3. **What bumps what.**
   - **Major:** anything that breaks a user's script or a cross-run comparison. That covers
     removed, renamed, retyped or reinterpreted columns, tables, manifest keys, flag tokens, config
     keys or CLI options; an output-schema major; and a database pin whose old database no longer
     works.
   - **Minor:** a tool version or image key change, a database pin bump, a new column, table,
     flag, option or module, a changed default, any change to a default run's output bytes on the
     golden and TE7 runs, or a calibration rerun that moves a gate number.
   - **Patch:** everything else, which means bug, crash, docs and performance fixes whose goldens
     and TE7 stay byte-identical.

   The output `schema_version` is `MAJOR.MINOR`. Its major equals the package major, and a new
   column bumps its minor. A static test diffs the table schemas against the previous release tag.
4. **Release candidates and gates.**
   - `X.Y.ZrcN` is cut on main by the release-prep PR. Rcs go to TestPyPI only, and finals go
     to PyPI.
   - Gate reruns scale with the bump (the table in M11-14 item 4). Calibration, pilot and
     benchmark evidence from commit C is reused at R only when `reusable(C, R)` holds, for every
     bump below major. The predicate holds when every path of
     `git diff --name-only --no-renames C R` (deletions and both sides of a move included) lies
     under a `NON_STATS_PATHS` entry, matched on whole path components.
     - The only exception is a diff to `src/giap/__init__.py` that is a single `__version__`
       hunk.
     - A diff touching `scripts/release/readiness.py`, which defines the list, is never reusable.
   - `NON_STATS_PATHS` holds only paths that cannot change a statistic:
     - docs and paperwork;
     - `.github/`, `.claude/`, `scripts/dev/`, `scripts/docs/`;
     - `src/giap/report/` and `src/giap/viz/`.
   - A disjoint never list is a constant too. It names `preflight/`, `runner/`, `setupdb/`,
     `cli.py`, `data/`, `config/`, `rng.py`, `io/`, `workflow/`, `pyproject.toml`, `uv.lock`, the
     simulator, the calibration tests and `ci/`. Changing either list needs `stats` review.
   - Khoa approves every release, patch included (§14.8 item 5).
5. **Flow and publishing rights.**
   - A merged change to a tool's inputs triggers `containers.yml`. The tool is built only if its
     key is new, and only when the repository variable `GHCR_PUBLISH` is `true` and Khoa approves
     the `ghcr` environment. Khoa alone creates both of these settings.
   - PRs that change `.github/workflows` (M11-14 14b, 14g, 14h, 14p) are CI changes and merge
     only with Khoa's approval.
   - The build is followed by a digest PR, then the FarmShare smoke, then the merge.
   - `release.yml` verify recomputes the keys. It also checks that every GHCR pin resolves
     anonymously.
   - Agents never push tags, never dispatch with `dry_run=false`, never set `GHCR_PUBLISH` and
     never approve `ghcr`. The first publish is Khoa's, after the GIAP → GIPPA rename, so no image
     is ever published under the old name.
6. **Compatibility.** The run manifest already records the package version, the tool versions
   and digests, and the database versions; M11-14 adds `operon_source`. `compare` warns across
   minor versions, tool versions or digests (printed as the version plus 12 hex of the digest),
   database versions and operon sources. It refuses to compare across output-schema majors unless
   `--allow-cross-major` is given. `output_schema_version` starts at `1.0` with 1.0.0. Before then
   no rule applies.

## Evidence (before → after: metrics, tests, simulations)
None yet. This is a policy ADR, and its tests are M11-14 AC 1–11, all offline.

Before:

- 0 GHCR images;
- 6 conda-only tools plus fgsea (7) under the container profiles;
- a mutable `:<tool version>` tag;
- any `src/` change reruns `calibration_full`.

After 14j:

- 0 conda-only tools on the pip path (UniOP is built locally by `setup`, not published);
- one immutable tag per build;
- a docs-or-CLI-only patch reuses the calibration evidence.

## Consequences
- 1.1, 1.2 and so on rebuild only the tools whose inputs changed. Older wheels keep resolving
  their own digests, so rollback is `pip install <pkg>==1.x`.
- The M11-06 `:<sha7>` and `:<tool version>` tags are dropped before anything is published, so no
  published tag is ever withdrawn.
- A stale pin (inputs changed, no new digest) blocks a release through readiness H4 and
  `image_key.py --check`.
- Releases reuse evidence only for diffs inside `NON_STATS_PATHS` (plus the version line). Any
  other diff pays the calibration, pilot and benchmark reruns, by design (R3).
- The expected pip path runs UniOP (Khoa, 2026-10-25): pip install → setup (incl. the local
  UniOP build) → preflight passes with UniOP. The distance-model fallback (`enabled: auto`)
  covers only failed builds (offline, no build rights). It becomes the default only after the
  statistician's attack and the ADR-0022 operating-point addendum (M11-14 14s1, 14s2, 14d2).
  Until then, a user whose build failed sets `enabled: auto` or `false` explicitly.
- On acceptance: M11-06 is amended (tags, publishing guard), and ADR-0008's image section gains
  "tags superseded by ADR-0060".
