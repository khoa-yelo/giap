# ADR-0060: Release versioning and immutable, content-addressed tool images
- Status: proposed
- Date: 2026-10-25
- Spec sections: §3 (Containers row), §6.3, §9.2 Versioning row, §12 (`containers.yml`, `release.yml`, Releases), §14.6, §14.8 item 5, §3.10; R5, R6, R11; ADR-0008, ADR-0059; M11-06, M11-13, M11-14
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
   - `<key>` is the first 12 hex of sha256 over the sorted `path\tsha256` lines of the tool's
     build inputs. The inputs are its `workflow/containers/<tool>/` files, `workflow/envs/<tool>.yml`
     and its lock files.
   - Dockerfiles install from the locks.
   - A tag is pushed once and never overwritten. A tool whose key already exists in the registry
     is not rebuilt, and its digest is reused.
   - `containers.config` pins bare digest references (ADR-0059). `workflow/conf/images.json`
     records key, tag, digest, processes, source commit and licence.
   - Both files ship in the wheel, so a release pins its exact images forever.
2. **No deletion, no retagging.** Agents never delete, retag or change the visibility of a package
   or a version. Every digest named by a released `containers.config` is kept forever. Pruning
   unreleased digests is Khoa's choice alone. The default is never; if Khoa prunes, the floor is 180
   days, keeping the latest digest per tool. UniOP (unlicensed, R5/R11) is never published. It is
   built only locally, on the user's request.
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
   - Gate reruns scale with the bump (the table in M11-14 item 4). A patch reuses calibration,
     pilot and benchmark evidence from commit C when `git diff C R` lies inside
     `NON_STATS_PATHS`.
   - Khoa approves every release, patch included (§14.8 item 5).
5. **Flow and publishing rights.**
   - A merged change to a tool's inputs triggers `containers.yml`. The tool is built only if its
     key is new, and only when the repository variable `GHCR_PUBLISH` is `true` and Khoa approves
     the `ghcr` environment.
   - The build is followed by a digest PR, then the FarmShare smoke, then the merge.
   - `release.yml` verify recomputes the keys. It also checks that every GHCR pin resolves
     anonymously.
   - Agents never push tags, never dispatch with `dry_run=false`, never set `GHCR_PUBLISH` and
     never approve `ghcr`. The first publish is Khoa's, after the GIAP → GIPPA rename, so no image
     is ever published under the old name.
6. **Compatibility.** The run manifest already records the package version, the tool digests and
   the database versions. `compare` warns across minor versions, tool digests and database
   versions. It refuses to compare across output-schema majors unless `--allow-cross-major` is
   given.

## Evidence (before → after: metrics, tests, simulations)
None yet. This is a policy ADR, and its tests are M11-14 AC 1–11, all offline.

Before:

- 0 GHCR images;
- 6 conda-only tools under the container profiles;
- a mutable `:<tool version>` tag;
- any `src/` change reruns `calibration_full`.

After 14j:

- 0 conda-only tools on the pip path (UniOP is optional);
- one immutable tag per build;
- a docs-or-CLI-only patch reuses the calibration evidence.

## Consequences
- 1.1, 1.2 and so on rebuild only the tools whose inputs changed. Older wheels keep resolving
  their own digests, so rollback is `pip install <pkg>==1.x`.
- The M11-06 `:<sha7>` and `:<tool version>` tags are dropped before anything is published, so no
  published tag is ever withdrawn.
- A stale pin (inputs changed, no new digest) blocks a release through readiness H4 and
  `image_key.py --check`.
- Patch releases become cheap only for changes outside `NON_STATS_PATHS`. Anything else pays the
  calibration rerun, by design (R3).
- If Khoa answers Q-K13 by refusing the distance-model default for UniOP, the pip path loses
  "preflight passes" for core runs until a licence arrives.
- On acceptance: M11-06 is amended (tags, publishing guard), and ADR-0008's image section gains
  "tags superseded by ADR-0060".
