# M11-04: `release.yml` (verify, build, and — only when Khoa dispatches it — publish and tag)

Milestone M11 · **pass 1 (MVP)** · depends on: M11-01 (`check_dist.py`), M11-03 (`readiness.py
verify`) · estimate ≈ 190 code lines + ≈ 170 test lines · not `stats` · roles: implementer
(`.github/workflows/release.yml`, `scripts/release/notes.py`), test-author
(`tests/unit/test_m11_04_*.py`) · Khoa: Q-K1–Q-K4 (account set-up; dispatch).

## Spec sections implemented

- §12 `release.yml`: "manual (version input) | verify gates, changelog and the version; build sdist
  and wheel; publish to PyPI (trusted publishing); create the tag and GitHub release (the Zenodo DOI
  comes from the GitHub integration). Bioconda recipe PR as a manual follow-up".
- §3.10: "Never push tags"; "Tags and releases come from `release.yml` (workflow_dispatch) or from
  Khoa". §14.8 item 5 "Each release". R12 (no tokens: trusted publishing needs none).

## Definition

Trigger: `workflow_dispatch` only, inputs:

| input | type | default | meaning |
|---|---|---|---|
| `version` | string, required | — | e.g. `1.0.0`; must match `^\d+\.\d+\.\d+(rc\d+)?$` |
| `dry_run` | boolean | **true** | verify + build + smoke + notes; no upload, no tag |
| `target` | choice `testpypi` / `pypi` | `testpypi` | index to upload to when `dry_run` is false |

Top-level `permissions: contents: read`. Jobs:

1. **`verify`** (any actor): checkout with `fetch-depth: 0` at `github.sha`; refuse unless
   `github.ref == 'refs/heads/main'`; `giap.__version__ == inputs.version`; tag `v<version>` absent
   (`git ls-remote --tags origin`); `readiness.py verify validation/release/<version>/readiness.json`
   exits 0 **and** its `overall` is `pass` — for `dry_run: true` a `pending`/`fail` overall is
   printed as a warning and the job continues (so pass 1 can rehearse); for `dry_run: false, target: pypi` it is
   fatal. For `target: testpypi` (the rehearsal) a missing or non-pass readiness is a warning, and
   an `rc` version reads the base version's file (`1.0.0rc1` → `validation/release/1.0.0/`).
2. **`build`** (`needs: verify`): `uv build`; `check_dist.py dist --version <version>`; the M11-01
   clean-venv smoke; `scripts/release/notes.py <version>` extracts the `## GIAP <version> (…)`
   section of `CHANGELOG.md` into `notes.md` (exit 1 if absent, except with `--draft`, which `dry_run: true` passes: then the towncrier draft for the version is used, headed `DRAFT`; for `rc` versions the base version's
   section, prefixed "Release candidate"); upload `dist/` and `notes.md` as artefacts.
3. **`publish`** (`needs: build`, `if: !inputs.dry_run && github.actor == 'khoa-yelo'`):
   `environment: ${{ inputs.target }}` (Khoa is the required reviewer, Q-K2);
   `permissions: id-token: write`; `pypa/gh-action-pypi-publish` pinned by SHA, `repository-url`
   set for TestPyPI. Uploads exactly the artefacts `build` produced.
4. **`github-release`** (`needs: publish`, same `if`, and `inputs.target == 'pypi'`):
   `permissions: contents: write`; `gh release create v<version> dist/* --target <github.sha>
   --title "GIAP <version>" --notes-file notes.md` (`--prerelease` for `rc`). This is the only place
   a tag is created. Zenodo archives it through its webhook (Q-K3).
5. **summary**: a step summary table: version, sha, readiness overall, artefact names + sha256,
   what was (not) published.

Third-party actions pinned by commit SHA with the tag in a comment (Q-L5). No secrets are read.

## Examples

| # | Dispatch | Expected |
|---|---|---|
| 1 | lead: `version=0.0.1 dry_run=true` on main after M11-01–03 (no readiness file, `CHANGELOG.md` unbuilt) | verify warns `readiness missing: validation/release/0.0.1/readiness.json` and continues; build passes; `notes.py --draft` finds no `## GIAP 0.0.1` section and uses the `towncrier build --draft --version 0.0.1` text, headed `DRAFT`; run green; summary `dry run: nothing published` |
| 2 | `version=0.0.1 dry_run=false target=pypi` by Khoa, same state | verify fails: `readiness missing` is fatal for a PyPI publish |
| 3 | `version=1.0.0 dry_run=true` while `__version__ = "0.0.1"` | verify fails: `version input 1.0.0 != giap 0.0.1` |
| 4 | `version=1.0.0 dry_run=false target=pypi` by an actor other than `khoa-yelo` | verify and build run; publish and github-release are **skipped**; summary: `publish skipped: actor <x> is not khoa-yelo` |
| 5 | Khoa: `version=1.0.0rc1 dry_run=false target=testpypi` | environment approval prompt; upload to TestPyPI; no tag, no GitHub release |
| 6 | Khoa: `version=1.0.0 dry_run=false target=pypi`, readiness `pass` | upload to PyPI; tag `v1.0.0` on `github.sha`; GitHub release with `notes.md` and both artefacts |
| 7 | as 6 but readiness `pending` | verify fails before build: `readiness overall pending (R3, R9)` |
| 8 | dispatched on a branch `m11/x` | verify fails: `release runs only on main` |
| 9 | tag `v1.0.0` already exists | verify fails: `tag v1.0.0 exists` |

## Edge cases

- **Re-run after a failed publish** (PyPI accepted the wheel, the release step failed): PyPI rejects
  re-uploads; `publish` sets `skip-existing: true` only when the existing file's sha256 equals the
  built one (checked against `https://pypi.org/pypi/giap/<v>/json` before upload); otherwise fail.
- **Readiness for a different commit**: handled by `verify`'s commit rule (M11-03).
- **`rc` versions**: PEP 440 `1.0.0rc1`; tag `v1.0.0rc1` only if Khoa ever targets PyPI with an rc
  (not the default).
- **Concurrent dispatches**: `concurrency: release` with `cancel-in-progress: false`.
- **Empty notes section**: `notes.py` exits 1 (a release with no notes is a mistake).
- Genomic edge cases: not applicable.

## Acceptance criteria

1. A static test parses `release.yml` with PyYAML and asserts: only `workflow_dispatch`; the three
   inputs and defaults; top-level `contents: read`; `id-token: write` only on `publish`;
   `contents: write` only on `github-release`; both guarded by `!inputs.dry_run` and the actor;
   `environment` on `publish`; every `uses:` outside `actions/*` pinned to a 40-hex SHA; no
   `secrets.` reference; no `git push` and no `git tag` anywhere in the file.
2. `notes.py` unit tests: section extraction (first, middle, last section), missing section exit 1,
   draft fallback only with `--draft`, `rc` prefix.
3. One real `dry_run=true` dispatch on main by the lead after merge (`gh workflow run release.yml
   -f version=<current> -f dry_run=true`) is green; its run URL goes into PROGRESS. **This is the
   end of pass 1** together with M11-05.
4. No agent ever dispatches with `dry_run=false` (stated in AGENTS.md by the lead and in M11-05).

## Estimate

Code ≈ 190 (workflow ≈ 140, `notes.py` ≈ 50). Tests ≈ 170.

## Questions for the lead (reversible defaults in use)

- Q-K1 pending publishers on PyPI and TestPyPI (Khoa's accounts). Q-K2 actor + environment guard.
  Q-K3 Zenodo. Q-K4 version string. Q-L5 SHA pins.
- May the lead dispatch `dry_run=true` runs? Default: yes (no tag, no upload; §3.10 forbids pushing
  tags, which a dry run cannot do).
