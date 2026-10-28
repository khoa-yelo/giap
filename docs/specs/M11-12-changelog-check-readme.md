# M11-12: the changelog-fragment check and the 1.0 README

Milestone M11 · **pass 2** · depends on: M11-02 (fragment naming rules) · estimate ≈ 120 code lines
+ ≈ 130 test lines + ≈ 80 prose lines · not `stats` · roles: implementer
(`scripts/ci/check_changelog.py`, its step in `ci.yml` › guards), docs (`README.md`), test-author
(`tests/unit/test_m11_12_*.py`) · Khoa: none.

## Spec sections implemented

- §14.2: "Every PR links an issue (or a `docs/backlog.md` entry) and the spec section it implements,
  and includes a changelog fragment." Today 2 of 317 merged PRs have one.
- §12 Releases "changelog fragments (towncrier)". §4 `README.md`.

## Definition

1. **`scripts/ci/check_changelog.py BASE`** (stdlib + git; same style as `check_roles.py`):
   - changed files = `git diff --name-only BASE...HEAD`;
   - *needs a fragment* when any changed path is under `src/` or `workflow/` (excluding
     `workflow/tests/`);
   - *has a fragment* when the diff adds at least one file under `changelog.d/` whose name matches
     M11-02's rule (`<PR number>.<type>.md` or `+<slug>.<type>.md`) and is non-empty;
   - exempt when the PR carries the label `no-changelog` (read from `GITHUB_EVENT_PATH`'s
     `pull_request.labels`; absent outside Actions → not exempt) or every commit subject in the
     range has conventional type `test`, `ci`, `chore` or `refactor`, with or without a scope
     (`chore(release): 1.0.0` counts; no user-visible change);
   - exit 0 / 1 with `changelog fragment missing: changed src/… ; add changelog.d/<PR>.<type>.md or
     label no-changelog`.
   Severity: fail (Q-L9). The lead announces it in AGENTS.md (one line) when it lands.
2. **README for 1.0**: status line from "under development (milestone M0)" to the released state
   (version badge from PyPI, DOI badge after Zenodo, docs link); "Install" (pip, Bioconda), a
   four-line quickstart pointing at the site, "Citing GIPPA" (points at `CITATION.cff`),
   "Licence" (unchanged), "Development on FarmShare" moved under a "Developers" heading. No numbers
   that are not in committed evidence.

## Examples

| # | Input (PR diff) | Expected |
|---|---|---|
| 1 | `src/gippa/cli.py` + `changelog.d/341.feature.md` | exit 0 |
| 2 | `src/gippa/cli.py` only | exit 1 with the message |
| 3 | `docs/user-guide/run.md` only | exit 0 (no code changed) |
| 4 | `workflow/tests/x.nf.test` only | exit 0 |
| 5 | `src/gippa/cli.py` + label `no-changelog` | exit 0, prints `exempt: label no-changelog` |
| 6 | `src/gippa/x.py`, all commits `refactor: …` | exit 0, prints `exempt: refactor/test/ci/chore only` |
| 7 | `src/gippa/x.py` + an empty `changelog.d/342.bugfix.md` | exit 1 `empty fragment` |
| 8 | `src/gippa/x.py` + `changelog.d/notes.txt` | exit 1 (name does not match) |

## Edge cases

- **Fragment deleted** in the PR (a release-prep build removes all fragments): the release-prep PR
  touches `src/gippa/__init__.py` (version) — exempt by its `chore(release):` subject (rule above).
- **Merge commits from main** in the range: `BASE...HEAD` (three dots) compares with the merge base,
  so main's changes do not count.
- **Renamed files** count as changed at their new path.
- Genomic edge cases: not applicable.

## Acceptance criteria

1. Examples 1–8 are unit tests on throwaway git repos in `tmp_path` (the `check_roles.py` test
   pattern) with a fake event JSON for labels.
2. The CI step runs only on `pull_request` and is listed in `guards`.
3. README renders on GitHub with no broken relative links (a test resolves every relative link to an
   existing file).

## Estimate

Code ≈ 120 (script ≈ 100, ci.yml ≈ 5, README not counted). Tests ≈ 130. Prose ≈ 80.

## Questions for the lead (reversible defaults in use)

- Q-L9 fail vs. warn. Default: fail, with the label and the conventional-type exemptions.
- Should `docs:` commits that change user-guide pages need a fragment (`doc` type)? Default: no
  (they do not touch `src/` or `workflow/`).
