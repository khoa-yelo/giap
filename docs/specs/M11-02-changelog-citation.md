# M11-02: the 1.0 changelog content and a complete citation file

Milestone M11 · **pass 1 (MVP)** · depends on: none (merge after M11-01: both touch `pyproject.toml`,
different tables) · estimate ≈ 30 code lines + ≈ 160 test lines + ≈ 150 prose lines · not `stats` ·
roles: docs (`changelog.d/*.md`, `CITATION.cff`), implementer (`[tool.towncrier]` in
`pyproject.toml`), test-author (`tests/unit/test_m11_02_*.py`) · Khoa: Q-K7 (authors).

## Spec sections implemented

- §13 row M11 "changelog, citation file"; §12 `release.yml` "verify … changelog and the version";
  §12 **Releases** "changelog fragments (towncrier) … a Zenodo DOI each time".
- §4 `CHANGELOG.md`, `CITATION.cff`; §14.1 docs role owns `CHANGELOG.md`; §14.2 "Every PR …
  includes a changelog fragment".
- R12 (no personal email addresses in the citation file).

## What exists

`towncrier build --draft --version 1.0.0 --date 2026-11-01` on main renders today:

```markdown
# giap 1.0.0 (2026-11-01)

## Features

- `validation/fetch_regulondb.py` downloads … truth sets of GIAP_SPEC.md §15.4. (m6-11-regulondb)
- Project scaffold: FarmShare environment … the simulator skeleton. (m0-scaffold)
```

Three problems: 2 fragments for 317 PRs; the fragment *name* is printed as an issue reference
(`(m0-scaffold)`); the order is alphabetical by name, so M6 precedes M0.

## Definition

1. **Rendering rules** (`[tool.towncrier]`): `name = "GIAP"`, `title_format = "## GIAP {version} ({project_date})"`
   (one level below the file's `# Changelog`). Milestone
   fragments are named `+<slug>.<type>.md` (towncrier orphan prefix) so no reference is printed;
   PR fragments named `<PR number>.<type>.md` render `([#N](https://github.com/khoa-yelo/giap/pull/N))`
   via `issue_format = "[#{issue}](https://github.com/khoa-yelo/giap/pull/{issue})"`. Types, in
   order: `feature` "Features", `bugfix` "Fixes", `doc` "Documentation", `removal` "Removals and
   breaking changes", `misc` "Other" (shown with content). The two existing fragments are renamed
   to `+m00-scaffold.feature.md` and `+m06-regulondb.feature.md`.
2. **Retroactive fragments** (Q-L3): one `+m<k>-<slug>.feature.md` per milestone M1–M10 (M0 and the
   RegulonDB one exist), each 3–8 lines, user-facing (what a user can now do; command names, output
   files), written from `PROGRESS.md` › Done and the milestone index, no internal issue numbers,
   no metrics that are not in committed evidence. Plus one `+known-limitations.misc.md` listing
   what §16 defers that users would expect (several samples, V4/V9, RO-Crate, full validation) and
   the open non-gate issues at release time (Q-K12; M11-13 updates it).
   Slugs sort in milestone order: `m00-…`, `m01-…`, …, `m10-…` (two digits), so towncrier's
   alphabetical order is the milestone order.
3. **`CITATION.cff`** (cff-version 1.2.0): `type: software`; `title`; `abstract` (≤ 120 words, the
   README's first paragraph); `authors` (Q-K7 default: Khoa Hoang only, `given-names`/`family-names`,
   no email); `license: MIT`; `repository-code: https://github.com/khoa-yelo/giap`; `url` (the docs
   site, Q-K9; omitted until Pages exists); `keywords` (proteogenomics, bacteria, genome context,
   operons, target-decoy, Nextflow); `version` = `giap.__version__`; `date-released` present only
   from the release-prep PR on (M11-13); `identifiers` (the Zenodo concept DOI) added after the first
   release (M11-13 post-release PR).
4. **Agreement test.** `CITATION.cff › version` equals `giap.__version__`; if `date-released` is
   present it is a valid ISO date and `CHANGELOG.md` has a `## GIAP <version> (<that date>)`
   heading.
5. **Schema check.** The CFF 1.2.0 JSON Schema (`schema.json` from
   `github.com/citation-file-format/citation-file-format`, tag `1.2.0`, CC-BY-4.0) is committed under
   `tests/data/schemas/cff-1.2.0.schema.json` with a `PROVENANCE.md` (URL, tag, sha256, licence,
   date) and validated with `jsonschema` (already in the venv). No network in tests.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `towncrier build --draft --version 1.0.0 --date 2026-11-01` after this issue | starts `## GIAP 1.0.0 (2026-11-01)`, then a Features heading one level below it (the test pins the exact level towncrier renders for that title) with 12 bullets in M0 → M10 order (one per milestone plus the RegulonDB one, which sorts after M6's), no `(m0-scaffold)` suffix; then an Other section with the known-limitations bullet |
| 2 | a later PR adds `331.bugfix.md` "Fix X." | the draft gains a Fixes section › `- Fix X. ([#331](https://github.com/khoa-yelo/giap/pull/331))` |
| 3 | `CITATION.cff` with `version: 0.0.1`, `__version__ = "0.0.1"` | agreement test passes |
| 4 | `__version__` bumped to `1.0.0`, CFF still `0.0.1` | agreement test fails: `CITATION.cff version 0.0.1 != giap 1.0.0` |
| 5 | CFF with `date-released: 2026-11-01` but no `## GIAP 1.0.0 (2026-11-01)` in `CHANGELOG.md` | agreement test fails naming the missing heading |
| 6 | CFF with an `email:` key under an author | a test fails (R12) |

## Edge cases

- **Fragment with an unknown type** (`x.feat.md`) — towncrier ignores it silently; a test lists
  `changelog.d/` and fails on any file not matching `^(\+[a-z0-9-]+|\d+)\.(feature|bugfix|doc|removal|misc)\.md$`
  (README files excepted: none exist).
- **Empty fragment** — test fails (a blank bullet).
- **Ties in order** (two fragments of one milestone) — alphabetical by slug; acceptable.
- **`CHANGELOG.md` already has the version** (re-running a build) — towncrier refuses unless
  `--yes`; the release prep (M11-13) builds once.
- Genomic edge cases: not applicable.

## Acceptance criteria

1. Example 1's draft output is reproduced by a test that runs towncrier in `tmp_path` on a copy of
   `changelog.d/` and `pyproject.toml` (towncrier is in the venv; no network).
2. Every milestone M0–M10 has at least one `+m<kk>-*.feature.md` fragment (M6 has two: operons and RegulonDB).
3. `CITATION.cff` validates against the committed CFF 1.2.0 schema and has every item-3 key that
   is due now (`date-released`, `identifiers`, `url` are optional until M11-13).
4. Tests for examples 3–6 and the fragment-name rule pass.
5. No change to `CHANGELOG.md` itself in this PR (it is written once, by the release prep).

## Estimate

Code ≈ 30 (pyproject table). Tests ≈ 160. Prose ≈ 150 (12 fragments, CFF).

## Questions for the lead (reversible defaults in use)

- Q-K7 (Khoa) authors/ORCID/affiliation/AI acknowledgement. Default as in item 3.
- Q-L3 one fragment per milestone. Alternative: one per merged PR (317 lines of noise).
- Vendoring the CFF schema (CC-BY-4.0, attribution in PROVENANCE.md) vs. checking keys by hand.
  Default: vendor; if the lead prefers no third-party file, a hand-written required-key check.
