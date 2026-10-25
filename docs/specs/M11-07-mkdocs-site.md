# M11-07: the MkDocs site, generated references and `docs.yml`

Milestone M11 · **pass 2** · depends on: M11-01 (soft) · estimate ≈ 260 code lines + ≈ 180 test
lines · not `stats` · roles: implementer (`mkdocs.yml`, `scripts/docs/gen_reference.py`,
`.github/workflows/docs.yml`, the `docs` CI step), docs (`docs/index.md`, nav stubs, link fixes in
`docs/adr/`), statistician (approves any edit inside `docs/methods/`), test-author
(`tests/unit/test_m11_07_*.py`) · Khoa: Q-K9 (enable Pages).

## Spec sections implemented

- §13 row M11 "MkDocs site"; §4 `docs/` = "MkDocs: user guide, methods/ (one note per statistical
  method), adr/, specs/"; `pyproject` extra `docs: mkdocs-material` (§13.1 step 5).
- Appendix A.5 docs agent: "Document every CLI option and output table from the code and schemas,
  not from memory." → the references are **generated**, and CI fails when they drift.
- §6.4: `--executor local` "is documented under 'Developer guide', not in the user guide".
- §9.5 schemas; §5.1 CLI; §5.2 configuration defaults.

## Definition

1. **`mkdocs.yml`** at the repository root: `site_name: GIAP`, `docs_dir: docs`, `theme: material`
   (no external fonts: `font: false`, so the site makes no Google Fonts request),
   `exclude_docs: specs/  backlog.md` (Q-L6), `strict` builds in CI, `validation` of links and
   anchors at `warn` (strict turns warnings into errors), `markdown_extensions`: `tables`, `toc`
   (permalink), `admonition`, `pymdownx.superfences`, `pymdownx.arithmatex: {generic: true}` with
   MathJax 3 from jsDelivr in `extra_javascript` (Q-L7).
2. **Nav** (pages that M11-08 fills are created here as stubs with a one-line purpose):

   ```text
   Home (index.md)
   Get started (docs/get-started/): install, setup-databases, quickstart
   User guide (docs/user-guide/): inputs, run, outputs, interpreting-results, plots-and-report, compare, makedb,
               convert, licences, platforms (M11-09), runtimes (M10-12's runtimes.md)
   Tutorials: ecoli-k12, nissle-pilot
   Reference: cli (generated), configuration (generated), tables (generated), flags (generated)
   Methods: index (generated list) + the 29 docs/methods notes
   Decisions: index (generated list) + docs/adr/ADR-*.md (TEMPLATE.md excluded)
   Developer guide: developer (local executor, tests, FarmShare), release (M11-05)
   Changelog (CHANGELOG.md, included via a page that embeds it), Citation (CITATION.cff rendered)
   ```
3. **`scripts/docs/gen_reference.py [--check]`** writes, from code only:
   - `docs/reference/cli.md`: every command and option from the Click tree (`giap.cli`), in command
     order, with type, default, choices and help text; hidden commands (`_step`, `version` if
     hidden) excluded;
   - `docs/reference/configuration.md`: every key of `defaults.yaml` (dotted path, default, type and
     description from `config.schema.json`), grouped by top-level section;
   - `docs/reference/tables.md`: one section per resource in `schemas/datapackage.json`: path,
     description, coordinate system, then a field table (name, type, unit, description,
     constraints);
   - `docs/reference/flags.md`: the §9.4 flag vocabulary from the schema enum;
   - `docs/methods/index.md` and `docs/adr/index.md`: title lists (first `#` heading of each file;
     ADR status line).
   Output is deterministic (sorted where the source has no order; LF; trailing newline).
   `--check` rewrites into a temp dir and exits 1 with a unified diff when the committed files
   differ (the pattern of `write_config_schema.py --check`).
4. **CI**: in `ci.yml` › lint, `python scripts/docs/gen_reference.py --check` and
   `mkdocs build --strict --site-dir /tmp/site`.
   **`docs.yml`**: on push to main (paths `docs/**`, `mkdocs.yml`, `src/giap/cli.py`,
   `src/giap/config/**`, `src/giap/schemas/**`, `CHANGELOG.md`) and dispatch: build strict, then
   `actions/upload-pages-artifact` + `actions/deploy-pages` (SHA-pinned), `permissions: pages: write,
   id-token: write` on the deploy job, `environment: github-pages`. The deploy job runs only when
   `vars.PAGES_ENABLED == 'true'` (Khoa sets the variable when enabling Pages; Q-K9).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `gen_reference.py` on main | `docs/reference/cli.md` has sections `giap compare`, `giap convert`, `giap makedb`, `giap preflight`, `giap report`, `giap run`, `giap setup`; under `giap run` a row `--genome` `FILE` required; `--profile` default `farmshare,apptainer` |
| 2 | a PR adds option `--foo` to `giap run` without regenerating | CI lint fails: `docs/reference/cli.md is stale` + diff with the `--foo` row |
| 3 | `configuration.md` | contains a row `seed` with the value in `defaults.yaml` and its schema description |
| 4 | `tables.md` › `genes.tsv` | coordinate line "1-based closed" (from the schema), field `locus_tag` first |
| 5 | `mkdocs build --strict` with a methods note containing `[0, 1](1 − q)` (parsed as a link) | build fails with the link warning; the fix is an escaped bracket or code span, made by the docs role and approved by the statistician |
| 6 | `docs/specs/M3-10-placement-fdr.md` | not in the built site (`/tmp/site/specs/` absent) |
| 7 | push to main with `PAGES_ENABLED` unset | build job green; deploy job skipped |

## Edge cases

- **Links from docs to repo files** (`../src/giap/…`): outside `docs_dir` they break `--strict`;
  rewrite as absolute GitHub URLs on `main` (docs role). ADRs that cite specs (`docs/specs/…`, now
  excluded) link to the GitHub URL instead.
- **Click options with dynamic defaults** (e.g. thread count from the environment): printed as the
  help text says, never the build machine's value; a test sets env vars and checks the output is
  unchanged.
- **Unicode** (§, ≤, en dashes) in help strings: written as-is in UTF-8.
- **Schema fields without a description**: written as `—` and listed by `--check` as a warning (not
  an error), so missing descriptions become visible without blocking.
- **Empty methods/ADR directories**: an index with "No entries" (only in tests).
- Genomic edge cases: not applicable (the coordinate system is copied from the schema).

## Acceptance criteria

1. `mkdocs build --strict` passes locally (in the allocation, venv `mkdocs`) and in CI.
2. `gen_reference.py` unit tests: a small Click group and a small datapackage in `tmp_path`
   produce the expected Markdown (golden strings in the test); `--check` exit 0/1 behaviour;
   determinism (two runs byte-identical).
3. Every `giap` subcommand and every non-hidden option appears in `cli.md`; every `defaults.yaml`
   key in `configuration.md`; every datapackage resource in `tables.md` (tests walk the real
   sources).
4. `docs.yml` static test: triggers, the `PAGES_ENABLED` guard, permissions only on deploy, pins.
5. No page in the built site requests a third-party host other than jsDelivr (MathJax): a test
   greps `/tmp/site/**/*.html` for `https?://` in `src=`/`href=` of `<script>`/`<link>`.

## Estimate

Code ≈ 260 (`gen_reference.py` ≈ 190, `mkdocs.yml` ≈ 45, `docs.yml` ≈ 25). Tests ≈ 180. Prose: stubs.

## Questions for the lead (reversible defaults in use)

Q-K9 (Pages), Q-L6 (ADRs in, specs out), Q-L7 (MathJax via jsDelivr). Also: put the generated
pages in git (default, so reviewers see diffs) rather than generating at build time.
