# M11-18: publish the Bakta JSON and GBFF under `annotation/tools/bakta/`

Milestone M11 · pass 2 · depends on: M2-07 (`BAKTA` process, ANNOTATE subworkflow), M2-20
(`annotation/tools/<tool>/`), M11-08 (user guide) · estimate ≈ 15 workflow lines + docs, one PR ·
roles: implementer (`workflow/modules/local/bakta/main.nf`, `docs/user-guide/`,
`docs/get-started/from-spectra.md`), test-author (`workflow/tests/publish_bakta.nf.test`,
`tests/unit/test_m11_18_publish_bakta.py`) · source: docs/backlog.md, "From-spectra guide gaps",
"The Bakta JSON is not published".

## Spec sections implemented

- §7.1.8 Outputs of Module 1: "`annotation/tools/<tool>/`: native outputs, untouched."
- §9.1 output tree (`docs/user-guide/outputs.md`: `annotation/tools/  native outputs of each tool`).
- §6.3 process conventions: `publishDir` copies final outputs into `OUT/` (M2-06).

## The gap

`docs/user-guide/outputs.md` lists `annotation/tools/` ("bakta/, uniop/, genomad/, …"), but no
tool process has a `publishDir`. Only `GIPPA_STEP` publishes, so runs have no
`annotation/tools/`. The pilot's pass-3 `out_a2` is one example. `gippa makedb` needs
`--annotation` (exit 2 without it), so the from-spectra routes tell users to run Bakta by hand,
even after a `gippa run` already ran Bakta on the same genome.

## Published paths

When the `BAKTA` process runs (no `annotation.precomputed.bakta_json`):

| File | Source (BAKTA task) |
|---|---|
| `OUT/annotation/tools/bakta/<sample>.json` | `bakta/<sample>.json` (emit `json`) |
| `OUT/annotation/tools/bakta/<sample>.gbff` | `bakta/<sample>.gbff` (emit `gbff`) |

- `<sample>` is the BAKTA `meta_id`, the Bakta `--prefix`. ANNOTATE sets it to the genome input's
  `simpleName`. Under `gippa run` the genome is the normalised `work/genome.norm.fa`, so the files
  are `annotation/tools/bakta/genome.json` and `genome.gbff`.
- Only the two declared outputs are published. `versions.yml` is never published. Bakta also
  writes `.gff3`, `.faa`, `.tsv` and other files, which are not declared outputs and are not
  published. The JSON has every feature, the GBFF serves antiSMASH and GenBank users, and copying
  the rest would only add size. A later item can add files by declaring them as outputs.
- **Mode `copy`**, the same as `GIPPA_STEP` (`publishDir "${params.outdir}", mode: 'copy'`) and
  `real.config`. A link into `work/` would break when `work/` is cleaned. The published files are
  regular files that are byte-identical to the task outputs (§7.1.8 "untouched": no gzip, no
  rename).
- Put the directive on the process, as `GIPPA_STEP` does. For example:
  `publishDir "${params.outdir}/annotation/tools", mode: 'copy', saveAs: { fn -> fn == 'versions.yml' ? null : fn }`.
  The task already writes the outputs as `bakta/<sample>.*`, so this gives
  `annotation/tools/bakta/<sample>.*`. A `withName: BAKTA` block in a config file is also
  acceptable, as long as it applies in every profile.
- The local executor (`--executor local`) needs `precomputed.bakta_json` (`REQUIRED_ROOTS`) and
  never runs Bakta, so it publishes nothing here. No change to `gippa.runner.local.publish`.

## `precomputed.bakta_json`: not published (recommendation)

When `annotation.precomputed.bakta_json` is set, BAKTA does not run and `OUT/annotation/tools/bakta/`
is **not** created. The precomputed JSON (and the GBFF next to it) is not copied into `OUT/`.

- `annotation/tools/` holds the native outputs of tools that **ran** (M2-20 example 1: "contains
  only tools that ran (none with full precomputed)"). A copied user file there would suggest that
  GIPPA produced it.
- The user already has the file and gives the same path to `makedb --annotation`. Its SHA-256 is
  in `run_manifest.json` (input checksums), so provenance does not depend on a copy.
- TE7 (M10-06) classes `annotation/tools/**` as `external`. Copying an input there adds no
  information.

The reversible alternative is to copy it to `annotation/tools/bakta/` with its own name. That
would make `OUT/` self-contained, but it duplicates the file and blurs provenance. It is not
recommended. If Khoa prefers it, the second nf-test below changes.

## Documentation

- `docs/user-guide/outputs.md`: the `tools/` line lists only tools that are published. After
  this item that is Bakta only: `bakta/<sample>.json`, `<sample>.gbff`, written when Bakta ran.
  `uniop/` and `genomad/` are removed from the line, because neither publishes yet (backlog:
  publish the other tools' native outputs). An undefined "…" is not enough.
- `docs/user-guide/makedb.md` and `docs/get-started/from-spectra.md`: after a `gippa run`, use
  `gippa makedb --annotation OUT/annotation/tools/bakta/genome.json`. Running Bakta by hand is
  still documented for users who build the database before any run.
- Changelog fragment `changelog.d/<PR>.added.md`.

## Acceptance criteria

| # | Check | Test |
|---|---|---|
| 1 | Stub run of `main.nf` (genome + replicons, nothing precomputed): `OUT/annotation/tools/bakta/genome.json` and `genome.gbff` exist, are regular files (not symlinks), and the JSON is byte-identical to `workflow/assets/stub/bakta.json` (the BAKTA stub output) | `publish_bakta.nf.test` test 1 |
| 2 | Same run: `OUT/annotation/tools/bakta/` holds exactly `genome.gbff` and `genome.json`. No `versions.yml` anywhere under `annotation/tools/` | test 1 |
| 3 | `precomputed.bakta_json` set (with a GBFF next to it): BAKTA does not run and `OUT/annotation/tools/bakta/` does not exist | `publish_bakta.nf.test` test 2 (regression guard: it already passes before the change) |
| 4 | Every tool that the `tools/` line of `outputs.md` names (`<tool>/`) has a `publishDir` to `annotation/tools` (on its module or in a `withName` config block), and the line names `bakta/` | `test_m11_18_publish_bakta.py` |
| 5 | `makedb.md` and `from-spectra.md` give `annotation/tools/bakta/` as the `--annotation` source after a run | `test_m11_18_publish_bakta.py` |
| 6 | Existing nf-test suites still pass (`--profile test,stub`). Publishing adds files only under `annotation/tools/` | CI `nf-test-stub` |
