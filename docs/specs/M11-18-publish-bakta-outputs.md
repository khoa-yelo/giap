# M11-18: publish the Bakta JSON and GBFF under `annotation/tools/bakta/`

Milestone M11 · pass 2 · depends on: M2-07 (`BAKTA` process, ANNOTATE subworkflow), M2-20
(`annotation/tools/<tool>/`), M11-08 (user guide) · estimate ≈ 15 workflow lines + docs, one PR ·
roles: implementer (`workflow/modules/local/bakta/main.nf`, `workflow/subworkflows/local/annotate.nf`, `docs/user-guide/`,
`docs/get-started/from-spectra.md`), test-author (`workflow/tests/publish_bakta.nf.test`,
`tests/unit/test_m11_18_publish_bakta.py`, `tests/unit/test_m11_18_bakta_tool_dir.py`) · source: docs/backlog.md, "From-spectra guide gaps",
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
- **Route: the M2-20 `tool_dirs` input of `context_layers`** (lead decision 2026-10-28, below).
  BAKTA declares a third output, `path 'tools/bakta', emit: tool_dir`. In the script and the stub
  it runs `mkdir -p tools/bakta` and hard-links `bakta/<sample>.json` and `bakta/<sample>.gbff`
  into it (`ln`, with a `cp` fallback when `ln` fails across filesystems). Only those two files
  go in: no symlinks, no `bakta/*`, no `cp -r`. ANNOTATE passes `kv('tool_dirs', BAKTA.out.tool_dir)`
  to the `context_layers` step, only in the branch where BAKTA runs. `copy_native`
  (`gippa.annotate.layers`) byte-copies the directory to `annotation/tools/bakta/`, keyed by the
  directory name, and CONTEXT_LAYERS publishes it with the rest of `annotation/` (mode `copy`).
  The published files are regular files, byte-identical to the task outputs (§7.1.8 "untouched":
  no gzip, no rename).
- BAKTA has **no `publishDir`**, neither on the process nor in a `withName` block.

### Lead decision 2026-10-28: `tool_dirs`, not a BAKTA `publishDir`

The first design put `publishDir "${params.outdir}/annotation/tools", mode: 'copy'` on BAKTA.
That cannot work. CONTEXT_LAYERS (a `GIPPA_STEP` with out item `annotation=annotation`) always
runs after BAKTA and publishes the whole `annotation` output to `OUT/annotation`. Nextflow's
default `overwrite` deletes an existing destination before it copies, so
`OUT/annotation/tools/bakta/` is removed again. This was checked with a two-process toy pipeline on
Nextflow 26.10.0. In stub runs, `annotation` is a touched file and replaces the directory in the
same way. Publishing `annotation/` file by file would need a different `GIPPA_STEP` output shape:
Nextflow does not publish files below a directory that is itself a declared output. The
`tool_dirs` route already exists (M2-20, `copy_native`; its annotate.nf wiring was backlog item
13 there) and keeps `annotation/` as the only publisher of `OUT/annotation/`.

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
- Changelog fragment `changelog.d/468.feature.md` (towncrier types: feature, bugfix, doc, removal, misc).

## Acceptance criteria

| # | Check | Test |
|---|---|---|
| 1 | Stub run of `main.nf` (genome + replicons, nothing precomputed): the BAKTA task has `tools/bakta/` with exactly `genome.gbff` and `genome.json`, regular files (not symlinks), byte-identical to `bakta/genome.*`; the JSON equals `workflow/assets/stub/bakta.json`. In stub runs the published `annotation` is a touched file, so the check is on the work dir | `publish_bakta.nf.test` test 1 |
| 2 | Same run: the CONTEXT_LAYERS task has one `tool_dirs` input, the BAKTA `tools/bakta` directory (named `bakta`). No Bakta file is published by BAKTA itself | test 1 |
| 2b | `context_layers` with `tool_dirs` = a recorded Bakta 1.12.2 `tools/bakta` writes `annotation/tools/bakta/genome.{json,gbff}` byte-identical, and nothing else in `annotation/` changes | `test_m11_18_bakta_tool_dir.py` |
| 3 | `precomputed.bakta_json` set (with a GBFF next to it): BAKTA does not run and `OUT/annotation/tools/bakta/` does not exist | `publish_bakta.nf.test` test 2 (regression guard: it already passes before the change) |
| 4 | Every tool that the `tools/` line of `outputs.md` names (`<tool>/`) reaches `annotation/tools` (a `publishDir` on its module or in a `withName` block, or a `tools/<tool>` output wired as `tool_dirs`), and the line names `bakta/`. BAKTA has no `publishDir`, emits `tools/bakta` with only the JSON and GBFF, and annotate.nf wires it in the BAKTA branch before the `context_layers` stepIn | `test_m11_18_publish_bakta.py` |
| 5 | `makedb.md` and `from-spectra.md` give `annotation/tools/bakta/` as the `--annotation` source after a run | `test_m11_18_publish_bakta.py` |
| 6 | Existing nf-test suites still pass (`--profile test,stub`). The tool dir adds files only under `annotation/tools/` | CI `nf-test-stub` |
