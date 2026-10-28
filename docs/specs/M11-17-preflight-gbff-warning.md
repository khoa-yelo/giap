# M11-17: preflight warns when a precomputed Bakta JSON has no GBFF next to it (bgc module)

Milestone M11 · pass 2 · depends on: M1 preflight (`gippa.preflight`, ADR-0005 item 3), M2-18
(antiSMASH on the Bakta GBFF), MVP-03 (manifest `warnings`) · estimate ≈ 60 code lines + ≈ 250 test
lines, one PR · roles: implementer (`src/gippa/preflight/`), test-author
(`tests/unit/test_m11_17_*.py`) · source: backlog "Preflight GBFF warning" (2026-10-27, lead).

## Spec sections implemented

- §7.0 Preflight ("every evaluated rule yields one check"; `warn` never stops a run; `--strict`
  turns warnings into exit 3, M1-12), §5.3 (a Bakta JSON given as `--annotation` is the
  precomputed Bakta layer), §7.1.1 (antiSMASH runs on the Bakta GBFF), §9.2 manifest `warnings`.

## The gap

With `annotation.precomputed.bakta_json` (or `--annotation <x>.json[.gz]`) Bakta does not run, so
there is no Bakta GBFF for antiSMASH. `workflow/subworkflows/local/annotate.nf` then looks for
`<stem>.gbff`, else `<stem>.gbff.gz`, next to the JSON (`<stem>` = the JSON path with a final
`.json` or `.json.gz` removed, any case). When neither exists, antiSMASH is skipped and Nextflow
only prints `WARN: layer_missing:antismash: …` on its console. The run then finishes with an
empty bgc layer, and nothing in `preflight.json`, `preflight.md` or `run_manifest.json` says why.
The user learns it after hours of compute, not before.

## Rule

`rule_id = "inputs.bakta_gbff"`, `input = "inputs"`. New module `gippa/preflight/gbff.py`:

```python
RULE_ID = "inputs.bakta_gbff"
def expected_gbff(bakta_json: Path) -> tuple[Path, Path]: ...       # pure, no I/O
def bgc_gbff_check(cfg: GippaConfig, annotation: Path | None) -> Check | None: ...
```

**The Bakta JSON** of the run is `cfg.annotation.precomputed.bakta_json` when set, else
`annotation` (`PreflightInputs.annotation`) when its name ends in `.json` or `.json.gz` (any case;
the same test as `gippa.cli._is_bakta_json`). `gippa run` puts `--annotation x.json` into the
config; `gippa preflight` does not, so the second source makes both commands agree.

**`expected_gbff(p)`** returns `(Path(s + ".gbff"), Path(s + ".gbff.gz"))` with
`s = re.sub(r"(?i)\.json(\.gz)?$", "", str(p))`, the Nextflow rule exactly. The GBFF suffixes are
lower case (Nextflow checks exactly these names). A relative path stays relative (resolved
against the working directory, as Nextflow resolves it against the launch directory).

| Condition (first match wins) | Level | `value` | `message` |
|---|---|---|---|
| no Bakta JSON (none configured; `annotation` absent or a GFF3/GBFF) | — | — | no check: the rule is not applicable (Bakta runs and writes its GBFF) |
| `"bgc" not in cfg.modules` | `info` | `null` | `not evaluated: the bgc module is off` |
| `cfg.annotation.precomputed.regions` is set | `info` | `null` | `not evaluated: annotation.precomputed.regions replaces antiSMASH` |
| `<stem>.gbff` is a file | `ok` | see below, `gbff` = that path | `""` |
| else `<stem>.gbff.gz` is a file | `ok` | see below, `gbff` = that path | `""` |
| otherwise | `warn` | see below, `gbff` = `null` | the wording below |

`value` (ok / warn) = `{"bakta_json": str(json), "expected": [str(gbff), str(gbff_gz)], "gbff": str | null}`,
paths as given (not resolved, not shortened to a file name: the user must find the file).

`rule` column text: `bgc module on and a precomputed Bakta JSON <stem>.json[.gz] -> <stem>.gbff or
<stem>.gbff.gz must exist next to it, else warn (antiSMASH is skipped: layer_missing)`.

**Warning wording** (`{json}`, `{gbff}`, `{gbff_gz}` = the paths of `value`):

> no GBFF next to the precomputed Bakta JSON {json}: expected {gbff} (or {gbff_gz}); antiSMASH will
> be skipped (layer_missing:antismash) and the bgc module will report no biosynthetic gene
> clusters. Put Bakta's GBFF output there, or leave bgc out of --modules.

The level is only `ok`, `info` or `warn`; `exit_code` is 0. It **never errors and never changes
results**: `NormalisedInputs` is unchanged and no step reads the check. Under `--strict` a `warn`
gives exit 3, as every preflight warning does.

## Where it runs and where it surfaces

1. `run_preflight` calls `bgc_gbff_check(cfg, inp.annotation)` and, when it returns a check, puts
   it directly after the `inputs.*` checks (after `inputs.annotation_kind` when present) and
   before `config.translation_table`.
2. **Preflight table**: `preflight.json` `checks[]` and `preflight.md`, through the existing
   writers (no new writer code). The rule id matches the schema pattern `^[a-z]+\.[a-z_]+$`.
3. **Manifest**: `run_manifest.json` `warnings` contains `preflight:inputs.bakta_gbff` through the
   existing `gippa.manifest._warnings` (only a test).
4. The Nextflow console `WARN: layer_missing:antismash` stays as it is (no workflow change).

## Examples (input → expected output)

`expected_gbff`:

1. `/d/x.json` → `(/d/x.gbff, /d/x.gbff.gz)`.
2. `/d/x.json.gz` → `(/d/x.gbff, /d/x.gbff.gz)`.
3. `/d/sample.bakta.JSON.GZ` → `(/d/sample.bakta.gbff, /d/sample.bakta.gbff.gz)` (only the final
   suffix is removed, any case).
4. `rel/x.Json` → `(rel/x.gbff, rel/x.gbff.gz)` (stays relative).

`bgc_gbff_check` (default modules `core,mobile,bgc,function`; `tmp/` a fresh directory):

5. `precomputed.bakta_json = tmp/x.json`, `tmp/x.gbff` exists → `ok`, `value.gbff = "tmp/x.gbff"`.
6. Same, only `tmp/x.gbff.gz` exists → `ok`, `value.gbff = "tmp/x.gbff.gz"`.
7. Both exist → `ok`, `value.gbff` is the plain `.gbff` (Nextflow's order).
8. Neither exists → `warn`; the message contains `tmp/x.gbff` and `tmp/x.gbff.gz` (full paths) and
   `layer_missing:antismash`; `value.gbff = null`.
9. Example 8 with `modules = core,function` → `info`, `not evaluated: the bgc module is off`.
10. Example 8 with `annotation.precomputed.regions = tmp/regions.tsv` → `info`.
11. No precomputed JSON and `annotation = None` → `None`; with `annotation = tmp/a.gff3` → `None`.
12. No precomputed JSON, `annotation = tmp/x.json.gz` and no GBFF → `warn` naming `tmp/x.gbff`.
13. `precomputed.bakta_json = tmp/a.json` (with `tmp/a.gbff`) and `annotation = tmp/b.json`
    (without) → `ok`: the configured JSON wins (as `--set` wins in `gippa run`).
14. A GBFF with another case (`tmp/x.GBFF`) does not count on a case-sensitive file system: `warn`.

## Acceptance criteria

1. Examples 1–14 hold. The level is never `error` and `exit_code` is always 0.
2. `run_preflight` has exactly one `inputs.bakta_gbff` check when a Bakta JSON is configured or
   given, none otherwise, placed as in "Where it runs" item 1; it equals `bgc_gbff_check(cfg,
   inp.annotation)`.
3. A non-strict preflight with the warning exits 0 and its `NormalisedInputs` equals the one of the
   same run with a GBFF present (only `checks` differ); `exit_code(report, strict=True)` is 3.
4. `preflight.json` validates against `preflight.schema.json` with the new check.
5. `build_manifest` with a preflight that has the warning lists `preflight:inputs.bakta_gbff` in
   `warnings`; at `ok` or `info` it does not.
6. `gippa preflight --annotation x.json --profile none` without a GBFF prints the warning (with
   the expected path) in `preflight.md` and exits 0; with `--strict` it exits 3; with
   `x.gbff` next to it, it exits 0 and the warning is absent.
7. Property: for any file name stem and any case of `.json` / `.json.gz`, `expected_gbff` keeps the
   directory and the stem and appends exactly `.gbff` / `.gbff.gz`.
8. Tests use no network, no `$GIPPA_DB_DIR` and no `$GIPPA_ROOT`. Lint, format and mypy are clean.

No statistical claim, so no calibration test.

## Test plan

- `tests/unit/test_m11_17_gbff_rule.py`: examples 1–14, AC 1, AC 7 (Hypothesis).
- `tests/unit/test_m11_17_gbff_surfaces.py`: AC 2–6 (`run_preflight` with a stubbed `EnvProbe` as
  in the M11-14w tests; manifest through `helpers.manifest_data`; the CLI through `CliRunner`).

## Questions for the lead

- Q1: the message and `value` hold the full path as given, while other checks keep file names
  only (`_scrub`). The user needs the directory to fix it. **Proposed: full path.**
- Q2: under the local executor antiSMASH never runs (external step), so the warning could be moot
  there. **Proposed: warn regardless of executor**; preflight does not know the executor.
- Q3: Nextflow's `file(p).exists()` also accepts a directory named `<stem>.gbff`; the rule asks for
  a file. **Proposed: file** (a directory would fail in antiSMASH anyway).
