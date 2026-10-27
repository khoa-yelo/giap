# GIAP → GIPPA rename plan: measured inventory (ADR-0064)

Measured on origin/main 45bc3599 (2026-10-27) with `git grep` over tracked files, matching
case-insensitively unless noted. All counts are occurrences. Nothing was changed. The method,
verification, evidence and rollout are in
[ADR-0064](../adr/ADR-0064-rename-giap-to-gippa.md). This page is the input for
`scripts/dev/rename_gippa.py` and its keep list.

## 1. Headline numbers

| | Count |
|---|---|
| tracked files | 2,920 |
| files containing `giap` | 2,243 |
| occurrences of `giap` / `GIAP` / `Giap` (exact case) | 22,181 in total |
| `GIAP_SPEC` tokens (stay) | 1,171 in 1,060 files |
| historical records that stay (docs/adr, docs/specs, PROGRESS.md, docs/backlog.md) | 3,140 |
| recorded tool data and PROVENANCE files that stay | 729 |
| rename scope, excluding historical, recorded, GIAP_SPEC.md and .github/ | ≈ 15,600 lines in ≈ 1,810 files |
| rename scope in src/ and workflow/ | 2,002 lines in 326 files |
| import lines `from/import giap…` | 5,591 |
| import lines `from/import giapsim…` | 659 |
| `src/giap` path references | 713 |
| tracked paths containing `giap` | 357 |
| third-party names containing `giap` | **none found**: every hit is ours |

## 2. By area (files with a hit / total occurrences, then by exact case)

| Area | Files | Total | giap | GIAP | Giap |
|---|---|---|---|---|---|
| src/ | 295 | 2,024 | 1,392 | 380 | 252 |
| tests/ | 1,335 | 13,637 | 9,224 | 2,984 | 1,429 |
| workflow/ | 108 | 510 | 280 | 230 | 0 |
| scripts/ | 88 | 1,147 | 568 | 577 | 2 |
| docs/ | 312 | 3,737 | 2,460 | 1,225 | 52 |
| benchmarks/ | 36 | 129 | 113 | 12 | 4 |
| validation/ | 25 | 539 | 269 | 253 | 17 |
| changelog.d/ (47 fragments) | 20 | 42 | 29 | 13 | 0 |
| .github/ (**Khoa**) | 7 | 41 | 29 | 12 | 0 |
| .claude/ | 3 | 4 | 0 | 4 (all `GIAP_SPEC`) | 0 |
| ci/ | 2 | 5 | 3 | 2 | 0 |
| root files | — | GIAP_SPEC.md 228 (stays), PROGRESS.md 50, README.md 22, pyproject.toml 18, AGENTS.md 12, NOTICE.md 12, CITATION.cff 4, CLAUDE.md 2, mkdocs.yml 1, .gitignore 1, uv.lock 1 | | | |

Sub-areas:

| Sub-area | Files | Occurrences | of which `GIAP_SPEC` |
|---|---|---|---|
| src/giap/schemas | 45 | 69 | 43 |
| src/giap/config | 2 | 14 | 2 |
| src/giap/data | 8 | 21 | 8 |
| workflow/conf | 6 | 30 | 5 |
| workflow/envs | 15 | 20 | 2 |
| workflow/containers | 8 | 51 | 6 |
| scripts/slurm | 42 | 797 | 17 |
| docs/user-guide | 11 | 120 | 3 |
| docs/methods | 32 | 249 | 33 |
| docs/specs (stays) | 214 | 2,805 | 38 |
| docs/adr (stays) | 42 | 153 | 9 |
| validation/pilot | 18 | 452 | 8 |
| tests/data/recorded (stays) | 58 | 662 | |
| tests/data/snapshot | 9 | 17 | |
| tests/golden | 8 | 80 | |

The most frequent tokens across the repository:

| Token | Count |
|---|---|
| `giap` | 11,406 |
| `giapsim` | 1,961 |
| `GiapConfig` | 1,740 |
| `GIAP` | 1,350 |
| `GIAP_SPEC` | 1,171 |
| `GIAP_ROOT` | 1,105 |
| `GIAP_DB_DIR` | 446 |
| `GIAP_REPO` | 355 |
| `GIAP_STEP` | 243 |
| `GIAP_DATA_DIR` | 163 |
| `giap_fgsea` | 119 |
| `giap_run` | 110 |
| `GIAP_PYTHON` | 89 |
| `giap_version` | 77 |
| `GIAPORF_` | 66 |
| `.giap` | 66 |

## 3. Paths to rename (`git mv`)

- `src/giap/` → `src/gippa/` (317 files).
- `tests/sim/giapsim/` → `tests/sim/gippasim/` (26 files).
- Single paths:
  - `benchmarks/test_giapsim.py`
  - `tests/unit/test_giap_run_sbatch_pins_code.py`
  - `tests/unit/test_giapsim_selenoproteins.py` and `tests/unit/test_giapsim_selenoproteins_exclusions.py`
  - `scripts/slurm/giap_run.sbatch` and `scripts/slurm/giap_setup.sbatch`
  - `workflow/bin/giap_fgsea.R`
  - `workflow/modules/local/giap_step/` (`main.nf`, `meta.yml`)
  - `workflow/tests/bin_giap/giap`
  - `workflow/tests/fake_giap.config`
- Stay: `docs/adr/ADR-0064-rename-giap-to-gippa.md` (the rename record). Since Q-K32
  (ADR-0064 amendment) `GIAP_SPEC.md` → `GIPPA_SPEC.md` and `docs/specs/M2-08-giap-setup-core.md`,
  `docs/specs/M2-09-giap-setup-resources.md` move too.
- In total, 357 tracked paths contain `giap`. 343 of them are in the two trees above.

## 4. Identifiers that users see

| Kind | Names |
|---|---|
| CLI | `giap` (pyproject `[project.scripts]`); `prog_name="giap"` (cli.py ×3); `giap version --json` key `"giap"`; the internal `giap _step` |
| distribution | pyproject `name = "giap"`; hatch `only-include src/giap`, `workflow` → `giap/workflow`; towncrier `package`/`name`/`title_format`/`issue_format` URL; `packaging/bioconda/meta.yaml.j2` |
| env vars that src/ and workflow/ read | `GIAP_ROOT`, `GIAP_DB_DIR`, `GIAP_REPO`, `GIAP_TOOLS_DIR`, `GIAP_FGSEA_PREFIX`, `GIAP_CONTAINER_RUNTIME` (set by `giap run`), `GIAP_SMOKE_DB_DIR`, `GIAP_HARNESS_REPO`, `GIAP_HARNESS_BAKTA_JSON`; Nextflow `System.getenv('GIAP_ROOT'|'GIAP_DB_DIR'|'GIAP_SMOKE_DB_DIR')` |
| other env vars | 63 distinct `GIAP_*` names in all (excluding `GIAP_SPEC`), mostly calibration, nightly and test harness, e.g. `GIAP_CALIBRATION_*`, `GIAP_SUITE`, `GIAP_WORKTREE`, `GIAP_NOW`, `GIAP_PYTHON`, `GIAP_UPDATE_GOLDEN` |
| config | no config key contains `giap`; only the schema title `GiapConfig` and the description text change |
| output fields | `giap_version` in run_manifest, preflight, makedb_manifest and the VCF writer; manifest tool row `name: "giap"`; datapackage `name: giap-run`; schema `$id` URLs (db_manifest, preflight, run_manifest) |
| output content | GFF3 column 2: `giap` (gff3.py region line, user_annotation) and `GIAP` (mapping/outputs.py, novel/outputs.py, uniop_adapter); VCF `##source=giap-<v>`; ORF IDs `GIAPORF_{class}_…` (makedb, novel/fdr.py parser, genetable/detection) |
| report text | `report.html.j2` (5), `compare.html.j2` (4), `tables.js` (1), `report/build.py` "not available in this GIAP version", `content.py` "GIAP version", `preflight/report.py` "GIAP <v>" |
| Nextflow | processes `GIAP_STEP`, `GIAP_MAPPING`, `GIAP_MAPPING_LOCI` (appear in trace and manifest process lists); staging file `.giap_inputs` |
| on-disk names | `.giap/capabilities.json` (probe_env; `.gitignore`); env marker `.giap-built`; `share/<tool>/GIAP_COMMIT` in envs and images; envs `giap-<tool>` (15) and `giap-wf` on `PATH`; `/tmp/giap-numba-$USER` |
| images | `ghcr.io/khoa-yelo/giap-<tool>` in update_digests.py, docs and tests; Dockerfile labels `org.giap.image.key`, `image.title giap-<tool>`; `/opt/giap/giap_fgsea.R` |
| Slurm | 44 distinct `giap-*` job names; `squeue | grep ' giap-'` in resume_status.sh and session_start.sh; scrontab `giap-nightly-launch` |
| URLs | `github.com/khoa-yelo/giap`: CITATION.cff, pyproject, cli.py, 3 schemas, user-guide/run.md, bioconda recipe, 3 tests |

## 5. Where a mechanical sed would be wrong

1. **`GIAP_SPEC`** (1,171). The file name stays, so this token goes on the keep list.
2. **Protein sequences.** GIAP is a valid tetrapeptide.
   - It occurs 33–38 times in each recorded Bakta proteome, in `.faa.gz`, `.gbff.gz` and `.json.gz`
     for K-12 and Nissle, and 34 times in `validation/pilot/db_a2/strain_db.faa.gz`.
   - Simulated and real runs can produce it at any time.
   - Never touch `.gz` or other binary files. The comparator compares sequences without
     normalisation (ADR-0064 item 8).
3. **Recorded tool data pinned by sha256.**
   - `tests/data/recorded/uniop/3ed14dfc/{k12,sim100k}/uniop_cds.gff` contain 300 `GIAP`
     source-column values.
   - Their sha256 is pinned in PROVENANCE.md, in `validation/calibrate_cutoffs.py`
     `K12_PINS.uniop_gff_sha256` (ADR-0022 cut-offs) and in
     `tests/unit/test_m6_12_calibrate_cutoffs.py`.
   - Keep all 58 files under tests/data/recorded and all 62 PROVENANCE files as recorded.
   - After the rename the adapter writes `GIPPA`. Any test that diffs fresh adapter output
     against a recording must normalise column 2. Today, `k12_cutoffs_world` regenerates its GFF
     with the adapter and re-pins it, so it is safe.
4. **Escaped strings.** `"\tgiap"`, `"\ngiap"` and `b"giap"` appear in gff3.py, test_gff3.py,
   makedb_genomes.py, novel_outputs_fixture.py and nightly.sbatch. A `\bgiap\b` regex misses
   them. Use substring replacement.
5. **Compound identifiers.** `NotAGiapRun`, `tGIAP`, `GIAPORF_`, `rev_GIAPORF_`, `GIAPNOVEL_`
   (tests) and `GIAPSIM` (a test constant). The case-preserving substring map handles them, but
   they must be checked in `--check`.
6. **The SVG hash salt `"giap"`** (`io/plot_files.py` HASHSALT, `viz/style.py`). Keep it
   (ADR-0064 item 7). Renaming it changes every SVG id.
7. **Snapshot fixtures.** These are `tests/data/snapshot/*/run_manifest.json` and
   `preflight.json` (key `giap_version`, tool row `giap`). Rename the keys so the reader still
   finds them. Their sha256 values are of input files and do not change.
8. **Sim outputs that carry the simulator's name**: GFF source `giapsim`, `source` column
   `giapsim`, `data-version: giapsim/synthetic`, `giapsim_git_sha`. Tests that pin sim-output
   sha256 values (about 30 files with 64-hex literals) may need re-pins. List each with its
   normalised-equal proof.
9. **Image keys.** They hash Dockerfiles, env ymls, locks and `giap_fgsea.R` by path. All keys
   change, which is expected (`images.json` is `{}`). Tests that pin a key value must be
   recomputed, not hand-edited.
10. **Env ymls and locks.** Changing `name: giap-<tool>` changes the sha256 behind
    `.giap-built`, so every prebuilt env reads as stale until rebuilt. This is intended (ADR-0064
    item 5). Do not hand-edit the markers.
11. **URLs.** `github.com/khoa-yelo/giap` redirects after the repo rename, but GitHub Pages does
    not. Rename the URLs in the same PR. Historical records keep the old URL.
12. **Historical prose.** Superseded by Q-K32 (ADR-0064 amendment): every document is renamed,
    including docs/adr, docs/specs, PROGRESS.md, docs/backlog.md and validation reports. Only
    the rename record (ADR-0064, its index line, this page) and PROVENANCE text stay as written.
13. **`scripts/release/readiness.py`.** NON_STATS_PATHS and NEVER_PATHS list `src/giap/...`, and
    `VERSION_FILE` is `src/giap/__init__.py`. The edit is mechanical but needs `stats` review
    (Q-L15), and the readiness unit tests (disjointness, examples 13c) follow.
14. **Paths outside git** (no sed): `$GIAP_ROOT/micromamba/envs/giap-*`, Nextflow `work/` and
    `-resume` caches, local Apptainer images, `local/logs/giap-*`, the installed scrontab and
    the shared venv's editable `giap` (`uv sync` in the main checkout after the merge).

## 6. `.github/` changes for Khoa (agents never edit .github/)

The `.github/` part has 41 occurrences in 7 files:

- `ci.yml`:
  - `mypy src/giap`;
  - the installed-wheel smoke (`giap --version`, `giap version --json`, `giap preflight`);
  - `site-packages/giap/workflow`;
  - `workflow/bin/giap_fgsea.R`;
  - `GIAP_DB_DIR: /tmp/giap-db-none`.
- `containers.yml`:
  - the `workflow/bin/giap_fgsea.R` path filter and the `changed=` patterns;
  - `from giap.annotate.tools import TOOLS`;
  - the label `org.giap.image.key`.
  - The prefix is already derived from pyproject `[project] name` (M11-14 §6).
- `docker-smoke.yml`: `GIAP_SMOKE_DB_DIR` (×4).
- `docs.yml`: path filters `src/giap/cli.py`, `src/giap/config/**`, `src/giap/schemas/**`.
- `golden-r.yml`: env `giap-golden-r`.
- `release.yml`:
  - `import giap; giap.__version__`;
  - `/tmp/pkg/bin/giap`;
  - the PyPI and TestPyPI JSON URLs `…/pypi/giap/$VERSION/json`;
  - the release title `GIAP $VERSION`.
- `ISSUE_TEMPLATE/feature.md`: `GIAP_SPEC.md` (stays).

The repository and publishing settings are Khoa's alone:

- the repository rename;
- PyPI and TestPyPI trusted publishers for project `gippa`, repo `khoa-yelo/gippa`, workflow
  `release.yml`;
- the `ghcr` environment and `GHCR_PUBLISH` (M11-14 14j);
- reinstalling the scrontab.

## 7. Script outline (for the rename PR, not this one)

```
rename_gippa.py [--apply | --check]
  keep paths (data only since Q-K32): ADR-0064 and this page (the rename record),
               tests/data/recorded/**, tests/data/uniop/**, **/PROVENANCE.md,
               validation/pilot/{REPORT.md,metrics/**,db_a2/**,inputs/**,diann_A*.cfg},
               the rename tooling, .github/** (Khoa)
  keep tokens: the HASHSALT/svg.hashsalt "giap" literals, .giap-complete, .giap_setup.lock,
               the checkout path, the ADR-0064 file name and index line
  rewrite:     every document; GIAP_SPEC.md → GIPPA_SPEC.md; the old expansion →
               "genome integrated prokaryotic proteomics analysis"
  skip:        binary files and *.gz
  map:         giap→gippa, GIAP→GIPPA, Giap→Gippa (substring, case-exact), then git mv paths
  regenerate:  config.schema.json, docs reference pages (gen_reference.py), uv.lock
  --check:     prints the remaining hits by keep category; exit 1 on any uncategorised hit
               or any giap left in the docs outside the rename record and kept data tokens
```
