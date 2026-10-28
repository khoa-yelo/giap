# M11-01: installable packages (wheel and sdist that carry the workflow)

Milestone M11 · **pass 1 (MVP)** · depends on: none · estimate ≈ 190 code lines + ≈ 230 test lines ·
not `stats` · roles: implementer (`pyproject.toml`, `src/gippa/__init__.py`, `scripts/release/check_dist.py`,
the `package` job in `.github/workflows/ci.yml`), test-author (`tests/unit/test_m11_01_*.py`) · Khoa: none.

## Spec sections implemented

- §13 row M11 "packages"; §12 `release.yml` "build sdist and wheel"; §12 **Releases** "semantic
  versioning … PyPI and Bioconda packages".
- §4: `pyproject.toml` "src layout, `gippa` console script, extras: dev, docs, viz"; `workflow/` is
  part of the product (`gippa run` launches `workflow/main.nf`, §6.3).
- R5 (no GPL/AGPL/unlicensed code shipped), R11 (no file > 5 MB; never redistribute RegulonDB,
  raw MS, DIA-NN, Bakta/geNOMAD DBs, UniOP code), R12 (no secrets).

## The gap

`[tool.hatch.build.targets.wheel] packages = ["src/gippa"]` ships only the Python package.
`gippa.runner.nextflow.workflow_main()` already looks for `gippa/workflow/main.nf` inside the installed
package, but nothing puts it there, so `pip install gippa && gippa run …` cannot find a workflow.
The version is written by hand in three places (`pyproject.toml`, `src/gippa/__init__.py`,
`CITATION.cff`).

## Definition

1. **Single-source version.** `pyproject.toml` declares `dynamic = ["version"]` and
   `[tool.hatch.version] path = "src/gippa/__init__.py"`. `__version__` stays the only literal.
   (`CITATION.cff` is checked against it in M11-02, not generated.)
2. **Workflow in the wheel.** `workflow/` is force-included as `gippa/workflow/` **except**
   `workflow/tests/`, `.nf-test*`, `.nextflow*`, `work/`. File modes are kept: every file under
   `gippa/workflow/bin/` and `gippa/workflow/assets/ps/ps` is executable after `pip install`
   (Nextflow puts `bin/` on the task PATH and needs `+x`).
3. **sdist content** (Q-L1 default): `src/`, `workflow/` (same exclusions), `pyproject.toml`,
   `README.md`, `LICENSE`, `NOTICE.md`, `CHANGELOG.md`, `CITATION.cff`. Nothing else.
4. **`scripts/release/check_dist.py DIST_DIR [--version V]`** opens every `*.whl` and `*.tar.gz` in
   DIST_DIR (stdlib `zipfile`/`tarfile`; never installs) and checks:
   - *required* (wheel): `gippa/__init__.py`, `gippa/cli.py`, `gippa/py.typed`,
     `gippa/config/defaults.yaml`, `gippa/schemas/datapackage.json`, `gippa/data/calibration.yaml`,
     `gippa/data/db_sources.yaml`, `gippa/workflow/main.nf`, `gippa/workflow/nextflow.config`, every
     `workflow/envs/*.yml` and `workflow/conf/*.config` of the checkout, the `gippa` console entry
     point in `entry_points.txt`, `License-Expression: MIT` (or `License: MIT`) in METADATA, and
     `LICENSE` in the dist-info licence files;
   - *required* (sdist): the item-3 files;
   - *forbidden* (both): any path under `local/`, `tests/`, `validation/`, `.github/`, `.claude/`,
     `docs/`, `benchmarks/`, `workflow/tests/`; suffixes `.sif`, `.raw`, `.d`, `.mzML`, `.parquet`;
     `__pycache__`/`.pyc`; any file whose basename matches `*regulondb*` (case-insensitive) or
     lives under a `uniop/` directory other than `gippa/workflow/modules/local/uniop/` and
     `gippa/workflow/containers/uniop/`; any file > 5 MB (R11);
   - *version*: with `--version V`, both file names carry V and METADATA `Version: V`; without it,
     V = `gippa.__version__` of the checkout;
   - exit 0 and one `ok …` line per artefact; exit 1 and one `FAIL <artefact>: <rule>: <path>` line
     per violation (all violations, not the first); exit 2 when DIST_DIR has no wheel or no sdist.
5. **CI `package` job** (ubuntu, Python 3.12, `needs: lint`): `uv build` → `check_dist.py dist`
   → a **clean** venv outside the checkout (`uv venv /tmp/pkg && uv pip install dist/*.whl`) →
   from `/tmp` (not the checkout): `gippa --version`; `gippa version --json`; `gippa preflight` on
   `tests/data/genomes/` tiny genome (path passed in); `python -c` asserting
   `workflow_main()` is under the venv's `site-packages/gippa/workflow/` and exists;
   with Nextflow installed as in `nf-test-stub`, `nextflow config <that dir> -profile test,stub` and
   `-profile farmshare,apptainer` both exit 0 (Q-L2).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `uv build` on this branch | `dist/gippa-0.0.1-py3-none-any.whl`, `dist/gippa-0.0.1.tar.gz`; `check_dist.py dist` exits 0 with two `ok` lines |
| 2 | the same wheel, `unzip -l` | contains `gippa/workflow/main.nf`, `gippa/workflow/modules/local/bakta/main.nf`, `gippa/workflow/envs/bakta.yml`; no `gippa/workflow/tests/` |
| 3 | wheel installed in a clean venv, `ls -l site-packages/gippa/workflow/bin/gippa_fgsea.R` | mode has `x` for the owner |
| 4 | a wheel with an extra `gippa/data/huge.bin` of 6 MB (fixture made in `tmp_path`) | exit 1, `FAIL gippa-…whl: file > 5 MB: gippa/data/huge.bin` |
| 5 | an sdist that contains `tests/data/x.tsv.gz` | exit 1, `FAIL gippa-….tar.gz: forbidden path: gippa-0.0.1/tests/data/x.tsv.gz` |
| 6 | `check_dist.py dist --version 1.0.0` on 0.0.1 artefacts | exit 1, one `version` FAIL per artefact |
| 7 | empty DIST_DIR | exit 2, `no wheel and no sdist in <dir>` |
| 8 | `__version__ = "1.0.0"` edited in `src/gippa/__init__.py` only, then `uv build` | artefacts named `1.0.0` (the pyproject has no literal version) |

## Edge cases

- **Checkout vs installed lookup.** In the checkout `workflow_main()` keeps preferring
  `<repo>/workflow/main.nf` (editable installs, nightly). The smoke runs from `/tmp`, so the
  checkout candidate (`parents[3]` = `…/lib/python3.12`) does not exist and the package path wins.
  A test fakes both layouts in `tmp_path`.
- **Read-only site-packages.** Nothing in the workflow writes into `projectDir` (`.nextflow/` and
  `work/` go to the launch and work dirs). A test greps `workflow/` for writes to `projectDir`
  other than the `assets/ps/ps` read-only bind.
- **Symlinks** under `workflow/` (none today) are rejected by `check_dist.py` (`symlink`), since
  wheels cannot carry them portably.
- **Generated files** in the tree (`__pycache__`, `.nextflow.log`) are excluded through `.gitignore`
  (hatch honours it) and caught by the forbidden list if they slip in.
- **Coordinates, strands, circular/linear genomes:** not applicable (no genomic logic).
- **Empty input:** see example 7.

## Acceptance criteria

1. `pyproject.toml` has no literal `version =` in `[project]`; `uv build` yields artefacts named
   with `gippa.__version__` (example 8, run in CI only: building needs hatchling from PyPI).
2. `check_dist.py` implements every rule of item 4; unit tests build small fake wheels/sdists in
   `tmp_path` (no network, no hatchling) and cover examples 4–7, the symlink rule and the
   "all violations reported" rule (two violations → two FAIL lines).
3. CI `package` job passes on the PR and runs the four smoke commands from a directory outside the
   checkout; the job is listed under the `summary` job's `needs`.
4. A unit test asserts the hatch config excludes `workflow/tests` and includes `workflow/bin` (reads
   `pyproject.toml`, no build).
5. `scripts/ci/check_imports.py src` still passes; no new runtime dependency.

## Test plan

`tests/unit/test_m11_01_check_dist.py` (fake artefacts), `tests/unit/test_m11_01_pyproject.py`
(static config), `tests/unit/test_m11_01_workflow_main.py` (both layouts). The real build runs in
the CI job only.

## Estimate

Code ≈ 190: `check_dist.py` ≈ 140, pyproject ≈ 15, ci.yml job ≈ 35. Tests ≈ 230.

## Questions for the lead (reversible defaults in use)

- Q-L1 sdist contents (index). Default: no tests in the sdist. Bioconda (M11-10) tests with
  `gippa --version` and imports only.
- Q-L2 installed smoke without `gippa run` (the CLI rejects `test`/`stub` profiles). Default:
  `nextflow config` parse + `workflow_main()` path; no CLI change.
