# M11-10: the Bioconda recipe (rendered from `pyproject.toml`)

Milestone M11 · **pass 2** · depends on: M11-01 (packaging), M11-04 (the PyPI sdist exists after the
release) · estimate ≈ 130 code lines + ≈ 110 test lines · not `stats` · roles: implementer
(`scripts/release/bioconda_recipe.py`, `packaging/bioconda/meta.yaml.j2`), test-author
(`tests/unit/test_m11_10_*.py`) · Khoa: opens the `bioconda/bioconda-recipes` PR (their account).

## Spec sections implemented

- §12 `release.yml`: "Bioconda recipe PR as a manual follow-up"; §12 Releases "PyPI and Bioconda
  packages"; §13 row M11 "packages".
- §3.2 conda channels: conda-forge and bioconda only, never `defaults`.
- R6 (the recipe's source is the PyPI sdist with its sha256).

## Definition

1. **Template** `packaging/bioconda/meta.yaml.j2` (Jinja2, already a runtime dependency): `package:
   name: giap`, `version`; `source: url: https://pypi.org/packages/source/g/giap/giap-{{ version
   }}.tar.gz`, `sha256`; `build: noarch: python`, `number: 0`, `script: {{ PYTHON }} -m pip install
   . --no-deps --no-build-isolation -vv`, `entry_points: giap = giap.cli:main`, `run_exports` none;
   `requirements: host: python >=3.11, pip, hatchling`; `run:` python ≥ 3.11, every
   `[project].dependencies` entry mapped to its conda name, plus `nextflow` and `openjdk >=17`
   (Q-L10); `test: commands: giap --version`, `giap --help`; `imports: giap, giap.cli`;
   `about: home`, `license: MIT`, `license_file: LICENSE`, `summary`, `doc_url` (site), `dev_url`;
   `extra: recipe-maintainers: [khoa-yelo]`, `identifiers: doi:<Zenodo DOI>` when known.
2. **Name map** (PyPI → conda), explicit in the script: `pyahocorasick` → `pyahocorasick`,
   `codon-bias` → `codon-bias` (bioconda), `scikit-learn` → `scikit-learn`, `pydantic>=2` →
   `pydantic >=2`, `biopython` → `biopython`, others identical. A dependency not in the map and not
   identical by rule fails the render (so a new dependency cannot be forgotten).
3. **`scripts/release/bioconda_recipe.py --version V --sha256 H [--doi D] --out DIR`** renders
   `DIR/giap/meta.yaml`. `--sha256` may be `auto`, which reads it from
   `https://pypi.org/pypi/giap/V/json` (network; used only after the release, never in tests).
4. **Hand-off**: the rendered recipe goes into the post-release PROGRESS entry and to Khoa
   (M11-05 step 16), who forks `bioconda-recipes`, adds `recipes/giap/meta.yaml` and opens the PR.
   The agent does not open PRs in other organisations' repositories.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `--version 1.0.0 --sha256 <64 hex>` on main's pyproject | `meta.yaml` with 21 runtime deps + python + nextflow + openjdk under `run:`, `noarch: python`, the URL `…/giap-1.0.0.tar.gz` |
| 2 | pyproject gains `"newlib>=1"` | render fails: `no conda name for newlib (add it to the map)` |
| 3 | `--sha256 xyz` | exit 1, `sha256 must be 64 hex characters` |
| 4 | `--doi 10.5281/zenodo.1234567` | `extra: identifiers: [doi:10.5281/zenodo.1234567]` |
| 5 | version specifiers `pydantic>=2`, `click>=8.1` | rendered as `pydantic >=2`, `click >=8.1` (conda spacing) |

## Edge cases

- **Extras** (`dev`, `docs`) are not part of the recipe.
- **`codon-bias` exists only on bioconda** (checked 2026-10-24): fine for a bioconda recipe; noted
  in the template comment.
- **Environment markers** in a dependency (`; python_version < "3.12"`): cannot become conda
  selectors in a `noarch: python` recipe, so the render fails with a message (none exist today).
- **Version pins duplicated** between pyproject and recipe: the test compares the two, so they
  cannot drift.
- Genomic edge cases: not applicable.

## Acceptance criteria

1. Rendering tests for examples 1–5 (no network; `auto` is not exercised).
2. A test that every `[project].dependencies` name appears in the rendered `run:` list with the same
   lower bound.
3. The rendered YAML parses with PyYAML and has the keys of item 1.
4. Optional, after the release: `conda build`/`bioconda-utils lint` is **not** run here (no
   installs, R14); the bioconda CI does that on Khoa's PR.

## Estimate

Code ≈ 130 (script ≈ 90, template ≈ 40). Tests ≈ 110.

## Questions for the lead (reversible defaults in use)

- Q-L10 `nextflow` + `openjdk` as run deps. Alternative: leave Nextflow to the user (smaller env,
  but `giap run` then fails until they install it).
- Maintainer handle on the recipe: `khoa-yelo` (Khoa's GitHub login).
