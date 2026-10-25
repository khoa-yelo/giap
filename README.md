# giap

**GIAP** (genome-integrated analysis of proteomics data) takes one assembled bacterial genome plus
proteomics evidence and asks: *what does genome context explain about which proteins are detected,
and how abundant they are?* It places peptides and proteins on the genome with its own error
control, builds one gene table keyed by Bakta locus tags, tests genome context (regions,
autocorrelation, abundance islands, operons, function, sequence determinants), calls REF/ALT
residues and novel ORFs, and writes standard tables, plots and an HTML report.

[![PyPI version](https://img.shields.io/pypi/v/giap.svg)](https://pypi.org/project/giap/)

Documentation: [khoa-yelo.github.io/giap](https://khoa-yelo.github.io/giap/). Changes between releases are listed in
[`CHANGELOG.md`](CHANGELOG.md).

## Install

```bash
pip install giap                                     # Python package and the `giap` command
micromamba install -c conda-forge -c bioconda giap   # or from Bioconda
```

GIAP runs its steps through Nextflow; see the documentation for the container and Conda profiles.

## Quickstart

```bash
giap setup --db-dir db                                   # databases and models, once
giap preflight --genome genome.fna --peptides peptides.tsv   # check the inputs in seconds
giap run --genome genome.fna --peptides peptides.tsv --db-dir db --outdir out --profile docker
open out/report/report.html
```

Inputs, options and outputs are described in the
[documentation](https://khoa-yelo.github.io/giap/).

## Citing GIAP

If you use GIAP, please cite it as described in [`CITATION.cff`](CITATION.cff) (GitHub's
"Cite this repository" button reads the same file).

## Developers

The build specification is [`GIAP_SPEC.md`](GIAP_SPEC.md); progress is tracked in
[`PROGRESS.md`](PROGRESS.md). Every PR that changes `src/` or `workflow/` adds a fragment under
[`changelog.d/`](changelog.d/) (checked in CI by `scripts/ci/check_changelog.py`).

### Development on FarmShare

```bash
source scripts/farmshare_env.sh          # all data and environments live in ./local (gitignored)
bash scripts/bootstrap_farmshare.sh      # once: uv, Nextflow, nf-test, gh
uv sync --frozen --all-extras
uv run pytest -m "not slow and not calibration_full and not e2e" -n auto
```

See [`AGENTS.md`](AGENTS.md) for test layers and conventions.

## Licence

MIT. Third-party tools keep their own licences and run only as separate programs; see
[`NOTICE.md`](NOTICE.md).
