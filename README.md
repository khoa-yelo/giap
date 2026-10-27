# GIPPA

**GIPPA** (genome integrated prokaryotic proteomics analysis) takes one assembled bacterial genome
plus proteomics evidence and asks: *what does genome context explain about which proteins are
detected, and how abundant they are?* It places peptides and proteins on the genome with its own
error control, builds one gene table keyed by Bakta locus tags, tests genome context (regions,
autocorrelation, abundance islands, operons, function, sequence determinants), calls REF/ALT
residues and novel ORFs, and writes standard tables, plots and an HTML report.

[![PyPI version](https://img.shields.io/pypi/v/gippa.svg)](https://pypi.org/project/gippa/)

Documentation: [khoa-yelo.github.io/gippa](https://khoa-yelo.github.io/gippa/) (from the first
release; until then the same pages are under [`docs/`](docs/index.md)). Changes between releases
are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Status

Pre-1.0 (`0.0.1`); nothing is published yet. The pipeline runs end to end on Stanford FarmShare
from a checkout (Slurm, Apptainer and prebuilt environments), on the synthetic demo below and on
the E. coli Nissle 1917 pilot (PXD078459; [`validation/pilot/README.md`](validation/pilot/README.md)).
Pending:

- the first release: PyPI, the GHCR images of four tools (Promotech, Promoter Calculator,
  TerminatorNet, OSTIR), the hosted docs and Bioconda;
- the final pilot report;
- genetic code 4 (Mycoplasma and relatives) is planned after 1.0; only code 11 is supported now.

Off FarmShare, `gippa preflight`, `convert`, `compare` and `report` work anywhere; until the first
image publish, `gippa run` under `--profile docker` or `apptainer` alone stops in preflight for the
four tools above ([Install](docs/get-started/install.md)).

## Install

From the first release:

```bash
pip install gippa                                     # Python package and the `gippa` command
micromamba install -c conda-forge -c bioconda gippa   # or from Bioconda
```

Until then, from a checkout (Python >= 3.11):

```bash
git clone https://github.com/khoa-yelo/gippa.git && cd gippa
python3 -m venv .venv && . .venv/bin/activate && pip install .
```

`gippa run` also needs Nextflow (with Java 17+) and Docker or Apptainer; see
[Install](docs/get-started/install.md) for the prerequisites, profiles and the FarmShare setup.

## Quickstart

```bash
pip install gippa                         # or: conda install -c bioconda gippa
gippa demo                                # about 2 minutes; prints the report path
gippa setup                               # databases and images for the detected runtime
gippa run --genome g.fna --peptides p.tsv --outdir out
```

`gippa demo` runs a small bundled synthetic genome with no Nextflow, containers, databases or
network. `gippa setup` and `gippa run` detect the Nextflow profile (pass `--profile` to choose).
Until the first release, install from a checkout as above; the full run then works on FarmShare.
The guided walk-through, with example inputs and outputs, is
[Quick start](https://khoa-yelo.github.io/gippa/get-started/quickstart/)
([source](docs/get-started/quickstart.md)). With raw spectra and no peptide table yet, see
[From raw spectra to a GIPPA run](docs/get-started/from-spectra.md); for converting search-engine
output, [`gippa convert`](docs/user-guide/convert.md).

## Citing GIPPA

If you use GIPPA, please cite it as described in [`CITATION.cff`](CITATION.cff) (GitHub's
"Cite this repository" button reads the same file).

## Developers

The build specification is [`GIPPA_SPEC.md`](GIPPA_SPEC.md); progress is tracked in
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

geNOMAD (module `mobile`) is for non-commercial use only (`gippa setup --no-genomad`,
`--modules core,bgc,function`); KEGG (`--kegg`) is free for academic use only; see
[Licences](docs/user-guide/licences.md).
