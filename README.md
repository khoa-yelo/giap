# GIPPA

**GIPPA** (genome integrated prokaryotic proteomics analysis) takes one assembled bacterial genome
plus proteomics evidence and asks: *what does genome context explain about which proteins are
detected, and how abundant they are?* It places peptides and proteins on the genome with its own
error control, builds one gene table keyed by Bakta locus tags, tests genome context (regions,
autocorrelation, abundance islands, operons, function, sequence determinants), calls REF/ALT
residues and novel ORFs, and writes standard tables, plots and an HTML report.

[![PyPI version](https://img.shields.io/pypi/v/giap.svg)](https://pypi.org/project/giap/)

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
- the final pilot report, and the code rename: until it lands, the package, the command and the
  environment variables are `giap` (`giap run`, `GIAP_DB_DIR`);
- genetic code 4 (Mycoplasma and relatives) is planned after 1.0; only code 11 is supported now.

Off FarmShare, `giap preflight`, `convert`, `compare` and `report` work anywhere; until the first
image publish, `giap run` under `--profile docker` or `apptainer` alone stops in preflight for the
four tools above ([Install](docs/get-started/install.md)).

## Install

From the first release:

```bash
pip install giap                                     # Python package and the `giap` command
micromamba install -c conda-forge -c bioconda giap   # or from Bioconda
```

Until then, from a checkout (Python >= 3.11):

```bash
git clone https://github.com/khoa-yelo/gippa.git && cd gippa
python3 -m venv .venv && . .venv/bin/activate && pip install .
```

`giap run` also needs Nextflow (with Java 17+) and Docker or Apptainer; see
[Install](docs/get-started/install.md) for the prerequisites, profiles and the FarmShare setup.

## Quickstart

```bash
giap setup --db-dir db --antismash --go --kegg --cog --contaminants --images --profile apptainer
GIAP_DB_DIR=db giap preflight --genome tests/data/genomes/sim100k.fna.gz --peptides docs/tutorials/data/sim100k_peptides.tsv.gz --profile apptainer
giap run --genome tests/data/genomes/sim100k.fna.gz --peptides docs/tutorials/data/sim100k_peptides.tsv.gz --db-dir db --outdir out --profile apptainer
xdg-open out/report/report.html
```

Run it from a repository checkout (`--profile docker` works the same way): the genome and the
peptide table are a small synthetic demo made by giapsim, the test simulator
([provenance](docs/tutorials/data/PROVENANCE.md)). Today the full run works on FarmShare, where
`--profile` defaults to `farmshare,apptainer`. The walk-through is
[Quickstart](https://khoa-yelo.github.io/gippa/get-started/quickstart/)
([source](docs/get-started/quickstart.md)). For your own data, convert the search-engine output
with [`giap convert`](docs/user-guide/convert.md); inputs, options and outputs are in
[Inputs](docs/user-guide/inputs.md), the [CLI reference](docs/reference/cli.md) and
[Outputs](docs/user-guide/outputs.md).

## Citing GIPPA

If you use GIPPA, please cite it as described in [`CITATION.cff`](CITATION.cff) (GitHub's
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

geNOMAD (module `mobile`) is for non-commercial use only (`giap setup --no-genomad`,
`--modules core,bgc,function`); KEGG (`--kegg`) is free for academic use only; see
[Licences](docs/user-guide/licences.md).
