# giap

**GIAP** (genome-integrated analysis of proteomics data) takes one assembled bacterial genome plus
proteomics evidence and asks: *what does genome context explain about which proteins are detected,
and how abundant they are?* It places peptides and proteins on the genome with its own error
control, builds one gene table keyed by Bakta locus tags, tests genome context (regions,
autocorrelation, abundance islands, operons, function, sequence determinants), calls REF/ALT
residues and novel ORFs, and writes standard tables, plots and an HTML report.

Status: under development (milestone M0, scaffold). The build specification is
[`GIAP_SPEC.md`](GIAP_SPEC.md); progress is tracked in [`PROGRESS.md`](PROGRESS.md).

## Development on FarmShare

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
