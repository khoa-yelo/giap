# M9-25 snapshot fixtures (`tests/data/snapshot/`)

Frozen inputs of the plot snapshot tests (`tests/snapshot/`, docs/specs/M9-25-snapshot-baselines.md
item 1). Synthetic only: built from giapsim (`tests/sim/giapsim`), no real genome, peptide or
database content. The function database is the synthetic one of `tests/helpers/function_db.py`
(made-up GO/KEGG/COG IDs); `$GIAP_DB_DIR` points at it during the build, never at the real one.

## How they are built

```
source scripts/farmshare_env.sh
PYTHONPATH=$PWD/src:$PWD/tests/sim:$PWD/tests MPLBACKEND=Agg \
  $UV_PROJECT_ENVIRONMENT/bin/python scripts/dev/make_snapshot_fixture.py
```

About 6 minutes on one CPU. `scripts/dev/make_snapshot_fixture.py`:

1. simulates giapsim `SIM100K`'s signals (plasmid, rRNA/tRNA, CRISPR, oriC, a 12-kb prophage, a
   10-kb BGC, pseudogenes, sORFs, promoter/RBS/terminator layers) at 300 kb, seed 20261014, with a
   SNV-only strain (1 SNV/kb), de novo peptides with their paired reference peptides, 40 planted
   novel ORFs and 3 N-terminal extensions (the M7 chain settings of `tests/helpers/m7_chain.py`);
2. adds synthetic COG/GO/KEGG `db_xrefs` to the reference Bakta JSON;
3. plants, in the peptide inputs: an abundance island (26 consecutive chr1 genes from CDS 120,
   intensities x 2^12), low detection in the prophage and the BGC (every second gene's spectra
   dropped); and in the terminator layer: four confidence-1 calls inside operons (flagged pairs);
4. runs `giap run --executor local --modules core,mobile,bgc,function --genome-relation related`
   twice (`run_a` seed 7; `run_b` seed 11 with every third spectrum dropped) with
   `report.max_locus_plots=2`, `context.islands.detrend_span_frac=0.5` (the planted island is
   wider than the default 10 % span on a 300-kb genome), `determinants.permutations=19`,
   `stats.bootstrap=199`, then `giap compare run_a run_b`;
5. keeps the `OUT/` files `RunView` reads (minus the VCF, pileup, FASTA/GFF3 and every Parquet
   file: no plot reads them, and R11 allows no committed Parquet here), `logs/steps` and
   `plots/support`, writes `resolved_config.yaml` (defaults +
   the settings above; no input paths), replaces the scratch path by `WORK` and time stamps by
   `2026-01-01T00:00:00+00:00` and step timings and free disk by 0 (so a rebuild gives the same
   bytes), and gzips TSVs over 16 kB (in place for schema tables, which are
   read by magic bytes; as `.tsv.gz` for support files);
6. copies that tree to `sim_empty/` with every TSV cut to its header (support files included).

| directory | content |
|---|---|
| `sim_plots/OUT/` | `run_a`'s `OUT/` subset (every plot non-empty) |
| `sim_plots/CMP/` | `compare_{detection,abundance,stability}.tsv`; `runs/{run_a,run_b}/` with `tables/genes.tsv` and `preflight/preflight.json` |
| `sim_empty/` | the same layout, header-only (empty states) |

Regenerate only in a PR that says why; then regenerate the baselines and the structure snapshot
(below) and show before/after images in the PR (GIAP_SPEC.md §10.1).

## Baselines and structure snapshot

From the repository root (with `PYTHONPATH=$PWD/src:$PWD/tests/sim:$PWD/tests`, the shared venv):

```
pytest tests/snapshot/test_plots.py -k test_plot_snapshot --mpl-generate-path=tests/snapshot/baseline
(cd tests && python -m helpers.snapshot_fixture --write-structure)
```

The first regenerates every baseline (rerun it after any plot style change, before a packet);
the second rewrites `tests/snapshot/report_structure.json`. Baselines are PNG at 100 dpi, compared at RMS <= 2 by the `snapshot` CI job (ubuntu, Python
3.12) only. Khoa approves them through `scripts/dev/approval_packet.py`.
