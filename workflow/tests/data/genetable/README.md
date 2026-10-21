# Fixtures for genetable.nf.test and context.nf.test (MVP-01)

Stub runs only: `GIAP_STEP`'s stub never reads an input, so every file here is a placeholder that
only has to exist and to carry the name the real pipeline gives it. What the names pin:

* `xrefs.parquet`, `annotation_meta.json`, `tus.parquet`, `layers.parquet` — the M4-06 contract
  tables ANNOTATE emits as `contract` (C3, C4, C5, C6); the gene model (C1, C2, features) comes
  from `../mapping/`.
* `regions.parquet`, `region_genes.parquet` — C7/C8, ANNOTATE's `regions` emission.
* `mapping/` — MAPPING's loci directory (§7.2.9): `gene_table` takes the directory, while
  `context_stats` takes `peptide_hits.tsv`, `protein_hits.tsv` and `unmapped.tsv` out of it.
* `genetable/genes.parquet` — the final gene table, deliberately named like the C2 table
  `../mapping/genes.parquet`: `context_stats` stages both, which is why `GIAP_STEP` stages each
  input into its own `in<n>/` directory.
* `stratum0_classes.parquet` (ADR-0019), `peptides.norm.parquet`, `abundance.norm.parquet`,
  `preflight.json`, `cog-20.def.tab` — the optional inputs `giap run` passes when it has them.
