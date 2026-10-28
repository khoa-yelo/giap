# Fixtures for operons.nf.test (M6-08)

Stub runs only: `GIPPA_STEP`'s stub never reads an input, so each file here is a placeholder. It
only has to exist and carry the name the real pipeline gives it.

* `pairs_raw.parquet`: the M2-12 `operon_pairs` output. It is the only required member of
  ANNOTATE's `layers_for_pairs` emission.
* `promoter_windows.parquet`, `terminators.parquet`, `tx_rates.parquet`: the three optional
  layers that `annotate.nf` emits for the pair evidence (M6-08 lead decision 1).

The gene table, the gene model and the contract tables come from `../genetable/` and
`../mapping/`, as they do in context.nf.test.
