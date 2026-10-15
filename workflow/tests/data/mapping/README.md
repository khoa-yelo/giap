Stub-mode inputs for workflow/tests/mapping.nf.test (M3-17 example 4). Contents are never read in
stub mode: genome.fa / proteins.long.faa copy data/miniprot; the Parquet files are one-row or
empty tables with the M1-12 (`peptides.norm.parquet`) and M2-04 (`genes`, `features`,
`replicons`) columns.
contaminants.fasta (first record of tests/data/contaminants_mini.faa) and reference_peptides.tsv
(one peptide) only show the optional inputs are wired (M3-17 decision 5); never read in stub mode.
step_registry.json = the giap.steps.mapping py steps (kept equal by
tests/unit/test_mapping_nf_registry.py).
db/ is a test $GIAP_DB_DIR whose manifest lists the contaminant library (M3-15 layout), for
mapping_nf_main.nf.test.
