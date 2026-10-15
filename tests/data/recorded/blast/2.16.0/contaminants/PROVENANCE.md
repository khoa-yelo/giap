# BLAST+ 2.16.0 blastp recorded run: synthetic queries vs the contaminant fixture library (M3-15)

- Tool: BLAST+ 2.16.0 (`blastp: 2.16.0+`, `Package: blast 2.16.0, build Mar 28 2025 16:30:45`; NCBI, public domain)
- Image (pinned in workflow/conf/containers.config, selectors `BLAST` and `TBLASTN`):
  `docker://quay.io/biocontainers/blast@sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c`
  (quay.io/biocontainers/blast:2.16.0--h66d330f_5), run with Apptainer 1.4.4 from `$GIAP_ROOT/containers/`
- Date: 2026-10-14
- Host: iron-01 (Stanford FarmShare, Slurm allocation), 2 threads requested (BLAST ignores
  `-num_threads` with `-subject` and warns on stderr)

## Inputs
- Library: `tests/data/contaminants_mini.faa` (trypsin P00761, KRT1 P04264, KRT10 P13645,
  albumin P02769), copied as `contaminants_mini.faa`;
  sha256 `eaec55339063a909598888997db9ad78a1ac61400cebbf2a83537c3e7526f586`
- Queries: `proteins.faa`, written by `make_proteins.py` (this directory, seed 15):
  `PYTHONPATH=src python make_proteins.py tests/data/contaminants_mini.faa proteins.faa`;
  per library entry two `HIT_<ACC>_<k>` (substitutions to 93-99 % identity, full length) and one
  `LOWCOV_<ACC>_1` (first 40 % of the entry + 60 % random residues), plus `RAND_1..4`;
  sha256 `0477c64ce7d66fc53c5945f7190833f7e414b2be40ec23e1fb75ec69fecc3629`

## Command (`giap.mapping.contaminants.blastp_cmd(proteins.faa, contaminants_mini.faa, 2)`)
```
blastp -query proteins.faa -subject contaminants_mini.faa -num_threads 2 -seg no -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen nident positive" > blastp.tsv
```

## Output
- `blastp.tsv`: 48 HSPs; sha256
  `0125008777e9cfe2fef4e563d4f7fdec5384b81490256ac9bc79c2cba949822d`.
  Subject IDs are reported as `sp|ACC|NAME` (no database build, no `-parse_seqids`);
  `match_contaminant_proteins` resolves them to the accession like the protein FASTA reader.
  With the default thresholds (0.90 identity, 0.50 coverage) all eight `HIT_` queries match their
  own entry (identity 0.931-0.987; `HIT_P13645_1` covers 259/380 = 0.68); `LOWCOV_` (coverage
  0.40, or 0.52 identity for KRT10) and `RAND_` queries match nothing.
