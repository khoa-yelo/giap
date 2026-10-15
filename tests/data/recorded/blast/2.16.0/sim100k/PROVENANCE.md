# BLAST+ 2.16.0 tblastn recorded run: giapsim SIM100K leftover proteins (M3-13 example 8)

- Tool: BLAST+ 2.16.0 (`tblastn: 2.16.0+`, `Package: blast 2.16.0, build Mar 28 2025 16:30:45`; NCBI, public domain)
- Image (pinned in workflow/conf/containers.config, selectors `BLAST` and `TBLASTN`):
  `docker://quay.io/biocontainers/blast@sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c`
  (quay.io/biocontainers/blast:2.16.0--h66d330f_5), run with Apptainer 1.4.4 from `$GIAP_ROOT/containers/`
- Date: 2026-10-14
- Host: iron-01 (Stanford FarmShare, Slurm allocation), 2 threads

## Inputs
- Genome: `simulate(SIM100K, SIM100K_SEED, outdir).genome_fasta` (giapsim, tests/sim), not stored here;
  sha256 `dd10216ee21cb0ac13940b2f7d5a26bf0fdc8c2af54af8baf26b380fc65a8d69` (`genome.fa`)
- Proteins: `proteins.faa`, written by `make_proteins.py` (this directory; seed 13): the 82 CDS/sORF
  translations > 40 aa of the simulated `bakta_like.json` (IDs are the locus tags, terminal `*`
  removed, 60 residues per line), each with random substitutions for a target identity drawn
  uniformly in [0.82, 0.93] (initial residue kept);
  sha256 `7af33b79ca2642098b34059d381f9d6d92a217a71a42ee7c66633d7dabb727b2`

## Commands (`-db_gencode 11 -seg no -evalue 1e-5` = default `mapping.proteins.tblastn_args`)
```
makeblastdb -in genome.fa -dbtype nucl -parse_seqids -out db/genome
tblastn -db db/genome -query proteins.faa -num_threads 2 -db_gencode 11 -seg no -evalue 1e-5 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen nident positive" > tblastn.tsv
```

## Output
- `tblastn.tsv`: 82 HSPs (one per protein); sha256
  `a4f52649b9e83813338bf8c29764fa629fc7f6e4b8c04ce7dd1659e7332a58f4`.
  After `merge_hsps` + `assign_alignments` (default thresholds) all 82 proteins are assigned
  `divergent` to their own gene (identity 0.818-0.928, coverage >= 0.969).
