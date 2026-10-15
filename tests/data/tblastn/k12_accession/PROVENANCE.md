# Accession-style contig ID, real tblastn run (M3-13 review decision 1)

Guards against `makeblastdb -parse_seqids`: with it BLAST+ reports an accession-style subject ID
as `ref|NC_000913.3|`, which never matches GIAP's contig ID; without it the first header token is
reported unchanged.

- Tool: BLAST+ 2.16.0 (`tblastn: 2.16.0+`), NCBI public domain
- Image (pinned in workflow/conf/containers.config):
  `docker://quay.io/biocontainers/blast@sha256:9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c`,
  Apptainer from `$GIAP_ROOT/containers/`
- Date: 2026-10-14; host iron-01 (FarmShare Slurm allocation), 2 threads

## Inputs
- `genome.fa`: one contig `NC_000913.3` (header `>NC_000913.3 Escherichia coli K-12 MG1655 excerpt
  757001-760000`), the 3,000 nt `[757000, 760000)` (0-based half-open) of E. coli K-12 MG1655
  U00096.3 (the code reviewer's K-12 probe genome, sha256
  `a6dfe06d30ddb353d5ed84a5888d149629c35fe9ea176ecc5df09a718e2bf438`);
  sha256 `4fb005381faacaee53f06e0c6f9601df292b10323574ace2e8710d4803aee052`
- `proteins.faa`: UniProt P37005 (228 aa) and P0AD86 (21 aa), copied from the reviewer's probe
  queries; sha256 `038e79c74a02f27dabe596b45e5992bef0f07890c6d45e7496c797267efa5918`
  - P37005: plus strand, excerpt `[941, 1625)`; P0AD86: plus strand, excerpt `[1841, 1904)`

## Commands
```
makeblastdb -in genome.fa -dbtype nucl -out db/genome
tblastn -db db/genome -query proteins.faa -num_threads 2 -db_gencode 11 -seg no -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen nident positive"
```
Contrast (negative example, not what GIAP runs): the same with
`makeblastdb -in genome.fa -dbtype nucl -parse_seqids -out dbp/genome` and `-db dbp/genome`.

## Output
- `tblastn.tsv` (2 HSPs, subject `NC_000913.3`); sha256
  `696d0db598617fb9d55d061758facdc1ab2fa75296f47213e8b970c9e95994ad`
- `tblastn_parse_seqids.tsv` (contrast; subject `ref|NC_000913.3|`); sha256
  `dba9a4984ea93952f52ff8604833d9a31cea085da11a3e1005aae6e426f500e9`
