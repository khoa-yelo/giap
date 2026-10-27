# Bundled sequences of giapsim (issue M2-05, lead decision Q1; rules R6/R12)

Real *Escherichia coli* str. K-12 substr. MG1655 rRNA and tRNA genes, inserted verbatim by
`SimConfig(real_rrna_trna=True)` so that Bakta's own rRNA/tRNA searches find them.

| | |
|---|---|
| Source | NCBI Nucleotide (RefSeq), official E-utilities `efetch` |
| Accession.version | NC_000913.3 (*E. coli* str. K-12 substr. MG1655, complete genome) |
| Retrieved | 2026-10-14 |
| Licence | NCBI places no restrictions on the use or distribution of GenBank/RefSeq data (https://www.ncbi.nlm.nih.gov/home/about/policies/) |

URLs (gene coordinates taken from the RefSeq feature table of the same record):

- feature table: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=ft&retmode=text`
  (sha256 of the file as retrieved: `a49fcaf434ea042f66d139f78879fbc241e74657067e524a129c78b8edb03a82`)
- each gene: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text&seq_start=<start>&seq_stop=<stop>&strand=<1|2>`
  with the 1-based closed coordinates in each FASTA header (`c` = complement strand, `strand=2`).

Files (sequences unchanged, upper case, 70 nt per line; headers rewritten to
`<gene> NC_000913.3:<coordinates> <annotation from the feature table>`):

| File | Content | sha256 |
|---|---|---|
| `ecoli_rrnB.fa` | rrnB operon rRNA genes: rrsB (16S, 1,542 nt, ends `GATCACCTCCTTA`), rrlB (23S, 2,904 nt), rrfB (5S, 120 nt) | `661c5d10031aa1f65a173fdbf117100be788e83761f2d439045fcd73a37c8d22` |
| `ecoli_trna.fa` | 16 tRNA genes, one per amino acid listed (alaV, argX, aspT, gltT, glyV, hisR, ileT, leuT, lysT, metT, pheV, proM, thrU, trpT, tyrT, valT); anticodon (DNA alphabet) and amino acid from the feature table notes | `3c8c30dc01870ca6f5b2487e47c9d9afb5b9c0a459d727ae913f123900efaf69` |

Everything else giapsim writes (genes, spacers, ITS, CRISPR repeats, the oriC-like region) is
synthetic and drawn from the simulator's RNG.
