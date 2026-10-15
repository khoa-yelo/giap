# miniprot 0.18 recorded run: giapsim SIM100K (M3-12 review decision 2)

- Tool: miniprot 0.18-r281 (MIT)
- Image (pinned in workflow/conf/containers.config):
  `docker://quay.io/biocontainers/miniprot@sha256:2eb53fea53743b56b1c10ac4b202469b92db959ffbb8bafd83156fdbde80e22f`
  (quay.io/biocontainers/miniprot:0.18--h577a1d6_0), run with Apptainer 1.4.4 from `$GIAP_ROOT/containers/`
- Date: 2026-10-14
- Host: iron-01 (Stanford FarmShare, Slurm allocation), 2 threads

## Inputs
- Genome: `simulate(SIM100K, SIM100K_SEED, outdir).genome_fasta` (giapsim, tests/sim), not stored here;
  sha256 `dd10216ee21cb0ac13940b2f7d5a26bf0fdc8c2af54af8baf26b380fc65a8d69` (`genome.fasta`)
- Proteins: `proteins.faa`, the 82 CDS/sORF translations > 40 aa of the simulated `bakta_like.json`
  (IDs are the locus tags, terminal `*` removed, 60 residues per line);
  sha256 `620cb22f254531101c2acce103d58f82a0caed15a686338f809b8fcb8806128e`

## Commands (GIAP_SPEC.md §7.1.1: -T 11 only at indexing)
```
miniprot -t 2 -T 11 -d genome.mpi genome.fasta
miniprot -t 2 -S -N 5 --outs=0.95 --gff genome.mpi proteins.faa > miniprot.gff
```
(`-S -N 5 --outs=0.95` = default `mapping.proteins.miniprot_args`.)

## Output
- `miniprot.gff`: 82 mRNA; sha256 `156d05256e56a3e448a1a535801d6a7c57a54cc8fc4ab5fe63f75db4ae14d8a2`.
  10 alignments open with an alternative start (`cs:Z:*gtgM` x7, `*ttgM` x3).
