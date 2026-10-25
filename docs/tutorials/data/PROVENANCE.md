# Demo peptide table: provenance

`sim100k_peptides.tsv.gz` is the quickstart's demo peptide table. It is **synthetic**: giapsim, the
GIAP test simulator (`tests/sim/giapsim`), generated every row. It contains no real mass
spectrometry data, no DIA-NN output and no pilot data.

- Genome: the table belongs to `tests/data/genomes/sim100k.fna.gz` (the SIM100K genome, one 100 kb
  circular chromosome and one 8 kb circular plasmid; see `tests/data/genomes/PROVENANCE.md`). The
  generator checks that the genome it simulates is byte-identical to that file.
- Simulator: giapsim at repository commit `f80fb217` (last change to `tests/sim/giapsim`:
  `e2e847bd`), Python 3.12.3, NumPy 2.5.4.
- Configuration: `SIM100K` with `proteomics=ProteomicsConfig(n_runs=1)`, all other proteomics
  settings at their defaults (tryptic peptides of 7–30 aa, 1 % wrong peptides, unscored decoys
  that name `rev_` proteins, `search_db = same`).
- Seed: `SIM100K_SEED` = 20261014.
- Command, from a repository checkout (no network access):

  ```bash
  python scripts/docs/make_demo_peptides.py           # writes the table
  python scripts/docs/make_demo_peptides.py --check   # exits 1 if the committed table differs
  ```

- Gzip: Python `gzip.compress(..., compresslevel=9, mtime=0)`; `--check` compares the
  uncompressed bytes.
- Contents: 1,685 rows (1,668 target and 17 decoy rows) of one run, `run1`, in the generic peptide
  table format (see [Inputs › Peptide table](../../user-guide/inputs.md#peptide-table-peptidestsv)):
  `peptide`, `protein_ids`, `q_value`, `pep`, `score`, `is_decoy`, `intensity`, `run`, `charge`,
  `spectrum_id`, `source`, `search_db`.
- Size: 68,421 bytes.

## Checksums (sha256)

```text
b7b7f7b11eda77383d25dfe1558ce1beb36dd430850628f74fb507c1a928654d  sim100k_peptides.tsv (uncompressed)
4894177927b371ad472d5cabb3c12155df965346f22a063f12993a31fc09da3b  sim100k_peptides.tsv.gz
```
