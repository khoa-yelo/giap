# UniProt K-12 proteome UP000000625 excerpt (M3-22 TE8 CI test; GIAP_SPEC.md §11.4, §11.6, §15.3)

- Command: `python -P tests/helpers/uniprot_excerpt.py $GIAP_DATA_DIR/proteomes/UP000000625.fasta tests/data/recorded/uniprot/2026_03/UP000000625/k12_excerpt.fasta`
  (deterministic: entries sorted by accession, every 110th plus P0CE47 tufA and P00722 lacZ;
  43 records written verbatim)
- Source (official, R6): `https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000000625%29`,
  UniProt release 2026_03, 4403 entries, downloaded 2026-10-15T01:11:53Z by
  `scripts/slurm/fetch_resources.sbatch` (row `proteome_UP000000625` of `$GIAP_DB_DIR/prefetch_manifest.tsv`)
- Inputs: `$GIAP_DATA_DIR/proteomes/UP000000625.fasta` sha256 `a8052a607a98d2444c9d5dd685bd1370150df8c5570c3116811544bbf3d093a2` (not committed, R11)
- Date: 2026-10-20 (excerpt built)
- Host: iron-01 (FarmShare, inside the session Slurm allocation; no network used)
- Licence: UniProt data are CC BY 4.0 (The UniProt Consortium); the excerpt is unmodified.
- The full proteome is never committed; the weekly TE8 run fetches it with `validation/fetch_uniprot_k12.py`.

## Files (sha256sum)

```
deede156effbd4e6ad5771aa532235461d75dddea263d2b0f9fd487d878601ab  k12_excerpt.fasta
```
