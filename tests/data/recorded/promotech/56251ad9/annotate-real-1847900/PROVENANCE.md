# Promotech RF-HOT targeted (`-s`) run: real ANNOTATE review run, Slurm job 1847900

Copied unmodified (test-author, 2026-10-15, M2-13) from the Nextflow work directory
`local/scratch/annotate-real/1847900/work/35/04a2f2d14dbd7caf4804a96cde321d/` of the reviewer's real
`ANNOTATE` run (`-profile farmshare`, conda env `giap-promotech`):

- `windows_000.fa`: the input chunk (86 records of exactly 40 nt; headers are locus tags written by
  the pre-M2-13 stub window writer, not `w<index>`).
  sha256 `230bcde33234ee197c64f1bbe8211051b087c90431baddd65ffd3c336c27a637`
- `sequences_predictions.csv`: Promotech's output `promotech_000/sequences_predictions.csv`
  (tab-separated despite the name; header `CHROM<TAB>SEQ<TAB>PRED`; one row per input record).
  sha256 `f1b1025f287172c463e8cff1774805b13f9a0dbcb4c7a2e09ec310a6adb4813c`

Tool: Promotech commit `56251ad9b883ef831b4753fc623d5ec970fe65e0` (`versions.yml`,
`$PROMOTECH_HOME/GIAP_COMMIT`), model `RF-HOT.model` from `RF-HOT.zip` sha256
`fcd96ffe5a43b579af965839ee2bb896de94de3601978d88bcabd32c4fb5093f`.

Command (`.command.sh`; cwd holds `models/` linking the env models plus `RF-HOT.model`):

```
python $PROMOTECH_HOME/promotech.py -s -m RF-HOT -f windows_000.fa -o promotech_000/
```

Reader in the pinned source: `sequences/process_sequences.py` parses the FASTA with
`Bio.SeqIO.parse(..., "fasta")` (record id = first whitespace-delimited token) and writes
`pd.DataFrame({"CHROM", "SEQ", "PRED"}).to_csv(..., index=None, sep='\t')`.
