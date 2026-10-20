# Promotech RF-HOT targeted (`-s`) run on K-12 MG1655: excerpt (M2-21)

Excerpt (test-author, 2026-10-16, M2-21) written by `../m2_21_make_promotech_excerpt.py` from the
real recording `$GIAP_ROOT/recordings/promotech/k12_mg1655/` (Slurm job `1849513_6`, node
`barley-03`, 17 CPUs, wall 236 s; `record_tools.sbatch`, `nextflow run workflow/tests/real/promotech.nf
-profile farmshare_local,apptainer`). Inputs were written by GIAP
(`record_tools_inputs.py promotech-windows k12_mg1655.json genome.norm.fa .`, giap `fa19c16`) from the
recorded full-DB Bakta 1.12.2 K-12 JSON: 950,228 windows of 40 nt, headers `w0 … w950227`, in ten
chunks `windows_000.fa … windows_009.fa`.

Tool: Promotech commit `56251ad9b883ef831b4753fc623d5ec970fe65e0` (`versions.txt`), model
`RF-HOT.model` from `RF-HOT.zip` sha256
`fcd96ffe5a43b579af965839ee2bb896de94de3601978d88bcabd32c4fb5093f`; env locks
`promotech.explicit.txt` sha256 `5a670f656cef63cc1d6a84220a79feea52dc47b89e10a95e31a65021bf7ffabb`,
`promotech.pip.txt` sha256 `6ae104da67d523a2b1a00bfd63970a326aa279ba5f11a32877cda3244dc3ddff`.

Command per task (`.command.sh`):

```
python $PROMOTECH_HOME/promotech.py -s -m RF-HOT -f windows_000.fa -o promotech_000/
python $PROMOTECH_HOME/promotech.py -s -m RF-HOT -f windows_009.fa -o promotech_009/
```

Kept (lines copied byte-for-byte; prediction header line kept):

| excerpt file | source | source sha256 | excerpt sha256 |
|---|---|---|---|
| `windows_000.fa` (records `w0 … w99`) | `windows/windows_000.fa` | `d3d44c45b3c90c0779c7c4b127adbbcc64f43ad24c7313e6c2a0370e42e784c5` | `8a41c419294a38c81b45e3842bce680733f348de2fbd2da05687f05a7ec17782` |
| `windows_000.sequences_predictions.csv` | `promotech_000/sequences_predictions.csv` | `d91d1c29fd265873b92085afa82105fe6d845ea51c1b5c73628fe181791b5477` | `0a8a5c0c05e8651b6b7b4442bbd979cc2fb348b779e4ad2767f0417c038b23dd` |
| `windows_009.fa` (records `w950128 … w950227`) | `windows/windows_009.fa` | `d992255745677f3b0a6a17aa122c053ecbde649eae72ae34447cbdef5687f5c2` | `53dc81dad58c79a4bbec64290418ec33d628148210e0419c71fc6f2b731cb531` |
| `windows_009.sequences_predictions.csv` | `promotech_009/sequences_predictions.csv` | `5af0a34feb78ef86485e7cb95dae9f70ff52a082c5216b9ca4ec8b5ace53dd30` | `61074b04f0ed71f2aa63ca68036d1a81183ca411db0d5ce0f1848c8a575bf0b5` |

Promotech scores each 40-nt record independently in `-s` mode, so the kept rows are exactly the
tool's output for the kept records.

Note: `$GIAP_ROOT/recordings/promotech/sim100k/` holds a genome-wide (`-pg`/`-g`) run, not the
targeted `-s` mode GIAP uses; its `genome_predictions.csv` reports `chrom` = `chr1,plasmid_1` for
every row (the two replicon IDs joined), so it is not used as a contract fixture.
