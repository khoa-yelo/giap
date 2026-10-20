# Recorded UniOP run on K-12 MG1655: excerpt (M2-21)

- Excerpt (test-author, 2026-10-16, M2-21) written by `../m2_21_make_uniop_excerpt.py` from the real
  recording `$GIAP_ROOT/recordings/uniop/k12_mg1655/` (Slurm job `1849513_8`, node `barley-03`,
  17 CPUs, wall 7 s; `record_tools.sbatch`, `nextflow run workflow/tests/real/uniop.nf -profile
  farmshare_local,apptainer`).
- Tool: UniOP, git commit `3ed14dfc2a5d69b993d99a26e50c1f8e9723c087` (pin 3ed14dfc; ADR-0008), run
  as a separate process. Env locks: `uniop.explicit.txt` sha256
  `a16e1731537f901edfff5c1b7604b40e05de40ba8a062c7d4647ef54b7ff463b`, `uniop.pip.txt` sha256
  `3a2d491cf84fd08caf79b73d9ca48ea06442975cb41f3f5cdc17a3f8b1b5eaff`.
- Input written by GIAP (`record_tools_inputs.py uniop k12_mg1655.json genome.norm.fa .`, giap
  `fa19c16`) from the recorded full-DB Bakta 1.12.2 K-12 JSON: 4,372 genes on `U00096.3`, in start
  order. Command (`.command.sh`):

      python $UNIOP_HOME/src/UniOP -f uniop_cds.gff --gff_source prodigal_gff -t out

- Kept: the first 300 genes (`uniop_cds.gff` header + 300 CDS lines, `uniop_map.tsv` header + 300
  rows) and `uniop.pred` title, header and rows `1 2 … 299 300` (row `k` pairs gene `k` with gene
  `k+1`, both in the excerpt). Lines copied byte-for-byte.

| file | source sha256 | excerpt sha256 |
|---|---|---|
| `uniop_cds.gff` | `ea5870894e26c649b4e7520c65cf57e1138e1786688db7663a7b3e3eea062906` | `651143d405778e4e30918df2fd93cf27bd04bd5b3168dddd5d44f546fce1bbd2` |
| `uniop_map.tsv` | `7c1586c182aae9a771d637bc2a60e9a89f7e18583a2befd55e0a8ebcbfca45b1` | `094a7523c694f3fca4232fc752bcf999e3b3c850126960c8f8e40b44493c65da` |
| `uniop.pred` | `7d852a93a436a9bcef2cf4cf38d608bcd5d59054d98bd1b1ddc21196690b5ee0` | `8d8a1ce7891e4b675fafec2b7970cd2b254d697c2e2a14e54a1562e4171358eb` |
