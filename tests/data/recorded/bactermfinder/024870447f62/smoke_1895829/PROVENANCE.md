# BacTermFinder smoke run on FarmShare (Slurm job 1895829): head of the result table

Recorded by the test-author on 2026-10-26 from the job log
`local/logs/giap-btf-smoke.1895829.out` (sha256
`1bf47822bac6bf3cd0304cfaa3ef3f2dde25aef373e6007c0942380393e3dd6a`). The job's working directory
was a `mktemp` dir deleted on exit, so only the 5 lines it printed (`head -5`) survive. They are
copied byte for byte.

- Command: `micromamba run -n giap-bactermfinder python genome_scan.py tiny.fa 10 out 1000 > log.out`
  (cwd = a copy of `genome_scan.py`, `iLearnPlus/` and the four `*.h5` models).
- Commit: `024870447f6288ebc9d29b314e957a6dfdfeb5f6` (`GIAP_COMMIT`). Env from `workflow/envs/bactermfinder.yml`
  (sha256 `cf9677bba5289b37a7b6db06c4a4e0396bc4ae16ab120a466b84487e0e706702`). tensorflow 2.8.0,
  keras 2.8.0, numpy 1.23.0, pandas 1.5.0, biopython 1.81.
- Host: rye-02, Slurm job 1895829, 2026-10-26T01:01:58 to 01:40:35, wall 810 s.
- Inputs: `tiny.fa` is the job script's `>tiny CP054306.1 first 3000 bp` + the first 3000 bp of
  the bundled `ref-genome-test.fasta`, i.e. all of its 630 bp. It is rebuilt here with the job's
  own Python lines. The full table has 106 rows (53 windows × 2 strands, step 10).
- `SampleName` format and coordinates: see `../ref_genome_test/PROVENANCE.md`.

## Files
```
0f06593fc1e6ca678dfad92b42cb5e30446a8df2f45006c54ca71f89915493ad  outtiny.fa_mean.head.csv
22c20ebe68374bc5792584a89d891391f7e5febc551ed4d34fc9927617d4c34f  tiny.fa
```
