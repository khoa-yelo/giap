# terminatornet 721803095459 on k12 (M2-10, GIAP_SPEC.md §11.4)

- Command: `python /home/classes/bios/270/khoa/giap/local/micromamba/envs/giap-terminatornet/share/terminatornet/TerminatorNet.py -f genome.fa -s -n 17 -o terminatornet.tsv` (cwd holds model.pkl -> $TERMINATORNET_HOME/model.pkl)
- Env lock: `workflow/envs/locks/terminatornet.explicit.txt` sha256 `d437c79205d1b313c5d060d0bff8c3fa676c6289de7c20011bbc598c8ba27c16`; tool commit `72180309545988dba63b8bd5db5d0b09bf4dba70`
- Date: 2026-10-15T06:41:54Z (run finished, UTC)
- Host: wheat-04 (FarmShare; SLURM_JOB_ID not recorded; 17 CPUs; wall time 3603 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/terminatornet/k12_mg1655`
- Inputs:
  - `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full/k12_mg1655/preflight/work/preflight/genome.norm.fa`: sha256 `918f35c7bfc8a3c8d824ed281a4185903025a0361242219b6b0e0ea6ecaa14f2`
- `terms.tsv.gz` = `terminatornet.tsv` verbatim (sha256 `673078256ba58086e91d2398d7fee9d5bfb49e3b3d539221191eb1923810cd4d`)

## Files (sha256sum)

```
9be61c75c0d7835f7d8c8f3c3948b95adfa34445003ec570ff9bbb2ca1e3a2ff  terms.tsv.gz
```
