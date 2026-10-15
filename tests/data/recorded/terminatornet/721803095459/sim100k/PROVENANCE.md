# terminatornet 721803095459 on sim100k (M2-10, GIAP_SPEC.md §11.4)

- Command: `python /home/classes/bios/270/khoa/giap/local/micromamba/envs/giap-terminatornet/share/terminatornet/TerminatorNet.py -f genome.fa -s -n 17 -o terminatornet.tsv` (cwd holds model.pkl -> $TERMINATORNET_HOME/model.pkl)
- Env lock: `workflow/envs/locks/terminatornet.explicit.txt` sha256 `d437c79205d1b313c5d060d0bff8c3fa676c6289de7c20011bbc598c8ba27c16`; tool commit `72180309545988dba63b8bd5db5d0b09bf4dba70`
- Date: 2026-10-15T05:41:51Z (run finished, UTC)
- Host: wheat-04 (FarmShare; SLURM_JOB_ID not recorded; 17 CPUs; wall time 80 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/terminatornet/sim100k`
- Inputs:
  - `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full/sim100k/preflight/work/preflight/genome.norm.fa`: sha256 `369d1da812bd925afe260cad42a51eb1181d4f358f4d94db6ee38e05516f2a7f`
- `terms.tsv.gz` = `terminatornet.tsv` verbatim (sha256 `d0898b107e1f9d08f814c2e68c3bb502952aaa9b7fda336ae7d463f4f517c8f0`)

## Files (sha256sum)

```
ea70da2edfcdd46a4e19ff9b1d0e4976ee12d6a00b6825730f83686b0b9857fa  terms.tsv.gz
```
