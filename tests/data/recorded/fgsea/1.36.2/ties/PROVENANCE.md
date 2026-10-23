# ties: real fgsea 1.36.2 output (giap_fgsea.R)

- Input: the T5 tie ranking: 50 TUs, 20 tied at the top
- Command (in this directory, giap-fgsea env at /home/classes/bios/270/khoa/giap/local/micromamba/envs/giap-fgsea): `giap_fgsea.R --stats stats.tsv --sets sets.tsv --jobs jobs.tsv --sample-size 101 --eps 0.0 --n-perm-simple 1000 --seed 20261022 --out results.tsv --versions versions.txt`
- Env lock workflow/envs/locks/fgsea.explicit.txt sha256 d91eaaca1705fd828bef422dc736f2da8b15368fd82b8b22b5f08c5a293c6895
- Date 2026-10-22T20:51:03Z, host wheat-01, Slurm job 1868708, giap f407a71
- Recorded by scripts/dev/record_fgsea.py (scripts/slurm/record_fgsea.sbatch).
