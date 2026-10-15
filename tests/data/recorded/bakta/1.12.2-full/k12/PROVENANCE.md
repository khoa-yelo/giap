# bakta 1.12.2-full on k12 (M2-10, GIAP_SPEC.md §11.4)

- Command: `bakta --db /home/classes/bios/270/khoa/giap/local/db/bakta/db --output /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/k12_mg1655/bakta --prefix k12_mg1655 --threads 16 --keep-contig-headers --replicons /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/k12_mg1655/preflight/work/preflight/replicons.tsv /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/k12_mg1655/preflight/work/preflight/genome.norm.fa`
- Image: `sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` (versions.txt: bakta 1.12.2)
- Database: Bakta DB full v6.0 (2025-02-24, doi 10.5281/zenodo.14916843)
- Date: 2026-10-15T16:57:17Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID 1849096; 16 CPUs; wall time 592 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/k12_mg1655`
- Inputs:
  - `tests/data/genomes/k12_mg1655.fna.gz`: sha256 `4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c`
  - `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/k12_mg1655/preflight/work/preflight/genome.norm.fa`: sha256 `918f35c7bfc8a3c8d824ed281a4185903025a0361242219b6b0e0ea6ecaa14f2`
- Bakta outputs from `bakta/k12_mg1655.*`; raw JSON sha256 `ae0ccaf186b098ecd704b78c76cde8e464f8717030bd5890b38397f1c0e9b169`
- JSON keys stripped (R11, M2-10 Q1): sequences.nt; the genome is tests/data/genomes/<sample>.fna.gz
- Preflight (`giap preflight`) wrote replicons.tsv and contig_map.tsv from `k12_mg1655.fna.gz`

## Files (sha256sum)

```
3c9f31147afe19dfcc8b34f03f709855ab403cf20c8786705181020ab17c1205  contig_map.tsv
79ee60bd1d06377b6a18690abe3c15f9c057fc0d024670041da18a019eabc161  k12_mg1655.faa.gz
3acfbf6eb75c15332145247255b90c82976e48b8cb0434a4333b974445a648bd  k12_mg1655.gbff.gz
9b8811a8af33fb075296e327826b3490741ed2a4e02962b733374f0460ec99b4  k12_mg1655.gff3.gz
c08c6466edf99d16cbe79679a8971b6e06e081cfc6adb9cc6f1ea9d51d3f083d  k12_mg1655.json.gz
9340127d4875e2b19c3787a1cd709adc3c7bb3391fb04308ea9af0626067fe2f  replicons.tsv
```
