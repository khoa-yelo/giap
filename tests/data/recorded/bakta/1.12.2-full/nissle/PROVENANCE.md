# bakta 1.12.2-full on nissle (M2-10, GIAP_SPEC.md §11.4)

- Command: `bakta --db /home/classes/bios/270/khoa/giap/local/db/bakta/db --output /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/nissle1917/bakta --prefix nissle1917 --threads 16 --keep-contig-headers --replicons /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/nissle1917/preflight/work/preflight/replicons.tsv /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/nissle1917/preflight/work/preflight/genome.norm.fa`
- Image: `sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` (versions.txt: bakta 1.12.2)
- Database: Bakta DB full v6.0 (2025-02-24, doi 10.5281/zenodo.14916843)
- Date: 2026-10-15T16:57:19Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID 1849097; 16 CPUs; wall time 595 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/nissle1917`
- Inputs:
  - `tests/data/genomes/nissle1917.fna.gz`: sha256 `fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee`
  - `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/nissle1917/preflight/work/preflight/genome.norm.fa`: sha256 `f64e01e8cb859c36221a97f730e6e40067c5abbdf5a861c62982c91d78485807`
- Bakta outputs from `bakta/nissle1917.*`; raw JSON sha256 `ac3e6190c4e761c3b15e9e91c20aca0efd4d76c20d0356b5c72303719665053c`
- JSON keys stripped (R11, M2-10 Q1): sequences.nt; the genome is tests/data/genomes/<sample>.fna.gz
- Preflight (`giap preflight`) wrote replicons.tsv and contig_map.tsv from `nissle1917.fna.gz`

## Files (sha256sum)

```
01fa56250e55344c67d375e47d1c6230469b1bbc87e0e49fd3597ab6236ba3e5  contig_map.tsv
936099687ebb5266abf828170716c94bf1fdac8f9d963405bbcbc7635c6f25aa  nissle1917.faa.gz
8e5cd1b9faa1a7791affeb2f93c5e0651dcb88e6e8318b16d874ed9e21cf99a7  nissle1917.gbff.gz
8b90bb44f06c786515150f94a47b60f27bb7367e7c1e508ef59019ee8e21b8a7  nissle1917.gff3.gz
a012798400eaa71587dc7954473a201ef046b4907bdfa31936a2927a1cb6c97c  nissle1917.json.gz
2567150e5d545ae8f0b20e052cb232b122ee040722c5b397e3751cb92c2d82f5  replicons.tsv
```
