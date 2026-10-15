# bakta 1.12.2-light on sim100k (M2-10, GIAP_SPEC.md §11.4)

- Command: `bakta --db /home/classes/bios/270/khoa/giap/local/db/bakta/db-light --output /home/classes/bios/270/khoa/giap/local/recordings/bakta_light_circ/sim100k/bakta --prefix sim100k --threads 16 --keep-contig-headers --replicons /home/classes/bios/270/khoa/giap/local/recordings/bakta_light_circ/sim100k/preflight/work/preflight/replicons.tsv /home/classes/bios/270/khoa/giap/local/recordings/bakta_light_circ/sim100k/preflight/work/preflight/genome.norm.fa`
- Image: `sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` (versions.txt: bakta 1.12.2)
- Database: Bakta DB light v6.0 (2025-02-24, doi 10.5281/zenodo.14916843)
- Date: 2026-10-15T17:03:14Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID 1849098; 16 CPUs; wall time 34 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/bakta_light_circ/sim100k`
- Inputs:
  - `tests/data/genomes/sim100k.fna.gz`: sha256 `f08b406acbc55bb9c5403b2ab04d5785df6383d147ebf328e78a25e00b3b0ad1`
  - `/home/classes/bios/270/khoa/giap/local/recordings/bakta_light_circ/sim100k/preflight/work/preflight/genome.norm.fa`: sha256 `369d1da812bd925afe260cad42a51eb1181d4f358f4d94db6ee38e05516f2a7f`
- Bakta outputs from `bakta/sim100k.*`; raw JSON sha256 `31768946d24932fd522b3d7548a103e7ffda6243a3e7ddc0d1a6d9f372afdb9b`
- Preflight (`giap preflight`) wrote replicons.tsv and contig_map.tsv from `sim100k.fna.gz`

## Files (sha256sum)

```
001f3f55a5316894f3ebfc14468693f0a1064dbb2bff5a18c794c6f78c116077  contig_map.tsv
c57ccafedbef1a689e02b46bd0b051f531d8fe382ec0db654c6b83a6e6a455e2  replicons.tsv
779baf71942913ba055538be2fc2948218b9b59b3ccf7edb35dbc426de120473  sim100k.faa.gz
713d90687e4ab2c1ae1a89a424b0c4e77155849bbb84051567adfb32477b64b7  sim100k.gbff.gz
13908bb393f5270efc4447893006675a9cadb92ddcfb040ee0b9360b7a8f41b3  sim100k.gff3.gz
0f8d7b72d993ba256f1be6dd8a8217cec8e4ef712128bf2b5fdd161b98acb702  sim100k.json.gz
7b401ba328d11fafe1909b127dfd48b4fd6bd25fde70ec85cc6af9c3dd0e41ba  sim100k.tsv.gz
```
