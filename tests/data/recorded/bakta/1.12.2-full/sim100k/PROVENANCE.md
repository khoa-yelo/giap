# bakta 1.12.2-full on sim100k (M2-10, GIAP_SPEC.md §11.4)

- Command: `bakta --db /home/classes/bios/270/khoa/giap/local/db/bakta/db --output /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/sim100k/bakta --prefix sim100k --threads 16 --keep-contig-headers --replicons /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/sim100k/preflight/work/preflight/replicons.tsv /home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/sim100k/preflight/work/preflight/genome.norm.fa`
- Image: `sha256:1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18` (versions.txt: bakta 1.12.2)
- Database: Bakta DB full v6.0 (2025-02-24, doi 10.5281/zenodo.14916843)
- Date: 2026-10-15T17:02:36Z (run finished, UTC)
- Host: barley-01 (FarmShare; SLURM_JOB_ID 1849096; 16 CPUs; wall time 317 s)
- Recording: `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/sim100k`
- Inputs:
  - `tests/data/genomes/sim100k.fna.gz`: sha256 `f08b406acbc55bb9c5403b2ab04d5785df6383d147ebf328e78a25e00b3b0ad1`
  - `/home/classes/bios/270/khoa/giap/local/recordings/bakta_full_circ/sim100k/preflight/work/preflight/genome.norm.fa`: sha256 `369d1da812bd925afe260cad42a51eb1181d4f358f4d94db6ee38e05516f2a7f`
- Bakta outputs from `bakta/sim100k.*`; raw JSON sha256 `a7b40a6824eedea71032b3b5e61102f536a1a8824a26f3e24105921a690151cb`
- Preflight (`giap preflight`) wrote replicons.tsv and contig_map.tsv from `sim100k.fna.gz`

## Files (sha256sum)

```
001f3f55a5316894f3ebfc14468693f0a1064dbb2bff5a18c794c6f78c116077  contig_map.tsv
c57ccafedbef1a689e02b46bd0b051f531d8fe382ec0db654c6b83a6e6a455e2  replicons.tsv
779baf71942913ba055538be2fc2948218b9b59b3ccf7edb35dbc426de120473  sim100k.faa.gz
b162cd737eee5fe5b905f0a520090cc7a29c11979b1690b00c4cb09a56d28d6b  sim100k.gbff.gz
06d072108311898c101dbe65b3ff1f421c04f6dd299026bb1f15892edf887ec9  sim100k.gff3.gz
5c3d3d31993176a0d434b67676857aa231b90d131d58d897c25942fe1bd3dfb4  sim100k.json.gz
21562adac2ea5532b723c22cc9673986601140c30a23860b0b71b69cb475ce8a  sim100k.tsv.gz
```
