# Test genomes (M2-10; GIAP_SPEC.md §11.4)

- Command: `python scripts/fetch_test_genomes.py --out tests/data/genomes --compare /home/classes/bios/270/khoa/giap/local/data/genomes`
- Date: 2026-10-15T16:46:20Z
- Host: iron-01
- Source: NCBI Datasets v2 (`https://api.ncbi.nlm.nih.gov`), assemblies:
- `k12_mg1655.fna.gz`: GCA_000005845.2
- `nissle1917.fna.gz`: GCA_043228065.1
- Cross-check: NCBI FTP downloads of `scripts/slurm/fetch_resources.sbatch`, fetched
  20261014 (2026-10-14) into `$GIAP_DATA_DIR/genomes/<name>/` (md5 verified).
- Headers are NCBI's with ` [topology=circular]` appended to every replicon (chromosomes and plasmids
  of these complete assemblies are circular; lead decision 1). Sequences are unchanged.
- `sim100k.fna.gz`: giapsim `simulate(SIM100K, SIM100K_SEED)` genome FASTA, verbatim.
- Gzipped without name or mtime (`gzip -n -9` equivalent).

## Sequences (sha256 of the upper-case sequence, no header or newlines)

- `k12_mg1655`: sequences identical to the NCBI FTP file `GCA_000005845.2_ASM584v2_genomic.fna.gz`.
  - `U00096.3` 4,641,652 bp, sequence sha256 `e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0`
- `nissle1917`: sequences identical to the NCBI FTP file `GCA_043228065.1_ASM4322806v1_genomic.fna.gz`.
  - `CP171242.1` 5,051,730 bp, sequence sha256 `14d7914bdd9073ccc04beab477fc4b0341d17ce3180be59fb46b9e4ea370d3fb`
  - `CP171243.1` 4,523 bp, sequence sha256 `40c11dbad5ab475d2b9f87cccf56190c705c018840ea3eb807064e90bc8d5c34`
  - `CP171244.1` 3,173 bp, sequence sha256 `7278bc2cf0fd2a8af41f3d7f3e6751a5d0038651526e2fcc1a133a6c11ae7ee8`
  - `CP171245.1` 5,514 bp, sequence sha256 `89b32fa4e38e2127d0c664bc37b3e1c041522fa060995617275d1aa44dc677ce`

## Files (sha256sum)

```
4bfd7163f6f0de6394a2131cbfdbe75100daa55cbee227ec3e403a44640ac03c  k12_mg1655.fna.gz
fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee  nissle1917.fna.gz
f08b406acbc55bb9c5403b2ab04d5785df6383d147ebf328e78a25e00b3b0ad1  sim100k.fna.gz
```
