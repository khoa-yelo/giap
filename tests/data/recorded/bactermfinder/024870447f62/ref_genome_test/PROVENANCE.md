# BacTermFinder genome scan of its bundled test genome: excerpt (fix-bactermfinder-parser)

Excerpt by the test-author on 2026-10-26. Source: the output that BacTermFinder ships with its
repository at the pinned commit `024870447f6288ebc9d29b314e957a6dfdfeb5f6` (ADR-0008; the env's
`share/bactermfinder/GIAP_COMMIT`). Upstream commit `127b404e` ("added a sample file input and
outputs", 2024-12-04) added it, and it is unchanged at the pin (`git status` clean). The file was
read from `$GIAP_ROOT/micromamba/envs/giap-bactermfinder/share/bactermfinder/` (env from
`workflow/envs/bactermfinder.yml`, sha256 `cf9677bba5289b37a7b6db06c4a4e0396bc4ae16ab120a466b84487e0e706702`;
no lock file exists yet, M11-14b).

- Command (README example at the pin): `python genome_scan.py ref-genome-test.fasta 3 out 10000 > log.out`
- Host: upstream authors' machine (not recorded upstream). The same env reran the script on FarmShare
  in Slurm job 1895829 (node rye-02), see `../smoke_1895829/`. Its first two rows equal rows 1-2
  here to within 6e-8.
- Inputs: `ref-genome-test.fasta` (CP054306.1, 630 bp; first 630 bp of Synechocystis sp. PCC 7338)
- Source table: `outref-genome-test-CP054306.1-Synechocystis-PCC-7338.fasta_mean.csv`, 354 rows,
  sha256 `a40499e57401eefd217afc54b096b2f375c862fd809740f1a15b51a9cb066b60`.
- `genome_scan.py` sha256 `697c4896f8c2fc0dd69430c1d9003f94f73ef232f6b952183c14bb329ec1352d`.

The excerpt keeps source lines 1 (header), 2-77 (windows 1-38), 158-159 (window 79, the highest
`probability_mean`, 0.1909) and 314-355 (windows 157-177, the last window `528_629`). That is 120
rows copied byte for byte (`sed -n '1p;2,77p;158,159p;314,355p'`).

## Format (genome_scan.py at the pin)
`SampleName,probability_binary,probability_ENAC,probability_PS2,probability_NCP,probability_mean`.
`SampleName` = `<u_id>_<seq_id>_<start>_<end>_<strand>` (lines 81-83). `u_id` is a 1-based window
counter that restarts for each FASTA record. `seq_id` is the Biopython record id, which may contain `_`
(e.g. `NZ_CP0001.1`), so split off `u_id` at the first `_` and the other three fields at the
last three `_` (from the right). Window `i` covers `seq[i:i + 101]`, stepping by 3, with
`start = i` and `end = i + 101` (lines 44-53). So `start`/`end` are **0-based half-open**, the
same as GIAP internally (§6.5). Each window appears twice, `+` then `-`, with the same
coordinates (the `-` row scores the reverse complement). Upstream's README script writes them to
BED unchanged. Every window is reported. The README suggests 0.13 / 0.3 / 0.47 cut-offs on
`probability_mean`.

## Files
```
bd8a44678ca35031042d2850f1f6c053298145d4b86d633705895c2fd01fdb9c  outref-genome-test.fasta_mean.excerpt.csv
87b7d8022a7e2379cca414981f774282735414a656e7bb22516c294d0404671a  ref-genome-test.fasta
```
