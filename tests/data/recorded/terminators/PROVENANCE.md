# Recorded TerminatorNet and TransTermHP outputs (M2-15 test fixtures)

Copied verbatim (whole files or line excerpts) from the M2-10 / record_tools recordings under
`$GIAP_ROOT/recordings/` (FarmShare, `$GIAP_ROOT=/home/classes/bios/270/khoa/giap/local`) on
2026-10-15. Nothing was edited; excerpts are concatenated line ranges of the original files.

## TerminatorNet (commit 72180309545988dba63b8bd5db5d0b09bf4dba70, env `giap-terminatornet`)

Command (from `command.txt`): `python $TERMINATORNET_HOME/TerminatorNet.py -f genome.fa -s -n 17 -o terminatornet.tsv`
(cwd holds `model.pkl`). Real format: TAB-separated, no header, 11 columns
`name TERM_i A_tract hairpin U_tract structure GH pred start stop strand`; 1-based closed,
minus-strand calls already on forward coordinates (`hairpin == revcomp(genome[start-1:stop])`,
verified on all 17 071 calls of the three recordings); numbering restarts per strand;
a contig/strand without calls gives the 2-column line `<contig>\tNO TERMINATORS FOUND`.

| file | source | lines |
|---|---|---|
| `terminatornet/sim100k/terminatornet.tsv` (+ `command.txt`, `versions.txt`) | `recordings/terminatornet/sim100k/` whole (sha256 d0898b10…c4f8c0); genome = giapsim `SIM100K`/`SIM100K_SEED` (`chr1` 100 000 bp, `plasmid_1` 8 000 bp) | 133 calls + `plasmid_1 NO TERMINATORS FOUND` (minus strand) |
| `terminatornet/k12_excerpt.tsv` | `recordings/terminatornet/k12_mg1655/terminatornet.tsv` (sha256 67307825…7cd4d) lines 1-3, 4165-4169, 8293-8294 | first/last `+` calls, first `-` calls near the origin, last `-` calls (U00096.3, 4 641 652 bp) |
| `terminatornet/nissle_plasmids.tsv` | `recordings/terminatornet/nissle1917/terminatornet.tsv` (sha256 c319d52f…e071e), every line of CP171243.1, CP171244.1, CP171245.1 | 26 calls on 3 plasmids, both strands |

## TransTermHP 2.09 (image sha256:fc6fe040…dfc2d03)

Command: `transterm -p /usr/local/data/expterm.dat genome.fa genes.coords`; `genes.coords` written by
`scripts/dev/record_tools_inputs.py coords <bakta_full JSON>` (`<locus> <start1> <end1> <contig>`,
minus-strand genes with start > end). The recorded `genes.coords` of all three samples equal, line
for line, the spec rule applied to `parse_bakta_json(<bakta_full JSON>).model.genes` (CDS, sORF and
pseudogene, model order).

| file | source | notes |
|---|---|---|
| `transtermhp/sim100k/{transtermhp.tt,genes.coords,command.txt,versions.txt}` | `recordings/transtermhp/sim100k/` whole (tt sha256 4855b9ca…0cd4d, coords 5903eadd…7bdf87f) | 30 TERM on chr1, none on plasmid_1 |
| `transtermhp/sim100k/genes.tsv` | derived: `parse_bakta_json(recordings/bakta_full/sim100k/bakta/sim100k.json)` (sha256 d28c8781…799cf) genes, 0-based half-open | input of `write_transterm_coords` |
| `transtermhp/k12_excerpt.tt` | `recordings/transtermhp/k12_mg1655/transtermhp.tt` (sha256 0fe1eb3d…42272) lines 1-32, 458-461, 598-601, 4029-4032, 8070-8088 | header, TERM 1-3 near the origin, loc codes `GF`/`RH`/`Tf`, integer hp `-13`/`-7`, TERM 1845-1847 at the end; one `SEQUENCE U00096.3` |
| `transtermhp/nissle_plasmids.tt` | `recordings/transtermhp/nissle1917/transtermhp.tt` (sha256 5aa8a5f1…cbd202) lines 1-19, 8318-end | header + three plasmid sections (CP171243.1: 5 TERM, CP171244.1: 0, CP171245.1: 1) |
