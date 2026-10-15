# tests/data/recorded/makedb (M3-23 test fixture)

Excerpt of a real Bakta 1.12.2 run (full DB 6.0) recorded on FarmShare on 2026-10-14 under
`$GIAP_ROOT/recordings/bakta_full/sim100k/` (command in `command.txt` there; input = the giapsim
`sim100k` genome, a simulated 100-kb chromosome plus the 8-kb circular `plasmid_1`).

- `sim100k_plasmid_1.bakta.json`: `bakta/sim100k.json`
  (SHA-256 `d28c87810c2e3feed57c643c8845824b077fe252d8afb250ada86c42a99799cf`) reduced to the
  `plasmid_1` replicon: its `sequences` record (full 8,000 nt, circular) and all of its
  `features` (two CDS, KABHIH_0098 with a TTG start on the minus strand and KABHIH_0099), plus
  the top-level `genome` and `version` objects. Nothing else was edited.
- `sim100k_plasmid_1.bakta.faa`: the KABHIH_0098 and KABHIH_0099 records copied verbatim from
  `bakta/sim100k.faa` (SHA-256 `1cec03b4d95bb5bfd986277179d4546df3391d4229d8bb90b9d11aa068103b52`),
  i.e. Bakta's own protein FASTA for these loci; `giap makedb` must write the same IDs, products
  and sequences.

Extracted by the test author (Agent-Role: test-author) on 2026-10-15. The genome is simulated, so
there is no third-party sequence content; Bakta output itself is GPL-3 tooling output (data, not
code).
