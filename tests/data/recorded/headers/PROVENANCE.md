# Real output headers of the annotation tools (M2-07 lead decision 7 after the PR #65 review)

One file per tool output, `<tool>.<file>.txt`: the header line(s) of the real tool output, LF line
endings. The stub blocks in `workflow/modules/local/<tool>/main.nf` must write exactly these
headers (tests/unit/test_nf_stub_headers.py, workflow/tests/annotate.nf.test). Only the header is
recorded; no tool code is vendored.

Sources (real runs from the PR #65 review, `local/scratch/review-pr65`, 2026-10-14; inputs from
giapsim SIM100K, `mkinputs.py` there):

- `genomad.{virus,plasmid}_{summary,genes}.tsv.txt`: geNOMAD 1.12.0 (pinned container,
  `genomad end-to-end --cleanup --splits 8`), `genomad/genome_summary/genome_<k>.tsv`, first line.
  `virus_genes` / `plasmid_genes` have 20 columns.
- `ostir.ostir.csv.txt`: OSTIR 1.1.3 (`giap-ostir` env), `ostir.csv`, first line (11 columns,
  `RBS_distance_bp` after `expression`).
- `promoter_calculator.csv.txt`: Promoter Calculator 1.2.4 (`giap-promoter_calculator` env),
  `tu0001.csv`, first line (27 columns); the tool writes CRLF, the CR is stripped here. When the
  region has no promoter the tool writes no CSV at all ("No promotors found.").
- `promotech.sequences_predictions.csv.txt`: no real run yet (the review run failed on the missing
  PyQt5). Header from the pinned source, `share/promotech/sequences/process_sequences.py` lines
  60-64 (`pd.DataFrame({"CHROM", "SEQ", "PRED"})`, `to_csv(..., sep='\t', index=None)`) written to
  `<out>/sequences_predictions.csv`. Replace with the recorded header once the real run passes.
- `uniop.uniop.pred.txt`: first two lines of `tests/data/recorded/uniop/3ed14dfc/sim100k/uniop.pred`.
