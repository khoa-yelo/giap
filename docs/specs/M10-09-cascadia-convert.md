# M10-09: `gippa convert --from cascadia` (Cascadia 0.0.7 `.ssl`) and its contract test

Milestone M10 · **pass 2** · depends on: M10-08 (the `.ssl` format facts below; the real 1,000-row
sample exists only after M10-08's smoke run gives `verdict = pass`) · estimate ≈ 140 code lines
(index: 200) + ≈ 330 test lines, one PR · **not `stats`** (no estimator, no threshold; see
"Statistician notes") · roles: test-author (CI tests, the FarmShare contract test), implementer
(`src/gippa/io/converters/cascadia.py`, the registry, the CLI branch), lead (the
`PROVENANCE.md` section after the FarmShare check) · Khoa: none needed (his 2026-10-23 data rule is
applied as already decided for DIA-NN; question Q8).

The CI half (converter, CLI, synthetic tests) needs no real data and may start now. Only the
FarmShare evidence (criteria 12–14) waits for M10-08's smoke run. If M10-08 says "fails on
timsTOF", the converter still lands (§5.4 lists Cascadia), the contract test skips forever, and A3
moves to M10-11a/11b.

## Spec sections implemented

- §5.4 converters: "Each converter is a pure function from input files to the generic table, with
  a contract test against a recorded sample file (§11.4). Until a real sample exists, the
  converter is written from the tool's documentation and its contract test is marked
  `pending_sample`; the pilot (§15) supplies DIA-NN and Cascadia samples." Row "Cascadia 0.0.x |
  its result file (mzTab-like) | as Casanovo; `source` = denovo". The Casanovo row it points to:
  "`peptide` ← stripped `sequence`; `score` ← `search_engine_score[1]`; `spectrum_id` ←
  `spectra_ref`; `run` ← the `ms_run` location; `charge` ← `charge`; `source` = denovo". The real
  0.0.7 file is not mzTab-like but a BiblioSpec-style `.ssl` table (M10-08); "as Casanovo" is
  applied column by column below.
- §5.4 stripping: "removes bracketed or parenthesised modifications (`(UniMod:35)`, `[+15.995]`,
  `[UNIMOD:35]`, `+15.995`, `n[42]`) and lower-case residues, and keeps only uppercase letters"
  (`gippa.io.converters.strip.strip_modifications`, unchanged).
- §5.3 generic peptide table: columns and types; "`score` | float | engine or de novo score;
  higher is better"; "`spectrum_id` | scan or spectrum identifier, unique within `run` (G2)";
  "`q_value` | precursor or peptide q-value"; "I/L are kept as written; equivalence is applied
  only in matching".
- §5.1 `gippa convert --from cascadia INPUT… --out peptides.tsv`.
- §7.2.5 m = 0, de novo input: π̂_wrong is "1 − the lower 95 % bound of the de novo precision for
  the peptide's score bin"; "De novo precision … For score bin b (deciles of `score`)". §7.3.1:
  "it passes the engine threshold (`q_value ≤ input.q_value_max` **when a q-value exists**)". These
  two decide the score → q question (decision D3).
- §7.5.4 pairing de novo rows with `--reference-peptides` by `spectrum_id` (decision D2).
- §6.6 determinism and "No timestamps or absolute paths in `tables/`" (decision D1).
- §11.4 recorded outputs and `PROVENANCE.md`; §15.2 step 7, the Cascadia half, as amended by
  Khoa's 2026-10-23 rule (M10-02 lead decisions 4 and 6; M10-03 lead decision 2; M10-08 Q4):
  recorded samples and converted tables stay on FarmShare.
- §6.5 coordinates: **not applicable**. The converter reads and writes no genomic coordinate;
  circular/linear topology and strand do not arise (see "Edge cases").
- R2 (tests first), R5 (Cascadia is never imported: the converter reads its output file only),
  R11 (nothing row-level committed), R14 (no new Slurm job; the contract test runs in seconds).

## What Cascadia 0.0.7 writes (read from the pinned source)

Source: the PyPI sdist `cascadia-0.0.7.tar.gz`, sha256 `6e9fcc50…ab347` (the hash M10-08
records; re-downloaded and re-checked 2026-10-24 on `iron-13`, read only, nothing installed or
executed). Writer: `cascadia/utils.py::write_results`; vocabulary:
`cascadia/depthcharge/primitives.py::MSKB_TO_UNIMOD`.

| # | Fact (0.0.7 code) | Consequence for the converter |
|---|---|---|
| W1 | One header line, written literally: `file scan charge sequence score-type score retention-time start-time end-time` (tab-separated), then one line per kept prediction; `\n` line ends; no quoting | read with `gippa.io.tsv.read_tsv` (UTF-8, BOM dropped, `\r` tolerated, ragged rows → `ValueError`, `NA`/empty → missing); gzip accepted by magic bytes |
| W2 | `file` = the spectrum-file argument exactly as given on the Cascadia command line (M10-08: an absolute path, `$GIPPA_DATA_DIR/pxd078459/mzml/<RUN>.mzML`); one value per output file | `run` ← the file's stem (D1) |
| W3 | `scan` = `str(idx)`, the 0-based index of the prediction in Cascadia's prediction list, **not** a scan number. Rows that fail the filters are skipped, so the values are strictly increasing with gaps; each value occurs at most once per file | `spectrum_id` ← `cascadia:<scan>` (D2) |
| W4 | `charge` = `str(int(charge))` | `charge` (Int64), as written |
| W5 | `sequence`: I already folded to L (`replace_isoleucine_with_leucine=True`). Modifications in ProForma-name form are converted to rounded mass shifts by `convert_mods_to_shifts` (`'+' + str(round(mass, 2))`): `M[Oxidation]` → `M[+15.99]`, `N[Deamidated]`/`Q[Deamidated]` → `N[+0.98]`/`Q[+0.98]`, `C[Carbamidomethyl]` → `C[+57.02]`, N-terminal `[Acetyl]-` → `[+42.01]-`, `[Carbamyl]-` → `[+43.01]-`, `[Ammonia-loss]-` → `[-17.03]-`; `[+25.980265]-` (carbamyl + NH₃ loss, not in Unimod) stays as is. A name Unimod cannot resolve is left as the name. A `-` anywhere but right after a leading N-terminal token makes Cascadia skip the row | `modified_peptide` ← `sequence` verbatim; `peptide` ← `strip_modifications(sequence)`; no I/L folding by the converter |
| W6 | `score-type` = the literal `UNKNOWN` on every row | any other value is an error (D5): it is not a Cascadia 0.0.7 file |
| W7 | `score` = `str(np.exp(log_confidence))`, so in (0, 1]; a row is written only if `score > score_threshold` (**strict**), `-t` default 0.8. NaN fails `>` and is never written | `score` ← `score` (float). No `q_value`, no `pep` (D3). Rows at or below 0.8 are kept and counted (D4) |
| W8 | `retention-time`, `start-time`, `end-time`: minutes when the largest RT of the file exceeds 500 (in the seconds Cascadia derives), else seconds. The pilot runs are 539.9 s long, next to that limit | **never read**. §5.3 has no RT column, so the unit ambiguity cannot reach any GIPPA table (D6) |
| W9 | No protein IDs, no decoys, no intensities, no q-values | those columns are absent from the output |

The 0.0.7 docs (`file_formats.md`) describe the same nine columns; where docs and code differ the
code wins (M10-08).

## Decisions taken here (each reversible; the matching lead question has the same default)

- **D1 `run` = the stem of `file`.** Take the part after the last `/`, remove one trailing `.gz`
  (case-insensitive), then one trailing `.mzML` (case-insensitive). `/data/mzml/WT_1.mzML` →
  `WT_1`; `WT_2.mzml.gz` → `WT_2`; `run3` → `run3`. Reason: the verbatim value is an absolute
  FarmShare path, which §6.6 keeps out of `tables/` and which would differ between hosts for the
  same data. (For DIA-NN, M10-03 Q3 kept `Run` verbatim because DIA-NN already writes a stem.)
  An empty or missing `file`, or one whose stem is empty, is the error `bad_value:file` (with
  rows; lead, 2026-10-24). (Q1)
- **D2 `spectrum_id` = `cascadia:<scan>`.** W3: `scan` is a prediction index, not a scan number.
  Written bare, it could pair by accident with a reference search that uses real scan numbers in
  the same run (§7.5.4, preflight C17), and would then give wrong ALT error rates. The prefix keeps
  it unique within the run (G2), keeps the validator's duplicate checks working, and stops any
  false pairing. With A2 as the reference (DIA-NN writes no `spectrum_id`) nothing pairs either
  way, so the pilot is unaffected. A `scan` that is not a non-negative integer is an error. (Q2)
- **D3 Score → q/confidence: `score` only; `q_value` and `pep` are not written.** Cascadia gives
  no q-value and no posterior. §7.3.1 applies the engine threshold only "when a q-value exists",
  and §7.2.5 calibrates de novo rows from the genome by score deciles. So GIPPA needs `score` as a
  ranking, and nothing else. Writing `1 − score` as `pep`, or any q, would claim a calibration that
  Cascadia does not provide (see "Statistician notes"). It would also flip
  `steps.mapping.resolve_sources` to `database` on a table without a `source` column. (Q3)
- **D4 No score filter in the converter.** Filtering is not a converter's job (M10-03 F10: one
  output row per input row). A Cascadia run with the default threshold has no row ≤ 0.8. Rows
  with `score ≤ 0.8` (a run with a lowered `-t`, or a hand-made file) are kept, and an info issue
  `score_at_or_below_default_cutoff` counts them, so the deviation is visible. (Q4)
- **D5 Strict input checks.** Each is an error issue: the CLI exits 1 and writes nothing.
  - a missing required column (`file`, `scan`, `charge`, `sequence`, `score-type`, `score`);
  - a `score-type` other than `UNKNOWN`;
  - a missing, non-numeric or NaN `score`, or a score outside [0, 1] (`inf` included);
  - the cases in D1 and D2.

  The validator only checks that `score` is a float, and accepts `nan` text in `score` as
  missing. A de novo row without a score cannot be binned (bin −1, so π̂ = 1). It would pass
  silently, although Cascadia 0.0.7 can never write one. (Q5)

  Lead decision 2026-10-24: missing charge/sequence are left to the validator; Cascadia 0.0.7
  never writes them.
- **D6 RT columns ignored.** They are neither required nor parsed. A file without them converts.
  (Q11)
- **D7 CLI.** With `--from cascadia`:
  - `--source denovo` is accepted (a no-op); today the CLI rejects it for every non-generic
    engine, which is wrong for a de novo engine;
  - `--source database` exits 1 ("`--source database` is invalid with `--from cascadia`");
  - `--search-db` exits 1. The converter returns the error issue `search_db_invalid`: a de novo
    row has no search database, and `same` would skip placement FDR for exact hits (§7.2.5
    step 5);
  - `--abundance-out` exits 1 through the existing `no_abundance` error.

  `casanovo` and `instanovo` keep their exit 5 ("not implemented"). (Q7)

## The converter (`src/gippa/io/converters/cascadia.py`, implementer)

```python
def convert_cascadia(
    paths: Sequence[Path], *, search_db: str | None, cfg: InputConfig
) -> ConvertResult: ...
```

The signature is that of `convert_diann`, so the CLI's non-generic call path is unchanged. Added
to `CONVERTERS` as `"cascadia"` and to `__all__`. The function is pure apart from reading `paths`.

1. `search_db is not None` → error `search_db_invalid`; return an empty table.
2. Per path, in argument order: `read_tsv(path)`. `ValueError`/`OSError`/`UnicodeDecodeError` →
   error `unreadable` (a zero-byte file, a ragged row with its line number, a bad encoding). The
   header is matched by exact name; order is free; columns other than the nine of W1 → one
   `unknown_column` warn per file. A required column missing → error
   `missing_column:<name>`. A header-only file → info `empty_input` naming the file (it is not an
   error by itself; see step 8).
3. **Run** (D1) per row (no stem → error `bad_value:file` with rows); across inputs, a run that comes from more than one input file → error
   `run_in_several_inputs`, naming both paths. This covers the same `.ssl` passed twice, a re-run
   on the same mzML, and two mzMLs with the same stem in different directories.
4. **spectrum_id** (D2): `scan` must fullmatch `[0-9]+`, else error `bad_type:scan` with rows.
5. **score-type** (W6, D5): every non-missing value must be `UNKNOWN`, else error
   `unsupported_score_type` with rows and the first offending value. Missing → the same error.
6. **score** (D5): parse with `float()`. Missing, unparseable or NaN → error `bad_type:score`;
   finite but outside [0, 1], or infinite → error `out_of_range:score`. Rows with
   `score ≤ 0.8` → info `score_at_or_below_default_cutoff` (D4). The constant
   `CASCADIA_DEFAULT_CUTOFF = 0.8` is named in the module with its source (W7). It is not a
   config key: GIPPA never filters on it.
7. **charge**: through the peptide-table parser (`_typed`, Int64); non-integral → error
   `bad_type:charge`. The validator adds `out_of_range:charge` for < 1.
8. **Peptides:** `modified_peptide` ← `sequence` (string, verbatim); `peptide` ←
   `strip_modifications(sequence)` for every row whose `sequence` is not already `[A-Z]+`; an
   empty result becomes missing. An info issue `stripped` counts the stripped rows, as in the
   generic converter. No I→L or L→I change.
9. `source` = `denovo` on every row. `q_value`, `pep`, `protein_ids`, `is_decoy`, `intensity` and
   `search_db` are not produced.
10. Any error issue so far → return `ConvertResult(empty, None, issues)`. Otherwise
    `validate_peptides(df, cfg)` adds its issues (as `convert_diann` does): `empty_table` (all
    inputs header-only), `short_peptide`/`long_peptide`/`non_letters` warnings,
    `duplicate_spectrum_row` (warn), `duplicate_spectrum_id` (error). The returned table keeps
    **every** converted row; preflight drops what it drops.
11. Output: `finish_peptides(df, cols)` with `cols` = `peptide`, `modified_peptide`, `score`, `run`,
    `charge`, `spectrum_id`, `source` (that is §5.3 order). The sort is stable (mergesort) by
    `run`, `peptide`, `modified_peptide`, `charge`; ties keep the argument order, then the file
    order. `abundance` is `None`.

CLI (`src/gippa/cli.py::convert`): the D7 branch, a few lines; the rest of the command is
unchanged. Floats are written `%.6g` by `write_tsv` (§6.6), and `.gz` outputs have a fixed gzip
header, so the output bytes depend only on the input bytes.

## Examples

Rows are hand-made in the 0.0.7 shape (W1–W7), never copied from a real run. `⇥` marks a tab in
the input tables only; outputs are shown as TSV lines. Unless stated, `file` is
`/x/mzml/WT_1.mzML` and the RT columns are `4.5⇥4.4⇥4.6`. Default config.

| # | Input | Expected |
|---|---|---|
| 1 | `WT_1.ssl`: header (W1); `…⇥0⇥2⇥PEPTLDEK⇥UNKNOWN⇥0.95123⇥…`; `…⇥3⇥3⇥M[+15.99]PEPTLDEKR⇥UNKNOWN⇥0.88120001⇥…` | exit 0; `peptides.tsv` = <br>`peptide  modified_peptide  score  run  charge  spectrum_id  source`<br>`MPEPTLDEKR  M[+15.99]PEPTLDEKR  0.8812  WT_1  3  cascadia:3  denovo`<br>`PEPTLDEK  PEPTLDEK  0.95123  WT_1  2  cascadia:0  denovo`<br>(sorted by peptide; `%.6g` score; no `q_value`/`pep` column); info `stripped` (1 row) |
| 2 | Sequences `[+42.01]-SAMPLLER`, `[-17.03]-QLLDPEKR`, `AC[+57.02]DN[+0.98]LLK`, `[+25.980265]-GLLDPEKR`, `Q[+0.98]VLLDPEK`, `M[Oxidation]VLLDPEK` (an unresolved name, W5) | `peptide` = `SAMPLLER`, `QLLDPEKR`, `ACDNLLK`, `GLLDPEKR`, `QVLLDPEK`, `MVLLDPEK`; each `modified_peptide` verbatim |
| 3 | Hand-made `PEPTIDEK` (an I, which 0.0.7 never writes) and `PEPTLDEK`, same charge, scans 4 and 5 | two rows: `peptide` `PEPTIDEK` and `PEPTLDEK` (kept as written, §5.3); `modified_peptide` likewise. The matcher's I ≡ L makes them one genome query later; the converter merges nothing |
| 4 | Scores `0.8`, `0.5`, `0.80000001`, `1.0` on four rows | all four kept; info `score_at_or_below_default_cutoff` with rows 1 and 2 (2 rows). `0.80000001` is **not** counted (it is > 0.8) but is written `0.8` by `%.6g`; `1.0` is written `1` |
| 5 | Score cells `nan`, `NA`, `abc`, empty, `inf`, `1.2`, `-0.1` (one per row, plus one valid row) | exit 1, no output; `bad_type:score` rows 1–4; `out_of_range:score` rows 5–7 |
| 6 | `score-type` = `PERCOLATOR QVALUE` on one row | exit 1; `unsupported_score_type` (row 1, value `PERCOLATOR QVALUE`) |
| 7 | Scan `0` twice, same sequence and charge; then scan `0` twice with `PEPTLDEK` and `SAMPLLER` | first file: exit 0 (a warning is not an error), output has both rows, warn `duplicate_spectrum_row`. Second file: exit 1, `duplicate_spectrum_id` |
| 8 | Two inputs, `WT_2.ssl` then `WT_1.ssl` (`file` `/x/mzml/WT_2.mzML` and `/x/mzml/WT_1.mzML`), one row each, both `PEPTLDEK` charge 2 scan 0 | two rows sorted by run: `PEPTLDEK … WT_1 2 cascadia:0 denovo` then `… WT_2 …`. `spectrum_id` repeats across runs, which is legal (unique within `run`) |
| 9 | Two inputs whose `file` values are `/a/WT_1.mzML` and `/b/WT_1.mzML`; then the same `WT_1.ssl` given twice | both: exit 1, `run_in_several_inputs` naming the two paths |
| 10 | Empty-input cases: (a) a header-only `WT_1.ssl`; (b) (a) plus `WT_2.ssl` with one row; (c) a zero-byte file; (d) a file without the `score` column | (a) exit 1, `empty_table` (plus info `empty_input`); (b) exit 0, one `WT_2` row, info `empty_input` naming `WT_1.ssl`; (c) exit 1, `unreadable` "empty TSV (no header)"; (d) exit 1, `missing_column:score` |
| 11 | Charge `0`; charge `2.5`; charge `2.0`; `PEPK` (4 aa, score 0.9); a 41-mer | `0` → exit 1 (`out_of_range:charge`); `2.5` → exit 1 (`bad_type:charge`). On their own, `2.0` → `2`; `PEPK` and the 41-mer → exit 0, rows written, warnings `short_peptide`/`long_peptide` (preflight drops them later) |
| 12 | CLI: `gippa convert --from cascadia WT_1.ssl --out p.tsv` with `--source database`; `--search-db same`; `--abundance-out a.tsv`; `--source denovo`; and the example 1 input as `WT_1.ssl.gz` with CRLF line ends and a BOM | exit 1 for the first three (stderr names the flag; `search_db_invalid`; `no_abundance`), nothing written; `--source denovo` → exit 0 and the example 1 table; the gzip/CRLF/BOM file → byte-identical output to example 1. A second run of any exit-0 case → byte-identical output |
| 13 | A file without the three RT columns; a file with an extra column `ion-mobility` | both convert; the second gives warn `unknown_column` naming `ion-mobility`; RT values are never read, so `retention-time = abc` converts too |

## Edge cases

- **Modifications** (W5): every mass-shift and name form strips to the residues (example 2).
  `strip_modifications` deletes bracket groups first, then bare signed masses, then every
  non-`A–Z` character (`-` included). A bare `+15.995` outside brackets (not 0.0.7 output) strips
  too. M10-08's `evaluate` uses its own stripping rule. On every example 2 sequence the two must
  give the same peptide; criterion 4 checks this.
- **I/L:** never folded or unfolded by the converter (example 3). §5.3 aggregation groups
  `PEPTIDEK` and `PEPTLDEK` as two peptides, and the matcher's `i_equals_l` treats them as one
  query sequence. Real 0.0.7 output has no `I`. The contract test asserts this (C6) and treats an
  `I` as a finding: either the tokenizer setting changed, or the file is not 0.0.7 output.
- **Charge:** integers ≥ 1. One prediction is written per (precursor, charge), so one peptide
  appears on several rows with different charges and `spectrum_id`s. They aggregate in preflight
  (§5.3), not here.
- **Empty input:** example 10. A header-only `.ssl` is what Cascadia writes when nothing passes
  0.8. On its own it is an `empty_table` error, matching M10-08, which treats it as a "fail"
  verdict, so M10-10 never converts it. Among several inputs it is an info note.
- **Duplicate scans:** impossible in one 0.0.7 file (W3). Example 7 shows what the validator does
  if they occur. Across files, duplicates can only come from the same run twice, which step 3
  stops first (example 9).
- **Scores at exactly 0.8:** kept and counted (example 4). A default 0.0.7 run cannot produce one
  (strict `>`). After `%.6g`, values in (0.8, 0.8000005) read back as `0.8`. That affects only
  readers that re-apply 0.8 to the written table; GIPPA does not (D4). A count of such rows is
  recorded on FarmShare by the contract test (C4), not committed.
- **NaN and non-finite scores:** errors (example 5, D5). `NA`/empty cells become missing in
  `read_tsv`, so they are `bad_type:score` too.
- **Ties:** in sort keys, the stable sort keeps argument then file order (example 3 has distinct
  keys; two rows with the same run, peptide, modified form and charge, but different scans, keep
  file order). In score: the converter does not rank. Ties matter only for §7.2.5's decile edges,
  which merge on ties (`score_bins`).
- **RT unit:** never read (D6, example 13). The 539.9-s gradient next to the 500 heuristic (W8)
  therefore cannot change any output.
- **Run names:** the A3 run labels are the mzML stems (D1). M10-08 names the mzML `<RUN>.mzML`,
  so they are `WT_1`…`WT_3`. If those differ from DIA-NN's `Run` values, nothing in GIPPA joins on
  run across arms (`gippa compare` joins on locus tags). M10-10 records both lists.
- **Paths with odd characters:** a `file` value with a Windows `\` separator is not split on it;
  0.0.7 runs on Linux only. A `file` with a tab cannot occur (it would be a ragged row).
- **Large files:** a full run is ≈ 10⁴–10⁵ rows (M10-08 load estimate ≈ 143 k augmented spectra
  before the 0.8 filter). `read_tsv` holds the whole text in memory: ≈ 100 B per row, so well under
  1 GB. The A3 conversion job is M10-10's (2 CPUs, 8 GB, as M10-03's).
- **Circular vs linear, minus strand, coordinates:** not applicable. The converter handles
  sequences and identifiers only. Topology and strand enter at mapping (§7.2.1, M3), which this
  issue does not touch.

## What is committed (reconciling the index wording)

The index row says "contract test on a 1,000-row sample under `tests/data/recorded/cascadia/0.0.7/`".
Under Khoa's 2026-10-23 rule (row-level search output is never committed), that splits into three
parts:

| Item | Where | Committed |
|---|---|---|
| The 1,000-row sample `WT_1.ssl.gz`: the header plus the first 1,000 data rows of `WT_1.ssl` in file order (index decision 6), `gzip -n` | `$GIPPA_ROOT/pilot/recorded/cascadia/0.0.7/` (the `tests/data/recorded/` layout mirrored under `$GIPPA_ROOT/pilot/recorded/`; M10-08 deliverable 13 writes it) | **no** |
| `PROVENANCE.md`: Cascadia version, env lock sha256, weights sha256, route, job IDs, input archive sha256, the W1 header from the source (M10-08). **M10-09 appends a section** with the sample's sha256 and data-row count (1,000, or fewer if the run had fewer), the selection rule, the gippa commit the contract test ran at, and its result (passed/skipped/failed per check C1–C8, no values) | `tests/data/recorded/cascadia/0.0.7/PROVENANCE.md` | **yes**: the only file in that directory |
| The converter, the CLI branch, CI tests whose `.ssl` fixtures are **written in `tmp_path` from text in the test file**, and the FarmShare contract test | `src/`, `tests/unit/`, `tests/contract/` | yes |
| Converted A3 tables (`a3.tsv.gz`) | `$GIPPA_ROOT/pilot/inputs/` (M10-10) | **no** (as M10-03 Q-K1's default for DIA-NN) |

No `.ssl` file is committed anywhere, not even a synthetic one. M10-08 criterion 12 already fails
on any tracked `.ssl`, so CI fixtures are built at test time.

If M10-08 did not write the sample (for example, its 08b landed without that step), the lead
makes it in a `dev` job (1 CPU, 2 GB, 10 min):
`head -n 1001 "$GIPPA_ROOT/pilot/A3_cascadia/WT_1/WT_1.ssl" | gzip -n > "$GIPPA_ROOT/pilot/recorded/cascadia/0.0.7/WT_1.ssl.gz"`.
The job ID goes into PROGRESS.

## The contract test (`tests/contract/test_m10_09_cascadia_contract.py`, test-author)

- **Data:** `SAMPLE = $GIPPA_ROOT/pilot/recorded/cascadia/0.0.7/WT_1.ssl.gz`. Optional, `slow`:
  `FULL = $GIPPA_ROOT/pilot/A3_cascadia/WT_1/WT_1.ssl`.
- **Skip, never fail, when absent:** `GIPPA_ROOT` unset (CI) → skip "`GIPPA_ROOT` is not set: the
  Cascadia sample lives on FarmShare only (Khoa 2026-10-23)"; file missing → skip naming the path
  and "M10-08 smoke run not done or verdict fail". The test asserts that the sample path is
  outside the checkout. No module-level read of `$GIPPA_ROOT`. No `pending_sample` marker (the
  sample is not expected in the repository); the M10-03 pattern.
- **Checks on the raw sample** (each its own test function; the raw file is read with the stdlib
  `csv` module, independently of the converter):
  - **C1 header:** exactly the nine W1 names, in W1 order.
  - **C2 score-type:** `UNKNOWN` on every row.
  - **C3 score:** finite floats with 0.8 < s ≤ 1. Default threshold: `denovo_pilot.sbatch` passes
    no `-t` (M10-08 deliverable 8).
  - **C4 scan:** non-negative integers, strictly increasing in file order, so unique. The number of
    scores in (0.8, 0.8000005) is printed for the PROVENANCE note, not asserted.
  - **C5 file:** one distinct value; its D1 stem is `WT_1`.
  - **C6 sequence:** `strip_modifications(sequence)` fullmatches `[A-Z]+` and contains no `I`.
    Every bracket token belongs to the W5 set {`[+15.99]`, `[+0.98]`, `[+57.02]`, `[+42.01]`,
    `[+43.01]`, `[-17.03]`, `[+25.980265]`}, and `-` appears only right after a leading token. A
    name token such as `[Oxidation]` fails, naming M10-08's offline Unimod launcher (its edge case
    "modified sequences without mass shifts").
  - **C7 charge:** integers ≥ 1.
- **Checks on the conversion:**
  - **C8:** `convert_cascadia([SAMPLE], search_db=None, cfg=InputConfig())`:
    - no error issue;
    - one output row per sample row; the multiset of (`run`, `spectrum_id`, `modified_peptide`,
      `charge`, `score`) equals the raw (stem, `cascadia:`+scan, sequence, charge, float(score));
    - every `source` is `denovo`; no `q_value` or `pep` column.

    The CLI on the sample (`--out t/a3.tsv.gz`) exits 0, a second call is byte-identical, and
    `read_peptide_table` + `validate_peptides` on the written file give no error issue.
- **No golden output file** is read or written (it would be derived data and would test the
  converter against itself).
- **`slow` variant:** C1–C8 on `FULL`, skipped when absent. It runs on FarmShare in M10-10's
  conversion job; it is not in the fast suite.

## Acceptance criteria

On CI (every commit; no Cascadia, no real data, no `$GIPPA_ROOT`):

1. `CONVERTERS["cascadia"] is convert_cascadia`; `gippa convert --from cascadia` no longer exits 5;
   `--from casanovo` and `--from instanovo` still exit 5.
2. Examples 1–13 hold exactly: the output bytes for exit-0 cases, the exit codes, and the issue
   codes with severities and 1-based row lists. Fixtures are written in `tmp_path`.
3. Determinism: each exit-0 example gives byte-identical files on two runs; `.gz` output is
   byte-identical across runs; inputs given in a different argument order give the same bytes
   when runs are distinct.
4. The converter's `peptide` equals M10-08's stripping rule (strip, then fold I → L) on every
   example 2 sequence once both are folded. If `validation/pilot/a3_feasibility.py` is not merged
   yet, the test encodes M10-08's four-step rule inline and cites the spec.
5. The output table has no `q_value`, `pep`, `protein_ids`, `is_decoy`, `intensity` or
   `search_db` column, and its columns are in §5.3 order. Read back with `read_peptide_table` and
   passed to `steps.mapping.resolve_sources`, every row is `denovo`, even with `source` dropped
   from the frame (no q/pep present).
6. No `spectrum_id` of the output equals its bare `scan` value; each starts with `cascadia:`.
7. `src/gippa/io/converters/cascadia.py` imports nothing from Cascadia, pyteomics or torch
   (`scripts/ci/check_imports.py` passes). It reads no `retention-time`, `start-time` or
   `end-time` value (example 13 proves it).
8. The existing tests in `tests/unit/test_cli_convert.py`, `test_convert_generic.py` and
   `test_convert_diann.py` pass unchanged, with one exception. The `"cascadia"` case of
   `tests/unit/test_cli_convert.py::test_other_engines_not_implemented` (which expects exit 5) is
   removed from its parametrize list by the lead in the implementation step (role guard: the
   test-author does not edit existing tests), citing this spec. `casanovo` and `instanovo` stay.
9. Guard: `git ls-files` lists no `*.ssl` or `*.ssl.gz`, and `tests/data/recorded/cascadia/` holds
   only `0.0.7/PROVENANCE.md`.
10. The contract test collects on CI and skips with the reasons above, never errors.
11. Lint and types: `ruff check`, `ruff format --check`, `mypy src/gippa` clean.

On FarmShare (evidence: PR description, PROGRESS, PROVENANCE; not CI):

12. With `GIPPA_ROOT` set, `python -m pytest -p no:cacheprovider
    tests/contract/test_m10_09_cascadia_contract.py` (venv Python, `PYTHONPATH` on the worktree)
    passes C1–C8 on the real sample. Its output is pasted into the PR, without any row.
13. Any failed check becomes a converter change with its own CI test on a hand-made fixture that
    shows the condition and never copies rows; or, for C3/C6 (wrong threshold, unshifted names), a
    finding for M10-08.
14. `PROVENANCE.md` has the M10-09 section (sample sha256, row count, commit, per-check result).

## Test plan

Owner: the test-author, tests first (R2), in new files only.

- `tests/unit/test_m10_09_convert_cascadia.py`: criteria 2–7 for the function (examples 1–11 and
  13), with a small `write_ssl(tmp_path, rows, header=W1)` helper.
- `tests/unit/test_m10_09_cli_cascadia.py`: criterion 1 and example 12 through `CliRunner`
  (gzip, CRLF and BOM; the flag errors; byte identity).
- `tests/unit/test_m10_09_guard.py`: criterion 9.
- `tests/contract/test_m10_09_cascadia_contract.py`: criteria 10 and 12; the `slow` full-file
  variant.

Local runs: targeted files only, with the venv Python and `PYTHONPATH=$PWD/src:$PWD/tests/sim:$PWD/tests`,
`-n 2` at most, never `uv run`. All four files run in seconds. The full suite goes to Slurm
(`test_branch.sbatch`). No new sbatch script; no GPU.

## Estimate

| Piece | Code lines |
|---|---|
| `cascadia.py` (read, run/scan/score checks, mapping, issues) | ≈ 125 |
| `CONVERTERS` / `__all__` | 2 |
| CLI branch (D7) | ≈ 10 |
| **Total** | **≈ 140** (index: 200; there is no sbatch or summary script, since the conversion job is M10-10's) |

Tests ≈ 330 lines. One PR. Agent time ≈ ½ day plus the FarmShare check after M10-08's smoke run.

## Statistician notes

This issue makes one statistical claim and changes no estimator: **Cascadia's `score` is not
treated as a calibrated probability**, so no `q_value` or `pep` is derived from it.

- The score is `exp(mean token log-probability)` of an autoregressive decoder trained on Orbitrap
  DIA. It has not been calibrated on any data set and is out of domain on timsTOF. Using
  `1 − score` as a PEP would put an uncalibrated number into π̂_wrong for m = 0 (§7.2.5 step 4
  uses mean PEP only for database input) and into the §7.3.1 engine threshold. §7.2.5 already
  calibrates de novo rows from the genome (target versus shuffled-decoy exact matches per score
  decile, with `c_ℓ` from A2 and a bootstrap CI). That method needs only a ranking.
- Decile binning (`mapping.denovo_precision.score_bins`, quantiles of `score`) is invariant to
  any strictly increasing transform of the score. Mapping the raw `score` loses nothing against,
  say, its logarithm.
- `%.6g` keeps 6 significant digits. Rows whose scores differ only beyond that become ties after
  writing, and tied scores share a bin (`score_bins` merges coinciding edges). With ≥ 10³ distinct
  scores in (0.8, 1] the effect on decile edges is ≤ 5 × 10⁻⁷, far below any bin width.
- The truncation at 0.8 means A3's deciles span (0.8, 1] only. The precision estimate is
  conditional on Cascadia's own cut-off. M10-10's report must say so, and must not compare the
  bins with those of another engine.

No `stats` label: the converter computes no statistic and sets no threshold (the 0.8 is
Cascadia's CLI default, recorded, never applied by GIPPA).

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| Q1 | `run`: the stem of `file` (D1), or the value verbatim (§5.4's "as Casanovo", whose `ms_run` location is a path)? | stem; the verbatim path is an absolute FarmShare path (§6.6) |
| Q2 | `spectrum_id`: `cascadia:<scan>` (D2), `scan` bare, or no `spectrum_id` at all? | `cascadia:<scan>`: unique within run, no false pairing with scan-numbered references |
| Q3 | Score → q: write only `score` (no `q_value`, no `pep`)? | yes (D3, statistician notes) |
| Q4 | Rows with `score ≤ 0.8` (only possible with a lowered `-t` or a hand-made file): keep them with an info count, drop them, or fail? | keep, with info `score_at_or_below_default_cutoff` |
| Q5 | NaN, missing, non-finite or out-of-[0, 1] scores: an error (exit 1), or warn and drop the rows? | error: 0.0.7 cannot write them, so the file is not Cascadia output |
| Q6 | `score-type` other than `UNKNOWN` (another tool's `.ssl`, e.g. a q-value where lower is better): error? | error `unsupported_score_type` |
| Q7 | With `--from cascadia`: accept `--source denovo`, reject `--source database` and `--search-db` (exit 1)? | yes (D7); `--source denovo` was wrongly rejected for every non-generic engine |
| Q8 | The converted A3 table is "derived" like the DIA-NN tables (M10-03 Q-K1 default): keep it under `$GIPPA_ROOT/pilot/inputs/`, uncommitted? | yes; M10-10 converts and records the sha256 and row count. Not asked of Khoa separately: the answer to Q-K1 covers both |
| Q9 | Who writes the 1,000-row sample: M10-08 (deliverable 13), with the lead's `head \| gzip -n` dev job as the fallback? | M10-08; fallback as stated |
| Q10 | `PROVENANCE.md` gets an M10-09 section with the sample sha256, the row count and per-check pass/fail. Is the count of `%.6g` boundary scores (C4) also allowed there, or does it stay in PROGRESS? | sha256, row count and pass/fail in PROVENANCE; the boundary count only in the PR text and PROGRESS |
| Q11 | RT columns: ignore entirely (D6), or validate them as numeric, which would make a corrupt RT fail an otherwise good file? | ignore |
| Q12 | The index estimate was 200 lines; this design is ≈ 140, because the A3 conversion job and its summary belong to M10-10. Confirm that split? | yes; M10-10 reuses `diann_convert.sbatch`'s pattern for `cascadia_convert.sbatch` |
| Q13 | `scan` → `spectrum_id` relies on W3 (unique per file). If a later Cascadia version reuses indices across output files of one run, the `run_in_several_inputs` rule already stops a merge. Is any other cross-version protection wanted (a `cascadia_version` check)? | no: the `.ssl` carries no version; the version is pinned in the env (M10-08) and recorded in PROVENANCE |

## Lead decisions (2026-10-24)

Q1–Q13 all take their defaults. The lead removes the `"cascadia"` case of `tests/unit/test_cli_convert.py::test_other_engines_not_implemented` during implementation (criterion 8, role guard).
