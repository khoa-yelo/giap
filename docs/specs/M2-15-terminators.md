# M2-15: Terminators: TerminatorNet, TransTermHP, cross-check, BacTermFinder (opt-in), `terminator_between`

Milestone M2 · depends on: M2-04 · estimate ≈ 390 code lines (tests not counted) · capabilities: none
(TerminatorNet fixture comes from M2-10; TransTermHP fixture recorded with the harness, `needs_conda`)

## Spec sections implemented

- §7.1.1 TerminatorNet: "`python TerminatorNet.py -f genome.fa -s -n N -o terms.tsv` | intrinsic
  terminators, confidence `pred`, hairpin ΔG". Gotchas: "Its argument parser matches flags by
  substring. Pass only the documented flags, use staged file names without spaces, and always pass
  `-s`." "Output lines are `name, TERM_i, A_tract, hairpin, U_tract, structure, GH, pred, start,
  stop[, strand]`, 1-based, with minus-strand coordinates already converted to the forward strand."
  "Skip lines reading `NO TERMINATORS FOUND`." TransTermHP: "`transterm -p expterm.dat genome.fa
  genes.coords` | terminators with confidence". BacTermFinder: "per its README | extra terminator layer, source-tagged".
- §7.1.5: "**TerminatorNet** (default) runs on the whole genome with `-s`." "**Terminator between two
  same-strand adjacent genes:** a terminator on that strand whose interval lies within [upstream gene
  stop − 50 bp, downstream gene start], in transcription direction, with confidence ≥
  `terminatornet_yes`." "**Calibration.** Until RegulonDB calibration (§15.4) the cut-off is a
  placeholder flagged `uncalibrated` in the manifest." "**TransTermHP cross-check.** Report the share
  of TerminatorNet calls within ±20 bp of a TransTermHP call in `annotation/terminators_crosscheck.tsv`."
  "**`--terminators transtermhp`** makes TransTermHP the main tool, with confidence = CONF/100."
  "**BacTermFinder calls** (opt-in) are added to `terminators.bed` with `source=BacTermFinder` but are
  not used in Module 4 flags in version 1."
- G12; §5.2 `terminators: {tool, cross_check_transtermhp: true, bactermfinder: false, terminatornet_yes: 0.5}`;
  §9.3 `genes.tsv` `terminator_conf_downstream`; Module 1 test "a TerminatorNet `-s` and coordinate test".

## Public interface

`src/giap/annotate/terminators.py`, `terminatornet_cmd.py`, `transtermhp_cmd.py`, `bactermfinder_cmd.py`;
steps `transterm_input`, `terminator_layer`.
```python
Source = Literal["TerminatorNet", "TransTermHP", "BacTermFinder"]
@dataclass(frozen=True, slots=True)
class Terminator: interval: Interval; confidence: float; source: Source; hairpin_dg: float | None = None

def terminatornet_command(genome_fa: Path, out_tsv: Path, threads: int, home: Path) -> list[str]
def parse_terminatornet(path: Path, replicons: Mapping[str, Replicon]) -> list[Terminator]
def write_transterm_coords(path: Path, genes: Sequence[Gene]) -> None
def transterm_command(genome_fa: Path, coords: Path, expterm: Path) -> list[str]
def parse_transterm_tt(path: Path, replicons) -> list[Terminator]          # confidence = CONF / 100
def parse_bactermfinder(path: Path, replicons) -> list[Terminator]         # format fixed by the recorded fixture
def crosscheck(primary: Sequence[Terminator], other: Sequence[Terminator], tol_bp=20) -> list[CrossRow]
def terminator_between(up: Gene, down: Gene, terms, replicon, slack_bp=50, yes=0.5) -> float | None
def terminator_conf_downstream(gene: Gene, next_gene: Gene | None, terms, replicon, slack_bp=50) -> float | None
def terminators_bed(terms, replicons) -> list[BedRecord]
```
Rules:
- TerminatorNet: genome staged as `genome.fa` (no spaces); command `python <home>/TerminatorNet.py
  -f genome.fa -s -n N -o terms.tsv` — exactly these flags. Parser: skip blank, header and
  `NO TERMINATORS FOUND` lines; `name` = contig; `start`, `stop` 1-based forward → `Interval(contig,
  start−1, stop, strand)`; `confidence = pred`; `hairpin_dg = GH`; strand from the 11th column; if
  missing (no `-s`), raise (we always pass `-s`).
- `.coords`: one line per CDS/sORF/pseudogene `"<locus_tag> <start1> <end1> <contig>"`, minus-strand
  genes written with start > end (TransTermHP convention).
- `.tt` lines `TERM <n> <a> - <b> <strand> <loc> <conf> <hp> <tail> | …` → interval
  `[min(a,b)−1, max(a,b))`, confidence `conf/100`; the current contig from `SEQUENCE <name>` lines.
- `crosscheck`: a primary call matches when an `other` call on the same contig and strand lies within
  `tol_bp` (gap between intervals ≤ 20; 0 if overlapping). Output rows per contig + `all`.
- Between window in transcription direction: `+` pair (up `[a,E)`, down `[S,…)`): terminator inside
  `[E − 50, S)`; `−` pair (up = higher coordinates `[s_up, e_up)`, down `[…, e_down)`): inside
  `[e_down, s_up + 50)`. Returns the max confidence ≥ `yes`, else `None`. Only primary-tool calls
  (not BacTermFinder) count.
- `terminator_conf_downstream(g)`: max confidence of primary calls in the window after `g` up to the
  next same-strand gene start (or 300 bp if none, Q3), no threshold.
- `terminators.bed`: BED6, name `<source>_<n>` (e.g. `TerminatorNet_1`), score `round(1000 × conf)`.

`annotation/terminators_crosscheck.tsv` columns (Q2): `contig`, `n_terminatornet`, `n_transtermhp`,
`n_matched`, `share_matched`; last row `contig = all`.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `terminatornet_command(Path("genome.fa"), Path("terms.tsv"), 8, H)` | `["python", f"{H}/TerminatorNet.py", "-f", "genome.fa", "-s", "-n", "8", "-o", "terms.tsv"]` |
| 2 | line `chr1,TERM_3,AAAAAA,GCCGCC…,TTTTTTT,((((…)))),-12.4,0.87,1201,1240,-` | `Terminator(Interval("chr1",1200,1240,"-"), 0.87, "TerminatorNet", -12.4)` |
| 3 | file with only `NO TERMINATORS FOUND` | `[]` |
| 4 | genes `GSIM_00005 (99,400,+)`, `GSIM_00010 (1000,1300,−)` | `.coords` lines `GSIM_00005 100 400 chr1` and `GSIM_00010 1300 1001 chr1` |
| 5 | `.tt` line `  TERM 1   4342 - 4366   +   F   93 -11.5 -3.8 | …` under `SEQUENCE chr1` | `Interval("chr1",4341,4366,"+")`, confidence 0.93 |
| 6 | primary `(1000,1040,+)`; TransTermHP `(1060,1090,+)` (gap 20) / `(1061,1090,+)` (gap 21) | matched / not matched |
| 7 | + pair up `(100,400)`, down `(420,700)`; terminator `(360,400,+)` conf 0.8 / `(345,400,+)` / `(360,400,−)` | `0.8` / `None` (starts before 350) / `None` (wrong strand) |
| 8 | same with conf 0.4, `yes = 0.5` | `None` |
| 9 | `--terminators transtermhp` | primary = TransTermHP calls; cross-check table not written |
| 10 | BacTermFinder call | present in `terminators.bed` as `BacTermFinder_<n>`; ignored by `terminator_between` |

## Edge cases

- Terminators crossing the origin (circular) keep `end > length`; BED splits them.
- Overlapping adjacent genes (`S < E − 50`): empty window → `None`.
- Contig names with characters TerminatorNet may mangle cannot occur (safe IDs, M1-02); a test
  feeds a contig `chr_1.a-b` through a recorded-format line.
- `cross_check_transtermhp: false` → TransTermHP not run, cross-check file absent (manifest notes it).
- Manifest flag `uncalibrated` while `calibration.yaml` says `terminatornet_yes.calibrated: false`.

## Acceptance criteria

- [ ] `tests/unit/test_terminators.py`: examples 1–10.
- [ ] `tests/contract/test_terminatornet_contract.py` (sim100k, k12 from M2-10) and
  `test_transtermhp_contract.py`, `test_bactermfinder_contract.py`; `pending_sample` until recorded.
- [ ] `test_nf_cmd_sync.py` covers TerminatorNet (`-s` present) and TransTermHP.

## Questions for the lead

1. Where does `expterm.dat` come from? **Proposed default:** the file shipped in the bioconda
   `transtermhp=2.09` package (`$CONDA_PREFIX/data/expterm.dat`), path resolved in the env.
2. `terminators_crosscheck.tsv` columns are not in §9.3. **Proposed default:** the 5 columns above.
3. `terminator_conf_downstream` window when no next same-strand gene. **Proposed default:** 300 bp
   past the stop (or the contig end on linear contigs).

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions 1–5 are accepted (tab-separated TerminatorNet without header; two-letter
   `.tt` location codes and integer hp scores; terminators.parquet = primary-tool calls (+ BacTermFinder)
   only; step keys matched by name; cross-check file may be header-only; per-source BED numbering; windows
   wrap at the origin on circular contigs).
2. The manifest `uncalibrated` flag and step-level `cross_check_transtermhp: false` are tested when manifest
   and optional-output support land (backlog; operon_meta → manifest wiring).

## Lead decisions after code review (PR #75, changes requested)
3. nf-sync test accepts a variable after `-p` (e.g. `$EXPTERM`) and checks the positional order
   `genome.fa genes.coords` (it must pass against #65's TRANSTERMHP script).
4. The TerminatorNet parser never skips data lines by prefix: only an exact header line (if any) is
   skipped; contig IDs starting with "name" parse normally. Test with `name_contig`.
5. BacTermFinder (provisional): an empty output directory → no calls (`[]`); exactly one CSV → parse it;
   several CSVs → ValueError naming them (no arbitrary pick). Stays provisional until a recording exists.
6. `terminator_between(up, down)` raises ValueError unless both genes share contig and strand.
7. Tests kill the surviving mutants: `−L` shift in `_max_in` at the origin; downstream window slack with
   no next gene; BED score rounding (`round`, not `int`); NaN `share_matched` for a contig without
   primary calls; cross-check contig order = FASTA order; the `bactermfinder` gate;
   `cross_check_transtermhp: false`; layer sorted.
