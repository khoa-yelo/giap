# M2-14: Promoter Calculator: one region per call, parser, `tx_rates`, TU transcription features

Milestone M2 · depends on: M2-12 (TUs, `write_tu_regions`) · estimate ≈ 200 code lines (tests not
counted) · capabilities: none (recording needs `needs_conda`)

## Spec sections implemented

- §7.1.1 Promoter strength: "Promoter Calculator 1.2.4 (PyPI `promotercalculator`, GPL) | own env |
  `promoter-calculator -i region.fa -o region.csv -j 1`, once per TU upstream region; regions run in
  parallel (G11) | σ70 Tx rate, TSS, strand". Gotcha: "It ignores FASTA headers and joins all lines
  into one sequence, so pass exactly one region per call."
- §7.1.4: "**Promoter Calculator.** Run once per TU upstream region. Keep promoters on the gene's
  strand whose TSS lies in the region. Features: `tx_rate_max`, `tx_rate_sum`, `tss_distance_bp`."
- §9.3 `genes.tsv`: `tx_rate_max_tu`, `tx_rate_sum_tu`; `pairs.tsv`: `tx_rate_between` (M4 uses the
  helper defined here). §5.2 `promoter_calculator: {enabled: true}`.
- Module 1 tests: "a Promoter Calculator one-region-per-call test".

## Public interface

`src/giap/annotate/promoter_calculator.py`, `src/giap/annotate/promoter_calculator_cmd.py`; step `tx_layer`.
```python
def region_fasta_text(tu: Tu, genome: Mapping[str, str]) -> str
    # ">{tu_id}\n{seq}\n": exactly one record, one sequence line, 5'→3' in TU orientation (revcomp for −)
def promoter_calculator_command(region_fa: Path, out_csv: Path) -> list[str]   # ["promoter-calculator","-i",…,"-o",…,"-j","1"]

@dataclass(frozen=True, slots=True)
class PcHit: tss: int; strand: str; tx_rate: float      # tss 1-based in region coordinates
def parse_pc_csv(path: Path) -> list[PcHit]

def tx_layer(tus: Sequence[Tu], hits: Mapping[str, Sequence[PcHit]], replicons) -> pd.DataFrame
    # columns: contig, strand, tss (0-based genome), tx_rate, tu_id — only hits with strand "+" in region
    # orientation (= TU strand) and TSS inside the region
def tu_tx_features(tus, tx: pd.DataFrame) -> pd.DataFrame
    # tu_id, tx_rate_max, tx_rate_sum, tss_distance_bp (distance from the max-rate TSS to the first start codon)
def tx_rate_between(up: Gene, down: Gene, tx: pd.DataFrame, replicon) -> float | None
    # max tx_rate of TSSs on the pair's strand lying between up's end and down's start codon (transcription order)
```
Region → genome: `+` TU with region `[a, b)`: genome TSS `= a + tss − 1`; `−` TU: `= b − tss`.
`tss_distance_bp` = start codon position − TSS in transcription direction (positive = upstream).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `region_fasta_text` for `TU_1` (+, region `[850, 1050)`) | `">TU_1\n" + genome[850:1050] + "\n"`; 200 nt on one line |
| 2 | same for a − TU, region `[5550, 5800)` | `">TU_3\n" + revcomp(genome[5550:5800]) + "\n"` |
| 3 | `promoter_calculator_command(r.fa, r.csv)` | `["promoter-calculator","-i","r.fa","-o","r.csv","-j","1"]` |
| 4 | TU_1 hit `tss=101, strand="+", tx_rate=1234.5` | `tx_layer` row `contig=chr1, strand=+, tss=950, tx_rate=1234.5, tu_id=TU_1`; first start codon S = 1000 → `tss_distance_bp = 50` |
| 5 | TU_3 (−) hit `tss=51, "+", 80.0` | genome TSS `5800 − 51 = 5749`; start codon first base 5599 → `tss_distance_bp = 150` |
| 6 | hit with strand `-` (opposite to the TU) | dropped |
| 7 | TU with hits 10.0, 30.0, 5.0 (all +) | `tx_rate_max = 30.0`, `tx_rate_sum = 45.0`, `tss_distance_bp` from the 30.0 hit |
| 8 | TU with no kept hits | `tx_rate_max = 0`, `tx_rate_sum = 0`, `tss_distance_bp = NA` (Q1) |
| 9 | a region file with two records | writer never produces it (assert); test: `region_fasta_text` output has exactly one `>` and two lines |

## Edge cases

- Regions crossing the origin are written as one contiguous sequence (wrapped slice).
- Regions shorter than the tool's minimum (e.g. heavily clipped) are still run; an empty CSV is a
  valid result with no hits.
- `promoter_calculator.enabled: false` → `tx_layer` not planned; `tx_rate_*` columns `NA` and run
  flag `layer_missing:tx_rates` (M2-20).
- Hits whose TSS lies outside `[1, len(region)]` are dropped.

## Acceptance criteria

- [ ] `tests/unit/test_promoter_calculator.py`: examples 1–9 with hand-written CSVs.
- [ ] `tests/contract/test_promoter_calculator_contract.py` on recorded CSVs for 5 sim100k TU regions;
  `pending_sample` until recorded.
- [ ] `test_nf_cmd_sync.py` extended with this tool (M2-07).

## Questions for the lead

1. Value for a TU without any kept promoter. **Proposed default:** `tx_rate_max = tx_rate_sum = 0`
   (no predicted σ70 promoter is information), `tss_distance_bp = NA`.
2. Promoter Calculator CSV column names are not in the spec. **Proposed default:** fixed by the
   recorded fixture; unit tests use a doc-derived header (`TSS`, `strand`, `Tx_rate`) until then.
