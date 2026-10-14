# M1-09: Preflight framework and genome rules (part A)

Milestone M1 · depends on: M1-01, M1-02 · estimate ≈ 400 changed lines (≈ 170 code, ≈ 230 tests)

## Spec sections implemented

- §7.0 purpose: "Validate inputs and environment in seconds … Errors stop the run (exit 1). Warnings
  become flags carried into every table and the report." · Tests: "`T-PRE-*` covers every rule with one
  passing and one failing example." · "A complete genome with two small plasmids raises no contiguity
  warning."
- §7.0 genome rows (quoted in the table below), §2.1 decision 12 ("Warn on N50 < 100 kb or largest
  contig < 1 Mb. The median contig length is recorded for information only."), §5.2 `preflight.*`.
- §7.1.2: "`replicons.tsv` for Bakta is built by preflight: the largest contig is the chromosome if it
  is circular or ≥ 1 Mb; other circular contigs < 500 kb are plasmids; everything else is a contig."
- §5.3: contigs are circular when tagged, or all of them when `annotation.bakta.complete: true`.

## Public interface

`src/giap/preflight/checks.py`
```python
Level = Literal["ok", "info", "warn", "error"]
InputName = Literal["inputs", "genome", "proteins", "peptides", "abundance", "environment"]

@dataclass(frozen=True, slots=True)
class Check:
    rule_id: str            # stable, e.g. "genome.contig_count"; T-PRE-n in tests maps 1:1
    input: InputName
    level: Level            # "ok" when the rule passed
    rule: str               # human rule with the threshold, e.g. "contig count > 10 → warn"
    value: JsonValue        # the measured value (number, string, list or object)
    message: str = ""
    exit_code: int = 0      # 2 only for environment errors (M1-11)

def worst(checks: Iterable[Check]) -> Level: ...
```
`src/giap/preflight/genome.py` (pure; no file I/O)
```python
@dataclass(frozen=True, slots=True)
class GenomeSummary:
    ids: tuple[str, ...]; original_ids: tuple[str, ...]; lengths: tuple[int, ...]
    circular_tags: tuple[bool, ...]; gc: int; acgt: int
    n_runs: tuple[tuple[str, int, int], ...]        # (contig, start0, end0), length ≥ n_run_min_bp
    ambiguous: tuple[tuple[str, int, str], ...]     # (contig, pos0, code) for IUPAC ≠ ACGTN
    invalid: Mapping[str, int]; empty_records: tuple[str, ...]

def summarize_genome(records: Sequence[GenomeRecord], n_run_min_bp: int) -> GenomeSummary: ...
def n50(lengths: Sequence[int]) -> int: ...
def resolve_topology(s: GenomeSummary, complete: bool | Literal["auto"]) -> tuple[Topology, ...]: ...
def classify_replicons(s: GenomeSummary, topo: Sequence[Topology], small_bp: int,
                       chromosome_min_bp: int = 1_000_000) -> list[Replicon]: ...
def genome_checks(s: GenomeSummary, cfg: GiapConfig, parse_error: str | None = None) -> list[Check]: ...
```
`n50`: sort descending; the length at which the cumulative sum first reaches ≥ half the total.
GC = (G + C) / (A + C + G + T), ambiguous codes excluded. Positions in `value` are 1-based closed
(JSON is human-facing, like TSV); lists are capped at 1,000 entries with a total count.
`classify_replicons` kinds: `chromosome`, `plasmid`, else `unknown` (Bakta type `contig`). Ties for
"largest" go to the first contig in FASTA order.

## Rules, with one passing and one failing example each

Tests use a `PreflightConfig` with small thresholds where noted, so sequences stay tiny.

| Test | `rule_id` | Rule (spec) | Level if failing | Passing example | Failing example |
|---|---|---|---|---|---|
| T-PRE-01 | `genome.parse` | "≥ 1 record, IUPAC nucleotides only" | error | `>c\nACGTN\n` | empty file; and `>c\nACGU\n` (value `{"U": 1}`); and record with empty sequence |
| T-PRE-02 | `genome.ids_unique` | "unique first tokens" (after safe mapping) | error | `>a`, `>b` | `>a x`, `>a y`; and `>a\|b`, `>a:b` (both map to `a_b`) |
| T-PRE-03 | `genome.ids_renamed` | "unsafe characters mapped …, map saved" | warn | `>NZ_CP1.1` | `>chr:1` → value `[{"from":"chr:1","to":"chr_1"}]` |
| T-PRE-04 | `genome.contig_count` | "> 10" | warn | 10 contigs | 11 contigs (value 11) |
| T-PRE-05 | `genome.total_length` | "< 100 kb" | warn | total 100,000 bp | total 99,999 bp |
| T-PRE-06 | `genome.contiguity` | "N50 < 100 kb or largest contig < 1 Mb; median length recorded" | warn | lengths `[4_600_000, 60_000, 5_000]` (complete genome + 2 small plasmids): ok, value `{"n50":4600000,"largest":4600000,"median":60000}` | lengths `[900_000, 800_000]`: warn (largest < 1 Mb); lengths `[1_200_000] + [50_000]*30`: warn (N50 50,000) |
| T-PRE-07 | `genome.gaps` | "runs of ≥ 10 N and other IUPAC codes, with positions" | warn | `ACGTNNNNNNNNNACGT` (9 N) | `ACG` + 10 N + `T` → `[{"contig":"c","start":4,"end":13}]`; `ACRT` → `[{"contig":"c","pos":3,"code":"R"}]` |
| T-PRE-08 | `genome.gc` | "outside 25–75 %" | warn | `ACGT` (0.5); `AAAC` (0.25, boundary) | `AAAAAAAC` (0.125); `GGGGGGGA` (0.875) |
| T-PRE-09 | `genome.topology` | "`circular=true` or `[topology=circular]` → Bakta replicon table" | info | always `info`; value = replicon table | ok-case: `>c circular=true` → `circular`; `complete: true` → all circular; untagged → `linear` |
| T-PRE-10 | `genome.small_replicons` | "contigs < 500 kb listed as putative plasmids or fragments" | info | one contig of 600 kb → `ok`, value `[]` | contigs 4.6 Mb + 60 kb → `info`, value `["p1"]` |

Replicon classification examples (`small_bp=500_000`):
- `chr` 4.6 Mb linear, `p1` 60 kb circular → `chromosome`, `plasmid`.
- `c1` 900 kb circular, `c2` 50 kb linear → `chromosome` (largest and circular), `unknown`.
- `c1` 900 kb linear, `c2` 600 kb circular → `unknown`, `unknown` (not largest; ≥ 500 kb).

## Edge cases

- When `genome.parse` is an error, the other genome rules are still evaluated when the summary
  exists (to report everything at once), except if there are zero records.
- `n50([])` raises `ValueError`; `n50([5]) == 5`; `n50([3, 3, 2]) == 3`.
- Lower-case `n` counts as `N` (reader upper-cases). An N run crossing a line break is one run.
- An N run at the origin of a circular contig is reported as two runs (no wrap in reporting).

## Acceptance criteria

- [ ] `tests/unit/preflight/test_genome_rules.py` has `test_T_PRE_01_pass`, `test_T_PRE_01_fail`, …
  through `T_PRE_10`, plus the three classification examples and the `n50` edge cases.
- [ ] Every `Check` produced has a `rule_id` from the table; unknown IDs fail a test.
- [ ] Rules are pure functions of `GenomeSummary` + config (no file access; checked by importing
  only `giap.preflight.genome` in tests with in-memory records).

## Questions for the lead

1. Should passing rules appear in `preflight.json`? **Proposed default:** yes, with `level: "ok"`, so
   the report shows every rule evaluated and tests can assert the pass case.
2. The chromosome threshold "≥ 1 Mb" (§7.1.2) has no config key. **Proposed default:** reuse
   `preflight.min_largest_contig_warn` (1,000,000) instead of a new key, avoiding a `defaults.yaml` ADR.
3. §7.0 says warnings "become flags carried into every table"; §9.4 has no preflight flag names.
   **Proposed default:** M1 records warnings in `preflight.json`; the flag vocabulary is extended (via
   the schema) when M4 builds `genes.tsv`.
