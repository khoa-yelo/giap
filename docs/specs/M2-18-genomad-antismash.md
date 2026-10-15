# M2-18: geNOMAD and antiSMASH: command builders and parsers

Milestone M2 · depends on: M2-03 (gene model, locus tags) · estimate ≈ 300 code lines (tests not
counted) · capabilities: none (fixtures from M2-10; live runs in M2-21 need `needs_genomad_db`, `needs_apptainer`)

## Spec sections implemented

- §7.1.1: "Prophages, plasmids | geNOMAD 1.12.0, DB v1.9 (non-commercial) | conda `genomad=1.12.0` |
  `genomad end-to-end --cleanup --splits 8 --threads N genome.fa out DB/genomad_db` | provirus
  coordinates, plasmid scores, virus gene annotations". "BGCs | antiSMASH 8.0.4 (AGPL) | container
  `antismash/standalone` … or conda `antismash=8.0.4` | `antismash --genefinding-tool none --minimal
  --cpus N --output-dir out bakta.gbff` | regions, products, core biosynthetic genes".
- §6.1 "geNOMAD [mobile] (parallel with Bakta)"; §6.2 `genomad` "genome → summary TSVs" (mobile),
  `antismash` "Bakta GBFF → antiSMASH JSON" (bgc).
- §7.1.7 uses: prophage "geNOMAD proviruses … geNOMAD markers"; plasmid "geNOMAD plasmid score"; bgc
  "antiSMASH regions | core biosynthetic genes flagged".
- R5: antiSMASH AGPL, geNOMAD non-commercial (user-guide note); never imported.

## Public interface

`src/giap/annotate/genomad.py`, `genomad_cmd.py`, `antismash.py`, `antismash_cmd.py`.
```python
def genomad_command(genome: Path, outdir: Path, db: Path, threads: int) -> list[str]
@dataclass(frozen=True, slots=True)
class Provirus: contig: str; interval: Interval | None; virus_score: float; fdr: float | None
                n_hallmarks: int; taxonomy: str          # interval None = whole contig classified viral
@dataclass(frozen=True, slots=True)
class PlasmidCall: contig: str; plasmid_score: float; fdr: float | None; n_hallmarks: int
                   conjugation_genes: tuple[str, ...]; amr_genes: tuple[str, ...]
@dataclass(frozen=True, slots=True)
class VirusGene: contig: str; interval: Interval; marker: str | None; description: str | None
                 virus_hallmark: bool
@dataclass(frozen=True, slots=True)
class GenomadResult: proviruses: tuple[Provirus, ...]; plasmids: tuple[PlasmidCall, ...]; virus_genes: tuple[VirusGene, ...]
def parse_genomad(outdir: Path, prefix: str, replicons: Mapping[str, Replicon]) -> GenomadResult
    # reads <prefix>_summary/<prefix>_{virus_summary,plasmid_summary,virus_genes}.tsv

def antismash_command(gbff: Path, outdir: Path, threads: int) -> list[str]
@dataclass(frozen=True, slots=True)
class BgcRegion: contig: str; interval: Interval; products: tuple[str, ...]; core_genes: tuple[str, ...]
def parse_antismash_json(path: Path, replicons: Mapping[str, Replicon]) -> list[BgcRegion]
```
geNOMAD rules: `seq_name` `<contig>|provirus_<a>_<b>` → `Interval(contig, a−1, b, ".")` (1-based
closed); `seq_name` without `|provirus_` → whole-contig virus (`interval=None`); virus gene names
`<seq_name>_<k>` → contig from the part before `|`, coordinates from `start`, `end`, `strand`
(1/−1) columns, offset by `a−1` for provirus genes (geNOMAD gene coordinates are relative to the
provirus sequence; Q2); `NA` → `None`. Contig IDs are safe IDs (M1-02), so `|` cannot occur in them.
antiSMASH rules: `records[].id` = contig; `records[].areas[]` `start`, `end` (0-based, end
exclusive) and `products`; `start > end` on a circular record → origin-crossing (`end + length`);
core genes = CDS features in the area whose qualifier `gene_kind` contains `biosynthetic` → their
`locus_tag` qualifier (Bakta locus tags, since the input is Bakta's GBFF).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `genomad_command(Path("genome.fa"), Path("out"), Path("/db/genomad_db"), 16)` | `["genomad","end-to-end","--cleanup","--splits","8","--threads","16","genome.fa","out","/db/genomad_db"]` |
| 2 | virus_summary row `chr1|provirus_12001_43000`, virus_score 0.97, fdr NA, n_hallmarks 9, taxonomy `Viruses;Duplodnaviria;…` | `Provirus("chr1", Interval("chr1",12000,43000,"."), 0.97, None, 9, "Viruses;Duplodnaviria;…")` |
| 3 | plasmid_summary row `p1`, plasmid_score 0.93, conjugation_genes `MOBP;T4SS_virB4` | `PlasmidCall("p1", 0.93, None, …, ("MOBP","T4SS_virB4"), ())` |
| 4 | virus_genes row `chr1|provirus_12001_43000_5`, start 201, end 1100, strand −1, marker `VV00123`, description `phage integrase` | `VirusGene("chr1", Interval("chr1",12200,13100,"-"), "VV00123", "phage integrase", …)` |
| 5 | headers only in all three TSVs (sim100k) | `GenomadResult((), (), ())` |
| 6 | `antismash_command(Path("bakta.gbff"), Path("out"), 8)` | `["antismash","--genefinding-tool","none","--minimal","--cpus","8","--output-dir","out","bakta.gbff"]` |
| 7 | area `{start: 1000, end: 25000, products: ["NRPS"]}` on `chr1`, CDS `ESC_1234` with `gene_kind: ["biosynthetic"]` inside | `BgcRegion("chr1", Interval("chr1",1000,25000,"."), ("NRPS",), ("ESC_1234",))` |
| 8 | `{"records": []}` | `[]` |

## Edge cases

- Provirus coordinates beyond the contig length → `ValueError` (corrupt input).
- A record id in antiSMASH not in the genome → `ValueError` naming it.
- Multiple products are kept in antiSMASH order; duplicate products collapsed preserving order.
- geNOMAD `--cleanup` removes intermediates; only `*_summary/` is kept in `annotation/tools/genomad/`.

## Acceptance criteria

- [ ] `tests/unit/test_genomad.py` (examples 1–5) and `tests/unit/test_antismash.py` (6–8) on
  hand-written files.
- [ ] `tests/contract/test_genomad_contract.py` and `test_antismash_contract.py` on the M2-10
  recordings (sim100k, k12, nissle); Nissle yields ≥ 1 provirus and ≥ 1 BGC region, K-12 ≥ 1 provirus;
  `pending_sample` until recorded.
- [ ] `test_nf_cmd_sync.py` covers both tools.

## Questions for the lead

1. Whole-contig virus calls (no `|provirus_`). **Proposed default:** kept with `interval=None`; M2-19
   turns them into a prophage region spanning the contig with note `whole_contig`.
2. Are geNOMAD gene coordinates relative to the provirus? **Proposed default:** treat them as
   relative for provirus genes (offset by `a−1`) and absolute otherwise; the K-12 contract test
   checks that every virus gene of a provirus lies inside it, which decides the rule.

## Lead decisions (after the test-author's pass)
1. The test-author's resolutions are accepted (geNOMAD gene coordinates absolute — no a−1 offset; core
   genes = exactly `biosynthetic`; CDS counted only when wholly inside the area; fuzzy and `join{}`
   locations parse; `start > end` on a linear contig and unknown contigs raise; `;`-split lists, `NA` →
   `()`; parsers read `genomad/genome_summary/` and `antismash/<name>.json`).
2. The ANTISMASH process and `antismash_cmd` pass `--databases <db_dir>/antismash` (the databases
   `giap setup` installed and verified; the DB directory is bound read-only since #78), as the recorded
   real runs did. The cmd-sync pair covers it.
3. No py step in M2-18; M2-19's `region_layer` registers the parsing step and must match annotate.nf's
   `stepIn('region_layer', …)` keys.

## Lead decisions after code review (PR #84, changes requested — tests/API)
4. `antismash_command(..., databases)` makes `databases` REQUIRED (real runs always need it); spec
   example 6 is amended accordingly; the cmd-sync pair passes it and the antiSMASH sync checks include
   `--databases` and `/antismash`.
5. Tests kill the surviving mutants: a biosynthetic CDS crossing the area's right edge is not a core gene
   ("wholly inside"); a `join{}` CDS with one part outside the area is excluded (`all`, not `any`); an
   area with `end > length` on a linear record raises ValueError; `provirus_0_…` naming at position 0.
6. tests/data/recorded/antismash/PROVENANCE.md is updated (main.nf now passes `--databases`).
