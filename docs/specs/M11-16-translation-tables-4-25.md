# M11-16: genetic codes 4 and 25 alongside 11 (`annotation.translation_table`)

Milestone M11 · pass 2 · depends on: #413 (preflight refuses `annotation.translation_table != 11`),
M11-14w (`genome.genetic_code` warning, `gippa/preflight/genetic_code.py`), M3-01 (six-frame
index), M7-16 (consumable stop cells, index format 2, ADR-0058), M2-03 (Bakta argv), M3-12/M3-13
(miniprot, tblastn), M2-17 (codon usage), M7 novel ORFs and variants, M3-23 (`makedb`) ·
estimate ≈ 1,190 code lines + ≈ 1,750 test lines in 7 PRs (table at the end) · **`stats`**: the
reverse direction of the `genome.genetic_code` rule (item 2.3) and the new calibration arms
(item 9) need the statistician's review · roles: implementer (`src/gippa/model/`,
`src/gippa/preflight/`, `src/gippa/mapping/`, `src/gippa/annotate/`, `src/gippa/genetable/`,
`src/gippa/novel/`, `src/gippa/variants/`, `src/gippa/makedb/`, `src/gippa/steps/`,
`src/gippa/manifest.py`, `src/gippa/report/`, `workflow/`, and `tests/sim/gippasim/`, which is
implementer tooling), test-author (`tests/unit/test_m11_16*_*.py`,
`tests/integration/test_m11_16*_*.py`, `workflow/tests/*.nf.test`, calibration arms),
validation (M. hyorhinis re-run, Slurm), docs · decision: **Khoa, 2026-10-28: support genetic
codes 4 and 25 alongside 11. This supersedes Q-K30** (2026-10-27: code 11 only for 1.0, code 4
after 1.0).

## Spec sections implemented

- GIPPA_SPEC.md §7.2.1 says "translate all six frames of every contig with table 11". §4 sets the
  config default `translation_table: 11`. §12.1 lists "table-11 translation" as a unit test. This
  issue lets the table be 4 or 25. The spec is read-only, so **ADR-0067 "Genetic codes 4, 11 and
  25"** (written in PR 16a) records the change and overrides those sentences. Table 11 stays the
  default.
- §7.0 preflight ("every evaluated rule yields one check"; `warn` never stops a run; `--strict`
  turns warnings into exit 3), §7.1.1 (Bakta, miniprot, BLAST+ pins), §7.2.1–§7.2.3, §7.5 variants,
  §7.9.2 `makedb`, §9.2 manifest, §6.5 (0-based half-open), §6.6 (randomness), R1, R3, R6, R14.
- ADR-0058 / M7-16 (stop cells, `sec_site` / `pyl_site`, recoded sites, index format 2).
- The backlog entry "Genetic code 4 support (post-1.0)" (docs/backlog.md) lists the places that
  assume table 11. This issue replaces it (16g removes the entry and points to this file).

The three tables (NCBI; Biopython `CodonTable.unambiguous_dna_by_id`, checked 2026-10-28):

| Table | NCBI name | Stops | TGA | Bakta 1.12.2 `--translation-table` | codon-bias `genetic_code` |
|---|---|---|---|---|---|
| 11 | Bacterial, Archaeal and Plant Plastid | TAA, TAG, TGA | stop (Sec when recoded) | yes (default) | yes |
| 4 | Mold/Protozoan/Coelenterate Mitochondrial; Mycoplasma; Spiroplasma | TAA, TAG | **W** (Trp) | yes | yes (column `4` of its `genetic_code_ncbi.csv`) |
| 25 | Candidate Division SR1 and Gracilibacteria | TAA, TAG | **G** (Gly) | yes | yes (column `25`) |

No other codon differs between these tables. So the whole change is "TGA is a stop under 11, W
under 4, G under 25". The other 63 codons translate the same way in all three tables.

## What exists today (main 1f7203d3, 2026-10-28)

Every place that assumes table 11. Found with grep over `src/`, `workflow/` and `tests/sim/`
for `translate11`, `STOP`, `stop_codons`, `STOP_CODONS`, `db_gencode`, `genetic_code`,
`translation_table`, `TGA`, `sec_site`, `Selenocysteine`, `-T 11`.

| Place | Today | Change (item) |
|---|---|---|
| `config/model.py:162`, `defaults.yaml:29`, `config.schema.json:185` | `translation_table: int = 11`, no allowed set | allowed set {4, 11, 25} checked in preflight (1) |
| `preflight/run.py:219` `_translation_table_check` | error unless 11 ("code 4 support is planned") | error unless in {4, 11, 25} (2.1) |
| `preflight/genetic_code.py` (M11-14w) | `info` when the table is not 11; the warning says "GIPPA 1.0 supports genetic code 11 only" | new wording at 11; a reverse check at 4/25 (2.3) |
| `annotate/bakta_cmd.py` `bakta_command`, `workflow/modules/local/bakta/main.nf` | no `--translation-table` | pass it (3.1) |
| `mapping/miniprot.py:31`, `workflow/modules/local/miniprot/main.nf:22` | `-T 11` hard-coded at indexing | `-T <table>` (3.2) |
| `defaults.yaml:64`, `model.py:188`, `workflow/conf/base.config:44` | `tblastn_args: "-db_gencode 11 -seg no -evalue 1e-5"` passed verbatim | `-db_gencode` derived from the table (3.3) |
| `mapping/sixframe.py` `_TABLE11_TCAG`, `translate11`, `translate11_text`, `_INDEX_FORMAT = 2` | one table-11 LUT; the index does not record its table | per-table LUTs; format 3 records the table; a mismatch is refused (4) |
| `mapping/readthrough.py:38-40` `STOP_CODONS`, `_RECODED` | {TAA, TAG, TGA}; U → TGA, O → TAG | from the table; U only under 11 (5.1) |
| `mapping/loci.py:69-71` `_PLUS_STOPS`, `_MINUS_STOPS`, `_RECODED_FLAG` | table 11 | from the table; `sec_site` only under 11 (5.2) |
| `steps/mapping.py:246-280` (`sixframe_index`, `map_*`, loci) | no table passed | pass `cfg.annotation.translation_table` (5.3) |
| `mapping/decoys.py`, `chance.py`, `qgram.py`, `tolerant*.py` | work on the residue codes (`STOP` = 19) | no change: they only see `*` where the index put one (5.4) |
| `annotate/codon_usage.py:58` `_code`, `steps/codon.py:67` | raises for any table but 11 | accept 4 and 25 (6.1) |
| `genetable/covariates.py:128-160` | already table-aware (Biopython) | tests only (6.2) |
| `novel/candidates.py:70` `_stops(table)` | already table-aware | no change |
| `novel/outputs.py:125-126, 290` | `_stops(11)`, `translate11_text` | the run's table (6.3) |
| `novel/evidence.py:33` `_STOPS` | {TAA, TAG, TGA} | the run's table (6.3) |
| `variants/flags.py:187`, `vcf_records.py:64-81, 135`, `align.py:97`, `pileup.py:69, 248` | `translate11_text`; U → TGA zero-change ALT | the run's table; U only under 11 (6.4) |
| `makedb/build.py:102-363` | `translate11_text`, table-11 six-frame ORFs | the configured table; recorded in `makedb_manifest.json` (6.5) |
| `annotate/bakta.py` `parse_bakta_json` | ignores the JSON `genome.translation_table` | checks it against the config (3.4) |
| `manifest.py:379` `_annotation`, `run_manifest.schema.json:114` | no table | `annotation.translation_table` (7) |
| `report/content.py:102` | shows the M11-14w warning | also shows the table (7) |
| `tests/sim/gippasim/codons.py` `CODE_11`, `STOPS`, `SENSE_USAGE`; `simulate.py:88` `stop_codons` (TGA 0.30); `layout.py:62` `_STOPS`; `readthrough.py`, `strain.py`, `proteomics.py`, `features.py`, `scenario_truth.py`, `real_genome.py`, `function_sim.py`, `selenoproteins.py` | table 11 everywhere | `SimConfig.genetic_code` (8) |
| `workflow/modules/local/genomad/main.nf`, antiSMASH, OSTIR, promoters, terminators, UniOP | geNOMAD picks its own code (prodigal-gv); antiSMASH reads Bakta's GBFF `/transl_table`; the rest are nucleotide-only | no change |

## The gap

M. hyorhinis SK76 (GCF_000313635.1, PR #411) uses code 4. It ran at table 11, and the code-11
output in `$GIPPA_ROOT/more_datasets/mhyo/out` shows what goes wrong. The spec-writer recomputed
these figures on 2026-10-28 from `genes.parquet`, `mapping/peptide_hits.tsv` and
`mapping/unmapped.tsv`:

- Bakta called 1,474 features: 1,200 CDSs and 274 pseudogenes. RefSeq has about 750.
- **801 of the 1,200 CDSs end in TGA** (TAA 340, TAG 59).
- The median protein is 90 aa over all features and 80 aa over CDSs.
- 390 pairs of adjacent CDSs sit on the same strand and in the same frame, ≤ 300 bp apart, where
  the upstream CDS ends in TGA. These are genes cut at a TGA codon.
- `peptide_hits.tsv` flags: 975 `nterm_extension_evidence`, 321 `strain_difference` and 316
  `readthrough`.
- 8,163 organism peptides have no hit (`reason = no_hit`, `label = unmapped`), and 1,524 of them
  contain W.
- The M11-14w ratio is 0.227.

The backlog entry says "all 628 unmapped organism peptides contain W" and "median 90 aa". The
median agrees with the all-feature figure. The 628 count comes from a narrower peptide set that
the backlog does not define. 16g re-derives it (see Q-L8).

Users can now set `annotation.translation_table: 4` or `25`. Every step must then translate,
search and flag with that table. A six-frame index or a precomputed annotation built for another
table must be refused, never silently reused.

## Definition

### 1. The table module (new, pure): `gippa/model/genetic_code.py`

- `SUPPORTED_TABLES: tuple[int, ...] = (4, 11, 25)`.
- `codon_table(table: int) -> Mapping[str, str]`: the 64 upper-case ACGT codons → one letter,
  with `*` for stops. It is built from the NCBI TCAG strings:
  - 11: `FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG`
  - 4: `FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG`
  - 25: `FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG`

  A test checks each one against Biopython. It is not imported from Biopython at run time, so
  numba and the LUT builders do not depend on it.
- `stop_codons(table) -> frozenset[str]`: `{TAA, TAG, TGA}` for 11 and `{TAA, TAG}` for 4 and 25.
- `minus_stop_codons(table) -> frozenset[str]`: the reverse complements, read on the plus string
  (`{TTA, CTA, TCA}` and `{TTA, CTA}`).
- `recoded_codons(table) -> Mapping[str, str]`: the Bakta letter → the stop codon it must sit over.
  `{"U": "TGA", "O": "TAG"}` for 11, `{"O": "TAG"}` for 4 and 25. TGA is a sense codon under 4 and
  25, so Sec cannot be recoded over it. Pyl over TAG stays possible.
- `table_name(table) -> str`: `"Bacterial, Archaeal and Plant Plastid"` (11), `"Mycoplasma and
  Spiroplasma"` (4), `"Candidate Division SR1 and Gracilibacteria"` (25).
- Every function raises `ValueError(f"unsupported translation table {table}; GIPPA supports 4, 11 and 25")`
  for any other value. Results are cached (`functools.cache`), and the module has no I/O.

Every other item takes the table from `cfg.annotation.translation_table` and calls this module.
No module keeps its own stop set: `novel.candidates._stops` becomes a call to `stop_codons`.

### 2. Config and preflight

**2.1 `config.translation_table`** (`_translation_table_check`, same rule id):

| Value | Level | Message |
|---|---|---|
| 4, 11 or 25 | `ok` | `""` |
| anything else | `error`, exit `INVALID_INPUT` (as #413) | `annotation.translation_table={t}: GIPPA supports genetic codes 4, 11 and 25 (the tables Bakta accepts)` |

The `rule` text becomes `annotation.translation_table in {4, 11, 25} → ok; else error`. The
pydantic model stays `int`. The schema gets `"enum": [4, 11, 25]` in its description only, not as
a constraint, so a bad value still reaches the preflight table and is not a pydantic traceback
(decision D7). `value` is the table.

**2.2 `config.tool_genetic_code`** (new, one check, `input = "inputs"`). GIPPA sets the code of
every tool from the one key, so a tool option that sets a code itself is a conflict:

- `annotation.bakta.extra_args` contains `--translation-table` (as its own token or as
  `--translation-table=N`), whatever its value;
- `mapping.proteins.miniprot_args` contains `-T` (as its own token, or `-T<N>`). `-T` is an
  indexing option that GIPPA sets (§7.1.1 gotcha);
- `mapping.proteins.tblastn_args` contains `-db_gencode N` (or `-db_gencode=N`) with `N != table`,
  or contains `-db_gencode` more than once.

Any of these gives `error` (`INVALID_INPUT`). The message names each offending key and tells the
user to remove the option and set `annotation.translation_table`. `value` = the list of offending
keys, or `[]`. If there is no conflict, the check is `ok`. `-db_gencode N` with `N == table` is
allowed, so configs that today spell out `-db_gencode 11` keep working at 11.

**2.3 `genome.genetic_code` (M11-14w), now evaluated at every supported table.** The scan,
`bp11`, `bp4`, `n4`, `MIN_CODONS = 300`, `MIN_LONG_ORFS = 20`, the threshold 0.55 and the exact
integer comparison stay as they are. The decision now depends on the table:

| Condition | Level | Message |
|---|---|---|
| table not in {4, 11, 25} | `info` | `not evaluated: unsupported annotation.translation_table` (`config.translation_table` already errors) |
| no records, or `n4 < 20` | `info` | unchanged (`not evaluated: fewer than 20 open reading frames of at least 300 codons`) |
| table 11 and `20 * bp11 < 11 * bp4` (ratio < 0.55) | `warn` | W11 below |
| table 4 or 25 and `20 * bp11 >= 11 * bp4` (ratio ≥ 0.55) | `warn` | W4 below |
| otherwise | `ok` | `""` |

The genome alone cannot tell code 4 from code 25: both read TGA as a sense codon. So W11 names
both codes, and the rule never warns "4 instead of 25" (edge cases).

`value` gains `"translation_table": table`. The other keys are unchanged. The `rule` text becomes:
`long ORFs (>= 300 codons, six frames): bp with TGA as stop / bp with TGA as sense; table 11:
< 0.55 -> warn; tables 4, 25: >= 0.55 -> warn (needs >= 20 such ORFs)`.

The ratio is printed truncated as `{r // 100}.{r % 100:02d}` with `r = (100 * bp11) // bp4`.
Today's `0.{r:02d}` prints `0.100` when `bp11 == bp4`, and that can now happen in W4.
Truncation keeps W11 from printing "0.55 < 0.55" and keeps W4 from printing a value below 0.55.

**W11** (table 11; `{q4}`, `{q11}` whole percentages as today):

> This genome probably uses a genetic code in which TGA is not a stop: code 4 (TGA = Trp;
> Mycoplasma, Spiroplasma, Mesoplasma, Ureaplasma and other Mollicutes) or code 25 (TGA = Gly;
> candidate division SR1, Gracilibacteria). Long open reading frames cover {q4} % of the genome
> when TGA is read as a codon but only {q11} % when it is a stop (ratio {ratio} < 0.55). This run
> uses annotation.translation_table: 11, so genes will be cut at TGA, peptides with W or G at TGA
> codons will not map, and strain_difference, readthrough and nterm_extension_evidence flags will
> be artefacts. Set annotation.translation_table to 4 or 25 from the organism's taxonomy; the
> genome sequence alone cannot tell the two apart. This warning does not change any result.

**W4** (table 4 or 25; `{table}`, `{name}` = `table_name(table)`):

> This genome probably uses genetic code 11 (TGA = stop). Long open reading frames cover {q11} %
> of the genome when TGA is a stop and {q4} % when it is read as a codon (ratio {ratio} >= 0.55;
> genomes that use code 4 or 25 are below 0.55). This run uses annotation.translation_table:
> {table} ({name}), so genes that end in TGA will be read through into the next gene and
> peptides may be placed in the wrong frame. Most bacteria and archaea use code 11; check the
> organism's genetic code. This warning does not change any result.

The level stays `warn` (decision D5) and the warning never changes results. `--strict` turns it
into exit 3, as M11-14w Q2 decided. The report's `SUMMARY_WARNINGS` alert shows W11 and W4 the
same way, with no change to `report/content.py` beyond item 7.

**Why the same 0.55 cut works for W4** (stats review, decision D4). The M11-14w evidence
(61 official NCBI assemblies plus draft stress tests) has every code-4/25 genome at ≤ 0.470
(≤ 0.477 as drafts) and every code-11 genome at ≥ 0.679 (≥ 0.641 as drafts). Read the other way,
a code-4/25 run whose genome scores ≥ 0.55 is a code-11 genome. On that panel this gives 0
false W4 on the 22 code-4/25 assemblies and 0 missed code-11 genomes. No new constant is needed.

### 3. Tool plumbing

**3.1 Bakta.** `bakta_command(..., translation_table: int)` adds `--translation-table
<t>` only when the table is not 11, placed before `--replicons` (decision D2, Khoa 2026-10-28: no
flag at 11, Bakta's own default, so code-11 caches stay valid). The Nextflow `BAKTA` module adds
the same flag, under the same condition, at the same place. The
existing argv parity test (Python argv vs the module script) covers both. The `--regions` route
(user GFF3/GBFF annotation) is unchanged: Bakta re-calls the CDSs with the table.

**3.2 miniprot.** `miniprot_index_cmd(genome, mpi, threads, table)` gives `-T <table>`, and the
module gives `-T ${params.annotation?.translation_table ?: 11}`. miniprot 0.18 accepts NCBI tables
1–33 except 7–8 and 17–20. The 16d contract test runs `miniprot -T 4` and `-T 25` on a 2-gene
fixture in the pinned env, on Slurm.

**3.3 tblastn `-db_gencode`.** The effective arguments are
`effective_tblastn_args(args: str, table: int) -> str` (pure, in `mapping/tblastn.py`):
- if `args` has no `-db_gencode`, append ` -db_gencode <table>` (after the user's tokens; the
  result is stripped of surrounding spaces);
- if it has `-db_gencode <table>`, return `args` unchanged (no duplicate; BLAST+ rejects a repeated
  argument);
- if it has a different value, raise `ValueError`. Preflight 2.2 has already refused such a run;
  the raise is defence in depth.

`tblastn_cmd` uses it. The Nextflow side mirrors it in `workflow/conf/base.config`
(`withName: 'TBLASTN'`, a closure over `params.mapping.proteins.tblastn_args` and
`params.annotation.translation_table`). A parity test runs both on the same case table (examples
7a–7c). The default `mapping.proteins.tblastn_args` becomes `"-seg no -evalue 1e-5"`. This changes
a value in the protected `src/gippa/config/defaults.yaml`, so it needs Khoa's approval, ADR-0067
accepted-by-khoa and the `Gate-Change-Approved: ADR-0067` trailer (decision D1, with a fallback).
The stale comment `# -T 11 is applied at indexing time` (defaults.yaml:63) becomes `# -T
<annotation.translation_table> is applied at indexing time` in the same change.

**3.4 Precomputed Bakta JSON (`--annotation x.bakta.json`, `gippa makedb --annotation`).**
`parse_bakta_json` reads `genome.translation_table` when the JSON has it. If the value differs from
the config, it raises `AnnotationTableMismatch(ValueError)`:
`annotation JSON was made with translation table {j}; this run uses {t}: re-run Bakta with
--translation-table {t} or set annotation.translation_table: {j}`. The step then exits
`INVALID_INPUT`. If the key is absent (older or hand-made JSON), the config table is assumed, and
the parser logs `annotation.table_unrecorded` at `info`. `annotation_meta.json` gains
`translation_table` (the config value, which equals the JSON value after the check). Preflight does
not parse the JSON (decision D9).

### 4. Six-frame index

- `sixframe.py` builds `_CODON_LUT` / `_CODON_ASCII` per table from `codon_table(table)`, cached.
  The numba kernel is unchanged: it already takes the LUT as an argument.
- `translate(nt: bytes, table: int = 11)` and `translate_text(nt: str, table: int = 11, offset: int = 0)`
  are the new names. `translate11` and `translate11_text` stay as aliases with `table=11`
  (decision D13), so call sites that are not table-specific do not change. Every call site in
  items 6.3–6.5 passes the table explicitly.
- `build_sixframe(genome, replicons, wrap_nt=WRAP_NT, table=11)`. `SixFrameIndex` gains
  `translation_table: int = 11`. Under 4, TGA becomes `W` (code 17); under 25 it becomes `G`
  (code 5). The `STOP` code (19) appears only at TAA/TAG under 4 and 25.
- **Index format 3.** `_INDEX_FORMAT = 3`, and `meta` gains `"translation_table": t`.
  `save_index` writes `index.translation_table`. `load_index(path, *, table: int)` and
  `load_stop_cells(path, *, table: int)` take the run's table as a required keyword, so no caller
  can forget it:
  - format 3 with `meta.translation_table == table`: load;
  - format 3 with a different table: `ValueError(f"{path}: six-frame index built for translation
    table {m}; this run uses {table}; rebuild the six-frame index")`;
  - format 2 (M7-16; it has no table key and was always built with table 11): load when `table == 11`.
    Otherwise raise the mismatch error with `m = 11` (decision D3);
  - format 1: unchanged ("has no stop cells; rebuild").
- Nextflow `-resume` already re-runs `sixframe_index` when the table changes, because
  `annotation` is in every step's config slice (`steps/registry.py:183`). The format check covers
  files passed by hand (`gippa step`, copied work dirs).
- **Consumable-stop cells.** `annotated_stops(gm, genome, recoded=False, *, table=11)` uses
  `stop_codons(table)` for terminal sites and `recoded_codons(table)` for recoded sites (5.1).
  `frame_cells` and `mirror_cells` are unchanged. Their `aa[p] == STOP` check now also catches a
  site that should not exist under the table.

### 5. Mapping

**5.1 `mapping/readthrough.py`.** `STOP_CODONS` becomes `stop_codons(table)`. `_RECODED` becomes
`recoded_codons(table)`. The `StopSite.codon` Literal stays `TAA | TAG | TGA` (TGA occurs only
under 11). Rule 21(e) under 4/25: a `U` in an annotated protein has no recoding codon in the
table. It gives no site and is counted under `recoded_skipped.codon_mismatch`, once per residue,
the same as a `U` over a non-TGA codon today (decision D6). An `O` over TAG is a `recoded` site in
all three tables. A CDS whose last codon is TGA under 4/25 is counted under `not_a_stop` (rule 1)
and gives no terminal site.

**5.2 `mapping/loci.py`.** `_PLUS_STOPS` and `_MINUS_STOPS` become `stop_codons(table)` and
`minus_stop_codons(table)`, reached through a `table` argument on the classifier (`_classify`,
`_nterm_plus`, `_nterm_minus`, `LocusClassifier` and `annotate_loci`). The N-terminal extension
upstream walk then stops only at the table's stops. `_RECODED_FLAG = {"TGA": "sec_site", "TAG":
"pyl_site"}` stays a codon → flag map. Under 4/25 no TGA site is `recoded` (5.1), so `sec_site` is
never written. `pyl_site` still is. No vocabulary or schema change. `context/summary.py`
`FLAG_GROUPS` keeps `sec_site`, whose count is 0 under 4/25.

**5.3 `steps/mapping.py`.** `sixframe_index`, `map_exact`, `map_tolerant`, the loci/placement step
and the protein routes pass `ctx.cfg.annotation.translation_table` to `build_sixframe`,
`annotated_stops`, `load_index`, `load_stop_cells` and the classifier. The `mapping.stop_sites`
event `data` gains `"translation_table": t`.

**5.4 No change** to `decoys.py` (reverse and shuffle mirror the target frames, and stops keep
their positions), `chance.py` (residue frequencies come from the index, so they follow the
table), `qgram.py` or `tolerant*.py`. Under 4/25 there are fewer `*` per frame, so stop-free
segments are longer and the shuffle windows (`shuffle_window_aa = 50`) are unchanged. Item 9
checks the decoy balance and chance-hit calibration under 4 and 25.

### 6. Other places

**6.1 Codon usage** (`annotate/codon_usage.py`, `steps/codon.py`). `_code(table)` accepts 4, 11
and 25 through `codon_table`. codon-bias is already called with `genetic_code=table`. Effects:
- under 4, Trp is a 2-codon family (TGA, TGG) and gets synonymous shares and a rare-codon status;
- under 25, Gly is a 5-codon family (GGN + TGA);
- the eligibility rule ("no internal stop in codons 2..n−1") uses the table's stops, so a code-4 CDS
  with an internal TGA is eligible.

tAI uses the anticodon counts from Bakta's tRNAs as today. The Mycoplasma tRNA-Trp(UCA) reads
TGA and TGG. Whether codon-bias's `TrnaAdaptationIndex(genetic_code=4)` gives TGA a weight > 0
from anticodon `TCA` is checked in 16e (AC 11). The module docstring drops "only translation
table 11 is accepted".

**6.2 Covariates** (`genetable/covariates.py`). Already table-aware: `protein_sequence` translates
`nt` with Biopython's table when Bakta gave no `aa`. 16e adds tests at 4 and 25 only.

**6.3 Novel ORFs.** In `novel/outputs.py`, `_walk_up(sv, t_end, _stops(11))` uses the run's stops,
and the two `translate11_text` calls (lines 126, 290) use `translate_text(..., table)`.
`novel/evidence.py` `_STOPS` (the internal-stop check) becomes `stop_codons(table)`, passed from
the step's config. `candidates.reconstruct` is unchanged: it already uses
`cfg.annotation.translation_table`. The novel FASTA (`novel_orfs.faa`) of a code-4 ORF has `W` at
TGA. `annotation.start_codons` stays a config key and is not derived from the table (edge cases).

**6.4 Variants.** `flags.py` (block mass equivalence), `vcf_records.py` (`AA_REF`, `_codons_of`
and `minimal_alt_codon`) and `align.py` (locus translation) use `translate_text(..., table)`.
`_codons_of(aa, table)` then gives `W → (TGA, TGG)` under 4 and `G → (GGA, GGC, GGG, GGT, TGA)`
under 25. `pileup.RECODED_STOP` and `vcf_records._RECODED` become `recoded_codons(table)`. Under
4/25 a `U` is never a zero-change ALT, so it is `rare`-skipped like every other U (M7-16 rule 11).
`O` over TAG is unchanged. A strain SNV TGA → TGG is a synonymous W → W under 4, not a stop loss.

**6.5 `gippa makedb`.** `MakedbOptions.translation_table: int = 11` is filled from the resolved
config (`--config` / `--set annotation.translation_table=4`; no new CLI flag). Six-frame ORFs,
`nterm_extension` entries and `_protein` use it. The annotation JSON table is checked (3.4), and
`makedb_manifest.json` gains `"translation_table"`. A `--patch` CDS with an internal `*` under the
table is rejected as today. Under 4/25, TGA is no longer `*`, so it is not rejected for a TGA.

### 7. Manifest and report

- `run_manifest.json` `annotation` gains `translation_table` (integer). It is read from
  `annotation_meta.json`, or from the resolved config when the annotation block is null. The
  schema gets an explicit `translation_table: {"type": "integer", "enum": [4, 11, 25]}` property
  next to the string-valued `additionalProperties`. This bumps `OUTPUT_SCHEMA_VERSION` by one
  minor version (additive, M11-14h policy). The config snapshot in the manifest already holds the
  table and is unchanged.
- **Report.** The Methods section gains one sentence for every run: `Genes were called and all
  six frames translated with NCBI genetic code {t} ({name}).` When `t != 11`, the `summary`
  section also opens with `<p class="note" data-genetic-code="{t}">Genetic code {t} ({name}):
  TGA is read as {W|G}.</p>`, after any `SUMMARY_WARNINGS` alert. The Methods sentence changes the
  code-11 snapshots once; the note does not appear at 11 (decision D8).

### 8. Simulator (gippasim, implementer tooling; tests by the test-author)

- `SimConfig.genetic_code: int = 11`, in {4, 11, 25}, else `ValueError`. **With 11, every output
  is byte-identical to today** for every existing scenario and seed. This is the first test of PR
  16c, and the 16c fixtures prove it.
- `gippasim.codons`: `code(table) -> dict[str, str]` and `stops(table) -> frozenset[str]`.
  `CODE_11` and `STOPS` stay as the table-11 values. Every module that reads them (`features`,
  `layout`, `readthrough`, `strain`, `proteomics`, `scenario_truth`, `real_genome`,
  `function_sim`, `signals`, `simulate`, `blocks`) takes the table from the config.
- Sense-codon usage: `sense_usage(table)`. For 11 it is today's `SENSE_USAGE`. For 4, TGA joins Trp
  and takes `SimConfig.tga_sense_share` of Trp's total (default 0.6), with TGG keeping the rest.
  For 25, TGA joins Gly with `tga_sense_share` of Gly's total (default 0.05). The total per-mille
  weight of each amino acid is unchanged, so protein composition is the same in all three tables.
  `tga_sense_share` must be in (0, 1). It is ignored at 11.
- Stop weights: under 4/25 the default `stop_codons` drops TGA and renormalises TAA/TAG
  (0.61/0.09 → 0.871/0.129). A user `stop_codons` that contains TGA under 4/25 raises.
- `n_selenoproteins > 0` with `genetic_code != 11` raises (Sec needs TGA as a stop). `n_readthrough`
  works: terminal stops are TAA/TAG and the readthrough residue map (`TAA`/`TAG` → Q) is unchanged.
- The simulated Bakta JSON writes `genome.translation_table` only when `genetic_code != 11`, so the
  code-11 fixtures stay byte-identical. The truth tables are unchanged in format.
- `real_genome.py` (peptides from a real genome) gains the same `genetic_code` argument.
- The M11-14w fixture generator (`n_genes = 300`, 400 codons, …) gains code 25 (TGA → Gly, rate
  `r`). The code-4 path keeps its seeds and values.

### 9. Calibration (R3: no tolerance or seed changes; `tolerances.yaml` is not touched)

There are no new suites. The existing suites whose statistic depends on where stops fall gain a
`genetic_code` parameter with arms 11 (today, unchanged and still first), 4 and 25. The arms
reuse the suite's existing tolerance keys, replicate counts and seeds:

| Suite | Why the table matters | Arms in the fast tier | Arms in `calibration_full` |
|---|---|---|---|
| `test_exact_recall.py` (≥ `mapping_recall_min` 0.999) | W/G at TGA must map | 4 | 4, 25 |
| `test_tolerant_recall.py` | fewer stops, longer frames | 4 | 4, 25 |
| `test_chance_hits.py` | longer stop-free frames → more chance hits; the chance model must stay calibrated | 4 | 4, 25 |
| `test_placement_fdr_sim.py` | target-decoy balance with fewer stops | none | 4, 25 |
| `test_readthrough_decoy_balance.py` | cells only at TAA/TAG | none | 4, 25 |
| `test_novel_class_fdr_calibration.py` | ORF lengths and the upstream walk | none | 4, 25 |
| `test_sim_layout_calibration.py` (the simulator's own checks) | sense usage and stop weights at 4/25 | 4, 25 | 4, 25 |

An arm that fails is a finding and is reported, never loosened. Arms 4 and 25 reuse the seeds of
arm 11, which keeps the protected seed list unchanged. The statistician reviews the
chance-hit and placement arms before 16f merges (`stats`). The fast-tier cost is 3 suites × 1 arm
(about +6 min on the 8-CPU allocation, by the suites' current durations). It runs on Slurm via
`test_branch.sbatch`, never in the session.

### 10. Real-data acceptance (validation agent; Slurm, R14)

Re-run M. hyorhinis SK76 at code 4 with the #411 settings:

```
sbatch -J gippa-run-mhyo-code4 scripts/slurm/gippa_run.sbatch --genome $GIPPA_ROOT/more_datasets/mhyo/raw/genome/GCF_000313635.1_ASM31363v1_genomic.fna.gz \
  --peptides $GIPPA_ROOT/more_datasets/mhyo/inputs/peptides.tsv.gz --db-dir $GIPPA_DB_DIR --threads 6 --memory 24.GB \
  --set annotation.bakta.complete=true --set annotation.translation_table=4 --outdir $GIPPA_ROOT/more_datasets/mhyo/out_code4
```

Use `executor { queueSize = 4 }` in `OUT/nextflow.config`, as #411's later runs did. Record the job
ID in PROGRESS.md. Fetch the RefSeq annotation of GCF_000313635.1 (`_genomic.gff.gz`) from the NCBI
FTP and check it against `md5checksums.txt` (R6). Write `validation/datasets/mhyo/CODE4.md` with
these comparisons of `out` (code 11) and `out_code4`:

| Measure | Code 11 (recomputed above) | Expected at code 4 |
|---|---|---|
| preflight `genome.genetic_code` | `warn` (W11), ratio 0.227 | `ok`, ratio 0.227, `translation_table: 4` |
| Bakta CDS / pseudogene | 1,200 / 274 | CDS within ±10 % of RefSeq's protein-coding genes; pseudogenes ≤ 25 % of 274 |
| CDSs ending in TGA | 801 | 0 |
| median CDS protein | 80 aa | ≥ 250 aa |
| RefSeq CDSs with a Bakta CDS on the same strand and stop codon | (report) | ≥ 90 % |
| the 390 TGA-split adjacent pairs | 390 | ≥ 90 % lie inside one code-4 CDS (same strand, frame) |
| organism peptides `no_hit` / those with W | 8,163 / 1,524 | W-containing `no_hit` ≤ 20 % of 1,524 |
| distinct mapped peptides (target rows) | 10,214 | ≥ 10,214 (report the gain) |
| flags `readthrough` / `strain_difference` / `nterm_extension_evidence` | 316 / 321 / 975 | each ≤ 25 % of the code-11 count |
| `sec_site` | 0 | 0 |
| `trna_counts.tsv` anticodon `TCA` (tRNA-Trp UCA) | absent | report whether Bakta calls it |
| codon usage: Trp family in `codon.parquet` synonymous shares | n/a | TGA and TGG both present |

The "expected" column is the validation expectation, not a CI gate (decision D12). A miss is
reported as a finding with its cause. In the same session, a code-11 control re-runs the
GIPPA-relevant steps of E. coli K-12 (`$GIPPA_ROOT/data/genomes/k12_mg1655/`). Its tables must be
byte-identical to the last main run, apart from the documented changes (AC 1).

## Examples (input → expected output)

Coordinates are 0-based half-open. "+0" is plus strand, frame offset 0.

1. **Translation.** `translate_text("ATGTGGTGATAA", table)`: 11 → `"MW**"`, 4 → `"MWW*"`,
   25 → `"MWG*"`. `translate_text("ATGTGN", 4)` → `"MX"` (ambiguous codon). The trailing partial
   codon is dropped in every table.
2. **Six frames, linear, both strands.** Contig `c` (linear, 15 nt) = `ATGGCCTGAGCCTAA`. Frame +0:
   11 → `MA*A*`, 4 → `MAWA*`, 25 → `MAGA*`. The same contig given as its reverse complement
   (`TTAGGCTCAGGCCAT`) gives the same strings in frame −0, and `to_interval(-0 frame, 2, 3)` is
   `[6, 9) −` (the TCA on the plus string). In table 11 that cell is `*`; in table 4 it is `W`.
3. **Exact mapping of a W peptide.** Contig `g` (linear, 30 nt) =
   `ATG AAA GCC TGA GCC GAA GCC TGG AAA TAA` (no spaces). The peptide `AWAEAWK` (7 aa):
   - table 4: one target row `[6, 27) +`, `cds_inframe` when a CDS `[0, 30) +` is annotated;
   - table 11: no exact row, because frame +0 reads `MKA*AEAWK*`;
   - table 25: no exact row for `AWAEAWK`; `AGAEAWK` maps at `[6, 27) +` instead.

   The reverse complement of `g` as the contig gives `[3, 24) −` for the same peptide under 4.
   On a circular contig with the gene across the origin (start 20), the hit is `[26, 47) +`,
   unwrapped, `origin_spanning` (§6.5).
4. **Index refused across tables.** `save_index(p, build_sixframe(G, R, table=4), …)` then
   `load_index(p, table=11)` raises `ValueError` "built for translation table 4; this run uses 11;
   rebuild the six-frame index". `load_index(p, table=4)` loads. A format-2 file (M7-16) loads
   with `table=11` and raises with `table=4` ("built for translation table 11").
5. **Stop sites (M7-16 example 21 under each table).** Contig `c` linear, 60 nt. G1 `cds +`
   `[3, 33)` = `ATG GCT GCT GCT TGA GCT GCT GCT GCT TAA`.
   - `aa = "MAAAUAAAA"`, table 11, `recoded=True`: sites `(c,+,15,(G1,),TGA,recoded)` and
     `(c,+,30,(G1,),TAA,terminal)`; skips `{length_mismatch 0, codon_mismatch 0}` (as today).
   - Same `aa`, table 4: one site `(c,+,30,(G1,),TAA,terminal)`; skips `{length_mismatch 0,
     codon_mismatch 1}`. A peptide over codon 15 reads `W` in the index; it is never `sec_site`.
   - `aa = "MAAAWAAAA"` (what Bakta writes at code 4), table 4: one terminal site; skips 0. Table
     25 with `aa = "MAAAGAAAA"`: the same.
   - G2 `cds +` `[33, 45)` = `ATG GCT GCT TGA`, table 4: no site, and `not_a_stop` = 1.
6. **Pyl survives.** Under 4 and 25, a CDS with `O` over TAG inside gives a `recoded` TAG site, and
   a peptide consuming it gets `pyl_site`, as under 11.
7. **tblastn arguments** (`effective_tblastn_args`; the Nextflow closure gives the same string):
   - a. `("-seg no -evalue 1e-5", 4)` → `"-seg no -evalue 1e-5 -db_gencode 4"`;
   - b. `("-db_gencode 25 -seg yes -evalue 1e-7", 25)` → unchanged;
   - c. `("-db_gencode 11 -seg no", 4)` → `ValueError`. Preflight gives `config.tool_genetic_code`
     `error`, `value = ["mapping.proteins.tblastn_args"]`;
   - d. `("", 11)` → `"-db_gencode 11"`.
8. **Bakta and miniprot argv.** At table 25, `bakta_command(...)` contains `"--translation-table",
   "25"` right before `"--replicons"`, and `miniprot_index_cmd(...)` contains `"-T", "25"`. With
   `annotation.bakta.extra_args: ["--translation-table", "4"]`, preflight gives
   `config.tool_genetic_code` `error`, `value = ["annotation.bakta.extra_args"]`, exit
   `INVALID_INPUT`.
9. **Preflight `genome.genetic_code` at each table** (genomes and ratios from the M11-14w table):
   - M. hyorhinis at 11 → `warn` W11, `ratio 0.22`;
   - M. hyorhinis at 4 → `ok`;
   - E. coli K-12 at 4 → `warn` W4, `ratio 0.88`, `(Mycoplasma and Spiroplasma)`;
   - E. coli K-12 at 11 → `ok`;
   - the synthetic code-25 fixture at 11 → `warn`, and at 25 → `ok`;
   - `_decide(550, 1000, 20, table=4)` → `warn` and `_decide(549, 1000, 20, table=4)` → `ok`, the
     complement of 11 (`ok` and `warn`);
   - `_decide(1000, 1000, 20, table=25)` → `warn`, with `ratio 1.00` printed (not `0.100`).
10. **Codon usage.** A gene `ATG TGG TGA TGA AAA TAA`; codons 2..n−1 are `TGG TGA TGA AAA`.
    - table 4: eligible; `synonymous_shares` gives Trp `{TGA: 2/3, TGG: 1/3}`;
    - table 25: eligible; TGA counts in the Gly family, and TGG is a 1-codon Trp family with no share;
    - table 11: not eligible (internal stop), as today.
11. **Novel ORF protein.** ORF nt `ATG AAA TGA GCC TAA`: under 4 the FASTA protein is `MKWA`, under
    25 `MKGA`, and under 11 the ORF ends at the TGA (`MK`).
12. **Variants codon sets.** `_codons_of("W", 4)` = `("TGA", "TGG")`; `_codons_of("W", 11)` =
    `("TGG",)`; `_codons_of("G", 25)` = `("GGA", "GGC", "GGG", "GGT", "TGA")`. Under 4, a peptide `U`
    over a TGA codon is `rare`-skipped, never a zero-change ALT.
13. **Precomputed JSON mismatch.** A Bakta JSON with `genome.translation_table: 11` and a run at 4:
    `parse_bakta` exits `INVALID_INPUT` with the 3.4 message. The same JSON at 11 parses, and
    `annotation_meta.json` has `translation_table: 11`. A JSON without the key parses at either
    table and logs `annotation.table_unrecorded`.

## Edge cases

| Case | Expected |
|---|---|
| table 4/25, **circular** replicon, TGA codon across the origin (bases n−2, n−1, 0) | read as W/G in the wrapped frame. An origin-spanning peptide over it maps (unwrapped end), and the `origin_spanning` flag is as today |
| table 4/25, **linear** contig ending in `…TG` + `A` missing | partial codon dropped; no wrap; open frame as today |
| **minus strand** | stops on the plus string are `minus_stop_codons(table)`. `TCA` on the plus string is a stop only under 11, so the `loci` minus-strand N-terminal walk and `novel.evidence` both use it. Example 2 covers it |
| **empty input**: genome with no records / all N | preflight `genome.parse` errors as today; `genome.genetic_code` `info` with `ratio: null` at every table; `build_sixframe({}, {}, table=4)` gives an index with no frames, saved and loaded with its table |
| no peptides, no proteins | the table still reaches Bakta, codon usage and the manifest; mapping steps are skipped as today |
| **ties**: ratio exactly 0.55 | table 11 → `ok`, tables 4/25 → `warn`. Each genome is "code-11-like" iff ratio ≥ 0.55, so the two decisions never disagree on one genome |
| ratio 1.00 (`bp11 == bp4`) at table 4/25 | `warn`, prints `1.00` |
| `-db_gencode` twice in `tblastn_args` | `config.tool_genetic_code` error, even if both values equal the table |
| `annotation.translation_table` given as `"4"` through `--set` | pydantic coerces it to 4, as for any int key; `4.5` or `"four"` is a config error as today |
| **selenoproteins under 11** | unchanged: `U` over TGA is a recoded site, the `sec_site` flag, and the variants zero-change ALT (M7-16) |
| **selenoproteins under 4/25** | no Sec is possible (TGA is sense). Bakta at code 4/25 writes W/G there, never `U`. A precomputed protein with `U` gives `codon_mismatch` (5.1), no site and no flag. Variants `rare`-skip the U |
| **Pyl (O over TAG)** | the same in all three tables |
| **alternative starts** | `annotation.start_codons` (default ATG, GTG, TTG) is not derived from the table. NCBI lists more starts for code 4 (TTA, CTG, ATT, ATC, ATA), but Bakta's gene caller (Prodigal rules) uses ATG/GTG/TTG at every table. A user may extend `start_codons`; the initiator is read as M at every table, as today |
| **six-frame index reused across tables** | refused with the rebuild message (example 4). `-resume` re-runs the step anyway because the config slice changes |
| a code-25 genome run at table 4 (or the reverse) | not detectable from the genome: TGA is sense in both. Peptides with G (or W) at TGA codons map only through the tolerant route with one difference (`strain_difference` under `related`). Documented in the user guide; possible post-mapping corroboration goes to the backlog (D11) |
| a plasmid or phage that uses a different code from its chromosome | not supported: one table per run (§4 config). Documented |
| user GFF3/GBFF `--annotation` with `/transl_table=4` | Bakta re-calls CDSs through `--regions` with the configured table. The GBFF qualifier is not read (no change) |
| geNOMAD | picks its own code (its `genetic_code` column); GIPPA does not pass the table, and the column is not compared |
| antiSMASH | reads Bakta's GBFF, which carries `/transl_table=<t>`; no change |
| strain SNV TGA → TGG under 4 | synonymous (W → W), no `stop_codon` flag; under 11 it is a stop loss (M7-16 rule 20) |
| `gippa demo` | table 11; unchanged |
| a format-3 index without `translation_table` in meta (hand-edited) | `ValueError` "unsupported sixframe_index format" (the key is required in format 3) |

## Acceptance criteria

1. **Table 11 is unchanged.** On the gippasim pilot (`tests/integration/test_m10_05_sim_pilot.py`
   inputs) at table 11, every table under `tables/` and `mapping/` is byte-identical to the
   output before M11-16. The documented changes are (the Bakta argv is unchanged at 11, D2):
   `sixframe.npz` meta (format 3, `translation_table`), `run_manifest.json`
   `annotation.translation_table` and `output_schema_version`, `preflight.json` (the
   `genome.genetic_code` `rule` text and `value.translation_table`), and the report's Methods
   sentence. The test lists them explicitly.
2. `config.translation_table` is `ok` for 4, 11 and 25, and `error` with exit `INVALID_INPUT` for
   every other value tested (1, 2, 12, 0, −1, 26). The message names the supported set.
3. `config.tool_genetic_code` gives exactly the levels and `value` of examples 7c and 8, and the
   rows of 2.2 one by one. A config with no conflict is `ok`.
4. `genome.genetic_code` matches example 9 and the tie rows. `_decide` at 4/25 is the exact
   complement of 11 on the evaluated domain (`n4 ≥ 20`), on Hypothesis-generated
   `(bp11 ≤ bp4, n4)`. The W11 and W4 texts match the spec word for word (placeholders filled).
   The level is never `error`.
5. `codon_table(t)` equals Biopython's `CodonTable.unambiguous_dna_by_id[t]` (forward table plus
   stops) for t ∈ {4, 11, 25}. Every function in item 1 raises for 1, 2, 12 and 26.
6. `translate_text` and `build_sixframe` agree with a brute-force per-codon reference that uses
   `codon_table(t)`, on Hypothesis sequences of 0–600 nt over `ACGTN`, for each table and both
   strands. Circular wrap is covered.
7. Index format 3 round-trips at each table. The cross-table load raises (example 4), format 2 is
   accepted only at 11, and the message contains "rebuild the six-frame index".
8. Examples 3, 5 and 6 hold through the real `sixframe_index` → `map_exact` → loci steps (step-level
   test with a tiny genome), including the `mapping.stop_sites` event `data.translation_table`.
9. `effective_tblastn_args` and the Nextflow closure give the same string on examples 7a–7d plus 6
   more cases (an nf-test parity table). The Bakta and miniprot argv parity tests pass at 4, 11
   and 25.
10. `parse_bakta_json` follows example 13. `makedb` at table 4 writes ORFs with W at TGA and
    `makedb_manifest.json` `translation_table: 4`, and refuses a JSON recorded at 11.
11. Codon usage follows example 10. With `trna_counts = {"TCA": 1, "CCA": 1, …}` at table 4,
    TGA's tAI weight is > 0. If codon-bias gives 0, that is reported as a finding with the
    library's behaviour, not patched in GIPPA (D10).
12. Novel ORFs follow example 11. Variants follow example 12, and a code-4 sim strain run has no
    `stop_codon` flag on TGA → TGG sites.
13. `run_manifest.json` validates against the bumped schema with `annotation.translation_table` at
    each table. The report at 4 has the `data-genetic-code="4"` note and the Methods sentence. At
    11 it has the Methods sentence and no note.
14. gippasim: `SimConfig()` output (all existing scenarios, seeds 1–3) is byte-identical to main.
    At `genetic_code=4`, no CDS has an internal TAA/TAG; TGA appears in CDSs at the expected share
    (±3 SE over 20 seeds); no CDS ends in TGA. The same holds at 25 with Gly. `n_selenoproteins=1,
    genetic_code=4` raises.
15. Calibration: every arm of item 9 passes with the existing tolerances. `tolerances.yaml`,
    `ci/gates.yaml` and the seed lists are unchanged in the diff (checked by `ci/protected.txt`).
16. The M. hyorhinis code-4 run completes (exit 0), and `CODE4.md` reports every row of the item 10
    table. Misses are explained.
17. Lint, format and mypy are clean. Tests use no network and no `$GIPPA_DB_DIR`; the optional
    `slow` real-genome tests skip when the files are absent.

## Test plan

- `tests/unit/test_m11_16a_genetic_code.py`: AC 5, item 1 errors.
- `tests/unit/test_m11_16a_preflight.py`: AC 2–4, examples 7c, 8 and 9 (the `_decide` rows with
  the M11-14w fixture generator at codes 4, 11 and 25 and GC 0.28/0.40/0.50/0.70).
- `tests/unit/test_m11_16b_sixframe.py`, `tests/property/test_m11_16b_translate_props.py`: examples
  1, 2 and 4, AC 6–7.
- `tests/unit/test_m11_16b_stops.py`: examples 5–6, `loci` with stop sets on both strands, and the
  `sec_site` / `pyl_site` flags at each table; step-level AC 8 (`tests/unit/steps/`).
- `tests/unit/test_m11_16c_sim.py`: AC 14.
- `tests/unit/test_m11_16d_tools.py`, `workflow/tests/bakta.nf.test`, `miniprot.nf.test`,
  `tblastn.nf.test` (parity table): examples 7–8 and 13, AC 9–10 (parse part), and AC 13
  (manifest and report). The contract test (`tests/contract/`, Slurm, pinned envs) runs miniprot
  `-T 4/25` and tblastn `-db_gencode 4` on a 2-gene fixture whose protein has W at a TGA.
- `tests/unit/test_m11_16e_downstream.py`: examples 10–12, AC 10 (makedb), 11 and 12.
- `tests/calibration/…`: the arms of item 9 (16f).
- `tests/integration/test_m11_16_table11_identity.py`: AC 1 (Slurm, `test_branch.sbatch`).

**Existing tests that pin today's behaviour.** By the role guard, the lead edits these, or the
test-author's PR lists them for the lead:
- `tests/unit/preflight/test_translation_table_11_only.py`: 4 now passes, and the refusal set
  becomes {1, 2, 12, …};
- `tests/unit/test_m11_14w_check.py`: example 9 (table 4 is now evaluated), `MSG_TABLE`, and the
  wording test (W11);
- `tests/unit/test_tblastn.py`, `test_tblastn_nf.py`, `workflow/tests/tblastn.nf.test`,
  `docker_smoke.nf.test`: the default args;
- `tests/unit/test_miniprot.py`, `test_miniprot_nf.py`, `workflow/tests/miniprot.nf.test`:
  `-T 11` → `-T <t>`;
- `tests/unit/test_sixframe.py`, `test_sixframe_review.py`, `tests/helpers/readthrough_index.py`,
  `tests/unit/test_mapping_step_readthrough.py`, `tests/unit/test_readthrough_sites.py`,
  `tests/integration/test_mapping_readthrough.py`: format 2 → 3, and the `table=` keyword on
  `load_index`;
- `tests/unit/test_codon_usage_review.py`: "only 11" raises;
- `tests/unit/test_m11_07_gen_reference.py`: config reference regenerated;
- report snapshots: the Methods sentence.

## Estimate and PR split (each ≤ ~400 code lines; tests in their own commit first, ADR-0007)

Feature-loop order per PR: the test-author commits failing tests and pushes; the implementer
makes them pass and pushes; review and fixes follow; then merge through the queue. MVP-first:
16a → 16b → 16d gives an end-to-end code-4 run (Bakta, index, mapping) before 16e/16f polish the
downstream modules.

| PR | Scope | Code | Tests | Depends on | Khoa |
|---|---|---|---|---|---|
| 16a | ADR-0067; `gippa/model/genetic_code.py`; preflight 2.1–2.3 (both warnings, the tool-args check); schema description; the lead updates the #413 and M11-14w tests | 150 | 300 | none | ADR acceptance (with D1) |
| 16b | six-frame per-table LUTs, `translate`/`translate_text`, `SixFrameIndex.translation_table`, format 3 + refusal; `readthrough` and `loci` stop sets and recoded codons; `steps/mapping` wiring | 230 | 320 | 16a | no |
| 16c | gippasim `genetic_code`, `tga_sense_share`, per-table usage, stops and JSON; M11-14w fixture code 25 (implementer tooling; test-author tests) | 220 | 180 | 16a (module names only; gippasim imports nothing from gippa) | no |
| 16d | Bakta `--translation-table`, miniprot `-T`, `effective_tblastn_args` + Nextflow closure, the defaults.yaml tblastn change (`Gate-Change-Approved: ADR-0067`), `parse_bakta_json` table check, `annotation_meta`, manifest field + schema bump, report note and Methods sentence | 230 | 320 | 16a, 16b | D1 (defaults.yaml) |
| 16e | codon usage, covariates tests, novel outputs/evidence, variants (flags, vcf_records, align, pileup), makedb | 240 | 330 | 16b | no |
| 16f | calibration arms (item 9); statistician review; table-11 identity integration test (AC 1) | 0 | 300 | 16c, 16d, 16e | no |
| 16g | docs: user guide "Genetic code" section (inputs/config), `docs/methods/` note on the 4-vs-25 limit, config reference regenerated, backlog entry removed; validation run and `CODE4.md` (item 10) | 0 (≈ 250 prose) | 0 | all | no |

Total ≈ 1,070 code lines + ≈ 1,750 test lines, plus ≈ 120 Nextflow/config lines in 16d (in its
230), plus prose. 16b and 16c can run in parallel, and so can 16d and 16e. Each PR adds a
`changelog.d` fragment (16a and 16d `feature`, the others `misc`).

## Lead decisions needed

Each has a recommended default. Under "don't block on Khoa", the default is used unless the lead
says otherwise. D1 needs Khoa because it touches a protected file.

- **D1 (Khoa: protected `defaults.yaml`) — the tblastn default.** Change the default
  `tblastn_args` to `"-seg no -evalue 1e-5"` and derive `-db_gencode` from the table (3.3).
  **Decided (Khoa, 2026-10-28): yes (option A)**, approved with ADR-0067; 16d makes the change
  with `Gate-Change-Approved: ADR-0067`. *Fallback (no longer needed):* keep the
  default string. `effective_tblastn_args` then replaces a `-db_gencode` value with the table, but
  only when the string equals the shipped default exactly. Any other explicit mismatch stays a
  preflight error. This is reversible, and it is removed once D1 is approved.
- **D2 — Bakta `--translation-table` always, or only when ≠ 11.** Always passing it changes the
  Bakta task hash once, so every cached Bakta run (the most expensive step) is re-run after the
  upgrade. ~~Recommended: always.~~ **Decided (Khoa, 2026-10-28): only when ≠ 11**; table 11 is the
  default when no flag is given, and code-11 caches stay valid.
- **D3 — format-2 indexes.** **Recommended: accept format 2 as table 11** (it can only have been
  built at 11) and refuse it at 4/25. The alternative is to refuse format 2 always, like format 1,
  but format 1 was refused because it lacked data and format 2 lacks nothing at 11.
- **D4 (`stats`) — the reverse warning W4 at the same 0.55 cut.** **Recommended: yes**, with the
  statistician confirming on the M11-14w 61-assembly panel that no code-4/25 assembly is ≥ 0.55
  (as drafts too). If the statistician wants a margin, a separate constant (for example W4 only
  at ≥ 0.60) is the alternative. It gives no false W4 on that panel either.
- **D5 — W11 level.** Keep `warn` (M11-14w Q2: `--strict` → exit 3). Making it `error` now that
  codes 4 and 25 are supported would stop users with a code-11 genome that scores low; the panel has
  none, but pseudogene-rich genomes reach 0.64. **Recommended: warn.**
- **D6 — a `U` under 4/25 in rule 21(e).** **Recommended: count under the existing
  `codon_mismatch`.** This means no event-schema change, and "U must sit over the table's
  recoding codon" covers it. The alternative is a new `recoded_skipped.not_recoded_in_table` key,
  which is clearer in logs.
- **D7 — validation layer.** **Recommended: preflight only** (2.1), as #413 does, with the schema
  description listing the values. A pydantic `Literal[4, 11, 25]` would fail at config load with a
  traceback-style message before the preflight table.
- **D8 — report.** **Recommended:** the Methods sentence on every run (one snapshot update), and
  the summary note only when ≠ 11. The alternatives are both only when ≠ 11 (no snapshot change
  at 11), or both always.
- **D9 — the precomputed JSON table check.** **Recommended: in `parse_bakta` and `makedb` only.**
  Preflight would have to parse a JSON of up to about 100 MB. The alternative is a streaming read
  of the `genome` block in preflight, which fails earlier but needs a partial JSON parser.
- **D10 — tAI under code 4.** If codon-bias gives TGA no weight from anticodon TCA, **recommended:**
  record the finding in `docs/methods/` and the backlog and keep codon-bias as is (R5, no
  patching). The tAI column is then documented as approximate at code 4.
- **D11 — telling code 4 from 25.** Out of scope. **Recommended:** a backlog item for a
  post-mapping corroboration: W versus G enrichment among tolerant one-difference hits at TGA
  codons. The user guide says the choice comes from taxonomy.
- **D12 — the item 10 thresholds.** These are validation expectations, not CI gates.
  **Recommended: report-only.** A miss is a finding for the lead, not a merge block for 16a–16f.
- **D13 — `translate11*` names.** **Recommended: keep them as aliases** (`table=11`) to limit the
  diff, and pass the table explicitly at every table-dependent call site. The alternative is to
  rename every call site and delete the aliases (≈ 40 more lines, cleaner).
- **D14 — calibration tiering.** **Recommended:** the table in item 9: code 4 in the fast tier for
  exact, tolerant and chance-hit recall, and all arms in `calibration_full`. If the fast tier is
  over budget, move the table-4 tolerant arm to full.
- **D15 — ADR number.** **Recommended: ADR-0067** (next free on 2026-10-28; no open PR uses it).
  Its status is `accepted-by-khoa` once Khoa approves D1. Until then, the parts other than D1 rest
  on Khoa's 2026-10-28 decision recorded in it.

## Gaps found while writing (questions for the lead)

- Q-L1: `gippasim` sense-codon defaults for code 4 (`tga_sense_share = 0.6` of Trp) and code 25
  (0.05 of Gly) are placeholders. Official per-genome codon usage (for example the RefSeq CDSs of
  M. hyorhinis and of Ca. Absconditicoccus praedator, the code-25 genome in the M11-14w panel) would
  ground them. The validation agent can compute them from official NCBI files in 16g, and 16c can
  update the defaults before 16f records any calibration result.
- Q-L2: `annotation.start_codons` stays independent of the table (edge cases). Should a code-4
  run warn when a user adds NCBI code-4 starts that Bakta never uses? **Proposed: no.**
- Q-L3: the VCF header carries no genetic code. A `##gippa_translation_table=<t>` line would make
  `AA_REF` unambiguous for downstream tools, but it changes the VCF snapshot at 11.
  **Proposed: add it only when ≠ 11**, in 16e.
- Q-L4: miniprot's accepted tables are taken from its documentation (1–33 except 7–8 and 17–20).
  The 16d contract test confirms 4 and 25 on the pinned 0.18 build. If it fails, the protein
  route at 4/25 needs a decision (skip miniprot, keep tblastn).
- Q-L5: Bakta's tRNA step (tRNAscan-SE) may not report tRNA-Trp(UCA) at code 4 (the code-11 run has
  no `TCA` anticodon). The item 10 run reports it. If it is missing, tAI at code 4 misses one
  tRNA (see D10).
- Q-L6: `gippa compare` across two runs with different tables. **Proposed:** a compare warning
  "runs use different genetic codes" from the manifests; it is not in this issue (the backlog
  gets it unless the lead adds it to 16d).
- Q-L7: the PR #411 branch `validation/more-datasets` is still open and its README says "GIAP 1.0
  assumes table 11". 16g updates that sentence there, or after #411 merges.
- Q-L8: the backlog's "628 unmapped organism peptides, all with W" does not match the 1,524 of
  8,163 recomputed above (a different peptide set). 16g re-derives the backlog number from the
  #411 notes, or replaces it with the recomputed one.
