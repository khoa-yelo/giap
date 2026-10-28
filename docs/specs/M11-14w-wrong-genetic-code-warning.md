# M11-14w: warn when a default-11 genome probably uses genetic code 4

Milestone M11 · pass 2 · depends on: #413 (preflight refuses `annotation.translation_table != 11`),
M1 preflight (`gippa.preflight`, ADR-0005 item 3), MVP-03 (manifest `warnings`), M9-16 (report
summary section) · estimate ≈ 120 code lines + ≈ 220 test lines, one PR · **`stats`**: the
statistician reviewed the rule and set the 0.55 threshold (#419 review; a classification
cut-off with a false-warning rate) · roles: implementer (`src/gippa/preflight/`, `src/gippa/report/`), test-author
(`tests/unit/test_m11_14w_*.py`) · decision: Khoa, 2026-10-27, Q-K30.

## Spec sections implemented

- §7.0 Preflight ("every evaluated rule yields one check"; `warn` never stops a run; `--strict`
  turns warnings into exit 3, M1-12), §6.5 (0-based half-open), §9.2 manifest `warnings`, R1, R3.
- Q-K30 (Khoa, 2026-10-27): GIPPA 1.0 supports genetic code 11 only; code 4 is a later release
  (backlog "Genetic code 4 support (post-1.0)"). For 1.0, warn when a genome left at the default
  11 probably uses code 4 (Mycoplasma, Spiroplasma, Mesoplasma, Ureaplasma: TGA = Trp).

## The gap

M. hyorhinis SK76 (GCF_000313635.1, PR #411) ran at table 11 with no warning. Bakta called 1,474
features against about 750 in RefSeq, because genes were cut at TGA. The median protein was 90 aa,
against about 287 in the other datasets. All 628 unmapped organism peptides contain W. There were
246 `strain_difference` and 354 `nterm_extension_evidence` flags, and all were artefacts. #413
only catches a table that the user sets explicitly. A user who keeps the default gets wrong
results with no warning.

## Signal (genome only, so it runs in preflight before Bakta)

In a code-4 genome, genes are long open reading frames when TGA is read as a sense codon. With TGA
read as a stop, they break into short pieces. Trp is about 1 % of codons, and most Trp codons are
TGA, so the chance that a 300-codon stretch has no TGA is about e^-3, or 5 %. In a code-11 genome,
genes contain no in-frame TGA, so the long ORFs barely change between the two readings. The rule
compares the two readings directly.

Definitions, for every genome record, both strands and all three frames per strand:
- A **run** is a maximal stretch of whole codons that contains no stop of the given stop set and
  no codon with a non-ACGT base. The run's coordinates exclude the stop codon. A run may end at the
  end of the record ("open"), and open runs count. Every record is treated as linear. A circular
  replicon loses at most one wrap-around run per frame (six per record), which is negligible and
  keeps the scan topology-free. The sequence is upper-cased first.
- **Minus strand**: scan the reverse complement. A run at `[a, b)` there maps to `[n-b, n-a)` on
  the record (n = record length).
- `bp11` = the number of genome positions covered by at least one run of **≥ 300 codons** under
  `S11 = {TAA, TAG, TGA}`. Each position is counted once (union over the six frames, per record,
  summed over records).
- `bp4` = the same under `S4 = {TAA, TAG}`. `n4` = the number of runs of ≥ 300 codons under `S4`.
- `ratio = bp11 / bp4`.

**Rule** (`rule_id = "genome.genetic_code"`, `input = "genome"`):

| Condition | Level | Message |
|---|---|---|
| `annotation.translation_table != 11` | `info` | `not evaluated: annotation.translation_table is not 11` (`config.translation_table` already errors) |
| no records, or `n4 < 20` | `info` | `not evaluated: fewer than 20 open reading frames of at least 300 codons` |
| `20 * bp11 < 11 * bp4` (ratio < 0.55, exact integer comparison) | `warn` | the wording below |
| otherwise | `ok` | `""` |

`value` (JSON) = `{"min_codons": 300, "bp_tga_stop": bp11, "bp_tga_sense": bp4, "n_orfs_tga_sense": n4, "ratio": ratio | null, "threshold": 0.55}`.
`ratio` is `null` when `bp4 == 0`. Every `S11` run lies inside an `S4` run, so `bp11 <= bp4` and
`ratio` is in [0, 1]. When the rule is evaluated (`n4 >= 20`), `bp4 >= 20 * 900 > 0`, so `ratio` is
never `null` in an `ok` or `warn` check. Floats are rounded by the existing `_jsonable` writer. The
decision always uses the exact integer comparison, never the rounded float.

`rule` column text: `long ORFs (>= 300 codons, six frames): bp with TGA as stop / bp with TGA as
sense < 0.55 -> warn (needs >= 20 such ORFs)`.

**Warning wording** (`message`; `{q11}` = bp11/genome bp and `{q4}` = bp4/genome bp, as whole
percentages; `{ratio}` = `(100 * bp11) // bp4` shown as `0.NN`, truncated and not rounded, so a
warning never prints `0.55 < 0.55`):

> This genome probably uses genetic code 4 (TGA = Trp; Mycoplasma, Spiroplasma, Mesoplasma,
> Ureaplasma) or another code in which TGA is not a stop. Long open reading frames cover {q4} % of
> the genome when TGA is read as a codon but only {q11} % when it is a stop (ratio {ratio} < 0.55).
> GIPPA 1.0 supports genetic code 11 only: genes will be cut at TGA, proteins and peptides that
> contain W will not map, and strain_difference and nterm_extension_evidence flags will be
> artefacts. This warning does not change any result.

### Thresholds and why

| Genome (local copies of the official NCBI files) | GC | bp11/G | bp4/G | ratio | n4 |
|---|---|---|---|---|---|
| M. hyorhinis SK76 GCF_000313635.1 (code 4) | 0.26 | 0.154 | 0.679 | **0.227** | 352 |
| E. coli K-12 GCA_000005845.2 | 0.51 | 0.655 | 0.744 | 0.881 | 3,187 |
| G. sulfurreducens GCF_000007985.2 | 0.61 | 0.790 | 0.961 | 0.823 | 6,993 |
| M. smegmatis GCF_000015005.1 | 0.67 | 0.862 | 0.998 | 0.864 | 12,873 |
| synthetic code 4 (GC 0.28, 1 % TGA-Trp) | 0.28 | 0.086 | 0.973 | 0.088 | 300 |
| synthetic code 4 (GC 0.40, 0.5 % TGA-Trp; M. pneumoniae-like) | 0.40 | 0.315 | 0.975 | 0.323 | 302 |
| synthetic code 11 (GC 0.28 / 0.50 / 0.70) | | | | 0.997 / 0.997 / 0.991 | 300–556 |

Extended evidence (statistician, #419 review, 2026-10-27): 61 official NCBI assemblies (60 RefSeq
and the GenBank E. coli K-12 copy above; md5 checked against `md5checksums.txt`; accessions, md5s and the full table are in the #419
review comment). Code as annotated by RefSeq (`transl_table` in the GBFF) or as known for the
lineage.

| Class (assemblies) | ratio range | nearest the cut |
|---|---|---|
| code 4 (21: Mycoplasma, Mycoplasmoides, Mycoplasmopsis, Mesoplasma, Spiroplasma, Ureaplasma, Hodgkinia MAG, Nasuia; 4 drafts of 6–165 contigs) | 0.201–0.470 | Nasuia 0.470 (GC 0.16, 113 kb), M. pneumoniae draft 0.384, M. genitalium 0.370, M. pneumoniae M129 0.358 |
| code 25 (1: Ca. Absconditicoccus praedator) | 0.006 | |
| code 11, low GC ≤ 0.35 (12, incl. Carsonella GC 0.17, Sulcia, Buchnera, Phytoplasma, Acholeplasma, Borrelia with 22 replicons) | 0.959–0.999 | M. jannaschii 0.959 |
| code 11, GC 0.35–0.60 (8, incl. pseudogene-rich Sodalis, M. leprae, S. symbiotica) | 0.679–0.939 | Sodalis glossinidius 0.679, M. leprae 0.717 |
| code 11, GC > 0.60 (19, incl. 3 drafts of 107–167 contigs) | 0.786–0.977 | M. tuberculosis draft 0.786, H37Rv 0.798 |
| draft stress test (10 genomes cut into exponential contigs, mean 1–50 kb; 100 N every 10 kb) | code 4 ≤ 0.477, code 11 ≥ 0.641 | Nasuia 0.477; Sodalis 0.641 at mean 1 kb |
| synthetic fixture, 5 seeds each (spec generator) | code 4 GC 0.28: 0.083–0.106; GC 0.40: 0.254–0.364; code 11: 0.988–0.999 | |

No assembly is misclassified at 0.55 (nor at 0.5). Lower-casing the input and inserting N gaps do
not change the decision.

The spec-writer computed these on 2026-10-27 with a throw-away script that implements the
definitions above. The synthetic genomes have 300 genes of 400 codons, 60 bp random spacers, half
on each strand, seed 7.
- **300 codons.** At 100 codons the code-4 genome stays at 0.68, because a 100-codon stretch has
  no TGA about 40 % of the time. At 300 codons the classes separate. Random ORFs of ≥ 300 codons
  are rare at any GC under `S11`. GC-rich antisense "shadow" ORFs inflate `bp4`, but they cover
  frames that already sit opposite genes, so the ratio stays near 0.8 or above in complete
  high-GC genomes (M. tuberculosis 0.798, M. smegmatis 0.864, Streptomyces coelicolor 0.931).
- **0.55** (statistician, #419; replaces the drafted 0.5). On 61 official NCBI assemblies
  (table above), every code-4/25 genome is ≤ 0.470 and every code-11 genome is ≥ 0.679. Three code-4
  genomes sit above 0.35: M. pneumoniae 0.358 (a 165-contig draft: 0.384), M. genitalium 0.370 and
  the 113 kb endosymbiont Candidatus Nasuia 0.470, only 0.03 below 0.5. The low code-11 tail comes
  from pseudogene-rich genomes (Sodalis glossinidius 0.679, M. leprae 0.717) and high-GC
  genomes (M. tuberculosis 0.798). Cutting complete genomes into draft-like contigs (exponential
  lengths, mean 1 kb, contigs < 200 bp dropped) moves code 4 to at most 0.477 and code 11 to at
  least 0.641 (Sodalis). 0.55 is the round value closest to the centre of that gap (0.56): it
  leaves 0.07 to the highest code-4 value and 0.09 to the lowest code-11 value (0.13 on complete
  assemblies). A missed code-4 genome gives silently wrong results; a false warning only prints a
  message (or exits 3 under `--strict`), so the cut is not moved further toward code 11.
- **n4 ≥ 20.** A ratio built on fewer than 20 long ORFs moves by 5 % or more per ORF. A typical
  plasmid or phage is below 20, so the check is not evaluated there; the warning is aimed at
  chromosomes. The smallest genomes tested pass the gate: Nasuia (113 kb) n4 = 41, Carsonella
  (160 kb) n4 = 75, Sulcia (191 kb) n4 = 94.

Constants, not config keys (no schema change): `MIN_CODONS = 300`, `RATIO_WARN = 0.55` (decided as `20 * bp11 < 11 * bp4`),
`MIN_LONG_ORFS = 20`, in `gippa/preflight/genetic_code.py`.

Rejected as the primary signal (they need annotation or mapping, which come after the expensive
steps): the median CDS length, the TGA-terminated CDS followed in frame by a same-strand CDS, and
W enrichment among unmapped peptides. They stay as possible post-annotation corroboration in the
backlog entry.

## Where it runs and where it surfaces

1. `gippa.preflight.genetic_code.genetic_code_check(seqs: Sequence[str], table: int) -> Check`.
   The function is pure (no I/O, no randomness, numpy only, no new dependency). `run_preflight`
   calls it on the records it has already read (`read_genome_fasta`) and puts the check directly
   after `*checks_genome`. Helper `long_orf_cover(seq, stops, min_codons) -> tuple[int, int]`
   returns `(bp covered, number of runs)` for one record, so tests can use a small `min_codons`.
2. **Preflight table**: `preflight.json` `checks[]` and the `preflight.md` table, through the
   existing writers. There is no new writer code.
3. **Manifest**: `run_manifest.json` `warnings` contains `preflight:genome.genetic_code`, through
   the existing `gippa.manifest._warnings`. There is no new manifest code, only a test.
4. **Report**: when the run's `preflight.json` has `genome.genetic_code` at `warn`, the `summary`
   section starts with `<div class="alert" data-warning="genome.genetic_code">` holding the
   message (escaped). It is read through `RunView.preflight()`. When there is no warning, nothing
   is shown, and the HTML is byte-identical to the HTML without this feature.
5. It **never errors and never changes results**. The level is only `ok`, `info` or `warn`, and
   `exit_code` is 0. `NormalisedInputs` (including `cfg.annotation.translation_table`) is
   unchanged. No step reads the check. Under `--strict`, a `warn` gives exit 3, as every
   preflight warning does (M1-12; see Q2).
6. The check also runs with `--annotation` (a precomputed Bakta JSON). The downstream steps still
   assume table 11, so the warning still applies.

## Examples (input → expected output)

`long_orf_cover` with an explicit `min_codons`. Only the named frame is listed; the other frames
are fixed by the definitions, and the test computes them with a brute-force reference.

1. `s = "ATG" + "GCC"*8 + "TGA" + "GCC"*8 + "TAA"` (57 nt), + strand frame 0, `min_codons=10`:
   under `S11`, runs of 9 codons `[0,27)` and 8 codons `[30,54)`, so neither qualifies. Under
   `S4`, one run of 18 codons `[0,54)` that qualifies (the TGA is read through and TAA ends it).
2. Minus strand: `revcomp(s)` gives the same frame on the − strand. The `S4` run is `[3,57)` in
   record coordinates, and there is none under `S11`.
3. `"ATG" + "GCC"*4 + "NNN" + "GCC"*4 + "TAA"`, frame 0, `min_codons=4`: runs `[0,15)` and
   `[18,30)` under both sets. The NNN codon breaks the run and is not covered.
4. Open end: `"ATGA"`, frame 0, `min_codons=1`: one open run `[0,3)` under `S4`; under `S11`
   also `[0,3)` (ATG; the trailing `A` is not a whole codon).
5. Empty record `""` gives `(0, 0)`.

`genetic_code_check` (default constants):

6. Synthetic code-4 genome (test fixture, below): `warn`, `value["ratio"] < 0.55`,
   `n_orfs_tga_sense ≥ 20`, and the message starts with "This genome probably uses genetic code 4".
7. Synthetic code-11 genome (same generator, no TGA inside genes): `ok`, `ratio > 0.9`.
8. Decision boundary on the integers: `bp11 = 550, bp4 = 1000` (ratio exactly 0.55) gives `ok`;
   `bp11 = 549, bp4 = 1000` gives `warn`, and its message shows `ratio 0.54` (truncated).
   These are tested through a small internal `_decide(bp11, bp4, n4)`. With `n4 = 19` the result is `info` whatever the ratio, and with
   `n4 = 20` the rule is evaluated.
9. `table = 4` gives `info` ("not evaluated: annotation.translation_table is not 11"), whatever
   the sequence.
10. A 50 kb plasmid-only input with fewer than 20 long ORFs gives `info`. A single poly-A record of
    10 kb has no stops in any frame, so it has six open runs and `n4 = 6`; that gives `info`, not
    `ok`.

## Edge cases

| Case | Expected |
|---|---|
| no records / parse error | `info` (not evaluated); `genome.parse` already errors |
| all-N or only non-ACGT | `n4 = 0` → `info`; `ratio = null` |
| tiny genome (< 20 long ORFs) | `info` |
| circular vs linear | identical result (topology ignored); test: same sequence, `complete=true` and `false` |
| minus-strand-only genes | same decision as the plus-strand mirror (test: the revcomp of the whole genome gives the same `bp11`, `bp4`, `n4`) |
| lowercase sequence | identical to upper case |
| many contigs | per-record scan, summed. Runs never cross records |
| ratio exactly 0.55 | `ok` (strict `<` on integers) |
| selenoprotein TGA (code 11) | a few genes; the ratio is not affected materially (no test) |
| `--annotation` given | evaluated as usual |
| translation_table ≠ 11 | `info`; the run already stops on `config.translation_table` |

## Acceptance criteria

1. `genome.genetic_code` appears exactly once in every `preflight.json` that has a genome, with
   `input = "genome"` and `level ∈ {ok, info, warn}`. It is never `error`, and `exit_code` is 0.
2. Examples 1–10 hold. `long_orf_cover` agrees with a brute-force reference (pure Python, codon by
   codon) on Hypothesis-generated sequences of 0–600 nt over `ACGTN`, for `min_codons ∈ {1..5}`
   and both stop sets.
3. The synthetic code-4 fixture warns and the code-11 fixture does not, for 5 fixture seeds each
   and at GC 0.28 and 0.40 (code 4) and GC 0.28, 0.50 and 0.70 (code 11).
4. The result is deterministic: the same input gives the same `Check`. Reverse-complementing every
   record, or reordering the records, leaves `value` unchanged.
5. A non-strict `run_preflight` on the code-4 fixture has the same exit code and an equal
   `NormalisedInputs` to the same run with the check monkeypatched to `ok`. Only `checks` differ.
6. `build_manifest` with a preflight that has the warning lists `preflight:genome.genetic_code`
   in `warnings`.
7. The report on a run view whose preflight has the warning contains the `data-warning` alert in
   `summary`. Without the warning, the HTML has no such element (snapshot unchanged).
8. Performance: on a 10 Mb synthetic genome on one core, the check takes ≤ 5 s
   (`benchmarks/`, reported only; not a CI gate).
9. Tests use no network, no `$GIPPA_DB_DIR` and no `$GIPPA_ROOT` (except the optional `slow`
   validation test below, which skips when the files are absent). Lint, format and mypy are clean.

## Test plan

- `tests/unit/test_m11_14w_long_orfs.py`: examples 1–5, the brute-force Hypothesis property,
  strand symmetry, lowercase input, and the N break.
- `tests/unit/test_m11_14w_check.py`: examples 6–10, the `_decide` boundary, the `info` paths,
  AC 3 seeds, AC 4 and AC 5 (`run_preflight` with `probe` stubbed, as the existing preflight tests
  do).
- Fixture generator (in the test module or `tests/sim/gippasim/`, implementer's choice per AGENTS
  "gippasim is tooling"): `n_genes = 300`, 400 codons each (ATG, 398 sense codons, TAA), drawn from
  base frequencies at GC `g` excluding TAA/TAG/TGA. For code 4, each sense codon is replaced by
  TGA with probability `r` (0.01, or 0.005 at GC 0.40). Genes are separated by 60 bp random
  spacers and half are reverse-complemented. Use `numpy.random.default_rng(seed)`.
- `tests/unit/test_m11_14w_surfaces.py`: AC 6 (manifest) and AC 7 (report, with a minimal
  `RunView` fixture as in the M9-16 tests).
- Optional validation (`@pytest.mark.slow`, skipped when the files are absent; never in CI): the
  four real genomes in the table, read from `$GIPPA_ROOT/more_datasets/*/raw/genome/` and
  `$GIPPA_ROOT/data/genomes/k12_mg1655/`, with the ratios to ±0.01.

## Questions for the lead

- Q1 (`stats`, validation): the real-genome evidence is one code-4 and three code-11 genomes, all
  above GC 0.25 apart from M. hyorhinis. Before merge, should the validation agent compute the
  ratio on official NCBI copies of low-GC code-11 genomes (Staphylococcus aureus, Clostridioides
  difficile, Buchnera APS, Candidatus Carsonella at 17 % GC) and on more code-4 genomes
  (M. pneumoniae at 40 % GC, Spiroplasma, Ureaplasma)? **Lead: yes. Done in the #419 statistician review: three code-4 genomes fell in [0.35, 0.65],
  none was misclassified, and the threshold moved from 0.5 to 0.55 (see Thresholds).** Proposed was: yes; adjust 0.5 only if a
  genome falls between 0.35 and 0.65.
- Q2: under `--strict`, this warning makes preflight exit 3, like every warning. Should it be
  exempt? **Proposed: no.** A strict user should not run a code-4 genome at table 11.
- Q3: the report shows only this warning. A general "Run warnings" list (every manifest
  `warnings` entry) would also cover `decoy_unscored` (backlog #359). **Proposed: this warning
  only for 1.0; the general list goes to the backlog.**
- Q4: there are no config keys for the thresholds, so a false warning cannot be silenced. That is
  acceptable for a warning. **Proposed: keep them as constants.**
