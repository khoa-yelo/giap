# M10-07: §15.5 metric 8, A1's free entrapment (`validation/pilot/entrapment.py`)

Milestone M10 · **pass 2** · depends on: M10-05 (`metrics.py`, `metrics.json`, the report
template; spec on `m10/spec-05`), M10-03 (the converted A1 table `a1.tsv.gz`), M10-02 (the A1
search settings) · estimate ≈ 230 code lines + ≈ 320 test lines + a methods note · **`stats`**
(the statistician reviews this spec before any test is written) · roles: validation
(`validation/pilot/entrapment.py`, `scripts/slurm/entrapment.sbatch`,
`docs/methods/free-entrapment.md`), test-author (`tests/unit/test_m10_07_*.py`), statistician
(review, the methods note's review) · Khoa: none (the M10-05 Q-K1 rule covers the outputs).

## Spec sections implemented

- §15.5 metric 8, verbatim: "Peptides absent from Nissle's translation (no exact match after the
  I/L and mass-equivalent rules) cannot be true if the genome is complete. With N_E such peptides,
  N_T peptides present, and r the in-silico ratio of K-12-only to shared tryptic peptides
  (computed with the same exclusions, within the searchable length, charge and m/z range), the
  combined estimator of Wen 2025 gives FDP̂ = N_E (1 + 1/r) / (N_T + N_E). These entrapment
  peptides are near-homologs, easier to mismatch than random sequences, so the estimate errs high.
  Only wild-type runs count, because engineered strains may carry K-12-derived constructs.
  Implemented in `validation/pilot/entrapment.py`."
- §8.7: validation measures the real false discovery proportion of the inherited identifications
  (A1's free entrapment).
- §15.7 hypothesis 2 ("A1's genome-inconsistent FDP is at least its entrapment FDP (full
  validation, FDRBench)") stays `deferred`; this issue only places the two numbers side by side
  (question L4).
- R1 and §16: "FDRBench entrapment" is a Later item. This issue therefore adds **no** entrapment
  database and **no** DIA-NN search (see "Does this need a new DIA-NN search?").
- R4: the estimate is labelled as what it is (an estimated upper bound in expectation, with the
  lower bound beside it); its sampling uncertainty is `unestimated` unless the statistician says
  otherwise (question S3).
- Khoa's 2026-10-23 rule and M10-05 Q-K1: nothing derived from the DIA-NN reports is committed.

## The estimator (source, confirmed)

Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U. "Assessment of false discovery rate
control in tandem mass spectrometry analysis using entrapment." *Nature Methods* 22 (2025),
https://doi.org/10.1038/s41592-025-02719-x (open access, CC BY 4.0; the DOI the spec's reference
list gives). The formulas below were read on 2026-10-24 from the publisher's article page that
the DOI resolves to. They are **quoted, not reconstructed from memory**:

- **Equation (1), the "combined" method.** It estimates the FDP among the target + entrapment
  discoveries: FDP̂(T ∪ E) = N_E (1 + 1/r) / (N_T + N_E). Here r is "the effective ratio of the
  entrapment to original target database size", and N_T and N_E are "the number of original
  target and entrapment discoveries". Supplementary Note 2.2 proves that, "under an assumption
  analogous to the equal-chance assumption that TDC relies on", (1) "provides an estimated upper
  bound, that is, on average it overestimates the true FDP". With r = 1 it reduces to Elias and
  Gygi's concatenated estimate. This is §15.5's formula, unchanged.
- **Equation (2), the lower bound.** N_E / (N_T + N_E). Assuming that every entrapment discovery
  is false, it "represents a lower bound on the FDP" (Supplementary Note 2.1). The paper names
  using it as evidence of control "incorrect"; used as evidence of *failure*, it is correct.
- **Equation (3), the "sample" method.** N_E (1/r) / N_T. The paper calls it "inherently flawed".
  It is **not** computed.
- **The paired method (4)** needs one-to-one shuffled pairs (r = 1), which a free entrapment does
  not have. It is not applicable.
- **Interpretation (paper, Fig. 1b).** An upper bound at or below the threshold is empirical
  evidence of control, a lower bound above it is evidence of failure, and anything in between is
  inconclusive. The paper also notes that shuffled entrapments cannot model errors due to homologs,
  and that foreign species raise the question of evolutionary distance. Our entrapment set is a
  related strain's peptides, which is the near-homolog case §15.5 calls "errs high".

The paper's averaged upper bound is a statement about expectation. A single estimate can fall below
the true FDP. The report says "estimated upper bound", never "upper bound".

**Mapping to the pilot.** The A1 search database is UniProt K-12 UP000000625 + contaminants.
Against Nissle spectra, its peptides split into the "original targets" (K-12 peptides present in
Nissle's translation) and the "entrapment" (K-12-only peptides). The entrapment is free because it
is already inside A1's database. DIA-NN controlled its FDR over the whole K-12 database, so the
combined method's quantity (the FDP among T ∪ E discoveries) is exactly the FDP of A1's
accepted set.

## Does this need a new DIA-NN search?

**No.** Metric 8 is computed from A1's existing search (M10-02 job 1879955) and a pure in-silico
digest. An extra A2 search with an entrapment-extended library (about 15 min on Slurm) would be a
*constructed* entrapment experiment (shuffled K-12 or Nissle sequences, or a foreign proteome).
That is FDRBench-style validation, listed in §16 as Later ("FDRBench entrapment"), so R1 forbids it
in version 1. If the lead or Khoa wants it anyway (question L1), it would be built as follows:

- **Database:** a shuffled entrapment from the A2 database `strain_db.faa`. r = 1 for the paired
  method; the paper's protocol is to shuffle each peptide with its C-terminal residue fixed, check
  it is distinct from every target, and retry up to 20 times.
- **Tool:** FDRBench, from Noble-Lab's GitHub releases (an allowed §3.2 source; licence to be
  read at the pinned tag).
- **Approvals:** an ADR, and Khoa, because §16 is involved.

A foreign proteome is the weaker choice because of its evolutionary-distance problem. If used, it
would come from `rest.uniprot.org`.

## What exists (main 1f567901; M10-05 spec on `m10/spec-05`, M10-02 on `m10/diann-searches`)

| Piece | State | Used for |
|---|---|---|
| A1 search (M10-02) | DIA-NN 2.3.1 Academia, two-step library-free, job 1879955. `diann_A1.cfg`: `--cut K*,R*`, `--missed-cleavages 1`, `--min-pep-len 7`, `--max-pep-len 30`, `--unimod4` (fixed carbamidomethyl C), `--var-mods 1`, `--var-mod UniMod:35,15.994915,M`, `--min-pr-charge 1`, `--max-pr-charge 4`, `--min-pr-mz 300`, `--max-pr-mz 1800`, `--min-fr-mz 200`, `--max-fr-mz 1800`, `--peptidoforms`, `--qvalue 0.01`, `--report-decoys`; **no** `--met-excision`, **no** MBR (`--reanalyse` was left out for this metric: README) | the searchable space (below) |
| `giap.mapping.digest.digest(protein, enzyme, missed, min_len, max_len)` | `trypsin_p` = DIA-NN's `K*,R*`; keeps the initiator M; 0-based half-open products | the in-silico digest of K-12 and contaminants |
| `giap.mapping.sixframe.build_sixframe`, `giap.mapping.exact.match_exact` | table-11 six-frame index, I → L, 120-nt origin overlap on circular replicons, `X` never matches; Aho–Corasick, every locus | "present" (exact rule) |
| `giap.mapping.tolerant.match_tolerant` | ≥ 10 aa; budget 1 (10–11 aa) or 2 (≥ 12 aa); rows carry `diff_types` (`substitution`, `swap`, `block`, `near_block`, `frameshift`) | "present" (mass-equivalent rule) |
| `giap.mapping.masses` | residue masses (Appendix B.5; `RESIDUE_MASS`, `MOD_DELTA` with `Oxidation` 15.994915 and `Carbamidomethyl` 57.021464) | precursor m/z |
| M10-05 `metrics.py` | `metrics.json` (`schema giap-pilot-metrics/1`), `render`, the template with metric 8 under "Pending in pass 2", the accepted-row rule (target, `q_value ≤ 0.01`, NA not accepted), `metrics.json` `m2` genome-inconsistent count | the A1 accepted set; the hook for `m8`; the cross-check |
| `validation/pilot/inputs/pxd078459_files.tsv` (committed) | the `run` column `WT_1`–`WT_3`, the three selected untreated wild-type runs | the "only wild-type runs" check |
| `validation/pilot/db_a2/strain_db.faa.gz` (committed) | Nissle's Bakta CDS + sORF proteins, initiator M as Bakta writes it | the non-ATG N-terminal rule (L2) |
| `tests/data/genomes/nissle1917.fna.gz` (committed, circular-tagged) | the Nissle genome | the six-frame translation |

### The searchable space (derived from `diann_A1.cfg`)

The space is distinct tryptic peptides of 7–30 aa (`trypsin_p`, ≤ 1 missed cleavage, initiator M
kept), built from the 20 standard residues. Precursors are those peptides with carbamidomethyl on
every C and 0 or 1 Met oxidation (`--var-mods 1`; each oxidised M position is its own peptidoform,
as with `--peptidoforms`), at charge 1–4, with 300 ≤ m/z ≤ 1800. m/z = (M + z·1.007276)/z, where
M is the sum of residue masses + 18.010565.

**At peptide level the charge and m/z bounds remove nothing.** Every 7–30-aa peptide has a mass
between 7 × G + H₂O = 417.16 Da and 30 × W + H₂O = 5600.39 Da. Charges 1, 2, 3 and 4 cover masses
of [298.99, 1798.99], [597.99, 3597.99], [896.98, 5396.98] and [1195.97, 7195.97] Da, and their
union contains that whole interval. `entrapment.py` still applies the bounds (they are parameters,
so a settings change is handled), and a test pins this fact. At **precursor** level, the bounds
decide how many charge states each peptide has (examples 1, 2).

## Inputs (on FarmShare; none committed)

| Name | Default | Role |
|---|---|---|
| `PEP_A1` | `$GIAP_ROOT/pilot/inputs/a1.tsv.gz` (M10-03) | A1 discoveries |
| `K12_FASTA` | `$GIAP_DATA_DIR/proteomes/UP000000625.fasta` | the target + entrapment database |
| `CONTAMINANTS` | `$GIAP_DB_DIR/contaminants/contaminants.fasta` (the file A1 searched; its path is in `command_step1.txt`) | exclusions |
| `GENOME` | `tests/data/genomes/nissle1917.fna.gz` | Nissle's translation |
| `STRAIN_DB` | `validation/pilot/db_a2/strain_db.faa.gz` | the N-terminal rule (L2) |
| `RUNS_MANIFEST` | `validation/pilot/inputs/pxd078459_files.tsv` | wild-type runs |
| `OUT_A1` (optional) | the M10-04 A1 run directory (as M10-05 resolves it) | cross-checks (exact route, genome-inconsistent count) |
| `DIANN_CMD` (optional) | `$GIAP_ROOT/pilot/A1/command_step1.txt` | asserts the digest and precursor settings equal the defaults below (exit 2 otherwise) |

Before computing, `entrapment.py` checks three things:

1. The sha256 of `K12_FASTA` and `CONTAMINANTS` equals the one M10-02's `provenance.json` records,
   when that file is given (`--provenance`). These are the files the search used.
2. Every `run` value of `PEP_A1` is in the manifest's `run` column. A non-wild-type run means the
   q-values were computed jointly with it, and its rows cannot simply be dropped.
3. `GENOME` carries its topology tags.

Any failure exits 2, names the file, and writes nothing.

## Definitions

Notation: `norm(s)` = `s` with I → L. Any residue outside the 20 standard letters (U, X, B, Z, O)
makes the peptide **non-standard**.

**Digest D(F)** of a FASTA F: `{norm(p) | p ∈ digest(protein, "trypsin_p", 1, 7, 30) for each
protein of F, p standard}`. A set of distinct normalised sequences. The I/L collapse matches what
DIA-NN can tell apart by mass (question S5).

**Universe.** U = D(K12_FASTA) \ D(CONTAMINANTS). A K-12 peptide that a contaminant protein also
yields is excluded on both sides, because it could be a true contaminant identification. This is
"the same exclusions" in §15.5.

**Present in Nissle's translation**, written `present(p)`. It is **one function**, applied
unchanged to the in-silico universe and to the discoveries. `p` is present when either:

1. **exact:** `match_exact` returns a target row for `p` on Nissle's six-frame index (I → L; every
   frame, both strands, origin wrap on circular replicons); or `p` is a substring of a
   `STRAIN_DB` protein after I → L (L2: this catches N-terminal products of non-ATG starts, where
   the genome reads V/L/I but the protein, and DIA-NN, read M); or
2. **mass-equivalent:** `match_tolerant`, run with the GIAP defaults (`mapping.tolerant`, ≥ 10 aa),
   returns a target row whose `diff_types` is non-empty and contains only `block`. Those are the
   exact-composition exchanges of §7.2.4, such as N ↔ GG and Q ↔ AG. `near_block`, `swap`,
   `substitution` and `frameshift` rows do **not** make a peptide present (question S2).

Peptides of 7–9 aa can only be present through rule 1, because the tolerant route does not
search them. The asymmetry is the same on both sides, since `present` is one function.

**Classes.** For p ∈ U: `shared` if `present(p)`, else `k12_only`.

**Ratio.** `r = |k12_only| / |shared|` at peptide level (§15.5's definition, the primary).
At precursor level (secondary, S1), `r_prec = Σ_{k12_only} n_prec(p) / Σ_shared n_prec(p)`.
Here `n_prec(p)` counts the precursors of p in the searchable space: for each form (unmodified,
plus one form per M position oxidised), the charges z ∈ {1..4} with 300 ≤ m/z ≤ 1800. Masses
come from `giap.mapping.masses`, with carbamidomethyl added to every C. A literal I or L in the
DIA-NN sequence does not change the mass.

**Discoveries.** These are the A1 accepted target rows, using M10-05's rule: `is_decoy` false and
`q_value ≤ 0.01` (inclusive, the written value; NA is not accepted). Rows from every run are used,
and all runs are wild type by the precondition. At peptide level the unit is `norm(peptide)`,
distinct. Each discovery falls in exactly one bucket:

| Bucket | Rule | Counted in |
|---|---|---|
| `excluded_contaminant` | `norm(peptide)` ∈ D(CONTAMINANTS) | reported only |
| `not_in_universe` | not in D(K12_FASTA) (for example a Met-excised or semi-specific product) | reported; > 0 is a `warnings` entry `digest_mismatch` (the in-silico space does not describe what DIA-NN searched) |
| `target` (N_T) | in U and `shared` | N_T |
| `entrapment` (N_E) | in U and `k12_only` | N_E |

The precursor level (secondary) uses distinct (`modified_peptide`, `charge`) of the accepted rows.
Each precursor takes the bucket of its stripped peptide.

**Estimates**, with α = 0.01, the q cut-off the discoveries were selected at:

- `fdp_combined = N_E (1 + 1/r) / (N_T + N_E)`. It is computed in exact rational arithmetic as
  `N_E · (|k12_only| + |shared|) / (|k12_only| · (N_T + N_E))` and converted to float once, so the
  boundary at α is exact (example 6).
- `fdp_lower = N_E / (N_T + N_E)`.
- `assessment`: `consistent_with_control` if `fdp_combined ≤ α`; `evidence_of_failure` if
  `fdp_lower > α`; else `inconclusive`. These are the paper's three readings. The label is
  descriptive: metric 8 gates nothing (§15.6).
- `ci`: the literal `unestimated` (R4), with the reason "no interval for the combined estimator
  is defined in the source; dependence between peptides of one protein" (question S3).
- `caveats`: three fixed strings, rendered into the report:
  - "near-homolog entrapment: the estimate errs high (§15.5)";
  - "an estimated upper bound in expectation (Wen 2025, Suppl. Note 2.2), not a bound on this
    sample";
  - "assumes the Nissle genome is complete".

**Levels.** The primary is the peptide level (§15.5's N and r). The secondary is the precursor
level (DIA-NN's controlled unit, and the level Wen 2025 applied to DIA-NN). There are also
per-run numbers (descriptive): discoveries are the accepted rows of the run, and r is unchanged.

## Outputs

All under `METRICS_DIR` (default `$GIAP_ROOT/pilot/metrics/`). `--out` inside the checkout exits 2.

| File | Content | Row-level DIA-NN data? |
|---|---|---|
| `entrapment.json` | `{"schema": "giap-pilot-entrapment/1", "spec": <SHA>, "inputs": {name: {path, sha256}}, "settings": {enzyme, missed, min_len, max_len, charges, mz, var_mods, tolerant}, "universe": {"k12_peptides", "excluded_contaminant", "non_standard", "shared", "shared_exact", "shared_mass_equivalent", "shared_strain_db_only", "k12_only", "r", "shared_precursors", "k12_only_precursors", "r_prec"}, "peptide": {N_T, N_E, excluded_contaminant, not_in_universe, fdp_combined, fdp_lower, assessment, ci}, "precursor": {…same…}, "per_run": {run: {…peptide fields…}}, "cross_checks": {…}, "caveats": […], "warnings": […]}` | no (counts and ratios) |
| `entrapment_discoveries.tsv` | one row per discovered peptide: `peptide`, `bucket`, `present_rule` (`exact`, `strain_db`, `mass_equivalent`, or empty) | **yes**: never committed |

JSON follows M10-05's byte rules: sorted keys, two-space indent, `repr` floats, final newline, no
timestamp. The `universe` block is a pure function of public inputs (UniProt, NCBI genome,
contaminants), and the M10-05 Q-K1 answer governs it like everything else.

**Hook into M10-05.** `metrics.py compute` gets an optional `--entrapment METRICS_DIR/entrapment.json`
and copies the `peptide`, `precursor` and `universe` blocks, the caveats and the cross-checks
under `m8`. The template moves metric 8 from "Pending in pass 2" to its own section, with
`fdp_combined`, `fdp_lower`, `assessment`, r, both levels and the caveat sentences. The H2 line
stays `deferred` and shows the A1 genome-inconsistent share and `fdp_combined` side by side, with
no verdict (L4).

## Cross-checks (warnings, not exit codes, unless stated)

1. **Against the run (when `OUT_A1` is given).** Every discovered peptide with a `route = exact`
   target row in A1's `peptide_hits.tsv` is `present` (rule 1). A violation exits 1: two exact
   matchers disagree on one genome.
2. **Against the genome-inconsistent count (when `OUT_A1` is given).** It is M10-05's `m2`, the
   count of A1 non-contaminant query peptides without an exact placement. N_E ≤ that count, because
   an entrapment peptide has no exact match and contaminant handling here only removes more. A
   violation is a `warnings` entry. The difference is reported in two parts: the peptides that
   rule 1's `STRAIN_DB` clause or rule 2 made present, and the placements that are exact but
   rejected.
3. **Digest agreement.** `not_in_universe = 0` is expected, because DIA-NN and GIAP use the same
   digest (M10-02 README "Digest flags = GIAP's digest"). If it is not 0, `warnings` gets
   `digest_mismatch` with the count.
4. **Settings (when `DIANN_CMD` is given).** The flags in the table above equal `settings`.
   Otherwise exit 2.

## Methods note (`docs/methods/free-entrapment.md`, validation role; statistician reviews)

The note must contain:

- the citation and DOI, and equations (1) and (2) as quoted above, with what the paper proves
  about each;
- why the sample method (3) is not used;
- the mapping of T and E onto shared and K-12-only peptides;
- the `present` rules and why one function serves both sides;
- the searchable space and the peptide-level m/z fact;
- the assessment rule;
- the bias table: near-homologs (high), 7–9-aa variants not tolerantly matched (high), equal-chance
  violations in the low direction (false matches attracted to abundant shared proteins: low,
  unquantified), an incomplete genome (high), K-12-derived constructs (excluded by "wild type
  only");
- why no CI is given (S3).

It is written in the house style of `docs/methods/target-decoy.md`, with no pilot numbers in it.

## Examples

Hand-made sequences and counts. None comes from the real reports.

**Shared fixture (examples 1–2).**

- `K12_FASTA` holds three proteins:
  - P1 `MSAGLLKVEPTIDEKAAGGWWR`;
  - P2 `AGNVLTEPWR`;
  - P3 `QQWWLLDR`.
- `CONTAMINANTS` holds C1 `QQWWLLDRPEPTIDES`.
- Nissle's six-frame translation contains `MSAGLLKVEPTLDEKAAGGYWR` and `AGGGVLTEPWR`, and nothing
  else that matches these peptides.
- `STRAIN_DB` adds nothing.

| # | Input | Expected |
|---|---|---|
| 1 | **Digest, universe, r.** The shared fixture | P1 digests to `MSAGLLK` [0,7), `VEPTLDEK` [7,15) (I → L), `AAGGWWR` [15,22), `MSAGLLKVEPTLDEK` [0,15) and `VEPTLDEKAAGGWWR` [7,22); P2 to `AGNVLTEPWR`; P3 to `QQWWLLDR`. D(C1) = {`QQWWLLDR`, `PEPTLDES`, `QQWWLLDRPEPTLDES`} (trypsin_p cuts R‖P). U has 6 peptides (`QQWWLLDR` is `excluded_contaminant`). `shared`: `MSAGLLK`, `VEPTLDEK`, `MSAGLLKVEPTLDEK` (exact) and `AGNVLTEPWR` (mass-equivalent: N ↔ GG, `diff_types = block`, 10 aa, budget 1). `k12_only`: `AAGGWWR` (7 aa, so no tolerant search) and `VEPTLDEKAAGGWWR` (tolerant hit W → Y, `diff_types = substitution`, not present). **r = 2/4 = 0.5.** Precursors per peptide: `MSAGLLK` 4 (z1 719.41201, z2 360.20964; ox z1 735.40693, z2 368.20710), `VEPTLDEK` 3 (z1–z3), `MSAGLLKVEPTLDEK` 8, `AGNVLTEPWR` 3, `AAGGWWR` 2, `VEPTLDEKAAGGWWR` 4. **r_prec = 6/18 = 0.3333…** |
| 1b | Same, but P2 = `GNVTEPWR` (8 aa) and the translation has `GGGVTEPWR` | `GNVTEPWR` is `k12_only`: 8 aa is below the tolerant route's 10, so rule 2 cannot apply. r = 3/3 = 1.0 |
| 2 | **Discoveries.** A1 rows (run, peptide, modified, z, q, decoy; `modified` equals `peptide` where only one sequence is shown). 1: WT_1 `MSAGLLK`/`MSAGLLK`/2/0.001/F. 2: WT_2 `MSAGLLK`/`M(UniMod:35)SAGLLK`/2/0.003/F. 3: WT_1 `VEPTIDEK`/2/0.004/F. 4: WT_3 `AAGGWWR`/2/0.008/F. 5: WT_1 `VEPTIDEKAAGGWWR`/3/0.02/F. 6: WT_2 `QQWWLLDR`/2/0.001/F (`Cont_…`). 7: WT_1 `AGNVLTEPWR`/2/**0.01**/F. 8: WT_2 `KEDITPEPR`/2/0.001/**T** | Peptide level: N_T = 3 (`MSAGLLK`, `VEPTLDEK`, `AGNVLTEPWR`; q = 0.01 is accepted), N_E = 1 (`AAGGWWR`), `excluded_contaminant` 1, `not_in_universe` 0. Row 5 is not accepted and row 8 is a decoy. `fdp_combined` = 1 · 3 / 4 = **0.75**, `fdp_lower` = 0.25, `evidence_of_failure`. Precursor level: N_T = 4 (rows 1, 2, 3, 7), N_E = 1, r_prec = 1/3, `fdp_combined` = 1 · 4 / 5 = 0.8, `fdp_lower` 0.2. Per run: WT_3 has N_T 0, N_E 1, `fdp_combined` 3.0 (above 1: reported as computed, flagged `above_one`), `fdp_lower` 1.0 |
| 3 | **Estimator, inconclusive.** N_T = 4,950, N_E = 20, \|k12_only\| = 25,000, \|shared\| = 100,000 (r = 0.25) | `fdp_combined` = 20 · 125,000 / (25,000 · 4,970) = 100/4,970 = 0.0201207; `fdp_lower` = 20/4,970 = 0.0040241; `inconclusive` |
| 4 | **Consistent with control.** N_T = 9,900, N_E = 10, r = 0.5 (k12_only 50,000, shared 100,000) | `fdp_combined` = 30/9,910 = 0.0030272; `fdp_lower` = 0.0010091; `consistent_with_control` |
| 5 | **Failure.** N_T = 900, N_E = 100, r = 1 | `fdp_combined` 0.2; `fdp_lower` 0.1 > 0.01; `evidence_of_failure` |
| 6 | **Boundary.** N_T = 2,990, N_E = 10, k12_only 1, shared 2 (r = 0.5) | `fdp_combined` = 10 · 3 / (1 · 3,000) = 0.01 exactly (rational arithmetic); `consistent_with_control` (≤ is inclusive). With N_T = 2,989: 30/2,999 = 0.0100033, `inconclusive` |
| 7 | **No entrapment discoveries.** N_T = 5,000, N_E = 0, r = 0.4 | `fdp_combined` 0.0, `fdp_lower` 0.0, `consistent_with_control`, `warnings` `n_entrapment_zero` (the report says "no K-12-only peptide was accepted"; R4: an FDP estimate of 0 with no interval) |
| 8 | **Digest mismatch.** As example 2, plus an accepted target row `SAGLLKVEPTIDEK` (a Met-excised N-terminal product that GIAP's digest never makes) | `not_in_universe` 1, `warnings` `digest_mismatch`; N_T and N_E unchanged |
| 9 | **Non-ATG start (L2).** As example 1, plus K-12 protein P4 `MTEKLLAGGPR` (digest: `MTEKLLAGGPR` 11 aa; `MTEK` and `LLAGGPR` → only `LLAGGPR`, 7 aa, is ≥ 7). The genome reads `VTEKLLAGGPR` (GTG start), and `STRAIN_DB` has `MTEKLLAGGPR…` | `LLAGGPR`: exact, `shared`. `MTEKLLAGGPR`: not exact in the six-frame translation, but a substring of a `STRAIN_DB` protein, so `shared` (`present_rule = strain_db`, counted in `shared_strain_db_only`). Without `STRAIN_DB` (L2 answered "no") it is `k12_only`, because rule 2 sees M ↔ V as a substitution |
| 10 | **Wild type only.** `PEP_A1` has a row with run `WT_ROS_1` | exit 2 naming the run; nothing written |

## Edge cases

- **Circular and linear contigs.** A peptide whose only occurrence crosses the origin of a
  circular replicon is present: the index has the 120-nt overlap. On a contig without a circular
  tag, the same text split across the end is absent. Fixture: a 2-contig genome, one
  circular-tagged with `AGGGVLTEPWR` across its origin, one linear with the same split. Expected:
  `shared`, and `k12_only` when only the linear one holds it.
- **Minus strand.** A peptide that occurs only in a reverse-complement frame is present: both
  strands are indexed. Fixture: example 1's translation placed on the minus strand gives the same
  classes.
- **Empty input.** No accepted target row: N_T = N_E = 0, both FDPs `null` with
  `reason: empty_discoveries`, `assessment: not_evaluable`. An empty `K12_FASTA` exits 2.
- **r undefined.** `|shared| = 0`: r and every estimate are `null` with `reason:
  no_shared_peptides`. `|k12_only| = 0`: r = 0, so 1/r is undefined. `fdp_combined` is `null`
  with `reason: no_entrapment_space`, even when N_E = 0, and `fdp_lower` is still reported (0 when
  N_E = 0). An entrapment discovery with `|k12_only| = 0` cannot happen because N_E ⊆ k12_only;
  an assertion guards it.
- **`fdp_combined` above 1.** With small N, N_E (1 + 1/r) can exceed N_T + N_E (example 2, per
  run). It is reported unclipped and flagged `above_one`; the paper's estimator is not clipped.
- **Ties.** Only distinct counts are taken, so no ranking is involved. One peptide in several runs,
  charges or peptidoforms is one peptide (and several precursors). Two K-12 proteins sharing a
  peptide give one universe element. A peptide that is both exact and mass-equivalent takes
  `present_rule` in the order `exact`, `strain_db`, `mass_equivalent`.
- **I/L.** A K-12 peptide with I and a Nissle text with L at that position is `shared`. Two K-12
  peptides that differ only by I/L are one element of U, and two discoveries that differ only by
  I/L are one discovery (S5).
- **Non-standard residues.** K-12 selenoproteins carry U. Any digest product with U, X, B, Z or O
  is dropped from U and counted in `non_standard`. DIA-NN's handling of U is not documented in the
  README, so the precursor count would be wrong either way. Such a product cannot be a discovery
  of the converted table's peptides; if one is, it counts as `not_in_universe`.
- **q at the cut-off and NA q.** q = 0.01 is accepted (example 2). Target rows with NA q are not
  discoveries, as in M10-05 metric 1. The `q_crossings` superset of M10-05 is not repeated here.
- **Contaminant-labelled but K-12-yielded peptides.** The bucket is decided by sequence (D(C)),
  not by `protein_ids`. A row labelled `Cont_…` whose sequence is not in D(C) is impossible with
  the same digest and is counted under `digest_mismatch`.
- **Decoys.** They are never discoveries and never in U.
- **Several M.** Each single-oxidation position is its own precursor (`MMGGGGGK`: 3 forms × 2
  charges = 6 precursors).

## Acceptance criteria

On CI (synthetic inputs only, built in `tmp_path`; no DIA-NN data, no network):

1. Examples 1–10 give exactly the values above: counts equal, r, r_prec and the FDPs equal to
   1e-12, assessment labels equal, and m/z values equal to 1e-5.
2. **One classification function.** A spy on `present` shows that the universe and the
   discoveries go through the same callable with the same index. A discovery classified `shared`
   is in U's `shared` set, and vice versa, for every example.
3. **Formula.** `fdp_combined` equals `N_E·(k+s)/(k·(N_T+N_E))` from `fractions.Fraction` for 200
   seeded random (N_T, N_E, k, s) with k, s ≥ 1. Example 6 is `consistent_with_control` and its
   N_T − 1 variant is `inconclusive`.
4. **Searchable space.** For every length 7–30 and the residues G and W at both extremes, the
   peptide-level filter keeps the peptide. A test with `--max-pr-charge 1` and `--max-pr-mz 600`
   drops `VEPTLDEK` (z1 930.48), which shows the bounds are applied and not hard-coded away.
   n_prec matches example 1 for each peptide.
5. **Digest equals the search.** `settings` equals the flags parsed from a synthetic
   `command_step1.txt` built from `diann_A1.cfg`. Changing `--missed-cleavages` to 2 in that file
   exits 2.
6. **Present rules.** The edge-case fixtures hold: origin-crossing on a circular contig is shared,
   the linear split is not, and minus strand is shared. `near_block`, `swap`, `substitution` and
   `frameshift` tolerant rows do not make a peptide present; a `block` row does.
7. **Preconditions.** A non-wild-type run, a FASTA hash differing from `provenance.json`, a missing
   input, a genome without topology tags, or `--out` under the checkout each exits 2 and writes
   nothing. Cross-check 1 failing exits 1. Cross-checks 2 and 3 failing give `warnings` and exit 0.
8. **Determinism.** Two `compute` calls give byte-identical `entrapment.json` and
   `entrapment_discoveries.tsv`.
9. **Hook.** `metrics.py compute --entrapment` puts the blocks under `m8`, and `render` resolves
   every metric-8 placeholder. Without `--entrapment`, metric 8 renders as pending, as before.
10. **Guard.** `git ls-files` lists no `entrapment.json` and no `entrapment_discoveries.tsv`. This
    extends M10-05's guard test.
11. **Methods note.** `docs/methods/free-entrapment.md` exists, cites
    `10.1038/s41592-025-02719-x`, and contains equations (1) and (2) as written above. A text test
    checks the DOI and the two formula strings.

On FarmShare (evidence, not CI):

12. `entrapment.py compute` on the real inputs runs in `scripts/slurm/entrapment.sbatch` and exits
    0. The sbatch request is `dev`, 4 CPUs, 16 GB, 1 h, job name `giap-entrapment`; the expected
    run is minutes, dominated by the tolerant search of the K-12-only peptides of ≥ 10 aa. The job
    ID, r, `fdp_combined`, `fdp_lower`, the assessment and the warnings go to PROGRESS › Measured
    numbers as aggregates (under the M10-05 Q-K1 default, not committed elsewhere). A
    `digest_mismatch` warning is a lead finding before the number is reported.

## Test plan

The test-author writes the tests first (R2), in new files only:

- `tests/unit/test_m10_07_entrapment.py`: criteria 1–6;
- `tests/unit/test_m10_07_cli.py`: criteria 7–9 (reusing M10-05's synthetic run-directory
  builder);
- `tests/unit/test_m10_07_guard.py`: criteria 10–11.

The tests use small synthetic genomes and FASTAs built in `tmp_path`. Local runs are targeted
(venv Python, `PYTHONPATH=$PWD/src:$PWD/tests/sim:$PWD`). No integration test is needed, because
every piece of GIAP that `entrapment.py` calls is already integration-tested.

**Attack (statistician, optional; S4).** A simulation of the combined estimator under the free-
entrapment structure:

- **Setup.** Draw a synthetic K-12-like universe with share r/(1+r) entrapment, and draw false
  discoveries with equal chance, then with a homolog attraction factor λ > 1 towards entrapment.
- **Check.** The mean of `fdp_combined` is ≥ the true FDP at λ = 1 (the paper's averaged bound,
  within Monte-Carlo error), and the excess grows with λ.

It documents "errs high". If the statistician wants it, it is
`tests/calibration/attacks/test_attack_free_entrapment.py`, with seeds from `tolerances.yaml`.

## Estimate

| Piece | Code lines |
|---|---|
| digest universe, exclusions, non-standard filter | 30 |
| `present` (six-frame exact, strain-db substring, tolerant block-only) and the index build | 45 |
| precursor enumeration and m/z bounds, r and r_prec | 35 |
| discoveries, buckets, per run | 30 |
| estimator (Fraction), assessment, caveats | 20 |
| CLI, preconditions, settings parse, cross-checks, writers | 50 |
| `metrics.py` hook and template section (M10-05 code) | 15 |
| `scripts/slurm/entrapment.sbatch` | 25 |
| **Total** | **≈ 250** (index: 220; +30 for the precursor level and the sbatch script) |

Tests come to about 320 lines, plus the methods note (prose). One PR. It depends on M10-05's
`metrics.py` being merged, or stacked on its branch.

## Questions

### For the statistician (review before tests)

| # | Question | Default |
|---|---|---|
| S1 | Level. §15.5 defines N and r on peptides; DIA-NN controls precursor-level FDR at 1 %, and Wen 2025 assessed DIA-NN at precursor level. Which is primary? | peptide level primary (the spec's text); precursor level secondary with r_prec. The assessment label uses α = 0.01 at both levels, and the report notes that A1's peptide-level FDR is not controlled at exactly 1 % |
| S2 | "Mass-equivalent rules". Which exchanges make a K-12 peptide "present"? | the §7.2.4 `block` exchanges only, found by the tolerant route (≥ 10 aa, its budget). `near_block` (Δ ≤ 0.05 Da, e.g. K/Q-like) excluded. Alternative: also `near_block`, which would move peptides from entrapment to target and lower r and N_E |
| S3 | R4: an interval for `fdp_combined`? The source gives none; peptides of one protein are dependent | `ci: unestimated` with the reason; the lower bound is shown beside it. Alternative: a protein-block bootstrap over discoveries, with r fixed |
| S4 | The equal-chance assumption. Near-homologs bias high, but false matches could also be attracted to abundant shared proteins (bias low). Do you want the attack above, or a pre-registered secondary estimate restricted to entrapment peptides from K-12 proteins with no shared peptide ("foreign-like" stratum)? | neither by default; the bias table in the methods note states both directions in words |
| S5 | Distinctness: I/L-collapsed sequences for both U and discoveries (what DIA-NN can separate by mass), or literal sequences (what DIA-NN lists as distinct library entries)? | I/L-collapsed on both sides |
| S6 | The assessment labels (paper Fig. 1b readings) with inclusive `≤ α` for `consistent_with_control` and strict `> α` for `evidence_of_failure`. Is that acceptable? | yes |

### For the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| L1 | A constructed entrapment (shuffled `strain_db.faa` via FDRBench, one extra ~15-min DIA-NN search) as a cross-check? | **no**: §16 lists FDRBench entrapment as Later (R1). If wanted, an ADR and Khoa (§16) |
| L2 | "Nissle's translation": the six-frame genome translation only, or also the Bakta proteins (`strain_db.faa`, M at non-ATG starts)? With the six-frame only, every true N-terminal peptide of a non-ATG-start gene that K-12 shares would count as entrapment | six-frame ∪ `strain_db.faa` substrings (exact rule only) |
| L3 | Contaminant exclusion by sequence (D(contaminants)) rather than by `protein_ids` (M10-05's rule) | by sequence, because the in-silico side has no `protein_ids`. The difference shows up in cross-check 2 |
| L4 | H2 (§15.7). It names "entrapment FDP (full validation, FDRBench)". Should this free estimate stand in for it in the pilot? | no: H2 stays `deferred`; the report shows the genome-inconsistent share and `fdp_combined` side by side, without a verdict |
| L5 | Where metric 8 lives. A separate `entrapment.json` plus a `--entrapment` hook in M10-05's `metrics.py`, or code inside `metrics.py`? | separate script and file (§15.5 names `entrapment.py`), with a small hook in `metrics.py` (≈ 15 lines in this PR) |
| L6 | Run on Slurm (`dev`, 4 CPUs, 16 GB, 1 h) or in the session allocation? | Slurm: the tolerant search over ~10⁵ K-12-only peptides can exceed a few minutes on the shared 8-CPU session |
| L7 | Should the in-silico DIA-NN library (`A1.predicted.speclib`, on FarmShare) be used to check that the precursor space (charges 1–4, m/z bounds) matches DIA-NN's? The library is binary, and parsing it needs DIA-NN or a format reader | no library check in this issue. The peptide-level estimate does not depend on it (the m/z fact above); r_prec is labelled "from the configured bounds, not the library" |

## Statistician decisions (2026-10-24; binding)

Reviewer: statistician. Verdict: **approved, with the binding changes below**. Where this section
conflicts with the text above, this section wins; the "Supersedes" list names the lines it
replaces. The test-author writes tests against the body as amended here.

### What was verified

- **Source (R6).** Wen et al. 2025, *Nat. Methods*, doi:10.1038/s41592-025-02719-x, read from the
  publisher's article page (nature.com, open access) on 2026-10-24. Equation (1) is
  `N_E(1 + 1/r)/(N_T + N_E)`, with r "the effective ratio of the entrapment to original target
  database size", proved in Supplementary Note 2.2 to be "an estimated upper bound, that is, on
  average it overestimates the true FDP" under an equal-chance-like assumption. Equation (2) is
  `N_E/(N_T + N_E)`, a lower bound "assuming that any entrapment discovery is indeed false"
  (Suppl. Note 2.1). The sample method (3) is "inherently flawed". The three Fig. 1b readings are
  as the spec states. The paper evaluated DIA-NN at the **precursor** level ("each discovery is a
  (possibly modified) peptide and a corresponding charge state"). It also warns that a foreign
  entrapment may coincide with a source of contamination. The spec's quotations are accurate.
  The supplementary notes themselves were not read; nothing below depends on them beyond the
  statements quoted from the main text.
- **Algebra.** `1 + 1/r = (k + s)/k` with k = |k12_only| and s = |shared|. The rational form is
  therefore exact.
- **Examples 1–9 and the `MMGGGGGK` edge case** were recomputed with `giap.mapping.digest` and
  `giap.mapping.masses`. Every digest product, bucket, r, r_prec, precursor count, m/z value
  (719.41201, 360.20964, 735.40693, 368.20710, VEPTLDEK z1 930.47783) and FDP is correct.
  Examples 3–7 are also correct: 0.0201207 and 0.0040241; 0.0030272 and 0.0010091; 0.2 and 0.1;
  exactly 0.01; 0.0100033.
- **The m/z claim.** It holds. The mass range of a 7–30-aa peptide is [417.160785, 5600.389865]
  Da (all-G and all-W). Carbamidomethyl-C (160.03) and oxidised M (147.04) both lie inside the
  G–W residue range, so they cannot leave it. Charge 1 covers [298.99, 1798.99] Da and charge 4
  covers [1195.97, 7195.97] Da. These two intervals already overlap and contain the whole range,
  so every peptide in the searchable space has at least one in-bounds charge. The claim depends
  on `--min-pr-charge 1`, `--max-pr-charge 4` and the 300–1800 bounds, which is why criterion 4
  keeps its "bounds applied" test. It says nothing about any implicit pruning of charges inside
  DIA-NN (L7).
- **The q column.** M10-03 writes `q_value` from `Global.Q.Value`, with no fallback to the
  run-level `Q.Value`. That is an experiment-wide **precursor** q-value. The set of distinct
  precursors with `q_value ≤ 0.01`, pooled over the three runs, is therefore the unit whose FDR
  DIA-NN claims to control at 1 %. Peptide-level sets and per-run sets carry no such claim.
- **The search structure.** A1's "two step" is a library prediction followed by **one** search
  against the full predicted library (`diann_A1.cfg`). The competed space is the full in-silico
  digest, so the in-silico r is the right r. DIA-NN does not know which library entries are T and
  which are E, and its precursor q-values use no protein information. The search therefore
  treats T and E symmetrically, which is the precondition for the paper's argument.
- **Erratum (statistician, 2026-10-24, after the implementation review).** The previous bullet is
  wrong in two places, according to the DIA-NN README at 2fb4df29, the revision that documents
  2.3.1.
  - **Protein information.** Precursor q-values use protein information unless `--no-refine-q`
    is set, and A1 does not set it. T and E are therefore not treated symmetrically. The
    direction is low, and the methods note's bias table carries it. `consistent_with_control`
    has less support for A1 as a result. `fdp_lower` and `evidence_of_failure` are unaffected.
  - **Rows also pass the run-level filter.** The main report is already filtered at run-level
    `Q.Value` ≤ `--qvalue`, so accepted rows pass both the run-level and the global filter. The
    S1 per-run `assessment_reason` should therefore read "runs were selected by the global and
    the run-level q-value together; no FDR claim attaches to a run's intersection". That string
    is a code constant, and changing it is an implementer change. The precursor-level
    assessment stands.

### S1: level — **changed**

Exact text:

> The §15.5 metric-8 value is the **peptide-level** `fdp_combined`, with r at peptide level
> (primary number, unchanged). The **assessment** (the paper's three readings against α) is
> issued **only at the precursor level, pooled over runs**: distinct (`norm(modified_peptide)`,
> `charge`) of the accepted rows. This is the unit `Global.Q.Value` controls at α = 0.01. The
> peptide block carries `assessment: "not_assessed"` with `assessment_reason:
> "peptide-level FDR is not controlled by the search (q_value is DIA-NN's precursor-level
> Global.Q.Value)"`. Every `per_run` block carries `assessment: "not_assessed"` with
> `assessment_reason: "runs were selected by the global q-value, not a run-level q-value"`.
> Per-run numbers are descriptive.

Why: Fig. 1b compares an FDP bound with the FDR the tool claims **for that list**. DIA-NN claims
nothing for the distinct peptides derived from its precursor list. False precursors spread over
distinct peptides, while true ones cluster several precursors per peptide, so a peptide-level FDP
above 1 % is not evidence that DIA-NN failed. The spec's default would print
`evidence_of_failure` for exactly that case (example 2), which R4 forbids.

Example 2 changes accordingly. Peptide block: `assessment: not_assessed`, and the numbers are
unchanged (0.75, 0.25). Precursor block: `evidence_of_failure` (fdp_lower 0.2 > 0.01). `per_run`:
`not_assessed` everywhere, and the numbers are unchanged.

The pure `assess(combined, lower, alpha)` function keeps the three labels, and examples 3–7 test
it directly.

### S2: mass-equivalent rules — **default approved, plus one sensitivity**

> "Present" by rule 2 requires a tolerant target row whose `diff_types` is non-empty and ⊆
> {`block`}: the §7.2.4 "mass-equivalent block exchange" list, Δ < `BLOCK_DELTA_DA`. `swap`,
> `near_block`, `substitution` and `frameshift` rows do not make a peptide present. This follows
> §7.2.4's own vocabulary, where only the block list is called mass-equivalent.

Required addition, a pre-registered descriptive sensitivity `sensitivity.present_broad` in the
peptide block. Rule 2 is widened to rows with `diff_types` ⊆ {`block`, `swap`, `near_block`}:
adjacent swaps are isobaric, and near-blocks are within 0.05 Da (§7.2.4, §8.5 attack 6). The block
reports `shared`, `k12_only`, `r`, `N_T`, `N_E`, `fdp_combined`, `fdp_lower` and `ucb95` (S3),
with no assessment. It reuses the tolerant rows already computed, so no new search is needed.

A known limit goes in the methods note's bias table, direction high and negligible: the tolerant
route reports only the minimal-cost loci, and breaks ties within a locus in favour of
`substitution`. A block-only locus that sits behind a cheaper non-block locus is therefore not
seen.

New example **1c** (swap). This is example 1, but the translation holds `AGNVLTEPRW` instead of
`AGGGVLTEPWR`. Primary: `AGNVLTEPWR` has a `swap` row (10 aa, budget 1), so it is `k12_only`;
shared 3, k12_only 3, **r = 1.0**. `present_broad`: it becomes shared; shared 4, k12_only 2,
**r = 0.5**.

### S3: interval — **changed (an upper confidence bound is added)**

> `fdp_lower` gets no interval. It is not an estimate: if every entrapment discovery is false, the
> realised FDP of this list is ≥ `fdp_lower`, with no sampling error. The condition is that the
> genome is complete and identical to the sample's strain, and that no K-12 material is in the
> sample. The report calls it "realised lower bound (conditional on …)".
>
> `fdp_combined` gets a one-sided 95 % **upper confidence bound** `ucb95`. The model is the
> paper's equal-chance assumption in its plainest form, with false matches independent: given M
> false discoveries among the N = N_T + N_E, N_E ~ Binomial(M, π) with π = k/(k + s) = r/(1 + r).
> Then `fdp_combined` = N_E/(π N), and E[N_E/π] = M. The bound inverts the one-sided test:
>
> `M_max = max{ M ∈ {N_E, …, N} : BinomCDF(N_E; M, π) > 0.05 }`, `ucb95 = M_max / N`
>
> CDF(N_E; N_E, π) = 1, so M_max ≥ N_E. Use `scipy.stats.binom.cdf` and a monotone scan or
> bisection over M. `ucb95` is `null` whenever `fdp_combined` is `null`, with the same reason.
> At the precursor level, π = K/(K + S) from the precursor sums.
>
> The JSON field `ci` becomes `{"ucb95": <float|null>, "method":
> "binomial_inversion_equal_chance", "assumes": "equal chance and independent false matches;
> not valid under clustered (protein-level) errors"}`. No two-sided interval is given.

Why: a bare `0.0` (example 7) or a per-run 3.0 says nothing about precision. The binomial bound
costs about 10 lines and has nominal coverage under the stated model (the attack checks this). A
protein-block bootstrap was rejected because the discoveries are a single realised list; the
resampling unit for false matches is not the protein of the discovered peptides.

`ucb95` values, added to the examples (tolerance 1e-12, since they are rationals M_max/N):

| Example | π | N_E, N | M_max | ucb95 |
|---|---|---|---|---|
| 2, peptide | 1/3 | 1, 4 | 4 | 1.0 |
| 2, precursor | 1/4 | 1, 5 | 5 | 1.0 |
| 2, per run WT_3, peptide | 1/3 | 1, 1 | 1 | 1.0 |
| 3 | 1/5 | 20, 4,970 | 140 | 140/4,970 = 0.0281690 |
| 4 | 1/3 | 10, 9,910 | 47 | 47/9,910 = 0.0047427 |
| 5 | 1/2 | 100, 1,000 | 225 | 0.225 |
| 6 | 1/3 | 10, 3,000 | 47 | 47/3,000 = 0.0156667 (N_T 2,989: 47/2,999 = 0.0156719) |
| 7 | 2/7 | 0, 5,000 | 8 | 8/5,000 = 0.0016 |

The labels stay on the point estimates, as in the paper (S6). `ucb95` is reported beside them.
The report states it in words: "with 95 % confidence, under equal chance, the FDP is at most …".

### S4: equal-chance bias — **attack required; a stratum replaces the "foreign-like" proposal**

Does free entrapment work? Yes, with these caveats in the bias table:

- Every K-12-only peptide is false only if the genome is complete and is the sample's strain, and
  no K-12 material (carryover, lab contamination) is present. K-12 is a common lab strain, and the
  paper names contamination as the foreign-entrapment risk. If any of these fails, the direction
  is high, and the lower bound stops being a bound.
- Neighbour-driven errors concentrate on E, because every K-12-only peptide near a Nissle variant
  has its neighbour's signal in the sample. Direction: high (§15.5).
- Six-frame chance hits move short K-12-only peptides into `shared`. About 0.08 expected exact
  chance hits per 7-mer and 0.005 per 8-mer (§7.2.5 table). Direction: slightly lower N_E, so a
  more conservative lower bound; negligible for the combined estimate.
- E and T may differ in length and charge composition (variant-bearing peptides are longer).
  Under equal chance **per precursor**, a peptide's chance of a false discovery scales with
  n_prec, so r at peptide level may be mis-weighted. Direction: either.
- A false-match attraction towards T (for example interference from abundant co-eluting shared
  proteins) biases low. Unquantified.
- The report has no protein-level filter, and must not get one. A `PG.Q.Value` filter would
  remove K-12-only proteins first and bias low. The methods note states, from the DIA-NN 2.3.1
  documentation, which filters the main report applies. A protein-level filter is a lead finding.

Required additions:

1. **`sensitivity.r_prec_weighted`** (peptide block): `fdp_combined_rprec = N_E (1 + 1/r_prec) /
   (N_T + N_E)`, plus `mean_nprec_shared` and `mean_nprec_k12_only`. Example 2: 1 · 4/4 =
   **1.0**; means 18/4 = 4.5 and 6/2 = 3.0.
2. **`sensitivity.far`** (peptide block), pre-registered and descriptive, with no assessment.
   E_far = `k12_only` peptides of ≥ 10 aa with **no tolerant target row of any type** (no Nissle
   near-homolog within the §7.2.4 budget). T_10 = `shared` peptides of ≥ 10 aa.
   r_far = |E_far| / |T_10|. Discoveries are restricted the same way:
   `fdp_combined_far = N_E_far (1 + 1/r_far) / (N_T_10 + N_E_far)`, plus `fdp_lower_far` and
   `ucb95_far`. It also reports the breakdown of `k12_only` and of N_E into `lt10` (not searched),
   `near_homolog` (≥ 1 tolerant target row) and `far`. This replaces the "K-12 proteins with no
   shared peptide" stratum: a protein with no shared peptide can still be a close homolog, whereas
   "no tolerant row" directly measures the neighbour that drives the high bias, and costs nothing
   (the rows exist). Example 1 / 2: `k12_only` lt10 1 (`AAGGWWR`), near_homolog 1
   (`VEPTLDEKAAGGWWR`, `substitution`), far 0. N_E lt10 1, near_homolog 0, far 0. |E_far| = 0, so
   `fdp_combined_far` is `null` (`no_entrapment_space`), `fdp_lower_far` = 0.0, N_T_10 = 1
   (`AGNVLTEPWR`).
3. **Attack** `tests/calibration/attacks/test_attack_free_entrapment.py`, marker `calibration`
   (fast). It is written by the statistician against the pure API below, after the test-author's
   files are committed. Under §8.5 it gates the merge of the implementation PR. Seed
   `seeds.calibration_base` + 707, k = `se_multiplier` (no `tolerances.yaml` change), R = 2,000
   replicates (a module constant; pure numpy, seconds). The simulated universe keeps the real
   structure: peptides grouped into proteins (sizes 1–30, geometric), s = 4,000 shared peptides,
   k = 1,000 K-12-only peptides (r = 0.25), and a true-present subset of 2,000 shared peptides.
   Each replicate draws true discoveries (each present peptide with probability 0.9) and
   M ~ Poisson(30) false matches. The scenarios:
   - **A, edge equal chance:** each false match lands in E with probability π. Check that
     mean(`fdp_combined`) ≥ mean(true FDP) − k·SE.
   - **B, uniform over null elements** (E ∪ absent shared): the same check.
   - **C, homolog attraction λ ∈ {1, 2, 4}** (E weight × λ): the mean excess
     `fdp_combined − FDP` is increasing in λ.
   - **D, attraction to T** (λ = 0.5): mean(`fdp_combined`) < mean(FDP). This confirms the
     low-direction row of the bias table.
   - **E, clustered:** false matches arrive in whole proteins, each cluster in E with
     probability π. The mean check of A passes; `ucb95` coverage is not gated, which is why its
     `assumes` string says so.

   Further checks, in every scenario and replicate:
   - `fdp_lower` ≤ true FDP (deterministic, since E is false by construction);
   - no replicate with true FDP ≤ α is labelled `evidence_of_failure`;
   - in A and B, the coverage P(FDP ≤ `ucb95`) ≥ 0.95 − k·√(0.95·0.05/R).

   As written (2026-10-24): the fast tier uses R = 2,000. A `calibration_full` tier uses
   R = 20,000 (Slurm, about 45 s). Both run under `calsuite("compare")`.

   There are two teeth checks, both in A:
   - `fdp_lower` used as an estimate fails the mean check;
   - `ucb95(level=0.5)` fails the coverage check.

   The coverage of C, D and E is recorded but not gated. Against a scratch reference, E covered
   about 0.72 and D about 0.75.

**Pure API** (in `validation/pilot/entrapment.py`, used by the attack and the unit tests):
`combined(n_t, n_e, n_k12_only, n_shared) -> Fraction | None`, `lower(n_t, n_e) -> Fraction |
None`, `ucb95(n_t, n_e, n_k12_only, n_shared, level=0.05) -> float | None`, `assess(combined,
lower, alpha=Fraction(1, 100)) -> str`. The `None` cases follow "Edge cases".

### S5: distinctness — **default approved**

> I/L-collapsed on both sides and at both levels. At the precursor level the unit is
> (`norm(modified_peptide)`, `charge`), and `n_prec` is counted on collapsed universe elements.

Two K-12 library entries that differ only by I/L would give a collapsed element two chances of a
false match. This is rare, and the methods note mentions it.

### S6: labels and boundaries — **approved**

> `consistent_with_control` iff `combined ≤ α`; `evidence_of_failure` iff `lower > α` (checked
> first; the two cannot both hold, since lower ≤ combined); else `inconclusive`. Both comparisons
> are made in `fractions.Fraction`, with α = `Fraction(1, 100)` (never the float 0.01). The labels
> are issued only where S1 allows.

### Testing families (§8.1) and calibration

Metric 8 is a descriptive estimate. It makes no p or q, and it is not in any §8.1 family: the
precursor row "audited with entrapment in validation" is exactly this audit. The assessment label
is a reading, not a test, so no level criterion applies. The calibration evidence is the attack
above: the bias direction, the deterministic lower bound, and `ucb95` coverage.

### Other required changes

- **Mapping wording.** Replace "is exactly the FDP of A1's accepted set" with "is the FDP of
  A1's accepted set after removing contaminant-sequence and `not_in_universe` discoveries".
- **Data rule (Khoa, 2026-10-24).** Aggregates may be committed; row-level DIA-NN content may
  not. `entrapment_discoveries.tsv` is row-level and is never committed. `entrapment.json` holds
  aggregates only. Committing it is the lead's decision at report time, not this PR's. Criterion
  10 becomes:
  - `git ls-files` lists no `entrapment_discoveries.tsv`;
  - a content guard: in the CLI test, no discovered peptide string, and no `modified_peptide`,
    occurs anywhere in the `entrapment.json` text.
- **Caveat strings.** The third caveat becomes "assumes the Nissle genome is complete, is the
  sample's strain, and that no K-12 material is in the sample". A fourth is added: "assessment is
  made at precursor level, the unit DIA-NN's global q controls".
- **Methods note.** It adds the S3 bound with its assumption and the S4 bias rows above. It
  states that `fdp_lower` is a realised, conditional bound, and why the peptide level carries no
  assessment.
- **Acceptance criteria**, additions in the same style:
  - **S-1:** examples 1c, the `ucb95` table, and the S4 sensitivity values hold, with
    `assessment` per S1.
  - **S-2:** `ucb95` equals a brute-force scan over M for 200 seeded (N_T, N_E, k, s).
  - **S-3 (reworded 2026-10-24, statistician; the earlier "`ucb95` ≥ `fdp_combined`" was false).**
    Write N = N_T + N_E, π = k/(k + s), and M_max = N · `ucb95`. Note that N · `fdp_combined` =
    N_E/π = N_E (k + s)/k. Whenever both are defined:

    > M_max ≥ min(⌊N_E (k + s)/k⌋, N), that is,
    > `ucb95` ≥ min(⌊N · fdp_combined⌋ / N, 1) > min(`fdp_combined`, 1) − 1/N.

    Proof: let M = min(⌊N_E/π⌋, N). Then M ≥ N_E and Mπ ≤ N_E. Every median of
    Binomial(M, π) is at most ⌈Mπ⌉ (Kaas and Buhrman 1980), and ⌈Mπ⌉ ≤ N_E. So
    BinomCDF(N_E; M, π) ≥ 1/2 > `level` for any `level` < 1/2, and M lies in the scan set.
    The bound is therefore never below the point estimate by a full grid step 1/N, and equals 1
    whenever `fdp_combined` ≥ 1.

    The weaker form `ucb95` ≥ `fdp_combined` is not a property of this bound, and that is **not
    a defect**. `ucb95` bounds the realised FDP M/N, which lies on the grid {0, 1/N, …, 1}, and it
    has exact one-sided coverage under the model (it inverts an exact test). `fdp_combined`
    estimates E[M]/N off that grid and is unclipped. They part in two places:
    - **Above one.** If `fdp_combined` > 1, then `ucb95` = 1 (the spec's per-run WT_3 row: 3.0
      against 1.0).
    - **π near 1.** Binomial(M, π) is then almost degenerate, so the bound stops at ⌊N_E/π⌋/N.
      An example: N_T 1, N_E 1, k 39, s 1 gives 0.5 against 0.513. On the grid N ≤ 1,000,
      N_E ≤ 40, k, s ≤ 40, `ucb95` < min(`fdp_combined`, 1) occurs only at π ≥ 0.975. The pilot
      (r of order 0.1–1, π ≤ 0.5) is far from that.

    Rounding `ucb95` up to `fdp_combined` would add nothing to its coverage and would break
    S-2's definition, so it is not done.

    Tested as follows:
    - the inequality above, for every π (exact, in integers);
    - `ucb95` ≥ min(`fdp_combined`, 1) for π ≤ 0.9, as a regression property (empirical; the
      grid above, no proof);
    - both corners, pinned.

    With N_E = 0, `ucb95` ≥ 0 = `fdp_combined` trivially.
- **Estimate:** about +45 code lines (ucb95 10, the three sensitivities 30, assessment scoping
  5) and about +60 test lines. The attack (about 150 lines) is the statistician's and does not
  count against the PR.

### Supersedes

- Line "R4: … its sampling uncertainty is `unestimated` unless the statistician says otherwise"
  → S3.
- Estimates: the `assessment` bullet → S1 and S6 (scope); the `ci` bullet → S3; the `caveats`
  bullet → "Other required changes".
- Levels paragraph → S1.
- Example 2's peptide-level and per-run labels → S1.
- Criterion 10 → "Other required changes".
- "Attack (statistician, optional; S4)" → required, as specified in S4.

## Lead decisions after the code review (rev-1007, 2026-10-24; binding)
These supersede earlier text that differs.
- **B1:** cross-check 1 considers only discoveries in the T and E buckets. Unclassified discoveries are counted and reported, never an exit-1 failure.
- **B2:** the sbatch default output path must be git-ignored (`git check-ignore`). The default run must not exit 2.
- **B3:** merge main, including M10-06's `expand()`/INPUTS hooks, and re-apply the M10-07 hook on top.
- **S1:** an unexpected exception exits 4, with a traceback; the docstring lists the exit code. Exit 1 stays reserved for cross-check 1.
- **S2:** provenance matches each input by basename role, not by set membership.
- **S3:** cross-check 1 requires `exact`; `strain_db` is not agreement.
- **S4:** a missing `detection_summary.tsv` gets its own warning code, not cross-check 2's.
- **S5:** `OUT_A1` defaults to the M10-04 A1 directory. If it is set explicitly and missing, the run exits 2. The sbatch never skips the cross-checks silently.
- **S6:** done in the methods note (lead, 3953908d).
- **S7:** pin the mutation survivors listed in the review.
- **S8:** reuse the M10-05 metrics helpers and the M10-03 matching helpers instead of copies. Target ≤ about 400 counted lines in entrapment.py.
- **Statistician (2026-10-24):**
  - S1: m8 carries the warnings and the assumption string;
  - the erratum: precursor q-values use protein information, so the per-run reason string changes (see "Statistician decisions").
