# M7-18: The variants step on ambiguous and rare residues (B, J, O, U, X, Z)

Milestone M7 · track V (variants) · **pass 2 (hardening; a crash fix found by the M7-15 K-12
real-tool run)** · depends on: merged M7-01 … M7-05, M7-14 · touches: `variants/pileup.py`,
`steps/variants.py` · estimate ≈ 75 counted code lines (+ ≈ 220 test lines) · not statistical
(counting rules; no threshold, null or FDR changes)

## The bug

The M7-15 K-12 run (TE8 peptides, `farmshare,apptainer`) fails in the `variants` step with

```
ValueError: no table-11 codon encodes 'X'
```

raised by `minimal_alt_codon` in `src/giap/variants/vcf_records.py`. `flags.nt_class` calls the
same function and fails the same way. 13 of the 151,424 K-12 TE8 peptides contain B, J, O, U, X
or Z.

The cause is a chain of four steps. Each one is correct on its own.

1. **Preflight.** The peptide table accepts any uppercase letters A–Z (§5, `peptides.tsv`). For
   proteins, preflight counts X, B, Z and J and only warns (§5; §7.0 table: "20 aa + U; X, B, Z, J
   counted … Warn"). Nothing removes these letters.
2. **Mapping.** `encode_query` maps every letter outside the 19 residue codes (B J O U X Z) to
   `NOMATCH` (`mapping/sixframe.py`). The exact route can never match such a letter: it decodes to
   `#`. The M3-06 DP can consume it, but only by a **substitution** (cost 1), and never inside a
   swap, block or near-block (the DP guards `q0, q1 < nres`). So a peptide of 10 aa or more with at
   most k such letters is placed by the tolerant route. Each of those letters is one
   `substitution` difference.
3. **Alignment (M7-01).** The traceback records that residue as
   `ResidueOp("substitution", ref=<genome residue>, alt="X")`.
4. **Pileup (M7-02 rule 4).** Every substitution becomes ALT `(alt, substitution)`. So an X
   becomes ALT `X`, and `minimal_alt_codon("…", "X")` raises when the VCF record or the
   `snv`/`multi_nt` flag is built.

I checked the following against main at 94bb0029, with a 12-codon `+` ORF `AKCREFGHILWD` at
[0, 36):

| Peptide | Op recomputed by `align_hit` | `_category` today | Outcome today |
|---|---|---|---|
| `AKXREFGHILWD` | substitution C→X at codon 6 | ALT `X` | **crash** |
| `XKCREFGHILWD` | substitution A→X at codon 0 | ALT `X` | **crash** |
| `AKCREFGHILWX` | substitution D→X at codon 33 | ALT `X` | **crash** |
| `AKUREFGHILWD` | substitution C→U at codon 6 | ALT `U` | **crash** |
| `AKCREFGHJLWD` | substitution I→J at codon 24 | ALT `J` | **silent false ALT**: `minimal_alt_codon("ATT", "J")` returns `ATT` with 0 changes, which gives a VCF record with ALT = REF |
| `KAXREFGHILWD` (n_diff 2) | swap AK→KA at codons (0, 3), then substitution C→X at codon 6 | ALT `AK>KA` and ALT `X` | **crash** |

The J row is a second bug: M7-02 rule 5 uses `J` as the *internal* token for "the peptide read I
or L". A literal `J` from the input arrives with the same letter. Today it becomes an ALT even over
an I or L codon.

U over a TGA codon (selenocysteine) and O over a TAG codon (pyrrolysine) **cannot reach the
pileup** with the merged matcher. No DP move consumes a target `*`, and `align._residues` raises on
a `*` in the locus (M7-16 "The gap"). The rule below covers them anyway (lead decision 9). They
become reachable if M7-16 option B is accepted, and unit tests can build such alignments by hand.

## Spec sections implemented

- §7.5 preamble: "whether the peptides support the genome's residue (REF), a different residue
  (ALT), or both."
- §7.5.2 steps 2–4: the pileup, calls and confidence, unchanged for standard residues.
- §7.5.5: `pileup.tsv.gz`, `variants.tsv`, `variants.vcf`. Their columns, schemas and headers are
  unchanged.
- §7.2.1 ("ambiguous codons → `X`, which never matches"), §7.2.4 and M3-06 ("A query `NOMATCH`
  residue can only be consumed by a substitution"). These are the mapping convention that this
  issue makes the variants step follow.
- §5 / §7.0: the protein alphabet (20 aa + U; X, B, Z, J counted; warn), and `peptides.tsv`
  (uppercase A–Z).
- M7-15 lead decision 9 (docs/specs/M7-15-integration-hardening.md, in the m7-integ worktree):
  - an ambiguous residue (B, J, X, Z) is never an ALT residue;
  - U over a TGA codon and O over a TAG codon count as matches (selenocysteine and pyrrolysine
    readthrough), not ALTs;
  - any other U or O is never an ALT;
  - each skipped site is counted in one event.

## The mapping convention, and why the variants step skips rather than matches

The mapper gives **no** partial credit to any of B, J, O, U, X or Z. All six are `NOMATCH`: an X
never matches anything, a B never matches D or N, a Z never matches E or Q, and **a J never
matches I or L**. Each one costs a full substitution in the placement cost and the placement FDR.
I/L folding (I → L in the index and the queries) applies only to the letters I and L.

The variants step must not invent a second convention. If B over D counted as REF here, a peptide
would be scored as a mismatch when it is placed and as a match when it is piled up. So, at a codon
where the peptide's residue is one of these six letters (except the two readthrough cases), the
peptide gives **no residue evidence**. That position is neither REF nor IL nor ALT, and it is
**skipped**. This is the conservative reading of "ambiguous". It is also what lead decision 9 asks
for ("never an ALT"). Whether B over D/N should count as REF is a mapping question (Q2): changing it
would change Module 2 first, and this step would follow.

The two readthrough cases are the only exception. Here the genome codon is a stop, which the
mapper never consumes today. The residue is the standard decoding of that stop codon (UGA → Sec,
UAG → Pyl), so lead decision 9 counts it as REF.

## Public interface

No new public function, no new column and no schema change.

`src/giap/variants/pileup.py` (pure):

```python
NONSTANDARD: frozenset[str] = frozenset("BJOUXZ")     # peptide letters that are never ALT
READTHROUGH: Mapping[str, str] = {"U": "TGA", "O": "TAG"}   # peptide residue -> gene-strand codon
NONSTANDARD_EVENT = "variants.nonstandard_residues"   # defined in steps/variants.py
```

- `pileup(...)`: same signature. The returned frame gains `attrs["nonstandard"]`, a plain dict
  (rule 6). The internal `_alts` column and the other attrs are unchanged.
- `eligible(...)`: same signature. Rule 3 adds one keep-condition.
- `steps/variants.py::variants_step` emits `NONSTANDARD_EVENT` (rule 7). `run_variants` and
  `run_variants_qc` keep their signatures and return values.
- `vcf_records.minimal_alt_codon` and `flags.nt_class`: unchanged. `minimal_alt_codon` still
  raises on an unencodable residue, which is now a contract violation (rule 5).

## Rules

The "peptide residue" of an op is `op.alt` as written in the input. I and L are kept distinct, and
the letters are taken **before** the internal I/L → `J` token of M7-02 rule 5 is applied.

1. **Mapping (Module 2): no change.** B, J, O, U, X and Z stay `NOMATCH`. They are consumed only by
   a substitution, never by a swap, block or near-block. A target `*` or `X` is never consumed. The
   placement, its cost and its q-value are what M3 computed. Consequence: a placed peptide
   contains at most k = `mapping.tolerant` budget such letters (1 at 10–11 aa, 2 at ≥ 12 aa), and
   none at all if it was placed by the exact route.
2. **Alignment (M7-01): no change.** `align_hit` keeps recording the letter as a `substitution` op
   with `alt` = the raw letter. `_residues` keeps raising on a genome `*` or `X` in the locus (a
   readthrough alignment needs M7-16).
3. **Eligibility (M7-02 rule 1): one added condition.** Under `genome_relation = same`, a
   `database` peptide with `search_db = other` keeps its `match`/`il` ops, as today. It also keeps
   every **readthrough-shaped** op: a `substitution` with `ref == "*"` and `alt in {"U", "O"}`. The
   pileup then decides REF or skip with the codon bases (rule 4b). Every other op is stripped
   exactly as today. That includes substitutions to B/J/X/Z, which therefore never reach the
   pileup and are **not** counted in the event (§7.5.1 had already excluded them).
4. **Pileup (M7-02 rule 4): classify each op at its first codon `codons[0]`, in this order.**
   - a. If `op.kind == "substitution"` and `op.alt` is a single letter in `NONSTANDARD`, apply b
     to d. Otherwise the op is classified exactly as today: `match` → REF, `il` or I↔L → IL,
     substitution → ALT (I/L → `J`), swap/block → ALT token. This path does not change by a
     single byte.
   - b. **Readthrough match.** If `op.alt` is `U` and the gene-strand codon bases
     (`_codon_nt(seq, rep, codon, strand)`, so minus strands are reverse-complemented and an
     origin-crossing codon is read across the origin) are `TGA`, **or** `op.alt` is `O` and the
     bases are `TAG`, then the peptide supports **REF** at that codon. It counts in `n_ref`
     exactly like a `match`.
   - c. **Skipped (rare).** Any other `U` or `O` (over a sense codon, over TAA, U over TAG, O over
     TGA) is skipped.
   - d. **Skipped (ambiguous).** `B`, `J`, `X` and `Z` are always skipped, over any codon. That
     includes J over I or L, B over D or N, and Z over E or Q.

   A **skipped** op contributes no category: no REF, no IL, no ALT, no `cats` entry and no
   `peps` entry for that `(codon, folded peptide)`. Its codon still counts in the
   frameshift-spanning `cover` set, because the peptide spans that codon contiguously in this frame
   (M7-02 rule 8; Q4). A codon whose only ops are skipped has **no pileup row**: it is not
   "covered by" any residue evidence. A codon where some peptides are skipped and others are not
   gets the usual row from the others, with the skipped peptides absent from `n_ref`, `n_il`,
   `alt_counts`, AD and DP.

   Defensive: the DP cannot produce a swap, block or near-block op whose peptide side contains a
   `NONSTANDARD` letter. If one is met, the op is skipped too (rule 4d) and counted under its
   `op.alt` string. It does not raise.
5. **ALT construction, call and confidence: no change.** Because rules 4b–d run first, every ALT
   residue that reaches `sites`, `flags_for`, `minimal_alt_codon` and the VCF is a standard
   residue, the token `J` (folded I/L) or a block/swap token over standard residues. The call
   (`REF`/`ALT`/`mixed`/`no_call`) and the confidence (`high` ≥ `variants.high_min_peptides`)
   follow M7-02 rules 6–7 on the remaining counts. A readthrough peptide makes a codon `REF` (or
   `mixed` with an ALT peptide). It never gets the `stop_codon` flag, which needs an ALT site (Q1).
6. **Counting (`attrs["nonstandard"]`).** A **skipped site** is one `(codon key, folded peptide)`
   pair whose op at that codon was skipped by rule 4c or 4d. Duplicates through two loci of one
   peptide (unwrapped and wrapped) collapse by the M7-02 rule 3 key. A **readthrough site** is one
   pair counted REF by rule 4b. The dict holds JSON-native values, with keys sorted:

   | Key | Value |
   |---|---|
   | `skipped_sites` | int: number of skipped `(codon, folded peptide)` pairs |
   | `skipped_codons` | int: distinct codon keys with ≥ 1 skipped pair |
   | `skipped_peptides` | int: distinct folded peptides with ≥ 1 skipped pair |
   | `by_residue` | dict `{letter: n_pairs}` over skipped pairs, keys sorted (for example `{"J": 1, "X": 3}`) |
   | `by_rule` | dict `{"ambiguous": n, "rare": n}`: rule 4d versus 4c pairs |
   | `readthrough_sites` | int: rule 4b pairs (`U` at TGA plus `O` at TAG) |
   | `examples` | list of ≤ 10 original peptide strings with a skipped or readthrough pair, sorted lexicographically |

   With nothing skipped and no readthrough, it is the same dict with all zeros, empty `by_residue`,
   zero `by_rule` entries and an empty `examples` list.
7. **The event.** After `run_variants_qc`, `variants_step` emits **exactly one** event
   `variants.nonstandard_residues` at INFO (Q3), with `data = pile.attrs["nonstandard"]`. It does
   so only when `skipped_sites + readthrough_sites > 0`. The message is
   `"variants: %d site(s) with an ambiguous or rare peptide residue skipped (%s), %d readthrough
   match(es)"`, with `by_residue`. With no such sites there is no event, so the event log of an
   all-standard run is unchanged. `variants.summary` is unchanged.
8. **VCF writing and the variants table.** No code change. A skipped peptide is not in any AD or
   DP (`_dp` = n_ref + n_il + Σ n_alt over the pileup row). A codon with only skipped peptides
   has no pileup row, no `variants.tsv` row and no record. A readthrough-only codon is a `REF`
   pileup row with `ref_aa = "*"`. It has no `variants.tsv` row and no VCF record.
9. **Standard residues: byte-identical.** For inputs whose eligible peptides contain none of
   B, J, O, U, X or Z, `pileup.tsv.gz`, `variants.tsv`, `variants.vcf`, `variants_qc` and the
   event log must be byte-identical to main. That covers every committed golden and the M7-15
   goldens once they merge.
10. **Distinct-peptide folding: no change.** `_fold` folds I → L only. A literal `J` is not folded
    to `L` (mapping does not treat J as I/L either; Q5).

## Examples

Linear contig `c1`, 66 nt. Gene G1 is `+` [0, 36) and translates to `AKCREFGHILWD`. Its codons are
A@0 K@3 C@6 R@9 E@12 F@15 G@18 H@21 I@24 L@27 W@30 D@33. All peptides are eligible
(q ≤ 0.01, unique) and come from de novo input with `genome_relation = same`, unless stated. P0
is `AKCREFGHILWD` itself (exact, REF at every codon).

| # | Peptides | Expected |
|---|---|---|
| 1 | P1 `AKXREFGHILWD` (n_diff 1) | **Regression.** No crash. Codon 6 has no row (only P1 covers it, and P1 is skipped there). The other 11 codons are `REF` with `n_ref` 1. `variants.tsv` and VCF are header-only. Event data: `skipped_sites` 1, `skipped_codons` 1, `skipped_peptides` 1, `by_residue` `{"X": 1}`, `by_rule` `{"ambiguous": 1, "rare": 0}`, `readthrough_sites` 0, `examples` `["AKXREFGHILWD"]` |
| 2 | P0 and P1 `AKXREFGHILWD` | codon 6: row `ref_aa` C, `n_ref` 1 (P0), `n_il` 0, `alt_counts` NA, call `REF`. DP would be 1 (P1 is not counted). Event `skipped_sites` 1 |
| 3 | P1 `AKXREFGHILWD`, P2 `AKDREFGHILWD` (C→D) | codon 6: `alt_counts` `D:1`, call `ALT`, `candidate`, `n_alt` 1. One VCF record: AD `0,1`, DP 1 (not 2). Same output as with P2 alone, except for the event |
| 4 | P1 `AKCREFGHJLWD` (J at the I codon 24) | codon 24 has no row (today: a false ALT `J` with ALT codon = REF codon). `by_residue` `{"J": 1}` |
| 5 | P1 `AKCREFGHILWD` with V at codon 24 in the genome (`GTT`); peptides P2 `…J…` and P3 `…I…` at codon 24 | P3 gives ALT token `J` (I/L fold, M7-02 rule 5), `n_alt` 1, `candidate`. P2 is skipped. So J-the-letter and J-the-token are separate: only P3 supports the ALT |
| 6 | P1 `AKUREFGHILWD` (U over the C codon `TGT`) | codon 6 has no row. `by_rule` `{"ambiguous": 0, "rare": 1}`, `by_residue` `{"U": 1}` |
| 7 | Hand-built `HitAlignment` (unit test only; unreachable until M7-16): codon 6 rewritten to `TGA`; op `substitution` ref `*` alt `U` at codon 6 | codon 6: row `ref_aa` `*`, `n_ref` 1, call `REF`. No site and no record. `readthrough_sites` 1, `skipped_sites` 0. The event is emitted (readthrough > 0) |
| 8 | Same as 7, with alt `O` over `TGA`; and alt `U` over `TAG` | both skipped (rule 4c), `by_rule.rare` 1 each, no row |
| 9 | P1 `KAXREFGHILWD` (n_diff 2: swap + X) | codon 0: ALT token `AK>KA`, `alt_type` swap, `candidate`. Codon 3 gets no support from P1 (M7-02 decision 4). Codon 6 skipped. One VCF record at POS 1 (6-nt REF `GCTAAA`, ALT `AAAGCA`: minimal K codon from `GCT` is `AAA`, ties to lexicographic, then minimal A codon from `AAA` is `GCA`). `skipped_sites` 1 |
| 10 | Minus strand: G2 `−` [100, 136) on `c1′` (circular, 200 nt), translating to the same `AKCREFGHILWD`; P1 `AKXREFGHILWD` | the X sits at the 3rd residue, whose codon key is 136 − 9 = **127** (lowest coordinate). No row at 127; the rest REF. The same counts as example 1 |

## Edge cases

| Situation | Behaviour |
|---|---|
| X at the peptide's **first** residue (`XKCREFGHILWD`) | codon 0 skipped (no row unless another peptide covers it). Codons 3 … 33 `REF`. Event counts 1 |
| X at the **last** residue (`AKCREFGHILWX`) | codon 33 skipped. Same as above |
| a peptide that is **all X** (or has more nonstandard letters than k) | never placed by Module 2 (cost = length > k), so never in `peptide_hits.tsv`. A hand-built alignment whose every op is skipped gives no pileup rows and `skipped_sites` = its length. Header-only outputs if it is the only peptide |
| `XX` + 10 standard (12 aa, k = 2) | both codons skipped. `skipped_sites` 2, `skipped_peptides` 1, `skipped_codons` 2 |
| adjacent swap next to X | the swap is a normal ALT. X is skipped (example 9). X can never be **inside** the swap (DP guard) |
| X in a **de novo** read vs a **database** peptide | the same rule for both. With `genome_relation = same` and a database peptide searched with `search_db = other`, the X op is stripped by `eligible` before the pileup and is **not** counted (rule 3). The readthrough-shaped op survives `eligible` |
| `genome_relation = related` | the same rule. A skipped op never gets `strain_difference` because it is never a site |
| frameshift peptide with X in the prefix or suffix | the X op is skipped. The frameshift event (`fs±1`) is unchanged. The skipped codon still counts toward the spanning `cover` (rule 4) |
| skipped codon crossing the origin of a circular replicon | the key is wrapped (M7-02 rule 3). A wrapped and an unwrapped hit of the same peptide count as **one** skipped site |
| linear contig | no difference apart from wrap |
| two hits of one folded peptide at one codon, one skipped and one `match` | REF counts. The pair is still counted once in `skipped_sites` (the skip happened) (Q6) |
| ties | no new ranking. `alt_counts` order and `variant_id` order are unchanged (rule 5). `examples` is sorted lexicographically and capped at 10. `by_residue` keys are sorted |
| empty input (no eligible peptides) | header-only outputs as today. `attrs["nonstandard"]` is all zeros. No event |
| lowercase letters | impossible: `peptides.tsv` is uppercase A–Z, and other rows are dropped in preflight |
| genome `X` (ambiguous codon) under a peptide residue | unchanged: impossible for an accepted hit (M7-01 degenerate table), raises if met |

## Acceptance criteria

- [ ] **Regression (the crash):** `run_variants` on a small fixture genome (c1 above), with a
  `peptide_hits` row for `AKXREFGHILWD` (n_diff 1, `diff_types` `substitution`) and a `peptides`
  row for it, returns without error. On main at 94bb0029 the same call raises
  `ValueError: no table-11 codon encodes 'X'`. Check: codon 6 is absent from the pileup, and
  `variants.tsv` and the records are empty.
- [ ] The same regression for `U` (example 6), for `B` and `Z` (over D/N and E/Q codons
  respectively), and for `O` over a sense codon: no crash, no row at the codon.
- [ ] Example 4: J over an I codon gives no ALT and no VCF record with ALT == REF. A property test
  checks that no VCF record anywhere has an ALT allele equal to its REF.
- [ ] Examples 1–10 and every edge-case row, as unit tests in `tests/unit/variants/` (pileup on
  hand-built `HitAlignment`s; examples 1–4, 6 and 9 also end to end through `run_variants`).
- [ ] Example 7 (U at TGA → REF) and example 8 (O at TGA, U at TAG → skipped), both on `+` and on
  `−` (gene-strand codon bases), and once with the codon crossing the origin.
- [ ] `eligible` under `genome_relation = same`, `database`, `search_db = other`: the readthrough-
  shaped op is kept, and a substitution to X is stripped and not counted.
- [ ] `attrs["nonstandard"]` matches rule 6 exactly (keys, types, sort, the `examples` cap of 10).
  It is all zeros on all-standard input.
- [ ] `variants_step` emits exactly one `variants.nonstandard_residues` event with that data when
  any site is skipped or read through, and **none** otherwise (`caplog` / the events JSONL).
- [ ] **No change for standard residues:** a property test checks that on random all-standard
  alignments the pileup, sites and VCF records are equal with and without the new code path. In
  practice that means equal to fixtures recorded on main. Every committed golden (and the M7-15
  variants goldens, once merged) passes **without** `GIAP_UPDATE_GOLDEN`.
- [ ] Property: adding a peptide whose only difference from an existing one is a nonstandard
  letter never changes any call, confidence, `n_ref`, `n_il`, `n_alt`, AD or DP.
- [ ] The M7-15 K-12 real-tool run (Slurm, `giap_run.sbatch`, recorded in PROGRESS.md by the lead)
  passes the `variants` step. The event reports a `skipped_peptides` count ≤ 13 (only the
  eligible, placed ones among the 13 TE8 peptides).
- [ ] Lint, format, mypy clean. ≤ 150 counted code lines.

## Test plan

- **Unit (`tests/unit/variants/test_pileup_nonstandard.py`, new):** examples 1–10 and the edge
  table, on hand-built `HitAlignment`s through `pileup` and `sites`, plus the `eligible` cases.
  The fixture is the c1/c1′ genome above, built inline: no committed FASTA.
- **Regression (`tests/unit/variants/test_variants_nonstandard_regression.py`, new):** the
  `run_variants` end-to-end rows (X, U, B, Z, O, J). Each test names the crash in its docstring.
  The test-author confirms they are red on main before the implementation.
- **Properties (`hypothesis`, in the same files):** standard-only invariance; adding a
  nonstandard twin is a no-op on counts; no VCF ALT == REF.
- **Step:** one `variants_step` test (existing step-test harness) asserting the event is emitted
  once and not emitted on an all-standard input.
- **Goldens:** rerun the golden suite unchanged on Slurm (`test_branch.sbatch`; the fast suite
  needs no local run on the shared 8-CPU session).
- **Real tools:** the K-12 `farmshare,apptainer` rerun belongs to the M7-15 real-tool arm. This
  issue only needs the variants step to pass there.
- No calibration change. No golden regeneration is allowed: a golden diff is a failure of
  rule 9.

## Line estimate (counted: `src`, `workflow`, `scripts`, added + deleted)

| Part | Lines |
|---|---|
| `variants/pileup.py`: constants, the rule-4 classifier (readthrough / skip), the counters, `attrs["nonstandard"]` | 50 |
| `variants/pileup.py`: `eligible` keep-condition | 5 |
| `steps/variants.py`: event constant and emission | 12 |
| docstrings (module headers of both files cite M7-18) | 8 |
| **Total** | **≈ 75** (target ≤ 150) |

Tests ≈ 220 lines, not gated. One PR.

## Questions for the lead

1. **Spec conflict: §7.5.1/§7.5.3 versus lead decision 9.** §7.5.1 row 1 says "Sec sites … appear
   as ALT calls", and §7.5.3 G8 `stop_codon` says "The REF residue is a stop (e.g., Sec read
   through UGA)". Decision 9 makes U over TGA a REF match instead. **Proposed reversible
   default:** follow decision 9 (REF). The `stop_codon` flag then fires only for a *standard*
   residue at a stop, for example W or C read through UGA. The case is unreachable until M7-16
   anyway. Because it reads against GIAP_SPEC.md, it probably needs an ADR, or a note in the M7-16
   ADR, if Khoa keeps it. The alternative is to record U at TGA as an ALT `U` with `stop_codon`
   and teach `minimal_alt_codon` that U = TGA and O = TAG. That is also about 10 lines.
2. **B/Z/J as ambiguous matches?** Should B over D or N (and Z over E/Q, J over I/L) count as REF
   (or IL for J)? **Proposed default:** no, skip them, as the mapper does (all `NOMATCH`, each a
   full substitution). Changing it belongs in Module 2 first (query codes, the q-gram index and
   the placement FDR's chance-hit table B.2). The variants step would follow automatically. Filed
   as a backlog item, not this issue.
3. **Event level.** INFO (proposed) or WARNING? Preflight already warns about the alphabet for
   proteins, and a WARNING may set the run status to `warn` (as `layer_missing:*` does).
   **Proposed default:** INFO. Should preflight also count B/J/O/U/X/Z in **peptide** tables, as
   it does for proteins? Proposed: yes, as a separate small preflight issue (not here).
4. **Do skipped codons count toward the frameshift-spanning `cover`?** **Proposed default:** yes.
   The peptide spans the codon contiguously, which is the evidence that rule 8 of M7-02 needs. It
   changes nothing on standard input.
5. **Fold a literal J into L for distinct-peptide counting?** Two peptides that differ only by
   J/L, or J/I, are one spectrum's ambiguity, yet would count as two distinct peptides elsewhere
   in the pileup. **Proposed default:** no change (fold I → L only), matching mapping, where `J`
   ≠ `L`. Revisit if the pilot shows J in real data. K-12 TE8 has at most 13 such peptides in
   total.
6. **Counting a skip that another hit of the same peptide resolves.** **Proposed default:** count
   it (the event reports what was skipped, not what was lost). The alternative counts only pairs
   left with no category.

## Lead decisions (2026-10-23)
1. **Q1 (U/O over a stop):** M7-15 decision 9's "U over TGA / O over TAG is a REF match" conflicts
   with §7.5.1 ("Sec sites … appear as ALT calls") and §7.5.3 (`stop_codon`). Neither case can
   reach the pileup today, because the DP never consumes a target `*`, and the spec is read-only.
   So M7-18 doesn't settle the question: **every U and O is skipped as `rare`**, and
   `readthrough_sites` is dropped from the event. Sec/Pyl readthrough calls go to the M7-16 ADR
   (option B), which will add the rule when that path becomes live. The `eligible` keep-condition
   for `ref == '*'` is not added. This supersedes the U/O-over-stop examples and that part of
   the estimate.
2. **Q2–Q6:** the proposed defaults are taken:
   - B, J and Z are skipped, as in the mapper;
   - the event is logged at INFO;
   - skipped codons count toward frameshift-spanning coverage;
   - a literal J is not folded into L;
   - a skip that another hit resolves still counts.

   Q2's ambiguous-match idea and Q3's preflight peptide count go to the backlog.
3. **Literal J over an I codon** (the silent REF = ALT bug): it is in scope. A J op is skipped
   before I/L → `J` folding is applied, and a regression test pins that no VCF record has
   ALT == REF.
