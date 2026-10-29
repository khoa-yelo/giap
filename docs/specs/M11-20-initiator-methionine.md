# M11-20: initiator M at annotated non-ATG starts (exact route)

Milestone M11 · pass 2 · depends on: M3-01 (six-frame index), M3-04 (exact route), M3-14 (loci
and flags), M7-16 (consumable stop cells, the model for this issue; ADR-0058), M11-16 (index
format 3, tables 4/11/25; ADR-0067), M7-01/M7-02 (variants alignment and pileup) · estimate
≈ 190 code lines + ≈ 420 test lines, one test PR then one implementation PR (table at the end) ·
**`stats`**: the decoy search changes (credited windows at mirrored cells), so the placement FDR
and the de novo precision inputs change; the statistician signs off on § Statistics · roles:
implementer (`src/gippa/mapping/`, `src/gippa/variants/`, `src/gippa/steps/mapping.py`,
`src/gippa/steps/variants.py`, `src/gippa/schemas/`, `validation/proxies/te_proxies.py`),
test-author (`tests/unit/test_m11_20_*.py`, `tests/property/test_m11_20_*.py`, the updated
expectations listed in § Tests whose expectations change, the calibration arm), statistician (§
Statistics, the arm) · source: M11-16f statistician review (PR #481), backlog "Initiator M at
non-ATG starts" (on #481), M7-13 decision 14 (on #481).

## Spec sections implemented

- GIPPA_SPEC.md §7.2.1 (six-frame index; storage), §7.2.2 (route 1, exact), §7.2.5 (placement FDR:
  decoys are searched "identically"), §7.2.6 (locus class, uniqueness), §7.5.2 (pileup), §9.4
  (Placement flag vocabulary), §11.6 TE1/TE2/TE8 and their in-silico proxies, §6.5 (0-based
  half-open, unwrapped origin-crossing intervals), §6.6 (no new randomness).
- **Spec conflict, needs an ADR (Khoa accepts).** §2.2 improvement 6 "Recoding and start catalogue
  (Sec, prfB, dnaX, non-ATG starts)" and §15 "No: … non-ATG starts, Met removal" are declined, and
  the TE1 proxy wording ("the exceptions being N-terminal peptides at non-ATG starts") accepts the
  miss. This issue is not a catalogue: it uses only the run's own Bakta annotation, as ADR-0058
  did for annotated stops and M3-12 decision 1 did for miniprot (`*gtgM` at query position 0
  counts as a match). The proposed **ADR-0068 "Initiator M at annotated non-ATG starts"**
  (written in the implementation PR) overrides §7.2.1/§7.2.2 (credit rule, index format 4), §9.4
  (new flag `initiator_m`) and the TE1/TE8 proxy wording (the exception no longer exists for
  annotated genes). Q-K1.

## What exists today (main 81a661bc)

- `build_sixframe` translates every codon with the run's table, so an annotated GTG start reads
  `V` and TTG reads `L` (`mapping/sixframe.py`; tables 4, 11 and 25 agree on GTG/TTG).
- `match_exact` (`mapping/exact.py`) is a plain Aho–Corasick scan. Bakta writes `M` as the first
  residue of every complete CDS whatever its start codon, so the N-terminal tryptic peptide of a
  GTG/TTG gene (when the Met is retained) has no exact hit at its own gene.
- What happens to that peptide now:
  1. 7–9 aa: below `mapping.tolerant.min_len` (10): `unmapped` `no_hit` — **lost**.
  2. ≥ 10 aa: tolerant route at its own gene, `n_diff` 1 (substitution M/V or M/L), placement
     stratum 1 (π̂ = 1), `strain_difference` under `genome_relation = related`, not counted
     under `same`; variants: ALT `M` at the start codon, flag `start_codon` (G8, M7-05).
  3. Any length, when the same M-initial peptide also occurs exactly at a real ATG elsewhere (a
     chance copy): the exact route wins, the tolerant route never runs for it ("only failures
     move on"), and the copy is its only locus. If that locus is intergenic, out of frame or
     antisense, `novel/candidates.py` takes it as a unique candidate: an **unestimated false
     novel ORF** (M7-13 decision 14: 5 ORFs in R = 300 at code 25, counted as
     `nonatg_nterm_recurrence`), or a false placement/assignment
     (`engine_protein_mismatch` with same-strain input).
- Scale: GIPPA_SPEC §16 sim defaults put 10 % of starts at GTG/TTG; real bacteria 10–15 %. Met is
  removed in vivo when residue 2 ∈ {A, C, G, P, S, T, V}; then the observed peptide has no M and
  maps today. The gap is the Met-retaining half of the non-ATG genes. TE1's Nissle proxy lists 431
  such distinct N-terminal products of 164,181 peptides (`tests/integration/test_te1_proxy.py`).
- Precedents in the codebase: `readthrough.StopSite` / `frame_cells` / `mirror_cells` (annotated
  sites → per-frame cells, mirrored into decoys, stored in the index, ADR-0058); miniprot's
  `alt_start` (M3-12 decision 1); `variants.flags.start_codons` (annotated first codons).

## Options

| # | Option | Fixes recall | Fixes false novel / placement | Index / FDR impact | Size |
|---|---|---|---|---|---|
| **A** (recommended) | **Credit M at annotated non-ATG start cells, exact route only.** A query whose first residue is M matches a window whose first position is an annotated non-ATG start cell when its other residues match exactly; decoys get the mirrored cells (the M7-16 model) | yes, all lengths 7–40 | yes: the peptide gets a `cds_inframe` row at its own gene, so `candidate_rows` excludes it (M7-06 rule 1) and a chance copy makes it `shared` | format 4 (start sites + cells); decoys searched identically; ≈ 0.2 % more chance windows for M-initial peptides | ≈ 190 code lines |
| A′ | A, plus credit in the tolerant route (position 0 free at a cell) | also cost-2 cases | same | the tolerant DP and q-gram seeds change; decoy tolerant search too | ≈ +150, touches M3-06 kernels |
| A″ | Credit M at **every** in-frame start codon (GTG/TTG anywhere), not only annotated ones; also covers novel ORFs with non-ATG starts | yes, also novel ORFs | partly; adds false loci | ≈ 2/64 of all positions become cells: the first-residue match probability for M-initial peptides roughly doubles (f_M ≈ 0.025 vs f_GTG + f_TTG ≈ 0.03); E_i must change; most M-initial tryptic peptides are internal (after K/R) and would gain chance loci at internal GTG/TTG | large; a statistics change |
| B | Annotated-protein targets: search Bakta's protein strings as extra targets (and their reverses as decoys) | yes | yes | N grows ≈ 15 % with fully redundant positions; every CDS peptide hits twice (dedupe protein ↔ genome); new decoy space; new index layout | large |
| C | Backlog option 2: `candidate_rows` drops a peptide whose engine `protein_ids` name an annotated protein | **no** | consequence 2 only, database input with same-strain IDs only (not de novo, not cross-strain IDs) | none | small, but a patch on the symptom |
| D | Status quo: keep the TE1 exception and the decision-14 diagnostic | no | no | none | 0 |

**Recommendation: A.** It is the smallest change that is correct for both consequences, it uses
only the annotation (like ADR-0058 and miniprot's `alt_start`), and its decoy treatment is the
already-reviewed M7-16 cell mirroring. The tolerant route (A′) and novel-ORF starts (A″) are
follow-ups (backlog). C and D leave the recall loss. B costs an index redesign for the same gain.

## Definition (option A)

### 1. Start sites

`gippa.mapping.initiator` (new, pure):

```python
@dataclass(frozen=True, slots=True)
class StartSite:
    contig: str
    strand: Literal["+", "-"]
    codon_start: int              # lowest coordinate of the first codon, wrapped into [0, L)
    locus_tags: tuple[str, ...]   # sorted; features sharing the codon
    codon: str                    # gene-strand codon as read, e.g. "GTG"

@dataclass(frozen=True, slots=True)
class StartCells:
    sites: tuple[StartSite, ...]
    target: tuple[np.ndarray, ...]                    # per frame id, sorted int64 aa positions
    decoys: Mapping[DecoyType, tuple[np.ndarray, ...]]

def annotated_starts(gm, genome, start_codons, *, table) -> tuple[tuple[StartSite, ...], dict[str, int]]
def start_frame_cells(index, sites) -> tuple[np.ndarray, ...]
def mirror_start_cells(kind, frames, cells) -> tuple[np.ndarray, ...]
def initiator_rows(rows, sites, lengths, circular) -> list[bool]
```

A **start site** is the first codon (transcription direction: `start` on `+`, `end − 3` on `−`)
of a `cds` or `sorf` feature (never `pseudogene`) when all hold:

1. the gene's Bakta `aa` is non-empty and begins with `M`;
2. the genome codon there, read on the gene strand (wrapping on a circular replicon), is in
   `annotation.start_codons` (upper-cased; default ATG, GTG, TTG) **and is not `ATG`**;
3. on a linear contig the codon lies inside `[0, L)`.

Features failing a rule are counted, not raised: `atg` (rule 2, ATG; the normal case),
`not_a_start` (rule 2, another codon, including a codon with a non-ACGT base), `aa_not_m`
(rule 1), `off_contig` (rule 3). Features sharing a key `(contig, strand, codon_start)` give one
site with sorted `locus_tags`. Sites are sorted by (contig order, strand, codon_start), as
`annotated_stops`. ATG starts are **excluded** on purpose: they already read M, and a cell there
would add decoy windows with no target change (edge case "ATG start").

### 2. Cells

`start_frame_cells` maps each site to its frame position `p` exactly as `readthrough.frame_cells`
does (the `+` and `−` formulas, plus the copy at `codon_start + L` on a circular replicon when it
fits in the extended frame), and raises `ValueError` when the frame residue at `p` is not the
table's translation of `codon` (it can never be `*` or `X` for a start). `mirror_start_cells`:
`n − 1 − p` in a `reverse` decoy frame of `n` residues, `p` in a `shuffle` decoy (the cell is a
position; the residue the shuffle put there does not matter). No randomness.

### 3. The credit rule (exact route only)

`match_exact(..., starts: StartCells | None = None)`. With `starts`, in addition to today's hits:

> a query `q` (after I → L) of length `k` with `q[0] == "M"` matches frame `f` at window
> `[p, p + k)` when `p` is a cell of `f` (target: `starts.target[f]`; decoy kind `t`:
> `starts.decoys[t][f]`) and `q[1:] == frame[p + 1 : p + k]`.

- The credited row is a normal exact row: `route = exact`, `n_diff = 0`, `diff_types = ""`,
  `blocks` NA; target or decoy (`decoy_type`) as today; the same `3k ≤ L` guard, mirroring, and
  deduplication on the (mirrored) interval as `_scan`, so an origin copy is one row and a window
  found both plainly and by credit is one row.
- Only query position 0 is credited (an initiator is the protein's first residue). An M at any
  other query position never matches a cell.
- Queries whose first residue is not M are unaffected. `starts=None` (or no sites) gives today's
  output exactly.
- The hit table schema (`HIT_COLUMNS`) is unchanged; whether a row was credited is derived later
  (rule 5), as `readthrough` is derived by `consumed_sites`.
- Every `match_exact` call in `steps/mapping.py` passes the index's start cells: `map_exact` and
  `_reference_exact` (de novo `c_ℓ`). `_denovo` reads `hits_exact`, so `m_bℓ` and `r_bℓ` include
  credited rows.

### 4. Tolerant route

Unchanged (no credit). A credited peptide has a target exact row, so it is not in `left` and its
target tolerant search does not run; its decoy tolerant search runs as today. A peptide with no
exact hit (e.g. M-initial with one more difference) goes to the tolerant route as today, which
charges the M/V substitution (cost 2 for a real SNV plus the initiator) and still yields the
ALT `M` + `start_codon` variant. A′ is the follow-up.

### 5. Flag `initiator_m`

New §9.4 Placement flag, added to the `peptide_hits.tsv` `flags` pattern (additive
`OUTPUT_SCHEMA_VERSION` minor bump, M11-14h policy). `annotate_peptide_hits(..., initiator=…)`
sets it on a row when `initiator_rows` is true: `search = target`, `route = exact`, the raw
peptide's first letter is `M`, and the row's 5′ codon (`start` on `+`, `end − 3` on `−`, wrapped)
is a start site key. It does not change locus class, gene, uniqueness or any other flag. A
credited row at the gene's start is in frame inside that gene, so it is `cds_inframe` (or `sorf`)
by M3-14 as is.

### 6. Variants (`variants/align.py`, `variants/pileup.py`)

`align_hit` today raises for an `n_diff = 0` row whose letters differ beyond I/L, so it must learn
the rule: `align_hits(..., starts: frozenset[(contig, strand, codon_start)])`. At an exact row
whose residue 0 is `M` over a start-site codon, residue 0 is op `initiator` (new `OpKind`) and the
rest is the zip. `pileup` gives an `initiator` op **no category** at its codon (it is not
evidence about the genome residue: the initiator tRNA reads GTG as M by design) but the codon stays
in `cover`; the count goes to `attrs["initiator"]`. So a credited peptide never makes REF V, ALT
M or `mixed` at the start codon, and supports REF at its other codons. Q-L3.

### 7. Index format 4 and cache keys

- `_INDEX_FORMAT = 4`. `meta` gains `"start_sites"` (list of `asdict(StartSite)`); arrays
  `target_init`, `target_init_bounds`, `<kind>_init`, `<kind>_init_bounds` (packed exactly as the
  stop cells, `_pack_cells`). `save_index(path, index, decoys, cells, starts=None)`; `None`
  writes empty site and cell lists. New `load_start_cells(path, *, table) -> StartCells`.
- Loading: format 4 as above. Formats 2 and 3 load (unchanged rules, D3) with **no start cells**
  and `load_start_cells` logs one WARNING event `mapping.index_without_start_cells` ("built before
  M11-20; rebuild the six-frame index for initiator-M credit"). Format 1 is refused as today. The
  committed fixtures `sixframe_format1.npz` / `sixframe_format2.npz` and their tests stay valid.
  Q-L1.
- `sixframe_index` already takes `genes` and `genome` as inputs and `annotation` (holding
  `start_codons`) is in every step's config slice (`steps/registry.py`), so **no new input, no new
  config key and no task-hash change**: Nextflow `-resume` re-runs `sixframe_index` exactly when
  the genes, genome, peptides or config change, as today. A cached index from an older GIPPA is
  reused only together with its cached downstream tasks (the task hash does not include the
  GIPPA version, as for every code change); a mixed state (hand-run `gippa step`) gets the
  warning above.
- `sixframe_index` logs `mapping.start_sites` with `data = {sites, cells_target, cells_decoy:
  {kind: n}, skipped: {atg, not_a_start, aa_not_m, off_contig}, translation_table}`.
- No config toggle (MVP). The A/B calibration arm calls `match_exact(starts=None)` directly. Q-L2.

### 8. Genetic codes 4 and 25

No table-specific rule: GTG → V and TTG → L in 4, 11 and 25; the start set is
`annotation.start_codons` at every table (M11-16 "alternative starts"). The `start_frame_cells`
check uses the table's own letter, so a user-extended set (e.g. `ATT` under code 4, read `I`,
stored `L`) works. TGA inside the peptide reads W (4) or G (25) as today.

## Statistics (for the statistician; `stats`)

- **S1 Decoy symmetry.** Target and decoys get the same number of cells; a credited window needs
  the same thing in both (one cell, then `k − 1` exact residues), so under the null of a wrong
  peptide the target and decoy credited-hit rates are equal. Edge loss differs slightly: a target
  cell near a frame's 3′ end and a reverse-decoy cell near its 5′ end lose windows; with ≤ a few
  hundred cells in ≈ 10⁷ positions this is negligible. Please confirm.
- **S2 Size.** Extra chance windows for an M-initial wrong peptide of length k: `C · Π_{i≥1} f_i`
  against `N · f_M · Π_{i≥1} f_i` today, ratio `C / (N f_M)`. K-12: C ≈ 450 sites (≈ 10 % of
  4,300 CDS, ×1–2 copies only at the origin), N ≈ 9.3 × 10⁶, f_M ≈ 0.025 → ≈ 0.2 %.
- **S3 Ranking score.** E_i (`chance.expected_hits`) is unchanged: it is a ranking score and the
  decoys calibrate (§7.2.5). Credited hits are stratum m = 0; π̂ rules unchanged; `same_db` exact
  rows keep `placement_q = 0`.
- **S4 De novo precision.** `m_bℓ`, `r_bℓ` (from `hits_exact`) and `c_ℓ` (`_reference_exact`) all
  use the credit, so the estimator stays internally consistent; `c_ℓ` rises slightly (reference
  N-terminal peptides of non-ATG genes now match).
- **S5 Novel FDR.** Code unchanged. M7-13 decision 14's `nonatg_nterm_recurrence` cases leave the
  candidate set (their peptide gets a `cds_inframe` row). The diagnostic stays and is expected
  to be 0; retiring the amendment is the statistician's call after the arm. Q-S2.
- **S6 Arm.** See § Tests, calibration. Q-S1.

## Edge cases

| Case | Behaviour |
|---|---|
| **ATG start** | no site (counted `atg`); the peptide matches as today; no extra decoy windows; the row has no `initiator_m` |
| **GTG/TTG start on `+`** | site at `start`; credited row `[start, start + 3k)` `+` |
| **minus strand** (TTG on `−`, `CAA` on `+`) | site at `end − 3`; cell from the `−` formula; credited row `[end − 3k, end)` `−`, frame label from `frame_label` |
| **start codon on a circular replicon, peptide across the origin** | cell at `p`, window extends into the wrap; row unwrapped (`end > L`), flags `initiator_m;origin_spanning` |
| **start codon itself across the origin** (`+` gene at `L − 1`: bases L−1, 0, 1) | `codon_start = L − 1`; codon read with wrap; cell in the frame of offset `(L − 1) % 3`; row `[L − 1, L − 1 + 3k)` |
| origin copy of a cell (`codon_start + L` fits in the extended frame) | both cells credited; both windows map to one interval: one row |
| **linear contig** | no wrap; a site at `[0, 3)` on `+` or `[L − 3, L)` (first codon of a `−` gene) is valid; a first codon beyond the contig is `off_contig` |
| **internal GTG that is not a start** | no site: an M-initial peptide there does not match (reads V) |
| **M at query position ≥ 1 over a cell** (e.g. an N-terminal-extension peptide spanning the start) | no credit at that position; the peptide reading V there matches as today (`nterm_extension`) |
| query starting with V/L at a cell | matches as today (the cell is not consumed or replaced) |
| **tolerant hits** | no credit; an M-initial peptide that also needs a real difference is charged the M/V substitution as today (ALT M, `start_codon`) |
| **chance ATG copy elsewhere (ties)** | two target rows: credited `cds_inframe` at the gene and plain exact at the copy; `uniqueness = shared` (or `multi_identical` if the copy is in a gene with the identical protein), `multi` when the classes differ; the peptide leaves the novel candidate set (M7-06 rule 1) |
| two features sharing one start codon | one site, `locus_tags` sorted |
| two identical non-ATG paralogues (IS copies) | two credited rows, `multi_identical` |
| gene `aa` empty or not starting with M (5′-partial, Bakta `start_type = Edge`) | no site (`aa_not_m`) |
| codon with N | `not_a_start`; the frame holds `X` there |
| user removes GTG from `annotation.start_codons` | GTG-start genes get no site (`not_a_start`); Q-L4 |
| **codes 4 and 25** | same rule; the cell's residue is the table's letter; TGA inside a peptide is W / G |
| **empty input**: no genes, or no non-ATG starts | `start_sites = []`, every cell array empty; output identical to `starts=None` |
| empty input: no peptides | `hits_frame()` as today; the index still records the sites |
| shuffle decoy | the cell keeps its position; credited decoy window `[p, p + k)` regardless of the shuffled residue at `p` |
| reverse decoy cell `n − 1 − p` too close to the frame end for `k` | no window (no error) |
| format-3 index loaded by the new code | no start cells, WARNING `mapping.index_without_start_cells` |
| `same_db` input | credited rows are exact: `placement_q = 0`, `same_db` |
| `genome_relation = same` | credited rows count (they are exact); before, the tolerant row did not |

## Examples (input → expected output)

Coordinates 0-based half-open; table 11 unless stated; `mapping.exact.min_len` 7; frame ids are
`SixFrameIndex.frames` order (`+0 +1 +2 −0 −1 −2` per contig). Examples 1–8 were checked against
`build_sixframe` / `build_decoys` on main 81a661bc.

**Genome c1** (circular, 120 nt): `C×12`, then gene **G1** `+` `[12, 39)` =
`GTG GCT GAA CTG TGG CAT AAA GGC TAA` (Bakta `aa` = `MAELWHKG`), then `C×81`. Frame 0 (c1 `+0`,
80 aa with the wrap) reads `PPPPVAELWHKG*PPP…PPPPVAELWHKG*PPP…`; start site `(c1, +, 12)`, cells
{4, 44} (44 = the wrap copy).

| # | Input | Expected (after M11-20) | Today |
|---|---|---|---|
| 1 | c1, query `MAELWHK` (database, reverse decoy) | `annotated_starts` → one site `(c1, +, 12, ("G1",), "GTG")`, skipped all 0. One target row: frame_id 0, aa `[4, 11)`, `[12, 33)` `+`, frame 1, route exact, n_diff 0. Annotated: `cds_inframe`, gene G1, `unique`, flags `initiator_m` | no target row; `unmapped` `no_hit` (7 aa < tolerant `min_len`) |
| 2 | c1 with `[60, 81)` replaced by `ATG GCT GAA CTG TGG CAT AAA` (a chance copy), query `MAELWHK` | two target rows: `[12, 33)` `cds_inframe` G1 `initiator_m`, and `[60, 81)` `intergenic` (no `initiator_m`); both `shared` and `multi`; `candidate_rows` returns no row for the peptide | one row `[60, 81)` `intergenic` `unique` → a novel candidate (the decision-14 case) |
| 3 | example 1 with G1's first codon `ATG` | no site (`skipped.atg = 1`); one exact row `[12, 33)` as today, flags without `initiator_m`; every table byte-identical to today except `sixframe.npz` meta | same row |
| 4 | c1, query `PPMAELWHK` (M at position 2) | no target row | no target row |
| 5 | c1, query `PPVAELWHK` | one row `[6, 33)` `+`, locus class as today (M3-14), no `initiator_m` | same |
| 6 | **minus strand.** c2 linear 60 nt, gene G2 `−` `[6, 30)`; on `−` it reads `TTG AAA CAG TGG TAT TTT CGT TAG` (`aa` `MKQWYFR`; `+` strand `[27, 30)` = `CAA`); query `MKQWYFR` | site `(c2, −, 27, ("G2",), "TTG")`; cell 10 in frame c2 `−0`; one row aa `[10, 17)`, `[9, 30)` `−`, frame −1, `cds_inframe` G2, `initiator_m` | no row (`no_hit`) |
| 7 | **origin.** c3 circular 90 nt, gene G3 `+` `[87, 111)` unwrapped, `TTG GCT GAA CTG TGG CAT AAA TAA` (`aa` `MAELWHK`); query `MAELWHK` | site `(c3, +, 87)`; cell 29 in c3 `+0`; row `[87, 108)` `+`, flags `initiator_m;origin_spanning` | no row |
| 7b | **start codon across the origin.** c4 circular 90 nt, G4 `+` `[89, 113)`, codon = bases 89, 0, 1 = `GTG`; query `MAELWHK` | site `(c4, +, 89)`; cell 29 in c4 `+2`; row `[89, 110)` `+`, `initiator_m;origin_spanning` | no row |
| 8 | **reverse decoy.** c1, query `MPPPPPP` with decoy type `reverse` | decoy cells {75, 35} in frame 0's reverse (`n = 80`); window at 35 reads `VPPPPPP` → one **decoy** row: frame_id 0, aa `[35, 42)`, mirrored target interval `[114, 135)` `+`; the window at 75 runs off the frame (no row); no target row | no decoy row |
| 9 | **tolerant.** A GTG gene whose N-terminal tryptic peptide is 12 aa, query = that peptide with one more substitution at residue 6 | no exact row; tolerant row at the gene, n_diff 2 (M/V + the substitution), as today; variants ALT `M` with `start_codon` at the start codon | same |
| 10 | **variants.** Example 1's G1, `search_db = same`, `genome_relation = same`, plus `MAELWHKG`-covering peptides | pileup: no category from the credited peptide at codon 12 (`attrs["initiator"] = 1`), REF at codons 15…30; no `variants.tsv` row at codon 12 | (for a ≥ 10-aa N-terminal peptide) tolerant row, ALT `M` `start_codon` at codon 12 |
| 11 | **code 4.** Gene `+` `GTG GCT GAA TGA CTG CAT AAA …` (`aa` `MAEWLHK…`), table 4, query `MAEWLHK` | credited row, `initiator_m`; at table 25 the `aa` and query are `MAEGLHK` and the row is the same | no row |
| 12 | **shuffle decoy.** c1, de novo query `MAELWHK` | target row of example 1; decoy cells {4, 44} (positions kept); a decoy row exists iff the shuffled residues after a cell spell `AELWHK` | target: none |
| 13 | **empty.** A gene model with only ATG starts, or no genes | `start_sites = []`, `cells_target = 0`; `match_exact(starts=…)` equals `match_exact(starts=None)` row for row | same |

## Acceptance criteria

1. `annotated_starts` returns exactly the sites and skip counts of examples 1, 3, 6, 7, 7b and
   13, and merges two features sharing a codon into one site with sorted `locus_tags`.
2. `start_frame_cells` gives the cells of examples 1, 6, 7, 7b (both origin copies when they fit)
   and raises `ValueError` on a site whose frame residue is not the codon's translation.
   `mirror_start_cells` gives `n − 1 − p` (reverse) and `p` (shuffle) and draws no randomness.
3. `match_exact(..., starts=…)` gives the rows of examples 1, 2, 4–8 and 12; with `starts=None`
   or empty sites its output is byte-identical (parquet) to main's for the same inputs.
4. `annotate_peptide_hits` sets `initiator_m` exactly on credited target rows (examples 1, 2, 6,
   7, 7b), never on decoy, tolerant, ATG-start or position-≥ 1 rows; classes and uniqueness as in
   examples 1, 2 and the paralogue edge case.
5. `peptide_hits.tsv` with `initiator_m` validates against the bumped schema; `OUTPUT_SCHEMA_VERSION`
   has one additive minor bump.
6. Index format 4 round-trips sites and cells at tables 4, 11 and 25; format 3 and the committed
   format-2 fixture load with no start cells and the WARNING event; format 1 is refused as today.
7. `sixframe_index` logs `mapping.start_sites` with the counts of § 7.
8. Variants: example 10's pileup and sites; `align_hit` no longer raises on a credited row and
   still raises on any other `n_diff = 0` letter mismatch.
9. **No non-ATG starts → no change.** On a gippasim run with `start_codons=(("ATG", 1.0),)` every
   table under `tables/` and `mapping/` is byte-identical to main; only `sixframe.npz` differs
   (format 4, empty sites).
10. **Recall.** `tests/calibration/test_exact_recall.py` with N-terminal peptides of non-ATG genes
    included (the `off == 0 and start_codon != "ATG"` exclusion removed): recall ≥ 0.999 at the
    existing tolerance; the non-ATG N-terminal subset alone: 100 % at their own gene.
11. **TE1 proxy.** `te1_proxy` (Nissle, `validation/proxies/te_proxies.py`) without the
    `non_atg_nterm` exclusion: the 431 products are `exact_inframe`; the share ≥ 0.999 holds
    with them in the denominator. The exclusion is kept only as a reported diagnostic count
    (`non_atg_nterm`, now "mapped by initiator credit").
12. **Calibration arm** (§ Tests) passes at the existing tolerances.
13. Lint, format, mypy clean; no new dependency; no randomness.

## Tests

**Unit** (`tests/unit/test_m11_20_initiator.py`): examples 1–13; the skip counters; format 4
round-trip and the format-2/3 warning path; `initiator_rows`; `align_hit`/`pileup` example 10.

**Property** (`tests/property/test_m11_20_initiator_props.py`, hypothesis, random genomes ≤ 20
kb, tables 4/11/25, both strands, circular and linear):

- P1 completeness: for a random gene with a GTG/TTG start and its `aa` (M + translation), every
  prefix of `aa` of length 7–40 that stays inside the gene maps exactly with a row at
  `[start, start + 3k)` (or the `−` equivalent) flagged `initiator_m`.
- P2 brute force: `match_exact(starts)` rows = a brute-force scan that tries, at every cell, `M`
  plus the next `k − 1` residues, unioned with the plain scan (target and both decoy kinds).
- P3 no-op: with no non-ATG sites the output equals `starts=None`.
- P4 decoy mirror: for `reverse`, a credited decoy window maps (`mirror_aa`) to the target
  interval ending at the cell's codon; the number of decoy cells equals the number of target
  cells per frame.
- P5 strand symmetry: reverse-complementing the genome (and the gene model) swaps the strand of
  every credited row and keeps the peptide set.

**Calibration** (`tests/calibration/test_m11_20_initiator_calibration.py`, `calibration`, fast
tier R = 50 at 100 kb; `calibration_full` R = 300 at 1 Mb, as §8.4 placement FDR):

- Arm I1, placement FDR with credit: gippasim with start codons ATG 0.5 / GTG 0.35 / TTG 0.15
  (to load the cells) and the existing wrong-peptide models, with a share of the wrong peptides
  forced M-initial (half); stratum-0 FDP ≤ q + 2 SE at q ∈ {0.01, 0.05}, database and de novo.
  Also record T₀, D₀ with and without credit (`starts=None`) and the share of credited decoy rows
  (expected ≈ C / (N f_M), S2).
- Arm I2, novel: re-run the M7-13 code-25 G4 arm (the one that showed 5 ORFs in R = 300):
  `nonatg_nterm_recurrence` = 0 in every replicate; per-class FDP at level as before.
- The statistician may widen or narrow these (Q-S1). No tolerance changes (R3).

## Tests whose expectations change (test-author updates, listed so nothing breaks silently)

- `tests/integration/test_te1_proxy.py`: 431 `non_atg_nterm` now map; counts and wording.
- `tests/integration/test_te8_proxy.py`: the 4 `non_atg_nterm` K-12 products map.
- `tests/integration/test_te2_proxy.py::test_te2_proxy_on_gtg_starts_flags_every_high_site_start_codon`:
  `high_start_codon > 0` no longer holds (the N-terminal peptides map exactly; no ALT M). The test
  becomes "no ALT M at a credited start; `violations` 0".
- `tests/integration/test_m7_04_variants_run.py` example 3: `variants.tsv` rows only at annotated
  non-ATG starts → no rows from N-terminal peptides that are credited (rows remain only for
  M-initial peptides that also needed a tolerant difference).
- `tests/integration/test_m10_05_sim_pilot.py`: the planted non-ATG N-terminal product is no
  longer a TE1 miss.
- `tests/calibration/test_exact_recall.py`: the non-ATG exclusion goes (AC 10).
- M7-13 G2 `start_codon` high counts (recorded, not asserted) fall.
- `tests/helpers/m7_gate_novel_records.py` `nonatg_nterm_recurrence` (on #481): kept; expected 0.
- Golden files regenerated in the implementation PR (reviewed diff): `tests/data/golden/sim50k/`
  `mapping/` (and `m4/`, `m5/` if a credited gene changes detection), `sim50k_m7/` (`pileup`,
  `variants`, `novel_orfs`), `sim100k_real/m4` — only where the fixture has a Met-retaining
  N-terminal peptide of a non-ATG gene.

## What changes in outputs (for the K-12 byte-identity check)

Upstream of mapping nothing changes: annotation, codon usage, regions, promoters, terminators,
operons, RBS and every annotation-only table are byte-identical. `sixframe.npz` always changes
(format 4, `start_sites`, init cell arrays). With peptides, let **A** = the input peptides whose
first residue is M and that are an N-terminal product of a K-12 CDS/sORF start site (non-ATG,
`aa` M-initial). Expected:

| Output | Change |
|---|---|
| `mapping/peptide_hits.tsv`, `mappings.bed`, `mappings.gff3`, `mappings.probed` | new exact rows for A at their genes, flag `initiator_m`; for A, tolerant rows at the gene (n_diff 1, `strain_difference`) disappear; rows at chance ATG copies become `shared`/`multi`; rows of peptides not in A unchanged except `placement_q` (below) |
| `mapping/unmapped.tsv` | A peptides leave `no_hit` / `rejected_placement_fdr` |
| `mapping/mapping_fdr.tsv`, placement work files | T₀ rises by the placed A; D₀ by the credited decoy hits; stratum-1 counts lose A; so `placement_q` of other peptides may shift (`search_db = same` exact rows keep q = 0, so with same-strain input every non-A row is byte-identical) |
| gene table and everything downstream (detection, abundance, detectability, context stats, enrichment, determinants, report) | changes only through the genes of A gaining detection / unique peptides, and peptides of A at chance copies no longer counted elsewhere |
| `tables/variants.tsv`, `pileup.tsv.gz`, VCF | fewer `start_codon` ALT M rows at annotated starts (MVP-04 K-12 audit: 6 such rows); pileup loses those ALT entries |
| `novel_orfs.tsv`, annotation patch | ORFs supported only by A peptides at chance copies disappear; novel-ORF non-ATG starts unchanged (MVP-04's 9 rows stay) |
| `logs/events.jsonl` | new `mapping.start_sites` event |
| `run_manifest.json` | `output_schema_version` minor bump |

The check passes when every changed `peptide_hits` row's peptide is in A, or differs only in
`placement_q` (and, with `search_db = same`, there is none of the latter).

## Estimate and PR split (each ≤ ~400 code lines; tests first, ADR-0007)

| PR | Content | Code | Tests | Role |
|---|---|---|---|---|
| 20t | unit + property tests, calibration arm, updated expectations (tests fail until 20i) | 0 | ≈ 420 | test-author |
| 20i | `mapping/initiator.py`, `exact.py` credit, `sixframe.py` format 4, `loci.py` flag, schema + version, `variants/align.py` + `pileup.py`, `steps/mapping.py` + `steps/variants.py` wiring, `te_proxies.py`, ADR-0068, changelog, golden regeneration | ≈ 190 | — | implementer |

The variants change must ship with the credit (a credited row makes today's `align_hit` raise).

## Questions

For Khoa:

- **Q-K1.** Accept ADR-0068 overriding §2.2 improvement 6 / §15 "No" for this narrow case
  (annotation-only, exact route only)? Default if unanswered: hold the implementation PR; the
  spec and tests can merge. This is the only blocking question.

For the lead:

- **Q-L1.** Format 2/3 indexes: load with no start cells plus a warning (default; keeps the
  committed fixtures and silently-old behaviour visible in the log), or refuse with "rebuild"
  (strict, but edits the format-2 tests)?
- **Q-L2.** No config toggle (default). A `mapping.exact.initiator_m: true` key would let users
  reproduce pre-M11-20 output, at the cost of a schema change and a one-time full `-resume`
  invalidation.
- **Q-L3.** The credited codon in the pileup: no category (default) or REF?
- **Q-L4.** Site rule 2 uses `annotation.start_codons`. Alternative: any non-ATG first codon of a
  complete gene whose `aa` begins with M (trust Bakta entirely). Default: the configured set.
- **Q-L5.** Follow-ups to the backlog: A′ (tolerant credit) and A″ (novel-ORF non-ATG starts, which
  the mapping step cannot know); makedb's literal translation of non-ATG starts (backlog line
  1196) is untouched here.

For the statistician:

- **Q-S1.** Is arm I1 (forced M-initial wrong peptides with loaded GTG/TTG starts) plus arm I2
  enough, or is a full-tier R = 300 placement arm needed before merge?
- **Q-S2.** After I2 shows 0, retire M7-13 decision 14's exemption or keep it as a diagnostic?
- **Q-S3.** Confirm S1–S4: E_i unchanged; decoy cells mirrored as M7-16; de novo precision uses
  the credit in m, r and c alike.
