# M7-13c: The M7 gate evidence — full-tier runs, methods-note calibration logs, the PROGRESS gate record

Milestone M7 · track C (calibration) · **pass 2 — the gate, part 3 of 3** (docs/specs/M7-13-calibration-gate.md
§ Split) · depends on: § Dependencies · estimate ≈ 0 `src/` lines, ≈ 150 docs lines, ≈ 150 test
lines if Q1 is accepted, ≈ 70 allocated CPU-hours on Slurm · **stats** (the statistician runs and
signs off the evidence) · **Slurm**

Written 2026-10-23 against main 3be808dc, origin/m7/calibration-gate-variants 06724e53 (#255,
M7-13a), origin/m7/calibration-gate-novel ce802349 (M7-13b), origin/m7/variant-flags-adr49 67421421
and Khoa's ADR-0049 decision (local/handoff/adr0049_decision.md, 2026-10-23).

## Spec sections implemented

- §13 row M7 gate: "F1 ≥ 0.9 on simulated variants and novel ORFs; TE2 proxy passes".
- §8.4: "The full suite runs nightly on FarmShare as a Slurm job array, one suite per task"; the
  per-PR subset timings are recorded "in `benchmarks/baseline.json`".
- §8.5: "A failing attack blocks the merge" (attack 6 as restated by ADR-0049 decision 5).
- §14.7: `docs/methods/` holds "one note per statistical method with its simulation evidence; the
  statistician signs off"; `PROGRESS.md` holds "done … questions for Khoa".
- §14.9, definition of done for M7: "A methods note exists for each statistical method, signed off
  by the statistician", "The benchmark is recorded", "Coverage is ≥ 90 %", "`PROGRESS.md` is
  updated".
- M7-13 acceptance items 4–6 (full-tier results with job IDs; the gate listed check by check; the
  backlog closures). M7-13a lead decisions 6 and 7 (the task-13 cost measurement and the
  statistician's methods text both belong to this issue). M7-index decision 11 ("the M7 gate is
  claimed check by check … PROGRESS.md records any vacuous clause").

This issue makes no new statistical claim and changes no method. It runs the merged gate tests at
the full tier on main and writes down what they say.

## What changed since M7-13 was written

1. **ADR-0049, decided by Khoa on 2026-10-23.** The decision changes three checks before any
   evidence is collected:
   - **G1:** scored on *unflagged* `high` substitution calls against *identifiable* truth. A truth
     row leaves the recall denominator when
     `giap.variants.flags.mass_equivalent(ref_aa, alt_aa, variants.mass_equivalent_tol_da)` is
     true. Flagged calls are reported and not scored. `f1_min.variants_high` stays 0.90 (R3).
   - **A6a:** every isobaric non-swap ALT is flagged `mass_equivalent`. This now includes the left
     2→1 block reads.
   - **A6b:** every random adjacent-swap call is flagged. This replaces "no `high` call".
   - **A6c:** unchanged. It is still expected to fail.
   - The statistician's evidence for the decision (Slurm 1874972, scope (d)) puts G1 at 0.934 at
     100k and 0.933 at 1M. It is not gate evidence, because it ran on branch code before the
     flag PR.
2. **ADR-0048 is still pending Khoa.** The novel gate uses `_m7_gate.GATE_R` (10 / 300) under
   suite `novel_orfs`, which is array task 11. If Khoa accepts the ADR, the suite becomes
   `novel_orfs_gate`, task 15. R = 300 in both cases.
3. **M7-13b lead decisions 1–8.**
   - Attacks a, b, d and GC run in the fast tier only.
   - Attacks a and d run on the clustered-null arm. Their "pass" is a confirmed fail of the
     mutant.
   - G4-k3 is a plain test until the full tier rules on it (decision 4).
   - Calls made through the fallback are left out of the FDP and recorded per class.
4. **Task 13 exits non-zero by construction** while any check in it fails in both blocks. This
   follows from M7-13a decision 9 and statistician note N3: the merge does not know about
   xfails. So the evidence is read record by record, and the task's exit code is not the
   evidence (rule 4).

## Inputs

| # | Run | Command (from the main checkout at commit C, `GIAP_CALIBRATION_LOG_DIR=$GIAP_ROOT/logs/calibration/m7-gate-<C>`) | Size and R | Expected cost |
|---|---|---|---|---|
| E1 | variants gate, full | `sbatch --array=13 scripts/slurm/calibration_full.sbatch` (runs in-task shards 0/2 and 1/2 side by side, then `--merge`) | G1 `variant_f1` 1M, R = 100; G2, A6a–c and the reference-proteome arm at 100k, R = 100; M7-10's full tests | ≈ 40 min of shards plus ≈ 70–80 min of merge and block 2 on 4 CPUs, inside `--time=04:00:00` (M7-13a decisions 6–7) |
| E2 | novel gate, full, block 1 | ten jobs `sbatch --array=<T> scripts/slurm/calibration_full.sbatch --shard K/10`, K = 0…9, where T is 11 now and 15 after ADR-0048 | `novel_orf_power(…, "1M", n_per_class=NOVEL_PER_CLASS_GATE)`, R = 300; G3, the 24 G4 checks, G4-k3 | ≈ 45–90 min per shard (statistician's shard of 30 replicates: G3 record 2,458 s); 40 CPUs |
| E2b | novel gate, block 2 | only for checks that are `pending` after E2: `--shard K/10 --block 2` | R = 300, replicates 300–599 | as E2 |
| E2m | novel merge | `sbatch --array=<T> scripts/slurm/calibration_full.sbatch --merge`, after E2 (and E2b) | none | minutes, because no block 2 is left to run |
| E3 | fast tier, timings | the fast gate tests of both halves, through `scripts/slurm/test_branch.sbatch` on C | `variants` 20, `novel_orfs` gate 10 | for `benchmarks/baseline.json` |
| E4 | TE2 proxy | the M7-14 weekly Nissle run on C, or a hand-submitted one if no weekly run has used a C at or after the ADR-0049 flag merge | Nissle same-strain peptides | as Slurm 1871578 |

Before each submission, add up the CPUs of running and pending `giap-*` jobs (R14: ≤ 128 CPUs and
≤ 20 jobs). E2 alone takes 40 CPUs and 10 jobs. Record every job ID in PROGRESS.md when it is
submitted, not after it finishes.

The inputs read back are:
- the harness records `<suite>.jsonl` and `<suite>.merged.jsonl` (fields `suite, stream, key,
  criterion, level, observed, bound, se, R, block, status, guard, vacuous, start, stop, seed`);
- the diagnostics `variants_gate_diagnostics.jsonl`, `novel_orfs_gate_diagnostics.jsonl` and
  `novel_orfs_gate_replicates.jsonl`;
- the merge's `summary_markdown` from the job's `.out` file.

G2 and A6a do not go through the harness (M7-13a decision 7). Their evidence is the pytest
outcome in the task's `.out` file together with their diagnostics record.

## Outputs

### O1. The PROGRESS.md gate record

A block headed `### Gate M7 (M7-13c, <date>)` under PROGRESS.md › Done (with evidence). It has:

1. **One status line**, following rule 5: `Gate M7: met` / `partially met` / `not met`, with the
   main commit C and the job IDs E1–E4.
2. **A table with one row per check:** G1, G2, G3, G4 (24 rows: regime × class × q), G4-k3 (one
   row per class × q), A6a, A6b, A6c, attacks a, b, d and GC, and TE2. The columns are: check ·
   stream/key · criterion · R · observed (SE) · bound · block 1 / block 2 · status · source (job
   ID). A status is one of `pass`, `unconfirmed`, `vacuous`, `fail`,
   `known failure (ADR-NNNN)` and `teeth confirmed` (mutant-only attacks a and d).
3. **The vacuous clauses, each stated:**
   - every `vacuous` G4 row, with the class q floor sentence of M7-13 rule 6, word for word: "a
     class with fewer than about 1/α = 100 candidates cannot be called at α = 0.01";
   - G2's `stop_codon` exception, which stays unexercised while M7-16 is undecided (M7-index
     decision 5);
   - G2's `start_codon` exception, which is not exercised at 100k (M7-13 rule 10; M7-14
     exercises it with 80 sites).
4. **Costs:** E1's shard and merge wall times, including the block-2 rerun of G2 and A6a (M7-13a
   decision 6); the wall time of each E2 shard; the peak memory per task.
5. **Open items:** any `known failure`, a strict xfail that is still present, or an `unconfirmed`
   row, each with its ADR or backlog line. If the record leaves something open, the matching
   PROGRESS › Questions for Khoa row is added or updated.

### O2. `docs/methods/variants.md` — a new section `## Variant calls: the M7 gate (M7-13)`

The section has a `### Calibration` log and a `### Sign-off` line. The log gives:

- **The G1 estimand under ADR-0049:**
  - the identifiable-truth rule, with the function name;
  - the share of recall-eligible truth that the rule removes (about 26 % at 100k and 24 % at 1M,
    re-measured);
  - the share of `high` substitution calls that are flagged and therefore not scored (≈ 37 %).
- **G1 results:**
  - pooled P, Rc and F1 per block, with the jackknife SE;
  - R, size, seeds stream `variants.g1`, job ID and commit;
  - the recorded scopes (a)–(c) of the ADR-0049 evidence table, measured again on C.
- **The effective threshold** (statistician, #255 comment 6365982711). Under "fail only if both
  blocks fail" the pass rule has no α. Include the statistician's P(pass) table at the fast SE.
  At the full SE (≈ 0.005), a true F1 of 0.89 passes with P ≈ 0.045.
- **The indel mask:**
  - the masked prediction count;
  - the number of truth rows that fall inside the mask, i.e. the recall asymmetry;
  - the mask is c−3…c+5.
- **Rows outside G1:** the precision of the `high` block, swap and frameshift rows per type
  (ADR-0049 decision 6), recorded and not gated.
- **G2:**
  - 0 ALT rows on `same_strain_insilico`, and the teeth of the mutant on the A6 reads;
  - the rule-of-three bound 3/R on the per-replicate violation rate (statistician N1).
- **A6:**
  - A6a: the unflagged isobaric rows out of all isobaric rows;
  - A6b: the flagged share of random-swap calls;
  - A6c: the share of replicates with a `high` swap, plus the high tokens at the 30 drawn sites
    (statistician note 3);
  - shares carry Wilson 95 % intervals.
- **The reference-proteome arm:** the pooled mean rate per class, against the ≥ 0.95 degeneracy
  assertion.
- **The residual errors:** what the false positives are made of (unannotated start-codon M,
  reversed-stretch decompositions) and what the false negatives are made of (`candidate` tier).

### O3. `docs/methods/novel-orfs.md` — a new section `### Calibration (M7-13 gate)`

The section goes before `### Limitations`. It gives:

- **G3:**
  - pooled F1, P and Rc per block, and F1 per class;
  - start accuracy and the share of planted RBS starts that pass the threshold;
  - the realised planted counts per class, including the intergenic shortfall from
    `truth/novel_orfs.json`;
  - the F1 of the `denovo_subgroup` and `unestimated` arms.
- **G4:** a table of regime × class × q with mean FDP, SE, bound, block and status.
  - The vacuous cells come with the class q floor sentence (rule 6).
  - The fallback calls per class are recorded (M7-13b decision 3).
  - Every `unestimated` reason found is listed. Any reason in `UNESTIMATED_BUGS` is a failure,
    not a note.
- **G4-k3:** one row per class × q. If a row fails, name the method ADR on the per-ORF minimum q
  (M7-13 risk 3; this note's "ORF q = minimum peptide q" section).
- **The clustered-null arm:**
  - attack a's and attack d's mutant FDP against the bound at q 0.10 (the teeth);
  - the filtered G4 on the arm, recorded only;
  - M, and the stream `novel_orfs.cluster`.
- **Attack b's recall** at RBS quantiles 0.10 and 0.50, and the GC-poor attack's intergenic FDP.
  These are fast-tier numbers and are labelled as such.
- **A sign-off line.**

### O4. `benchmarks/baseline.json`

Add entries for the fast tier of each gate half (E3), with wall time, CPUs, node and commit (§8.4).
Add the full-tier wall times too, labelled as recorded and not budgeted.

### O5. Backlog

Close or re-point the backlog lines that M7-13's last acceptance item names (360 A6 part, 614,
665–673, 681–688, 778 first half), and the M7-index line 225 checkbox ("closes when 13 and 14 are
green"). Leave M7-13a decision 7's open backlog items open, but point them at this record. Those
items are the Wilson CI, the merge-side known-failures list, NFS appends and the indel mask across
the origin.

## Rules

1. **Main only.** The gate evidence comes from one full run of E1 and E2 on one main commit C,
   which contains every dependency. This follows M4 decision 7 and M5-11 Q2. Branch-snapshot runs
   (the statistician's `local/stat/m713b-full`, Slurm 1874972) are cited as the branch run, as
   in the M5-07 note's "Branch run / Main run" pairing. They never stand in for the main run.
2. **A fresh log directory per run.** Each run writes to `$GIAP_ROOT/logs/calibration/m7-gate-<C>`.
   Records from another commit, from the fast tier or from a partial attempt are never merged
   into it. A failed attempt keeps its own directory and its job ID in PROGRESS.md.
3. **Evidence is a full-R record.** A record counts only when `R` equals the full count (variants
   100, novel gate 300) and its replicate ranges cover `0…R−1`, plus `R…2R−1` when block 2 ran.
   A check also counts only when its status is not `shard`, `pending` or `skipped`.
4. **Read the records, not the exit code.** E1's exit 1 is expected when the only confirmed fails
   are listed known failures (A6c). Every other `fail` record is a gate failure, whatever the exit
   code. A `MERGE_ERROR` (exit 2) means there is no evidence: rerun, and record the attempt.
5. **The status line.**
   - **met:** every gated check (G1, G2, G3, every non-vacuous G4 cell, G4-k3, A6a, A6b, TE2) is
     `pass` or `unconfirmed`. Attacks a and d have `teeth confirmed`. Attacks b and GC pass.
     Every `vacuous` cell is one that rule 6 of M7-13 predicts (intergenic at q = 0.01, or a de
     novo class with a mean list size below 1/q). The only failure is a known failure whose
     acceptance by Khoa is written down (Q2).
   - **partially met:** each of the above holds except for checks that a pending ADR names.
   - **not met:** anything else.
   - The record never says "met" while a strict xfail other than a Khoa-accepted known failure
     remains in the gate tests.
6. **Determinism.** Seeds are fixed. If C's `src/giap/variants`, `src/giap/novel`, `tests/sim`
   and gate tests match the code of a branch run file for file, the main numbers must equal the
   branch numbers exactly. A difference is a bug, found before the record is written.
7. **No tolerance change (R3).** If a failure turns up here, the record shows it and an ADR goes
   to the lead. `tolerances.yaml` and the gate tests are not edited in this issue.
8. **Roles.** The statistician submits E1, E2 and E4, reads the records, writes O2 and O3 and
   signs them off. The lead writes O1 and O5. Whoever runs E3 writes O4. Docs only, so no
   test-author/implementer loop unless Q1 is accepted.

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | E1 merged G1 record `{stream: variants.g1, criterion: pooled_f1, level: 0.90, observed: 0.933, se: 0.004, R: 100, block: 1, status: pass}` | O1 row `G1 · variants.g1 · pooled F1 ≥ 0.90 · 100 · 0.933 (0.004) · 0.90 · 0.933 / — · pass · <E1 job>`; O2 states P and Rc and the unflagged/identifiable scope |
| 2 | G1 block 1 0.897 (fail), block 2 0.905 (pass) | status `unconfirmed`, counts as passed; the row shows both blocks; O2 cites the effective-threshold paragraph |
| 3 | G4 `picked_decoy`, intergenic, q 0.01: block-1 mean list targets 62 < 1/q = 100 | row `vacuous`, observed = mean FDP, bound —; the floor sentence appears in O1 and O3; the status line is not lowered (rule 6 predicts it) |
| 4 | G4 `picked_decoy`, antisense, q 0.05: mean list targets 14 < 20 | `vacuous` but not predicted (antisense was planted at 300) → status line `partially met` at best, plus a lead question: a broken scenario, not a property of the method |
| 5 | E1 exits 1; merged records: A6c `fail` (share 1.0 in both blocks), every other check passes | the exit code is noted as expected (rule 4); A6c row `known failure (ADR-0049 decision 5)`; the status line is `met` only if Q2's default holds |
| 6 | E1 exits 1; merged records: A6c `fail` **and** A6b `fail` (flagged share 0.97 < 1) | A6b is a gate failure, not hidden behind A6c → `not met`; a question goes to the lead with the unflagged swap rows |
| 7 | E2 merge: G4-k3 cds_outframe q 0.01 FDP 0.070 / 0.068 against bound 0.011 | `fail` → O3 names the per-ORF minimum q ADR; status line `partially met` once the ADR is proposed (M7-13b decision 4), `not met` before that |
| 8 | the log directory holds a G3 record with `R: 10` (fast) next to the full one | the fast record is ignored (rule 3); if it is the only one, there is no G3 evidence |
| 9 | ties: a G4 cell with mean FDP 0.0117 equal to its bound 0.0117; G3 F1 exactly 0.900 | both `pass` (≤ and ≥ are inclusive, §8.4 and `f1_min`) |
| 10 | E4: the latest Nissle TE2 run (1871578) predates the ADR-0049 flag merge | rerun on C; the O1 TE2 row cites the new job; 0 violations expected, since flags never change the tier (ADR-0049 decision 1) |

## Edge cases

| Situation | Behaviour |
|---|---|
| **circular replicon, origin-crossing keys** | none are planted (the M7-12 pool excludes them); O3 says so. The indel mask across the origin stays a backlog line (code review N5 on #255), cited in O2 |
| **linear contig** | the giapsim gate genomes are circular, so the open-3′ `stop ± 3` rule (M7-11 decision 13) is not exercised; O3 lists this under Limitations as not calibrated |
| **minus strand** | both strands are planted (antisense by construction); coordinates are lowest-base keys in every record. The evidence is not split by strand (Q4) |
| **empty arm** | an empty A6a or A6c arm, or a G1 guard miss, raises `RuntimeError` or gives a merged `fail` (M7-13a decisions 5, 7–8): never evidence, never shown as a known failure |
| **empty or missing shard** | a range gap in E2's records → merge error → no evidence; rerun that shard into the same directory |
| **G2 with 0 ALT rows** | `pass`, labelled trivial on this scenario, with the mutant's 45–107 violations per replicate as its teeth (M7-13a decision 4) |
| **a strict xfail XPASSes** | pytest fails the task; it means a fix landed. Remove the xfail in a test PR before the gate run, never in this issue |
| **ADR-0048 decided mid-run** | the run already submitted under task 11 stays valid: R is 300 both ways, and seeds come from `suite_seed(<stream>, i)` (keyed by the stream name, e.g. `novel_orfs.g3`). This holds as long as the ADR-0048 PR keeps the stream names; if it renames them, the seeds change and the run is redone. The record names the task number used |
| **ADR-0049 tests not yet merged** | G1, A6a and A6b would still be strict xfails against the old scope → do not start E1 (§ Dependencies) |
| **`unconfirmed` everywhere on one class** | counted as passing (ADR-0015); O1 lists the rows, so the pattern shows |
| **a timeout (`--time` 04:00:00)** | Slurm `TIMEOUT`, so no merged records; record the job ID and wall time, then rerun E1 with `--shard 0/2` / `1/2` as separate jobs and a separate `--merge` |
| **two shards on different nodes appending to one jsonl** | allowed: appends within one client are atomic, and a torn line fails closed at merge (statistician N4). E2 runs one job per shard, so each shard writes from its own node; a merge error from a torn line means a rerun |

## Acceptance criteria

- [ ] E1, E2 (+ E2b if needed) and E2m ran on one main commit C. Every job ID is in PROGRESS.md,
  and the merged records are under `$GIAP_ROOT/logs/calibration/m7-gate-<C>`.
- [ ] Every gated check in § O1.2 has a full-R record or a pytest outcome on C (rule 3). Q1's
  ledger test, if accepted, checks this from the logs.
- [ ] PROGRESS.md › Done holds the `Gate M7` block (O1). It has the status line under rule 5, the
  per-check table, the floor sentence word for word, the `stop_codon` and `start_codon` lines,
  the costs, and the open items.
- [ ] `docs/methods/variants.md` and `docs/methods/novel-orfs.md` carry the calibration sections
  O2 and O3, with numbers that match the merged records to the printed precision. Each has a
  statistician sign-off line that cites the jobs.
- [ ] `benchmarks/baseline.json` has the fast-tier entries (O4).
- [ ] TE2 has a passing run on C or later (E4).
- [ ] The §14.9 checklist for M7 is in the gate block. It covers:
  - coverage of `giap.variants` and `giap.novel` ≥ 90 %, from the latest main CI run;
  - the methods notes, signed off;
  - the schema, validated;
  - the benchmark;
  - the reviews.
  Any unchecked item is named.
- [ ] The backlog lines in O5 are closed or re-pointed.
- [ ] No change to `tolerances.yaml`, `src/` or existing tests in this PR. A docs-only PR is
  allowed under ≈ 400 lines.

## Estimate

- Docs: O1 ≈ 50 lines, O2 ≈ 50, O3 ≈ 45, O5 ≈ 5.
- `benchmarks/baseline.json` ≈ 20 lines.
- With Q1: `tests/integration/test_m7_gate.py` ≈ 150 test lines (test-author), modelled on
  `test_m5_gate.py`.
- Slurm:
  - E1 ≈ 8 CPU-h allocated (2 shards × 4 CPUs plus the merge);
  - E2 ≈ 10 × 4 CPUs × 1.5 h ≈ 60 CPU-h, plus E2b if block 2 is needed;
  - E3 and E4 are small.
- Wall time ≈ 1 day, most of it waiting on E2.

## Dependencies

All of these are merged to main before E1 and E2 are submitted:

1. **#255 M7-13a** (variants gate harness; queued).
2. **M7-13b** (`m7/calibration-gate-novel`), after the statistician's full tier on the branch
   rules on G4-k3 (decision 4). That full tier is the branch run of rule 1.
3. **#237 M7-11** (pure `giap.novel.fdr`, `docs/methods/novel-orfs.md`) and **#267 M7-11b**
   (step wiring and the giapsim makedb-ID fix, replacing the closed #256).
4. **The ADR-0049 PRs:**
   - the ADR set to accepted-by-khoa (lead queue: `docs/adr-0049-decision`);
   - the flag change (`m7/variant-flags-adr49`, 67421421);
   - the test-side change of G1's scope and of A6a's and A6b's criteria, with their strict xfails
     removed (Q3).
5. **ADR-0048 is not required.** Its outcome only renames the suite (§ Edge cases).

## Questions for the lead

1. **A gate ledger test.** M3, M4 and M5 each have `tests/integration/test_mN_gate.py`, which reads
   the logs and fails on a missing, short-R or failed record. M7-13's split says only "≈ 60 lines
   of docs". **Proposed default:** add `tests/integration/test_m7_gate.py` (test-author, ≈ 150
   lines) on the M5 model:
   - a list of required streams and keys;
   - the records at full R;
   - a `KNOWN_FAILURES` map from stream to ADR, holding A6c;
   - a `GIAP_M7_GATE=1` switch;
   - G2 and A6a read from their diagnostics.
   It goes in its own small PR before the record, so the acceptance criteria can be checked by a
   test.
2. **A6c under ADR-0049.** Khoa's decision 5 says "A6c: unchanged (expected fail)", while §8.5
   says "A failing attack blocks the merge". **Proposed default:** A6c stays a strict xfail, now
   naming ADR-0049 decision 5 as an *accepted* known failure. The status line can then read
   `met`, with A6c listed under open items. O2 also records the flagged share of systematic-swap
   calls (expected 1.0, as decision 3 flags every swap). The lead tells Khoa in the next batch,
   reversibly. The alternative is to restate A6c as "every systematic-swap call is flagged", like
   A6b; that test would pass trivially.
3. **Who carries the ADR-0049 test changes?** The new G1 scope, A6a with the left 2→1 block form,
   A6b as "every swap flagged", and the removal of the xfails all edit #255's tests. **Proposed
   default:** a separate test-author PR, "M7-13a2: ADR-0049 gate criteria", after #255 and the
   flag PR merge. It records the recomputed G1 from ADR-0049 decision 4 in the ADR. M7-13c starts
   only after it merges.
4. **Strand and size strata.** No check is split by strand. **Proposed default:** record nothing
   new. The keys are strand-aware, and both strands are planted. A per-strand recall row costs
   new diagnostics code, which this issue does not write (backlog, if wanted).
5. **Main run against reusing the branch run.** A full main rerun of E2 costs ≈ 60 CPU-h. It would
   repeat the statistician's branch run if the sources match. **Proposed default:** rerun on main
   (rule 1, the M4/M5 precedent), and use the branch run as the determinism cross-check (rule 6).
   If the cluster is busy, the lead may accept a file-for-file source match plus the branch
   records as an interim `partially met`.
6. **The statistician's methods text for variants.md** (M7-13a decision 7: "goes into M7-13c";
   #255 comment 6365982711: "the proposed text is in the lead's report"). I could not find it
   in the repository or under local/handoff. **Proposed default:** the lead puts it in
   `local/handoff/m713c_methods_draft.md`. If it is lost, the statistician writes O2 from the
   review comment's points 1–3 and notes N1–N2.
7. **The fast tier of attacks b and GC as gate evidence.** These attacks run in the fast tier only
   (M7-13b decision 1), so their R is 10, not 300. **Proposed default:** they count for rule 5,
   labelled "fast tier, R = 10". This follows that decision, and §8.5 asks for a passing attack,
   not a particular R.

## Lead decisions (2026-10-23)
1. **Q1 yes:** the test-author adds `tests/integration/test_m7_gate.py` (known-failures map naming
   the ADR for A6c; `GIAP_M7_GATE=1` switch) in a small PR before the evidence PR.
2. **Q2 yes:** A6c stays a strict xfail, recorded as a known failure accepted by Khoa (ADR-0049
   decision 5); the gate can read `met` with A6c listed as open. Record the flagged share of
   systematic swaps. Pending Khoa (reversible): confirm that this reading of decision 5 against §8.5
   is intended; asked in the next batch.
3. **Q3 yes:** separate test-author PR "M7-13a2: ADR-0049 gate criteria" after #255 and the flag PR
   merge; M7-13c waits for it.
4. **Q4 yes:** no strand split.
5. **Q5 yes:** rerun on main; the branch run is the determinism check. Interim `partially met` on a
   file-for-file source match only if the cluster is busy.
6. **Q6:** the statistician's variants.md text (#255 review) is saved at
   `local/handoff/m713c_methods_draft.md` (lead scratch, not committed). It describes G1 before
   ADR-0049; update the G1 scope sentence to the ADR-0049 wording.
7. **Q7 yes:** fast-tier-only attacks (b, GC) count, labelled fast tier.
8. **ADR-0048 seeds:** if accepting ADR-0048 renames the streams, redo E2.
9. **Ledger test (0a3a25e4): gaps 1–7 accepted as proposed.** A6b moves to a replicates check in
   M7-13a2 if that PR moves it out of the harness; the novel suite name comes from
   `GIAP_M7_NOVEL_SUITE` (default `novel_orfs`); the ledger test leaves the `GIAP_M7_GATE` switch in
   a follow-up once #255 and M7-13b merge; vacuous G4-k3 cells count as failures; attacks a/b/d/GC,
   TE2 and the reference-proteome arm stay in the PROGRESS record; G4-k3 gets a known-failure entry
   only with an ADR. Backlog: a commit field in calibration records.
