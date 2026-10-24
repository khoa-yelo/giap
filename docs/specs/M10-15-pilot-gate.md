# M10-15: the M10 gate: final `REPORT.md`, the review packet, Khoa's review

Milestone M10 · **pass 2, the gate** · depends on: M10-14 (freeze commit F, the freeze runs, the
freeze metrics directory, `reproduce.md` at `status=final`), and through it M10-05 to M10-12 ·
estimate ≈ 155 code lines (index: 60; see "Pieces") + ≈ 260 test lines, as two PRs (15a tooling,
15b evidence; question L1) · not `stats` (no new statistic; it renders and checks numbers already
computed) · roles: validation (`metrics.py` hooks, the final template, `report_notes.md`,
`review_packet.py`, the committed report), test-author (`tests/unit/test_m10_15_*.py`), lead
(PROGRESS gate block, the packet hand-off, the PR's `needs-khoa` label, transcribing Khoa's answer) ·
**Khoa: yes, §14.8 item 4.** The gate is a human gate: no agent, the lead included, marks it done.

Written 2026-10-24 against main 4b48786b and the M10 spec branches origin/m10/spec-05 (62d07aa5),
spec-06 (871fce56), spec-07, spec-08, spec-09, spec-12 and spec-14 (071d3947).

## Spec sections implemented

- §13 row M10 gate: "the pilot passes (§15.6) and **Khoa has reviewed** its report".
- §14.8 item 4, "The pilot report (M10)", one of "the only things the agent waits for". R10 says
  everything else the agent decides; this one it does not.
- §15.2 step 8: "then write `validation/pilot/REPORT.md` (§15.5–15.7)".
- §15.5 metrics 1–8, §15.6 (all five bullets, including "A failure is investigated and explained in
  the report before Khoa's review"), §15.7 hypotheses 1–4 ("reported where the pilot data allow").
- §14.9 definition of done, recorded for the milestone in the PROGRESS gate block.
- R3 (no threshold is moved to make a bullet pass), R4 (every claim carries an error estimate or
  `unestimated`), R8, R11 and Khoa's data rule of 2026-10-24 (M10-05 Q-K1, M10-14 lead decision on
  L2/K1): aggregates (counts, shares, CIs, ρ, Jaccards, runtimes, named example genes) and hashes
  may be committed in `REPORT.md` and `metrics.json`; **row-level DIA-NN content never** (no
  peptide or precursor sequence, no per-row table, no column list of a DIA-NN report).
- Index rows 15, the §15.6 evidence table and "Questions that need Khoa" item 1.

## What exists (on the M10 branches; none of it merged yet)

| Piece | Where | What M10-15 uses |
|---|---|---|
| `metrics.py compute/render`, `report_template.md`, `checklist()`, the notes rules (no 4+ digit number absent from `metrics.json`, no 7+ amino-acid string), `**missing: <path>**` and exit 1 for an unresolved placeholder | m10/spec-05 | the final render; the checklist rows become the report's status table |
| `te_real.py`, `te7.py`, the `--te-real/--te7` hook, `_bullet5` (`pass` / `reported (k not passed, explained)` / `pending (M10-06)`), render exit 1 for a non-pass TE without a `## Diagnosis: TE<n>` section naming a cause class | m10/spec-06 | bullet 5 and the TE section; "investigated and explained" is already enforced |
| `entrapment.py`, the `--entrapment` hook (`m8`), H2 shown `deferred` beside `fdp_combined` | m10/spec-07 | metric 8 |
| `a3_feasibility.json` (committed verdict `pass`/`fail`, `reason`) | m10/spec-08 | which A3 branch the report describes |
| `benchmark_standard.sbatch`, `benchmark_summary.py --write` → `benchmarks/baseline.json › standard_run.wheat.<cell>` (4 cells), ADR-M10-D, `docs/user-guide/runtimes.md` | m10/spec-12, spec-12a | bullet 4. M10-12 says bullet 4 "becomes `pass` when the 4 entries exist and ADR-M10-D is accepted (M10-15 reads it)"; **no branch implements that yet**: `checklist()` still returns `pending (M10-12)` |
| `reproduce.md` (header `freeze_commit`, `a3_branch`; Evidence table `step, label, job, state, exit, elapsed, maxrss, node, commit`), `reproduce_check.py parse`, `verify-<F7>.json`, freeze metrics at `$GIAP_ROOT/pilot/metrics/freeze-<F7>/` | m10/spec-14 | the source of every job ID and of F |
| `tests/unit/test_m10_05_guard.py` | m10/spec-05 | allows `metrics.json`, `te1.json`, `REPORT.md` under `validation/pilot/` (no peptide-like string); forbids every other name under `validation/pilot/metrics/` |
| M10-14's R11 audit rule 5 | m10/spec-14 | `validation/pilot/metrics/` holds only `a3_feasibility.json` "and whatever Khoa's Q-K1 answer later allows, added to the allow-list by the lead" |
| Review packets | `local/review/m9-*` (gitignored), M9-25's `local/review/m9-approval/` | the precedent: a gitignored directory, the path and a summary posted to Khoa, the PR labelled `needs-khoa`, merge only after their decision |
| Metric 6 and H4 | M10-10 (**not specced**) | gap G1 |

## Deliverables

| # | File | PR | Role | What |
|---|---|---|---|---|
| D1 | `validation/pilot/metrics.py` hooks | 15a | validation | bullet 4 from `baseline.json` + ADR status; bullet 1's A3 label; an `evidence` key copied from `reproduce.md`; one `evidence` string per checklist row; a `final` render mode (below) |
| D2 | `validation/pilot/report_template.md` (final) | 15a | validation | sections 1–10 below; no numbers |
| D3 | `validation/pilot/review_packet.py` | 15a | validation | builds `local/review/m10-pilot/r<k>/` after its checks pass (stdlib + `metrics.py` imports) |
| D4 | `tests/unit/test_m10_15_*.py` | 15a | test-author | fixtures only; criteria 1–8 |
| D5 | `validation/pilot/REPORT.md`, `validation/pilot/metrics/metrics.json`, `validation/pilot/metrics/te1.json`, `validation/pilot/report_notes.md` | 15b | validation | the rendered report at F and the inputs that reproduce it |
| D6 | `validation/pilot/README.md` › "Pilot report (M10-15)" | 15b | validation | F, the report's sha256, the packet round, a pointer to `REPORT.md` |
| D7 | `PROGRESS.md` › Done: `### Gate M10 (M10-15, <date>)`; › Questions for Khoa: the review row | 15b | lead | the gate record (below) |
| D8 | the edits to existing guard tests that D5 needs (allow-lists) | 15b | **lead** | role guard: the implementer and test-author never edit existing tests (question L7) |

## The final report (`validation/pilot/REPORT.md`)

Rendered by `metrics.py render --final` from `validation/pilot/metrics/metrics.json` and
`validation/pilot/report_notes.md`. It is the freeze metrics directory's `metrics.json`, copied
unchanged (byte for byte), so the committed report can be re-rendered on CI.

**Sections** (the template's order; each number followed by its source, as M10-05):

1. **Status.** The §15.6 table with four columns: bullet · status · evidence · explanation. The
   GIAP commit F, the spec SHA, the report's own statement "Submitted for Khoa's review (§14.8
   item 4)". The report never states that it was reviewed (see "Khoa's review").
2. **Data and arms.** PXD078459 runs, DIA-NN 2.3.1 Academia (ADR-0051), the databases, the A3
   branch: `Cascadia 0.0.7 on the same runs` or `Casanovo 5.2.1 on PXD053727 / UTI89 (fallback,
   M10-08 verdict <reason>)` or `not run (<reason>)`.
3. **Metrics 1–8**, A1/A2/A3 side by side where the metric applies; metric 6 (A2/A3 overlap) and
   metric 8 (A1 free entrapment) included.
4. **TE1 on real data** (M10-05's section, unchanged).
5. **TE2–TE5 and TE7** (M10-06's section) with every `## Diagnosis: TE<n>` from the notes.
6. **Hypotheses 1–4.** H1 and H3 verdicts as M10-05; H2 `deferred (full validation, FDRBench)`
   with A1's genome-inconsistent share and `fdp_combined` side by side and no verdict (M10-07 L4);
   H4 from metric 6 when A3 is the Cascadia branch, `unestimated` otherwise (index L9).
7. **Not evaluated.** One line per metric or hypothesis whose value is `null` with a `reason`
   (for example `m6: different_genome`), saying why the data do not allow it (§15.7 "reported
   where the pilot data allow").
8. **Runtimes.** Metric 7 (pilot measurements) and, separately, the benchmark medians from
   `baseline.json › standard_run` with node class, commit and ADR-M10-D, against §13.2's
   35 min (16–55) and 25 min estimates.
9. **Reproducibility.** F, the `reproduce.md` path and its `verify-<F7>.json` result, TE7's verdict
   at F, the benchmark commit (M10-12 L12: not re-run at F), the clean-up line.
10. **Deviations from the spec.** M10-05's list, plus every deferral the lead made at the freeze
    (M10-14 L3), the A3 branch when it is the fallback or absent, and any accepted ADR that
    changes a pilot setting.

The "Pending in pass 2" section of the draft template is removed. `render --final` rejects a
`metrics.json` in which any checklist status or TE status starts with `pending`.

## Checklist rows: status and evidence (D1)

`checklist()` gains an `evidence` string per row (file › key, job IDs from `metrics.json ›
evidence`, which `compute --runbook validation/pilot/reproduce.md` copies from the runbook's
Evidence table for the rows S8–S14 and S-bench).

| Bullet | Status rule (final) | Evidence string |
|---|---|---|
| 1 GIAP runs end to end on A1, A2 and A3 | `pass` when A1, A2 and A3 exit 0 at F; the label names the A3 branch: `pass (A3 = Cascadia)` or `pass (A3 = Casanovo fallback on UTI89, different strain)`. `fail (A3 not run: <reason>)` when no A3 run exists and the runbook's `a3_branch` is not `pending`. `fail (<arm> exit <n>)` otherwise | the Evidence rows S8, S9, S11 (job, state, exit, commit) |
| 2 TE1 | unchanged from M10-05 (`te1.passed`; ≥ 0.995 inclusive) | `metrics.json › te1.share_exact_inframe`, `te1.passed`, S14 job |
| 3 A1 and A2 compared | unchanged from M10-05 | `metrics.json › m4`, `m5`; S12 job |
| 4 Measured runtimes replace the estimates | `pass` when `baseline.json › standard_run.wheat` holds the 4 cells, each `complete`, **and** ADR-M10-D's `Status:` is `accepted` or `accepted-by-khoa`; `pending (ADR-M10-D <status>)` when the cells exist but the ADR is `proposed`; `fail (incomplete: <cells>)` when a cell is missing or not `complete` | `baseline.json › standard_run.wheat.<cell>.median_s`, ADR path and status, `docs/user-guide/runtimes.md` |
| 5 TE2–TE5 reported, failures explained | M10-06's `_bullet5`; both `pass` and `reported (k not passed, explained)` satisfy the bullet | `metrics.json › te_real.te<n>.status`; for each non-pass TE, the diagnosis heading and cause class |

New `compute` options: `--baseline benchmarks/baseline.json`, `--benchmark-adr docs/adr/ADR-NNNN-*.md`,
`--a3-feasibility validation/pilot/metrics/a3_feasibility.json`, `--runbook validation/pilot/reproduce.md`.
All are read-only, committed inputs; their sha256 go into `metrics.json › inputs` like the others.

## The review packet (`local/review/m10-pilot/r<k>/`, D3)

`python validation/pilot/review_packet.py --metrics-dir $GIAP_ROOT/pilot/metrics/freeze-<F7>
--runbook validation/pilot/reproduce.md --round <k> --out local/review/m10-pilot` (session, one
CPU, under a minute; it copies files and hashes them). Never committed (`/local/` is gitignored).

| Path in the round | Source |
|---|---|
| `index.md` | written by the script: the status table (copied from `REPORT.md` section 1), "what changed since r<k−1>" (the diff of `REPORT.md` against the previous round, when one exists), the questions for Khoa (the open K-rows of PROGRESS that M10 raised), how to answer (below) |
| `REPORT.md`, `metrics.json`, `report_notes.md` | the committed files at the PR head |
| `runs/A1/report.html`, `runs/A2/report.html`, `runs/A3/report.html`, `runs/compare.html` | the freeze outdirs (`$R/out_*/report/`) |
| `reproduce.md`, `verify-<F7>.json` | the runbook and S-check's result |
| `benchmark/summary.md` | `$GIAP_ROOT/bench/<bench_id>/summary.md` (M10-12) |
| `a3_feasibility.json` | the committed verdict |
| `adrs/` | copies of ADR-0051, ADR-M10-B, -C (when it exists), -D and any M10 waiver ADR, with their `Status:` lines listed in `index.md` |
| `MANIFEST.sha256` | `sha256  path` for every file above, sorted by path |

The per-arm `report.html` files hold per-gene content; that is why the packet stays under `local/`.
Khoa reads it on FarmShare; nothing from it is posted to GitHub except the path and the summary.

**Checks before anything is written** (exit 1 on a failed check, 2 on a precondition, nothing written):

1. `render --final` on the committed `metrics.json` and `report_notes.md` reproduces the committed
   `REPORT.md` byte for byte.
2. `metrics.json › giap_commit` = the runbook's `freeze_commit` = every `commit` of the Evidence
   rows S8–S14; the runbook is `status=final`.
3. No `**missing:` string and no status starting with `pending` in `REPORT.md`.
4. Every checklist row has a non-empty `evidence` string.
5. Every file in the table above exists (A3's `report.html` only when bullet 1 reports an A3 run).
6. Round `r<k>` does not exist yet; rounds are never overwritten.

## The gate record (D7) and Khoa's review

**PROGRESS › Done, `### Gate M10 (M10-15, <date>)`**, written by the lead in 15b:

1. **One status line**, from the rules below: `Gate M10: awaiting Khoa's review` /
   `not met: <bullets>` / `met`.
2. **The §15.6 table** (the five rows with status and evidence, as in `REPORT.md`), plus a sixth
   row `Khoa has reviewed the report (§14.8 item 4)` whose status is `awaiting` until their answer.
3. **The report's identity:** the sha256 of `validation/pilot/REPORT.md` (file content, so it
   survives the squash merge), F, the packet path and round.
4. **The §14.9 checklist for M10:** tests per layer (the M10 test files, CI run link), coverage of
   `src/giap` from the latest main CI (validation scripts are outside the 90 % gate and are listed
   with their test files), the methods note `docs/methods/free-entrapment.md` signed off by the
   statistician, schemas validated (`compare_*` tables, `metrics.json` schema string), the
   benchmark (M10-12), the reviews (PR links). An unchecked item is named.
5. **Open items:** every non-pass TE with its cause class; every deferral; every ADR still
   `proposed`; M10-13 (post-gate, not blocking).

**Status rules.**

- `not met: <bullets>`: any of bullets 1–4 is not `pass`, or bullet 5 is neither `pass` nor
  `reported (… explained)`. A failing bullet is never relabelled; a different criterion needs an
  ADR that Khoa accepts (R3, R13).
- `awaiting Khoa's review`: bullets 1–5 are satisfied and the packet round is posted.
- `met`: bullets 1–5 are satisfied **and** an "Answered by Khoa" entry (below) with verdict
  `accept` names the same `REPORT.md` sha256. Only then may the M10 milestone bar show done.

**How the review is requested.** The lead posts in chat to Khoa: the packet path, the report
sha256, a summary of at most 15 lines (the five statuses, the A3 branch, every non-pass TE with
its cause class, the open ADRs), and the question "Verdict on the M10 pilot report: accept,
accept with changes, or reject?". The same path and summary go into the 15b PR description, the
PR gets the `needs-khoa` label, and a PROGRESS › Questions for Khoa row is added ("M10 pilot report
review, round r<k>"). The PR is not merged before the verdict (M9-25's precedent).

**How the review is recorded.**

- The only valid source is **Khoa's own message**: a user turn in a session, or their edit of
  PROGRESS. A message from any agent, a subagent report, a PR comment, or silence is not a review.
  GitHub comments do not count, because the agents post through Khoa's `gh` login and the two
  cannot be told apart.
- The lead transcribes it into PROGRESS › Answered by Khoa as one entry:
  `- <date>: **M10 pilot report** (REPORT.md sha256 <64 hex>, packet r<k>): "<Khoa's words,
  verbatim>" → verdict accept | changes | reject`. The same quote goes into a PR comment.
- "Don't block on Khoa" and the reversible-default practice do **not** apply: a §14.8 checkpoint
  has no default. While no answer exists, the gate line stays `awaiting Khoa's review` and work
  continues elsewhere (M10-13, M11 preparation).
- After `accept`: the lead removes `needs-khoa`, closes the Questions row, sets the status line
  to `met`, and merges. After `changes`: the requested edits are made (prose in
  `report_notes.md` or the template; a method change re-opens the relevant issue and, if code
  changes, a new freeze per M10-14), the report is re-rendered, and a new round `r<k+1>` is posted.
  After `reject`: the gate stays `not met`, and the lead writes the follow-up issues from their
  reasons. After three rounds without `accept`, R9: a Questions-for-Khoa row asks how to proceed.

## Examples

Numbers are hand-made, not measurements.

| # | Input | Expected |
|---|---|---|
| 1 | Freeze at F; A1, A2, A3 (Cascadia) exit 0; TE1 0.997; compare tables valid; 4 benchmark cells `complete`, ADR-M10-D `accepted`; TE2–TE5 `pass` | checklist all `pass`, bullet 1 `pass (A3 = Cascadia)`; `render --final` exit 0; the packet `r1` is written with `MANIFEST.sha256`; the gate line reads `awaiting Khoa's review`, never `met` |
| 2 | As 1, but TE3 `fail` and `report_notes.md` has `## Diagnosis: TE3` with cause class `biology` | bullet 5 `reported (1 not passed, explained)`, satisfied; section 5 shows the diagnosis; the gate block lists TE3 under open items; status `awaiting Khoa's review` |
| 3 | As 2, but the notes have no TE3 diagnosis | `render` exit 1 naming TE3 (M10-06 rule); no `REPORT.md`, no packet; nothing is committed |
| 4 | TE1 `share_exact_inframe` 0.994 (`passed: false`), explained in the notes as a mapping defect (`pipeline_bug`) | bullet 2 `fail`; status `not met: TE1`. The packet is still built (Khoa reviews failures too), and the summary leads with the failure. Fixing the defect means a new freeze; the 0.995 threshold is not touched (R3) |
| 5 | M10-08 verdict `fail` (`low_genome_match`); M10-11b ran Casanovo on UTI89, `giap run` exit 0; compare of A1, A2 only | bullet 1 `pass (A3 = Casanovo fallback on UTI89, different strain)`; metric 6 `null`, `reason: different_genome`, listed in "Not evaluated"; H4 `unestimated`; section 10 names the fallback |
| 6 | M10-08 verdict `fail` (`route_unavailable`) and Khoa declined ADR-M10-C, so no A3 run exists; runbook `a3_branch=fallback`, A3 steps without Evidence | bullet 1 `fail (A3 not run: route_unavailable; ADR-M10-C declined)`; status `not met: bullet 1`; metric 6 and H4 under "Not evaluated"; the lead drafts a waiver ADR ("M10 gate without A3", Needs Khoa: yes) for the same review round |
| 7 | The template has `{{m6.A2_A3.jaccard}}` but `metrics.json` has no `m6` key (M10-10's hook never ran) | `**missing: m6.A2_A3.jaccard**`, `render` exit 1; the packet script refuses (check 1 or 3); this is a defect to fix, not a "Not evaluated" line |
| 8 | `m8.peptide` has `fdp_combined: null`, `assessment: not_evaluable`, `reason: empty_discoveries` | rendered as "not evaluable (empty_discoveries)" in metric 8 and in "Not evaluated"; exit 0 (a key that exists with a null value and a reason is allowed) |
| 9 | 4 cells present, `sim5m_core` `incomplete`; ADR-M10-D `proposed` | bullet 4 `fail (incomplete: sim5m_core)` (an incomplete cell wins over the ADR status); status `not met: bullet 4`; with the cell complete and the ADR still `proposed`: `pending (ADR-M10-D proposed)`, and `render --final` rejects it |
| 10 | Khoa answers in chat: "Report reads fine; reword the H1 paragraph so it does not sound like a test of biology." | entry with verdict `changes`; the H1 prose changes in the template; re-render; new sha256; packet `r2` with a diff of `REPORT.md` against `r1`; the gate line stays `awaiting Khoa's review` |
| 11 | A subagent's hand-off says "Khoa approved the pilot report" | not a review record. The lead writes nothing to Answered by Khoa; the gate line is unchanged |

## Edge cases

| Situation | Behaviour |
|---|---|
| **Circular and linear replicons** | Not computed here. Origin-crossing genes and hits are handled where the metrics are computed (TE1 shifts ±L in M10-05; M10-07 and M10-09 follow the genome's topology tags). The report states "circular-tagged chromosome and plasmids" from the committed genome. Under the fallback, the UTI89 genome must carry its own tags (M10-11a); `metrics.json › inputs` holds its hash |
| **Minus strand** | Not computed here; the strand logic lives in the reused M10-05/06 functions. The report carries no coordinates |
| **Empty input** | An arm with no accepted target peptide reaches the report as M10-05's `fail (no data)` and blocks bullet 2 or 3; it is never rendered as a pass. An empty or absent `report_notes.md` is allowed only when every TE passes. A `metrics.json` without `checklist` is a precondition failure (exit 2) |
| **Ties at a threshold** | TE1 exactly 0.995 passes (inclusive, M10-05). H3's `ci_low` exactly 0.7 holds. A benchmark median exactly at the §13.2 estimate is reported, not judged. Two accepted attempts of one arm at F: the latest passing one (M10-04), and the evidence string names it |
| **A failed TE** | Examples 2 to 4. TE2–TE5: explained, then reported, and the bullet is satisfied. TE1: bullet 2 fails. TE7 is not a §15.6 bullet: a TE7 `fail` (including `known_lapack`, M10-06 L9) is shown in sections 5 and 9 with its diagnosis and listed as an open item; it does not change the status line |
| **A3 skipped by the M10-08 verdict** | Fallback ran: example 5. No A3 at all: example 6. The runbook still at `a3_branch=pending`: bullet 1 stays `partial`, which `render --final` rejects; the report is not built |
| **A missing metric** | A key that is absent is a defect (example 7). A key with `null` and a `reason` is a finding (example 8), listed under "Not evaluated" |
| **A freeze bug found after the packet** | A new F′ (M10-14); every arm re-runs at F′; a new `metrics.json`, a new report sha256, a new round. A verdict on the old sha256 does not carry over |
| **`metrics.json` from another commit** | `giap_commit` ≠ `freeze_commit`: packet check 2 fails (exit 1). The lead re-runs S14 at F |
| **Main moves between packet and verdict** | The report is identified by its sha256 and by F, not by main's HEAD. Merging main into the 15b branch changes no committed report file. A change to `metrics.py` or the template on main is a conflict for the lead to resolve; the CI re-render test (criterion 2) catches a renderer change that alters the report |
| **Khoa answers with a question** (no verdict) | Answered in chat; recorded in PROGRESS as a question; no verdict is inferred |
| **Khoa accepts but asks for later work** (for example a full validation item) | verdict `accept`; the requests become backlog lines linked from the gate block |
| **Row-level content slips into the notes** | `render` rejects 7+ amino-acid strings and unknown 4+ digit numbers (M10-05). Criterion 6 adds the identifier-list check for `metrics.json` |

## Acceptance criteria

On CI (every commit; fixtures only, no `$GIAP_ROOT`):

1. **Checklist hooks.** `checklist()` gives the statuses of the table above on fixtures for:
   bullet 1 with Cascadia, fallback, no A3 (example 6) and `pending`; bullet 4 with 4 complete
   cells and an `accepted` ADR (`pass`), a `proposed` ADR (`pending (ADR-M10-D proposed)`), and a
   missing or incomplete cell (`fail (incomplete: …)`); every row carries a non-empty `evidence`
   string naming a file and key.
2. **The committed report reproduces.** After 15b, `render --final --metrics
   validation/pilot/metrics/metrics.json --notes validation/pilot/report_notes.md` equals the
   committed `validation/pilot/REPORT.md` byte for byte. (In 15a the test skips when the file is
   absent, and runs on a fixture report.)
3. **Final mode.** `render --final` exits 1 on a fixture with any `pending` checklist or TE status,
   with a `**missing:` placeholder, or with a non-pass TE that lacks a diagnosis; the draft mode keeps
   M10-05's behaviour.
4. **Template.** The final template has sections 1–10, no "Pending in pass 2", no number, and
   placeholders for metric 6, metric 8, H2 and H4; the M10-05 guard tests still pass on it.
5. **Packet.** `review_packet.py` on a sandbox (`tmp_path`) writes exactly the files of the packet
   table plus `MANIFEST.sha256`, whose hashes match; it exits 1 and writes nothing when the
   re-render differs, when `giap_commit` ≠ `freeze_commit`, or when the round exists; exit 2 when
   an input is missing; it never writes outside `--out`.
6. **Data rule on the committed files.** `REPORT.md`, `metrics.json`, `te1.json` and
   `report_notes.md` contain no 7+ amino-acid string, and `metrics.json` holds no array of more
   than 20 strings matching a locus tag or accession pattern (`[A-Z][A-Z0-9]+_\d{3,}`, UniProt
   accessions) (question L6). No `te_real.json`, `te7.json`, `entrapment.json`, exceptions file
   or run table is tracked.
7. **Gate record shape** (`test_m10_15_gate_record.py`, on `PROGRESS.md`): when a `### Gate M10`
   block exists, its status line is one of the three forms; a `met` line requires an Answered-by-
   Khoa entry with verdict `accept`, a quoted string, and a 64-hex sha256 equal to the sha256 of
   the committed `REPORT.md`. The test checks the shape only; it cannot check who wrote the quote,
   which is why the recording rule above binds the lead.
8. **Consistency.** `metrics.json › giap_commit` equals `reproduce.md`'s `freeze_commit`, and its
   `inputs` entries for committed files hash to their committed contents.

On FarmShare and in the PRs (evidence):

9. `metrics.py compute` (with `--baseline`, `--benchmark-adr`, `--a3-feasibility`, `--runbook`,
   `--te-real`, `--te7`, `--entrapment`, `--a3` when it applies) and `render --final` exit 0 at F
   into the freeze metrics directory; the files of D5 are copies of that directory's output.
10. `review_packet.py` wrote round `r1` (or later); the path, the report sha256 and the summary are
    in the 15b PR description and in chat; the PR carries `needs-khoa`.
11. PROGRESS holds the gate block (D7) with the status line by the rules, the §15.6 table with
    evidence, the report identity, the §14.9 checklist and the open items; a Questions-for-Khoa row
    for the review exists.
12. **The human gate.** The status line reads `met` only after an Answered-by-Khoa entry
    transcribed from Khoa's own message, with verdict `accept` and the report's sha256. No agent
    ticks this criterion; the lead ticks it only by quoting Khoa.

## Pieces

| Piece | PR | Role | Lines |
|---|---|---|---|
| `metrics.py`: bullet 4 from `baseline.json` and the ADR `Status:` line | 15a | validation | 25 code |
| `metrics.py`: bullet 1 A3 label from `a3_feasibility.json` and the runbook header | 15a | validation | 15 code |
| `metrics.py`: `--runbook` evidence copy (via `reproduce_check.parse`) and per-row `evidence` strings | 15a | validation | 25 code |
| `metrics.py`: `render --final` (pending, missing, diagnosis checks) | 15a | validation | 15 code |
| `review_packet.py` (checks 1–6, copies, manifest, `index.md`, round diff) | 15a | validation | 75 code |
| final `report_template.md` | 15a | validation | ≈ 90 docs |
| `tests/unit/test_m10_15_checklist.py`, `…_render_final.py`, `…_packet.py`, `…_committed.py`, `…_gate_record.py` | 15a (gate record and re-render tests skip until 15b) | test-author | ≈ 260 test |
| `REPORT.md`, `metrics.json`, `te1.json`, `report_notes.md` | 15b | validation | data and prose |
| README section, PROGRESS gate block and Questions row | 15b | validation, lead | ≈ 15 + 40 docs |
| allow-list edits in the existing M10-05 guard and M10-14 R11 audit tests | 15b | lead | ≈ 5 test |
| **Code total** | | | **≈ 155** (index 60: it counted neither the bullet-4 hook nobody else implements, nor a packet builder with checks) |

Wall time: 15a ≈ 1 day of agent time; 15b ≈ 1 hour after the freeze metrics exist, then Khoa's
review (unbounded; the gate waits, other work does not). Compute: one session CPU for minutes;
no Slurm job.

## Gaps found while writing this spec

- **G1. Metric 6 and H4 have no spec.** M10-10 is not specced, so the key names (`m6`, `h4`), the
  overlap definition (A2/A3 novel ORFs and variants by class and confidence) and its `reason`
  values are not fixed. The template's placeholders here are provisional (question L10).
- **G2. Bullet 4's hook is unowned.** M10-12 says "M10-15 reads it", but `checklist()` still
  hard-codes `pending (M10-12)`. This spec takes it (D1).
- **G3. Whether ADR-M10-D needs Khoa** is not settled: M10-12 says only "accepted" (question L2).
- **G4. Recording a human gate.** No earlier gate spec says how Khoa's answer is told apart from an
  agent's text, and GitHub cannot do it here (the shared `gh` login). The rule above is new
  (question L4).
- **G5. Allow-lists.** M10-05's guard allows only `metrics.json`, `te1.json`, `REPORT.md` under
  `validation/pilot/metrics/`; M10-14's R11 audit allows only `a3_feasibility.json` there. The two
  must agree once D5 lands, and only the lead may edit them (question L7).

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| L1 | ≈ 155 code lines against the index's 60. One PR, or 15a (tooling + tests, mergeable before the freeze) and 15b (evidence, `needs-khoa`)? | split; 15b holds no code beyond the lead's allow-list lines |
| L2 | Bullet 4: is ADR-M10-D `accepted` by the lead enough, or must it be `accepted-by-khoa`? It supersedes the §13.2 estimates (a spec number, R13). | `accepted` suffices for the bullet; the ADR is `Needs Khoa: yes` and goes into the same review round, so their verdict covers it |
| L3 | No A3 run at all (example 6): status `not met`, with a waiver ADR in the same round, or `awaiting review` with the gap noted? | `not met: bullet 1` plus a waiver ADR (Needs Khoa: yes); the gate can read `met` only if Khoa accepts both the ADR and the report |
| L4 | Recording the review: only Khoa's own message (user turn or their PROGRESS edit), quoted verbatim with the report's content sha256; GitHub comments excluded? | yes |
| L5 | Commit the notes as `validation/pilot/report_notes.md`, so CI can re-render the report byte for byte? (The guard forbids new names under `metrics/`.) | yes |
| L6 | Identifier lists in the committed `metrics.json` (warnings, TE7 `differing` paths, contig names): ban arrays of more than 20 locus tags or accessions as row-level content? | yes, criterion 6; table paths, column names and contig names stay |
| L7 | The two existing guard tests (M10-05 guard, M10-14 R11 audit) need matching allow-lists: the lead edits them in 15b (role guard)? | yes, the lead, in 15b, with the trailer `Agent-Role: lead` |
| L8 | Merge 15b only after Khoa's verdict (M9-25 precedent), or merge first and record the review afterwards? | after the verdict; the branch is kept current by merging main |
| L9 | Named example genes (Khoa allowed them): add any beyond the ones TE3/TE4 already name? | no new names in the template; at most five in `report_notes.md`, each tied to a TE or metric |
| L10 | Metric 6 / H4 key names before M10-10 is specced (G1): placeholders `m6.<class>.<confidence>` and `h4.verdict`, with M10-10 conforming, or wait for M10-10? | provisional placeholders now; M10-10's spec must match them or update this template in its own PR |
| L11 | `metrics.json` gains `evidence` and the A3 label: keep the schema string `giap-pilot-metrics/1` (additive keys) or bump to `/2`? | keep `/1`; keys are additive and M10-05's guard lists top-level keys, so the lead adds `evidence` to that list in 15b |
| L12 | §14.9 "coverage ≥ 90 %": the validation scripts sit outside `src/giap`. Report their test files only, or add them to a coverage run? | report only; no new coverage gate (R3: gates are not changed here) |

## Questions for Khoa (asked with the packet; nothing here has a default)

| # | Question |
|---|---|
| K1 | The gate itself: verdict on `validation/pilot/REPORT.md` (sha256 in the packet): accept, accept with changes, or reject (§14.8 item 4). |
| K2 | ADR-M10-D (benchmark definition; the measured runtimes replace §13.2's estimates): accept? |
| K3 | Conditional: a waiver ADR when A3 did not run (example 6), or the fallback's "different strain" caveat when it did (example 5): accept the pilot on that basis? |
| K4 | Conditional: every TE with cause class `definition` or `search` (M10-06): the deviation it states, for their decision (no definition or setting is changed by the agents). |

## Lead decisions (2026-10-24)

L1–L12 take the recommended defaults above. On G1, M10-10 is being specced now (m10/spec-10): 15a uses
the provisional `m6`/`h4` keys and renames them to M10-10's final names when that spec lands. On G5,
the lead edits the guard allow-lists in 15b (L7). K1–K4 go to Khoa with the packet; there is no default.

### Lead decisions on the 15a tests (2026-10-24)

- **The M10-05 guard's `TOP_KEYS`:** gains `m6`, `h4` and `evidence` in 15a, as a lead edit. It is not deferred to 15b.
- **Bullet 1:** without `--runbook` it keeps M10-05's plain `fail` / `pass` strings, so the existing M10-05
  test still holds. The labelled forms (`fail (A3 exit 1)`, …) apply only with `--runbook`.
- **`m6` naming:** follows M10-10 (`m6.overlap_a2_a3`, reason `different_strain`). Examples 5 and 7
  are read with those names.
- **A complete benchmark cell:** it is present in `baseline.json › standard_run` with a non-null
  `median_s`. M10-12 writes no `status` key.
- **A missing packet input:** exit 2 (criterion 5). A check that fails exits 1.
- **15a's implementation:** waits for M10-06 (#306) and M10-14a to merge, because it uses `diagnosis_problems`
  and `reproduce_check.parse`.

## Lead decisions after the 15a code review (2026-10-24)
These supersede earlier text that differs. Review: rev-1015.
- **B1, the runbook input:** `--runbook` is an input of `render --final` and of the packet, not of `compute`. `metrics.json` therefore never records the runbook. This avoids the circularity in which the runbook is still a draft inside S14, and the final runbook's Outputs then hashes `metrics.json`. Bullet 1 and the S8–S14 evidence are computed at render time, from the final runbook.
- **B2:** superseded Evidence rows (M10-14) are ignored by `gate_evidence` and by packet check 2.
- **S1:** "A3 not run" is a finding. M10-10 emits `reason: "a3_not_run"` when `--a3` is absent and keeps `missing:<table>` for an absent table. `--final` renders `a3_not_run` under "Not evaluated". M10-10's "byte for byte" sentence is read as "byte for byte apart from the keys M10-15 adds".
- **S2:** the bullet-1 label comes only from the runbook's `a3_branch`. `pending` gives `partial`, and a Cascadia branch with a `fail` feasibility verdict gives `fail (A3 feasibility fail)`.
- **S3:** the packet's `--out` must be outside the checkout, or a path the checkout git-ignores; otherwise it exits 2.
- **S4:** on an error the packet removes its `.r<k>.tmp`, never overwrites a round (it checks this explicitly, with its own message) and returns only the documented exit codes.
- **S5:** check 1's message carries the refusal reason.
- **S6:** checks 3 and 4 are folded into `final_problems` (the evidence check). The comparison of the freeze directory's `metrics.json` with the committed one stays, and gets a test.
- **S7:** an absent `--a3-feasibility` gives the feasibility block a `pending (M10-08 verdict)` reason, which `--final` refuses.
- **S8:** the gate-record checker:
  - rejects "awaiting" or "met" while a checklist bullet fails, unless a Khoa-accepted waiver ADR is cited;
  - accepts a quote containing `"` (escaped) or spanning several lines;
  - lets a later `changes` or `reject` on the same sha256 override an earlier `accept`;
  - checks every `### Gate M10` block.
- **S10:** the 15a PR title does not start with `M10-15`, so the progress bar doesn't count it as done before Khoa's review.
- **Nits:**
  - N1: a `rejected` or `superseded` ADR gives `fail (ADR-M10-D <status>)`.
  - N2: drop "; ADR-M10-C declined" from example 6's label.
  - N4: in draft mode the placeholders render the pending marker.
  - N5: the packet reuses `adr_status`.
  - N6: no `__pycache__` in the repo (`sys.dont_write_bytecode`, or `-B`).
  - N7: check 1 compares bytes.
  - N3 (a mistyped path under a finding reason): documented in the methods/README, no code change.
