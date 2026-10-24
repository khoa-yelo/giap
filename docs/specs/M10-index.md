# M10 "Pilot validation" (E. coli Nissle 1917, PXD078459): issue index

Spec: GIAP_SPEC.md §13 row M10. Deliverables: "§15: raw data download, DIA-NN (A1, A2) and de
novo (A3) searches as Slurm jobs, conversion, runs A1/A2/A3, `giap compare`, pilot report, the
16-thread benchmark". **Gate:** "the pilot passes (§15.6) and **Khoa has reviewed** its report"
(§14.8 checkpoint 4).

Sources read: §15 (whole: §15.1 data and arms, §15.2 steps 1–9 and the DIA-NN template, §15.3,
§15.5 metrics 1–8, §15.6 pass criteria, §15.7 hypotheses), §13 rows M10 and §13.2 (release
criteria: TE1 ≥ 99.5 %, mapping FDR ≤ 1 %, "a standard run on 16 threads and 64 GB takes about
35 min"), §11.6 TE1–TE8, §14.4 speed loop ("after M3, M8 and M10"; `sacct` MaxRSS, same node
class; < 5 % gain over 3 attempts), §3 (§3.1 partitions and QOS, §3.2 source table rows
"iProX / ProteomeXchange" and `regulondb.ccg.unam.mx`, §3.6 `needs_diann`/`needs_gpu`, §3.7, §3.8
job header and polling, §3.9 disk budget), §0.2 R5, R6, R9–R11, R13, R14, §5.4 converters
(DIA-NN 2.x, Casanovo 5.x, Cascadia 0.0.x; `pending_sample`), §7.9.3 `giap compare`, §14.8, §17
risks (Cascadia on timsTOF, acquisition mode, DIA-NN flags, download source down).
ADRs: ADR-0042 (`compare` replicate counts), ADR-0046 (`giap compare` layout: locus-tag join,
genome identity, exit 1 on a different genome), ADR-0009/0047 (picked FDR, decoy q filter, which
the DIA-NN `--report-decoys` output feeds). No accepted ADR covers DIA-NN 2.3.1, de novo search or
the benchmark. Prior specs: M3-24 (the pilot kit and its 17 lead decisions), M3-22 (TE1/TE8
proxies), M3-23 (`giap makedb`), M9-22/M9-23 (`giap compare`, V16), MVP-04 (the K-12 evidence run).

One issue = one PR; size counts **code** only (ADR-0007). Runs on real data are evidence; the
commands, job IDs and numbers go into `PROGRESS.md`, and only small derived files are committed.
The **validation** role owns `validation/` (§14.1); converters live in `src/giap/io/converters/`
(implementer); statistics-bearing scripts get the `stats` label.

## How M10 is sequenced (Khoa's standing direction: MVP first)

- **Pass 1 (MVP: the pilot end to end on A1 and A2), 01–05.** Download and acquisition-mode check
  (01), the two DIA-NN searches (02), conversion with the converter contract-tested on the real
  sample (03), `giap run` for A1 and A2 plus `giap compare out_a1 out_a2` (04), and a first metrics
  pass with TE1 and a draft `REPORT.md` (05). After 05, three of the five §15.6 bullets have
  evidence (A1/A2 end to end, TE1, A1/A2 compared).
- **Pass 2 (A3, the remaining pass criteria, hardening and the gate), 06–15.** Real-data TE2–TE5
  and TE7 (06), the free-entrapment estimator (07), A3 feasibility, converter and run (08–10), the
  conditional Casanovo fallback (11a/11b), the 16-thread benchmark (12), the speed loop (13),
  reproducibility and cleanup (14), and the final report with Khoa's review (15).

A3 is the only arm with a real chance of failing (§17 "Cascadia fails on timsTOF data"), so it
starts in pass 2 but its feasibility issue (08) may run in parallel with pass 1: it shares no code.

## How this design lands

| Commit | Files |
|---|---|
| this one | this index and the full spec of the first issue, `M10-01-pxd078459-download.md` |
| later | per-issue specs M10-02 … M10-15, written one wave ahead of the feature loop |

## Already built (checked on main, d29a2ea5)

| Piece | State |
|---|---|
| `validation/pilot/` (M3-24, #122) | `diann_A1.cfg`, `diann_A2.cfg` (two-step predict-then-search, no `--met-excision`, no MBR, decoys reported), `db_a2/` (4,692 proteins, 989 KB), `README.md` (version of record, flag table, run commands). **Khoa accepted the flag table on 2026-10-20** ("looks fine") |
| DIA-NN | **2.3.1 Academia** Apptainer image `/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif` (532 MB, sha256 `1946452d…6381` pinned in the README and in `tests/unit/test_m3_24_pilot_kit.py`). The spec names 2.7.0; no ADR records the change yet (see ADR-M10-A). **The file is mode `rwxrwxrwx`**: anyone can replace it, so every DIA-NN job checks its sha256 first (02) |
| `scripts/slurm/fetch_pxd078459.sbatch` | written in M3 with **unverified** defaults (`https://download.iprox.cn/dataset/PXD078459/WT_{1,2,3}.d.zip`). A HEAD request on 2026-10-23 returned **404** for that base URL; the ProteomeXchange record (`…/cgi/GetDataset?ID=PXD078459`) answers 301 to `…/ui?pxid=PXD078459` and `www.iprox.cn` answers 200. M10-01 replaces the defaults with the file list of the official record |
| `scripts/slurm/diann_pilot.sbatch`, `denovo_pilot.sbatch` | written and stub-tested; never run on real data. `denovo_pilot.sbatch` exits 2 until `workflow/envs/cascadia.yml` and the weights exist |
| `giap convert --from diann` | implemented (M1-07) from the documentation; `tests/contract/test_diann_contract.py` waits for a real sample. `--from cascadia`/`casanovo`/`instanovo` are CLI choices only: `CONVERTERS` holds `diann` and `generic`, so they exit 5 |
| `giap run`, `giap compare` | on main: `--reference-peptides`, `--genome-relation`, `--profile`, `--executor`, `--threads`, `--seed`; Nextflow writes `OUT/logs/trace.tsv` (metric 7). `giap compare` (M9-22) and V16 / `report/compare.html` (M9-23) are merged |
| Inputs on disk | `$GIAP_DATA_DIR/proteomes/UP000000625.fasta`, `$GIAP_DB_DIR/contaminants/contaminants.fasta` (201 KB), Nissle genome (`tests/data/genomes/nissle1917.fna.gz`, circular-tagged) and its recorded full-DB Bakta JSON (`tests/data/recorded/bakta/1.12.2-full/nissle/`), Bakta full DB, geNOMAD, Promotech, antiSMASH |
| TE proxies | `validation/proxies/te_proxies.py` (TE1, TE2, TE8 in silico); real-data TE1–TE5 do not exist |
| Missing | `validation/pilot/entrapment.py`, `scripts/slurm/benchmark_standard.sbatch`, `workflow/envs/cascadia.yml`, `casanovo.yml`, any `.d → mzML` route, `tests/data/recorded/diann/`, `…/cascadia/`, `validation/pilot/inputs/`, `REPORT.md` |
| `local/` | 100 GB used of the 400 GB budget (§3.9); class share 197 TB free |

## Issues

| Order | Issue | Pass | Scope | Depends on | Code lines | Stats | Khoa |
|---|---|---|---|---|---|---|---|
| 1 | [M10-01](M10-01-pxd078459-download.md) | 1 | PXD078459: file list from the official record, download job with checksums, acquisition-mode check of all three runs | none | 250 | no | only if the data are login-gated or the selection is ambiguous |
| 2 | M10-02 | 1 | DIA-NN A1 and A2 as Slurm jobs; image sha256 check; version, settings and runtimes recorded; 1,000-row report samples under `tests/data/recorded/diann/2.3.1/`; ADR-M10-A | 01 | 120 | no | ADR-M10-A acceptance; re-review if the acquisition mode forces a settings change |
| 3 | M10-03 | 1 | `giap convert --from diann` on both reports; contract test on the real sample (test-author); converter fixes the sample exposes; `validation/pilot/inputs/a{1,2}.tsv.gz` (+ `abundance` from `PG.MaxLFQ`) | 02 | 150 | no | no |
| 4 | [M10-04](M10-04-pilot-runs.md) | 1 | **MVP closer:** `giap run` A1 and A2 (`giap_run.sbatch`, `farmshare,apptainer`), `giap compare out_a1 out_a2 --labels A1,A2`; run directories, job IDs, exit codes, `report.html` per arm | 03; MVP-04 acceptance run (soft) | 80 | no | no |
| 5 | M10-05 | 1 | `validation/pilot/metrics.py`: §15.5 metrics 1–5 and 7 for A1/A2, **TE1 on real data**; `REPORT.md` draft with the §15.6 checklist | 04 | 360 | light | no |
| 6 | M10-06 | 2 | real-data TE2–TE5 (`validation/pilot/te_real.py`) and TE7 (a second A2 run, byte comparison); failures investigated in the report | 05 | 280 | no | no |
| 7 | M10-07 | 2 | §15.5 metric 8: `validation/pilot/entrapment.py` (Wen 2025 combined estimator, in-silico r under the searchable length, charge and m/z range), methods note | 05 | 220 | **yes** | no |
| 8 | M10-08 | 2 | A3 feasibility: `workflow/envs/cascadia.yml`, weights, the `.d → mzML` route (ADR-M10-B), one-run GPU smoke, the "fails on timsTOF" decision rule | 01 | 150 | no | if the conversion route needs a vendor-licence acceptance |
| 9 | M10-09 | 2 | `giap convert --from cascadia` (`src/giap/io/converters/cascadia.py`), contract test on a 1,000-row sample under `tests/data/recorded/cascadia/0.0.7/` | 08 | 200 | no | no |
| 10 | M10-10 | 2 | A3 run: 3 runs on the `gpu` partition, convert, `giap run` A3 with `--reference-peptides a2`, `giap compare out_a1 out_a2 out_a3`; metric 6 (A2/A3 overlap) and hypothesis 4 | 09, 05 | 180 | light | no |
| 11a | M10-11a | 2 (conditional) | fallback inputs: ADR-M10-C (PRIDE or the record's host as a source), PXD053727 3 DDA runs, UTI89 genome GCA_000013265.1 and its full-DB Bakta recording, `.raw → mzML`, `workflow/envs/casanovo.yml` and weights | 08 says "fails" | 200 | no | **yes**: ADR-M10-C adds a download source (§14.8 item 1) |
| 11b | M10-11b | 2 (conditional) | `giap convert --from casanovo`, contract test, Casanovo 5.2.1 run, `giap run` on UTI89; metric 6 with the "different strain" caveat | 11a | 260 | light | no |
| 12 | M10-12 | 2 | `scripts/slurm/benchmark_standard.sbatch` (16 CPU, 64 GB, one node class, median of 3; standard and core-only), measured runtimes replace the estimates (ADR-M10-D, `benchmarks/baseline.json`, user guide, PROGRESS) | 04 | 160 | no | no |
| 13 | M10-13 | 2 (post-gate) | speed loop after M10 (§14.4; performance role) on the 5-Mb simulated genome and K-12, stop rule < 5 % over 3 attempts; **does not block the gate** | 12 | per attempt ≤ 400 | no | no |
| 14 | M10-14 | 2 | reproducibility: `validation/pilot/reproduce.md` runbook (every command, input hash, job class), the report-freeze re-runs of A1/A2/A3 on the final main, `nextflow clean`, deletion of libraries and raw reports per §3.9, R11 audit of committed files | 06, 07, 10 (or 11b), 12 | 100 | no | no |
| 15 | M10-15 | 2 | **gate:** final `REPORT.md` (§15.5–15.7, every §15.6 bullet with evidence, failures explained), the review packet under `local/review/m10-pilot/`, Khoa's review | 14 | 60 | no | **yes** (§14.8 item 4) |

Total ≈ 2,800 code lines (≈ 960 pass 1), most of them validation scripts and their tests, plus an
estimated 2,500 test lines. The effort the spec gives (1 week) holds only if the Slurm work queues
promptly; the searches and the de novo runs dominate wall time, not code.

## Dependency DAG

```text
Pass 1   01 ─► 02 ─► 03 ─► 04 ─► 05
Pass 2   05 ─► 06 ─┐
         05 ─► 07 ─┤
         01 ─► 08 ─► 09 ─► 10 ─┤          (08 "fails": 08 ─► 11a ─► 11b ─┤)
         04 ─► 12 ─► 13 (post-gate, not on the gate path)
                                └─► 14 ─► 15 (Khoa reviews REPORT.md)
```

## Waves

| Wave | Issues | Notes |
|---|---|---|
| A (now) | 01; 08 in parallel (env build and weights do not need the raw data; its smoke run does) | 01's resolve job is a few minutes; the download is one job |
| B | 02 | two 16-CPU jobs at once |
| C | 03, then 04 | 04 runs the arms one after the other (R14, below) |
| D | 05; 09 | 05 is the end of the MVP pass |
| E | 06, 07, 10, 12 | 10's GPU jobs and 12's benchmark jobs are sized to fit together under the caps |
| F | 14, then 15 | 15 ends with Khoa's review; 13 starts after 12 and runs past the gate |

## What M10 depends on that is not merged (2026-10-23)

Pass 1 needs nothing unmerged: `giap convert --from diann`, `giap run`, `giap compare` and V16 are
on main. The **report-freeze runs** (M10-14) wait for the hardening below, so the reviewed report
shows the final methods; pass-1 runs are MVP evidence and are re-run then (DIA-NN is not re-run:
its inputs do not change).

| Work | Status | What M10 needs it for |
|---|---|---|
| MVP-04 K-12 evidence run (acceptance job not yet submitted, PROGRESS next action 2) | spec and code merged (#247, #258) | the first full real-genome run; lessons before 04 (soft) |
| M6-12 part C (cut-offs into `calibration.yaml`, ADR-0022, Khoa) | parts A, B, C0 merged | operon flags without `uncalibrated` in the frozen runs |
| M7-11b (#267), ADR-0049 flags (#269), M7-13b (#271), M7-15 (integration hardening) | open | metric 6, TE2, `novel_orfs.tsv` class FDR |
| M8-15 (#266), M8-16 (#265) | open | metric 5 "enriched categories" and the determinants section |
| M9-16 (#228), M9-16b, M9-17 (#272), M9-19/20 (#275, #276), M9-26 (#274), M9-25 (Khoa) | open | per-arm `report.html` content; not needed for `REPORT.md`, which is written from tables |

## Slurm budget (R14: ≤ 128 CPUs, ≤ 4 GPUs, ≤ 20 queued or running `giap-*` jobs at once)

Every number below is a request, not a measurement; measured values go to PROGRESS › Measured
numbers. Before each `sbatch` the lead sums the CPUs of running and pending `giap-*` jobs of every
agent (they share the cap).

| Issue | Job (script) | Partition / QOS | CPUs | GPUs | Mem | Time limit (expected) | Count |
|---|---|---|---|---|---|---|---|
| 01 | resolve the record (`fetch_pxd078459.sbatch --resolve`) | normal / dev | 1 | 0 | 2 GB | 30 min (2 min) | 1 |
| 01 | download (`fetch_pxd078459.sbatch`) | normal | 2 | 0 | 8 GB | 24 h (1–4 h for 5–10 GB) | 1 |
| 01 | acquisition mode (`acquisition_mode.py`) | dev | 1 | 0 | 4 GB | 30 min (< 2 min) | 1 |
| 02 | `diann_pilot.sbatch A1`, `… A2` | normal | 16 each | 0 | 64 GB | 24 h (2–8 h, unmeasured) | 2, concurrent: 32 CPUs |
| 02, 03 | report samples, `giap convert` | dev | 4 | 0 | 16 GB | 1 h | 1–2 |
| 04, 06, 10, 14 | `giap_run.sbatch` (head) + Nextflow tasks | normal | 2 + ≤ 8 × 16 | 0 | 8 GB + tasks | 2 d (35–60 min per arm, §13.2 estimate) | 1 head + ≤ 8 tasks per arm |
| 08 | `build_envs.sbatch cascadia` | normal | 4 | 0 | 16 GB | 2 h | 1 |
| 08, 10 | `.d → mzML` per run | normal | 8 | 0 | 32 GB | 4 h | 1 (08), array of 3 (10) |
| 08, 10 | `denovo_pilot.sbatch` per run | gpu | 8 | 1 × L40S | 64 GB | 24 h (unknown) | 1 (08), array of 3, `%3` (10) |
| 11a | PXD053727 download; UTI89 Bakta full (`record_bakta_full.sbatch`) | normal | 2; 16 | 0 | 8; 64 GB | 24 h; 1 h | 1; 1 |
| 11b | Casanovo per run | gpu | 8 | 1 | 64 GB | 24 h | array of 3 |
| 12 | `benchmark_standard.sbatch` (standard and core, 3 repeats each) | normal, `--constraint` one class | 16 each | 0 | 64 GB | 3 h (16–55 min) | 6, at most 3 at a time: 48 CPUs |

**One `giap run` can already exceed the CPU cap.** `workflow/conf/farmshare.config` allows
`queueSize = 8` tasks of up to 16 CPUs (128) and its comment counts only those; the head job adds
2, so a single arm can reach 130 CPUs, and two concurrent arms 260. M10 therefore runs the arms
**one at a time** and with `--threads 8` when any other `giap-*` job is running (8 × 8 + 2 = 66);
the profile fix (queueSize 7, or the head counted) is a backlog item for the lead (question L1).

## Licences and what may be committed (R5, R11)

| Item | Licence (to confirm in the named issue) | Committed | Never committed |
|---|---|---|---|
| DIA-NN 2.3.1 Academia | proprietary; free for academic use; no redistribution (Khoa accepted it, M3-24) | `diann_A{1,2}.cfg`, README, the image's path and sha256, 1,000-row samples of the report (outputs are tables, not binaries: §15.2 step 7), the converted generic tables | the image (`.sif`), any DIA-NN binary, spectral libraries (`*.speclib`, `*_empirical.parquet`: large, re-creatable), full `report.parquet` files |
| PXD078459 raw runs | the record's data licence, recorded in the README by 01 | `validation/pilot/inputs/pxd078459_files.tsv` (names, URLs, sizes, checksums), `acquisition.json` (metadata only) | any raw file or archive (`.d`, `.zip`, `.tdf`, `.tdf_bin`, `.mzML`); the downloaded PX XML record stays under `$GIAP_DATA_DIR` (question L4) |
| Cascadia 0.0.7 | expected Apache-2.0; read the LICENSE at the pinned tag (08); the weights may carry their own terms | `workflow/envs/cascadia.yml` + lock, the weights' URL and sha256, a 1,000-row result sample | the weights (`.ckpt`), mzML conversions, full results |
| `.d → mzML` converter | depends on the route (ADR-M10-B): ProteoWizard's vendor-reader image needs the Bruker licence accepted; open readers (e.g. timsconvert) may load Bruker's proprietary TDF-SDK library | the route, version, image digest or env | vendor libraries; converted spectra |
| Casanovo 5.2.1, ThermoRawFileParser (fallback only) | Apache-2.0 (both); Thermo's RawFileReader inside has its own terms; weights from GitHub releases | env files, weights' URL and sha256, sample | weights, raw data, spectra |
| UniProt UP000000625 | CC BY 4.0 | nothing new (bulk file > 5 MB stays in `$GIAP_DATA_DIR`, §11.4) | the FASTA |

All external tools run as separate processes in their own environments (R5); GIAP imports none
of them. Every committed file stays ≤ 5 MB, and `tests/data/` stays < 50 MB (R11).

## §15.6 pass criteria → where the evidence comes from

| §15.6 bullet | Issue | Test or check |
|---|---|---|
| GIAP runs end to end on A1, A2 and A3 | 04 (A1, A2), 10 or 11b (A3), re-run in 14 | exit 0 of each `giap_run.sbatch` job; every §9.1 table present; job IDs in PROGRESS |
| ≥ 99.5 % of A2 peptides map exactly (TE1) | 05 | `validation/pilot/metrics/te1.json` `share_exact_inframe ≥ 0.995` (the definition is fixed in 05 before the run, question L5) |
| A1 and A2 detections and abundances are compared (`giap compare`) | 04, 10 | `CMP/tables/compare_{detection,abundance,stability}.tsv` validate against their schemas; ρ with its CI in `REPORT.md` |
| Measured runtimes replace the estimates | 12 | `benchmarks/baseline.json` entries with node class, ADR-M10-D, user-guide table updated |
| TE2–TE5 are reported; a failure is investigated and explained | 06 | `validation/pilot/metrics/te_real.json` has one entry per TE with `pass`/`fail` and, for each fail, a `REPORT.md` section with the diagnosis |
| (gate) Khoa has reviewed the report | 15 | Khoa's answer recorded in PROGRESS › Answered by Khoa |

§15.5 metric 8 and the §15.7 hypotheses are reported, not gated (§15.7: "reported where the pilot
data allow"). Hypothesis 2 (FDRBench) is the full validation's and is marked "deferred" (§16).

## ADRs M10 owes (numbers: next free at merge time)

- **ADR-M10-A "DIA-NN 2.3.1 Academia instead of 2.7.0"** (02): the version, the image, the
  two-step workflow, and the template deviations (`--met-excision` out, `--reanalyse` and
  `--rt-profiling` out, extra m/z ranges and `--peptidoforms`). Khoa accepted the README on
  2026-10-20, but a dependency deviation from the spec is an ADR (R13, §14.8 item 1); this one
  records his answer.
- **ADR-M10-B "`.d → mzML` route for A3"** (08): a new tool and its licence.
- **ADR-M10-C "PXD053727 source"** (11a, conditional): §3.2 lists iProX/ProteomeXchange for
  PXD078459 only; the UTI89 runs are hosted where their record says (probably PRIDE). Khoa.
- **ADR-M10-D "Benchmark definition and first baseline"** (12): what a "standard run" is (input,
  modules, live Bakta or recorded JSON), the node class, and the measured numbers against the
  §13.2 estimate.
- Conditional: if the acquisition mode is not diaPASEF (01), the DIA-NN settings change; that is
  an amendment to ADR-M10-A, with Khoa's re-check.

## Risks (and the Blocked-log plan)

R6 governs every download: one later retry, then a row under PROGRESS › Blocked (date, what, host,
what is needed from Khoa), and work moves to issues that do not need the data. Never an
unofficial mirror or a re-hosted copy, never TLS verification off.

| Risk | Signal | Plan |
|---|---|---|
| PXD078459 files are not yet public (embargoed or under review), or iProX requires a login | the record lists no file URIs, or the URIs answer 401/403 or a login page | Blocked row "PXD078459 raw files login-gated / private"; ask Khoa (he may have access or know the authors). Pass 1 stops at 01; 08's env build and 12's benchmark (on simulated peptides) continue |
| iProX serves only plain `ftp://` | the record's URIs are `ftp://` | M3-24 decision 11: do not download over plain FTP; ask whether the same host's `https://` path is acceptable (question L2); Blocked row until answered |
| The M3 default URL is wrong | HEAD 404 on 2026-10-23 | 01 takes URLs from the record only |
| Checksum mismatch or truncated archive | record checksum differs, `unzip -t` fails | retry once, then Blocked; never accept a partial run |
| Acquisition mode is not diaPASEF (ddaPASEF, or mixed across runs) | 01's `acquisition.json` | stop before 02; amend the cfgs (DIA-NN handles DDA only partly), ADR-M10-A amendment, Khoa re-checks the flags (§17) |
| The DIA-NN image is replaced or altered (it is world-writable) | sha256 differs from the pin | the job refuses to run; ask Khoa for the authentic image; never fetch DIA-NN elsewhere |
| DIA-NN 2.3.1 rejects a cfg flag, or a column the converter expects is missing | job log; 03's contract test | fix the cfg (lead decision recorded in the README) or the converter; R9 after 3 attempts |
| Cascadia cannot read timsTOF data, or its output is degenerate | 08's smoke run and decision rule | the §15.1 fallback (11a/11b); ADR-M10-C goes to Khoa first |
| The `.d → mzML` route needs a vendor licence | ProteoWizard's vendor-reader image, Bruker SDK terms | Khoa accepts (§14.8 item 6 by analogy) or an open route is chosen; Blocked row until then |
| Cascadia or Casanovo weights unavailable | download fails | Blocked row; R6 |
| GPU queue waits | `PENDING (Resources)` on `gpu` | arrays capped at 3 GPUs; never above 4 |
| One `giap run` exceeds 128 CPUs (above) | `squeue` sum | arms one at a time, `--threads 8` when sharing |
| A3 on UTI89 cannot be compared with A1/A2 | `giap compare` exit 1 (different genome, ADR-0046) | the fallback reports A3 separately; §15.5 metric 6 and hypothesis 4 get the caveat |
| A TE fails on real data | 06 | investigated and explained in the report before Khoa's review (§15.6); no threshold changes (R3) |
| Disk | raw 5–10 GB, libraries and reports unknown, work directories tens of GB | 14 cleans; the probe warns above 400 GB |
| Session allocation expires during long jobs | session ends | job IDs in PROGRESS; scripts are resumable (download `-C -`, DIA-NN step 1 skip, `-resume`) |

## Feasibility on FarmShare

- 5–10 GB over HTTPS in one 2-CPU job is a few hours at most; the class share has 197 TB free.
- DIA-NN library-free on ~4.4–4.7 k proteins plus contaminants with 3 timsTOF runs fits 16 CPU /
  64 GB (the M3-24 sizing); wall time is unmeasured, hence the 24 h limit.
- Cascadia on one L40S per run: unmeasured; a 24 h limit per run and an array of 3.
- `giap run` on a 5-Mb genome: §13.2 expects 16–55 min on 16 threads; the arms run in sequence.

## Decisions taken in this index (for the lead to confirm)

1. **MVP first:** A1 and A2 reach `giap compare` and TE1 (01–05) before A3, the benchmark or any
   hardening. A3 feasibility (08) runs in parallel because it shares no code.
2. **Version of record is DIA-NN 2.3.1 Academia** (M3-24 decision 2) with ADR-M10-A; 2.7.0 is not
   awaited. If it arrives, only the cfgs, the README and the recorded sample's version change.
3. **Arms run one at a time** under the R14 cap; `--threads 8` when sharing.
4. **Report-freeze runs:** the arms are re-run once on the final main in 14, so the reviewed report
   matches the released methods; DIA-NN and the de novo search are not re-run.
5. **Committed derived data:** `validation/pilot/inputs/` holds the gzipped generic tables (§15.2
   step 7), the file manifest and `acquisition.json`; `validation/pilot/metrics/` holds the metric
   TSV/JSON files; nothing else from real data.
6. **Samples for contract tests** come from the real reports: the first 1,000 rows in file order
   of each report, every column, gzipped, with a `PROVENANCE.md` (tool version, image sha256, job ID,
   input checksums). The DIA-NN sample lives under `tests/data/recorded/diann/2.3.1/` (the version
   that produced it, as the README says), not `2.7.0/`.
7. **REPORT.md is written from tables** (`validation/pilot/metrics/`), never from screenshots or
   numbers typed by hand; each number cites the file it comes from.
8. **The speed loop (13) is post-gate:** the M10 gate needs measured runtimes (12), not a faster
   pipeline.

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| L1 | `farmshare.config` lets one run reach 130 CPUs (head + 8 × 16). Fix the profile (`queueSize = 7`) as a chore, or rely on M10's one-arm-at-a-time rule? | one arm at a time now; a backlog chore for the profile |
| L2 | If the record gives only `ftp://` URIs, may the job use `https://` on the **same host and path** when a HEAD returns 200 with the record's byte size? (Not a mirror; same official host.) | yes, recorded per file in `pxd078459_files.tsv`; else Blocked |
| L3 | More than three wild-type runs in the record: which three? | the first three by natural sort of file name, all candidates listed in the manifest (M10-01 Q2) |
| L4 | Commit the PX XML record itself? | no: only its URL, retrieval date and sha256; the derived manifest is committed |
| L5 | TE1's denominator on real data: distinct stripped target sequences at q ≤ 0.01 (run-level? global?), contaminant-bucket peptides out, N-terminal peptides of non-ATG genes out (as in the in-silico proxy)? | `Global.Q.Value` ≤ 0.01 as converted; contaminants out; non-ATG N-terminal products counted and reported separately but kept in the denominator (§11.6 names the exception only for the proxy) |
| L6 | "Fails on timsTOF" for Cascadia needs a rule before the smoke run. | fail when, on run 1, fewer than 500 distinct de novo peptides of 7–30 aa pass Cascadia's own score cut-off, or fewer than 20 % of the top 1,000 by score map exactly to the Nissle six-frame translation (random 7–30-mers almost never do) |
| L7 | Benchmark "standard run": live Bakta full DB from the FASTA, or the recorded JSON? Node class? | both: "standard" = live Bakta + default modules; "core" = recorded JSON + core modules; class `CPU_SKU:6148` (wheat, 8 nodes, the class of most measured numbers so far) |
| L8 | Which documents carry "measured runtimes replace the estimates"? (GIAP_SPEC.md cannot be edited, R13.) | `benchmarks/baseline.json`, the user-guide runtime table, PROGRESS › Measured numbers, and ADR-M10-D stating measured vs §13.2 |
| L9 | A3 fallback (UTI89): without a same-strain database search there are no `--reference-peptides`, so de novo precision c and hypothesis 4 are unestimated. Run DIA-NN on the DDA runs for a reference, or report c as `unestimated`? | `unestimated`, stated in the report (R4) |
| L10 | Should the pilot runs enable the function module (KEGG content; MVP-04 turned it off)? Metric 5 needs enriched categories. | on, with the KEGG retrieval date in the report (ADR-0017); no KEGG-derived content committed |

## Questions that need Khoa (not the lead)

1. The gate: review of `validation/pilot/REPORT.md` (M10-15).
2. ADR-M10-A (DIA-NN 2.3.1 as the version of record; formalises his 2026-10-20 answer).
3. Conditional: access if PXD078459 is login-gated or private; the vendor-licence acceptance if
   the A3 conversion route needs one (ADR-M10-B); ADR-M10-C if the A3 fallback is triggered.
4. Conditional: a re-check of the DIA-NN flags if the runs are not diaPASEF.

## Lead decisions (2026-10-23)
1. **L1–L10 and M10-01's Q1–Q5:** the proposed defaults are taken, all reversible:
   - one arm at a time, with a backlog chore to cap the farmshare profile at 128 CPUs including the head job;
   - https on the same host when HEAD gives the record's byte size;
   - the first three runs by natural sort;
   - the record XML is not committed;
   - TE1's denominator as proposed;
   - the Cascadia failure rule as proposed;
   - the standard and core benchmark variants on wheat;
   - runtimes go into baseline.json, the guide, PROGRESS and ADR-M10-D;
   - de novo precision is `unestimated` in the fallback;
   - the function module is on.
2. **The DIA-NN image:** the shared `.sif` belongs to another user and is world-writable. M10-02 checks its sha256 against the pin. On a match, it copies the image to `$GIAP_ROOT/images/diann_2.3.1.sif` with mode 0444, and every job checks and uses that copy. It is never committed.
3. **ADR-M10-A (DIA-NN 2.3.1 vs the spec's 2.7.0):** M10-02 drafts it. Pass 1 proceeds on 2.3.1 as Khoa accepted in M3-24; the ADR records it and goes to Khoa in the next batch.
4. **The wrong iProX URL in the M3-24 fetch script:** M10-01 replaces it, using names and URLs from the official record only (R6).
