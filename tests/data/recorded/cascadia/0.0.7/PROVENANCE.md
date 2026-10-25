# Cascadia 0.0.7: A3 de novo pilot, WT_1 (provenance only)

M10-08b, deliverable 13 (docs/specs/M10-08-a3-feasibility.md). The A3 verdict on WT_1 is a
**pass** (`validation/pilot/metrics/a3_feasibility.json`: `n_pass` 863, `n_top` 863,
`n_match` 549, `share` 0.6362, 95 % Wilson CI [0.6035, 0.6676]). This directory holds this note
only: no `.ssl`, no peptide, no data (acceptance item 12; Khoa, 2026-10-23, M10-02 decision 4;
spec Q4).

## Recorded sample (not committed)

- Path: `$GIAP_ROOT/pilot/recorded/cascadia/0.0.7/WT_1.ssl.gz` on FarmShare.
- Content: the `.ssl` header plus the top `n_top = min(1000, n_pass) = 863` rows (spec steps 1-7:
  strip, I -> L, length 7-30, `score > 0.8`, one row per folded sequence at its maximum score,
  sorted by score descending then folded sequence ascending). Fewer than 1,000 rows because
  `n_pass` = 863 < 1,000. The rank order equals the evaluator's `feasibility_top1000.tsv`.
- 864 lines (header + 863 rows), 34,258 bytes, gzip with mtime 0;
  sha256 `38873d141daadbbd43b0d72f1d7ae01878499f1decfe5fc70ce3336c6484cc36`.
- Source: `$GIAP_ROOT/pilot/A3_cascadia_cpu/WT_1/WT_1.ssl` (1,445 rows; sha256
  `601344b23f89c3caf8142b4a2569187f1a40a8af2e3017377375267734a3712c`).

## Tool and environment

- Cascadia **0.0.7**, from PyPI: wheel `cascadia-0.0.7-py3-none-any.whl`, sha256
  `50b5204349ab91d6bc34f9313ad18e9d2d305814d19e53cd1cd75d8208f092a3`. Licence Apache-2.0.
- Env `giap-cascadia`, built by job 1888993 from the committed spec and locks (commit `b94a2ec2`):
  - `workflow/envs/cascadia.yml` sha256 `11832ac20c7178c0d4714afd6f66bd6ab15e3497cebb353ddfab5dfcba90762b`
  - `workflow/envs/locks/cascadia.explicit.txt` sha256 `5c4c5d65ef246b8f4045cef63e981077e5baf05bdfed3f90b8b7b9635077235c`
  - `workflow/envs/locks/cascadia.pip.txt` sha256 `e6f8980ce1735236e48e3c18ecce02df5279cc762338cf91e3617d7c7c698abd`
- Weights `cascadia.ckpt` (stored as `$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt`), 584,332,141 bytes,
  sha256 `0f5e91413c716b22019699c8cb25655216f9c96d0f2aeb6ec42084bb65652192`
  (`validation/pilot/a3_assets.tsv`; host per ADR-0054).
- Unimod table (psims 1.4.0 vendored), sha256
  `7b303537f45f8f0d7ddd1e5b19cf59060474eafb18e33bef0b17fef078ef4c38`.
- Score cut-off 0.8 (Cascadia's CLI default, `cutoff_source` `cli_default`); batch size 32.

## Route and input

- Route **B2**: TIMSCONVERT at commit `3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b`
  (`route_version`), with pyTDFSDK `4808da81d3e37f82da96f392dbb446c3733edacb` (ADR-0054).
- Input: PXD078459 run WT_1, `MLA1914689_WT_P1_1_Y-D7_1_13119.d.zip`, 1,001,043,899 bytes,
  sha256 `58d09a41d6642f35f0cb25e88e4a715fbfe0c30f7b4f839c9fd30601bd99f013`
  (`validation/pilot/inputs/pxd078459_files.tsv`).
- mzML `$GIAP_DATA_DIR/pxd078459/mzml/WT_1.mzML`, sha256
  `691556c139e9ebfc31f64413541ea39eaa965f78c2b1b0bde1510fb145bcc75f`; `mzml-check` M1-M8 all pass.
  143,064 spectra, 4,471 batches.
- Genome for the match: `tests/data/genomes/nissle1917.fna.gz`, sha256
  `fe74801661e688a525eef7a7a6be6210659ebc9a2e61691efb1f10b65ea623ee`.

## Slurm jobs (checked with `sacct`, 2026-10-25)

| Job | Name | What | State |
|---|---|---|---|
| 1883978 | giap-build-envs | cascadia env, first build (CUDA solve on a CPU node) | FAILED |
| 1884085 | giap-build-envs | cascadia env rebuild | COMPLETED |
| 1888993 | giap-build-envs | cascadia env with `setuptools<81` (the env used) | COMPLETED |
| 1884175 | giap-fetch-a3-assets | weights and Unimod, `--resolve` | COMPLETED |
| 1884404 | giap-fetch-a3-assets | weights and Unimod, fetch | COMPLETED |
| 1884087 | giap-mzml-pilot | WT_1 `.d` -> mzML (B2) and `mzml-check` | COMPLETED |
| 1889426 | giap-cascadia-cpu | WT_1 Cascadia, CPU shard 0 (task 0 completed; tasks 1-5 cancelled, replaced by the work queue) | COMPLETED (task 0) |
| 1890599 | giap-cascadia-work | CPU work-queue array, tasks 0-7 | COMPLETED |
| 1890876 | giap-cascadia-work | CPU work-queue array, tasks 0-1 (tasks 2-5 cancelled before start) | COMPLETED (0-1) |
| 1890902 | giap-cascadia-merge | WT_1 chunk merge -> `WT_1.ssl` (23:04-23:05) | COMPLETED |
| 1890903 | giap-a3-evaluate | WT_1 `a3_feasibility.py evaluate` -> the verdict | COMPLETED |
| 1890952 | giap-cascadia-work | CPU work-queue array, tasks 0-3 | COMPLETED |
| 1891287 | giap-cascadia-work | CPU work-queue array, tasks 0-13 | COMPLETED |
| 1891944 | giap-cascadia-work | CPU work-queue array, tasks 0-13 | COMPLETED |

The verdict JSON's `job_ids` are 1883978, 1884085, 1884175, 1884087, 1884404, 1889426, 1890903.
The work queue is shared by WT_1, WT_2 and WT_3; the WT_1 chunks were all written by 23:04 on
2026-10-24 (before the merge 1890902). Arrays 1890952, 1891287 and 1891944 started after the
merge and served WT_2/WT_3 only.

## The `.ssl` header, as defined in the Cascadia source

Taken from the source code of the installed package, not from any data file:
`cascadia/utils.py`, function `write_results` (line 27), header written at **line 50**, in
`$MAMBA_ROOT_PREFIX/envs/giap-cascadia/lib/python3.10/site-packages/cascadia/utils.py`
(file sha256 `457723db67331f82bdb3d7770189892b13cd82df4acae4c14071dc1dcab5e53e`, matching the
wheel's `RECORD` entry `sha256=RXcj22czH4K9s9d3AYmJKxPNgt9KyuTBQHHcHcq15T4`):

```python
out.write("\t".join(['file',	'scan',	'charge',	'sequence',	'score-type',	'score',	'retention-time',	'start-time',	'end-time']) + '\n')
```

Tab-separated, nine columns, in this order:

| # | Column |
|---|---|
| 1 | `file` |
| 2 | `scan` |
| 3 | `charge` |
| 4 | `sequence` |
| 5 | `score-type` |
| 6 | `score` |
| 7 | `retention-time` |
| 8 | `start-time` |
| 9 | `end-time` |

The file is named `<-o value>.ssl` (`open(results_file + '.ssl', 'w')`, line 49). Rows are
written only when `conf > score_threshold` (strict; line 57).
