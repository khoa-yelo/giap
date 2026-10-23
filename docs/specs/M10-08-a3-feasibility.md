# M10-08: A3 de novo feasibility: Cascadia on the timsTOF pilot runs

Milestone M10 · **pass 2** (may run in parallel with pass 1: it shares no code) · depends on:
M10-01 (the raw runs and `acquisition.json`; the environment, weights and image steps do not need
them) · estimate ≈ 330 code lines (above the index's 150: see question Q7) + ≈ 400 test lines ·
not `stats` (one Wilson interval, reported only) · roles: implementer (env file, registry entry,
`scripts/`), validation (`validation/pilot/`) · **Khoa: needed** if the `.d → mzML` route loads a
vendor library (options B1, B2 below), and if the weights carry restrictive terms or sit on a
host outside §3.2.

## Spec sections implemented

- §15.1 arm A3: "Cascadia 0.0.7 on the same runs (exploratory), on the `gpu` partition
  (1 × L40S). If DIA de novo fails on timsTOF data: Casanovo 5.2.1 on 3 UTI89 DDA runs
  (PXD053727; genome GCA_000013265.1)". This issue decides **whether** the fallback is needed;
  the fallback itself is M10-11a/11b.
- §15.2 step 6 (`denovo_pilot.sbatch`, `--partition=gpu --gres=gpu:1`) and step 7 (a 1,000-row
  sample of the raw Cascadia result under `tests/data/recorded/cascadia/0.0.7/`).
- §17 row "Cascadia fails on timsTOF data | A3 empty or degenerate | Casanovo on 3 UTI89 DDA
  runs". Index question L6 makes "empty or degenerate" measurable (lead decision 1 took it).
- §3.1 (`gpu` partition: `oat-01..06`, 4 × L40S each, QOS `gpu` ≤ 4 GPUs; checked 2026-10-23:
  partition limit 2 days, GRES `gpu:4` without a type name), §3.2 sources, §3.8 job header, §3.9
  disk budget ("de novo outputs: measure"), R5 (separate process, own environment), R6, R11, R14.
- §5.4 row "Cascadia 0.0.x | its result file (mzTab-like) | as Casanovo; `source` = denovo": this
  issue records the real file; the converter is M10-09.
- §7.2.1–7.2.2 (six-frame translation, exact route): the genome check of the failure rule reuses
  `giap.mapping.sixframe` and `giap.mapping.exact` unchanged.

## What exists (main d29a2ea5 + m10/index)

| Piece | State |
|---|---|
| `scripts/slurm/denovo_pilot.sbatch` | M3-24: `gpu`, `--gres=gpu:1`, 8 CPU, 64 GB, 24 h, job name `giap-denovo`; env `giap-cascadia`; weights default `$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt`; input glob `*.d` under `$GIAP_DATA_DIR/pxd078459`; runs **all** runs in one call; the command line `cascadia sequence --model … --num_workers … --output …` is a **guess** from M3, never checked against the docs. Exits 2 while the env or the weights are missing |
| M3-24 tests on it (`tests/unit/test_m3_24_pilot_kit.py`) | partition `gpu`, GRES `gpu(:type)?:1`, name `giap-denovo*`, the text names `Cascadia`, `0.0.7`, `Casanovo`, `5.2.1`, `PXD053727`, `UTI89`; no `diann` in its code; **no network** (`curl`, `wget`, `aria2c`, `git clone`) in its code; outputs only under `$GIAP_ROOT`; the kit's CPU sum ≤ 128. **They must keep passing unchanged** |
| `scripts/build_envs.sh` / `build_envs.sbatch` | builds `workflow/envs/<tool>.yml` only for tools found in the registry (`giap.annotate.tools.TOOLS`, or the frozen `fgsea`); writes `workflow/envs/locks/<tool>.{explicit,pip}.txt`; `build_envs.sbatch` also pulls every digest in `workflow/conf/containers.config` |
| Six-frame and exact matching | `build_sixframe(genome, replicons)` (table 11, I folded into L, `*` at stops, `X` for ambiguous codons, circular replicons wrapped by `WRAP_NT` = 126 nt); `build_decoys(index, "shuffle", 50, rng)`; `match_exact(queries, index, decoys, decoy_of, min_len, max_len)`; `giap.rng.for_step(name, seed)`; `config.seed` default 20261014; `mapping.decoys.denovo_input = "shuffle"` |
| Nissle genome | `tests/data/genomes/nissle1917.fna.gz`: **four** circular-tagged records, 5,064,940 bp: CP171242.1 chromosome, CP171243.1 plasmid **p15A-ILA-dapA**, CP171244.1 pMUT1, CP171245.1 pMUT2 (question Q3); recorded full-DB Bakta JSON in `tests/data/recorded/bakta/1.12.2-full/nissle/` |
| Missing | `workflow/envs/cascadia.yml`, its locks, the registry entry, the weights, any `.d → mzML` route, `tests/data/recorded/cascadia/`, any ADR on de novo or mzML conversion (none of ADR-0001…0049 covers them) |

### What a HEAD request showed (2026-10-23, from `iron-13`; no download)

| URL | Answer | Meaning |
|---|---|---|
| `https://github.com/Noble-Lab/cascadia` | 200 | the official repository exists; `blob/main/LICENSE` 200 |
| `…/cascadia/releases/latest` | 302 to `/releases` | **no GitHub release**; tags `v0.0.7` and `0.0.7` answer 404 |
| `https://pypi.org/pypi/cascadia/0.0.7/json` | 200 | **0.0.7 exists on PyPI** (`0.0.8`, `0.1.0`: 404). That the PyPI project is the Noble Lab's must be checked from its metadata (`project_urls`) in the build job |
| `https://cascadia.readthedocs.io/` | 302 to `/en/latest/` (200) | the §18 docs; no versioned docs path checked |
| `anaconda.org/bioconda/cascadia` | 200 | **inconclusive**: anaconda.org answers 200 for a made-up name too |
| `hub.docker.com/v2/repositories/proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses` | 200 | the ProteoWizard vendor-reader image exists on Docker Hub (a §3.2 source) |
| `https://github.com/gtluu/timsconvert` | 200 | TIMSCONVERT exists; `blob/main/LICENSE` 404 (licence file name or branch unknown); `pypi.org/pypi/timsconvert` 404 |
| `https://github.com/MannLabs/alphatims` | 200 | AlphaTims (open TDF reader) exists |

Where the **weights** are published is not known: there is no GitHub release to hold them. The
build step reads the pinned version's README/docs for the official URL (question Q11).

## Deliverables

| # | File | Owner | What |
|---|---|---|---|
| 1 | `workflow/envs/cascadia.yml` + `workflow/envs/locks/cascadia.{explicit,pip}.txt` | implementer | conda-forge/bioconda only (`--override-channels`), Python and CUDA-enabled PyTorch at the versions `cascadia==0.0.7` declares; `cascadia==0.0.7` from PyPI in the `pip:` block. Header comment: source (PyPI 0.0.7, project URL checked to be `github.com/Noble-Lab/cascadia`), licence as read from the **0.0.7 distribution's** LICENSE file and metadata, date. The locks are written by `build_envs.sh` |
| 2 | registry entry `cascadia` | implementer | a `PILOT_TOOLS` mapping in `giap.annotate.tools` (version `0.0.7`, licence as found, `version_cmd = _pkg_version("cascadia")`, env `giap-cascadia`), consulted by `build_envs.sh`'s `registry()` after `TOOLS` and the frozen map. **Not in `TOOLS`**, so `required_tools`, `build_envs.sh all` and a run's tool versions never see it (question Q4) |
| 3 | `validation/pilot/a3_assets.tsv` | validation | committed manifest, one row per external artefact of A3: `asset` (`cascadia_weights`, `mzml_converter`), `version`, `url` (official source), `sha256`, `bytes`, `licence`, `path` (under `$GIAP_ROOT`), `retrieved` (UTC date). No binary is committed |
| 4 | `scripts/slurm/fetch_a3_assets.sbatch` | implementer | the **only** A3 script that reaches the network: downloads the weights row (https only, no TLS-off flag, `--continue-at -`), checks `sha256`, renames to `$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt`, mode 0444; pulls the converter image **by digest** with `apptainer pull` to `$GIAP_ROOT/images/<name>@<short-digest>.sif`, mode 0444. `--resolve` prints the sha256/bytes of what it fetched for the agent to commit into deliverable 3 (first fetch: TOFU from the official source, recorded; every later fetch must match) |
| 5 | `scripts/slurm/mzml_pilot.sbatch RUN [--slice SECONDS]` | implementer | `.d → mzML` for one run through the ADR-M10-B route (below); writes `$GIAP_DATA_DIR/pxd078459/mzml/<RUN>.mzML` (or `<RUN>.slice<SECONDS>.mzML`), then the mzML sanity check (deliverable 7, `--mzml-check`) |
| 6 | `scripts/slurm/denovo_pilot.sbatch` (changed) | implementer | `cascadia` branch takes **one run** (`RUN=WT_1` env or argument, default all, as now), reads `$GIAP_DATA_DIR/pxd078459/mzml/<RUN>*.mzML` (or `.d` under option B0), writes `$GIAP_ROOT/pilot/A3_cascadia/<RUN>/`; the command line comes from the pinned version's documentation (`cascadia --help` captured into the job log first). Checks the weights' sha256 against deliverable 3 before running (exit 2 on mismatch). Every M3-24 assertion still holds (no network, `gpu:1`, names) |
| 7 | `validation/pilot/a3_feasibility.py` | validation | the L6 rule (below); CLI `python -I validation/pilot/a3_feasibility.py RESULT --genome G --annotation J --out OUT_JSON [--top-tsv T] [--cutoff X]`, plus `--mzml-check MZML --acquisition validation/pilot/inputs/acquisition.json`. Imports GIAP (`giap.mapping`, `giap.annotate.bakta`, `giap.rng`, `giap.config`); stdlib/pandas reader for the Cascadia result; never imports Cascadia (R5) |
| 8 | `validation/pilot/metrics/a3_feasibility.json` | validation | the committed verdict (fields below) |
| 9 | `tests/data/recorded/cascadia/0.0.7/` | validation | **on a pass only** (question Q6): the first 1,000 rows in file order of run 1's result, every column, gzipped, plus `PROVENANCE.md` (Cascadia version, env lock sha256, weights sha256, converter route and image digest, job IDs, input `.d` checksum from `pxd078459_files.tsv`), as index decision 6 |
| 10 | `docs/adr/ADR-NNNN-a3-mzml-route.md` (ADR-M10-B) | lead (drafted here) | the route chosen from the options below; `Needs Khoa: yes` for B1/B2 |
| 11 | `validation/pilot/README.md` "A3 de novo" section | validation | versions, licences, assets, route, commands, job IDs, sizes, wall times, the verdict |

Nothing else is committed: no weights, no image, no mzML, no full Cascadia result, no top-1,000
table from real data (it stays under `$GIAP_ROOT`; index decision 5).

### Where things live under `$GIAP_ROOT` (never committed)

| Item | Path | Expected size |
|---|---|---|
| Cascadia env | `$MAMBA_ROOT_PREFIX/envs/giap-cascadia` | 5–8 GB (PyTorch + CUDA); `micromamba clean -a` after |
| Weights | `$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt` (0444; original file name in deliverable 3) | 0.1–1 GB (unknown) |
| Converter image | `$GIAP_ROOT/images/` (0444) | 2–5 GB (B1, Wine + ProteoWizard) |
| mzML | `$GIAP_DATA_DIR/pxd078459/mzml/` | measure: 2–10 GB per run, peak-picked, zlib |
| Cascadia result | `$GIAP_ROOT/pilot/A3_cascadia/<RUN>/` | small (tables) |
| Evaluator table | `$GIAP_ROOT/pilot/A3_cascadia/<RUN>/feasibility_top1000.tsv` | < 1 MB |

## The `.d → mzML` route and ADR-M10-B (draft, options)

Cascadia is trained on and documented for mzML DIA data (to confirm in its 0.0.7 docs). Bruker
`.d` holds `analysis.tdf` (SQLite, open) and `analysis.tdf_bin` (compressed peaks). Every
converter that writes calibrated m/z and mobility from `analysis.tdf_bin` either uses Bruker's
proprietary TDF-SDK library (`timsdata.dll` / `libtimsdata.so`, Bruker licence terms) or
re-implements it.

**Draft ADR-M10-B "`.d → mzML` route for arm A3"** (Status: proposed · Needs Khoa: yes, unless B0)

| Option | Tool, source, licence | Vendor library or licence acceptance | Fit | Verdict |
|---|---|---|---|---|
| **B0** Cascadia reads `.d` natively | Cascadia 0.0.7 itself | none | only if the 0.0.7 docs/`--help` list `.d` input | **preferred if it exists**; no ADR dependency beyond Cascadia |
| **B1** ProteoWizard `msconvert`, Docker Hub `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses`, pinned by digest, run with Apptainer (Wine inside) | ProteoWizard: Apache-2.0; the image bundles the vendor readers, incl. Bruker's TDF-SDK, under each vendor's licence | **yes**: the image name is an acceptance of every bundled vendor licence (Bruker's among them). **Khoa needed**: he accepts, as for DIA-NN (§14.8 item 6 by analogy) | the de-facto route the DIA tools document; diaPASEF supported with `--combineIonMobilitySpectra`; Wine under Apptainer needs a writable `WINEPREFIX`/`HOME` (`--writable-tmpfs`, binds; §17 Apptainer row) | **default proposal** |
| **B2** TIMSCONVERT (`gtluu/timsconvert`, bioconda) | licence to confirm (LICENSE path 404 on `main`); it bundles Bruker's `libtimsdata.so` | **yes**: Bruker TDF-SDK terms. **Khoa needed** | native Linux, no Wine; diaPASEF export to be confirmed in its docs | second choice |
| **B3** AlphaTims (Apache-2.0) reads TDF without the SDK, plus a GIAP-written mzML writer | none | no | new GIAP code that writes spectra (R1 scope; correctness risk) | rejected |
| **B4** skip Cascadia on timsTOF, go to the fallback | none | no | gives up A3 on the same spectra | only after B0–B2 fail (R9) |

Decision order (reversible default): B0 if documented; else B1 after Khoa accepts the vendor
licences; B2 if B1 cannot run under Apptainer after 3 attempts (R9); B4 only after both. Until
Khoa answers, PROGRESS › Blocked carries "A3 `.d → mzML`: vendor licence acceptance (ADR-M10-B)"
and the env build, weights and evaluator proceed.

**B1 command (proposal, to be checked against the ProteoWizard docs before the first run):**

```bash
apptainer exec --cleanenv --writable-tmpfs -B "$GIAP_DATA_DIR:$GIAP_DATA_DIR" "$PWIZ_SIF" \
  wine msconvert "$GIAP_DATA_DIR/pxd078459/WT_1.d" --mzML --zlib --64 \
  --combineIonMobilitySpectra --filter "peakPicking vendor msLevel=1-" \
  [--filter "scanTime [0,SECONDS]"] -o "$GIAP_DATA_DIR/pxd078459/mzml" --outfile WT_1.mzML
```

The options Cascadia's documentation gives for timsTOF (if any) override this proposal; the
command used is recorded in the README.

**mzML sanity check** (`a3_feasibility.py --mzml-check`, stdlib `xml.etree.iterparse`, streaming):
counts MS1 and MS2 spectra, the distinct isolation windows (target m/z ± offsets), and the RT
range. It passes when MS2 > 0, MS1 > 0, the number of distinct isolation windows is ≥ the
window-group count in `acquisition.json` for that run, and the RT range covers ≥ 90 % of the
gradient length there (100 % of the slice for `--slice`). A failed check is a **route** failure
(exit 1), never a Cascadia verdict.

## The smoke run (one GPU, one run)

Run 1 = `WT_1`, the first row of `validation/pilot/inputs/pxd078459_files.tsv`.

| Step | Script | Partition / QOS | CPUs | GPUs | Mem | Time limit (expected) |
|---|---|---|---|---|---|---|
| 1 env | `build_envs.sbatch cascadia` | normal | 4 | 0 | 16 GB | 2 h (20–40 min) |
| 2 assets | `fetch_a3_assets.sbatch` | normal | 1 | 0 | 4 GB | 2 h (10–30 min) |
| 3a slice | `mzml_pilot.sbatch WT_1 --slice 600` | normal | 8 | 0 | 32 GB | 1 h (≤ 15 min) |
| 3b plumbing | `RUN=WT_1 denovo_pilot.sbatch cascadia` on the slice (`--time=02:00:00` on the command line) | gpu | 8 | 1 | 64 GB | 2 h (unknown) |
| 4a full | `mzml_pilot.sbatch WT_1` | normal | 8 | 0 | 32 GB | 4 h (unknown) |
| 4b smoke | `RUN=WT_1 denovo_pilot.sbatch cascadia` | gpu | 8 | 1 | 64 GB | 24 h (unknown); one resubmission at `--time=47:00:00` after a TIMEOUT |
| 5 verdict | `a3_feasibility.py` (in a `dev` job: 2 CPU, 8 GB, 30 min) | dev | 2 | 0 | 8 GB | 30 min (< 5 min) |

Steps 1–2 need no raw data and run now (wave A). Steps 3–5 wait for M10-01's download. Steps
3b/4b use `--dependency=afterok:<conversion job>`. Peak use: 8 CPUs and **1 GPU** (R14 ≤ 4);
the lead sums all `giap-*` CPUs before each `sbatch`. GRES is `gpu:1` (the partition has no GPU
type names). Every job ID, `sacct` Elapsed/MaxRSS, node, `nvidia-smi` GPU name and memory, and
output sizes go to PROGRESS (In progress, then Measured numbers).

The plumbing run (3b) only proves the chain works (exit 0, a non-empty result); the verdict uses
the **full** run 1 (4b), as L6 says.

## The failure rule L6 (lead decision 1)

**Cascadia "fails on timsTOF"** when, on the complete result of run 1, **either**

- (a) fewer than **500** distinct de novo peptides of 7–30 aa pass the score cut-off, **or**
- (b) fewer than **20 %** of the top 1,000 of them match the Nissle six-frame translation exactly.

Otherwise it passes, and M10-09/M10-10 proceed. A fail triggers M10-11a (ADR-M10-C to Khoa).

### Computing (a)

1. Read every PSM row of the Cascadia result (the format of the pinned version: column names are
   located by header, listed in a small constant in the script; the format found is recorded in
   the README for M10-09).
2. **Strip** each predicted sequence: delete bracketed or parenthesised modification text
   (`[…]`, `(…)`), any signed mass (`+15.995`, `-17.027`) and any character outside `A–Z`;
   uppercase. This is the generic table's `peptide` (§5.3).
3. **Length filter:** keep 7 ≤ len ≤ 30.
4. **Score cut-off:** keep rows with `score ≥ cutoff`. `cutoff` = Cascadia 0.0.7's own default
   reporting threshold as documented or as its CLI default (recorded); if the tool documents none,
   every reported row passes and `cutoff = "none"` is recorded (question Q1). `--cutoff` overrides
   (recorded as `cutoff_source = "override"`). Higher scores are better (§5.3); a result whose
   documented score is lower-is-better is negated and the fact recorded.
5. **Distinct:** one entry per I/L-folded stripped sequence (`I → L`; the matcher folds the same
   way, so `PEPTIDEK` and `PEPTLDEK` are one peptide), keeping the **maximum** score (§5.3
   aggregation) and the first as-written form in file order.
6. `n_pass` = the number of distinct entries. (a) holds when `n_pass < 500`.

### Computing (b) with existing GIAP code

7. **Top 1,000:** sort the distinct entries by score descending, then by folded sequence
   ascending (deterministic ties); take the first `min(1000, n_pass)`. `n_top` is that number.
8. **Genome:** `tests/data/genomes/nissle1917.fna.gz` and the recorded full-DB Bakta JSON, loaded
   as `validation/proxies/te_proxies.py` does (`read_genome_fasta`, then
   `parse_bakta_json(json, genome=fasta).model.replicons`), so contig IDs, lengths and
   **circular** topologies are the ones `giap run` uses. All four records are searched
   (question Q3). The FASTA's sha256 is recorded.
9. **Index:** `index = build_sixframe(genome, replicons)` with the default `WRAP_NT`, so a
   peptide may span the origin of a circular replicon; both strands, three phases each, table 11.
10. **Decoys (diagnostic only):** `decoys = build_decoys(index, "shuffle", 50,
    for_step("a3_feasibility", 20261014))`; window and seed from the config defaults
    (`mapping.decoys.shuffle_window_aa`, `config.seed`); `shuffle` is the de novo decoy type
    (`mapping.decoys.denovo_input`).
11. **Match:** `hits = match_exact(top_as_written, index, {"shuffle": decoys},
    {p: "shuffle" for p in top_as_written}, min_len=7, max_len=30)`.
12. A peptide **matches** when `hits` has ≥ 1 row with `search == "target"` (any locus, any
    replicon, either strand; several loci count once). "Exactly" means the matcher's rules only:
    I ≡ L, no other mass equivalence (no N ↔ GG, Q ↔ AG, K ↔ Q), no stop read-through, `X`
    codons and query letters outside the residue alphabet (`B J O U X Z`) never match.
13. `share = n_match / n_top`. (b) holds when `share < 0.20`. With `n_top = 0`, `share` is `null`
    and (b) is not evaluated (only (a) decides).
14. Reported with it: the Wilson 95 % interval of `share` (R4; the decision uses the point
    estimate), the shuffled-decoy match share `r = n_decoy_match / n_top` (expected about 1 %:
    10.1 M residue positions in six frames over 20⁷ = 1.28 × 10⁹ seven-mers gives ≈ 0.8 % for a
    uniform 7-mer, far less for longer ones), matches per replicon, the share of the top 1,000
    that matches `contaminants.fasta` (`match_contaminant_peptides`; diagnostic: contaminants
    count as non-matches in (b)), and the length distribution.

Why the thresholds separate a working from a failed run: a working de novo run on a matched
genome puts most of its top-scored peptides on the genome; random or garbled 7–30-mers match at
about `r`. 20 % is ≥ 10 × the chance rate and well below what a working tool gives.

### `validation/pilot/metrics/a3_feasibility.json`

`schema_version` 1; `run`; `cascadia_version`; `weights_sha256`; `route` (B0–B2) and converter
version/digest; `mzml_sha256`; `genome_sha256`; `cutoff`, `cutoff_source`
(`documented|cli_default|none|override`); `n_rows`, `n_len_7_30`, `n_pass`, `n_top`, `n_match`,
`share`, `share_ci95` `[low, high]`, `decoy_share`, `contaminant_share`,
`matches_by_replicon` `{contig: n}`, `rule` `{min_peptides: 500, min_share: 0.20, top_n: 1000}`,
`fails` `{too_few_peptides: bool, low_genome_match: bool | null}`, `verdict` `pass|fail`,
`job_ids`, `date`. Keys sorted, 2-space indent, trailing newline (byte-stable for a fixed input).

Exit codes: **0** pass; **3** fail (the rule holds: the fallback decision); **1** error (unreadable
result, failed mzML check, missing column); **2** precondition (missing file, unknown option).

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | 2,400 distinct passing 7–30-mers; of the top 1,000, 640 match the target six-frame, 9 the shuffled decoy | `n_pass` 2400, `share` 0.64, CI ≈ [0.61, 0.67], `decoy_share` 0.009, `verdict` pass, exit 0; the 1,000-row sample is recorded |
| 2 | 420 distinct passing peptides, 300 of them match | `fails.too_few_peptides` true, `n_top` 420, `share` 0.714 still reported, `verdict` fail, exit 3 |
| 3 | 3,000 passing; 150 of the top 1,000 match | `low_genome_match` true, `share` 0.15, fail, exit 3 |
| 4 | exactly 500 passing; exactly 100 of them match (`n_top` 500, `share` 0.20) | both rules use "fewer than": pass, exit 0 |
| 5 | 1,200 passing; ranks 999–1,002 share one score, sequences `AAAK…`, `AAAL…`, `ACDE…`, `ADDE…` (folded) | ranks 999, 1,000 = the first two in folded-sequence order; the result is identical on reruns and on shuffled input rows |
| 6 | rows `PEPTIDEK` (score 0.91) and `PEPTLDEK` (0.95); `M[+15.995]PEPTIDEK`, `n[+42.011]SAMPLER` | one distinct `PEPTLDEK`-folded entry with score 0.95, kept as written `PEPTIDEK` (first in file order); modified rows strip to `MPEPTIDEK`, `SAMPLER` |
| 7 | a peptide whose codons are the last 12 and first 12 nt of the circular chromosome (8 aa) | matches (origin wrap), counted once, `matches_by_replicon["CP171242.1"]` += 1 |
| 8 | the reverse complement of a planted 30-nt stretch translated in phase 2 of the minus strand | matches |
| 9 | a 9-mer whose genomic stretch contains a TAA codon in frame (the peptide reads through it as `Q`) | no target match (the frame has `*` there) |
| 10 | a result with only the header line | `n_rows` 0, `n_pass` 0, `n_top` 0, `share` null, `too_few_peptides` true, `low_genome_match` null, fail, exit 3 |

Examples 7–9 run on a synthetic genome built at test time (a short circular contig plus a linear
contig), not on Nissle.

## Edge cases

- **Circular vs linear.** The origin wrap applies only to circular replicons (all four Nissle
  records are tagged circular). On a linear contig, the same origin-spanning 8-mer of example 7
  does not match (test on the synthetic linear contig).
- **Minus strand:** covered by `build_sixframe` (frames `-0, -1, -2`); the peptide must read in
  the minus strand's 5' → 3' direction.
- **Empty input:** header only (example 10); a zero-byte or missing result file is exit 2
  ("no Cascadia result for WT_1"), not a verdict; a job that died (non-zero exit, TIMEOUT twice)
  is exit 1 and R9, never a fail verdict. A verdict needs a completed run with a passing mzML
  check.
- **Ties:** in score at the top-1,000 boundary (example 5); between I/L variants of one peptide
  (example 6: maximum score, first as-written form).
- **Fewer than 1,000 but ≥ 500 passing:** `n_top = n_pass`; (b) on all of them.
- **Letters outside the residue alphabet** (`X`, `B`, `Z`, `U` in a prediction): stay in the
  counts, never match (reported as `n_nonstandard`).
- **N-terminal Met from a non-ATG start** (the genome has GTG/TTG translated `V`/`L` in the six
  frames): such a peptide does not match; rare, not corrected (conservative, it lowers `share`).
- **Contaminant peptides** (keratins, trypsin): count as non-matches; reported separately.
- **Duplicated sequences across replicons** (e.g. IS elements): one match per peptide.
- **Cascadia writes several candidates per spectrum:** only the rank-1 prediction per spectrum is
  read (if the format has a rank column), then the distinct step applies.
- **Score column missing or non-numeric:** exit 1 naming the column; no verdict.
- **A diaPASEF variant** (synchro/slice/midia-PASEF, M10-01 Q3): converted as diaPASEF; the
  window count in the mzML check comes from `acquisition.json`.
- **`dataset_mode` not diaPASEF:** if run 1 is ddaPASEF, Cascadia (a DIA tool) is not the right
  tool: the issue stops with a lead question, no verdict.
- **mzML too large for the disk budget:** `mzml_pilot.sbatch` checks free space ≥ 4 × the `.d`
  size before converting (exit 2 otherwise).
- **Weights sha256 differs from the committed one** (re-published file): exit 2; never run on an
  unverified file (R6).

## Acceptance criteria

Unit and stub tests (CI, every commit; no GPU, no Cascadia, no converter, no network):

1. `a3_feasibility.py` on synthetic Cascadia-format fixtures (written by the test-author under
   `tests/data/pilot/a3/`, a few kB) gives exactly the counts, verdicts and exit codes of examples
   1–6 and 10, including both boundary values of example 4 (499 → fail, 500 → pass;
   share 0.199 → fail, 0.200 → pass).
2. The genome check on a synthetic genome reproduces examples 7–9 and the linear-contig case,
   and its result equals a direct call of `match_exact` on the same index (no reimplementation:
   the module imports `build_sixframe`, `build_decoys` and `match_exact`).
3. The JSON has every field above, sorts its keys, and is byte-identical across two runs and
   across a row-shuffled input.
4. `--mzml-check` on small synthetic mzML files (a handful of spectra) passes and fails as the
   rules say (no MS2; fewer windows than `acquisition.json`; RT coverage below 90 %); a failing
   check exits 1.
5. `fetch_a3_assets.sbatch` against stubbed `curl` and `apptainer`: requests only the URLs in
   `a3_assets.tsv`, only `https://`, no TLS-off flag; a sha256 mismatch exits 1 and leaves no
   file at the final path; the image is pulled by `@sha256:` digest; final files are mode 0444
   under `$GIAP_ROOT`.
6. `mzml_pilot.sbatch` against a stubbed `apptainer`: the converter command carries the route's
   options and writes only under `$GIAP_DATA_DIR/pxd078459/mzml/`; too little disk (stubbed `df`)
   exits 2.
7. `denovo_pilot.sbatch` against stubbed `micromamba`, `nvidia-smi`: `RUN=WT_1` passes one input;
   a weights sha256 mismatch exits 2 before Cascadia starts; outputs go to
   `$GIAP_ROOT/pilot/A3_cascadia/WT_1/`. **Every M3-24 test passes unchanged.**
8. `build_envs.sh cascadia` resolves the `PILOT_TOOLS` entry (stubbed `micromamba`); `build_envs.sh
   all` does not list `cascadia`; `required_tools` never returns it.
9. `scripts/ci/check_imports.py` passes: nothing in `src/` or `validation/` imports `cascadia`.
10. No committed file under `validation/` or `tests/data/recorded/cascadia/` exceeds 5 MB;
    `git ls-files` lists no `.ckpt`, `.sif`, `.mzML`, `.d`, `.tdf`, `.tdf_bin`.

Evidence on FarmShare (PR description and PROGRESS, not CI):

11. The env builds (job ID; `cascadia` version command prints `0.0.7`; locks committed); the
    0.0.7 licence and the weights' licence are stated in the env header and the README.
12. The weights and the converter image are fetched from official sources with sha256/digest in
    `a3_assets.tsv`, or a Blocked row names what is missing.
13. ADR-M10-B is written with the chosen option; for B1/B2 Khoa's acceptance is recorded in
    PROGRESS › Answered by Khoa before the first conversion.
14. The plumbing and full smoke jobs' IDs, wall times, MaxRSS, GPU and node are in PROGRESS;
    `validation/pilot/metrics/a3_feasibility.json` holds the verdict; on a pass, the 1,000-row
    sample and its `PROVENANCE.md` are committed; on a fail, the lead opens M10-11a.

Link to §15.6: "GIAP runs end to end on A1, A2 and A3" needs an A3 input; this issue decides
whether it comes from Cascadia (M10-09/10) or from the Casanovo fallback (M10-11a/11b).

## Test plan (stubs only)

- **Unit** (`tests/unit/test_m10_08_a3_feasibility.py`): stripping, length filter, cut-off
  handling (`documented`, `none`, `override`), I/L-folded distinct with max score, ranking ties,
  both rules and their boundaries, JSON layout and determinism, exit codes; synthetic genome for
  circular/linear, minus strand, stop codons, multiple loci.
- **Unit** (`tests/unit/test_m10_08_mzml_check.py`): streamed mzML fixtures for the sanity check.
- **Script stubs** (`tests/unit/test_m10_08_scripts.py`, the M3-24 `tests/helpers/pilot_kit.py`
  sandbox pattern): fake `curl`, `apptainer`, `micromamba`, `nvidia-smi`, `df`, `sha256sum` on
  `PATH`; assert requested URLs, commands, output paths, exit codes, file modes.
- **Registry**: `PILOT_TOOLS` entry and `build_envs.sh` resolution with a stubbed `micromamba`.
- Not tested in CI: the GPU run, the conversion, the real verdict (evidence only). No
  `needs_gpu` test is added: a GPU test would only re-run the smoke job.

## Risks (what may be unreachable or fail)

| Risk | Signal | Plan |
|---|---|---|
| Weights not published at a §3.2 source (no GitHub release exists; may be Zenodo, a lab server or a cloud drive) | README of 0.0.7 names a host outside §3.2 | an ADR adding the host (§14.8 item 1, **Khoa**); Blocked row until then; never a re-hosted copy (R6) |
| The PyPI `cascadia` 0.0.7 is not the Noble Lab's package (name collision) | `project_urls` in the PyPI metadata | stop; Blocked; ask the lead (no GitHub tag exists to build from) |
| Cascadia 0.0.7's licence or the weights' terms restrict use | LICENSE in the 0.0.7 distribution, weights page | record; non-commercial → documented as for geNOMAD; anything stricter → Khoa |
| No PyTorch/CUDA build of the pinned deps solves on conda-forge, or CUDA mismatch on L40S | solve fails; `torch.cuda.is_available()` false in the plumbing job | pip PyTorch wheels with bundled CUDA inside the env; R9 after 3 attempts |
| Vendor licence needed for `.d → mzML` (B1, B2) | ADR-M10-B | **Khoa** accepts or not; Blocked row; B4 if he declines |
| ProteoWizard under Wine does not run in Apptainer (read-only image, `WINEPREFIX`) | conversion job errors | `--writable-tmpfs`, sandbox build under `$GIAP_ROOT`, binds (§17); then B2 (R9) |
| Docker Hub pull rate limit | HTTP 429 | retry later (§17) |
| Cascadia cannot use diaPASEF spectra (ion-mobility-collapsed windows differ from its training data) | the L6 rule | that is what the rule measures: fail → fallback |
| Full run 1 exceeds 24 h | TIMEOUT | one resubmission at 47 h (partition limit 2 d); then R9 and a lead question (slice-based verdict?) |
| mzML far larger than expected | disk check | peak picking and zlib already on; clean after the verdict (§3.9) |
| GPU queue waits | `PENDING (Resources)` | wait; never above 1 GPU for this issue |
| M10-01 blocked (data login-gated) | M10-01 Blocked row | steps 1–2 and the evaluator still land; the smoke run waits |

## Estimate

| Part | Code lines | Test lines |
|---|---|---|
| `cascadia.yml`, `PILOT_TOOLS` entry, `build_envs.sh` lookup | 45 | 50 |
| `fetch_a3_assets.sbatch`, `a3_assets.tsv` | 80 | 90 |
| `mzml_pilot.sbatch` | 55 | 50 |
| `denovo_pilot.sbatch` changes | 30 | 40 |
| `a3_feasibility.py` (reader, rule, genome check, mzML check, JSON) | 120 | 170 |
| Total | ≈ 330 | ≈ 400 |

Wall time: env and assets about 1 h of jobs; conversion and the full Cascadia run are
unmeasured (hours); agent time about 1 day plus queue waits and Khoa's answer on ADR-M10-B.

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| Q1 | L6 says "Cascadia's own score cut-off". If 0.0.7 documents no threshold, what counts as passing? | the documented threshold, else the CLI default; if neither exists, every reported row passes (`cutoff = "none"`) and rule (b) carries the weight; never a number GIAP invents |
| Q2 | "Distinct" peptides: as written, or I/L-folded? | I/L-folded (Cascadia cannot tell I from L; the matcher folds too) |
| Q3 | The committed Nissle FASTA has a fourth record, **CP171243.1 plasmid p15A-ILA-dapA**, not in §15.1's "chromosome, pMUT1 and pMUT2". Search it in the six-frame check? (It also bears on A2, TE1 and metric 8.) | yes, all four records, as `giap run` does; matches per replicon reported; the lead checks the assembly record and tells M10-05/07 |
| Q4 | Register Cascadia in a separate `PILOT_TOOLS` mapping (so runs and `build_envs.sh all` never see it), or in `TOOLS` with module `opt_in`? | `PILOT_TOOLS` |
| Q5 | ADR-M10-B route when Cascadia has no native `.d` input: B1 (ProteoWizard vendor image) or B2 (TIMSCONVERT)? Both need Khoa's vendor-licence acceptance. | B1, B2 if B1 fails under Apptainer; Khoa asked once for both in the next batch |
| Q6 | Record the 1,000-row Cascadia sample on a fail too? | no: only on a pass (M10-09 needs it then; on a fail it would be dead weight in `tests/data/`) |
| Q7 | The issue is ≈ 330 code lines, not 150. Split it into 08a (env, registry, assets, route, ADR-M10-B) and 08b (evaluator, smoke run, verdict)? | one PR (under the ~400-line cap); split only if review asks |
| Q8 | L6 is evaluated on run 1 only. Run 2 if run 1 is borderline (for example `share` in [0.15, 0.25])? | no: run 1 decides, as L6 says; M10-10 runs all three anyway on a pass |
| Q9 | §5.4 calls the Cascadia result "mzTab-like". If 0.0.7 writes another format (for example a BiblioSpec `.ssl` table), is an ADR needed? | no ADR: §5.4 says "its result file"; M10-08 records the real format, M10-09 maps it |
| Q10 | Run the 10-minute RT-slice plumbing job before the full run? | yes (one cheap GPU job catches env, CUDA and format errors before a 24 h job) |
| Q11 | Where the weights are published is unknown (no GitHub release). If the 0.0.7 docs point to a host outside §3.2, ADR + Khoa? | yes (§14.8 item 1); Blocked row meanwhile |

## Needs Khoa (batched by the lead)

1. **Conditional:** acceptance of the vendor licences bundled in the converter (ADR-M10-B,
   options B1/B2): the ProteoWizard vendor-reader image ("i-agree-to-the-vendor-licenses",
   Bruker TDF-SDK among them) or TIMSCONVERT's bundled Bruker `libtimsdata.so`.
2. **Conditional:** a new download source if the Cascadia weights are hosted outside §3.2, or
   their terms if they are more restrictive than non-commercial.

## Lead decisions (2026-10-23)
1. **Q1–Q11:** the proposed defaults are taken, all reversible:
   - the cut-off is the documented or CLI-default one, else `none`;
   - peptides are counted after folding I into L;
   - all four records are in the six-frame check, with matches reported per replicon;
   - the tool goes in `PILOT_TOOLS`;
   - the conversion route is B1, then B2;
   - the sample is recorded on a pass only;
   - one PR;
   - run 1 decides;
   - no ADR for the result format;
   - the slice runs first;
   - weights hosted outside §3.2 need an ADR and Khoa, with a Blocked entry meanwhile.
2. **Sequencing:** M10-08 waits behind pass 1 (M10-01 to M10-05) and Khoa's answer on the vendor licences (B1/B2). Only the environment build and the asset discovery (where the weights live, and the 0.0.7 licence) may start earlier. No raw data or vendor image is used before then.
3. **The fourth Nissle record:** CP171243.1, plasmid p15A-ILA-dapA, is in `tests/data/genomes/nissle1917.fna.gz`, but §15.1 does not name it. That affects A2, TE1 and metric 8. It goes to Khoa in the next batch, and M10-05 and M10-07 report per replicon meanwhile.
