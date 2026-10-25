# M10-08: A3 feasibility: Cascadia 0.0.7 on the diaPASEF pilot runs

Milestone M10 · **pass 2** (may run alongside pass 1: it shares no code) · depends on: M10-01
(merged: the three `.d` runs under `$GIAP_DATA_DIR/pxd078459/`, `validation/pilot/inputs/acquisition.json`
with `dataset_mode = diaPASEF`) · estimate ≈ 455 code lines + ≈ 450 test lines, proposed as two
PRs (08a ≈ 270, 08b ≈ 185; question Q1) · not `stats` (one Wilson interval, reported only) ·
roles: implementer (env files, registry, `scripts/`), validation (`validation/pilot/`),
test-author (`tests/`) · **Khoa: needed** for ADR-M10-B: accepting the Bruker TDF-SDK licence
on the record, and adding `drive.google.com` (where the Cascadia weights live) as a download
source.

This spec replaces the earlier draft on `m10/a3-feasibility` (16b08360, 8c57871f). It keeps that
draft's failure rule (index L6, lead decision 1) and its genome check. Four of that draft's lead
decisions no longer fit what the official sources say: Q1 (cut-off), Q5 (route order), Q6
(committed sample), Q10 (slice first). Questions Q2–Q5 below raise them again.

## Spec sections implemented

- §15.1 arm A3: "Cascadia 0.0.7 on the same runs (exploratory), on the `gpu` partition (1 ×
  L40S). If DIA de novo fails on timsTOF data: Casanovo 5.2.1 on 3 UTI89 DDA runs (PXD053727;
  genome GCA_000013265.1)". This issue decides **whether** the fallback is needed. The fallback
  itself is M10-11a/11b.
- §15.2 step 6 (`denovo_pilot.sbatch`, `--partition=gpu --gres=gpu:1`), and the Cascadia half
  of step 7, as amended by Khoa's 2026-10-23 rule: recorded samples stay on FarmShare (M10-02
  lead decision 4).
- §17 row "Cascadia fails on timsTOF data | A3 empty or degenerate | Casanovo on 3 UTI89 DDA
  runs". Index L6 (lead decision 1) makes "empty or degenerate" measurable.
- §3.1 `gpu` partition. Checked 2026-10-24 with `sinfo`: `oat-[01-06]`, `gpu:4` per node, limit
  2-00:00:00; QOS `gpu` `MaxTRESPU gres/gpu=4`. Also §3.2 (sources), §3.8 (job header), §3.9
  ("de novo outputs: measure"), R5, R6, R9, R11, R14, and §14.8 item 1 (a new dependency or
  source needs an ADR and Khoa).
- §5.4 row "Cascadia 0.0.x | its result file (mzTab-like)". The real format is the `.ssl` table
  (below). This issue records the format; the converter is M10-09.
- §7.2.1–7.2.2: the genome check reuses `giap.mapping.sixframe.build_sixframe`,
  `giap.mapping.decoys.build_decoys` and `giap.mapping.exact.match_exact`, unchanged.
- Khoa, 2026-10-24: a Bruker-licensed conversion route is allowed "as long as it's free and
  legal".

## What the official sources say (read 2026-10-24 UTC from `iron-13`; nothing executed)

### Sources consulted

| # | Source (official) | Reached | Used for |
|---|---|---|---|
| S1 | https://github.com/Noble-Lab/cascadia (API: repo, tree, tags, releases, issues #1, #5, #6, #12, #15, #16, PRs #9, #14, #18, #19) | yes | licence (Apache-2.0), no tags, no releases, open issues |
| S2 | https://pypi.org/pypi/cascadia/json and the 0.0.7 sdist `cascadia-0.0.7.tar.gz` (sha256 `6e9fcc50…ab347`; read, not installed) | yes | the code of the pinned version: CLI, reader, writer, LICENSE |
| S3 | https://cascadia.readthedocs.io/en/latest/ (`usage.html`), plus its source `docs/source/{getting_started,usage,file_formats}.md` in S1 | yes | input format, CLI, output format, weights link |
| S4 | Sanders et al., *Nat Methods* 22, 1447–1453 (2025), https://doi.org/10.1038/s41592-025-02718-y | **partly**: the full text is paywalled. Only the abstract and the Data and Code availability statements could be read | weights licence; evaluation datasets |
| S5 | bioRxiv 10.1101/2024.06.03.597251 (CC BY), JATS XML via https://api.biorxiv.org/details/biorxiv/10.1101/2024.06.03.597251 | yes. The first try got HTTP 429; the one later retry worked (R6) | instruments evaluated |
| S6 | ProteomeXchange XML of S4's datasets PXD056793, PXD053291, PXD042704 (`proteomecentral.proteomexchange.org`, a §3.2 host) | yes | instruments of the evaluation data |
| S7 | Google Drive folder linked from S3, https://drive.google.com/drive/folders/1UTrZIrCdUqYqscbqga_KdX8kc8ZjMMfr | listing only. Files `cascadia.ckpt`, `cascadia_astral_tuned.ckpt`, `cascadia_pretrained.ckpt`, `demo.mzML`; no sizes, no checksums, no licence file | weights source |
| S8 | https://proteowizard.sourceforge.io/licenses.html, `core_license.html`, `download.html`, `tools/msconvert.html`, `js/main.js` | yes | ProteoWizard licence and download routes |
| S9 | Vendor EULAs that S8 loads from https://github.com/ProteoWizard/pwiz/tree/master/pwiz/data/vendor_readers (`Bruker/EULA.Bruker.txt` and the Agilent, Sciex, Shimadzu, Thermo, Waters files) | yes | **the TDF-SDK licence text** (quoted below) |
| S10 | Docker Hub `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses` (S8's download link points to the old `chambm/…` name, whose page says "This Docker image has moved") | yes | image digest, licensing statement |
| S11 | `pwiz/data/vendor_readers/Bruker/{TdfData.cs,SpectrumList_Bruker.cs}`, `pwiz_tools/Commandline/msconvert-help.txt` (S9 repository) | yes | how msconvert writes diaPASEF spectra |
| S12 | https://github.com/gtluu/timsconvert (API, README, docs sources, `requirements.txt`, `setup.py`, `bin/cmd.py`, `timsconvert/{arguments,parse,write}.py` at tag `v2.0.0`) | yes | licence, install route, diaPASEF export |
| S13 | https://github.com/gtluu/pyTDFSDK at commit `4808da81…` (tree, README, `pyTDFSDK/classes.py`) | yes | the bundled `TDF-SDK/libtimsdata.so` |
| S14 | Bruker's own TDF-SDK download page, which S9 names (`https://www.bruker.com/service/support-upgrades/software-downloads/mass-spectrometry.html`) | **no**: it redirects to a login page (`/en.login.html?targetUrl=/protected/…`) | the TDF-SDK text was read from S9 instead. pyTDFSDK's bundled `EULA TDF-SDK.pdf` was not opened |
| S15 | https://www.unimod.org/unimod_help.html; https://www.unimod.org/xml/unimod_tables.xml (HTTPS answers 206 to a range request) | yes | Unimod terms |
| S16 | psims 1.4.0: PyPI JSON and the wheel's file list (read with `zipfile`, not installed); conda-forge `psims` 1.4.0 (Apache-2.0) | yes | a vendored `unimod_tables.xml.gz` from a §3.2 source |
| S17 | https://github.com/levitsky/pyteomics `pyteomics/mass/unimod.py`; PyPI pyteomics (5.0.1, Apache-2.0) | yes | why Cascadia reaches the network |
| S18 | quay.io `protio/cascadia` tags; https://github.com/mriffle/cascadia-docker (Dockerfile); https://github.com/mriffle/nf-skyline-dia-ms (`modules/cascadia.nf`, `container_images.config`). This is the workflow S4 names | yes | a second copy of the weights (question Q8) |

### Cascadia 0.0.7: what the pinned code does (S2, S3)

| Fact | Evidence | Consequence |
|---|---|---|
| Licence **Apache-2.0** | `LICENSE` in the 0.0.7 sdist; PyPI classifier; GitHub licence API | free; may be used and redistributed. GIAP does not redistribute it |
| Source: **PyPI only**. No git tag and no GitHub release; `pyproject.toml` on `main` still says 0.0.5 | S1 tags `[]`, releases `[]`; S2 releases 0.0.1–0.0.7 (0.0.7 uploaded 2025-02-14; wheel sha256 `50b52043…f092a3`, sdist `6e9fcc50…ab347`) | pin `cascadia==0.0.7` with the wheel's sha256 in the pip lock. PyPI `project_urls.Homepage` = `https://github.com/Noble-Lab/cascadia` (the provenance check of the earlier draft is done) |
| Requires `torch>=2.0,<2.1`, `lightning>=2.0,<2.1`, `pytorch-lightning>=1.9,<2.0`, `numpy<2.0`, `pyteomics>=4.6`, `lxml>=4.9.1`, `sqlalchemy>=2.0`, `rdkit`, `numba`, … | PyPI `requires_dist` of 0.0.7 | Python 3.10 (S3) with PyTorch 2.0.x and CUDA 11.8 (Ada / L40S is sm_89, which CUDA 11.8 supports). Unmerged PR #18 says the `lightning` + `pytorch-lightning` pair breaks fresh pip installs: risk table |
| **Input: mzML only** (`pyteomics.mzml.read`). No `.d` reader | `augment.py`; S3 `file_formats.md`: "Cascadia takes as input spectra in mzML format … converted … using msConvert … peakPicking filter at both the MS1 and MS2 levels using the vendor algorithm" | option B0 (native `.d`) does not exist: a conversion route is required |
| It reads, per MS2 spectrum, **one** `isolationWindow` (target m/z, lower and upper offsets) and the scan start time, and per MS1 spectrum the peaks. **Ion mobility is never read** | `augment.get_centers`, `extract_spectra` | the mzML must hold **one spectrum per diaPASEF isolation window per frame**, with its own window. A whole frame merged into one spectrum would break the precursor model (mzML check below) |
| CLI: `cascadia sequence SPECTRUM_FILE MODEL [-o OUTPREFIX] [-t 0.8] [-b 32] [-w 2] [-c 4] [-p mskb]`; **one** spectrum file, model positional | `cascadia.py` `sequence()` argparse (S3 `usage.md` shows `-m MODEL` and a `-f` flag that 0.0.7 does not have; the code wins) | the current `denovo_pilot.sbatch` call (`--model … --num_workers … --output … RUN…`) **fails** in argparse. One call per run |
| Writes a scratch directory `cascadia_<mm-dd-HH:MM:SS>/` in the **current directory**, with a text `temp.asf` and an `index.hdf5` | `sequence()` | `cd` into a per-run directory on the node's `$TMPDIR` (§3.1 local SSD) before the call |
| DataLoader `num_workers=4`, hard-coded | `sequence()` | `--cpus-per-task=8` is enough |
| **Score cut-off: `-t/--score_threshold`, default 0.8**. A row is written only when `exp(log-confidence) > threshold` (**strict**) | `utils.write_results` | this is "Cascadia's own score cut-off" in L6. It is a CLI default, not something GIAP picks |
| Output `<OUTPREFIX>.ssl`, tab-separated, columns `file, scan, charge, sequence, score-type, score, retention-time, start-time, end-time`; `score-type` = `UNKNOWN`; `scan` = the prediction's index, not a scan number | `utils.write_results`; S3 | the format M10-09 maps. One row per augmented spectrum × charge; no rank column |
| Sequences come out **I → L folded** (`replace_isoleucine_with_leucine=True`), with modifications as mass shifts (`M[+15.99]`; an N-terminal shift as `[+42.01]-…`) | tokenizer call; `convert_mods_to_shifts` | stripping rule below. I ≡ L already holds |
| Retention-time unit decided by a heuristic: minutes when the largest RT > 500, else seconds | `write_results` | **these runs are 539.9 s long** (`acquisition.json`), close to that limit. M10-09 must not assume a unit. The evaluator ignores RT |
| **At the end of every run it downloads `http://www.unimod.org/xml/unimod_tables.xml`** (pyteomics `mass.unimod.Unimod()` with no path). That is plain HTTP to a host outside §3.2, and with lxml ≥ 5.4 it fails ("Attempt to load network entity", issues #6 and #12, both open) | `write_results`; S17 lines 61, 772 | R6 forbids the download. Fix: run Cascadia through an offline launcher (deliverable 4) that points pyteomics at a local copy of the table **before** Cascadia starts. The copy comes from psims 1.4.0, which vendors `unimod_tables.xml.gz` and is on conda-forge, a §3.2 source (S16) |

### Does Cascadia support diaPASEF (timsTOF)?

**Not documented, not excluded, never tested by its authors.**

- The 0.0.7 docs and code name no instrument, no ion mobility, no timsTOF and no diaPASEF (S2,
  S3). An issue search on the repository for `timsTOF OR diaPASEF OR bruker OR tdf` returns
  one hit (#16), which is about the training file format, not timsTOF.
- The preprint (S5) evaluates only Orbitrap instruments: Fusion Lumos, the DeepNovo-DIA test
  set, and Orbitrap Astral. "timsTOF", "PASEF" and "mobility" appear 0 times in it.
- The published paper's datasets (S4 Data availability; S6) are CPTAC (Orbitrap), MSV000082368,
  PXD053291 (Astral, Q Exactive HF-X), PXD042704 (Astral, Fusion Lumos) and PXD056793 (Astral,
  Exploris 480, Fusion Lumos). None is timsTOF. The paper's full text was not reachable (S4).
- The model is trained on Orbitrap DIA. diaPASEF MS2 spectra come from a TOF at lower
  resolution, with mobility-resolved window groups. Cascadia's reader accepts them only as plain
  isolation-window spectra (the table above).

**Consequence (rule E0, stated explicitly).** Nothing says Cascadia **cannot** take diaPASEF data
once it is converted to per-window mzML. So the "fails on timsTOF" rule **does not fire on
documentation alone**, and the one-run smoke test goes ahead. It is cheap: one 9-minute run,
about 1 GPU-hour.

The rule **fires early**, before the evaluator runs, in exactly these cases. Each one ends
M10-08 with `verdict = fail`, and A3 moves to the **M10-11a/11b Casanovo fallback** (ADR-M10-C
goes to Khoa first):

| Code | Early-fire case | Evidence required |
|---|---|---|
| E1 `route_unavailable` | No `.d → mzML` route is clearly free and legal (the licence rule below), or Khoa declines the Bruker TDF-SDK terms | Khoa's answer in PROGRESS, or the licence finding |
| E2 `structurally_incompatible` | Both routes (B2, then B1) give mzML that fails the structural check on the points Cascadia needs (one isolation window per MS2 spectrum, MS1 present), after R9 (3 attempts per route) | the two `mzml_check.json` files |
| E3 `tool_incompatible` | Cascadia 0.0.7 exits non-zero on an mzML that **passed** the check, for a reason in its data handling, after R9. Environment, CUDA, Unimod, memory and time-limit failures do **not** count | job logs; the lead confirms E3 (question Q6) |

What does **not** fire the rule: missing weights (a Blocked row: Khoa, ADR-M10-B), a failed env
build (R9 escalation), GPU queue waits, or a missing M10-01 input. None of these is evidence
about timsTOF data.

## Licence findings

### The two conversion routes

Bruker `.d` stores peaks in `analysis.tdf_bin`. Both candidate routes decode it with Bruker's
proprietary **TDF-SDK** (`libtimsdata.so` / `timsdata.dll`). Its licence, read from S9 because
S14 is login-gated, says (quoted):

> "These license terms … shall be the legal basis for any noncommercial or commercial usage of
> Licensor's Software. … You may only use the Software as a whole, if you have downloaded it from
> Licensor's Webserver … However, certain individual Software Libraries might be conveyed to you
> by others. In both cases, you may only legally use Software or any of the Software Libraries
> conveyed to you after having accepted this Agreement as binding between Licensor and you."
>
> "2.1 You may install and use the Software on any of your computing devices."
>
> "3.2 You may not distribute the Software to others, breach of which terminates the rights
> granted under this Agreement. … However, you may distribute Own Software which requires the
> usage of or links to Software Libraries if you comply with the requirements set down in
> Section 4."
>
> "6.1 As Licensor delivers the software for free, …"
>
> §1: the Software is not "designed to perform any queries or operations which might require
> reliable results and/or which might be the basis for business and/or medical decisions."

| Question | B2: TIMSCONVERT 2.0.0 + pyTDFSDK (bundles `libtimsdata.so`) | B1: ProteoWizard `msconvert`, Docker Hub image `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses:4.0.26282-8321319` |
|---|---|---|
| Free of charge? | **Yes.** TIMSCONVERT and pyTDFSDK are Apache-2.0 (S12, S13); TDF-SDK: "delivers the software for free" (§6.1) | **Yes.** S10: "In general, ProteoWizard is free for both commercial and non-commercial use, but you may not charge a fee for use of ProteoWizard (with vendor support enabled)". The core is Apache-2.0 plus an addendum binding the user to every vendor library's terms (S8 `core_license.html`) |
| Licences we must accept | **One**: Bruker's TDF-SDK agreement. pyTDFSDK conveys the library ("certain individual Software Libraries might be conveyed to you by others") and carries the required notice: "Included software components: Copyright © 2022 by Bruker Daltonics GmbH & Co. KG" (S12, S13) | **Seven**: "as the image name implies, you must agree to the vendor licenses to use this image" (S10). That is Bruker (CompassXtract: "non-commercial usage", "beta"; Baf2Sql; TDF-SDK), Agilent, Sciex, Waters (each non-commercial only), Thermo (no commercial exploitation without consent) and Shimadzu (S9) |
| May we use it on FarmShare? | **Yes**, under §2.1 "any of your computing devices", reading the Stanford allocation GIAP runs under as such (question Q7). Academic research is allowed (any non-commercial or commercial use). §1's disclaimer fits an exploratory arm: A3 is never a basis for business or medical decisions | **Yes**, for academic non-commercial use under all seven EULAs, by the same reading. It runs under Wine inside Apptainer (risk) |
| Redistribution limits | We may not redistribute the SDK (§3.2), and §5.2 bars any combined distribution with GPL software. **GIAP redistributes nothing**: the env is built on FarmShare from the upstream sources, and nothing from it is committed (R11) | No redistribution of the image or of any vendor library (Thermo §3.3: end users are bound by an agreement "that prohibits them from redistributing this Software to others"; Waters: "You may not distribute, assign, rent, sublicense, “timeshare,” or transfer"). GIAP keeps the image under `$GIAP_ROOT/images/` and never commits it |
| Within Khoa's 2026-10-24 permission ("Bruker-licensed … free and legal")? | **Yes, exactly**: the only proprietary component is Bruker's | **Goes beyond it**: it also accepts six licences from other vendors, for libraries this pilot never loads |
| Source (§3.2, R6) | GitHub: release tag `v2.0.0` (commit `3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b`), which pins pyTDFSDK at commit `4808da81d3e37f82da96f392dbb446c3733edacb` (one commit behind tag `v0.4.1`). Not on PyPI or bioconda (both 404). A pinned GitHub commit has a precedent in UniOP (§17 "clone fails at the pinned commit"); ADR-M10-B records it | Docker Hub, by digest `sha256:34ed787cd4c3182c2d81adcbf7514031e594609ca460abc5d3d114dc2334dae8` (1.89 GB compressed, 2026-10-09). The cleanest §3.2 source |
| diaPASEF output fits Cascadia? | **Yes, in its code**: `parse.py` handles `ScanMode 9, MsMsType 9` frames and writes one spectrum per `DiaFrameMsMsWindows` row, with `isolation window target / lower / upper` (S12, S13); `--exclude_mobility` sums the window's scans into one centroided spectrum | **Yes, in its code**: with `--combineIonMobilitySpectra`, plain diaPASEF gives one combined spectrum per frame and window (`TdfData.cs` `BuildSpectrumIndex`; the whole-frame mode is only for diagonal PASEF or when requested) (S11). It is the route Cascadia's docs name |
| Runs natively on Linux? | yes (TDF-SDK `.so`) | under Wine (the image's documented use). ProteoWizard also offers a Linux .NET build that reads Bruker data, but only from its TeamCity server, which is not a §3.2 source (S8) |

**Route chosen: B2 (TIMSCONVERT 2.0.0 with pyTDFSDK's bundled TDF-SDK).** It is the route that is
**clearly** free and legal within Khoa's permission. It is free (Apache-2.0 plus a free SDK), its
only proprietary licence is Bruker's, it is used and not redistributed, it runs natively, and its
code exports diaPASEF exactly as Cascadia needs. It needs one Khoa action: accepting the TDF-SDK
agreement on the record (ADR-M10-B). The SDK says use is legal only "after having accepted this
Agreement".

**Second choice: B1**, only if B2 fails after R9 (environment, conversion or mzML check), and only
after Khoa also accepts the six non-Bruker EULAs that come with the image (question Q2).

**If neither is clearly free and legal** (Khoa declines the TDF-SDK terms, or a licence reading
changes), rule **E1 fires**:

- `verdict = fail`, `reason = route_unavailable`;
- PROGRESS › Blocked gets a row;
- A3 moves to M10-11a.

GIAP never falls back to an unofficial converter or a GIAP-written TDF decoder (option B3,
AlphaTims plus a GIAP mzML writer, stays rejected: R1 scope, and AlphaTims ships Bruker's library
too).

### The other A3 artefacts

| Artefact | Licence (source) | Free? Use on FarmShare? Redistribution | Where it comes from |
|---|---|---|---|
| Cascadia 0.0.7 | Apache-2.0 (S1, S2) | yes; yes; allowed, but we do not redistribute | PyPI (§3.2) |
| Cascadia weights `cascadia.ckpt` | **Apache-2.0** per S4's Code availability: "The open-source Cascadia software and associated model weights are available with an Apache license via GitHub". The Drive folder (S7) itself carries no licence file and no checksum | yes; yes; allowed, but we do not redistribute. The `.ckpt` is never committed (R11 size, and no need) | **Google Drive**, the folder the docs link (S3, S7). It is **not a §3.2 source**: ADR-M10-B adds it (§14.8 item 1, Khoa). No checksum is published, so the first download is trust-on-first-use and its sha256 is recorded. The paper's "via GitHub" does not hold: there is no GitHub release (S1) |
| Unimod table (`unimod_tables.xml`) | Unimod: "a public domain database, distributed under a copyleft licence" (S15) | yes; yes; copying allowed. We do not commit it | the copy psims 1.4.0 vendors (conda-forge/PyPI, Apache-2.0 package; file md5 `2d64d1e08013808d05982885117c02da` per psims' `record.json`) (S16). No new host. Fallback: `https://www.unimod.org/xml/unimod_tables.xml` (HTTPS works), which would need an ADR source row (question Q9) |
| pyteomics (Cascadia dependency) | Apache-2.0 (S17) | yes | conda-forge, pinned `<5` (question Q10) |

## Deliverables

| # | File | PR | Owner | What |
|---|---|---|---|---|
| 1 | `workflow/envs/cascadia.yml` + `locks/cascadia.{explicit,pip}.txt` | 08a | implementer | conda-forge only (`--override-channels`): `python=3.10`, `pytorch=2.0.*` CUDA 11.8 build, `numpy<2`, `pyteomics>=4.6,<5`, `psims=1.4.0`, `lxml`, `rdkit`, `numba`, `h5py`, `sqlalchemy>=2`. `pip:` block: `cascadia==0.0.7`. Header: source (PyPI 0.0.7, homepage `github.com/Noble-Lab/cascadia`), licence Apache-2.0, wheel sha256, date. Locks written by `build_envs.sh` |
| 2 | `workflow/envs/timsconvert.yml` + locks | 08a | implementer | conda-forge `python=3.11`, `pip`; `pip:` `timsconvert @ git+https://github.com/gtluu/timsconvert@3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b` (its `requirements.txt` pins pyTDFSDK at `4808da81…`). Header: route B2, the TDF-SDK licence, ADR-M10-B, "never redistribute" |
| 3 | `PILOT_TOOLS` in `giap.annotate.tools` (`cascadia`, `timsconvert`) and its lookup in `build_envs.sh` `registry()` | 08a | implementer | earlier lead decision 1 (Q4): **not** in `TOOLS`, so `required_tools`, `build_envs.sh all` and a run's tool versions never see them |
| 4 | `scripts/pilot/cascadia_offline.py` | 08a | implementer | runs **inside** `giap-cascadia` (`micromamba run -n giap-cascadia python -I scripts/pilot/cascadia_offline.py UNIMOD_XML sequence …`). It (a) refuses a `UNIMOD_XML` that is not an existing local file; (b) sets `pyteomics.mass.unimod._unimod_xml_download_url` to its absolute path; (c) hands the rest of `argv` to `cascadia.cascadia.main()`. Imports Cascadia only here, never from `src/` or `validation/` (R5 spirit; Cascadia is Apache-2.0) |
| 5 | `validation/pilot/a3_assets.tsv` | 08a, filled by the fetch | validation | one row per external A3 artefact: `asset` (`cascadia_weights`, `unimod_table`, `timsconvert`, `pytdfsdk`, `pwiz_image` if B1), `version`, `url` (official), `sha256`, `bytes`, `licence`, `path` (under `$GIAP_ROOT`), `retrieved` (UTC date). No binary is committed |
| 6 | `scripts/slurm/fetch_a3_assets.sbatch [--resolve]` | 08a | implementer | the **only** A3 script that reaches the network. Weights: HTTPS from the S7 file ID (`drive.usercontent.google.com/download?id=…&export=download&confirm=t`) **only once ADR-M10-B is accepted**, otherwise exit 2 "weights host not in §3.2: ADR-M10-B". It writes `$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt` (mode 0444; original name kept in deliverable 5). Unimod: gunzips psims' vendored file from the built env to `$GIAP_DB_DIR/cascadia/unimod_tables.xml` (0444), with no network. B1 only: `apptainer pull` by digest. `--resolve` prints the sha256 and bytes to commit into deliverable 5; later fetches must match |
| 7 | `scripts/slurm/mzml_pilot.sbatch RUN` | 08a | implementer | `.d → mzML` for **one** run. `MZML_ROUTE=timsconvert` (default) or `pwiz` (B1). Writes `$GIAP_DATA_DIR/pxd078459/mzml/<RUN>.mzML`, then runs deliverable 9's `mzml-check` |
| 8 | `scripts/slurm/denovo_pilot.sbatch` (changed) | 08b | implementer | Cascadia branch: `RUN=WT_1` (argument or env; default all three in sequence). Input `$GIAP_DATA_DIR/pxd078459/mzml/<RUN>.mzML`. Checks the weights' and Unimod table's sha256 against deliverable 5 (exit 2 on a mismatch), then `cd "$TMPDIR/<RUN>"` and runs `cascadia_offline.py UNIMOD sequence <mzML> <ckpt> -o $GIAP_ROOT/pilot/A3_cascadia/<RUN>/<RUN> -b 32` (threshold left at its default 0.8, recorded). It logs `cascadia --help` first, plus `nvidia-smi` and `python -c 'import torch; print(torch.cuda.is_available())'` (exit 1 when false). The Casanovo branch stays unchanged. **Every M3-24 assertion still holds** (no `curl`/`wget`/`git clone` in its code, `gpu:1`, `giap-denovo`, the fallback text) |
| 9 | `validation/pilot/a3_feasibility.py` | 08b (`mzml-check` in 08a) | validation | three subcommands: `mzml-check`, `evaluate`, `early-fail` (interfaces below). Imports GIAP (`giap.mapping`, `giap.annotate.bakta`, `giap.rng`, `giap.config`, `giap.stats.cutoffs.wilson_ci`). Reads the `.ssl` with the stdlib `csv` module; never imports Cascadia |
| 10 | `validation/pilot/metrics/a3_feasibility.json` | 08b | validation | the committed verdict: **aggregates only** (fields below) |
| 11 | `docs/adr/ADR-NNNN-a3-external-artefacts.md` (ADR-M10-B) | 08a | lead (drafted here) | Status: proposed · Needs Khoa: yes. Decides (1) route B2, then B1, then E1; (2) Khoa accepts the TDF-SDK agreement (and the seven B1 EULAs only if B1 is used); (3) `drive.google.com` as the source of the Cascadia weights; (4) the Unimod table from psims. The number is the next free one at merge (0051–0053 are taken on open branches) |
| 12 | `validation/pilot/README.md` "A3 de novo" section | 08b | validation | versions, licences (this spec's table), assets, the route, the commands, job IDs, sizes, wall times, the verdict |
| 13 | `tests/data/recorded/cascadia/0.0.7/PROVENANCE.md` (pass only) | 08b | validation | Cascadia version, the env lock's sha256, the weights' sha256, the route and its commit or digest, job IDs, the input `.d` archive sha256 (from `pxd078459_files.tsv`), the `.ssl` header **as defined in the Cascadia source** (not read from data). The 1,000-row sample itself goes to `$GIAP_ROOT/pilot/recorded/cascadia/0.0.7/WT_1.ssl.gz` and is **never committed** (Khoa, 2026-10-23, via M10-02 decision 4) |

Nothing else is committed: no weights, image, mzML, `.ssl`, top-1,000 list, per-peptide table, or
Unimod copy.

### Paths under `$GIAP_ROOT` (never committed)

| Item | Path | Expected size |
|---|---|---|
| Envs | `$MAMBA_ROOT_PREFIX/envs/giap-{cascadia,timsconvert}` | 5–8 GB (PyTorch + CUDA); 1–2 GB (PySide6 is pulled in by TIMSCONVERT). `micromamba clean -a` afterwards |
| Weights | `$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt` (0444) | unknown; ≈ 0.2–1.5 GB (a 512-dimension, 9-layer encoder–decoder) |
| Unimod table | `$GIAP_DB_DIR/cascadia/unimod_tables.xml` (0444) | ≈ 3 MB |
| mzML | `$GIAP_DATA_DIR/pxd078459/mzml/WT_1.mzML` | measure; ≈ 1–4 GB (1.0 GB archive, centroided, zlib) |
| Cascadia scratch | `$TMPDIR/WT_1/cascadia_*/` on the node's SSD (`temp.asf` text, `index.hdf5`), deleted by Cascadia and by the job trap | ≈ 2–5 GB, transient |
| Result | `$GIAP_ROOT/pilot/A3_cascadia/WT_1/WT_1.ssl` | small |
| Evaluator detail | `$GIAP_ROOT/pilot/A3_cascadia/WT_1/feasibility_top1000.tsv` | < 1 MB |
| Recorded sample | `$GIAP_ROOT/pilot/recorded/cascadia/0.0.7/WT_1.ssl.gz` | < 100 kB |

## Commands (B2, the route chosen)

```bash
# conversion (in mzml_pilot.sbatch; one run, one input, so TIMSCONVERT's Pool(cpu_count()-1) maps one job)
micromamba run -n giap-timsconvert timsconvert --input "$GIAP_DATA_DIR/pxd078459/WT_1.d" \
  --outdir "$GIAP_DATA_DIR/pxd078459/mzml" --mode centroid --compression zlib --exclude_mobility
# Cascadia (in denovo_pilot.sbatch)
mkdir -p "$TMPDIR/WT_1" && cd "$TMPDIR/WT_1"
micromamba run -n giap-cascadia python -I "$GIAP_REPO/scripts/pilot/cascadia_offline.py" \
  "$GIAP_DB_DIR/cascadia/unimod_tables.xml" sequence "$GIAP_DATA_DIR/pxd078459/mzml/WT_1.mzML" \
  "$GIAP_DB_DIR/cascadia/cascadia_0.0.7.ckpt" -o "$GIAP_ROOT/pilot/A3_cascadia/WT_1/WT_1" -b 32
```

B1 (only if Q2 allows it): `apptainer exec --cleanenv --writable-tmpfs -B "$GIAP_DATA_DIR" "$PWIZ_SIF"
wine msconvert WT_1.d --mzML --zlib --64 --combineIonMobilitySpectra --filter "peakPicking vendor
msLevel=1-" -o … --outfile WT_1.mzML`. This follows Cascadia's documented advice: vendor peak
picking at MS1 and MS2, otherwise defaults.

The commands actually used, with `--help` output captured in the job logs, go into the README.

## The mzML check (`a3_feasibility.py mzml-check MZML --acquisition validation/pilot/inputs/acquisition.json --run WT_1 --out J`)

It reads the file with streamed stdlib `xml.etree.ElementTree.iterparse`, decoding no binary
arrays except their lengths. It **passes** when every rule below holds; each failed rule is named
in `J`:

| # | Rule | WT_1 value it is checked against |
|---|---|---|
| M1 | MS1 spectra count = `msms_type_counts["0"]` | 734 |
| M2 | MS2 spectra count ≥ `msms_type_counts["9"]` and ≤ `msms_type_counts["9"] × dia_window_rows` | between 8,804 and 572,260 |
| M3 | every MS2 spectrum has exactly one precursor with `isolation window target m/z`, and lower and upper offsets > 0 | |
| M4 | no MS2 isolation width (lower + upper) exceeds half the acquisition m/z span. This catches whole-frame merging | span 343.83–1188.83 → width ≤ 422.5 |
| M5 | the number of distinct (target, lower, upper) triples, rounded to 0.01, lies in [`window_groups`, `dia_window_rows`] | [12, 65] |
| M6 | every spectrum is flagged centroid (`MS:1000127`) | |
| M7 | the largest scan start time ≥ 0.9 × `gradient_length_s` | ≥ 485.9 s |
| M8 | every MS2 spectrum has ≥ 1 peak (`defaultArrayLength`); the share of empty MS2 spectra is reported | |

A failed check is a **route** failure (exit 1), never a Cascadia verdict. It counts toward E2 only
as that table says.

## The failure rule L6 (lead decision 1), made exact for 0.0.7

Cascadia **fails on timsTOF** when, on the complete `.ssl` of run 1 (`WT_1`), **either**

- (a) fewer than **500** distinct de novo peptides of 7–30 aa pass the score cut-off, **or**
- (b) fewer than **20 %** of the top 1,000 of them match the Nissle six-frame translation exactly.

Otherwise it passes, and M10-09/M10-10 go ahead. A fail (or E1–E3) opens M10-11a.

### Computing (a) (`evaluate`)

1. Read every row of the `.ssl`. The header must hold `sequence` and `score` (located by name):
   otherwise exit 1. `score` must parse as a finite float: otherwise exit 1 naming the line.
2. **Strip** the sequence:
   - delete `[…]` and `(…)` with their contents;
   - delete signed masses (`[+-]\d+(\.\d+)?`) and `-`;
   - delete any character outside `A–Z`;
   - uppercase the rest.
   Then fold `I → L` (a no-op for 0.0.7 output, kept for robustness).
3. **Length:** keep 7 ≤ len ≤ 30.
4. **Cut-off:** keep `score > cutoff` (strict, as Cascadia applies it). `cutoff` = 0.8,
   `cutoff_source = "cli_default"`; `--cutoff X` overrides it (`"override"`). Rows at or below the
   cut-off are counted as `n_below_cutoff`, which is 0 for an unmodified 0.0.7 run.
5. **Distinct:** one entry per folded sequence, keeping the **maximum** score and the first
   as-written form in file order.
6. `n_pass` = the number of distinct entries. (a) holds when `n_pass < 500`.

### Computing (b) with the existing GIAP code

7. **Top:** sort by score descending, then folded sequence ascending; take the first
   `n_top = min(1000, n_pass)`.
8. **Genome:** `tests/data/genomes/nissle1917.fna.gz` and the recorded full-DB Bakta JSON, loaded
   as `validation/proxies/te_proxies.py` loads them. The contig IDs, lengths and **circular**
   topologies are those `giap run` uses. All four records are searched, CP171243.1 included
   (earlier lead decision 3), and matches are reported per replicon. The FASTA's sha256 is
   recorded.
9. **Index:** `build_sixframe(genome, replicons)` with the default `WRAP_NT` (origin-spanning
   peptides on circular replicons); both strands, three phases each, table 11.
10. **Decoys (diagnostic only):** `build_decoys(index, "shuffle", cfg.mapping.decoys.shuffle_window_aa,
    for_step("a3_feasibility", cfg.seed))`.
11. **Match:** `match_exact(top, index, {"shuffle": decoys}, {p: "shuffle" for p in top},
    min_len=7, max_len=30)`.
12. A peptide **matches** when it has ≥ 1 hit with `search == "target"`, on any locus, replicon
    or strand; several loci count once. Only the matcher's rules apply:
    - I ≡ L;
    - no other mass equivalence;
    - no read-through of a stop codon;
    - `X` codons and query letters `B J O U X Z` never match.
13. `share = n_match / n_top`. (b) holds when `share < 0.20`. With `n_top = 0`, `share` is
    `null`, (b) is not evaluated, and (a) decides.
14. Also reported (none of these decides):
    - the Wilson 95 % interval of `share` (R4);
    - `decoy_share` (shuffled-decoy matches / `n_top`);
    - `contaminant_share` (`match_contaminant_peptides` on `contaminants.fasta`; these count as
      non-matches in (b));
    - matches per replicon;
    - the length histogram;
    - `n_nonstandard`.

### `validation/pilot/metrics/a3_feasibility.json` (aggregates only)

The fields:

- `schema_version` 1, `run`, `cascadia_version`, `weights_sha256`;
- `route` (`B2` | `B1`) with `route_version` (commit or digest), `mzml_sha256`, `mzml_check`
  `{rule: bool}`;
- `genome_sha256`;
- `cutoff`, `cutoff_source`;
- the counts: `n_rows`, `n_below_cutoff`, `n_len_7_30`, `n_pass`, `n_top`, `n_match`,
  `n_nonstandard`;
- the shares: `share`, `share_ci95` `[low, high]`, `decoy_share`, `contaminant_share`;
- `matches_by_replicon` `{contig: n}`, `length_hist` `{len: n}`;
- `rule` `{min_peptides: 500, min_share: 0.20, top_n: 1000}`;
- `fails` `{too_few_peptides: bool | null, low_genome_match: bool | null}`;
- `verdict` (`pass` | `fail`), `reason` (`null` | `too_few_peptides` | `low_genome_match` |
  `both` | `route_unavailable` | `structurally_incompatible` | `tool_incompatible`);
- `evidence` (free text, early-fail only), `job_ids`, `date`.

Keys are sorted, the indent is 2 spaces, and the file ends with a newline, so the bytes are
stable for a fixed input. No sequence and no per-peptide value appears in the file.

`early-fail --reason {route_unavailable,structurally_incompatible,tool_incompatible} --evidence
TEXT --out J` writes the same schema, with `verdict = fail`, the counts `null` and both `fails`
`null`.

Exit codes:

- **0**: pass;
- **3**: fail (the fallback decision, including `early-fail`);
- **1**: error (unreadable result, missing column, failed mzML check);
- **2**: precondition (missing file, unknown option, empty file).

## The smoke run (one GPU, one run)

Run 1 = `WT_1` (`MLA1914689_WT_P1_1_Y-D7_1_13119.d`: 539.9 s, 734 MS1 frames, 8,804 diaPASEF
frames, 12 window groups, 65 window rows). The gradient is only 9 minutes, so the earlier draft's
600-s slice would be longer than the run. **There is no slice: the full run is the smoke test**
(question Q5).

Load estimate: 8,804 frames × about 5.4 windows each (65 rows / 12 groups) ≈ 48 k MS2 spectra.
That fits in one augmentation part (Cascadia splits at 50,000). Charges 2–4 give ≈ 143 k
augmented spectra. The documented demo does 3,868 in about 1 minute on a GPU, so expect ≈ 0.5–1 h.

| Step | Script | Partition / QOS | CPUs | GPUs | Mem | Time limit (expected) | Waits for |
|---|---|---|---|---|---|---|---|
| 1 envs | `build_envs.sbatch cascadia timsconvert` | normal | 4 | 0 | 16 GB | 2 h (20–40 min) | 08a merged |
| 2 assets | `fetch_a3_assets.sbatch` | normal | 1 | 0 | 4 GB | 1 h (≤ 15 min) | 1; ADR-M10-B (weights host) |
| 3 convert | `mzml_pilot.sbatch WT_1` | normal | 8 | 0 | 32 GB | 2 h (unknown; ≤ 30 min) | 1; Khoa's TDF-SDK acceptance |
| 4 smoke | `RUN=WT_1 sbatch scripts/slurm/denovo_pilot.sbatch cascadia` (`--time=04:00:00` on the command line) | gpu / gpu | 8 | **1** | 64 GB | 4 h (≈ 0.5–1 h) | `afterok` of 2 and 3 |
| 5 verdict | `a3_feasibility.py evaluate …` in a `dev` job | normal / dev | 2 | 0 | 8 GB | 30 min (< 5 min) | 4 |

- Peak use is 8 CPUs and **1 GPU** (R14: ≤ 4). GRES is `gpu:1`, since the partition has no type
  names. Before each `sbatch`, the lead sums all running and pending `giap-*` CPUs (≤ 128).
- A TIMEOUT at step 4 gets one resubmission at `--time=24:00:00`; after that, R9.
- 64 GB because Cascadia keeps every augmented spectrum of a part in memory: issue #15 was a
  host-memory kill on the demo.
- Job IDs, `sacct` Elapsed and MaxRSS, the node, the `nvidia-smi` GPU name and memory, and the
  output sizes go to PROGRESS (In progress, then Measured numbers).

## Examples

| # | Input | Expected output |
|---|---|---|
| 1 | `.ssl` with 2,400 distinct passing 7–30-mers; of the top 1,000, 640 match the target six-frame translation and 9 the shuffled decoy | `n_pass` 2400, `n_top` 1000, `share` 0.64, `share_ci95` [0.6098, 0.6692], `decoy_share` 0.009, `verdict` pass, `reason` null, exit 0 |
| 2 | 420 distinct passing peptides, 300 of which match | `fails.too_few_peptides` true, `low_genome_match` false, `n_top` 420, `share` 0.7143 (reported), `verdict` fail, `reason` too_few_peptides, exit 3 |
| 3 | 3,000 passing; 150 of the top 1,000 match | `low_genome_match` true, `share` 0.15, CI [0.1292, 0.1735], fail, `reason` low_genome_match, exit 3 |
| 4 | exactly 500 passing, exactly 100 match (`share` 0.20); then the same with 99 matching (0.198); then 499 passing | 500/100 → pass, exit 0 (both rules say "fewer than"); 500/99 → fail, low_genome_match; 499 → fail, too_few_peptides |
| 5 | rows `PEPTLDEK` 0.95, `PEPTIDEK` 0.91 (hand-made: Cascadia itself writes L), `M[+15.99]PEPTLDEK` 0.85, `[+42.01]-SAMPLER` 0.90, `PEPTLDEK` 0.80 | distinct `PEPTLDEK` (score 0.95, as written `PEPTLDEK`), `MPEPTLDEK`, `SAMPLER`; the 0.80 row is `n_below_cutoff` = 1 (strict `>`); `n_pass` 3 |
| 6 | 1,200 passing; ranks 999–1,002 share one score; folded sequences `AAAKLLLR`, `AAALLLLR`, `ACDEFGHK`, `ADDEFGHK` | ranks 999 and 1,000 are `AAAKLLLR` and `AAALLLLR`. The JSON is byte-identical on a rerun and on row-shuffled input |
| 7 | synthetic circular contig: a peptide whose codons are the last 12 nt and the first 12 nt (8 aa); the same sequence on a linear contig | circular: match, `matches_by_replicon[contig]` = 1; linear: no match |
| 8 | the reverse complement of a planted 30-nt stretch, read in phase 2 of the minus strand | match |
| 9 | a 9-mer whose genomic stretch has an in-frame TAA (the peptide reads it as `Q`) | no target match |
| 10 | `.ssl` with only the header line | `n_rows` 0, `n_pass` 0, `n_top` 0, `share` null, `too_few_peptides` true, `low_genome_match` null, fail, exit 3 |
| 11 | `mzml-check` on a synthetic mzML: 3 MS1 and 12 MS2 spectra, 4 distinct windows of width 25, acquisition stub `{"0": 3, "9": 3, window_groups 2, dia_window_rows 4, mz_range [400, 500], gradient 10 s}`, last RT 9.5 s | pass; `J` lists M1–M8 true. With one MS2 spectrum widened to 100 m/z (> 50, half the span): M4 false, exit 1 |
| 12 | `early-fail --reason route_unavailable --evidence "Khoa declined the TDF-SDK terms (PROGRESS 2026-…)"` | JSON with `verdict` fail, `reason` route_unavailable, counts null, exit 3 |

Examples 7–9 and 11 use synthetic files built at test time, not Nissle or real spectra.

## Edge cases

- **Circular vs linear:** the origin wrap applies only to circular replicons. All four Nissle
  records are circular-tagged. The linear case of example 7 does not match.
- **Minus strand:** frames `-0, -1, -2` come from `build_sixframe`. The peptide reads in the minus
  strand's 5' → 3' direction (example 8).
- **Empty input:**
  - A header-only `.ssl` is a verdict (example 10). Cascadia wrote no row above 0.8, which is the
    signal L6 measures.
  - A zero-byte or missing `.ssl` is exit 2 ("no Cascadia result for WT_1"), not a verdict.
  - An mzML with 0 MS2 spectra fails M2 (exit 1).
  - A job that died is exit 1, then R9; it never becomes a verdict.
- **Ties:**
  - in score at the top-1,000 boundary, broken by folded sequence (example 6);
  - between I/L forms of one peptide: the maximum score and the first as-written form
    (example 5);
  - a score exactly at the cut-off is excluded (strict `>`, example 5).
- **≥ 500 but < 1,000 passing:** `n_top = n_pass`, and (b) runs on all of them.
- **One peptide on several charges or spectra:** Cascadia writes one row per augmented spectrum
  × charge. The distinct step counts the peptide once.
- **Letters outside the residue alphabet** (`X`, `B`, `Z`, `U`): they stay in the counts and never
  match (`n_nonstandard`).
- **N-terminal Met on a non-ATG start:** GTG and TTG translate to V and L in the six frames, so
  such a peptide does not match. This is rare and not corrected; it is conservative.
- **Contaminants** (keratins, trypsin) count as non-matches and are reported separately.
- **Duplicated loci** (IS elements, rRNA operons) count as one match.
- **RT unit:** a 539.9-s run sits next to Cascadia's 500 heuristic. The evaluator never reads RT,
  and M10-09 must read the unit from the values.
- **Unimod:** if the job log shows any `unimod.org` URL, or the `.ssl` has modified sequences
  without mass shifts, the launcher did not take effect: exit 1, R9.
- **Weights sha256 differs from the committed one** (the file was re-uploaded to Drive): exit 2.
  Never run on an unverified file (R6). A changed file is a new asset row, approved by the lead.
- **CUDA not available on the node** (`torch.cuda.is_available()` false): exit 1 before Cascadia
  starts. Cascadia would otherwise fall back silently to the CPU and take days.
- **A diaPASEF variant** (synchro, slice or midia-PASEF): `acquisition.json` says plain diaPASEF
  for all three runs, so this does not arise. If a re-download says otherwise, stop with a lead
  question.
- **Disk:** `mzml_pilot.sbatch` checks that free space under `$GIAP_DATA_DIR` is ≥ 4 × the `.d`
  size before converting (exit 2 otherwise).

## Acceptance criteria

CI, on every commit: no GPU, no Cascadia, no converter, no network, stubs and synthetic files
only (under `tests/data/pilot/a3/`, a few kB).

1. `evaluate` on synthetic `.ssl` fixtures gives exactly the counts, shares, CI (4 decimals),
   verdicts, reasons and exit codes of examples 1–6 and 10, including every boundary of
   example 4.
2. The genome check on a synthetic genome reproduces examples 7–9 and the linear case. Its hits
   equal a direct `match_exact` call on the same index: the module imports `build_sixframe`,
   `build_decoys` and `match_exact`, and defines no matcher of its own.
3. The JSON has every field listed, sorted keys and a trailing newline. It is byte-identical
   across two runs and across row-shuffled input, and contains no amino-acid sequence (no
   string value matches `^[A-Z]{7,}$`).
4. `mzml-check` passes example 11 and fails each of M1–M8 on a fixture that breaks only that
   rule; every failure exits 1 and is named in `J`.
5. `early-fail` writes example 12's JSON for each of the three reasons, exits 3, and rejects any
   other reason with exit 2.
6. `cascadia_offline.py`, run with a stub `cascadia` package and a stub `pyteomics` on
   `PYTHONPATH`:
   - it sets the Unimod URL to the given absolute local path before `main()` runs;
   - it passes `argv` through unchanged;
   - an `http(s)://` argument or a missing file exits 2 without importing `cascadia`.
7. `fetch_a3_assets.sbatch`, with `curl`, `apptainer` and `sha256sum` stubbed:
   - it requests only the URLs in `a3_assets.tsv`, only over `https://`, with no TLS-off flag;
   - while ADR-M10-B is unaccepted (`A3_WEIGHTS_HOST_OK` unset), it exits 2 before any request;
   - a sha256 mismatch exits 1 and leaves no file at the final path;
   - final files are mode 0444 under `$GIAP_DB_DIR`.
8. `mzml_pilot.sbatch` with a stubbed `micromamba`:
   - the default route calls `timsconvert` with `--mode centroid --compression zlib
     --exclude_mobility` and one `--input`, and writes only under `$GIAP_DATA_DIR/pxd078459/mzml/`;
   - `MZML_ROUTE=pwiz` calls `wine msconvert` with `--combineIonMobilitySpectra` and the
     `peakPicking vendor` filter;
   - too little disk (stubbed `df`) exits 2.
9. `denovo_pilot.sbatch`, with `micromamba` and `nvidia-smi` stubbed:
   - `RUN=WT_1` passes exactly one mzML as the first positional argument and the checkpoint as
     the second, through `cascadia_offline.py`, with `-o …/A3_cascadia/WT_1/WT_1`;
   - it uses no `--model`, `--num_workers` or `--output` flag;
   - the working directory is under `$TMPDIR`;
   - a weights or Unimod sha256 mismatch exits 2 before Cascadia starts, and a stubbed
     `torch.cuda.is_available()` that returns false exits 1.
   **Every M3-24 test passes unchanged.**
10. `build_envs.sh cascadia` and `build_envs.sh timsconvert` resolve the `PILOT_TOOLS` entries
    (with `micromamba` stubbed). `build_envs.sh all` lists neither, and `required_tools` never
    returns them.
11. `scripts/ci/check_imports.py` passes, and no file under `src/` or `validation/` imports
    `cascadia`, `timsconvert` or `pyTDFSDK`.
12. `git ls-files` lists no `.ckpt`, `.sif`, `.mzML`, `.ssl`, `.d`, `.tdf`, `.tdf_bin`,
    `libtimsdata*` or `unimod_tables.xml*`. No committed file under `validation/` or
    `tests/data/recorded/cascadia/` is above 5 MB, and the recorded directory holds only
    `PROVENANCE.md`.

Evidence on FarmShare (in the PR description and PROGRESS, not CI):

13. ADR-M10-B is written. Khoa's acceptance of the TDF-SDK agreement and of the weights host is
    recorded in PROGRESS › Answered by Khoa before the first conversion or weights download, or
    E1 / a Blocked row is recorded instead.
14. Both envs build (job ID): `cascadia` reports 0.0.7, the TIMSCONVERT commit is recorded, and the
    locks are committed.
15. `a3_assets.tsv` holds the weights' and the Unimod table's sha256 and bytes, from the official
    sources.
16. WT_1's `mzml_check.json` passes (or the E2 path is documented). The smoke job's ID, wall
    time, MaxRSS, GPU and node are in PROGRESS. `a3_feasibility.json` holds the verdict.
    - On a pass: `PROVENANCE.md` is committed and the 1,000-row sample sits under `$GIAP_ROOT`.
    - On a fail: the lead opens M10-11a.

## Test plan (stubs only)

- `tests/unit/test_m10_08_a3_evaluate.py`:
  - stripping, length filter, strict cut-off and overrides;
  - I/L-folded distinct with the maximum score;
  - ranking ties, both rules and their boundaries;
  - JSON layout, determinism and the no-sequence check;
  - exit codes and `early-fail`;
  - a synthetic genome for circular/linear, the minus strand, stop codons and multiple loci.
- `tests/unit/test_m10_08_mzml_check.py`: small hand-written mzML fixtures, one per rule.
- `tests/unit/test_m10_08_scripts.py`:
  - the M3-24 `tests/helpers/pilot_kit.py` sandbox pattern, with fake `curl`, `apptainer`,
    `micromamba`, `nvidia-smi`, `df` and `sha256sum` on `PATH`;
  - checks the requested URLs, the argv, output paths, exit codes and file modes.
- `tests/unit/test_m10_08_launcher.py`: `cascadia_offline.py` against stub packages in a temp
  directory, run with `python -I` and an explicit `PYTHONPATH`.
- `tests/unit/test_m10_08_registry.py`: the `PILOT_TOOLS` entries and the `build_envs.sh` lookup.
- Not in CI: the conversion, the GPU run and the real verdict. They are evidence only, and no
  `needs_gpu` test is added. The real-row contract test is M10-09's: it reads `$GIAP_ROOT` and
  skips on CI.

## Risks

| Risk | Signal | Plan |
|---|---|---|
| Khoa declines the TDF-SDK terms, or does not answer | PROGRESS | E1 when he declines. While he has not answered, a Blocked row; the envs, evaluator, mzML check and launcher still land (08a/08b need no data) |
| Khoa declines `drive.google.com` | ADR-M10-B | Cascadia cannot run: the lead decides between M10-11 and Q8's image copy (question Q8). Not a timsTOF verdict |
| The Drive file changes or needs a confirmation page | sha256 mismatch; an HTML body | exit 1, retry once later, then Blocked (R6) |
| Cascadia 0.0.7's pins do not solve (`lightning` + `pytorch-lightning`, PR #18) | solve or pip error | conda-forge `pytorch=2.0.*`, `lightning=2.0.*`, `pytorch-lightning=1.9.*` in the conda part, then `pip install --no-deps cascadia==0.0.7`, with the lock recording it; R9 after 3 attempts |
| PyTorch 2.0 with CUDA 11.8 does not run on the `oat` nodes' driver | `torch.cuda.is_available()` false | a newer conda-forge CUDA 11.8 build; R9 |
| pyteomics 5 breaks Cascadia's mzML or Unimod calls | import or attribute errors | pin `pyteomics=4.7.5` (current when 0.0.7 shipped, Feb 2025), recorded (question Q10) |
| TIMSCONVERT's heavy pins (PySide6, `numpy==2.1.1`) fail to install, or `libtimsdata.so` fails on Ubuntu 24.04 | build or conversion error | separate env (done); R9, then B1 (needs Q2) |
| Converted spectra break M4/M5 (whole frames merged) | mzML check | the other route; E2 when both fail |
| Cascadia runs but output is garbage on TOF spectra (domain shift from Orbitrap training) | L6 | that is what L6 measures: a fail opens M10-11a |
| Full run exceeds 4 h | TIMEOUT | one resubmission at 24 h, then R9 |
| GPU queue waits | `PENDING (Resources)` | wait; never more than 1 GPU for this issue |

## Estimate

| Part | PR | Code lines | Test lines |
|---|---|---|---|
| `cascadia.yml`, `timsconvert.yml`, `PILOT_TOOLS` and the `build_envs.sh` lookup | 08a | 60 | 50 |
| `cascadia_offline.py` | 08a | 35 | 50 |
| `fetch_a3_assets.sbatch`, `a3_assets.tsv` | 08a | 75 | 80 |
| `mzml_pilot.sbatch` | 08a | 55 | 50 |
| `a3_feasibility.py mzml-check` | 08a | 45 | 60 |
| 08a total | | **≈ 270** | ≈ 290 |
| `a3_feasibility.py evaluate`, `early-fail`, JSON | 08b | 145 | 120 |
| `denovo_pilot.sbatch` changes | 08b | 40 | 40 |
| 08b total | | **≈ 185** | ≈ 160 |
| Total | | ≈ 455 (the index said 150) | ≈ 450 |

Wall time:

- jobs: ≈ 1 h for the envs and assets; the conversion and the smoke run are unmeasured, likely
  ≤ 2 h together;
- agent time: ≈ 1.5 days;
- plus queue waits and Khoa's answer on ADR-M10-B (the critical path).

## Questions for the lead (reversible defaults in use)

| # | Question | Default |
|---|---|---|
| Q1 | ≈ 455 code lines is over the ~400-line PR cap. Split into 08a (envs, registry, launcher, fetch, conversion, mzML check, ADR draft) and 08b (evaluator, `denovo_pilot` changes, smoke run, verdict)? This reverses the earlier "one PR" (Q7) | split; 08a needs no Khoa answer to merge |
| Q2 | Route order: B2 (TIMSCONVERT, one Bruker licence) first, then B1 (ProteoWizard image, seven vendor licences). This **reverses** the 2026-10-23 lead decision "B1, then B2", which came before Khoa's 2026-10-24 permission and before the licences were read. Is B1 allowed at all without asking Khoa about the six non-Bruker EULAs? | B2 first. B1 only after Khoa also accepts the six non-Bruker EULAs (asked in the same batch, so B1 is ready if B2 fails) |
| Q3 | Earlier Q1 ("cut-off if none is documented") is settled by the code: 0.8, strict, CLI default. Confirm that L6's "Cascadia's own score cut-off" is the 0.8 default and not a tuned value | 0.8, `cli_default`; never tuned on pilot data (R3) |
| Q4 | Earlier Q6 committed the 1,000-row sample on a pass. Khoa's 2026-10-23 rule (M10-02 decision 4) keeps samples on FarmShare. Apply it to Cascadia, committing only `PROVENANCE.md` (with the `.ssl` header taken from the tool's source code, not from data) and no `schema.tsv`? | yes; the sample goes to `$GIAP_ROOT/pilot/recorded/cascadia/0.0.7/` |
| Q5 | Earlier Q10 (a 10-minute slice first) is moot: the whole run is 9 minutes. Run the full WT_1 as the only smoke job? | yes |
| Q6 | E3 `tool_incompatible`: Cascadia crashes on a valid mzML for a data reason after R9. Count that as "fails on timsTOF" (go to M10-11), or escalate to Khoa? | count it as a fail, with the lead confirming the diagnosis in PROGRESS; Khoa is told in the next batch |
| Q7 | The TDF-SDK grants use "on any of your computing devices". GIAP reads FarmShare (a Stanford cluster we hold an allocation on) as such. Should ADR-M10-B state this reading for Khoa to confirm? | yes, stated in the ADR; Khoa's acceptance covers it |
| Q8 | The weights host (Google Drive) is not in §3.2. One alternative stays inside §3.2: `quay.io/protio/cascadia:0.0.7@sha256:a0227b5ed831bf541e9d5b9e4f9cf57f36dc856c92b87e9311a7b9a319b163db` (6.1 GB), built by paper co-author M. Riffle for the nf-skyline-dia-ms workflow that the paper names. It bakes in `/usr/local/bin/cascadia.ckpt`, but the Dockerfile does not say where the checkpoint came from, so R6 may count it as a re-hosted copy. Use it only as a sha256 cross-check after the Drive download? | Drive via ADR-M10-B (the docs' link); the image is not used. A cross-check only if the lead asks |
| Q9 | Unimod table: use the copy vendored in psims 1.4.0 (conda-forge, §3.2), or download `https://www.unimod.org/xml/unimod_tables.xml` (official, but not in §3.2, so an ADR row)? | the psims copy |
| Q10 | Pin `pyteomics` below 5 (`4.7.5`) in the Cascadia env to match the 0.0.7 era, or let the solver pick? | pin `4.7.5`; recorded in the env header |
| Q11 | Which checkpoint: `cascadia.ckpt` (the general model of the docs and paper), `cascadia_astral_tuned.ckpt` or `cascadia_pretrained.ckpt`? | `cascadia.ckpt` only. Trying a second checkpoint after a fail would tune on the pilot (R3 spirit) |
| Q12 | This spec reuses the file name of the draft on `m10/a3-feasibility`. Retire that branch (no PR) once this one merges? | yes: the lead closes it; its lead decisions 2 and 3 are carried here |

## Needs Khoa (batched by the lead; ADR-M10-B)

1. Accept Bruker's **TDF-SDK** licence agreement (text above; S9 copy, since Bruker's page S14 is
   login-gated) for route B2, on FarmShare, with no redistribution. Without it, E1 fires and A3
   moves to the Casanovo fallback.
2. Add **`drive.google.com`** (the Noble Lab folder the Cascadia docs link to) as the source of
   the Cascadia 0.0.7 weights (Apache-2.0 per the Nat Methods Code availability statement; no
   published checksum).
3. Conditional (Q2): accept the six other vendor EULAs bundled with the ProteoWizard image, only
   needed if B2 fails.

## Amendments after the PR #301 code review (lead, 2026-10-24)

These supersede the earlier text where they differ:
- **Unimod:** Cascadia 0.0.7 fetches Unimod over plain HTTP twice: at the end of a run via
  pyteomics, and at the start of every run via psims (`GenericModification(...).mass`). The
  launcher blocks both. It points pyteomics at the local table, registers a local psims resolver,
  and installs a Python socket guard. `cascadia.yml` pins `lxml>=5.4`.
- **cascadia.yml:** the conda channel is conda-forge only; the pip part installs
  `cascadia @ <PyPI wheel URL>#sha256=50b52043…f092a3` and `==`-pinned dependencies from PyPI.
  The env is therefore not "conda-forge only" overall (rows 1 and 2 of the deliverables).
- **Vendor licences:** the ProteoWizard route needs Bruker's plus **five** other vendor EULAs
  (six in all), not six others (§B1, Q2, Khoa item 3).
- **The weights fetch** rejects an HTML page or a body without the zip magic in both modes, and
  its `--time` is 02:00:00 so curl's worst case fits.
- **mzml-check M7:** the time unit is read from `unitAccession`, else `unitName`; no unit fails M7.

## Amendments (2026-10-24, 08b)

Lead decisions for 08b; they supersede the earlier text where they differ:
- **CLI:** `evaluate <RUN.ssl> --run <RUN> --out <J> [--genome] [--bakta-json] [--contaminants]
  [--cutoff]`. `--genome` and `--bakta-json` default to the Nissle files of step 8,
  `--contaminants` to `$GIAP_DB_DIR/contaminants/contaminants.fasta`. Exit codes: 0 pass, 3 fail
  on timsTOF, 2 precondition, 1 bad input. M10-10's `cascadia_convert.sbatch` already calls this
  interface. `early-fail` takes `--run` optionally (default `WT_1`).
- **Stripping (step 2):** after deleting the groups, signed masses and `-`, uppercase first, then
  delete every character outside A–Z; then fold I → L.
- **Counting:** `n_below_cutoff` counts after the length filter. `n_nonstandard` and `length_hist`
  are over the distinct passing set. A peptide on two replicons counts once in each
  (`matches_by_replicon`). A contaminant peptide counts as a non-match even if it also matches
  the genome (the conservative choice), in `n_match`, `share` and `matches_by_replicon`;
  `decoy_share` is diagnostic and not adjusted.
- **Provenance values** (`weights_sha256`, `route`, `route_version`, `mzml_sha256`,
  `mzml_check`, `job_ids`, `date`) come from the job (sbatch environment and logs), passed as
  options (`--weights-sha256`, `--route`, `--route-version`, `--mzml-sha256`, `--mzml-check J`,
  `--job-id` repeatable, `--date`, default today UTC), and are filled in by the validation run.
  The unit tests check only the keys. `--detail TSV` writes the per-peptide top list (under
  `$GIAP_ROOT` only, never committed).
- **Unimod URL:** `denovo_pilot.sbatch` tees each Cascadia call into
  `$GIAP_ROOT/pilot/A3_cascadia/<RUN>/cascadia.log` and exits 1 when that log shows a `unimod.org`
  URL. The Cascadia input directory is `MZML_DIR` (default `$GIAP_DATA_DIR/pxd078459/mzml`); the
  run list is `$RUN`, else the second argument, else `WT_1 WT_2 WT_3`.
- **Env build:** building `cascadia` on a CPU node needs `CONDA_OVERRIDE_CUDA=11.8` (the README's
  A3 section).

## Amendments after the 08b code review (lead, 2026-10-24)

These supersede the earlier text where they differ:
- **S1, resources:** `denovo_pilot.sbatch` asks for `--cpus-per-task=16 --mem=64000M`. The `gpu`
  partition has MaxMemPerCPU 4000M, so `--mem=64G` on 8 CPUs was silently raised to 17 CPUs; the
  CPUs asked for are now the CPUs held (R14 accounting).
- **S2, scratch:** the job makes one per-job directory
  `SCRATCH_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/giap-a3.XXXXXX")"` and runs each Cascadia call in
  `$SCRATCH_ROOT/<RUN>`; an EXIT trap removes `$SCRATCH_ROOT`. It never creates or removes
  `$TMPDIR/<RUN>` (deliverable 8, the scratch row and the command block above read accordingly).
  The scratch is still under `$TMPDIR`. The log's `unimod.org` check runs before the exit-code
  check, and an exit 0 without `<RUN>.ssl` is exit 1 (a header-only `.ssl` is still a result).
- **S3, named modifications:** the edge case "modified sequences without mass shifts → exit 1" is
  implemented: a `[...]` group whose content is not a signed number (`[+-]\d+(\.\d+)?`, e.g.
  `[Oxidation]`, `[UNIMOD:4]`, `[]`) is an `InputError` (exit 1, R9) naming the `.ssl` line.
  `(...)` groups are still stripped. Step 2's separate deletion of signed masses and `-` is
  subsumed by deleting every character outside A–Z after uppercasing.
- **N6, preconditions (exit 2):** a missing `--mzml-check` file (an unreadable one, or one without
  `rules`, is exit 1), a whitespace-only `.ssl` (no header; it counts as an empty file), and a
  `--detail` path inside the checkout that is not git-ignored.

### Lead decision: the CPU route (2026-10-25, #334 review item 6)
GIAP_SPEC.md §15.1 puts A3 on the `gpu` partition. At Khoa's request (2026-10-24) the three WT runs
were predicted on CPU instead, as a Slurm array (`cascadia_array.sbatch`, ≤ 3 h per task, R14), in
two forms:
- **Fixed blocks** (`predict`, used for the first jobs). Task 1889429_0 finished at 2:58:39 of a 3:00
  limit, so this mode has no time margin.
- **The work queue** (`work`, used for the rest of the run). Chunks lie on the global grid
  `[k·20, (k+1)·20)`; each is claimed atomically within a generation and saved atomically;
  `pending_chunks` skips batches already held by chunk or shard files, and `merge` refuses gaps,
  overlaps and records whose prediction or row counts do not match their block.

`write_results` runs once on the concatenated predictions, so the merged table equals an unsplit
**CPU** run of the same build (same row numbers, threshold 0.8, time width, model parameters). It
is **not** bit-for-bit a GPU run: float reductions differ, so a row near the strict `> 0.8` cut can
flip. The WT_1 verdict is far from its rule (863 peptides vs 500; share 0.636, CI 0.604–0.668, vs
0.20), so a handful of flipped rows cannot change it. Open, not needed for the verdict: claims are
keyed by chunk name, so releasing a reserved marker while a generation is planning could give two
overlapping claims (`merge` would refuse rather than write a bad table); claiming by grid cell
would close it.
