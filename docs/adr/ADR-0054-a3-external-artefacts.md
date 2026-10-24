# ADR-0054: External artefacts of the A3 de novo arm (ADR-M10-B): conversion route, TDF-SDK licence, Cascadia weights host, Unimod table
- Status: accepted-by-khoa (2026-10-24, in session: "go with your recommendation": Bruker TDF-SDK accepted for route B2; `cascadia.ckpt` allowed from the authors' Drive folder with its sha256 pinned; the ProteoWizard route not accepted)
- Date: 2026-10-24
- Spec sections: §15.1 (arm A3), §15.2 steps 6-7, §17, §3.2, §14.8 item 1; R5, R6, R9, R11, R14
- Needs Khoa: yes (a new download source, `drive.google.com`, and a proprietary licence, Bruker's TDF-SDK agreement)

## Context
Arm A3 runs Cascadia 0.0.7 (Apache-2.0, PyPI) on the PXD078459 diaPASEF runs
(docs/specs/M10-08-a3-feasibility.md). Three external artefacts are not covered by §3.2:

1. Cascadia reads mzML only, so the Bruker `.d` runs need a `.d → mzML` conversion. Both candidate
   routes decode `analysis.tdf_bin` with Bruker's proprietary TDF-SDK (`libtimsdata.so`):
   - **B2**: TIMSCONVERT 2.0.0 (GitHub tag `v2.0.0`, commit
     `3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b`, Apache-2.0) with pyTDFSDK at commit
     `4808da81d3e37f82da96f392dbb446c3733edacb` (Apache-2.0), which bundles the TDF-SDK. One
     proprietary licence: Bruker's TDF-SDK agreement.
   - **B1**: ProteoWizard `msconvert` from the Docker Hub image
     `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses@sha256:34ed787cd4c3182c2d81adcbf7514031e594609ca460abc5d3d114dc2334dae8`,
     under Wine. Six vendor licences: Bruker's, and five more (Agilent, Sciex, Shimadzu, Thermo,
     Waters) for libraries this pilot never loads.
2. The Cascadia weights (`cascadia.ckpt`, Apache-2.0 per the Nat Methods 2025 Code availability
   statement) are published only in the Google Drive folder that the Cascadia docs link,
   https://drive.google.com/drive/folders/1UTrZIrCdUqYqscbqga_KdX8kc8ZjMMfr. Drive is not a §3.2
   source, and no checksum is published.
3. Cascadia 0.0.7 fetches `http://www.unimod.org/xml/unimod_tables.xml` (plain HTTP, a host
   outside §3.2) twice per run: at the start, when its tokenizer resolves named modifications
   (pyteomics `proforma` -> psims `obo_cache.resolve(UNIMOD_OBO_URL)` -> `Unimod()`), and at the
   end, through pyteomics `mass.unimod.Unimod()`.

The TDF-SDK agreement (read from the copy in ProteoWizard's repository, since Bruker's own page is
login-gated) allows use "on any of your computing devices" (§2.1), forbids redistribution (§3.2),
is free of charge (§6.1) and binds the user only "after having accepted this Agreement".

## Decision
1. **Route order: B2, then B1, then E1.** B2 is the conversion route
   (`workflow/envs/timsconvert.yml`, `scripts/slurm/mzml_pilot.sbatch`). B1 is used only if B2
   fails after R9 **and** Khoa accepts the five non-Bruker EULAs. If neither route is available,
   rule E1 (`route_unavailable`) fires and A3 moves to the Casanovo fallback (M10-11a).
2. **Licences.** Khoa accepts Bruker's TDF-SDK licence agreement for route B2 (2026-10-24, "the
   Bruker TDF-SDK licence is accepted for route B2 (TIMSCONVERT + pyTDFSDK)"). GIAP reads the
   FarmShare allocation it runs under as "your computing devices" (§2.1; spec question Q7).
   Nothing from the env, the SDK or the converted data is redistributed or committed (R11).
   Khoa does **not** accept the ProteoWizard five-EULA route (2026-10-24): B1 stays unavailable,
   no image is pulled, and `fetch_a3_assets.sbatch` has no pull step.
3. **Weights host.** `drive.google.com` (file downloads through `drive.usercontent.google.com`,
   HTTPS) is the source of `cascadia.ckpt`, Drive file `14b7Y2-5ZqLY807MYtggoNocSfELdfHy9`
   (Khoa, 2026-10-24: "`cascadia.ckpt` may be fetched from the authors' Drive folder linked in the
   Cascadia docs, with its sha256 recorded and pinned"). The first download
   (`fetch_a3_assets.sbatch --resolve`) is trust on first use; its sha256 and size go into
   `validation/pilot/a3_assets.tsv`, and every later fetch and run must match it (R6). A changed
   file is a new asset row, approved by the lead. Only `cascadia.ckpt` is used (spec Q11). The
   `quay.io/protio/cascadia` copy is not used (spec Q8).
4. **Unimod table.** The `unimod_tables.xml.gz` that psims 1.4.0 (conda-forge, Apache-2.0)
   vendors is gunzipped from the built `giap-cascadia` env to
   `$GIAP_DB_DIR/cascadia/unimod_tables.xml` and pinned by sha256. Cascadia runs through
   `scripts/pilot/cascadia_offline.py`, which blocks both Unimod fetches before Cascadia is
   imported: it registers a psims resolver for the Unimod OBO URL that reads that file (start of
   run) and points pyteomics' `_unimod_xml_download_url` at it (end of run). It also refuses every
   Python-level socket connection, and `cascadia.yml` pins `lxml>=5.4`, which refuses network
   entities, so a missed call site fails instead of reaching unimod.org.

`fetch_a3_assets.sbatch` refuses to request the weights until this ADR is accepted
(`A3_WEIGHTS_HOST_OK=1`, set by the lead once the status changes).

## Evidence (before → after: metrics, tests, simulations)
- Licence texts and sources: docs/specs/M10-08-a3-feasibility.md, "What the official sources
  say" (S1-S18) and "Licence findings".
- Tests: `tests/unit/test_m10_08_scripts.py` (no request before acceptance, https only, sha256
  mismatch leaves no file, 0444 files), `tests/unit/test_m10_08_launcher.py` (both Unimod
  call sites read the local table, set before Cascadia is imported; Python sockets are refused), `tests/unit/test_m10_08_committed.py` (no weights, SDK, image,
  mzML or Unimod copy is tracked).

## Consequences
- Two pilot-only envs (`giap-cascadia`, `giap-timsconvert`; conda channel conda-forge only, pip
  requirements from PyPI, pinned) in
  `giap.annotate.tools.PILOT_TOOLS`; they never join `TOOLS`, `required_tools` or
  `build_envs.sh all`.
- `cascadia.ckpt` is a pickle (`torch.load`, `weights_only=False` on torch 2.0), so loading it
  runs code. It is trusted on first use, pinned by sha256 and loaded only inside the GPU job.
- The weights may disappear or change on Drive: a sha256 mismatch exits non-zero and the lead
  records a Blocked row; it is never evidence about timsTOF data.
- If B2 fails after R9, A3 needs Khoa again (B1's five EULAs) or goes to E1 and M10-11a.
