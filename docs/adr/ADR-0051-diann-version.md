# ADR-0051: DIA-NN 2.3.1 Academia (Apptainer image) instead of the spec's 2.7.0
- Status: proposed
- Date: 2026-10-23
- Spec sections: §15.1 (arms A1, A2: "DIA-NN 2.7.0"), §15.2 steps 3, 5 and 7 and the DIA-NN template, §3.3 (`GIAP_TOOLS_DIR`), §3.6 (`needs_diann`: "DIA-NN 2.7.0 under `$GIAP_TOOLS_DIR`"), §14.8 item 6, §0.1 item 4, §17 rows on DIA-NN; R5, R6, R11, R13
- Needs Khoa: yes. This changes a dependency version and the way it is installed (§14.8 items 1 and 6). It writes down Khoa's 2026-10-20 acceptance of the M3-24 README and does not replace it.

## Context
The spec names DIA-NN **2.7.0** in four places:
- the arms (§15.1);
- the install step: "Khoa: accept the DIA-NN 2.7.0 licence and unpack the Linux binary under
  `$GIAP_TOOLS_DIR/diann-2.7.0`" (§15.2 step 3, §0.1 item 4, §14.8 item 6);
- the test marker (§3.6);
- the recorded-sample path `tests/data/recorded/diann/2.7.0/` (§15.2 step 7).

What exists is different. Khoa supplied an Apptainer image of **2.3.1 Academia** and pointed at
their reference script `/farmshare/home/classes/bios/270/jordan/run_2step.sh` (M3-24 lead decision
1). The M3-24 pilot kit is written against it.

Every fact below was read from disk. Nothing was inferred about DIA-NN versions that are not on
disk.

From `validation/pilot/README.md` (M3-24, #122):
- version string `DIA-NN 2.3.1 Academia`, compiled 2025-12-05;
- the image `/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif` (532 MB), with sha256
  `1946452d56d2fb0f66e05ced473df88e53484e28040f1ecbd83c4a43405a6381` (pinned 2026-10-20) and the
  binary `/diann-2.3.1/diann-linux`;
- the module `apptainer/1.4.4`, binding `/home/classes`, `/scratch` and `/farmshare`;
- the section "Version actually used" says: "Khoa reviewed and accepted the flag table and this
  section on 2026-10-20 ('looks fine') … 2.3.1 Academia is the version of record unless he says
  otherwise." It adds: "If a 2.7.0 build arrives later, only `diann_A1.cfg`, `diann_A2.cfg` and
  this section change."

The README also documents the two-step workflow (M3-24 decision 4): "In DIA-NN 2.x, combining
`--fasta-search --predictor` with raw files in one call disables the FASTA digest". Step 1
therefore predicts a library from the FASTA files, and step 2 searches the runs against it. Khoa's
reference script carries the same note.

The cfgs (`validation/pilot/diann_A{1,2}.cfg`) deviate from the §15.2 template. All of this is
listed in the README:
- `--met-excision` is left out, because GIAP's digest keeps the initiator Met;
- `--lib ""` is not used;
- added: `--min-pr-mz 300 --max-pr-mz 1800 --min-fr-mz 200 --max-fr-mz 1800`, `--peptidoforms`
  and `--verbose 1`;
- added in step 2: `--matrices`, plus `--gen-spec-lib --out-lib <arm>_empirical.parquet`;
- Khoa's `--reanalyse` and `--rt-profiling` are not used, because match-between-runs weakens
  run-level FDR control.

Checked on 2026-10-23, in the session allocation:
- The shared image is 532,586,496 bytes, mode `-rwxrwxrwx`, owner `khoang99`, group `operator`,
  mtime 2026-07-09 18:20 −0700. Its directory is `drwxrwxrwx`.
- The sha256 still equals the pin.
- `apptainer inspect` gives the labels build date `Tuesday_6_January_2026_19:18:54_PST`,
  bootstrap `docker`, from `khoayelo189/diann:2.3.1`.
- `/diann-2.3.1/` holds `diann-linux` (25,747,216 bytes), `libtimsdata.so`,
  `ThermoFisher.CommonCore.*.dll`, `libtorch_cpu.so`, `models/*.pt` and `diann-stats.py`.
- `find / -maxdepth 3 -iname '*licen*'` inside the image finds no licence file.
- `$GIAP_TOOLS_DIR` (`local/tools`) is empty: no `diann-2.7.0`. `module avail diann` finds
  nothing.
- `run_2step.sh` collects `*.raw` files (Thermo) and a human FASTA. GIAP has not yet run this
  image on timsTOF `.d` data.

What is **not** known from anything on disk:
- whether a 2.7.0 build is available to Khoa;
- how 2.7.0 differs from 2.3.1 in flags, report columns or scoring;
- who built `khoayelo189/diann:2.3.1`, and from which DIA-NN release archive;
- the text of the Academia licence.

## Options
1. **Adopt 2.3.1 Academia as the version of record** (the M3-24 kit, index decision 2). The pilot
   runs now. The spec's 2.7.0 paths are read as their 2.3.1 equivalents:
   - the shared image, verified against the pin and copied read-only to
     `$GIAP_ROOT/images/diann_2.3.1.sif` (index lead decision 2, M10-02), in place of
     `$GIAP_TOOLS_DIR/diann-2.7.0`;
   - `tests/data/recorded/diann/2.3.1/` in place of `…/2.7.0/`.
2. **Wait for 2.7.0** under `$GIAP_TOOLS_DIR/diann-2.7.0`, installed by Khoa. This follows the spec
   to the letter. M10 pass 1 stops at M10-01 until then. The flag table would need Khoa's check
   again against the 2.7.0 documentation, and the effect on the gate's timeline is unknown.
3. **Run 2.3.1 now and re-run on 2.7.0 if it arrives.** Option 1 now. When a 2.7.0 build exists,
   only the cfgs, the README, the image or binary pin and the recorded sample's version directory
   change. DIA-NN, the conversion and the A1/A2 arms run again before the report freeze (M10-14),
   and both versions' numbers are kept in the pilot report.

## Decision
Proposed: **option 1, with option 3 kept open.** The pilot (M10-02 … M10-05) runs on 2.3.1
Academia from the read-only copy. The version string, the image sha256, the two-step workflow and
the template deviations above are the record. A 2.7.0 build is not awaited. If Khoa provides
one, option 3 applies, and this ADR is amended with the new pin and the re-run's job IDs.

Image handling (M10-02):
- every DIA-NN job checks the copy's sha256 against `validation/pilot/diann_image.sha256` and
  refuses a writable copy;
- the install job hashes the copy it made, not only the world-writable source;
- a mismatch is a Blocked row and a question to Khoa, never a re-pin or a download from
  elsewhere (R6);
- the image, and every file inside it, is never committed (R11).

## Questions for Khoa
1. Do you accept 2.3.1 Academia as the version of record for the pilot (option 1/3), or should the
   pilot wait for 2.7.0 (option 2)?
2. The copy is made by the agent from your image, under `$GIAP_ROOT/images/`, not under
   `$GIAP_TOOLS_DIR`. §14.8 item 6 reserves *installing* DIA-NN to you. Is a verified, read-only
   copy of the image you supplied acceptable?
3. The image's labels name `khoayelo189/diann:2.3.1` on Docker Hub as its source. Can you confirm
   that this image is yours, built from an official DIA-NN 2.3.1 release?
4. The image has no licence text. §15.2 step 7 allows committing a 1,000-row sample of each
   report, because DIA-NN output tables are not binaries. Does the Academia licence you accepted
   permit that?
5. If PXD078459 turns out not to be diaPASEF (M10-01), the settings change in an amendment to this
   ADR, with your check of the flags.

## Evidence
- `validation/pilot/README.md` ("DIA-NN: the container", "Version actually used", "The two-step
  DIA-NN workflow", "§15.2 template flags that the cfg files do not use");
- `docs/specs/M3-24-pilot-kit.md` lead decisions 1–7 and 9;
- `tests/unit/test_m3_24_pilot_kit.py`: the version string, the pin, the cfg-versus-digest checks,
  and the re-hash of the image on disk (`test_the_pinned_sha256_is_still_the_image_on_disk`);
- `docs/specs/M10-index.md`: decision 2 and lead decisions 2–3;
- `docs/specs/M10-02-diann-searches.md`: "What exists", the facts checked on 2026-10-23.

No search has run yet, so there are no before/after metrics. The first evidence will be the
M10-02 job IDs and `check.json` results.

## Consequences
- The spec text is not edited (R13). Once this ADR is accepted, it overrides the 2.7.0 mentions
  listed under "Spec sections" for the pilot.
- `needs_diann` and the capability probe still look for the jordan path or `diann-2.7.0`. Teaching
  them the copy is a backlog chore (M10-02 Q8).
- The contract-test sample lives under `tests/data/recorded/diann/2.3.1/`. A later 2.7.0 run adds
  `…/2.7.0/` and keeps the old one only if both are tested.
- The pilot report states "DIA-NN 2.3.1 Academia" for A1 and A2, and cites this ADR.
- Risk: 2.3.1 might differ from 2.7.0 in report columns. The converter's contract test on the
  real sample (M10-03) catches any column difference that matters to GIAP.

## Amendment (2026-10-24, M10-03, PR #288)
Khoa ruled on 2026-10-23 that no row-level DIA-NN content is committed, and on 2026-10-24 that
aggregates (counts, hashes, named example genes) may be. This replaces §15.2 step 7's committed
recorded sample for DIA-NN:
- The converted tables (`a1.tsv.gz`, `a2.tsv.gz` and their abundance tables) and any sample drawn
  from the reports stay under `$GIAP_ROOT/pilot/`, never in the repository.
- What is committed is the conversion record in `validation/pilot/README.md` › "Conversion
  (M10-03)": the job IDs, the sha256 of every input report and output table, and the row counts.
- The contract test (`tests/contract/test_m10_03_diann_contract.py`) runs on FarmShare against the
  `$GIAP_ROOT` copies; CI runs the converter's unit tests on synthetic fixtures. The
  `tests/data/recorded/diann/2.3.1/` path named above therefore holds no data rows.
