# M3-24: Pilot kit prepared (§15), DIA-NN request to Khoa

Milestone M3 · depends on: M3-23 (`giap makedb`); recorded Nissle Bakta full-DB JSON (M2-10 /
§15.2 step 1) · creates `validation/pilot/README.md`, `validation/pilot/diann_A1.cfg`,
`validation/pilot/diann_A2.cfg`, `validation/pilot/db_a2/strain_db.faa.gz` (+ `id_map.tsv`,
`makedb_manifest.json`), `scripts/slurm/diann_pilot.sbatch`, `scripts/slurm/denovo_pilot.sbatch`,
`scripts/slurm/fetch_pxd078459.sbatch`, `validation/fetch_uniprot_k12.py` (if not already from
M3-22) · **external tools:** DIA-NN 2.7.0 (installed by Khoa; licence), Cascadia 0.0.7 (GPU) —
prepared only, not run in M3 · estimate ≈ 250 code lines (scripts)

## Spec sections implemented

- §13 M3 deliverable: "**pilot kit (§15) prepared**, and Khoa asked to install DIA-NN".
- §15.2 steps 1–3 done now; steps 4–9 scripted now, run at M10:
  1. "Run `scripts/slurm/record_bakta_full.sbatch` (Bakta 1.12.2, full DB v6.0) on Nissle and K-12
     if M2 has not already done so. Commit the gzipped JSON/GFF3/GBFF under
     `tests/data/recorded/bakta/1.12.2-full/`"
  2. "Run `giap makedb --annotation <nissle JSON> --outdir validation/pilot/db_a2`. Commit
     `strain_db.faa.gz` and a `README.md` with the exact DIA-NN settings for A1 and A2
     (`validation/pilot/diann_A1.cfg`, `diann_A2.cfg`)"
  3. Khoa: "Accept the DIA-NN 2.7.0 licence and unpack the Linux binary under
     `$GIAP_TOOLS_DIR/diann-2.7.0`; check every flag of the template below against its
     documentation".
  4–6. Download PXD078459 (3 wild-type runs, checksums, confirm acquisition mode), DIA-NN for A1/A2
     (`normal`, 16 CPU, 64 GB, `--report-decoys`), Cascadia on the `gpu` partition with the
     Casanovo/UTI89 fallback.
- §15.1 arms A1 (UniProt K-12 UP000000625 + contaminants), A2 (Nissle Bakta proteins incl. sORFs +
  contaminants, from `giap makedb`), A3 (Cascadia). The DIA-NN template of §15.2, verbatim, is the
  content of the two cfg files (one flag per line; `--fasta` differs).
- R14 (job names `giap-*`, caps), §3.8 job header, §14.8 (Khoa checkpoints: the DIA-NN licence is
  Khoa's).

## Deliverables and examples

| # | Item | Expected |
|---|---|---|
| 1 | `diann_A1.cfg` | the §15.2 flags, `--fasta $GIAP_DATA_DIR/uniprot/UP000000625.fasta --fasta contaminants.fasta`, `--report-decoys`, `--threads 16`, `--out report.parquet` |
| 2 | `diann_A2.cfg` | identical except `--fasta validation/pilot/db_a2/strain_db.faa` (built with `--no-contaminants`, M3-23 Q2) |
| 3 | `diann_pilot.sbatch A1` | `#SBATCH --job-name=giap-diann-A1`, `--partition=normal`, `--cpus-per-task=16`, `--mem=64G`; refuses to start (exit 2, clear message) when `$GIAP_TOOLS_DIR/diann-2.7.0/diann` is absent |
| 4 | `denovo_pilot.sbatch` | `--partition=gpu --gres=gpu:1`, job name `giap-denovo`; fallback branch documented (Casanovo 5.2.1 on PXD053727 UTI89 DDA) |
| 5 | `fetch_pxd078459.sbatch` | downloads to `$GIAP_DATA_DIR/pxd078459/`, writes `SHA256SUMS`, never under the repository |
| 6 | `README.md` | the exact commands for steps 4–9, the three `giap run` lines of §15.2, and who does what |
| 7 | `PROGRESS.md` "Questions for Khoa" (by the lead, not this issue's PR) | a row asking Khoa to accept the DIA-NN 2.7.0 licence and unpack it under `$GIAP_TOOLS_DIR/diann-2.7.0`, and to check every flag of the template |

## Acceptance criteria

- [ ] `bash -n` passes on every new sbatch script; `scripts/ci/check_*` guards pass (job names,
  caps, no absolute paths outside `$GIAP_ROOT`).
- [ ] `tests/unit/test_pilot_kit.py`: the two cfg files contain exactly the §15.2 flag set (parsed
  and compared) and differ only in `--fasta`.
- [ ] `strain_db.faa.gz` was produced by `giap makedb` from the recorded Nissle JSON;
  `makedb_manifest.json` lists the JSON's SHA-256; size is a few MB (R11).
- [ ] The DIA-NN request appears in `PROGRESS.md` (lead) and the `needs_diann` marker skips cleanly
  until it is installed.

## Questions for the lead

1. The template has `--max-pep-len 30` while GIAP accepts up to 40 aa. **Proposed default:** keep
   the template as written (Khoa checks flags in step 3).
2. Is `contaminants.fasta` for DIA-NN the same Frankenfield library that `giap setup
   --contaminants` downloads? **Proposed default:** yes, the same file, copied (not regenerated)
   so GIAP and DIA-NN share one contaminant set.
