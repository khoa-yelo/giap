# M2-10: Recording harness, test genomes, wave-A recorded fixtures (Bakta light+full, TerminatorNet, geNOMAD, antiSMASH, BacTermFinder)

Milestone M2 · depends on: M2-01 (envs), M2-03 (`bakta_cmd`), M2-05 (SIM100K), M2-08 (DBs) ·
estimate ≈ 300 code lines (Python + sbatch; fixtures not counted) · capabilities: `needs_slurm`,
`needs_conda` (or `needs_apptainer`), `needs_bakta_db`, `needs_bakta_full`, `needs_genomad_db`,
`needs_net(api.ncbi.nlm.nih.gov)`

## Spec sections implemented

- §11.4: "**Recorded native outputs** live in `tests/data/recorded/<tool>/<version>/<input>/`, with a
  `PROVENANCE.md` giving the command, environment lock or digest, date and host." "The agent records
  them by running the tools on the 100-kb synthetic genome on FarmShare (a batch job, or the
  `record-fixtures` mode of the nightly Slurm suite)." "Re-recording happens only on a tool upgrade,
  which requires an ADR and a contract-test diff."
- §11.4 G17: "`scripts/slurm/record_bakta_full.sbatch` (16 CPU, 64 GB) runs Bakta 1.12.2 with the full
  DB v6.0 on: K-12 MG1655 (GCA_000005845.2); Nissle 1917 (GCA_043228065.1); the 100-kb synthetic
  genome. … commits the gzipped JSON, GFF3 and GBFF (a few MB each) under
  `tests/data/recorded/bakta/1.12.2-full/`. Until then … the full-DB contract test is marked `pending_sample`."
- §11.4: "The K-12 and Nissle FASTA files (gzipped, < 2 MB each) are committed under
  `tests/data/genomes/`, with a provenance note."
- R11: "No single file may exceed 5 MB, and `tests/data/` stays under 50 MB in total. Compress text
  fixtures with gzip." §3.7 row "Bakta with the full DB … `record_bakta_full.sbatch` (16 CPU, 64 GB)".

## Public interface

`scripts/record_fixture.py`
```
record_fixture.py TOOL --input NAME=PATH [--input …] --label INPUT_LABEL --out-root tests/data/recorded
                  [--keep GLOB …] [--strip-json-key KEY …] [--dry-run]
```
- Builds the command with the tool's `*_cmd.py` builder (M2-03, later M2-12…M2-18), runs it with
  `micromamba run -n giap-<tool>` in `$TMPDIR`, copies the kept files to
  `<out-root>/<tool>/<version>[-<db>]/<label>/`, gzips text (`gzip -n`, reproducible), and writes
  `PROVENANCE.md`: command line, env lock hash (`workflow/envs/locks/<tool>.explicit.txt` sha256) or
  image digest, DB name+version, input SHA-256s, date (UTC), host, node, `SLURM_JOB_ID`, wall time.
- Refuses (exit 1) a kept file > 5 MB after gzip, or a total `tests/data/` > 50 MB.
- `--strip-json-key sequences.nt` removes a key (recorded in PROVENANCE) to meet the size cap (Q1).

`scripts/fetch_test_genomes.py`: downloads GCA_000005845.2 and GCA_043228065.1 from
`api.ncbi.nlm.nih.gov` (NCBI Datasets), writes `tests/data/genomes/{k12_mg1655,nissle1917}.fna.gz`
(contig headers as NCBI gives them; circular flags added: `[topology=circular]`) + `PROVENANCE.md`.
Also writes `tests/data/genomes/sim100k.fna.gz` from `SIM100K`/`SIM100K_SEED` (M2-05).

`scripts/slurm/record_bakta_full.sbatch` (16 CPU, 64 GB, 12 h): preflight-normalises each genome
(M1-12 `genome.norm.fa`, `replicons.tsv`), then `record_fixture.py bakta` with the full DB for
`k12`, `nissle`, `sim100k`; keeps JSON, GFF3, GBFF (+ FAA, TSV for sim100k only).
`scripts/slurm/record_fixtures.sbatch` (`--export=TOOLS=a,b,…`; 16 CPU, 64 GB, 8 h): wave A =
`bakta` light on sim100k; `terminatornet` on sim100k and k12; `genomad` on sim100k, k12, nissle;
`antismash` on the full-DB GBFF of sim100k, k12, nissle; `bactermfinder` on sim100k.

Fixture layout produced (wave A):
```
tests/data/genomes/{k12_mg1655,nissle1917,sim100k}.fna.gz, PROVENANCE.md
tests/data/recorded/bakta/1.12.2-light/sim100k/{sim100k.json.gz,.gff3.gz,.gbff.gz,.faa.gz,.tsv.gz,PROVENANCE.md}
tests/data/recorded/bakta/1.12.2-full/{sim100k,k12,nissle}/…
tests/data/recorded/terminatornet/<git12>/{sim100k,k12}/terms.tsv.gz
tests/data/recorded/genomad/1.12.0-db1.9/{sim100k,k12,nissle}/<prefix>_summary/*.tsv.gz
tests/data/recorded/antismash/8.0.4/{sim100k,k12,nissle}/<prefix>.json.gz (trimmed, Q2)
tests/data/recorded/bactermfinder/<git12>/sim100k/…
```

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `record_fixture.py bakta --input genome=g.fa --label sim100k --dry-run` | prints the exact `bakta_command(...)` argv and the target dir `tests/data/recorded/bakta/1.12.2-light/sim100k/`; writes nothing |
| 2 | a kept file 6 MB after gzip | exit 1: `sim.json.gz is 6.0 MB > 5 MB (R11); consider --strip-json-key` |
| 3 | recording twice the same input | identical `.gz` bytes (gzip `-n`, no mtime); PROVENANCE differs only in date, node, job ID, wall time |
| 4 | after `record_bakta_full.sbatch` | `tests/contract/test_bakta_contract.py` and `test_te6_oric.py` lose `pending_sample` and pass |
| 5 | `du -sh tests/data` after wave A | < 50 MB; the PR lists each file's size |

## Edge cases

- NCBI genome headers contain spaces and dots; preflight's safe-ID mapping (M1-02) is applied
  before recording; `contig_map.tsv` is kept beside the fixture.
- A tool that finds nothing (e.g. geNOMAD on sim100k) still yields a valid, committed fixture with
  header-only tables; that is a contract case too.
- Nissle has two plasmids (pMUT1, pMUT2); the replicon table must list them as plasmids.
- Jobs are named `giap-record-*`; never more than 20 GIAP jobs at once (R14).

## Acceptance criteria

- [ ] `tests/unit/test_record_fixture.py`: examples 1–3 with a fake runner.
- [ ] Wave-A fixtures committed with `PROVENANCE.md`; examples 4–5 hold; contract tests of M2-03,
  M2-04 (TE6), M2-15 (TerminatorNet part), M2-18 switch from `pending_sample` to active.
- [ ] Measured wall time and MaxRSS per tool × genome (`sacct`) recorded in `PROGRESS.md` › Measured
  numbers with node class.
- [ ] `PROGRESS.md` "Done" cites the G17 evidence (job IDs, commit).

## Questions for the lead

1. K-12/Nissle full-DB Bakta JSON may exceed 5 MB gzipped. **Proposed default:** strip
   `sequences[].nt` (genome FASTA is committed separately; parser accepts `genome=`), recorded in
   PROVENANCE; if still > 5 MB, also strip per-feature `nt` (re-derived from the genome).
2. antiSMASH full JSON is large. **Proposed default:** keep only `records[].id`, `records[].areas`
   and the CDS entries of `records[].features` that lie inside an area (with their `locus_tag` and
   `gene_kind` qualifiers); drop the record sequence and all other modules. Trimming is done by
   `record_fixture.py --trim antismash` and documented in PROVENANCE.
3. Budget: JSON+GFF3+GBFF for three genomes ≈ 20 MB. **Proposed default:** GBFF only for sim100k
   (antiSMASH input for real genomes comes from the live run in the recording job).

## Lead decisions (after the test-author's pass)
1. Topology follows GIAP_SPEC.md (line ~770): only `circular=true`/`[topology=circular]` header tags (or
   `annotation.bakta.complete: true`) make a replicon circular. The committed fixture genomes
   (tests/data/genomes/*.fna.gz) carry `[topology=circular]` on the chromosome and the circular plasmids;
   the sequences equal NCBI's (sha256 of sequences recorded in PROVENANCE.md). The full-DB Bakta
   recordings are redone from these tagged genomes so replicons are circular and Nissle's plasmids typed
   `plasmid`. Users with plain NCBI FASTA get a preflight Info (spec) — documented in the user guide later.
2. GBFF and FAA are committed for all three genomes (each < 5 MB gzipped; G17), overriding lead decision 4
   (GBFF only for sim100k).
3. BacTermFinder (opt-in, env not built) and record_fixture.py examples 1–3 are deferred (backlog).
