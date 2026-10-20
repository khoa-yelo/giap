# M3-22: In-silico digestion and the TE1 / TE8 proxies

Milestone M3 · depends on: M3-17 (mapping chain), M3-01; recorded Nissle and K-12 Bakta 1.12.2
full-DB JSON and the committed genomes in `tests/data/genomes/` (both M2-10, G17, §15.2 step 1) ·
creates `src/giap/mapping/digest.py`, `validation/proxies/te_proxies.py`,
`validation/fetch_uniprot_k12.py`, an entry in `scripts/slurm/weekly.sbatch` · pure Python
(network for UniProt on FarmShare) · estimate ≈ 260 code lines

## Spec sections implemented

- §11.6 "TE1 | ≥ 99.5 % of A2 peptides map exactly, in frame, to annotated CDS | pilot (M10).
  In-silico proxy: in-silico tryptic peptides of Nissle Bakta proteins map ≥ 99.9 % exactly in
  frame, the exceptions being N-terminal peptides at non-ATG starts".
- "TE8 | cross-strain sanity: UniProt K-12 tryptic peptides on the K-12 genome map ≥ 99.5 %
  exactly; on the Nissle genome, the decoy-estimated placement FDR at defaults is ≤ 1 % and
  tolerant plus unmapped shares are reported | FarmShare (needs `rest.uniprot.org`)".
- §13 Gate M3 "TE1 and TE8 proxies pass". §12 weekly suite "in-silico peptide sets; TE6–TE8
  proxies". §11.4 "UniProt bulk downloads (> 5 MB) … Scripts download them into
  `$GIAP_DATA_DIR`." §15.3 UP000000625.
- `abundance.enzyme: trypsin_p` (the digest is reused by M4 for `n_observable_peptides` and
  `abundance.observable_len`).

## Public interface

```python
def digest(protein: str, enzyme: Literal["trypsin_p", "trypsin"] = "trypsin_p", missed: int = 1,
           min_len: int = 7, max_len: int = 30) -> list[tuple[int, int, str]]: ...  # (start, end, peptide), 0-based
def te1_proxy(genome: Path, bakta_json: Path) -> ProxyResult: ...   # share exact in frame, exceptions listed
def te8_proxy(k12_genome: Path, k12_json: Path, nissle_genome: Path, nissle_json: Path,
              uniprot_fasta: Path) -> ProxyResult: ...
```
`ProxyResult` holds counts, shares, the pass/fail flag, and up to 50 example exceptions; the
script writes `validation/proxies/<date>_te1.json` / `_te8.json` (dates only in validation
outputs, never in `mapping/`).

**TE1 rule.** Peptides = `digest` of every Bakta CDS/sORF protein of Nissle (missed ≤ 1, 7–30
aa), as a `database`, `search_db = same` input. Pass: share with `route = exact` and
`locus_class` ∈ {`cds_inframe`, `sorf`} at their own gene ≥ 0.999 after excluding the
N-terminal peptides of genes whose start codon is not ATG (Bakta writes `M`, the genome
translates `V`/`L`). Those exceptions are counted and listed.
**TE8 rule.** UniProt K-12 (UP000000625) digest as `database`, `search_db = other` input. On
K-12: exact share ≥ 0.995. On Nissle: the pooled decoy-estimated FDR among accepted placements
at defaults ≤ 0.01; tolerant and unmapped shares reported.

## Examples

| # | Input | Expected |
|---|---|---|
| 1 | `digest("MKAAKPRGGRLLLK")` (trypsin_p, missed 1, min_len 3 for the test) | cuts after K1, K4, R6, R9; products include `AAKPR` (missed cleavage at K4 kept with P) |
| 2 | same with `enzyme="trypsin"` | no cut at K4 (followed by P) |
| 3 | `digest(p, missed=0)` of a 300-aa protein | products tile the protein between cut sites; none > 30 or < 7 aa |
| 4 | TE1 on the 100-kb sim genome (all ATG starts) | share = 1.0 |
| 5 | TE1 on a sim genome with 8 % GTG starts | the only misses are N-terminal peptides of GTG/TTG genes |
| 6 | TE8 without network | test skipped by `needs_net(rest.uniprot.org)` |

## Acceptance criteria

- [ ] `tests/unit/test_digest.py`: examples 1–3; property: concatenating missed = 0 products with
  the short fragments reconstructs the protein.
- [ ] `tests/integration/test_te1_proxy.py` (`slow`): examples 4–5 on sims; the real Nissle TE1
  proxy runs when the recorded full-DB JSON is present (otherwise `pending_sample`), and passes.
- [ ] TE8 proxy runs in `weekly.sbatch` on FarmShare (needs `rest.uniprot.org`); its result JSON is
  linked from `PROGRESS.md` as gate evidence.
- [ ] `validation/fetch_uniprot_k12.py` downloads UP000000625 into `$GIAP_DATA_DIR/uniprot/` with a
  SHA-256 and retrieval date; never committed (R11).

## Questions for the lead

1. TE1 "exceptions being N-terminal peptides at non-ATG starts": also except peptides that start
   with the initiator Met removed (Met excision)? **Proposed default:** the digest keeps the
   protein as Bakta writes it (with `M`); peptides starting after Met excision are not generated,
   so only the non-ATG case arises.
2. If the recorded Nissle full-DB JSON is not ready when M3 closes: **Proposed default:** use the
   light-DB Bakta JSON for TE1 (CDS calls are identical between DBs for coordinates), and re-run
   once G17 fixtures exist.

## Lead decisions (2026-10-20)
1. Accept the default: the digest keeps the protein as Bakta writes it (with the initiator `M`); no
   Met-excised peptides are generated, so only the non-ATG-start exception arises for TE1.
2. The recorded full-DB Bakta JSONs exist (`tests/data/recorded/bakta/1.12.2-full/{nissle,k12}`): use
   the full-DB recording for TE1; the light-DB fallback is not needed.
3. TE8 (UniProt K-12) needs the network: it runs only in the weekly job (marker `e2e`/`slow` as the
   spec says); unit and CI tests use a small recorded excerpt with PROVENANCE.md, never the network.
