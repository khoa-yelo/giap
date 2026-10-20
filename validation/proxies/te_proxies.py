"""TE1 and TE8 in-silico proxies (GIAP_SPEC.md §11.6, §13 gate M3; issue M3-22).

Both proxies digest proteins with :func:`giap.mapping.digest.digest` (trypsin_p, missed <= 1,
7-30 aa) and send the distinct peptide sequences (lead decision 4) through the merged mapping
chain (``giap preflight`` normalisation, then the registered steps under the local executor) as a
``database`` input without q-values or PEPs, so every decoy hit weighs 1 (pure target-decoy).

* **TE1** (``search_db = same``): Nissle Bakta CDS/sORF proteins (pseudogene-flagged CDSs
  included) on the Nissle genome. A peptide passes when an accepted ``route = exact`` row has
  ``locus_class`` in {cds_inframe, sorf} at one of the genes it was digested from, or a class that
  outranks them (``rna_overlap``, ``pseudogene``) on a locus in frame inside such a gene
  (``exact_inframe_masked``). N-terminal products of strict non-ATG genes (first genome
  codon not ATG and the Bakta protein starts with M; decision 5) are counted as
  ``non_atg_nterm`` and leave the denominator. Pass: share >= 0.999.
* **TE8** (``search_db = other``): UniProt K-12 (UP000000625) peptides. ``k12_exact``: peptides
  with an exact target hit on the K-12 genome, after removing ``non_atg_nterm`` (decisions 6, 9:
  N-terminal in a UniProt protein and equal to an N-terminal product of a strict non-ATG K-12
  gene). On Nissle: accepted exact / tolerant placements, rejected placements, peptides without
  any target hit, and the pooled decoy-estimated FDR among accepted placements. Pass:
  ``k12_exact >= 0.995`` and FDR <= 0.01.

Usage (writes ``validation/proxies/<YYYY-MM-DD>_te1.json`` / ``_te8.json`` by default)::

    python validation/proxies/te_proxies.py te1 GENOME BAKTA_JSON [--out FILE]
    python validation/proxies/te_proxies.py te8 K12_GENOME K12_JSON NISSLE_GENOME NISSLE_JSON
        UNIPROT_FASTA [--out FILE]
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
import tempfile
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from giap.annotate.bakta import parse_bakta_json
from giap.config import GiapConfig, resolve_config
from giap.io.fasta import read_genome_fasta
from giap.io.protein_fasta import read_protein_fasta
from giap.io.tsv import write_tsv
from giap.mapping.chance import residue_freqs
from giap.mapping.digest import digest
from giap.mapping.placement_fdr import DECOY_TYPE_OF, STRATA, best_costs
from giap.mapping.sixframe import load_index
from giap.model.coords import Interval
from giap.model.genes import Gene
from giap.steps.mapping import select_queries

TE1_MIN, TE8_K12_MIN, TE8_FDR_MAX = 0.999, 0.995, 0.01
MAX_EXCEPTIONS = 50
INFRAME = frozenset({"cds_inframe", "sorf"})
MASKING = frozenset({"rna_overlap", "pseudogene"})  # classes ranked above sorf/cds_inframe
HERE = Path(__file__).resolve().parent


@dataclass(frozen=True)
class ProxyResult:
    """Counts, shares, the pass flag and up to 50 example exceptions of one proxy."""

    name: str
    counts: dict[str, int]
    shares: dict[str, float]
    passed: bool
    exceptions: list[dict[str, Any]] = field(default_factory=list)

    def to_json(self) -> dict[str, Any]:
        return {
            "proxy": self.name,
            "passed": self.passed,
            "counts": self.counts,
            "shares": self.shares,
            "exceptions": self.exceptions,
        }


# ---------------------------------------------------------------- pure helpers


def coding_genes(genes: Iterable[Gene]) -> list[Gene]:
    """Bakta CDS/sORF proteins: genes of type cds, sorf or pseudogene (a Bakta ``cds`` with a
    ``pseudogene`` record) that have a protein sequence."""
    return [g for g in genes if g.feature_type in ("cds", "sorf", "pseudogene") and g.aa]


def is_strict_non_atg(g: Gene) -> bool:
    """Lead decision 5: the first genome codon is not ATG and the protein starts with M."""
    return bool(g.nt) and str(g.nt)[:3].upper() != "ATG" and g.aa.startswith("M")


def nterm_products(protein: str) -> list[str]:
    """Digest products starting at position 0."""
    return [p for s, _, p in digest(protein) if s == 0]


def peptide_sources(proteins: Iterable[tuple[str, str]]) -> dict[str, set[str]]:
    """Distinct digest sequence -> IDs of the proteins it comes from."""
    out: dict[str, set[str]] = {}
    for pid, seq in proteins:
        for _, _, p in digest(seq):
            out.setdefault(p, set()).add(pid)
    return out


def non_atg_nterm(genes: Iterable[Gene]) -> dict[str, str]:
    """N-terminal product of a strict non-ATG gene -> its (smallest) locus tag."""
    out: dict[str, str] = {}
    for g in sorted(genes, key=lambda x: x.locus_tag):
        if is_strict_non_atg(g):
            for p in nterm_products(g.aa):
                out.setdefault(p, g.locus_tag)
    return out


def _exceptions(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for r in rows:
        if len(out) >= MAX_EXCEPTIONS:
            break
        out.append(r)
    return out


def inframe_inside(
    contig: str, start: int, end: int, strand: str, gene: Interval, circular_len: int | None
) -> bool:
    """``[start, end)`` (0-based, unwrapped) lies inside ``gene`` on its strand and frame; on a
    circular replicon the hit may be shifted by +-L (as in ``giap.mapping.loci``)."""
    if gene.contig != contig or gene.strand != strand:
        return False
    fs, fe = gene.start, gene.end
    if circular_len and fe <= fs:
        fe += circular_len
    shifts = (0, circular_len, -circular_len) if circular_len else (0,)
    for d in shifts:
        a, b = start + d, end + d
        if fs <= a and b <= fe and ((a - fs) % 3 == 0 if strand == "+" else (fe - b) % 3 == 0):
            return True
    return False


def te1_result(
    sources: Mapping[str, Collection[str]],
    excluded: Mapping[str, str],
    peptide_hits: pd.DataFrame,
    genes: Mapping[str, Interval],
    circular: Mapping[str, int | None],
) -> ProxyResult:
    """TE1 counts from the peptide sources, the non-ATG exclusions and ``peptide_hits.tsv``.

    A peptide is ``exact_inframe`` when an accepted exact row has ``locus_class`` cds_inframe or
    sorf at one of its own genes, or a higher-priority class (``rna_overlap``, ``pseudogene``,
    §7.2.6) on a locus that lies in frame inside one of its own genes; the latter are also
    counted as ``exact_inframe_masked`` (diagnostic). ``circular``: contig -> length when
    circular, else None.
    """
    good: set[str] = set()
    masked: set[str] = set()
    seen: dict[str, list[str]] = {}
    if len(peptide_hits):
        cols = ("peptide", "route", "locus_class", "gene", "contig", "start", "end", "strand")
        for pep, route, cls, gene, contig, s1, e1, strand in zip(
            *(peptide_hits[c] for c in cols), strict=True
        ):
            pep = str(pep)
            seen.setdefault(pep, []).append(f"{route}:{cls}:{gene}")
            own = sources.get(pep, ())
            if route != "exact" or pep in good:
                continue
            if cls in INFRAME and gene in own:
                good.add(pep)
            elif cls in MASKING and any(
                inframe_inside(str(contig), int(s1) - 1, int(e1), str(strand), genes[t],
                               circular.get(str(contig)))
                for t in own
                if t in genes
            ):  # fmt: skip
                masked.add(pep)
    masked -= good
    good |= masked
    peps = sorted(sources)
    n_ex = sum(p in excluded for p in peps)
    exact = sum(p in good for p in peps if p not in excluded)
    misses = [p for p in peps if p not in excluded and p not in good]
    denom = len(peps) - n_ex
    share = exact / denom if denom else 0.0
    rows = [
        {"kind": "miss", "locus_tag": min(sources[p]), "peptide": p,
         "hits": ";".join(seen.get(p, [])) or "unplaced"}
        for p in misses
    ]  # fmt: skip
    rows += [
        {"kind": "non_atg_nterm", "locus_tag": excluded[p], "peptide": p}
        for p in peps
        if p in excluded
    ]
    counts = {
        "peptides": len(peps),
        "non_atg_nterm": n_ex,
        "exact_inframe": exact,
        "misses": len(misses),
        "exact_inframe_masked": sum(p in masked for p in peps if p not in excluded),
    }
    return ProxyResult(
        "TE1", counts, {"exact_inframe": share}, bool(denom) and share >= TE1_MIN,
        _exceptions(rows),
    )  # fmt: skip


def pooled_decoy_fdr(costs: pd.DataFrame, per: pd.DataFrame) -> float:
    """Pooled decoy-estimated FDR among accepted placements (unit decoy weights).

    Per decoy type and stratum m with accepted targets, the accepted list is E <= e* (the largest
    accepted E); its estimate is (D_m(e*) + 1) / A_m, D counting decoys with best cost m and
    E <= e* (§7.2.5 with pi_wrong = 1, the input carries no q-value or PEP). Pooled:
    sum_m (D_m(e*) + 1) / sum_m A_m; 0 when nothing is accepted.
    """
    acc = per[per["accepted"].astype(bool)]
    if not len(acc):
        return 0.0
    c = costs.set_index("peptide")
    src = c["source"].astype(str)
    false = 0.0
    for source, dtype in DECOY_TYPE_OF.items():
        for m in STRATA:
            a = acc[(acc["decoy_type"] == dtype) & (acc["m_target"] == m)]
            if not len(a):
                continue
            e_star = float(np.nanmax(a["e_target"].to_numpy(dtype=float)))
            d = c[(src == source) & (c["m_decoy"] == m) & (c["e_decoy"] <= e_star)]
            false += len(d) + 1.0
    return min(1.0, false / len(acc))


# ---------------------------------------------------------------- the mapping chain (file I/O)


def _write_peptides(path: Path, sources: Mapping[str, Collection[str]], search_db: str) -> None:
    rows = [(p, ";".join(sorted(sources[p])), "database", search_db) for p in sorted(sources)]
    write_tsv(path, ["peptide", "protein_ids", "source", "search_db"], rows)


def run_chain(
    cfg: GiapConfig,
    genome: Path,
    bakta_json: Path,
    sources: Mapping[str, Collection[str]],
    search_db: Literal["same", "other"],
    out: Path,
    targets: Sequence[str],
) -> Path:
    """Preflight + the local executor up to ``targets``; returns the work directory."""
    from giap.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs
    from giap.runner.local import run_local
    from giap.steps import load_step_modules

    tsv = out / "peptides.tsv"
    out.mkdir(parents=True, exist_ok=True)
    _write_peptides(tsv, sources, search_db)
    pf = out / "pf"
    report, norm = run_preflight(
        PreflightInputs(genome, None, tsv, None, bakta_json), cfg, None, "preflight", pf
    )
    write_outputs(pf, report, norm)
    if exit_code(report, False) != 0:
        errors = [f"{c.rule_id}: {c.message}" for c in report.checks if c.level == "error"]
        raise RuntimeError(f"preflight failed for {genome}: {errors}")
    work = pf / "work" / "preflight"
    roots = {
        "genome": work / "genome.norm.fa",
        "peptides": work / "peptides.norm.parquet",
        "bakta_json": bakta_json,
    }
    load_step_modules()
    run_local(cfg, out / "work", targets, roots)
    return out / "work"


def read_coding_genes(bakta_json: Path, genome: Path) -> tuple[list[Gene], dict[str, int | None]]:
    """CDS/sORF genes of a Bakta JSON (nucleotides from ``genome`` when the JSON has none) and
    contig -> length for circular replicons (None for linear ones)."""
    seqs = {r.id: r.sequence for r in read_genome_fasta(genome)}
    model = parse_bakta_json(bakta_json, genome=seqs).model
    circ = {c: (r.length if r.topology == "circular" else None) for c, r in model.replicons.items()}
    return coding_genes(model.genes), circ


def _read_hits_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA"])


def te1_proxy(genome: Path, bakta_json: Path) -> ProxyResult:
    """TE1: Nissle Bakta CDS/sORF peptides map exactly in frame at their own gene (>= 0.999)."""
    genes, circular = read_coding_genes(bakta_json, genome)
    sources = peptide_sources((g.locus_tag, g.aa) for g in genes)
    excluded = non_atg_nterm(genes)
    cfg = resolve_config()
    with tempfile.TemporaryDirectory(prefix="giap-te1-") as tmp:
        work = run_chain(cfg, genome, bakta_json, sources, "same", Path(tmp), ["loci"])
        hits = _read_hits_tsv(work / "mapping" / "peptide_hits.tsv")
    return te1_result(sources, excluded, hits, {g.locus_tag: g.interval for g in genes}, circular)


def te8_result(
    sources: Mapping[str, Collection[str]],
    excluded: Collection[str],
    k12_exact_hits: Collection[str],
    per: pd.DataFrame,
    costs: pd.DataFrame,
) -> ProxyResult:
    """TE8 counts from the K-12 exact-hit set and the Nissle placement table and costs."""
    peps = sorted(sources)
    n = len(peps)
    ex = {p for p in peps if p in excluded}
    hit = set(k12_exact_hits)
    k12 = sum(p in hit for p in peps if p not in ex)
    denom = n - len(ex)
    acc = per[per["accepted"].astype(bool)]
    m = acc["m_target"].to_numpy()
    counts = {
        "peptides": n,
        "non_atg_nterm": len(ex),
        "k12_exact": k12,
        "nissle_exact": int((m == 0).sum()),
        "nissle_tolerant": int((m > 0).sum()),
        "nissle_rejected": int(len(per) - len(acc)),
        "nissle_unmapped": int(n - len(per)),
    }
    fdr = pooled_decoy_fdr(costs, per)
    shares = {
        "k12_exact": k12 / denom if denom else 0.0,
        "nissle_exact": counts["nissle_exact"] / n if n else 0.0,
        "nissle_tolerant": counts["nissle_tolerant"] / n if n else 0.0,
        "nissle_unmapped": counts["nissle_unmapped"] / n if n else 0.0,
        "nissle_placement_fdr": fdr,
    }
    rows = [
        {"kind": "k12_miss", "locus_tag": min(sources[p]), "peptide": p}
        for p in peps
        if p not in ex and p not in hit
    ]
    passed = bool(denom) and shares["k12_exact"] >= TE8_K12_MIN and fdr <= TE8_FDR_MAX
    return ProxyResult("TE8", counts, shares, passed, _exceptions(rows))


def te8_proxy(
    k12_genome: Path,
    k12_json: Path,
    nissle_genome: Path,
    nissle_json: Path,
    uniprot_fasta: Path,
) -> ProxyResult:
    """TE8: UniProt K-12 peptides on K-12 (exact share) and on Nissle (placement FDR)."""
    prots = [(r.id, r.sequence) for r in read_protein_fasta(uniprot_fasta)]
    sources = peptide_sources(prots)
    uni_nterm = {p for _, seq in prots for p in nterm_products(seq)}
    k12_genes, _ = read_coding_genes(k12_json, k12_genome)
    excluded = uni_nterm & set(non_atg_nterm(k12_genes))  # decision 9: gene-model reading
    cfg = resolve_config()
    with tempfile.TemporaryDirectory(prefix="giap-te8-") as tmp:
        k_work = run_chain(cfg, k12_genome, k12_json, sources, "other", Path(tmp) / "k12",
                           ["map_exact"])  # fmt: skip
        hx = pd.read_parquet(k_work / "hits_exact")
        k12_hit = set(hx.loc[hx["search"] == "target", "peptide"].astype(str))
        n_work = run_chain(cfg, nissle_genome, nissle_json, sources, "other",
                           Path(tmp) / "nissle", ["placement"])  # fmt: skip
        per = pd.read_parquet(n_work / "placement")
        costs = _costs(n_work)
    return te8_result(sources, excluded, k12_hit, per, costs)


def _costs(work: Path) -> pd.DataFrame:
    """Best target/decoy costs of the placement step's queries (``best_costs`` on its inputs)."""
    queries = select_queries(pd.read_parquet(work / "peptides"), resolve_config())
    parts = [pd.read_parquet(work / k) for k in ("hits_exact", "hits_tolerant")]
    hits = pd.concat([p for p in parts if len(p)] or parts[:1], ignore_index=True)
    index, _ = load_index(work / "sixframe")
    return best_costs(hits, queries, residue_freqs(index))


# ---------------------------------------------------------------- CLI


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="proxy", required=True)
    p1 = sub.add_parser("te1")
    p1.add_argument("genome", type=Path)
    p1.add_argument("bakta_json", type=Path)
    p8 = sub.add_parser("te8")
    for name in ("k12_genome", "k12_json", "nissle_genome", "nissle_json", "uniprot_fasta"):
        p8.add_argument(name, type=Path)
    for p in (p1, p8):
        p.add_argument("--out", type=Path, default=None)
    a = ap.parse_args(argv)
    if a.proxy == "te1":
        r = te1_proxy(a.genome, a.bakta_json)
    else:
        r = te8_proxy(a.k12_genome, a.k12_json, a.nissle_genome, a.nissle_json, a.uniprot_fasta)
    out = a.out or HERE / f"{dt.date.today().isoformat()}_{a.proxy}.json"
    out.write_text(json.dumps(r.to_json(), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{r.name}: passed={r.passed} counts={r.counts} shares={r.shares} -> {out}")
    return 0 if r.passed else 1


if __name__ == "__main__":
    sys.exit(main())
