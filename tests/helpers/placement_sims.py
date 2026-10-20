"""One genome-level placement replicate through the registered mapping steps (issue M3-21).

``run_replicate(kind, seed, tmp, ...)`` simulates a giapsim genome with M3-19 proteomics (and M3-18
strain / M3-20 de novo reads when the suite needs them), runs ``giap preflight`` and the mapping
chain (``helpers.mapping_steps.drive``, test ``$GIAP_DB_DIR``, no tools), and scores it with
``helpers.placement_truth.evaluate``. Kinds:

* ``database``: M3-19 database errors (``wrong_source``), ``search_db = other``.
* ``related_dbstrain``: as ``database`` with a strain, but every *correct* peptide is reported as
  the database (reference) strain's sequence at the same protein offset (what a search against the
  related proteome identifies); its origin is the sample locus of that offset (decision 2).
* ``denovo``: the M3-20 reads of the same spectra as de novo input (shuffle decoys) with the
  reference peptides; ``decoys`` overrides ``mapping.decoys.denovo_input`` (attack 5 documentation).
"""

from __future__ import annotations

import csv
import dataclasses
import shutil
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from giap.config import resolve_config
from helpers.mapping_steps import (
    drive,
    fresh_steps,
    module2_files,
    preflight,
    read_tsv_df,
    roots_of,
)
from helpers.placement_truth import (
    Replicate,
    Truth,
    _reference_cds,
    database_truth,
    denovo_truth,
    evaluate,
)


@dataclasses.dataclass
class Chain:
    res: Any
    avail: dict[str, Path]
    truth: Truth
    cfg: Any
    planted: dict[str, tuple[str, int, int, str]] = dataclasses.field(default_factory=dict)


def _db_strain_peptides(res: Any, out: Path) -> tuple[Path, dict[str, str]]:
    """``peptides.tsv`` with each correct peptide replaced by the reference strain's peptide at
    the same offset of the same locus; returns the new path and {database peptide: sample
    peptide}."""
    tr = read_tsv_df(res.truth["peptides"])
    cds = _reference_cds(res.reference_json)
    sub: dict[str, str] = {}
    for r in tr[(tr["is_wrong"] == "false") & (tr["origin"] == "gene")].itertuples(index=False):
        if r.peptide in sub or r.locus_tag not in cds:
            continue
        _c, s0, e0, strand, aa = cds[r.locus_tag]
        n = len(r.peptide)
        a = (int(r.start) - 1 - s0) // 3 if strand == "+" else (e0 - int(r.end)) // 3
        ref = aa[a : a + n] if a >= 0 else ""
        sub[r.peptide] = ref if len(ref) == n else r.peptide
    with Path(res.peptides).open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        if row["is_decoy"] != "true" and row["peptide"] in sub:
            row["peptide"] = sub[row["peptide"]]
    path = out / "peptides_dbstrain.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path, {v: k for k, v in sub.items()}


_COMP = str.maketrans("ACGT", "TGCA")
_ALT = {"W": "TGG", "A": "GCT"}  # replacement codon: W, or A when the residue is already W


def _revcomp(nt: str) -> str:
    return nt.translate(_COMP)[::-1]


def plant_paralogs(res: Any, n: int, seed: int) -> dict[str, tuple[str, int, int, str]]:
    """Copy ``n`` single-locus correct peptides (10-20 aa) into free intergenic space with one
    substituted middle codon, in both ``genome.fasta`` and the Bakta-like JSON: each gets a
    cost-1 non-origin locus (a near-paralog). Returns {peptide: copy interval (0-based)}."""
    import json

    import numpy as np

    rng = np.random.default_rng([seed, 21])
    fasta = Path(res.genome_fasta).read_text().splitlines()
    head, seq = fasta[0], "".join(fasta[1:])
    contig = head[1:].split()[0]
    feats = read_tsv_df(res.truth["features"])
    feats = feats[feats["contig"] == contig]
    occupied = sorted((int(a) - 1, int(b)) for a, b in zip(feats["start"], feats["end"],
                                                          strict=True))  # fmt: skip
    gaps, last = [], 0
    for a, b in occupied:
        if a - last >= 90:
            gaps.append((last + 15, a - 15))
        last = max(last, b)
    tr = read_tsv_df(res.truth["peptides"])
    cand = tr[
        (tr["is_wrong"] == "false")
        & (tr["n_loci_genome"] == "1")
        & (tr["contig"] == contig)
        & tr["peptide"].str.len().between(10, 20)
        & (tr["alt_start_nterm"] == "false")
    ]
    cand = cand.drop_duplicates("peptide")
    cand = cand[cand["end"].astype(int) <= len(seq)]
    order = rng.permutation(len(cand)).tolist()
    out: dict[str, tuple[str, int, int, str]] = {}
    s_list = list(seq)
    gi = rng.permutation(len(gaps)).tolist()
    for i, g in zip(order, gi, strict=False):
        if len(out) == n:
            break
        r = cand.iloc[i]
        s0, e0 = int(r["start"]) - 1, int(r["end"])
        nt = seq[s0:e0] if r["strand"] == "+" else _revcomp(seq[s0:e0])
        k = 3 * (len(r["peptide"]) // 2)
        nt = nt[:k] + _ALT["A" if r["peptide"][k // 3] == "W" else "W"] + nt[k + 3 :]
        a, b = gaps[g]
        if b - a < len(nt):
            continue
        put = nt if r["strand"] == "+" else _revcomp(nt)
        s_list[a : a + len(put)] = list(put)
        out[str(r["peptide"])] = (contig, a, a + len(put), str(r["strand"]))
    seq = "".join(s_list)
    Path(res.genome_fasta).write_text(
        head + "\n" + "".join(seq[j : j + 60] + "\n" for j in range(0, len(seq), 60))
    )
    doc = json.loads(Path(res.bakta_json).read_text())
    for rec in doc["sequences"]:
        if rec["id"] == contig:
            rec["nt" if "nt" in rec else "sequence"] = seq
    Path(res.bakta_json).write_text(json.dumps(doc, sort_keys=True, indent=1) + "\n")
    return out


def run_chain(
    kind: str,
    seed: int,
    tmp: Path,
    *,
    length: int = 100_000,
    proteomics: Any = None,
    strain: Any = None,
    denovo: Any = None,
    overrides: Sequence[str] = (),
    paralogs: int = 0,
) -> Chain:
    from giapsim.simulate import SimConfig, simulate

    cfg_sim = SimConfig(
        genome_length=length,
        proteomics=proteomics,
        strain=strain,
        denovo=denovo,
        write_precomputed=True,
    )
    res = simulate(cfg_sim, seed, tmp / "sim")
    planted = plant_paralogs(res, paralogs, seed) if paralogs else {}
    truth: Truth
    if kind == "denovo":
        truth = denovo_truth(res)
        res_in = dataclasses.replace(res, peptides=res.outdir / "peptides_denovo.tsv")
    elif kind == "related_dbstrain":
        truth = database_truth(res)
        path, back = _db_strain_peptides(res, tmp)
        for db_pep, sample_pep in back.items():
            if db_pep.replace("I", "L") != sample_pep.replace("I", "L"):
                truth.absent.add(db_pep)
                truth.add(
                    db_pep, truth.origins.get(sample_pep, []), truth.genes.get(sample_pep, ())
                )
        res_in = dataclasses.replace(res, peptides=path)
    else:
        truth = database_truth(res)
        res_in = res
    run = preflight(res_in, tmp / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={seed}", *overrides])
    roots = roots_of(run)
    if kind == "denovo":
        roots["reference_peptides"] = res.outdir / "reference_peptides.tsv"
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp / "work", roots)
    return Chain(res, avail, truth, cfg, planted)


def null_share_stratum0(ch: Chain) -> tuple[float, float]:
    """(π̂ of the database stratum-0 list, true share of its searched peptides without an exact
    origin: wrong or absent from the sample). Database input, ``search_db = other`` only."""
    from giap.steps.mapping import select_queries

    q = select_queries(pd.read_parquet(ch.avail["peptides"]), ch.cfg)
    q = q[q["source"] == "database"]
    peps = q["peptide"].astype(str)
    null = peps.isin(ch.truth.wrong | ch.truth.absent).to_numpy()
    f = pd.read_parquet(ch.avail["mapping_fdr"])
    row = f[(f["section"] == "placement") & (f["decoy_type"] == "reverse") & (f["n_diff"] == 0)
            & (f["length_class"] == "all")].iloc[0]  # fmt: skip
    return float(row["pi_wrong"]), float(null.mean()) if len(null) else 0.0


def hits_of(avail: Mapping[str, Path]) -> pd.DataFrame:
    parts = [pd.read_parquet(avail[k]) for k in ("hits_exact", "hits_tolerant")]
    return pd.concat([p for p in parts if len(p)] or parts[:1], ignore_index=True)


def planted_hits(ch: Chain) -> pd.DataFrame:
    """The planted near-paralogs as cost-1 target hits (0-based), the mutant's extra loci."""
    rows = [
        {"peptide": p, "search": "target", "contig": c, "start": a, "end": b, "strand": st,
         "n_diff": 1}
        for p, (c, a, b, st) in ch.planted.items()
    ]  # fmt: skip
    return pd.DataFrame(rows, columns=["peptide", "search", "contig", "start", "end", "strand",
                                       "n_diff"])  # fmt: skip


def score_chain(ch: Chain, levels: Sequence[float], tmp: Path | None = None) -> Replicate:
    per = pd.read_parquet(ch.avail["placement"])
    files = module2_files(ch.avail["mapping"].parent, ch.avail)
    ph = read_tsv_df(files["peptide_hits.tsv"])
    extra = planted_hits(ch)
    rep = evaluate(
        per, hits_of(ch.avail), ch.truth, levels, ph, ch.cfg.stats.alpha.placement, extra
    )
    if tmp is not None:
        shutil.rmtree(tmp, ignore_errors=True)
    return rep


def run_replicate(kind: str, seed: int, tmp: Path, levels: Sequence[float], **kw: Any) -> Replicate:
    """``run_chain`` + ``score_chain``; the replicate's directory is removed afterwards."""
    ch = run_chain(kind, seed, tmp, **kw)
    return score_chain(ch, levels, tmp)


def mapping_fdr(avail: Mapping[str, Path]) -> pd.DataFrame:
    return pd.read_parquet(avail["mapping_fdr"])
