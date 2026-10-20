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


def _db_strain_peptides(res: Any, out: Path) -> tuple[Path, dict[str, list[str]]]:
    """``peptides.tsv`` with each correct peptide replaced by the reference strain's peptide at
    the same offset of the same locus; returns the new path and {database peptide: [sample
    peptides]}."""
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
    back: dict[str, list[str]] = {}
    for sample_pep, db_pep in sub.items():
        back.setdefault(db_pep, []).append(sample_pep)
    return path, back


def inject_short_wrong(
    res: Any, out: Path, k: int, seed: int, q_max: float = 0.01
) -> tuple[Path, set[str]]:
    """Stressed stratum-0 arm (PR #113 review B4; statistician's b2sim ``--mode short``): add ``k``
    wrong database peptides, each a shuffled 7-mer (75 %) or 8-mer window of a random expressed
    protein that does not occur in any protein (I = L), as one confident row (q = 0.001). Every
    target row then gets PEP = the true wrong share among the peptides that pass ``q <= q_max``
    (constant, so the run collapse cannot bias it): π̂ = mean PEP is exact (delta = 0)."""
    import numpy as np

    rng = np.random.default_rng([seed, 31])
    with Path(res.peptides).open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    with Path(res.proteins_faa).open() as fh:
        text = "".join("|" if line.startswith(">") else line.strip() for line in fh)
    prots = [p for p in text.split("|") if p]
    hay = "|".join(p.replace("I", "L") for p in prots)
    tr = read_tsv_df(res.truth["peptides"])
    wrong = set(tr.loc[(tr["is_wrong"] == "true") & (tr["origin"] != "decoy"), "peptide"])
    present = {r["peptide"] for r in rows}
    usable = [p for p in prots if len(p) >= 9]
    inj: set[str] = set()
    tries = 0
    while len(inj) < k and tries < 100 * k:
        tries += 1
        prot = usable[int(rng.integers(len(usable)))]
        n = 7 if rng.random() < 0.75 else 8
        a = int(rng.integers(1, len(prot) - n))
        pep = "".join(rng.permutation(list(prot[a : a + n])).tolist())
        if set(pep) - set("ACDEFGHIKLMNPQRSTVWY") or pep in present or pep in inj:
            continue
        if pep.replace("I", "L") in hay:
            continue
        inj.add(pep)
    template = next(r for r in rows if r["is_decoy"] != "true")
    for j, pep in enumerate(sorted(inj)):
        rows.append({**template, "peptide": pep, "protein_ids": "INJECTED", "q_value": "0.001",
                     "score": "3", "spectrum_id": f"inj.{j + 1}"})  # fmt: skip
    best: dict[str, float] = {}
    for r in rows:
        if r["is_decoy"] != "true":
            best[r["peptide"]] = min(best.get(r["peptide"], 1.0), float(r["q_value"]))
    queries = [p for p, q in best.items() if q <= q_max]
    share = sum(p in wrong or p in inj for p in queries) / max(1, len(queries))
    for r in rows:
        if r["is_decoy"] != "true":
            r["pep"] = f"{share:.10g}"
    path = out / "peptides_injected.tsv"
    out.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return path, inj


def prefilter_to_genome(res: Any, out: Path, pep: float) -> Path:
    """ADR-0019 decision 4 (``input_prefiltered``): the input as an upstream pipeline that kept only
    genome-matching peptides would hand it over: every target row of a wrong identification
    (giapsim truth ``is_wrong``) is dropped, and every target row claims PEP = ``pep`` (giapsim's
    own PEP is ~0.01 at any pi_wrong after the q filter, too small for the warning to have power
    at 50 kb)."""
    with Path(res.peptides).open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    tr = read_tsv_df(res.truth["peptides"])
    wrong = set(tr.loc[tr["is_wrong"] == "true", "peptide"])
    keep = [r for r in rows if r["is_decoy"] == "true" or r["peptide"] not in wrong]
    for r in keep:
        if r["is_decoy"] != "true":
            r["pep"] = f"{pep:.10g}"
    path = out / "peptides_prefiltered.tsv"
    out.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(keep)
    return path


def dbstrain_input(res: Any, tmp: Path) -> tuple[Path, Truth, dict[str, list[str]]]:
    """Database-strain input of ``related_dbstrain`` and its truth (decision 2, review B1): every
    database peptide inherits the origins of the sample peptide(s) it replaces, including pairs
    equal under I = L (an L/I SNV: the mapper places them exactly at the origin); only pairs that
    differ under I = L are ``absent`` (no exact copy at the origin)."""
    truth = database_truth(res)
    path, back = _db_strain_peptides(res, tmp)
    for db_pep, samples in back.items():
        for sample_pep in samples:
            if db_pep == sample_pep:
                continue
            truth.add(db_pep, truth.origins.get(sample_pep, []), truth.genes.get(sample_pep, ()))
            if db_pep.replace("I", "L") != sample_pep.replace("I", "L"):
                truth.absent.add(db_pep)
    return path, truth, back


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
    inject_short: int = 0,
    prefilter: float | None = None,
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
        path, truth, _back = dbstrain_input(res, tmp)
        res_in = dataclasses.replace(res, peptides=path)
    else:
        truth = database_truth(res)
        res_in = res
        if inject_short:
            path, inj = inject_short_wrong(res, tmp, inject_short, seed)
            truth.wrong |= inj
            res_in = dataclasses.replace(res, peptides=path)
        if prefilter is not None:
            res_in = dataclasses.replace(res, peptides=prefilter_to_genome(res, tmp, prefilter))
            truth.wrong = set()
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
    f = pd.read_parquet(ch.avail["mapping_fdr"])
    for dt in ("reverse", "shuffle"):
        row = f[(f["section"] == "placement") & (f["decoy_type"] == dt) & (f["n_diff"] == 0)
                & (f["length_class"] == "all")]  # fmt: skip
        t = int(row["targets"].iloc[0]) if len(row) else 0
        # Ŵ0 = estimated false placements among all stratum-0 targets; V0 = realised (exercise)
        rep.extra[f"W0_{dt}"] = float(row["fdr_stratum"].iloc[0]) * t if t else 0.0
        rep.extra[f"V0_{dt}"] = float(rep.exercise[f"{dt}_m0"])
        rep.extra[f"T0_{dt}"] = float(t)
    if tmp is not None:
        shutil.rmtree(tmp, ignore_errors=True)
    return rep


def run_replicate(kind: str, seed: int, tmp: Path, levels: Sequence[float], **kw: Any) -> Replicate:
    """``run_chain`` + ``score_chain``; the replicate's directory is removed afterwards."""
    ch = run_chain(kind, seed, tmp, **kw)
    return score_chain(ch, levels, tmp)


def mapping_fdr(avail: Mapping[str, Path]) -> pd.DataFrame:
    return pd.read_parquet(avail["mapping_fdr"])
