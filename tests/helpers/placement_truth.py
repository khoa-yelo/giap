"""Truth oracle for the genome-level placement suites (issue M3-21; docs/specs/M3-21 lead
decisions 1-5). Independent of the code under test: it reads giapsim truth files and the mapping
work tables (hit Parquets, 0-based half-open; ``peptide_hits.tsv``, 1-based closed).

* **True origin** (decision 2): correct database peptide -> every truth row of the peptide (its
  ``locus_tag`` copies); de novo read -> the origins of its spectrum's ``true_peptide``, whatever
  the error; other-strain peptide -> the lifted ortholog in the sample (the reference protein's
  offset placed on the sample CDS; SNV-only strains, so coordinates agree); shuffled/random -> none.
* **Overlap**: same contig and strand, same codon phase (start for ``+``, end for ``-``), covering
  >= 50 % of the true interval (circular contigs compared modulo the length).
* **False placement** (decision 1): an accepted peptide none of whose accepted loci (target hits at
  its best cost m_T) overlaps a true origin.
* **Lists**: (decoy type, m_T); each accepted peptide belongs to exactly one, its best-cost list,
  so the pooled FDP (decision 5) is over the union with each peptide once.
* **Exercised** (decision 4, PR #113 review B2): mean >= 1 per replicate of the list's peptides
  (target hit at cost m, any q) that are false, i.e. whose best-cost loci miss every origin. A
  wrong string placed at its true origin (de novo read, other-strain peptide) is not a null
  candidate.
* **Unique guard** (decision 3): a non-false accepted peptide with a non-origin accepted locus must
  not be ``unique`` unless all its loci lie in a true gene.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

Iv = tuple[str, int, int, str]  # contig, start0, end0 (unwrapped), strand
LISTS: tuple[tuple[str, int], ...] = tuple(
    (d, m) for d in ("reverse", "shuffle") for m in (0, 1, 2)
)
POOLED = "pooled"


def list_name(lst: tuple[str, int]) -> str:
    return f"{lst[0]}_m{lst[1]}"


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA"])


@dataclass
class Truth:
    """Origins, true genes and wrong status per query peptide string."""

    lengths: Mapping[str, int]
    circular: Mapping[str, bool]
    origins: dict[str, list[Iv]] = field(default_factory=dict)
    genes: dict[str, set[str]] = field(default_factory=dict)
    wrong: set[str] = field(default_factory=set)
    absent: set[str] = field(default_factory=set)  # correct, but no exact copy at its origin

    def add(self, pep: str, ivs: Iterable[Iv], genes: Iterable[str]) -> None:
        self.origins.setdefault(pep, []).extend(ivs)
        self.genes.setdefault(pep, set()).update(g for g in genes if g)

    # ---------------------------------------------------------------- overlap

    def overlaps(self, hit: Iv, pep: str) -> bool:
        return any(self._ov(hit, o) for o in self.origins.get(pep, ()))

    def _ov(self, h: Iv, o: Iv) -> bool:
        c, s, e, st = h
        tc, ts, te, tst = o
        if c != tc or st != tst or te <= ts:
            return False
        n = self.lengths[c]
        shifts = (-n, 0, n) if self.circular.get(c, True) else (0,)
        for k in shifts:
            a, b = ts + k, te + k
            ov = min(e, b) - max(s, a)
            phase = (s - a) % 3 == 0 if st == "+" else (e - b) % 3 == 0
            if phase and 2 * ov >= (te - ts):
                return True
        return False


def _contigs(genome_fasta: Path) -> dict[str, int]:
    out: dict[str, int] = {}
    name, n = None, 0
    for line in Path(genome_fasta).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name is not None:
        out[name] = n
    return out


def _topology(bakta_json: Path) -> dict[str, bool]:
    doc = json.loads(Path(bakta_json).read_text())
    return {s["id"]: s.get("topology", "circular") == "circular" for s in doc["sequences"]}


def database_truth(res: object) -> Truth:
    """Truth of an M3-19 database table: correct rows' loci; wrong (shuffled/other-strain) rows."""
    lengths = _contigs(res.genome_fasta)  # type: ignore[attr-defined]
    t = Truth(lengths, _topology(res.bakta_json))  # type: ignore[attr-defined]
    tr = _tsv(res.truth["peptides"])  # type: ignore[attr-defined]
    for r in tr[(tr["is_wrong"] == "false") & (tr["origin"] == "gene")].itertuples(index=False):
        t.add(r.peptide, [(r.contig, int(r.start) - 1, int(r.end), r.strand)], [r.locus_tag])
    wrong = tr[(tr["is_wrong"] == "true") & (tr["origin"] != "decoy")]
    t.wrong |= set(wrong["peptide"])
    other = wrong[wrong["origin"] == "other_strain"]
    if len(other):
        src = _source_loci(res.peptides)  # type: ignore[attr-defined]
        for pep in other["peptide"]:
            for loc in src.get(pep, ()):
                t.add(pep, lifted_origin(res, pep, loc), [loc])
    return t


def _source_loci(peptides_tsv: Path) -> dict[str, set[str]]:
    df = _tsv(peptides_tsv)
    df = df[df["is_decoy"] != "true"]
    out: dict[str, set[str]] = {}
    for p, ids in zip(df["peptide"], df["protein_ids"], strict=True):
        out.setdefault(p, set()).update(x for x in str(ids).split(";") if x and x != "nan")
    return out


_REF_CACHE: dict[Path, dict[str, tuple[str, int, int, str, str]]] = {}


def _reference_cds(ref_json: Path) -> dict[str, tuple[str, int, int, str, str]]:
    """locus -> (contig, start0, end0 unwrapped, strand, aa) of the reference (lifted) CDS."""
    ref_json = Path(ref_json)
    if ref_json not in _REF_CACHE:
        doc = json.loads(ref_json.read_text())
        lens = {s["id"]: int(s["length"]) for s in doc["sequences"]}
        out = {}
        for f in doc["features"]:
            if not f.get("aa"):
                continue
            c = f.get("contig", f.get("sequence"))
            s0, stop = int(f["start"]) - 1, int(f["stop"])
            e0 = stop if stop > s0 else stop + lens[c]
            out[str(f["locus"])] = (c, s0, e0, f["strand"], str(f["aa"]))
        _REF_CACHE.clear()
        _REF_CACHE[ref_json] = out
    return _REF_CACHE[ref_json]


def lifted_origin(res: object, pep: str, locus: str) -> list[Iv]:
    """Sample intervals of ``pep`` in the reference protein ``locus`` (SNV-only strains)."""
    cds = _reference_cds(res.reference_json).get(locus)  # type: ignore[attr-defined]
    if cds is None:
        return []
    c, s0, e0, strand, aa = cds
    n = len(pep)
    out: list[Iv] = []
    a = aa.find(pep)
    while a >= 0:
        if strand == "+":
            s = s0 + 3 * a
            e = s + 3 * n
        else:
            e = e0 - 3 * a
            s = e - 3 * n
        out.append((c, s, e, strand))
        a = aa.find(pep, a + 1)
    return out


def denovo_truth(res: object) -> Truth:
    """Truth of the de novo reads: origins of each read's true peptide; wrong = never correct."""
    base = database_truth(res)
    t = Truth(base.lengths, base.circular)
    dn = _tsv(res.truth["denovo"])  # type: ignore[attr-defined]
    ok: dict[str, bool] = {}
    for r in dn.itertuples(index=False):
        t.add(r.peptide, base.origins.get(r.true_peptide, []), base.genes.get(r.true_peptide, ()))
        ok[r.peptide] = ok.get(r.peptide, False) or r.correct == "true"
    t.wrong = {p for p, good in ok.items() if not good}
    return t


# ---------------------------------------------------------------- evaluation


def best_loci(hits: pd.DataFrame, m_target: Mapping[str, int]) -> dict[str, list[Iv]]:
    """Target hits of each peptide at its best cost ``m_target`` (the accepted loci)."""
    h = hits[hits["search"] == "target"]
    h = h[h["peptide"].map(m_target).to_numpy() == h["n_diff"].to_numpy()]
    out: dict[str, list[Iv]] = {}
    for p, c, s, e, st in zip(h["peptide"], h["contig"], h["start"], h["end"], h["strand"],
                              strict=True):  # fmt: skip
        out.setdefault(str(p), []).append((str(c), int(s), int(e), str(st)))
    return out


@dataclass
class Replicate:
    """Per-replicate FDP per list and pooled at each level, plus side measures."""

    fdp: dict[str, np.ndarray]
    exercise: dict[str, int]  # list -> false peptides (best-cost loci miss every origin), any q
    accepted: dict[str, np.ndarray]
    loci_false_share: float  # accepted (peptide, locus) rows at alpha missing every origin
    guard: list[str]  # unique-guard violations (real output)
    guard_mutant: list[str]  # unique-guard violations of the cost <= m_T + 1 mutant
    fdp_mutant: np.ndarray  # pooled FDP at the levels under the mutant's loci
    extra: dict[str, float] = field(default_factory=dict)


def evaluate(
    per: pd.DataFrame,
    hits: pd.DataFrame,
    truth: Truth,
    levels: Sequence[float],
    peptide_hits: pd.DataFrame | None = None,
    alpha: float = 0.01,
    all_target_hits: pd.DataFrame | None = None,
) -> Replicate:
    """Score one replicate: ``per`` is the ``placement`` work table, ``hits`` exact + tolerant.

    ``all_target_hits``: further target loci (e.g. planted near-paralogs at cost 1) that the
    mutant may add; the tolerant matcher reports best-cost hits only.
    """
    per = per.drop_duplicates("peptide")
    m_t = dict(zip(per["peptide"], per["m_target"].astype(int), strict=True))
    loci = best_loci(hits, m_t)
    false = {p: not any(truth.overlaps(iv, p) for iv in loci.get(p, ())) for p in m_t}
    q = per["placement_q"].astype(float).to_numpy()
    peps = per["peptide"].to_numpy()
    key = [f"{d}_m{int(m)}" for d, m in zip(per["decoy_type"], per["m_target"], strict=True)]
    key_arr = np.array(key, dtype=object)
    is_false = np.array([false[p] for p in peps], dtype=bool)
    fdp: dict[str, np.ndarray] = {}
    acc: dict[str, np.ndarray] = {}
    for name in [*(list_name(x) for x in LISTS), POOLED]:
        sel = np.ones(len(peps), bool) if name == POOLED else key_arr == name
        f, a = np.zeros(len(levels)), np.zeros(len(levels))
        for j, lv in enumerate(levels):
            called = sel & ~np.isnan(q) & (q <= lv)
            n = int(called.sum())
            a[j] = n
            f[j] = (called & is_false).sum() / n if n else 0.0
        fdp[name], acc[name] = f, a
    exercise = {list_name(x): int((is_false & (key_arr == list_name(x))).sum()) for x in LISTS}
    rows = [(p, iv) for p in peps[~np.isnan(q) & (q <= alpha)] for iv in loci.get(p, ())]
    lf = sum(not truth.overlaps(iv, p) for p, iv in rows) / len(rows) if rows else 0.0
    g, gm, fm = [], [], np.zeros(len(levels))
    if peptide_hits is not None:
        g = unique_guard(peptide_hits, truth)
        pool = hits if all_target_hits is None else pd.concat([hits, all_target_hits])
        mut = mutant_rows(peptide_hits, pool, m_t)
        gm = unique_guard(mut, truth)
        fm = _mutant_fdp(per, mut, truth, levels, loci)
    return Replicate(fdp, exercise, acc, lf, g, gm, fm)


def _tsv_loci(ph: pd.DataFrame) -> dict[str, list[tuple[Iv, str | None, str]]]:
    out: dict[str, list[tuple[Iv, str | None, str]]] = {}
    for r in ph.itertuples(index=False):
        gene = None if pd.isna(r.gene) or r.gene == "" else str(r.gene)
        iv = (str(r.contig), int(r.start) - 1, int(r.end), str(r.strand))
        out.setdefault(str(r.peptide), []).append((iv, gene, str(r.uniqueness)))
    return out


def unique_guard(ph: pd.DataFrame, truth: Truth) -> list[str]:
    """Decision 3: peptides with an origin locus and a non-origin locus labelled ``unique`` whose
    loci are not all in a true gene (one line per violation)."""
    bad = []
    for pep, rows in _tsv_loci(ph).items():
        on = [truth.overlaps(iv, pep) for iv, _, _ in rows]
        if not any(on) or all(on):
            continue
        if any(u == "unique" for _, _, u in rows):
            genes = truth.genes.get(pep, set())
            if not all(gn is not None and gn in genes for _, gn, _ in rows):
                bad.append(f"{pep}: loci {[iv for iv, _, _ in rows]} true genes {sorted(genes)}")
    return bad


def mutant_rows(ph: pd.DataFrame, hits: pd.DataFrame, m_t: Mapping[str, int]) -> pd.DataFrame:
    """The decision-3 mutant: every accepted peptide also gets its target loci of cost <= m_T + 1
    (as rows of ``peptide_hits.tsv``, 1-based), keeping the peptide's uniqueness label."""
    lab = ph.drop_duplicates("peptide").set_index("peptide")["uniqueness"]
    h = hits[(hits["search"] == "target") & hits["peptide"].isin(lab.index)]
    h = h[h["n_diff"].to_numpy() <= h["peptide"].map(m_t).to_numpy() + 1]
    add = pd.DataFrame(
        {
            "peptide": h["peptide"].astype(str).to_numpy(),
            "contig": h["contig"].astype(str).to_numpy(),
            "start": (h["start"].astype(int) + 1).to_numpy(),
            "end": h["end"].astype(int).to_numpy(),
            "strand": h["strand"].astype(str).to_numpy(),
            "gene": None,
            "uniqueness": h["peptide"].map(lab).to_numpy(),
        }
    )
    base = ph[["peptide", "contig", "start", "end", "strand", "gene", "uniqueness"]]
    out = pd.concat([base, add], ignore_index=True)
    return out.drop_duplicates(["peptide", "contig", "start", "end", "strand"], keep="first")


def _mutant_fdp(
    per: pd.DataFrame,
    mut: pd.DataFrame,
    truth: Truth,
    levels: Sequence[float],
    loci: Mapping[str, list[Iv]],
) -> np.ndarray:
    """Pooled FDP when the mutant's extra loci count as accepted loci (peptides at q <= alpha
    get the mutant rows; the others keep their best-cost loci)."""
    mloci: dict[str, list[Iv]] = {}
    for p, rows in _tsv_loci(mut).items():
        mloci[p] = [iv for iv, _, _ in rows]
    q = per["placement_q"].astype(float).to_numpy()
    peps = per["peptide"].to_numpy()
    out = np.zeros(len(levels))
    for j, lv in enumerate(levels):
        called = [p for p, qq in zip(peps, q, strict=True) if qq == qq and qq <= lv]
        if called:
            n_false = sum(
                not any(truth.overlaps(iv, p) for iv in mloci.get(p, loci.get(p, ())))
                for p in called
            )
            out[j] = n_false / len(called)
    return out
