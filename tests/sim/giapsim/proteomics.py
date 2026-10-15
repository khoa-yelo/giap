"""Simulated proteomics: proteins, tryptic peptides and database-search errors (issue M3-19).

GIAP_SPEC.md §11.3 "Expression and proteomics". Expressed proteins are the *true* proteins of the
sample genome (non-pseudo CDS and sORFs from their true start, plus the unannotated novel ORFs);
each distinct Trypsin/P product is detected once by a logistic draw in its length and the
protein's log abundance and then identified in every run. Wrong identifications (a share
``pi_wrong`` of the distinct target peptides) come from a residue-shuffled proteome or from the
related strain's proteome and never occur, up to I/L, in a protein of the sample. Scores follow a
two-component normal mixture from which PEP and q-values are computed consistently; decoy rows
are reversed tryptic peptides with null scores.

All draws come from ``numpy.random.default_rng([seed, 3])``, so the genome and its truth are
unchanged by the proteomics options. Gene-level abundance is a plain log-normal per protein
(``gene_log_abundance``; M4 replaces it with the TU/region/AR(1) model).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np
from scipy import stats

from giapsim.codons import CODE_11, codons_of, revcomp, translate_cds
from giapsim.features import Feat

__all__ = [
    "Protein",
    "ProteinSimConfig",
    "ProteomicsConfig",
    "digest",
    "digest_spans",
    "expressed_proteins",
    "gene_log_abundance",
    "mixture_pep",
    "mixture_q",
    "simulate_proteomics",
]

AA20 = "ACDEFGHIKLMNPQRSTVWY"
_NA = "NA"
#: Mean and SD of the per-protein natural-log abundance (placeholder until M4).
LOG_ABUNDANCE = (1.0, 1.5)
INTENSITY_SD = 0.5
P_CHARGE_2 = 0.6
_STREAM = 3  # numpy SeedSequence key of the proteomics stream (1: features, 2: signals)


@dataclass(frozen=True)
class ProteinSimConfig:
    """Cross-strain protein copies in ``proteins.faa`` (for the miniprot/tblastn routes)."""

    identity: tuple[float, float] = (0.80, 0.99)
    fraction: float = 0.2  # share of annotated CDS proteins (>= 50 aa) that get a copy
    min_length_aa: int = 50

    def __post_init__(self) -> None:
        object.__setattr__(self, "identity", tuple(float(x) for x in self.identity))
        lo, hi = self.identity
        if not 0.0 < lo <= hi < 1.0:
            raise ValueError("identity must satisfy 0 < lo <= hi < 1")
        if not 0.0 <= self.fraction <= 1.0:
            raise ValueError("fraction must be in [0, 1]")


@dataclass(frozen=True)
class ProteomicsConfig:
    """Peptide simulation knobs (spec §11.3; docs/specs/M3-19-sim-peptides-db-errors.md)."""

    digest_len: tuple[int, int] = (7, 30)
    missed_cleavages: int = 1  # trypsin: after K/R, also before P (Trypsin/P)
    detect_logit: tuple[float, float, float] = (-2.0, 0.08, 0.9)  # intercept, per aa, log abund.
    pi_wrong: float = 0.01
    wrong_source: Literal["shuffled", "other_strain"] = "shuffled"
    score_mixture: tuple[float, float, float, float] = (3.0, 1.0, 0.0, 1.0)  # mu1, sd1, mu0, sd0
    n_decoy_rows: Literal["match_wrong"] | int = "match_wrong"
    search_db: Literal["same", "other"] = "same"
    n_runs: int = 3
    proteins: ProteinSimConfig | None = None

    def __post_init__(self) -> None:
        for name in ("digest_len", "detect_logit", "score_mixture"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        lo, hi = self.digest_len
        if not 1 <= lo <= hi:
            raise ValueError("digest_len must satisfy 1 <= lo <= hi")
        if self.missed_cleavages < 0:
            raise ValueError("missed_cleavages must be >= 0")
        if not 0.0 <= self.pi_wrong < 1.0:
            raise ValueError("pi_wrong must be in [0, 1)")
        if self.wrong_source not in ("shuffled", "other_strain"):
            raise ValueError("wrong_source must be 'shuffled' or 'other_strain'")
        if self.search_db not in ("same", "other"):
            raise ValueError("search_db must be 'same' or 'other'")
        if self.score_mixture[1] <= 0 or self.score_mixture[3] <= 0:
            raise ValueError("score_mixture SDs must be > 0")
        if self.n_decoy_rows != "match_wrong" and (
            not isinstance(self.n_decoy_rows, int) or self.n_decoy_rows < 0
        ):
            raise ValueError("n_decoy_rows must be 'match_wrong' or an int >= 0")
        if self.n_runs < 1:
            raise ValueError("n_runs must be >= 1")


# ---------------------------------------------------------------- digest


def digest_spans(
    protein: str, digest_len: tuple[int, int] = (7, 30), missed_cleavages: int = 1
) -> list[tuple[int, int]]:
    """Spans ``[a, b)`` of the Trypsin/P products of ``protein`` (cut after every K/R), in order
    of ``a`` then ``b``, with ``digest_len[0] <= b - a <= digest_len[1]``."""
    if not protein:
        return []
    cuts = [0, *(i + 1 for i, c in enumerate(protein[:-1]) if c in "KR"), len(protein)]
    lo, hi = digest_len
    out = []
    for i in range(len(cuts) - 1):
        for j in range(i + 1, min(i + 2 + missed_cleavages, len(cuts))):
            if lo <= cuts[j] - cuts[i] <= hi:
                out.append((cuts[i], cuts[j]))
    return out


def digest(
    protein: str, digest_len: tuple[int, int] = (7, 30), missed_cleavages: int = 1
) -> list[str]:
    """Distinct Trypsin/P products of ``protein`` in order of first occurrence."""
    spans = digest_spans(protein, digest_len, missed_cleavages)
    return list(dict.fromkeys(protein[a:b] for a, b in spans))


# ---------------------------------------------------------------- statistics


def mixture_pep(
    scores: np.ndarray, pi: float, mix: tuple[float, float, float, float]
) -> np.ndarray:
    """Posterior error probability pi f0 / (pi f0 + (1 - pi) f1) of each score."""
    mu1, sd1, mu0, sd0 = mix
    s = np.asarray(scores, dtype=float)
    if pi <= 0.0:
        return np.zeros_like(s)
    l0 = math.log(pi) + stats.norm.logpdf(s, mu0, sd0)
    l1 = math.log1p(-pi) + stats.norm.logpdf(s, mu1, sd1)
    out: np.ndarray = np.exp(l0 - np.logaddexp(l0, l1))
    return out


def mixture_q(scores: np.ndarray, pi: float, mix: tuple[float, float, float, float]) -> np.ndarray:
    """q-values: the mixture FDR pi S0 / (pi S0 + (1 - pi) S1) of the score threshold, made
    monotone (minimum over every threshold at or below the score)."""
    mu1, sd1, mu0, sd0 = mix
    s = np.asarray(scores, dtype=float)
    if pi <= 0.0 or s.size == 0:
        return np.zeros_like(s)
    l0 = math.log(pi) + stats.norm.logsf(s, mu0, sd0)
    l1 = math.log1p(-pi) + stats.norm.logsf(s, mu1, sd1)
    fdr: np.ndarray = np.exp(l0 - np.logaddexp(l0, l1))
    order = np.argsort(s, kind="stable")
    q: np.ndarray = np.empty_like(fdr)
    q[order] = np.minimum.accumulate(fdr[order])
    return q


# ---------------------------------------------------------------- proteins


@dataclass(frozen=True)
class Protein:
    """An expressed true protein; 0-based ``[start, end)`` unwrapped (end may exceed the contig)."""

    locus: str
    contig: str
    start: int
    end: int
    strand: str
    aa: str


def _slice(seq: str, s0: int, e0: int, strand: str) -> str:
    n = len(seq)
    nt = seq[s0:e0] if e0 <= n else seq[s0:] + seq[: e0 - n]
    return revcomp(nt) if strand == "-" else nt


def expressed_proteins(
    feats: Sequence[Feat], novel: Sequence[Sequence[Any]], seqs: dict[str, str]
) -> list[Protein]:
    """Non-pseudo CDS and sORFs translated from the true start, then the novel ORFs."""
    out = []
    for f in feats:
        if f.kind not in ("cds", "sorf") or f.pseudo_reason or f.extra.get("pseudo"):
            continue
        s0, e0 = f.span()
        aa = translate_cds(codons_of(_slice(seqs[f.contig], s0, e0, f.strand)))
        out.append(Protein(f.locus, f.contig, s0, e0, f.strand, aa))
    for orf_id, _cls, contig, s1, e0, strand, aa in novel:
        out.append(Protein(str(orf_id), str(contig), int(s1) - 1, int(e0), str(strand), str(aa)))
    return out


def gene_log_abundance(rng: np.random.Generator, n: int) -> np.ndarray:
    """Natural-log abundance of ``n`` proteins: i.i.d. normal (M4 replaces this model)."""
    out: np.ndarray = np.round(rng.normal(*LOG_ABUNDANCE, size=n), 6)
    return out


def _il(s: str) -> str:
    return s.replace("I", "L")


# ---------------------------------------------------------------- six-frame locus counts

_NT_CODE = np.full(256, 0, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _NT_CODE[_b] = _i
_CODON_AA = np.frombuffer(
    "".join(CODE_11[a + b + c] for a in "ACGT" for b in "ACGT" for c in "ACGT").encode(),
    dtype=np.uint8,
)
_K = 7  # prefix length of the six-frame index


def _translate_np(nt: bytes) -> bytes:
    m = len(nt) // 3
    c = _NT_CODE[np.frombuffer(nt[: 3 * m], dtype=np.uint8)].reshape(m, 3)
    out: bytes = _CODON_AA[16 * c[:, 0] + 4 * c[:, 1] + c[:, 2]].tobytes()
    return out


def _kmer_codes(aa: bytes, k: int) -> np.ndarray:
    a = np.frombuffer(aa, dtype=np.uint8).astype(np.int64) & 31
    m = len(a) - k + 1
    if m <= 0:
        return np.zeros(0, dtype=np.int64)
    code = np.zeros(m, dtype=np.int64)
    for j in range(k):
        code = (code << 5) | a[j : j + m]
    return code


def sixframe_counts(
    peptides: Sequence[str], seqs: dict[str, str], topology: dict[str, str]
) -> dict[str, int]:
    """Number of distinct loci (contig, plus-strand start, strand) whose six-frame translation is
    exactly each peptide; circular contigs are searched across the origin."""
    peps = sorted({p for p in peptides if len(p) >= _K})
    by_prefix: dict[int, list[str]] = {}
    for p in peps:
        by_prefix.setdefault(int(_kmer_codes(p[:_K].encode(), _K)[0]), []).append(p)
    keys = np.array(sorted(by_prefix), dtype=np.int64)
    loci: dict[str, set[tuple[str, int, str]]] = {p: set() for p in peps}
    max_len = max((len(p) for p in peps), default=0)
    for cid, seq in seqs.items():
        n = len(seq)
        circular = topology.get(cid, "circular") == "circular"
        ext = seq + seq[: 3 * max_len] if circular else seq
        for strand, src in (("+", ext), ("-", revcomp(ext))):
            raw = src.encode()
            for off in range(3):
                aa = _translate_np(raw[off:])
                hits = np.flatnonzero(np.isin(_kmer_codes(aa, _K), keys))
                text = aa.decode()
                for i in hits.tolist():
                    for p in by_prefix[(int(_kmer_codes(aa[i : i + _K], _K)[0]))]:
                        if text.startswith(p, i):
                            s = off + 3 * i if strand == "+" else len(ext) - off - 3 * (i + len(p))
                            loci[p].add((cid, s % n if circular else s, strand))
    out = {p: len(v) for p, v in loci.items()}
    for p in peptides:
        out.setdefault(p, 0)
    return out


# ---------------------------------------------------------------- simulation


@dataclass
class _Target:
    peptide: str
    wrong: bool
    origin: str  # gene | shuffled | other_strain
    loc: int  # index of the protein whose abundance sets the intensity
    ids: str  # protein_ids


def _detect(
    cfg: ProteomicsConfig, rng: np.random.Generator, prots: list[Protein], la: np.ndarray
) -> tuple[list[str], dict[str, list[tuple[int, int]]], dict[str, int]]:
    """Detected distinct peptides, every (protein, offset) occurrence and the first protein."""
    occ: dict[str, list[tuple[int, int]]] = {}
    first: dict[str, int] = {}
    for k, p in enumerate(prots):
        for a, b in digest_spans(p.aa, cfg.digest_len, cfg.missed_cleavages):
            pep = p.aa[a:b]
            if pep not in occ:
                occ[pep] = []
                first[pep] = k
            occ[pep].append((k, a))
    order = list(occ)
    b0, b1, b2 = cfg.detect_logit
    lengths = np.array([len(x) for x in order], dtype=float)
    abund = la[[first[x] for x in order]] if order else np.zeros(0)
    eta = b0 + b1 * lengths + b2 * abund
    prob = 1.0 / (1.0 + np.exp(-np.clip(eta, -700, 700)))
    hit = rng.random(len(order)) < prob
    detected = [x for x, h in zip(order, hit.tolist(), strict=True) if h]
    return detected, occ, first


def _wrong_peptides(
    cfg: ProteomicsConfig,
    rng: np.random.Generator,
    n: int,
    prots: list[Protein],
    reference: Sequence[Feat],
    hay: str,
) -> list[tuple[str, str]]:
    """``n`` distinct wrong peptides (peptide, source locus) absent from ``hay`` (I = L)."""
    out: dict[str, str] = {}
    if n <= 0:
        return []

    def ok(pep: str) -> bool:
        return pep not in out and _il(pep) not in hay

    if cfg.wrong_source == "other_strain":
        pool = sorted(
            {
                (pep, f.locus)
                for f in reference
                if f.kind in ("cds", "sorf") and f.aa
                for pep in digest(f.aa, cfg.digest_len, cfg.missed_cleavages)
            }
        )
        for i in rng.permutation(len(pool)).tolist():
            pep, loc = pool[i]
            if ok(pep):
                out[pep] = loc
                if len(out) == n:
                    break
        return list(out.items())
    usable = [p for p in prots if len(p.aa) >= cfg.digest_len[0]]
    tries = 0
    while usable and len(out) < n and tries < 100 * n + 1000:
        tries += 1
        p = usable[int(rng.integers(len(usable)))]
        shuffled = "".join(rng.permutation(np.array(list(p.aa))).tolist())
        prods = digest(shuffled, cfg.digest_len, cfg.missed_cleavages)
        if not prods:
            continue
        pep = prods[int(rng.integers(len(prods)))]
        if ok(pep):
            out[pep] = p.locus
    return list(out.items())


def _fmt(x: float) -> str:
    return f"{x:.10g}"


def _copies(
    pc: ProteinSimConfig, rng: np.random.Generator, annotated: list[tuple[str, str, str]]
) -> list[tuple[str, str]]:
    """Cross-strain copies (UniProt-like ID, sequence) at identity drawn from ``pc.identity``."""
    src = [
        (loc, aa) for loc, kind, aa in annotated if kind == "cds" and len(aa) >= pc.min_length_aa
    ]
    k = round(pc.fraction * len(src))
    if src and pc.fraction > 0:
        k = max(k, 1)
    picks = sorted(rng.choice(len(src), size=k, replace=False).tolist()) if k else []
    out, used = [], set()
    alnum = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
    for i in picks:
        aa = src[i][1]
        ident = float(rng.uniform(*pc.identity))
        n_sub = min(max(1, round((1 - ident) * len(aa))), int(0.2 * len(aa)))
        seq = list(aa)
        for j in rng.choice(np.arange(1, len(aa)), size=n_sub, replace=False).tolist():
            alt = [c for c in AA20 if c != seq[j]]
            seq[j] = alt[int(rng.integers(len(alt)))]
        while True:
            acc = (
                "OPQ"[int(rng.integers(3))]
                + str(int(rng.integers(10)))
                + "".join(alnum[int(x)] for x in rng.integers(36, size=3))
                + str(int(rng.integers(10)))
            )
            if acc not in used:
                used.add(acc)
                break
        out.append((f"tr|{acc}|{acc}_SIMUL", "".join(seq)))
    return out


def _write_fasta(path: Path, records: Sequence[tuple[str, str]]) -> None:
    parts = []
    for name, seq in records:
        parts.append(f">{name}\n" + "".join(seq[i : i + 60] + "\n" for i in range(0, len(seq), 60)))
    path.write_text("".join(parts))


def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    lines = ["\t".join(header)] + ["\t".join(r) for r in rows]
    path.write_text("\n".join(lines) + "\n")


PEPTIDE_HEADER = (
    "peptide", "protein_ids", "q_value", "pep", "score", "is_decoy", "intensity", "run", "charge",
    "spectrum_id", "source", "search_db",
)  # fmt: skip
TRUTH_HEADER = (
    "peptide", "is_wrong", "origin", "locus_tag", "contig", "start", "end", "strand",
    "n_loci_genome",
)  # fmt: skip


def simulate_proteomics(
    cfg: ProteomicsConfig,
    seed: int,
    out: Path,
    feats: Sequence[Feat],
    novel: Sequence[Sequence[Any]],
    seqs: dict[str, str],
    topology: dict[str, str],
    reference: Sequence[Feat] | None = None,
) -> dict[str, Path]:
    """Write ``peptides.tsv``, ``proteins.faa``, ``abundance.tsv`` and ``truth/peptides.tsv``."""
    if cfg.wrong_source == "other_strain" and reference is None:
        raise ValueError("wrong_source='other_strain' needs SimConfig.strain")
    rng = np.random.default_rng([seed, _STREAM])
    prots = expressed_proteins(feats, novel, seqs)
    la = gene_log_abundance(rng, len(prots))
    detected, occ, first = _detect(cfg, rng, prots, la)
    annotated = [(f.locus, f.kind, f.aa) for f in feats if f.kind in ("cds", "sorf") and f.aa]
    annotated_loci = {loc for loc, _, _ in annotated}
    hay = "|".join(_il(p.aa) for p in prots) + "|" + "|".join(_il(a) for _, _, a in annotated)
    pi = cfg.pi_wrong
    n_wrong = round(len(detected) * pi / (1 - pi))
    wrong = _wrong_peptides(cfg, rng, n_wrong, prots, reference or (), hay)
    origin = "other_strain" if cfg.wrong_source == "other_strain" else "shuffled"
    targets = [
        _Target(
            pep,
            False,
            "gene",
            first[pep],
            ";".join(dict.fromkeys(prots[k].locus for k, _ in occ[pep])),
        )
        for pep in detected
    ]
    targets += [
        _Target(pep, True, origin, int(rng.integers(len(prots))), loc) for pep, loc in wrong
    ]

    # target rows: every target peptide is identified once in every run
    n_t, n_runs = len(targets), cfg.n_runs
    mu1, sd1, mu0, sd0 = cfg.score_mixture
    is_wrong = np.array([t.wrong for t in targets], dtype=bool)
    z = rng.standard_normal((n_runs, n_t))
    scores = np.round(np.where(is_wrong, mu0 + sd0 * z, mu1 + sd1 * z), 6)
    abund = la[[t.loc for t in targets]] if n_t else np.zeros(0)
    intensity = np.exp(abund + INTENSITY_SD * rng.standard_normal((n_runs, n_t)))
    charge = np.where(rng.random((n_runs, n_t)) < P_CHARGE_2, 2, 3)
    peps = mixture_pep(scores.ravel(), pi, cfg.score_mixture).reshape(scores.shape)
    qs = mixture_q(scores.ravel(), pi, cfg.score_mixture).reshape(scores.shape)

    # decoy rows: reversed tryptic peptides of expressed proteins, null scores
    n_dec = int(is_wrong.sum()) * n_runs if cfg.n_decoy_rows == "match_wrong" else cfg.n_decoy_rows
    target_set = {t.peptide for t in targets}
    decoys: list[tuple[str, str, int, int]] = []  # peptide, protein_ids, run, charge
    usable = [k for k, p in enumerate(prots) if digest_spans(p.aa, cfg.digest_len)]
    tries = 0
    while usable and len(decoys) < n_dec and tries < 100 * n_dec + 1000:
        tries += 1
        p = prots[usable[int(rng.integers(len(usable)))]]
        spans = digest_spans(p.aa, cfg.digest_len, cfg.missed_cleavages)
        a, b = spans[int(rng.integers(len(spans)))]
        rev = p.aa[a:b][::-1]
        if rev in target_set:
            continue
        run = int(rng.integers(n_runs))
        decoys.append((rev, f"rev_{p.locus}", run, 2 if rng.random() < P_CHARGE_2 else 3))

    rows: list[list[str]] = []
    for r in range(n_runs):
        run_rows: list[list[str]] = []
        for i, t in enumerate(targets):
            run_rows.append([
                t.peptide, t.ids, _fmt(qs[r, i]), _fmt(peps[r, i]), _fmt(scores[r, i]), "false",
                f"{intensity[r, i]:.6g}", f"run{r + 1}", str(charge[r, i]),
            ])  # fmt: skip
        for pep, ids, run, ch in decoys:
            if run == r:
                run_rows.append([pep, ids, _NA, _NA, _NA, "true", _NA, f"run{r + 1}", str(ch)])
        for j, i in enumerate(rng.permutation(len(run_rows)).tolist()):
            rows.append([*run_rows[i], f"run{r + 1}.{j + 1}", "database", cfg.search_db])

    # truth
    counts = sixframe_counts(detected, seqs, topology)
    truth: list[list[str]] = []
    for pep in detected:
        n_loci = str(max(1, counts.get(pep, 0)))
        for k, a in occ[pep]:
            p = prots[k]
            n = len(seqs[p.contig])
            if p.strand == "+":
                s0 = p.start + 3 * a
                e0 = s0 + 3 * len(pep)
            else:
                e0 = p.end - 3 * a
                s0 = e0 - 3 * len(pep)
            if s0 >= n:
                s0, e0 = s0 - n, e0 - n
            truth.append(
                [pep, "false", "gene", p.locus, p.contig, str(s0 + 1), str(e0), p.strand, n_loci]
            )
    for pep, _loc in wrong:
        truth.append([pep, "true", origin, _NA, _NA, _NA, _NA, _NA, _NA])
    for pep in dict.fromkeys(d[0] for d in decoys):
        truth.append([pep, "true", "decoy", _NA, _NA, _NA, _NA, _NA, _NA])

    pseudo = {f.locus for f in feats if f.pseudo_reason or f.extra.get("pseudo")}
    faa = [(loc, aa) for loc, _, aa in annotated if loc not in pseudo]
    if cfg.proteins is not None:
        faa += [
            c
            for c in _copies(cfg.proteins, rng, [x for x in annotated if x[0] not in pseudo])
            if c[0] not in annotated_loci
        ]

    paths = {
        "peptides": out / "peptides.tsv",
        "proteins": out / "proteins.faa",
        "abundance": out / "abundance.tsv",
        "truth": out / "truth" / "peptides.tsv",
    }
    _write_tsv(paths["peptides"], PEPTIDE_HEADER, rows)
    _write_fasta(paths["proteins"], faa)
    _write_tsv(
        paths["abundance"],
        ("locus_tag", "log_abundance"),
        [[p.locus, f"{x:.6f}"] for p, x in zip(prots, la.tolist(), strict=True)],
    )
    _write_tsv(paths["truth"], TRUTH_HEADER, truth)
    return paths
