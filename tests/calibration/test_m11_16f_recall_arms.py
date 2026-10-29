"""Genetic-code arms 4 and 25 of the exact and tolerant recall suites (M11-16f; docs/specs/
M11-16-translation-tables-4-25.md item 9: "W/G at TGA must map", "fewer stops, longer frames";
decision D14; AC 15).

Arm 11 is ``test_exact_recall.py`` / ``test_tolerant_recall.py``, unchanged. Each arm here is
that suite's test with the genome simulated at genetic code ``t`` (``SimConfig.genetic_code``),
the six-frame index built at ``t`` and the oracle translating with gippasim's own table
(``gippasim.codons.code(t)``, independent of ``gippa.model.genetic_code``). Same seeds (exact:
101 / 102; tolerant: ``calibration_base`` + 700 + i), same threshold ``mapping_recall_min``, same
exact criteria of the tolerant suite. Non-vacuity: planted peptides whose origin holds an
in-frame TGA (read as W at 4, G at 25) are part of the planted set and counted separately; each
of them must be recovered (they are the peptides a table-11 index cannot map).

Fast tier: arm 4. ``calibration_full``: arms 4 and 25.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.genetic_code_arms import RECALL_FAST, RECALL_FULL, TGA_AA, in_frame_tga
from test_tolerant_recall import _brute_min, _overlaps, _plant, _tryptic

pytestmark = pytest.mark.calsuite("placement")

_RC = str.maketrans("ACGTN", "TGCAN")


def _revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def _world(length: int, seed: int, table: int, tmp: Path):  # type: ignore[no-untyped-def]
    from gippa.io.fasta import read_genome_fasta
    from gippa.mapping.sixframe import build_sixframe
    from gippa.model import Replicon
    from gippa.model.genes import gene_model_from_bakta_json
    from gippasim.simulate import SimConfig, simulate

    res = simulate(SimConfig(genome_length=length, genetic_code=table), seed=seed, outdir=tmp)
    (rec,) = read_genome_fasta(res.genome_fasta)
    rep = Replicon(rec.id, len(rec.sequence), "circular")
    model = gene_model_from_bakta_json(res.bakta_json, {rec.id: rep})
    idx = build_sixframe({rec.id: rec.sequence}, {rec.id: rep}, table=table)
    assert idx.translation_table == table
    return rec, rep, model, idx


def _origin_nt(genome: str, start: int, end: int, strand: str) -> str:
    n = len(genome)
    nt = "".join(genome[x % n] for x in range(start, end))
    return nt if strand == "+" else _revcomp(nt)


def _translate(nt: str, table: int) -> str:
    from gippasim.codons import code

    t = code(table)
    return "".join(t[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))


def _exact_arm(length: int, seed: int, table: int, tmp: Path, recall_min: float) -> None:
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.exact import match_exact
    from gippa.model import Interval
    from gippa.rng import for_step

    rec, rep, model, idx = _world(length, seed, table, tmp)
    planted: dict[str, set[Interval]] = {}
    tga: set[tuple[str, Interval]] = set()
    for g in model.genes:
        iv = g.interval
        for off, pep in _tryptic(g.aa):
            if not 7 <= len(pep) <= 40 or (off == 0 and g.start_codon != "ATG"):
                continue
            s = iv.start + 3 * off if iv.strand == "+" else iv.end - 3 * (off + len(pep))
            s %= rep.length
            want = Interval(rec.id, s, s + 3 * len(pep), iv.strand)
            nt = _origin_nt(rec.sequence, want.start, want.end, want.strand)
            assert _translate(nt, table) == pep  # oracle check at the arm's table
            planted.setdefault(pep, set()).add(want)
            if in_frame_tga(nt):
                assert TGA_AA[table] in pep
                tga.add((pep, want))
    assert len(planted) > 100
    assert len(tga) >= 20, f"only {len(tga)} planted peptides read TGA as sense: arm is vacuous"
    queries = sorted(planted)
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))}
    df = match_exact(queries, idx, decoys, dict.fromkeys(queries, "reverse"))
    found: dict[str, set[Interval]] = {}
    for r in df[df["search"] == "target"].itertuples(index=False):
        found.setdefault(r.peptide, set()).add(
            Interval(r.contig, int(r.start), int(r.end), r.strand)
        )
    n_planted = sum(len(v) for v in planted.values())
    n_found = sum(len(v & found.get(p, set())) for p, v in planted.items())
    missed_tga = sorted((p, iv) for p, iv in tga if iv not in found.get(p, set()))
    assert not missed_tga, (table, len(missed_tga), missed_tga[:10])
    recall = n_found / n_planted
    assert recall >= recall_min, (table, recall, n_found, n_planted)


@pytest.mark.calibration
@pytest.mark.parametrize("table", RECALL_FAST)
@pytest.mark.parametrize(("length", "seed"), [(50_000, 101), (100_000, 102)])
def test_planted_peptide_recall_code_arm(
    length: int, seed: int, table: int, tmp_path: Path, tolerances: dict[str, object]
) -> None:
    _exact_arm(length, seed, table, tmp_path, float(tolerances["mapping_recall_min"]))  # type: ignore[arg-type]


@pytest.mark.calibration_full
@pytest.mark.parametrize("table", RECALL_FULL)
@pytest.mark.parametrize(("length", "seed"), [(50_000, 101), (100_000, 102)])
def test_planted_peptide_recall_code_arm_full(
    length: int, seed: int, table: int, tmp_path: Path, tolerances: dict[str, object]
) -> None:
    _exact_arm(length, seed, table, tmp_path, float(tolerances["mapping_recall_min"]))  # type: ignore[arg-type]


def _origin_cost(genome: str, o, query: str, k: int, table: int) -> int:  # type: ignore[no-untyped-def]
    """``test_tolerant_recall._origin_cost`` with the origin translated at ``table``."""
    from brute_tolerant import INF
    from oracle_tolerant import align_window

    from gippa.mapping.sixframe import decode, encode_query, translate

    nt = _origin_nt(genome, o.start - 3 * k, o.end + 3 * k, o.strand)
    text = decode(translate(nt.encode("ascii"), table))
    q = decode(encode_query(query))
    m = (o.end - o.start) // 3
    best = INF
    for s in range(0, 2 * k + 1):
        for e in range(m, min(k + m + k, len(text)) + 1):
            if e - s < len(q) - k or e - s > len(q) + 2 * k:
                continue
            got = align_window(q, text, s, e, None, k)
            if got is not None:
                best = min(best, got[0])
    return best


def _tolerant_arm(length: int, i: int, table: int, tmp: Path, tol: dict[str, object]) -> None:
    from gippa.config.model import TolerantConfig
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.tolerant import match_tolerant
    from gippa.model import Interval
    from gippa.rng import for_step

    seed = int(tol["seeds"]["calibration_base"]) + 700 + i  # type: ignore[index]
    rec, rep, model, idx = _world(length, seed, table, tmp)
    rng = for_step("calibration_tolerant_recall", seed)
    planted: dict[str, tuple[int, Interval]] = {}
    clash: set[str] = set()
    types_seen: set[str] = set()
    for g in model.genes:
        iv = g.interval
        for off, pep in _tryptic(g.aa):
            if not 10 <= len(pep) <= 30 or off == 0 or "*" in pep:
                continue
            k = 1 if len(pep) <= 11 else 2
            n_ops = int(rng.integers(1, k + 1))
            query, kinds = _plant(pep, n_ops, rng)
            if not 10 <= len(query) <= 40 or (len(query) <= 11 and n_ops > 1):
                continue
            s = iv.start + 3 * off if iv.strand == "+" else iv.end - 3 * (off + len(pep))
            s %= rep.length
            origin = Interval(rec.id, s, s + 3 * len(pep), iv.strand)
            if query in planted:
                clash.add(query)
            planted[query] = (n_ops, origin)
            types_seen |= set(kinds)
    for q in clash:
        del planted[q]
    assert len(planted) > 200
    assert types_seen == {"substitution", "swap", "block", "near_block"}
    n_tga = sum(
        in_frame_tga(_origin_nt(rec.sequence, o.start, o.end, o.strand))
        for _, o in planted.values()
    )
    assert n_tga >= 20, f"only {n_tga} in-frame TGA codons under planted origins: arm is vacuous"
    queries = sorted(planted)
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))}
    df = match_tolerant(
        queries, idx, decoys, dict.fromkeys(queries, "reverse"), set(queries), TolerantConfig()
    )
    found: dict[str, set[tuple[str, int, int, str]]] = {}
    cost: dict[str, int] = {}
    for r in df[df["search"] == "target"].itertuples(index=False):
        found.setdefault(r.peptide, set()).add(
            (r.contig, int(r.start) % rep.length, int(r.end) - int(r.start), r.strand)
        )
        cost[r.peptide] = int(r.n_diff)
    missed, wrong_min, bad_plant, hit = [], [], [], 0
    for q, (n_ops, o) in planted.items():
        k_q = 1 if len(q) <= 11 else 2
        c_o = _origin_cost(rec.sequence, o, q, k_q, table)
        if c_o > n_ops:
            bad_plant.append((q, n_ops, c_o, o))
        c_min = _brute_min(q, idx.frames, k_q)
        if q not in cost or cost[q] != c_min:
            wrong_min.append((q, cost.get(q), c_min))
        if any(_overlaps(o, f, rep.length) for f in found.get(q, set())):
            hit += 1
        elif c_o == c_min:
            missed.append((q, n_ops, o, c_o, cost.get(q)))
    assert not bad_plant, bad_plant[:10]
    assert not wrong_min, wrong_min[:10]
    assert not missed, missed[:10]
    recall = hit / len(planted)
    assert recall >= float(tol["mapping_recall_min"]), (table, recall, hit, len(planted))  # type: ignore[arg-type]


@pytest.mark.calibration
@pytest.mark.calibration_nightly
@pytest.mark.parametrize("table", RECALL_FAST)
@pytest.mark.parametrize(("length", "i"), [(50_000, 0), (100_000, 1)])
def test_planted_tolerant_peptide_recall_code_arm(
    length: int, i: int, table: int, tmp_path: Path, tolerances: dict[str, object]
) -> None:
    _tolerant_arm(length, i, table, tmp_path, tolerances)


@pytest.mark.calibration_full
@pytest.mark.parametrize("table", RECALL_FULL)
@pytest.mark.parametrize(("length", "i"), [(50_000, 0), (100_000, 1)])
def test_planted_tolerant_peptide_recall_code_arm_full(
    length: int, i: int, table: int, tmp_path: Path, tolerances: dict[str, object]
) -> None:
    _tolerant_arm(length, i, table, tmp_path, tolerances)
