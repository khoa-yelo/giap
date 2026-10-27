"""Codon usage: CAI, tAI, ENC, Fop, GC3, codon-pair bias, ``at_2_10``, ``rare_codons_2_10``
(docs/specs/M2-17-codon-usage.md; GIPPA_SPEC.md §6.2, §7.1.1, §7.1.6, §9.3, §9.4).

Examples 1-10, the edge cases, the ``genetic_code=11`` spy (acceptance criterion 3), the
``codon_usage`` step (registry keys equal to the ANNOTATE call) and a K-12 excerpt of the
recorded Bakta JSON (``tests/data/codon_k12_excerpt.json.gz``) checked against the independent
oracle in ``tests/helpers/codon_oracle.py`` plus a few pinned numbers.

Ambiguities resolved here (test author):
* Q1/Q2 lead defaults: codons 2 .. n-1 for every index; empty reference → CAI/Fop ``None`` and
  run flag ``cai_reference_missing``.
* "< 10 codons" counts the full CDS (start and stop included): 30 nt is scored (8 codons left),
  27 nt is not. ``at_2_10`` / ``rare_2_10`` need ``len(nt) >= 30`` and ignore ``len % 3``.
* The CAI/Fop reference is the ribosomal-protein genes (``is_ribosomal_protein``, M2-04) whose
  ``scored_sequence`` is not ``None``; ``cai_reference`` is true exactly for those genes.
* Genome-wide shares (``rare_codons_2_10``) and the CPB reference use the scored sequences of
  every gene with ``scored_sequence`` not ``None``. In an observed family a codon never seen has
  share 0.0 (so it is rare); families never observed are absent from ``synonymous_shares``.
* ``codon.parquet``: ``locus_tag``, the 8 ``CodonFeatures`` fields in dataclass order, then
  ``cai_reference``; one row per gene of ``genes.parquet`` in that file's order; missing values
  are nulls (never NaN). ``trna_counts.tsv`` with no rows → ``tai`` null everywhere.
"""

from __future__ import annotations

import csv
import dataclasses
import importlib
import json
import math
import random
import re
import tempfile
from collections import Counter
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from helpers import codon_oracle as O
from helpers.codon_excerpt import load_excerpt, to_bakta, write_inputs

import gippasim
from gippa.config import GippaConfig
from gippa.model.coords import Interval
from gippa.model.genes import Gene

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "annotate.nf"
FIELDS = ("cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10", "rare_codons_2_10")
COLUMNS = ("locus_tag", *FIELDS, "cai_reference")
PARSE_OUT = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
REL = 1e-9


def cu() -> ModuleType:
    return importlib.import_module("gippa.annotate.codon_usage")


# ------------------------------------------------------------------ synthetic genes

_SENSE = [c for c in O.SENSE if c not in ("ATG", "TGG")] + ["ATG", "TGG"]
_POS = [0]


def make_gene(
    locus: str,
    nt: str,
    *,
    gene: str | None = None,
    product: str = "hypothetical protein",
    feature_type: str = "cds",
    partial: bool = False,
) -> Gene:
    start = _POS[0]
    _POS[0] += len(nt) + 100
    return Gene(
        locus_tag=locus,
        feature_type=feature_type,  # type: ignore[arg-type]
        interval=Interval("c1", start, start + len(nt), "+"),
        aa="M",
        nt=nt,
        start_codon=nt[:3],
        gene=gene,
        product=product,
        partial=partial,
    )


def random_cds(rnd: random.Random, n_codons: int, bias: dict[str, float] | None = None) -> str:
    pool = list(bias) if bias else O.SENSE
    weights = list(bias.values()) if bias else None
    body = rnd.choices(pool, weights=weights, k=n_codons - 2)
    return "ATG" + "".join(body) + rnd.choice(["TAA", "TGA", "TAG"])


def synthetic_set(seed: int = 3, n_rib: int = 8, n_other: int = 30) -> list[Gene]:
    """Ribosomal genes drawn from a skewed codon distribution, others uniformly."""
    rnd = random.Random(seed)
    skew = {c: (8.0 if k % 3 == 0 else 1.0) for k, c in enumerate(O.SENSE)}
    genes = [
        make_gene(
            f"R_{k:03d}",
            random_cds(rnd, rnd.randint(60, 200), skew),
            gene=f"rpl{chr(65 + k)}",
            product=f"50S ribosomal protein L{k + 1}",
        )
        for k in range(n_rib)
    ]
    genes += [
        make_gene(f"O_{k:03d}", random_cds(rnd, rnd.randint(40, 400))) for k in range(n_other)
    ]
    return genes


TRNA = {"CAT": 3, "CAG": 4, "GCC": 2, "TTC": 2, "GTC": 3, "TGC": 2, "GAA": 1, "CGT": 2}


def _ok(a: float | None, b: float) -> bool:
    return a is not None and math.isclose(a, b, rel_tol=REL, abs_tol=1e-12)


# ------------------------------------------------------------------ examples 1-6


def test_example_1_at_2_10_all_at() -> None:
    assert cu().at_2_10("ATG" + "AAA" * 9 + "TAA") == 1.0


def test_example_2_at_2_10_all_gc() -> None:
    assert cu().at_2_10("ATG" + "GCC" * 9 + "TAA") == 0.0


def test_example_3_eight_codons_gives_none() -> None:
    nt = "ATG" + "AAA" * 6 + "TAA"
    assert len(nt) == 24
    assert cu().at_2_10(nt) is None
    assert cu().rare_2_10(nt, frozenset({"AAA"})) is None


def test_at_2_10_uses_codons_2_to_10_only() -> None:
    m = cu()
    assert m.at_2_10("GGG" + "AAA" * 9 + "GGG" * 20 + "TAA") == 1.0  # codon 1 and 11+ ignored
    assert m.at_2_10("ATG" + "AAA" * 4 + "GGC" * 5 + "TAA") == pytest.approx(12 / 27)
    assert m.at_2_10("ATG" + "AAA" * 9) == 1.0  # exactly 10 codons (no stop) is enough
    assert m.at_2_10("ATG" + "AAA" * 8 + "AA") is None  # 29 nt


def test_at_2_10_ignores_ambiguous_bases() -> None:
    nt = "ATG" + "ANA" + "GCN" + "AAA" * 7 + "TAA"
    # window: A N A G C N + 21 A = A/T 23 over ACGT 25
    assert cu().at_2_10(nt) == pytest.approx(23 / 25)
    assert cu().at_2_10(nt) == pytest.approx(O.at_2_10(nt))


def test_example_4_gc3() -> None:
    assert cu().gc3("GCCGCGGCA") == pytest.approx(2 / 3)


def test_example_5_rare_codons_strict_cutoff() -> None:
    m = cu()
    leu = {"CTG": 50, "TTA": 13, "TTG": 13, "CTT": 11, "CTC": 10, "CTA": 3}
    seq = "".join(c * n for c, n in leu.items())
    shares = m.synonymous_shares([seq])
    assert shares["CTA"] == pytest.approx(0.03)
    assert shares["CTC"] == pytest.approx(0.10)
    assert sum(shares[c] for c in leu) == pytest.approx(1.0)
    rare = m.rare_codons(shares)
    assert isinstance(rare, frozenset)
    assert "CTA" in rare
    assert "CTC" not in rare  # 0.10 is not < 0.10
    assert rare == frozenset({"CTA"})
    assert m.rare_codons(shares, cutoff=0.11) == frozenset({"CTA", "CTC"})


def test_example_6_rare_2_10() -> None:
    nt = "ATG" + "CTA" + "CTG" + "CTA" + "AAA" * 6 + "CTA" + "TAA"  # codon 11 CTA not counted
    assert cu().rare_2_10(nt, frozenset({"CTA"})) == 2
    assert cu().rare_2_10(nt, frozenset()) == 0


# ------------------------------------------------------------------ shares and rare: edge cases


def test_single_codon_families_share_one_never_rare_and_stops_excluded() -> None:
    m = cu()
    shares = m.synonymous_shares(["ATGTGGTAACTGCTA", "ATGTGATAG"])
    assert shares["ATG"] == 1.0 and shares["TGG"] == 1.0
    assert not {"TAA", "TAG", "TGA"} & set(shares)
    rare = m.rare_codons(shares, cutoff=1.0)  # everything below 100 % is rare
    assert "ATG" not in rare and "TGG" not in rare
    assert shares["CTG"] == 0.5 and shares["CTA"] == 0.5
    assert shares["TTA"] == 0.0 and "TTA" in m.rare_codons(shares)  # unseen Leu codon
    assert "GCC" not in shares  # Ala never observed


def test_synonymous_shares_skips_codons_with_ambiguous_bases() -> None:
    m = cu()
    seqs = ["CTGCTNCTGNTGCTA", "CTRCTA"]
    shares = m.synonymous_shares(seqs)
    assert shares["CTG"] == pytest.approx(0.5) and shares["CTA"] == pytest.approx(0.5)
    assert all(re.fullmatch("[ACGT]{3}", c) for c in shares)
    assert shares == pytest.approx(O.shares(seqs))


def test_synonymous_shares_accepts_a_generator_and_table_11() -> None:
    m = cu()
    seqs = ["GCCGCAGCAATA", "ATAATT"]
    assert m.synonymous_shares(iter(seqs)) == pytest.approx(O.shares(seqs))
    assert m.synonymous_shares(seqs, table=11) == pytest.approx(O.shares(seqs))


# ------------------------------------------------------------------ scored_sequence


def test_scored_sequence_drops_start_and_stop() -> None:
    nt = "GTG" + "GCC" * 9 + "TGA"
    assert cu().scored_sequence(make_gene("A", nt)) == "GCC" * 9


@pytest.mark.parametrize(
    ("nt", "kw", "expect"),
    [
        ("ATG" + "AAA" * 8 + "TAA", {}, "AAA" * 8),  # 10 codons: scored
        ("ATG" + "AAA" * 7 + "TAA", {}, None),  # 9 codons
        ("ATG" + "AAA" * 9 + "TAAG", {}, None),  # length % 3 != 0
        ("ATG" + "AAA" * 20 + "TAA", {"feature_type": "pseudogene"}, None),
        ("ATG" + "AAA" * 20 + "TAA", {"partial": True}, None),
        ("ATG" + "AAA" * 9 + "TAA", {"feature_type": "sorf"}, "AAA" * 9),  # sORF >= 10 codons
    ],
)
def test_scored_sequence_eligibility(nt: str, kw: dict[str, Any], expect: str | None) -> None:
    assert cu().scored_sequence(make_gene("A", nt, **kw)) == expect


def test_scored_sequence_none_without_nt() -> None:
    g = dataclasses.replace(make_gene("A", "ATG" + "AAA" * 20 + "TAA"), nt=None)
    assert cu().scored_sequence(g) is None


# ------------------------------------------------------------------ codon_usage on genes


def test_codon_features_dataclass() -> None:
    cf = cu().CodonFeatures
    assert tuple(f.name for f in dataclasses.fields(cf)) == FIELDS
    x = cf(*([None] * 8))
    with pytest.raises(dataclasses.FrozenInstanceError):
        x.cai = 1.0  # type: ignore[misc]


def _expected(genes: Sequence[Gene], trna: dict[str, int]) -> dict[str, dict[str, Any]]:
    """Oracle features (ENC and CPB from codon-bias directly, the rest independent)."""
    from codonbias.scores import CodonPairBias, EffectiveNumberOfCodons

    from gippa.annotate.bakta_layers import is_ribosomal_protein

    sc = {g.locus_tag: O.scored(g.nt, g.feature_type, g.partial) for g in genes}
    eligible = [s for s in sc.values() if s]
    ref = [sc[g.locus_tag] for g in genes if is_ribosomal_protein(g) and sc[g.locus_tag]]
    rare = O.rare(O.shares(eligible))
    enc = EffectiveNumberOfCodons(genetic_code=11)
    cpb = CodonPairBias(ref_seq=eligible, genetic_code=11)
    cw = O.cai_weights(ref) if ref else None
    tw = O.tai_weights(trna) if trna else None
    out: dict[str, dict[str, Any]] = {}
    for g in genes:
        s, nt = sc[g.locus_tag], g.nt or ""
        row: dict[str, Any] = {
            "at_2_10": O.at_2_10(nt),
            "rare_codons_2_10": O.rare_2_10(nt, rare),
        }
        if s is None:
            row.update(dict.fromkeys(("cai", "tai", "enc", "fop", "gc3", "cpb")))
        else:
            row.update(
                cai=O.geomean(cw, s) if cw else None,
                fop=O.fop(s, ref) if ref else None,
                tai=O.geomean(tw, s) if tw else None,
                enc=float(enc.get_score(s)),
                gc3=O.gc3(s),
                cpb=float(cpb.get_score(s)),
            )
        out[g.locus_tag] = row
    return out


def _check(got: dict[str, Any], want: dict[str, dict[str, Any]]) -> None:
    assert set(got) == set(want)
    for lt, w in want.items():
        f = got[lt]
        for name in FIELDS:
            v, e = getattr(f, name), w[name]
            if e is None:
                assert v is None, (lt, name, v)
            elif name == "rare_codons_2_10":
                assert v == e and isinstance(v, int), (lt, name, v, e)
            else:
                assert _ok(v, e), (lt, name, v, e)


def test_codon_usage_matches_the_oracle_on_synthetic_genes() -> None:
    genes = synthetic_set()
    genes.append(make_gene("P_001", "ATG" + "AAA" * 30 + "TAA", feature_type="pseudogene"))
    genes.append(make_gene("S_001", "ATG" + "GCC" * 9 + "TAA", feature_type="sorf"))
    genes.append(make_gene("S_002", "ATG" + "GCC" * 6 + "TAA", feature_type="sorf"))
    feats, flags = cu().codon_usage(genes, TRNA)
    assert flags == ()
    _check(feats, _expected(genes, TRNA))
    # pseudogene: no score, but at_2_10 / rare_codons_2_10 (>= 10 codons)
    p = feats["P_001"]
    assert (p.cai, p.tai, p.enc, p.fop, p.gc3, p.cpb) == (None,) * 6
    assert p.at_2_10 == 1.0 and p.rare_codons_2_10 is not None
    assert feats["S_001"].cai is not None  # sORF with 11 codons scored
    assert all(getattr(feats["S_002"], n) is None for n in FIELDS)  # 8 codons


def test_scores_lie_in_their_ranges_and_are_never_nan() -> None:
    feats, _ = cu().codon_usage(synthetic_set(seed=11), TRNA)
    for f in feats.values():
        for n in FIELDS:
            v = getattr(f, n)
            assert v is None or (isinstance(v, (int, float)) and math.isfinite(v)), (n, v)
        assert 0 < f.cai <= 1 and 0 <= f.fop <= 1 and 0 < f.tai <= 1
        assert 20 <= f.enc <= 61 and 0 <= f.gc3 <= 1 and 0 <= f.at_2_10 <= 1
        assert 0 <= f.rare_codons_2_10 <= 9


def test_reference_genes_score_above_uniform_genes() -> None:
    """Simulated: genes drawn from the reference codon bias outrank uniform ones on CAI and Fop."""
    feats, _ = cu().codon_usage(synthetic_set(seed=5, n_rib=10, n_other=40), TRNA)
    rib = [f for k, f in feats.items() if k.startswith("R_")]
    oth = [f for k, f in feats.items() if k.startswith("O_")]
    assert min(f.cai for f in rib) > max(f.cai for f in oth)
    assert sorted(f.fop for f in rib)[5] > sorted(f.fop for f in oth)[20]


def test_example_7_no_ribosomal_reference() -> None:
    rnd = random.Random(7)
    genes = [make_gene(f"O_{k}", random_cds(rnd, 120)) for k in range(12)]
    feats, flags = cu().codon_usage(genes, TRNA)
    assert flags == ("cai_reference_missing",)
    for f in feats.values():
        assert f.cai is None and f.fop is None
        assert f.tai is not None and f.enc is not None and f.gc3 is not None
        assert f.cpb is not None


def test_example_7_partial_ribosomal_genes_do_not_form_a_reference() -> None:
    rnd = random.Random(8)
    genes = [make_gene(f"O_{k}", random_cds(rnd, 120)) for k in range(6)]
    genes.append(
        make_gene("R_1", random_cds(rnd, 120), gene="rplB", product="50S ribosomal protein L2",
                  partial=True)
    )  # fmt: skip
    genes.append(make_gene("R_2", random_cds(rnd, 120), gene="rpsA", feature_type="pseudogene"))
    feats, flags = cu().codon_usage(genes, TRNA)
    assert flags == ("cai_reference_missing",)
    assert all(f.cai is None and f.fop is None for f in feats.values())


def test_example_8_no_trna_counts() -> None:
    genes = synthetic_set(seed=9)
    feats, flags = cu().codon_usage(genes, {})
    assert "cai_reference_missing" not in flags
    assert all(f.tai is None for f in feats.values())
    assert all(f.cai is not None for f in feats.values())


def test_alternative_starts_and_stop_codons_do_not_change_any_feature() -> None:
    genes = synthetic_set(seed=13)
    swapped = []
    for k, g in enumerate(genes):
        assert g.nt is not None
        nt = ("GTG", "TTG", "ATG")[k % 3] + g.nt[3:-3] + ("TGA", "TAG", "TAA")[k % 3]
        swapped.append(dataclasses.replace(g, nt=nt, start_codon=nt[:3]))
    a, _ = cu().codon_usage(genes, TRNA)
    b, _ = cu().codon_usage(swapped, TRNA)
    assert a == b


def test_gene_order_does_not_matter() -> None:
    genes = synthetic_set(seed=21)
    a, fa = cu().codon_usage(genes, TRNA)
    b, fb = cu().codon_usage(list(reversed(genes)), TRNA)
    assert a == b and fa == fb


# ------------------------------------------------------------------ genetic_code=11 spy


@pytest.fixture
def spy(monkeypatch: pytest.MonkeyPatch) -> dict[str, list[tuple[tuple[Any, ...], dict[str, Any]]]]:
    import codonbias
    import codonbias.scores as scores

    m = cu()
    calls: dict[str, list[tuple[tuple[Any, ...], dict[str, Any]]]] = {}
    names = (
        "CodonAdaptationIndex",
        "TrnaAdaptationIndex",
        "EffectiveNumberOfCodons",
        "FrequencyOfOptimalCodons",
        "CodonPairBias",
    )
    for name in names:
        orig = getattr(scores, name)
        rec = calls.setdefault(name, [])

        def init(self: Any, *a: Any, _o: Any = orig, _r: Any = rec, **k: Any) -> None:
            snap = {
                key: (val.copy() if hasattr(val, "copy") and not isinstance(val, str) else val)
                for key, val in k.items()
            }
            _r.append((a, snap))
            _o.__init__(self, *a, **k)

        fake = type(name, (orig,), {"__init__": init})
        monkeypatch.setattr(scores, name, fake)
        if hasattr(codonbias, name):
            monkeypatch.setattr(codonbias, name, fake)
        if hasattr(m, name):
            monkeypatch.setattr(m, name, fake)
    return calls


def test_genetic_code_11_reaches_every_codon_bias_constructor(spy: dict[str, Any]) -> None:
    from gippa.annotate.bakta_layers import is_ribosomal_protein

    genes = synthetic_set(seed=4)
    genes.append(make_gene("P_1", "ATG" + "GCC" * 30 + "TAA", feature_type="pseudogene"))
    cu().codon_usage(genes, TRNA)
    for name, calls in spy.items():
        assert calls, f"{name} was never constructed"
        for args, kw in calls:
            assert kw.get("genetic_code") == 11, f"{name}: genetic_code={kw.get('genetic_code')}"
            assert not args, f"{name}: pass arguments by keyword ({args!r:.80})"

    sc = {g.locus_tag: O.scored(g.nt, g.feature_type, g.partial) for g in genes}
    ribo = sorted(sc[g.locus_tag] for g in genes if is_ribosomal_protein(g))
    eligible = sorted(s for s in sc.values() if s)
    assert ribo and len(eligible) == len(genes) - 1
    for name in ("CodonAdaptationIndex", "FrequencyOfOptimalCodons"):
        assert all(sorted(kw["ref_seq"]) == ribo for _, kw in spy[name]), name
    assert all(sorted(kw["ref_seq"]) == eligible for _, kw in spy["CodonPairBias"])
    for _, kw in spy["TrnaAdaptationIndex"]:
        assert kw.get("prokaryote") is True
        df = kw["tGCN"]
        assert {"anti_codon", "GCN"} <= set(df.columns)
        got = {str(a).upper(): int(n) for a, n in zip(df["anti_codon"], df["GCN"], strict=True)}
        assert got == TRNA


def test_no_codon_bias_model_without_reference_or_trnas(spy: dict[str, Any]) -> None:
    rnd = random.Random(2)
    cu().codon_usage([make_gene(f"O_{k}", random_cds(rnd, 90)) for k in range(5)], {})
    assert not spy["CodonAdaptationIndex"] and not spy["FrequencyOfOptimalCodons"]
    assert not spy["TrnaAdaptationIndex"]


# ------------------------------------------------------------------ K-12 excerpt (recorded data)


@pytest.fixture(scope="module")
def k12() -> tuple[list[Gene], dict[str, int]]:
    from gippa.annotate.bakta import parse_bakta_json
    from gippa.annotate.bakta_layers import trna_counts

    data, fasta = to_bakta(load_excerpt())
    head, seq = fasta.split("\n")[:2]
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "k12x.json"
        p.write_text(json.dumps(data))
        ann = parse_bakta_json(p, genome={head[1:]: seq})
    return list(ann.model.genes), trna_counts(ann.model.features)


# Independently computed (tests/helpers/codon_oracle.py) on the recorded K-12 MG1655 full-DB
# Bakta JSON; the excerpt keeps the whole ribosomal set and every tRNA, so these equal the
# whole-genome values. (CAI, Fop, tAI, GC3 over codons 2..n-1, at_2_10.)
PINNED = {
    "KOIAIB_03524": ("rplB", 0.803751, 0.735294, 0.298215, 0.463235, 0.629630),
    "KOIAIB_00982": ("rpsA", 0.840637, 0.766187, 0.318162, 0.523381, 0.666667),
    "KOIAIB_03546": ("tufA_E379K", 0.866198, 0.804071, 0.310152, 0.554707, 0.777778),
}


def test_k12_excerpt_matches_the_oracle(k12: tuple[list[Gene], dict[str, int]]) -> None:
    genes, trna = k12
    assert len(genes) == 173 and sum(trna.values()) == 87
    feats, flags = cu().codon_usage(genes, trna)
    assert flags == ()
    _check(feats, _expected(genes, trna))


def test_k12_pinned_values(k12: tuple[list[Gene], dict[str, int]]) -> None:
    genes, trna = k12
    feats, _ = cu().codon_usage(genes, trna)
    by = {g.locus_tag: g for g in genes}
    for lt, (name, cai, fop, tai, gc3, at) in PINNED.items():
        assert by[lt].gene == name
        f = feats[lt]
        got = (f.cai, f.fop, f.tai, f.gc3, f.at_2_10)
        assert got == pytest.approx((cai, fop, tai, gc3, at), abs=1e-6), (name, got)


def test_k12_alternative_starts_pseudogenes_and_sorfs(
    k12: tuple[list[Gene], dict[str, int]],
) -> None:
    genes, trna = k12
    feats, _ = cu().codon_usage(genes, trna)
    alt = [g for g in genes if g.start_codon in ("GTG", "TTG")]
    assert {g.gene for g in alt} >= {"rpsT", "nhaA", "carA"}
    for g in alt:
        assert g.nt is not None
        assert cu().scored_sequence(g) == g.nt[3:-3]
        assert feats[g.locus_tag].cai is not None
    pseudo = [g for g in genes if g.feature_type == "pseudogene"]
    assert len(pseudo) == 5
    for g in pseudo:
        f = feats[g.locus_tag]
        assert (f.cai, f.tai, f.enc, f.fop, f.gc3, f.cpb) == (None,) * 6
        assert f.at_2_10 is not None and f.rare_codons_2_10 is not None
    sorfs = [g for g in genes if g.feature_type == "sorf"]
    short = [g for g in sorfs if len(g.nt or "") < 30]
    ten = [g for g in sorfs if len(g.nt or "") == 30]
    assert short and ten
    assert all(feats[g.locus_tag].at_2_10 is None for g in short)
    assert all(feats[g.locus_tag].cai is not None for g in ten)


# ------------------------------------------------------------------ the codon_usage step


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


def test_step_declaration(steps: ModuleType) -> None:
    s = steps.get("codon_usage")
    assert (s.kind, s.module, s.precomputed) == ("py", "core", None)
    assert set(s.inputs) == {"genes", "trna_counts"}
    assert s.optional_inputs == ()
    assert s.outputs == ("codon",)
    assert steps.topo_order(["codon_usage", "parse_bakta"]) == ["parse_bakta", "codon_usage"]


def test_annotate_nf_calls_codon_usage_with_the_registry_keys() -> None:
    text = re.sub(r"(?m)//.*$", "", SUBWF.read_text(encoding="utf-8"))
    m = re.search(r"def\s+(\w+)\s*=\s*stepIn\(\s*'codon_usage'\s*,\s*\[(.*?)\]\s*\)", text, re.S)
    assert m, "annotate.nf: no stepIn('codon_usage', [...])"
    keys = set(re.findall(r"kv\(\s*'(\w+)'", m.group(2)))
    var_key = dict(re.findall(r"def\s+(\w+)\s*=\s*kv\(\s*'(\w+)'", text))
    keys |= {var_key[v] for v in re.findall(r"\b(\w+)\b", m.group(2)) if v in var_key}
    assert keys == {"genes", "trna_counts"}
    call = re.search(
        rf"CODON_USAGE\(\s*{m.group(1)}\.files\s*,\s*{m.group(1)}\.keys\s*,\s*\[(.*?)\]", text
    )
    assert call and re.findall(r"'(\w+)=([\w.]+)'", call.group(1)) == [("codon", "codon.parquet")]


def _parse(reg: ModuleType, js: Path, fa: Path, work: Path) -> dict[str, Path]:
    outs = {k: work / v for k, v in PARSE_OUT.items()}
    reg.run_step("parse_bakta", GippaConfig(seed=1), {"bakta_json": js, "genome": fa}, outs)
    return outs


def _codon(reg: ModuleType, genes: Path, trna: Path, out: Path, seed: int = 1) -> Path:
    reg.run_step(
        "codon_usage", GippaConfig(seed=seed), {"genes": genes, "trna_counts": trna},
        {"codon": out},
    )  # fmt: skip
    return out


@pytest.fixture
def excerpt_work(steps: ModuleType, tmp_path: Path) -> dict[str, Path]:
    js, fa = write_inputs(tmp_path / "in")
    return _parse(steps, js, fa, tmp_path / "work")


def test_step_writes_codon_parquet(
    steps: ModuleType, excerpt_work: dict[str, Path], tmp_path: Path,
    k12: tuple[list[Gene], dict[str, int]],
) -> None:  # fmt: skip
    from gippa.annotate.bakta_layers import is_ribosomal_protein

    path = _codon(
        steps, excerpt_work["genes"], excerpt_work["trna_counts"], tmp_path / "codon.parquet"
    )
    table = pq.read_table(path)
    assert tuple(table.column_names) == COLUMNS
    types = dict(zip(table.column_names, table.schema.types, strict=True))
    assert types["locus_tag"] == pa.string()
    assert all(types[c] == pa.float64() for c in FIELDS if c != "rare_codons_2_10")
    assert types["rare_codons_2_10"] == pa.int64() and types["cai_reference"] == pa.bool_()
    rows = table.to_pylist()
    order = pq.read_table(excerpt_work["genes"], columns=["locus_tag"]).column(0).to_pylist()
    assert [r["locus_tag"] for r in rows] == order

    genes, trna = k12
    want = _expected(genes, trna)
    by = {g.locus_tag: g for g in genes}
    for r in rows:
        w = want[r["locus_tag"]]
        for n in FIELDS:
            v = r[n]
            assert (v is None) == (w[n] is None), (r["locus_tag"], n, v)
            if v is not None:
                assert not (isinstance(v, float) and math.isnan(v))
                assert math.isclose(v, w[n], rel_tol=REL, abs_tol=1e-12), (r["locus_tag"], n)
        g = by[r["locus_tag"]]
        assert r["cai_reference"] is (
            is_ribosomal_protein(g) and O.scored(g.nt, g.feature_type, g.partial) is not None
        )
    assert sum(r["cai_reference"] for r in rows) == 59


def test_example_9_rplB_is_flagged_cai_reference(
    steps: ModuleType, excerpt_work: dict[str, Path], tmp_path: Path
) -> None:
    path = _codon(steps, excerpt_work["genes"], excerpt_work["trna_counts"], tmp_path / "c.parquet")
    rows = {r["locus_tag"]: r for r in pq.read_table(path).to_pylist()}
    assert rows["KOIAIB_03524"]["cai_reference"] is True  # rplB
    assert rows["KOIAIB_03546"]["cai_reference"] is False  # tufA: not ribosomal
    assert rows["KOIAIB_00374"]["cai_reference"] is False  # lacZ


def test_example_10_reruns_are_byte_identical(
    steps: ModuleType, excerpt_work: dict[str, Path], tmp_path: Path
) -> None:
    a = _codon(steps, excerpt_work["genes"], excerpt_work["trna_counts"], tmp_path / "a.parquet")
    b = _codon(
        steps, excerpt_work["genes"], excerpt_work["trna_counts"], tmp_path / "b.parquet", seed=99
    )
    assert a.read_bytes() == b.read_bytes()


def test_step_without_trnas_writes_null_tai(
    steps: ModuleType, excerpt_work: dict[str, Path], tmp_path: Path
) -> None:
    empty = tmp_path / "trna_counts.tsv"
    with excerpt_work["trna_counts"].open() as fh:
        header = fh.readline()
    empty.write_text(header)
    rows = pq.read_table(
        _codon(steps, excerpt_work["genes"], empty, tmp_path / "c.parquet")
    ).to_pylist()
    assert all(r["tai"] is None for r in rows)
    assert sum(r["cai"] is not None for r in rows) > 100


def test_step_on_sim100k_has_no_cai_reference(steps: ModuleType, tmp_path: Path) -> None:
    """gippasim's SIM100K has no ribosomal-protein genes (example 7 at step level)."""
    res = gippasim.simulate(gippasim.SIM100K, gippasim.SIM100K_SEED, tmp_path / "sim")
    outs = _parse(steps, res.bakta_json, res.genome_fasta, tmp_path / "work")
    with outs["trna_counts"].open(newline="") as fh:
        trna = {r["anticodon"]: int(r["n_copies"]) for r in csv.DictReader(fh, delimiter="\t")}
    rows = pq.read_table(
        _codon(steps, outs["genes"], outs["trna_counts"], tmp_path / "c.parquet")
    ).to_pylist()
    genes = pq.read_table(outs["genes"]).to_pylist()
    assert len(rows) == len(genes)
    assert all(r["cai"] is None and r["fop"] is None and r["cai_reference"] is False for r in rows)
    assert (any(r["tai"] is not None for r in rows)) is bool(trna)
    assert any(r["enc"] is not None for r in rows)
    kinds = Counter(g["feature_type"] for g in genes)
    assert kinds["cds"] > 0
