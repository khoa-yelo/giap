"""Golden codon usage on the recorded full-DB Bakta JSON (M2-17 acceptance criterion 2).

On E. coli K-12 MG1655: CAI of ``rplB``, ``rpsA`` and ``tufA`` lies in the top decile and the
median CAI of the ribosomal-protein genes exceeds the median of all genes. On sim100k, K-12 and
Nissle 1917 every gene's features equal the independent oracle (tests/helpers/codon_oracle.py;
ENC and CPB from codon-bias called directly). ``pending_sample`` until M2-10 commits the trimmed
recordings under ``tests/data/recorded/bakta/1.12.2-full/``; until then the run-time recordings
in ``$GIAP_ROOT/recordings/bakta_full/<genome>/bakta/<genome>.json`` are used (skip if absent).
In K-12's Bakta 1.12.2 annotation ``tufA`` carries the AMRFinder name ``tufA_E379K``.
"""

from __future__ import annotations

import importlib
import math
import os
import statistics
from pathlib import Path

import pytest
from helpers import codon_oracle as O

from giap.annotate.bakta_layers import is_ribosomal_protein, trna_counts

REPO = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("GIAP_ROOT", REPO / "local"))
COMMITTED = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full"
GENOMES = ("sim100k", "k12_mg1655", "nissle1917")

pytestmark = pytest.mark.pending_sample


def _json(name: str) -> Path:
    for p in (
        COMMITTED / f"{name}.json.gz",
        COMMITTED / f"{name}.json",
        ROOT / "recordings" / "bakta_full" / name / "bakta" / f"{name}.json",
    ):
        if p.is_file():
            return p
    pytest.skip(f"pending_sample: no recorded full-DB Bakta JSON for {name}")


def _run(name: str) -> tuple[list, dict[str, int], dict, tuple[str, ...]]:  # type: ignore[type-arg]
    from giap.annotate.bakta import parse_bakta_json

    ann = parse_bakta_json(_json(name))
    genes = list(ann.model.genes)
    trna = trna_counts(ann.model.features)
    cu = importlib.import_module("giap.annotate.codon_usage")
    feats, flags = cu.codon_usage(genes, trna)
    return genes, trna, feats, flags


@pytest.fixture(scope="module")
def k12() -> tuple[list, dict[str, int], dict, tuple[str, ...]]:  # type: ignore[type-arg]
    return _run("k12_mg1655")


def test_k12_reference_genes_in_the_top_cai_decile(k12: tuple) -> None:  # type: ignore[type-arg]
    genes, _, feats, flags = k12
    assert flags == ()
    cai = {lt: f.cai for lt, f in feats.items() if f.cai is not None}
    assert len(cai) > 4000
    q90 = sorted(cai.values())[int(0.9 * len(cai))]
    picks = {
        "rplB": [g for g in genes if g.gene == "rplB"],
        "rpsA": [g for g in genes if g.gene == "rpsA"],
        "tufA": [g for g in genes if (g.gene or "").startswith("tufA")],
    }
    for name, gs in picks.items():
        assert len(gs) == 1, (name, [g.locus_tag for g in gs])
        assert cai[gs[0].locus_tag] >= q90, (name, cai[gs[0].locus_tag], q90)


def test_k12_ribosomal_median_cai_above_genome_median(k12: tuple) -> None:  # type: ignore[type-arg]
    genes, _, feats, _ = k12
    rib = [feats[g.locus_tag].cai for g in genes if is_ribosomal_protein(g)]
    allc = [f.cai for f in feats.values() if f.cai is not None]
    assert len(rib) == 59 and None not in rib
    assert statistics.median(rib) > statistics.median(allc)


def test_k12_rare_codons(k12: tuple) -> None:  # type: ignore[type-arg]
    """K-12 genome-wide rare codons (< 10 % of their family; codons 2..n-1 of eligible genes)."""
    genes, _, _, _ = k12
    cu = importlib.import_module("giap.annotate.codon_usage")
    seqs = [s for g in genes if (s := cu.scored_sequence(g)) is not None]
    rare = cu.rare_codons(cu.synonymous_shares(seqs))
    assert rare == frozenset({"AGA", "AGG", "ATA", "CGA", "CGG", "CTA"})


@pytest.mark.parametrize("name", GENOMES)
def test_recorded_genome_matches_the_oracle(name: str) -> None:
    from codonbias.scores import CodonPairBias, EffectiveNumberOfCodons

    genes, trna, feats, flags = _run(name)
    sc = {g.locus_tag: O.scored(g.nt, g.feature_type, g.partial) for g in genes}
    eligible = [s for s in sc.values() if s]
    ref = [sc[g.locus_tag] for g in genes if is_ribosomal_protein(g) and sc[g.locus_tag]]
    assert flags == (() if ref else ("cai_reference_missing",))
    rare = O.rare(O.shares(eligible))
    cw = O.cai_weights(ref) if ref else None
    tw = O.tai_weights(trna) if trna else None
    enc = EffectiveNumberOfCodons(genetic_code=11)
    cpb = CodonPairBias(ref_seq=eligible, genetic_code=11)
    assert set(feats) == {g.locus_tag for g in genes}
    for g in genes:
        f, s, nt = feats[g.locus_tag], sc[g.locus_tag], g.nt or ""
        assert f.rare_codons_2_10 == O.rare_2_10(nt, rare), g.locus_tag
        a = O.at_2_10(nt)
        assert (f.at_2_10 is None) == (a is None) and (a is None or math.isclose(f.at_2_10, a))
        if s is None:
            assert (f.cai, f.tai, f.enc, f.fop, f.gc3, f.cpb) == (None,) * 6, g.locus_tag
            continue
        want = {
            "cai": O.geomean(cw, s) if cw else None,
            "fop": O.fop(s, ref) if ref else None,
            "tai": O.geomean(tw, s) if tw else None,
            "gc3": O.gc3(s),
            "enc": float(enc.get_score(s)),
            "cpb": float(cpb.get_score(s)),
        }
        for k, w in want.items():
            v = getattr(f, k)
            if w is None:
                assert v is None, (g.locus_tag, k)
            else:
                assert v is not None and math.isclose(v, w, rel_tol=1e-9, abs_tol=1e-12), (
                    g.locus_tag, k, v, w,
                )  # fmt: skip
