"""Review pins for ``parse_bakta`` (M2-04 "Lead decisions after code review", item 1).

docs/specs/M2-04-parse-bakta-step.md: genes/features ordered by FASTA replicon order then start;
"Putative 30S ribosomal protein" is not ribosomal (product match anchored at the start); a tRNA
without anticodon writes no anticodon row; ``nearest_genes`` ties keep input order; Parquet columns
carry no statistics; the 16S tail is upper-cased. Each test targets a mutant that survived review.
Item 2 (the writer moves to ``gippa.io``) names no public function, so it is not pinned here.
"""

from __future__ import annotations

import importlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow.parquet as pq
import pytest

from gippa.config import GippaConfig
from gippa.model import Interval
from gippa.model.coords import Replicon
from gippa.model.genes import Gene, NcFeature

OUT_FILES = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
PARQUETS = ("genes", "features", "replicons", "xrefs")

_rng = np.random.default_rng(4207)
# FASTA order "zchr" then "apl": neither alphabetical nor the JSON order
ZCHR = "".join(_rng.choice(list("ACGT"), 3_000))
APL = "".join(_rng.choice(list("ACGT"), 1_000))


@pytest.fixture
def reg() -> Any:
    from gippa.steps import registry

    mod = importlib.import_module("gippa.steps.annotation")
    if "parse_bakta" not in registry._REGISTRY:
        importlib.reload(mod)
    return registry


def _cds(locus: str, seq: str, start: int, stop: int, **kw: Any) -> dict[str, Any]:
    return {
        "type": "cds",
        "sequence": seq,
        "start": start,
        "stop": stop,
        "strand": "+",
        "locus": locus,
        "aa": "MKV",
        "product": "hypothetical protein",
        **kw,
    }


def _trna(fid: str, seq: str, start: int, anticodon: str | None) -> dict[str, Any]:
    f: dict[str, Any] = {
        "type": "tRNA",
        "sequence": seq,
        "start": start,
        "stop": start + 75,
        "strand": "+",
        "id": fid,
        "product": "tRNA-Leu" if anticodon else "tRNA-Xxx",
        "amino_acid": "Leu" if anticodon else None,
    }
    if anticodon:
        f["anti_codon"] = anticodon
    return f


def _data(features: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}},
        "sequences": [
            {"id": "apl", "type": "plasmid", "topology": "circular", "length": len(APL)},
            {"id": "zchr", "type": "chromosome", "topology": "circular", "length": len(ZCHR)},
        ],
        "features": features,
    }


def _run(reg: Any, tmp: Path, data: dict[str, Any], genome: str | None = None) -> dict[str, Path]:
    d = tmp / "in"
    d.mkdir(parents=True, exist_ok=True)
    js, fa = d / "bakta.json", d / "genome.norm.fa"
    js.write_text(json.dumps(data))
    fa.write_text(genome or f">zchr\n{ZCHR}\n>apl\n{APL}\n")
    outs = {k: tmp / "work" / v for k, v in OUT_FILES.items()}
    reg.run_step("parse_bakta", GippaConfig(seed=1), {"bakta_json": js, "genome": fa}, outs)
    return outs


def _gene(locus: str, start: int, end: int, *, product: str = "", gene: str | None = None) -> Gene:
    return Gene(
        locus_tag=locus,
        feature_type="cds",
        interval=Interval("chr1", start, end, "+"),
        aa="M",
        nt=None,
        start_codon=None,
        gene=gene,
        product=product,
    )


# ------------------------------------------------ order: FASTA replicon, then start


def test_genes_and_features_follow_fasta_then_start(reg: Any, tmp_path: Path) -> None:
    """Locus-tag, ID, JSON and alphabetical contig orders all disagree with the required order."""
    feats = [
        _cds("A_001", "apl", 101, 400),
        _cds("A_002", "zchr", 2001, 2300),
        _cds("A_003", "zchr", 101, 400),
        _cds("A_004", "apl", 11, 90),
        _trna("f_a", "zchr", 1501, "CAG"),
        _trna("f_m", "apl", 501, "GCC"),
        _trna("f_z", "zchr", 601, "CAG"),
    ]
    outs = _run(reg, tmp_path, _data(feats))
    genes = pq.read_table(outs["genes"]).to_pylist()
    assert [g["locus_tag"] for g in genes] == ["A_003", "A_002", "A_004", "A_001"]
    fts = pq.read_table(outs["features"]).to_pylist()
    assert [f["feature_id"] for f in fts] == ["f_z", "f_a", "f_m"]
    reps = pq.read_table(outs["replicons"]).to_pylist()
    assert [r["contig"] for r in reps] == ["zchr", "apl"]


# ------------------------------------------------------------------ ribosomal match anchored


@pytest.mark.parametrize(
    ("product", "bare"),
    [
        ("Putative 30S ribosomal protein S1", "30S ribosomal protein S1"),
        ("putative 50S ribosomal protein L2", "50S ribosomal protein L2"),
        ("Probable 30S ribosomal protein", "30S ribosomal protein"),
    ],
)
def test_putative_ribosomal_product_is_not_ribosomal(product: str, bare: str) -> None:
    from gippa.annotate.bakta_layers import is_ribosomal_protein

    assert is_ribosomal_protein(_gene("g1", 0, 300, product=product)) is False
    assert is_ribosomal_protein(_gene("g2", 0, 300, product=bare)) is True


def test_putative_ribosomal_not_flagged_in_genes_parquet(reg: Any, tmp_path: Path) -> None:
    feats = [
        _cds("R_001", "zchr", 101, 400, product="Putative 30S ribosomal protein S1"),
        _cds("R_002", "zchr", 501, 800, product="30S ribosomal protein S1"),
    ]
    rows = {
        r["locus_tag"]: r
        for r in pq.read_table(_run(reg, tmp_path, _data(feats))["genes"]).to_pylist()
    }
    assert rows["R_001"]["ribosomal_protein"] is False
    assert rows["R_002"]["ribosomal_protein"] is True


# ------------------------------------------------------------------ tRNA without anticodon


def test_trna_without_anticodon_writes_no_count_row(reg: Any, tmp_path: Path) -> None:
    from gippa.annotate.bakta_layers import trna_counts
    from gippa.io.tsv import read_tsv

    feats = [
        _trna("t1", "zchr", 101, "CAG"),
        _trna("t2", "zchr", 301, None),
        _trna("t3", "apl", 101, "GCC"),
    ]
    outs = _run(reg, tmp_path, _data(feats))
    header, rows = read_tsv(outs["trna_counts"])
    assert header == ["anticodon", "amino_acid", "n_copies"]
    assert rows == [["CAG", "Leu", "1"], ["GCC", "Leu", "1"]]
    fts = {f["feature_id"]: f for f in pq.read_table(outs["features"]).to_pylist()}
    assert fts["t2"]["anticodon"] is None  # the feature is kept, with no anticodon value

    nc = [
        NcFeature(id="x", type="trna", interval=Interval("chr1", 0, 76, "+"), anticodon=None),
        NcFeature(id="y", type="trna", interval=Interval("chr1", 100, 176, "+"), anticodon="CAG"),
    ]
    assert trna_counts(nc) == {"CAG": 1}


# ------------------------------------------------------------------ nearest_genes ties


def test_nearest_genes_ties_keep_input_order() -> None:
    from gippa.annotate.bakta_layers import nearest_genes

    rep = Replicon("chr1", 10_000, "circular", "chromosome")
    feat = NcFeature(id="ori", type="oric", interval=Interval("chr1", 5_000, 5_100, "+"))
    # all three at gap 100; locus tags in reverse lexical order, positions out of order
    genes = [
        _gene("z_up", 4_600, 4_900),
        _gene("y_down", 5_200, 5_500),
        _gene("x_up2", 4_700, 4_900),
        _gene("far", 0, 100),
    ]
    assert nearest_genes(feat, genes, rep, n=3) == [("z_up", 100), ("y_down", 100), ("x_up2", 100)]
    assert nearest_genes(feat, list(reversed(genes)), rep, n=2) == [("x_up2", 100), ("y_down", 100)]


# ------------------------------------------------------------------ parquet without statistics


def test_parquet_columns_carry_no_statistics(reg: Any, tmp_path: Path) -> None:
    feats = [
        _cds("S_001", "zchr", 101, 400, db_xrefs=["COG:COG0090", "GO:0003735"]),
        _trna("s_t1", "zchr", 1001, "CAG"),
    ]
    outs = _run(reg, tmp_path, _data(feats))
    for key in PARQUETS:
        md = pq.ParquetFile(outs[key]).metadata
        assert md.num_row_groups >= 1, key
        for i in range(md.num_row_groups):
            rg = md.row_group(i)
            for j in range(rg.num_columns):
                st = rg.column(j).statistics
                assert st is None or not st.has_min_max, (key, rg.column(j).path_in_schema)


# ------------------------------------------------------------------ 16S tail upper-cased


def test_anti_sd_upper_cases_lower_case_tail() -> None:
    from gippa.annotate.bakta_layers import anti_sd

    assert anti_sd(["acgt" * 20 + "gatcacctcctta"]) == ("ACCTCCTTA", ())
    assert anti_sd(
        [
            "ACGT" * 20 + "GATCACCTCCTTT",
            "acgt" * 20 + "gatcacctcctta",
            "ACGT" * 20 + "GATCACCTCCTTA",
        ]
    ) == ("ACCTCCTTA", ())


def test_lower_case_16s_nt_gives_upper_case_anti_sd(reg: Any, tmp_path: Path) -> None:
    s16 = ("acgt" * 150 + "gatcacctcctta").lower()
    feats = [
        {
            "type": "rRNA",
            "sequence": "zchr",
            "start": 1001,
            "stop": 1000 + len(s16),
            "strand": "+",
            "id": "rrs1",
            "gene": "rrs",
            "product": "16S ribosomal RNA",
            "nt": s16,
        }
    ]
    genome = f">zchr\n{ZCHR.lower()}\n>apl\n{APL.lower()}\n"
    meta = json.loads(_run(reg, tmp_path, _data(feats), genome)["annotation_meta"].read_text())
    assert meta["anti_sd"] == "ACCTCCTTA"
    assert "anti_sd_default" not in meta["flags"]
