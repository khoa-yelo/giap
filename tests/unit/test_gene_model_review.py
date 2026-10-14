"""Lead decisions from the M3-01 code review (consistent with docs/specs/M2-03-*.md).

1. Strand per kind: non-coding features (incl. CRISPR, gap) with strand ``?`` or missing get
   ``.``; a CDS/sORF with a missing or ``?`` strand is a ``ValueError``.
2. ``truncated``: ``null`` / empty -> ``partial=False``; present and non-empty -> ``True``.
3. ``oriV`` is a non-coding type (``"oriv"`` in ``NcType``).
4. A CDS without ``locus`` and a top-level JSON that is not an object raise ``ValueError``.
5. A string ``db_xrefs`` is a ``ValueError`` (must be a list).
6. ``GeneModel`` is documented as unhashable and does not alias the caller's replicons dict.
7. ``overlapping``: 100k queries on a 5-Mb / 5,000-gene model in < 1.5 s, and < 3 s with one
   1-Mb gap feature present.
"""

from __future__ import annotations

import json
import random
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any, get_args

import pytest

from giap.model import Interval, Replicon
from giap.model.genes import Gene, GeneModel, NcFeature, NcType, gene_model_from_bakta_json

CIRC = {"chr1": Replicon("chr1", 50_000, "circular", "chromosome")}


def _write(path: Path, data: Any) -> Path:
    path.write_text(json.dumps(data) + "\n")
    return path


def _doc(features: list[dict[str, Any]], reps: Mapping[str, Replicon] = CIRC) -> dict[str, Any]:
    return {
        "sequences": [
            {"id": r.id, "length": r.length, "topology": r.topology, "type": r.kind}
            for r in reps.values()
        ],
        "features": features,
    }


def _load(tmp_path: Path, features: list[dict[str, Any]]) -> GeneModel:
    return gene_model_from_bakta_json(_write(tmp_path / "b.json", _doc(features)), CIRC)


def _cds(**kw: Any) -> dict[str, Any]:
    f: dict[str, Any] = {
        "type": "cds",
        "contig": "chr1",
        "start": 101,
        "stop": 400,
        "strand": "+",
        "locus": "GSIM_00005",
        "aa": "M" + "A" * 98,
        "nt": "ATG" + "GCT" * 98 + "TAA",
        "start_type": "ATG",
    }
    f.update(kw)
    return f


def _nc(t: str, **kw: Any) -> dict[str, Any]:
    f: dict[str, Any] = {"type": t, "contig": "chr1", "start": 2001, "stop": 2500, "id": "f1"}
    f.update(kw)
    return f


# --- (1) strand normalisation per kind ------------------------------------------------------


@pytest.mark.parametrize("ftype", ["crispr", "gap", "ncRNA", "tRNA", "oriC"])
@pytest.mark.parametrize("strand", ["?", None])
def test_noncoding_unknown_or_missing_strand_is_dot(
    tmp_path: Path, ftype: str, strand: str | None
) -> None:
    f = _nc(ftype) if strand is None else _nc(ftype, strand=strand)
    (nc,) = _load(tmp_path, [f]).features
    assert nc.interval == Interval("chr1", 2000, 2500, ".")


@pytest.mark.parametrize("ftype", ["cds", "sorf"])
@pytest.mark.parametrize("strand", ["?", None])
def test_gene_missing_or_unknown_strand_is_error(
    tmp_path: Path, ftype: str, strand: str | None
) -> None:
    f = _cds(type=ftype)
    if strand is None:
        del f["strand"]
    else:
        f["strand"] = strand
    with pytest.raises(ValueError, match="GSIM_00005"):
        _load(tmp_path, [f])


# --- (2) truncated -> partial ---------------------------------------------------------------


@pytest.mark.parametrize("value", [None, "", []])
def test_truncated_null_or_empty_not_partial(tmp_path: Path, value: Any) -> None:
    (g,) = _load(tmp_path, [_cds(truncated=value)]).genes
    assert g.partial is False


@pytest.mark.parametrize("value", ["5'", "3'", "both"])
def test_truncated_nonempty_is_partial(tmp_path: Path, value: str) -> None:
    (g,) = _load(tmp_path, [_cds(truncated=value)]).genes
    assert g.partial is True


def test_no_truncated_key_not_partial(tmp_path: Path) -> None:
    (g,) = _load(tmp_path, [_cds()]).genes
    assert g.partial is False


# --- (3) oriV ---------------------------------------------------------------------------------


def test_oriv_in_nc_type() -> None:
    assert "oriv" in get_args(NcType)


def test_oriv_parsed_as_noncoding(tmp_path: Path) -> None:
    (nc,) = _load(tmp_path, [_nc("oriV", strand="?", id="ov1")]).features
    assert nc == NcFeature("ov1", "oriv", Interval("chr1", 2000, 2500, "."))


# --- (4) missing locus / non-object JSON ----------------------------------------------------


def test_cds_without_locus_is_value_error(tmp_path: Path) -> None:
    f = _cds()
    del f["locus"]
    with pytest.raises(ValueError, match="locus"):
        _load(tmp_path, [f])


@pytest.mark.parametrize("data", [[], [_cds()], "x", 3, None])
def test_top_level_not_object_is_value_error(tmp_path: Path, data: Any) -> None:
    with pytest.raises(ValueError):
        gene_model_from_bakta_json(_write(tmp_path / "b.json", data), CIRC)


# --- (5) db_xrefs must be a list ------------------------------------------------------------


def test_string_db_xrefs_is_value_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="db_xrefs"):
        _load(tmp_path, [_cds(db_xrefs="SO:0001217")])


def test_list_db_xrefs_ok(tmp_path: Path) -> None:
    (g,) = _load(tmp_path, [_cds(db_xrefs=["SO:0001217"])]).genes
    assert g.db_xrefs == ("SO:0001217",)


# --- (6) GeneModel: unhashable (documented), no aliasing of replicons ----------------------


def test_gene_model_documented_unhashable() -> None:
    assert GeneModel.__doc__ is not None
    assert "unhashable" in GeneModel.__doc__.lower()


def test_direct_construction_does_not_alias_replicons() -> None:
    reps = dict(CIRC)
    g = Gene("L1", "cds", Interval("chr1", 100, 400, "+"), "MA", None, "ATG")
    m = GeneModel(replicons=reps, genes=(g,), features=())
    reps["chr1"] = Replicon("chr1", 1_000, "linear", "unknown")
    reps["p9"] = Replicon("p9", 500, "circular", "plasmid")
    assert dict(m.replicons) == CIRC
    assert m.replicons["chr1"].length == 50_000
    # Still circular: the origin-wrapping query [49_990, 50_110) covers 0..110 and reaches L1
    # at 100..110; on the aliased 1-kb linear replicon it would be out of bounds instead.
    assert m.overlapping(Interval("chr1", 49_990, 50_110, ".")) == [g]


def test_reader_does_not_alias_replicons(tmp_path: Path) -> None:
    reps = dict(CIRC)
    m = gene_model_from_bakta_json(_write(tmp_path / "b.json", _doc([_cds()])), reps)
    reps.clear()
    assert dict(m.replicons) == CIRC


# --- (7) overlapping performance -----------------------------------------------------------

_LEN = 5_000_000
_N_GENES = 5_000
_N_QUERIES = 100_000


def _big_model(with_gap: bool) -> GeneModel:
    rng = random.Random(270)
    reps = {"chr1": Replicon("chr1", _LEN, "circular", "chromosome")}
    step = _LEN // _N_GENES
    genes = []
    for i in range(_N_GENES):
        s = i * step + rng.randrange(0, 100)
        genes.append(
            Gene(
                f"G{i:05d}",
                "cds",
                Interval("chr1", s, s + rng.randrange(300, 900), rng.choice("+-")),
                "M",
                None,
                "ATG",
            )
        )
    feats = (NcFeature("gap1", "gap", Interval("chr1", 2_000_000, 3_000_000, ".")),)
    return GeneModel(replicons=reps, genes=tuple(genes), features=feats if with_gap else ())


def _queries() -> list[Interval]:
    rng = random.Random(1)
    out = []
    for _ in range(_N_QUERIES):
        s = rng.randrange(0, _LEN)
        out.append(Interval("chr1", s, s + rng.randrange(1, 2_000), "."))
    return out


def _time_queries(model: GeneModel, queries: list[Interval]) -> float:
    t0 = time.perf_counter()
    for q in queries:
        model.overlapping(q)
    return time.perf_counter() - t0


def test_overlapping_100k_queries_fast() -> None:
    model, queries = _big_model(with_gap=False), _queries()
    elapsed = _time_queries(model, queries)
    assert elapsed < 1.5, f"100k overlapping queries took {elapsed:.2f} s (limit 1.5 s)"


def test_overlapping_100k_queries_fast_with_1mb_gap() -> None:
    model, queries = _big_model(with_gap=True), _queries()
    # sanity: the gap is found and genes inside it are still reported
    hits = model.overlapping(Interval("chr1", 2_500_000, 2_500_001, "."))
    assert any(isinstance(h, NcFeature) and h.id == "gap1" for h in hits)
    elapsed = _time_queries(model, queries)
    assert elapsed < 3.0, f"100k overlapping queries with a 1-Mb gap took {elapsed:.2f} s"
