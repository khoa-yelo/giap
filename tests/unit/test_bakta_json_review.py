"""Lead decisions after the M2-03 code review (PR #30; docs/specs/M2-03-gene-model-bakta-json.md).

1. ``start_type == "Edge"`` (partial start) -> ``start_codon = None``.
2. A JSON sequence's ``length`` must equal the genome sequence used (``replicons`` or JSON ``nt``);
   a mismatch is a ``BaktaParseError`` naming the sequence.
3. ``gene_model_from_bakta_json`` delegates to ``parse_bakta_json(path, replicons).model``.
4. Type validation: non-string ``nt``/``aa``, non-int ``start``/``stop`` (bools and floats too)
   and an ``expert`` that is not a list of dicts raise ``BaktaParseError``.
5. A JSON sequence absent from a given ``replicons`` set is a ``BaktaParseError``.
6. One focused test per surviving mutant listed in the review.
"""

from __future__ import annotations

import gzip
import json
import random
from pathlib import Path
from typing import Any

import pytest

from gippa.annotate.bakta import BaktaParseError, parse_bakta_json
from gippa.annotate.bakta_cmd import bakta_db_path
from gippa.model import Interval, Replicon
from gippa.model.genes import gene_model_from_bakta_json


def _rand_seq(n: int, seed: int) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


CHR1 = _rand_seq(10_000, 11)
P1 = _rand_seq(2_000, 12)
LIN = _rand_seq(1_500, 13)

REPS = {
    "chr1": Replicon("chr1", 10_000, "circular", "chromosome"),
    "p1": Replicon("p1", 2_000, "circular", "plasmid"),
    "lin1": Replicon("lin1", 1_500, "linear", "unknown"),
}


def _seq(sid: str, nt: str | None, kind: str, topology: str | None, **kw: Any) -> dict[str, Any]:
    rec: dict[str, Any] = {"id": sid, "type": kind, "length": len(nt) if nt else 0}
    if nt is not None:
        rec["nt"] = nt
    if topology is not None:
        rec["topology"] = topology
    rec.update(kw)
    return rec


def _seqs() -> list[dict[str, Any]]:
    return [
        _seq("chr1", CHR1, "chromosome", "circular"),
        _seq("p1", P1, "plasmid", "circular"),
        _seq("lin1", LIN, "contig", "linear"),
    ]


def _doc(
    features: list[dict[str, Any]],
    sequences: list[dict[str, Any]] | None = None,
    version: Any = None,
) -> dict[str, Any]:
    doc: dict[str, Any] = {
        "features": features,
        "sequences": _seqs() if sequences is None else sequences,
    }
    doc["version"] = (
        {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}}
        if version is None
        else version
    )
    return doc


def _cds(locus: str, start: int, stop: int, strand: str, contig: str = "chr1", **kw: Any) -> Any:
    f: dict[str, Any] = {
        "type": "cds",
        "sequence": contig,
        "start": start,
        "stop": stop,
        "strand": strand,
        "locus": locus,
        "product": "DNA polymerase",
        "aa": "MKV",
    }
    f.update(kw)
    return f


def _write(tmp_path: Path, doc: Any, name: str = "b.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(doc))
    return path


def _parse(tmp_path: Path, doc: Any, **kw: Any) -> Any:
    return parse_bakta_json(_write(tmp_path, doc), **kw)


def _gene(tmp_path: Path, f: dict[str, Any], **kw: Any) -> Any:
    (g,) = _parse(tmp_path, _doc([f]), **kw).model.genes
    return g


# --- (1) start_type "Edge" -> start_codon None --------------------------------------------------


def test_edge_start_type_has_no_start_codon(tmp_path: Path) -> None:
    g = _gene(tmp_path, _cds("L_1", 1, 300, "+", contig="lin1", start_type="Edge"))
    assert g.start_codon is None
    assert g.nt == LIN[0:300]  # nt is still the feature's sequence


def test_edge_start_type_with_nt_given(tmp_path: Path) -> None:
    nt = "AAA" * 99 + "TAA"
    g = _gene(tmp_path, _cds("L_1", 101, 400, "+", nt=nt, start_type="Edge"))
    assert g.start_codon is None


def test_edge_start_type_via_gene_model_reader(tmp_path: Path) -> None:
    path = _write(tmp_path, _doc([_cds("L_1", 101, 400, "+", start_type="Edge")]))
    (g,) = gene_model_from_bakta_json(path, REPS).genes
    assert g.start_codon is None


def test_real_start_type_kept(tmp_path: Path) -> None:
    g = _gene(tmp_path, _cds("L_1", 101, 400, "+", start_type="GTG"))
    assert g.start_codon == "GTG"


# --- (2) sequence length vs genome sequence -------------------------------------------------------


def test_json_length_differs_from_nt(tmp_path: Path) -> None:
    seqs = _seqs()
    seqs[1]["length"] = len(P1) + 1
    with pytest.raises(BaktaParseError, match="p1"):
        _parse(tmp_path, _doc([], sequences=seqs))


def test_json_length_differs_from_replicons(tmp_path: Path) -> None:
    seqs = _seqs()
    for s in seqs:
        del s["nt"]
    seqs[2]["length"] = 1_499
    genome = {"chr1": CHR1, "p1": P1, "lin1": LIN}
    with pytest.raises(BaktaParseError, match="lin1"):
        _parse(tmp_path, _doc([], sequences=seqs), replicons=REPS, genome=genome)


def test_matching_lengths_parse(tmp_path: Path) -> None:
    ann = _parse(tmp_path, _doc([]), replicons=REPS)
    assert dict(ann.model.replicons) == REPS


# --- (3) one parser -----------------------------------------------------------------------------


def _rich_features() -> list[dict[str, Any]]:
    return [
        _cds("L_1", 101, 400, "+"),  # no nt/start_type: sliced from the JSON nt
        _cds("L_2", 1001, 1300, "-", start_type="Edge"),
        _cds("L_3", 9901, 150, "+", edge=True, aa="MKV*"),
        _cds("L_4", 1, 300, "+", contig="lin1"),  # linear end -> partial
        _cds("L_5", 501, 800, "+", contig="p1", pseudogene={"causes": ["frameshift"]}),
        _cds("L_6", 3001, 3060, "+", type="sorf", product="hypothetical protein"),
        {"type": "tRNA", "sequence": "chr1", "start": 2001, "stop": 2076, "strand": "+",
         "id": "t1", "anti_codon": "cag", "amino_acid": "Leu", "product": "tRNA-Leu"},
        {"type": "crispr", "sequence": "chr1", "start": 5001, "stop": 5300, "id": "c1"},
        {"type": "mystery", "sequence": "chr1", "start": 1, "stop": 10},
    ]  # fmt: skip


def test_gene_model_reader_equals_parser(tmp_path: Path) -> None:
    path = _write(tmp_path, _doc(_rich_features()))
    assert gene_model_from_bakta_json(path, REPS) == parse_bakta_json(path, REPS).model


def test_gene_model_reader_derives_start_codon_and_partial(tmp_path: Path) -> None:
    path = _write(tmp_path, _doc(_rich_features()))
    genes = {g.locus_tag: g for g in gene_model_from_bakta_json(path, REPS).genes}
    assert genes["L_1"].nt == CHR1[100:400]
    assert genes["L_1"].start_codon == CHR1[100:103]
    assert genes["L_4"].partial is True
    assert genes["L_1"].partial is False


def test_gene_model_reader_raises_bakta_parse_error(tmp_path: Path) -> None:
    path = _write(tmp_path, _doc([_cds("X_1", 101, 400, "+"), _cds("X_1", 501, 800, "+")]))
    with pytest.raises(BaktaParseError, match="X_1"):
        gene_model_from_bakta_json(path, REPS)


def test_gene_model_reader_gzip(tmp_path: Path) -> None:
    path = tmp_path / "b.json.gz"
    path.write_bytes(gzip.compress(json.dumps(_doc(_rich_features())).encode()))
    assert gene_model_from_bakta_json(path, REPS) == parse_bakta_json(path, REPS).model


# --- (4) type validation ------------------------------------------------------------------------


@pytest.mark.parametrize("value", [123, ["ATG"], {"nt": "ATG"}, True])
def test_feature_nt_not_string(tmp_path: Path, value: Any) -> None:
    with pytest.raises(BaktaParseError, match="L_1"):
        _parse(tmp_path, _doc([_cds("L_1", 101, 400, "+", nt=value)]))


@pytest.mark.parametrize("value", [123, ["M", "K"], {"aa": "MK"}, True])
def test_feature_aa_not_string(tmp_path: Path, value: Any) -> None:
    with pytest.raises(BaktaParseError, match="L_1"):
        _parse(tmp_path, _doc([_cds("L_1", 101, 400, "+", aa=value)]))


@pytest.mark.parametrize("value", [123, ["ACGT"]])
def test_sequence_nt_not_string(tmp_path: Path, value: Any) -> None:
    seqs = [_seq("chr1", CHR1, "chromosome", "circular")]
    seqs[0]["nt"] = value
    with pytest.raises(BaktaParseError, match="chr1"):
        _parse(tmp_path, _doc([], sequences=seqs))


@pytest.mark.parametrize("key", ["start", "stop"])
@pytest.mark.parametrize("value", [True, False, 101.0, 400.5, "101", None])
def test_start_stop_must_be_int(tmp_path: Path, key: str, value: Any) -> None:
    f = _cds("L_1", 101, 400, "+")
    f[key] = value
    with pytest.raises(BaktaParseError, match="L_1"):
        _parse(tmp_path, _doc([f]))


@pytest.mark.parametrize("key", ["start", "stop"])
@pytest.mark.parametrize("value", [True, 2001.0, "2001"])
def test_nc_start_stop_must_be_int(tmp_path: Path, key: str, value: Any) -> None:
    f = {"type": "tRNA", "sequence": "chr1", "start": 2001, "stop": 2076, "strand": "+", "id": "t9"}
    f[key] = value
    with pytest.raises(BaktaParseError, match="t9"):
        _parse(tmp_path, _doc([f]))


@pytest.mark.parametrize(
    "value", ["amrfinder", {"type": "amrfinder"}, ["amrfinder"], [{"type": "a"}, 3], 7]
)
def test_expert_must_be_list_of_dicts(tmp_path: Path, value: Any) -> None:
    with pytest.raises(BaktaParseError, match="L_1"):
        _parse(tmp_path, _doc([_cds("L_1", 101, 400, "+", expert=value)]))


@pytest.mark.parametrize("value", [None, []])
def test_expert_null_or_empty_ok(tmp_path: Path, value: Any) -> None:
    assert _gene(tmp_path, _cds("L_1", 101, 400, "+", expert=value)).expert == ()


# --- (5) JSON sequence missing from replicons ----------------------------------------------------


def test_json_sequence_absent_from_replicons(tmp_path: Path) -> None:
    reps = {k: v for k, v in REPS.items() if k != "p1"}
    with pytest.raises(BaktaParseError, match="p1"):
        _parse(tmp_path, _doc([]), replicons=reps)


# --- (6) surviving mutants -----------------------------------------------------------------------


def test_mutant_pseudogene_key_alone_marks_pseudogene(tmp_path: Path) -> None:
    f = _cds("L_1", 101, 400, "+", pseudogene={"causes": ["internal_stop"]})
    assert "pseudo" not in f
    assert _gene(tmp_path, f).feature_type == "pseudogene"


def test_mutant_no_pseudo_markers_is_cds(tmp_path: Path) -> None:
    assert _gene(tmp_path, _cds("L_1", 101, 400, "+", pseudo=False)).feature_type == "cds"


@pytest.mark.parametrize(("aa", "expected"), [("MKV*", "MKV"), ("MKV", "MKV"), ("M*K", "M*K")])
def test_mutant_terminal_stop_strip(tmp_path: Path, aa: str, expected: str) -> None:
    assert _gene(tmp_path, _cds("L_1", 101, 400, "+", aa=aa)).aa == expected


@pytest.mark.parametrize(
    ("extra", "expected"),
    [
        ({"hypothetical": False, "product": "hypothetical protein"}, False),
        ({"hypothetical": True, "product": "DNA polymerase"}, True),
        ({"product": "hypothetical protein YbcX"}, True),
        ({"product": "putative hypothetical protein"}, False),
        ({"product": "DNA polymerase"}, False),
    ],
)
def test_mutant_hypothetical_precedence(
    tmp_path: Path, extra: dict[str, Any], expected: bool
) -> None:
    assert _gene(tmp_path, _cds("L_1", 101, 400, "+", **extra)).hypothetical is expected


def test_mutant_unknown_topology_is_error(tmp_path: Path) -> None:
    seqs = [_seq("chr1", CHR1, "chromosome", "spiral")]
    with pytest.raises(BaktaParseError, match="chr1"):
        _parse(tmp_path, _doc([], sequences=seqs))


def test_mutant_topology_case_insensitive(tmp_path: Path) -> None:
    seqs = [_seq("chr1", CHR1, "chromosome", "Circular")]
    ann = _parse(tmp_path, _doc([], sequences=seqs))
    assert ann.model.replicons["chr1"].topology == "circular"


@pytest.mark.parametrize("strand", ["x", "+-", 1])
def test_mutant_nc_invalid_strand_is_error(tmp_path: Path, strand: Any) -> None:
    f = {"type": "tRNA", "sequence": "chr1", "start": 2001, "stop": 2076, "strand": strand,
         "id": "t9"}  # fmt: skip
    with pytest.raises(BaktaParseError, match="t9"):
        _parse(tmp_path, _doc([f]))


@pytest.mark.parametrize("strand", ["?", ".", None])
def test_mutant_nc_unknown_strand_is_dot(tmp_path: Path, strand: Any) -> None:
    f = {"type": "tRNA", "sequence": "chr1", "start": 2001, "stop": 2076, "strand": strand,
         "id": "t9"}  # fmt: skip
    (nc,) = _parse(tmp_path, _doc([f])).model.features
    assert nc.interval == Interval("chr1", 2000, 2076, ".")


@pytest.mark.parametrize(
    ("version", "expected"),
    [
        ({"bakta": "1.12.2", "db": {"version": "6", "type": "custom"}}, ("1.12.2", "6", "unknown")),
        ({"bakta": "1.12.2", "db": {"version": "6", "type": "LIGHT"}}, ("1.12.2", "6", "light")),
        ({"bakta": "1.12.2", "db": {"version": "6.0"}}, ("1.12.2", "6.0", "unknown")),
        ({"bakta": "1.12.2", "db": "6.0"}, ("1.12.2", "", "unknown")),
        ({"bakta": "1.12.2"}, ("1.12.2", "", "unknown")),
        ("1.12.2", ("", "", "unknown")),
    ],
)  # fmt: skip
def test_mutant_db_type_fallback(
    tmp_path: Path, version: Any, expected: tuple[str, str, str]
) -> None:
    ann = _parse(tmp_path, _doc([], version=version))
    assert (ann.bakta_version, ann.db_version, ann.db_type) == expected


def test_mutant_genome_fill_from_replicons(tmp_path: Path) -> None:
    """A replicon missing from the JSON gets its sequence from ``genome``."""
    seqs = [_seq("chr1", CHR1, "chromosome", "circular")]
    genome = {"chr1": CHR1, "p1": P1, "lin1": LIN}
    doc = _doc([_cds("L_1", 501, 800, "-", contig="p1")], sequences=seqs)
    ann = _parse(tmp_path, doc, replicons=REPS, genome=genome)
    assert ann.sequences["p1"] == P1
    comp = str.maketrans("ACGT", "TGCA")
    assert ann.model.genes[0].nt == P1[500:800].translate(comp)[::-1]


def test_mutant_expert_dedup_and_source_fallback(tmp_path: Path) -> None:
    expert = [
        {"type": "amrfinder", "rank": 98},
        {"type": "amrfinder", "rank": 97},
        {"source": "pgap", "rank": 90},
        {"type": "", "source": "ips"},
        {"rank": 1},
    ]
    g = _gene(tmp_path, _cds("L_1", 101, 400, "+", expert=expert))
    assert g.expert == ("amrfinder", "pgap", "ips")


@pytest.mark.parametrize(
    ("setting", "expected"),
    [
        ("/dbs/bakta-light-6.0", "light"),
        ("/dbs/light/db", "full"),  # only the directory's own name counts
        ("/dbs/bakta-6.0", "full"),
    ],
)
def test_mutant_light_name_heuristic(tmp_path: Path, setting: str, expected: str) -> None:
    path, db_type = bakta_db_path(setting, tmp_path)
    assert path == Path(setting)
    assert db_type == expected


def test_mutant_feature_sequence_key_wins_over_contig(tmp_path: Path) -> None:
    f = _cds("L_1", 101, 400, "+", contig="chr1")
    f["contig"] = "p1"
    g = _gene(tmp_path, f)
    assert g.interval == Interval("chr1", 100, 400, "+")
    assert g.nt == CHR1[100:400]


def test_mutant_sequences_key_wins_over_contigs(tmp_path: Path) -> None:
    doc = _doc([], sequences=[_seq("chr1", CHR1, "chromosome", "circular")])
    doc["contigs"] = [_seq("p1", P1, "plasmid", "circular")]
    ann = _parse(tmp_path, doc)
    assert list(ann.model.replicons) == ["chr1"]
    assert dict(ann.sequences) == {"chr1": CHR1}
