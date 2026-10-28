"""Unit tests for the defensive Bakta 1.12 JSON parser (issue M2-03; GIPPA_SPEC.md §7.1.1, §7.1.2).

M2-03 extends ``gippa.model.genes`` (M3-01). What M3-01 already pins (type mapping, ``pseudo``,
``truncated``, ``hypothetical``, stop stripping, sorting, unknown contigs) is not re-tested here;
these tests cover what M2-03 adds: ``gippa.annotate.bakta.parse_bakta_json`` on the full Bakta 1.12
JSON shape (version/DB metadata, the replicon table, sORFs, ncRNA classes, CRISPR, gaps, oriC/oriV,
``db_xrefs``, ``edge: true`` origin-crossing features), ``BaktaParseError`` and the extra
``Gene``/``NcFeature`` fields. Fixtures are hand-written, minimal Bakta-shaped JSON (plain,
gzip). The simulator round trips live in ``test_bakta_json_sim.py``.
"""

from __future__ import annotations

import copy
import dataclasses
import gzip
import json
import random
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

from gippa.model import Interval, Replicon

# ----------------------------------------------------------------------------------------------
# fixtures


def _bakta() -> ModuleType:
    import gippa.annotate.bakta as mod

    return mod


def _rand_seq(n: int, seed: int) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


_COMP = str.maketrans("ACGTN", "TGCAN")


def _revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


CHR1 = _rand_seq(100_000, 1)  # the spec's chr1: circular, 100,000 bp
P1 = _rand_seq(5_000, 2)
LIN = _rand_seq(3_000, 3)

VERSION_FULL = {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}}


def _sequence(sid: str, nt: str, kind: str, topology: str | None) -> dict[str, Any]:
    """A Bakta 1.12 ``sequences[]`` record."""
    rec: dict[str, Any] = {
        "id": sid,
        "description": f"[{sid}]",
        "nt": nt,
        "length": len(nt),
        "complete": topology == "circular",
        "type": kind,
        "simple_id": f"contig_{sid}",
        "orig_id": sid,
        "orig_description": "",
    }
    if topology is not None:
        rec["topology"] = topology
    return rec


def _sequences() -> list[dict[str, Any]]:
    return [
        _sequence("chr1", CHR1, "chromosome", "circular"),
        _sequence("p1", P1, "plasmid", "circular"),
        _sequence("lin1", LIN, "contig", "linear"),
    ]


def _doc(
    features: list[dict[str, Any]],
    *,
    sequences: list[dict[str, Any]] | None = None,
    version: dict[str, Any] | None = VERSION_FULL,
) -> dict[str, Any]:
    """A Bakta 1.12-shaped top-level document."""
    doc: dict[str, Any] = {
        "genome": {
            "genus": "Escherichia",
            "species": "coli",
            "strain": "sim",
            "complete": True,
            "gram": "-",
            "translation_table": 11,
        },
        "stats": {"no_sequences": 3, "size": len(CHR1) + len(P1) + len(LIN), "gc": 0.5},
        "features": features,
        "sequences": _sequences() if sequences is None else sequences,
        "run": {"start": "2026-10-14 10:00:00", "end": "2026-10-14 10:05:00", "duration": "5"},
    }
    if version is not None:
        doc["version"] = version
    return doc


def _cds(
    locus: str, start: int, stop: int, strand: str, contig: str = "chr1", **kw: Any
) -> dict[str, Any]:
    """A Bakta 1.12 CDS record without ``nt``/``start_type`` (so they come from the sequence)."""
    f: dict[str, Any] = {
        "type": "cds",
        "sequence": contig,
        "start": start,
        "stop": stop,
        "strand": strand,
        "frame": 1,
        "gene": None,
        "product": "hypothetical protein",
        "hypothetical": True,
        "db_xrefs": [],
        "aa": "MKV",
        "aa_hexdigest": "0" * 32,
        "rbs_motif": None,
        "rbs_spacer": None,
        "genes": [],
        "locus": locus,
        "id": f"id_{locus}",
    }
    f.update(kw)
    return f


def _nc(
    ftype: str, fid: str, start: int, stop: int, strand: str | None, **kw: Any
) -> dict[str, Any]:
    f: dict[str, Any] = {
        "type": ftype,
        "sequence": kw.pop("contig", "chr1"),
        "start": start,
        "stop": stop,
        "id": fid,
    }
    if strand is not None:
        f["strand"] = strand
    f.update(kw)
    return f


def _write(tmp_path: Path, doc: Any, name: str = "x.json", gz: bool = False) -> Path:
    path = tmp_path / (name + (".gz" if gz else ""))
    raw = json.dumps(doc).encode()
    path.write_bytes(gzip.compress(raw) if gz else raw)
    return path


def _parse(tmp_path: Path, doc: Any, **kw: Any) -> Any:
    return _bakta().parse_bakta_json(_write(tmp_path, doc), **kw)


def _genes(ann: Any) -> dict[str, Any]:
    return {g.locus_tag: g for g in ann.model.genes}


def _one_gene(tmp_path: Path, feature: dict[str, Any], **kw: Any) -> Any:
    (g,) = _parse(tmp_path, _doc([feature]), **kw).model.genes
    return g


def _one_feature(tmp_path: Path, feature: dict[str, Any]) -> Any:
    (f,) = _parse(tmp_path, _doc([feature])).model.features
    return f


# ----------------------------------------------------------------------------------------------
# public interface


def test_bakta_parse_error_is_value_error() -> None:
    mod = _bakta()
    assert issubclass(mod.BaktaParseError, ValueError)


def test_bakta_annotation_fields() -> None:
    mod = _bakta()
    names = {f.name for f in dataclasses.fields(mod.BaktaAnnotation)}
    assert {
        "model",
        "sequences",
        "bakta_version",
        "db_version",
        "db_type",
        "skipped_types",
    } <= names


def test_model_is_m3_gene_model(tmp_path: Path) -> None:
    from gippa.model.genes import Gene, GeneModel

    ann = _parse(tmp_path, _doc([_cds("GSIM_00005", 101, 400, "+")]))
    assert isinstance(ann.model, GeneModel)
    assert all(isinstance(g, Gene) for g in ann.model.genes)


def test_new_gene_and_ncfeature_fields_have_defaults() -> None:
    from gippa.model.genes import Gene, NcFeature

    g = Gene("L_1", "cds", Interval("chr1", 0, 9, "+"), "MK", None, None)
    assert g.expert == ()
    nc = NcFeature("t1", "trna", Interval("chr1", 0, 76, "+"))
    assert (nc.gene, nc.product, nc.amino_acid, nc.anticodon, nc.nt) == (None, "", None, None, None)


def test_str_path_accepted(tmp_path: Path) -> None:
    path = _write(tmp_path, _doc([_cds("GSIM_00005", 101, 400, "+")]))
    ann = _bakta().parse_bakta_json(str(path))
    assert [g.locus_tag for g in ann.model.genes] == ["GSIM_00005"]


def test_gzip_json_same_as_plain(tmp_path: Path) -> None:
    doc = _doc(
        [
            _cds("GSIM_00005", 101, 400, "+"),
            _cds("GSIM_00010", 99901, 150, "-", edge=True),
            _nc("tRNA", "t1", 2001, 2076, "+", anti_codon="cag", amino_acid="Leu"),
        ]
    )
    mod = _bakta()
    plain = mod.parse_bakta_json(_write(tmp_path, doc, "plain.json"))
    gz = mod.parse_bakta_json(_write(tmp_path, doc, "packed.json", gz=True))
    assert gz.model.genes == plain.model.genes
    assert gz.model.features == plain.model.features
    assert dict(gz.sequences) == dict(plain.sequences)
    assert (gz.bakta_version, gz.db_version, gz.db_type) == ("1.12.2", "6.0", "full")


# ----------------------------------------------------------------------------------------------
# spec examples 1-9


def test_example1_cds_plus(tmp_path: Path) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00005", 101, 400, "+"))
    assert g.locus_tag == "GSIM_00005"
    assert g.feature_type == "cds"
    assert g.interval == Interval("chr1", 100, 400, "+")
    rep = Replicon("chr1", 100_000, "circular", "chromosome")
    assert not g.interval.crosses_origin(rep)


def test_example2_contig_key_identical_gene(tmp_path: Path) -> None:
    a = _cds("GSIM_00005", 101, 400, "+")
    b = copy.deepcopy(a)
    b["contig"] = b.pop("sequence")
    assert _one_gene(tmp_path, a) == _one_gene(tmp_path, b)


def test_example3_minus_strand_nt_and_start_codon_from_sequence(tmp_path: Path) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00010", 1001, 1300, "-"))
    assert g.interval == Interval("chr1", 1000, 1300, "-")
    expected = _revcomp(CHR1[1000:1300])
    assert g.nt == expected
    assert g.start_codon == expected[:3]


def test_example4_edge_feature_unwrapped(tmp_path: Path) -> None:
    ann = _parse(tmp_path, _doc([_cds("GSIM_00015", 99901, 150, "+", edge=True)]))
    (g,) = ann.model.genes
    assert g.interval == Interval("chr1", 99900, 100150, "+")
    assert g.interval.crosses_origin(ann.model.replicons["chr1"])
    expected = CHR1[99900:] + CHR1[:150]
    assert g.nt == expected
    assert g.start_codon == expected[:3]


def test_edge_feature_minus_strand_nt_origin_aware(tmp_path: Path) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00015", 99901, 150, "-", edge=True))
    assert g.interval == Interval("chr1", 99900, 100150, "-")
    assert g.nt == _revcomp(CHR1[99900:] + CHR1[:150])


def test_example5_edge_on_linear_contig_names_locus(tmp_path: Path) -> None:
    mod = _bakta()
    seqs = [_sequence("chr1", CHR1, "chromosome", "linear")]
    doc = _doc([_cds("GSIM_00015", 99901, 150, "+", edge=True)], sequences=seqs)
    with pytest.raises(mod.BaktaParseError, match="GSIM_00015"):
        _parse(tmp_path, doc)


def test_example6_pseudo_flag(tmp_path: Path) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00020", 101, 400, "+", pseudo=True))
    assert g.feature_type == "pseudogene"


def test_pseudogene_without_aa_gets_empty_string(tmp_path: Path) -> None:
    f = _cds("GSIM_00020", 101, 400, "+", pseudo=True)
    del f["aa"]
    del f["aa_hexdigest"]
    f["pseudogene"] = {"inference": "HMM", "causes": ["internal_stop"], "observations": []}
    g = _one_gene(tmp_path, f)
    assert g.feature_type == "pseudogene"
    assert g.aa == ""
    assert g.nt == CHR1[100:400]


def test_pseudogene_keeps_given_aa(tmp_path: Path) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00020", 101, 400, "+", pseudo=True, aa="MKVL"))
    assert g.aa == "MKVL"


def test_example7_trna(tmp_path: Path) -> None:
    f = _one_feature(
        tmp_path,
        _nc(
            "tRNA",
            "t1",
            2001,
            2087,
            "+",
            gene="leuW",
            product="tRNA-Leu",
            anti_codon="cag",
            amino_acid="Leu",
            locus="GSIM_00025",
        ),
    )
    assert f.type == "trna"
    assert f.anticodon == "CAG"
    assert f.amino_acid == "Leu"
    assert f.gene == "leuW"
    assert f.product == "tRNA-Leu"
    assert f.interval == Interval("chr1", 2000, 2087, "+")


def test_trna_anticodon_u_to_t(tmp_path: Path) -> None:
    f = _one_feature(
        tmp_path, _nc("tRNA", "t2", 2001, 2076, "-", anti_codon="ugc", amino_acid="Ala")
    )
    assert f.anticodon == "TGC"


def test_example8_versions(tmp_path: Path) -> None:
    version = {"bakta": "1.12.2", "db": {"version": "6.0", "type": "light"}}
    ann = _parse(tmp_path, _doc([], version=version))
    assert (ann.bakta_version, ann.db_version, ann.db_type) == ("1.12.2", "6.0", "light")


def test_db_type_full(tmp_path: Path) -> None:
    ann = _parse(tmp_path, _doc([]))
    assert (ann.bakta_version, ann.db_version, ann.db_type) == ("1.12.2", "6.0", "full")


def test_missing_version_block_is_unknown_db_type(tmp_path: Path) -> None:
    ann = _parse(tmp_path, _doc([_cds("GSIM_00005", 101, 400, "+")], version=None))
    assert ann.db_type == "unknown"
    assert len(ann.model.genes) == 1


def test_example9_duplicate_locus_tag(tmp_path: Path) -> None:
    mod = _bakta()
    doc = _doc([_cds("X_1", 101, 400, "+"), _cds("X_1", 1001, 1300, "-")])
    with pytest.raises(mod.BaktaParseError, match="duplicate locus tag X_1"):
        _parse(tmp_path, doc)


def test_duplicate_locus_across_cds_and_sorf(tmp_path: Path) -> None:
    mod = _bakta()
    sorf = _cds("X_1", 1001, 1060, "+", type="sorf")
    with pytest.raises(mod.BaktaParseError, match="X_1"):
        _parse(tmp_path, _doc([_cds("X_1", 101, 400, "+"), sorf]))


# ----------------------------------------------------------------------------------------------
# contig resolution and the replicon table


def test_feature_without_contig_key_is_error(tmp_path: Path) -> None:
    mod = _bakta()
    f = _cds("GSIM_00005", 101, 400, "+")
    del f["sequence"]
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, _doc([f]))


def test_replicon_table_from_sequences(tmp_path: Path) -> None:
    seqs = [*_sequences(), _sequence("x9", "ACGT" * 50, "contig", None)]
    ann = _parse(tmp_path, _doc([], sequences=seqs))
    assert list(ann.model.replicons) == ["chr1", "p1", "lin1", "x9"]
    assert ann.model.replicons["chr1"] == Replicon("chr1", 100_000, "circular", "chromosome")
    assert ann.model.replicons["p1"] == Replicon("p1", 5_000, "circular", "plasmid")
    assert ann.model.replicons["lin1"] == Replicon("lin1", 3_000, "linear", "unknown")
    # topology missing -> linear
    assert ann.model.replicons["x9"] == Replicon("x9", 200, "linear", "unknown")
    assert dict(ann.sequences) == {"chr1": CHR1, "p1": P1, "lin1": LIN, "x9": "ACGT" * 50}


def test_contigs_key_and_sequence_nucleotides_fallback(tmp_path: Path) -> None:
    rec = {"id": "chr1", "type": "chromosome", "topology": "circular", "sequence": CHR1}
    doc = _doc([_cds("GSIM_00005", 101, 400, "+")], sequences=[])
    del doc["sequences"]
    doc["contigs"] = [rec]
    ann = _parse(tmp_path, doc)
    assert ann.model.replicons["chr1"] == Replicon("chr1", 100_000, "circular", "chromosome")
    assert dict(ann.sequences) == {"chr1": CHR1}
    assert ann.model.genes[0].nt == CHR1[100:400]


def test_replicons_argument_wins(tmp_path: Path) -> None:
    reps = {
        "chr1": Replicon("chr1", 100_000, "circular", "chromosome"),
        "p1": Replicon("p1", 5_000, "linear", "unknown"),  # JSON says circular plasmid
        "lin1": Replicon("lin1", 3_000, "linear", "unknown"),
    }
    ann = _parse(tmp_path, _doc([]), replicons=reps)
    assert dict(ann.model.replicons) == reps


def test_replicons_argument_topology_governs_edge(tmp_path: Path) -> None:
    mod = _bakta()
    reps = {
        "chr1": Replicon("chr1", 100_000, "linear", "chromosome"),  # JSON says circular
        "p1": Replicon("p1", 5_000, "circular", "plasmid"),
        "lin1": Replicon("lin1", 3_000, "linear", "unknown"),
    }
    doc = _doc([_cds("GSIM_00015", 99901, 150, "+", edge=True)])
    with pytest.raises(mod.BaktaParseError, match="GSIM_00015"):
        _parse(tmp_path, doc, replicons=reps)


def test_replicons_argument_length_mismatch(tmp_path: Path) -> None:
    mod = _bakta()
    reps = {
        "chr1": Replicon("chr1", 99_999, "circular", "chromosome"),
        "p1": Replicon("p1", 5_000, "circular", "plasmid"),
        "lin1": Replicon("lin1", 3_000, "linear", "unknown"),
    }
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, _doc([]), replicons=reps)


def test_feature_on_contig_missing_from_replicons_argument(tmp_path: Path) -> None:
    mod = _bakta()
    reps = {"chr1": Replicon("chr1", 100_000, "circular", "chromosome")}
    seqs = [_sequence("chr1", CHR1, "chromosome", "circular")]
    doc = _doc([_cds("GSIM_00005", 101, 400, "+", contig="chrX")], sequences=seqs)
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, doc, replicons=reps)


def test_trimmed_sequences_need_genome(tmp_path: Path) -> None:
    mod = _bakta()
    seqs = _sequences()
    for s in seqs:
        del s["nt"]
    doc = _doc([_cds("GSIM_00010", 1001, 1300, "-")], sequences=seqs)
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, doc)
    genome = {"chr1": CHR1, "p1": P1, "lin1": LIN}
    ann = _parse(tmp_path, doc, genome=genome)
    assert dict(ann.sequences) == genome
    (g,) = ann.model.genes
    assert g.nt == _revcomp(CHR1[1000:1300])


# ----------------------------------------------------------------------------------------------
# gene records


def test_nt_and_start_type_from_feature_win(tmp_path: Path) -> None:
    nt = "GTG" + "AAA" * 98 + "TAA"
    g = _one_gene(tmp_path, _cds("GSIM_00005", 101, 400, "+", nt=nt, start_type="GTG"))
    assert g.nt == nt
    assert g.start_codon == "GTG"


def test_start_codon_from_nt_when_start_type_absent(tmp_path: Path) -> None:
    nt = "TTG" + "AAA" * 98 + "TAA"
    g = _one_gene(tmp_path, _cds("GSIM_00005", 101, 400, "+", nt=nt))
    assert g.start_codon == "TTG"


def test_plus_strand_nt_sliced_from_sequence(tmp_path: Path) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00005", 101, 400, "+"))
    assert g.nt == CHR1[100:400]
    assert g.start_codon == CHR1[100:103]


def test_sorf(tmp_path: Path) -> None:
    f = _cds(
        "GSIM_00030",
        3001,
        3090,
        "+",
        type="sorf",
        product="small protein",
        hypothetical=False,
        db_xrefs=["SO:0001217", "UniRef:UniRef100_A0A0A0"],
    )
    g = _one_gene(tmp_path, f)
    assert g.feature_type == "sorf"
    assert g.interval == Interval("chr1", 3000, 3090, "+")
    assert g.db_xrefs == ("SO:0001217", "UniRef:UniRef100_A0A0A0")


def test_db_xrefs_cog_go_ec_kegg_kept_in_order(tmp_path: Path) -> None:
    xrefs = [
        "COG:COG0587",
        "COG:L",
        "EC:2.7.7.7",
        "GO:0003887",
        "GO:0006260",
        "KEGG:K02337",
        "RefSeq:WP_000176290.1",
        "SO:0001217",
        "UniParc:UPI000000F8E5",
        "UniRef:UniRef50_P10443",
    ]
    f = _cds("GSIM_00005", 101, 400, "+", gene="dnaE", product="DNA polymerase III subunit alpha")
    f["hypothetical"] = False
    f["db_xrefs"] = xrefs
    g = _one_gene(tmp_path, f)
    assert g.db_xrefs == tuple(xrefs)
    assert g.gene == "dnaE"
    assert g.product == "DNA polymerase III subunit alpha"
    assert g.hypothetical is False


def test_amrfinder_expert(tmp_path: Path) -> None:
    expert = [
        {
            "type": "amrfinder",
            "rank": 98,
            "gene": "blaTEM-1",
            "product": "class A beta-lactamase TEM-1",
            "method": "EXACTP",
            "identity": 1.0,
            "query_cov": 1.0,
            "subject_cov": 1.0,
            "db_xrefs": ["NCBIProtein:WP_000027057.1"],
        }
    ]
    g = _one_gene(tmp_path, _cds("GSIM_00005", 101, 400, "+", expert=expert))
    assert "amrfinder" in g.expert
    plain = _one_gene(tmp_path, _cds("GSIM_00010", 101, 400, "+"))
    assert plain.expert == ()


@pytest.mark.parametrize(
    ("contig", "start", "stop", "expected"),
    [
        ("lin1", 1, 300, True),  # touches the linear start
        ("lin1", 2701, 3000, True),  # touches the linear end
        ("lin1", 1001, 1300, False),
        ("chr1", 1, 300, False),  # circular contigs have no ends
        ("chr1", 99701, 100_000, False),
    ],
)
def test_partial_at_linear_contig_ends(
    tmp_path: Path, contig: str, start: int, stop: int, expected: bool
) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00005", start, stop, "+", contig=contig))
    assert g.partial is expected


@pytest.mark.parametrize(("truncated", "expected"), [("5'", True), ("3'", True), ("", False)])
def test_truncated_value(tmp_path: Path, truncated: str, expected: bool) -> None:
    g = _one_gene(tmp_path, _cds("GSIM_00005", 1001, 1300, "+", truncated=truncated))
    assert g.partial is expected


@pytest.mark.parametrize("strand", ["?", None])
def test_gene_without_strand_is_error(tmp_path: Path, strand: str | None) -> None:
    mod = _bakta()
    f = _cds("GSIM_00005", 101, 400, "+")
    if strand is None:
        del f["strand"]
    else:
        f["strand"] = strand
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, _doc([f]))


def test_cds_without_locus_is_fatal(tmp_path: Path) -> None:
    mod = _bakta()
    f = _cds("GSIM_00005", 101, 400, "+")
    del f["locus"]
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, _doc([f]))


# ----------------------------------------------------------------------------------------------
# non-coding features


def test_rrna_keeps_nt_gene_product(tmp_path: Path) -> None:
    nt = CHR1[5000:6542]
    f = _one_feature(
        tmp_path,
        _nc(
            "rRNA",
            "r1",
            5001,
            6542,
            "+",
            gene="rrs",
            product="16S ribosomal RNA",
            nt=nt,
            db_xrefs=["RFAM:RF00177", "SO:0001000"],
            locus="GSIM_00035",
        ),
    )
    assert f.type == "rrna"
    assert f.gene == "rrs"
    assert f.product == "16S ribosomal RNA"
    assert f.nt == nt
    assert f.interval == Interval("chr1", 5000, 6542, "+")


@pytest.mark.parametrize(
    ("bakta_type", "extra", "expected"),
    [
        ("ncRNA", {"class": "riboswitch", "product": "FMN riboswitch", "gene": None}, "ncrna"),
        ("ncRNA", {"class": "sRNA", "product": "6S RNA", "gene": "ssrS"}, "ncrna"),
        ("ncRNA", {"class": "antisense", "product": "SymR", "gene": "symR"}, "ncrna"),
        ("ncRNA-region", {"class": "RNA_thermometer", "product": "ROSE"}, "ncrna_region"),
        ("ncRNA-region", {"class": "ribosomal_frameshift", "product": "PRF"}, "ncrna_region"),
        ("tmRNA", {"gene": "ssrA", "product": "transfer-messenger RNA"}, "tmrna"),
    ],
)
def test_ncrna_classes(
    tmp_path: Path, bakta_type: str, extra: dict[str, Any], expected: str
) -> None:
    f = _one_feature(tmp_path, _nc(bakta_type, "n1", 7001, 7100, "-", **extra))
    assert f.type == expected
    assert f.product == extra["product"]
    assert f.gene == extra.get("gene")
    assert f.interval == Interval("chr1", 7000, 7100, "-")


@pytest.mark.parametrize("strand", ["?", None])
def test_crispr_without_strand_is_unstranded(tmp_path: Path, strand: str | None) -> None:
    f = _one_feature(
        tmp_path,
        _nc(
            "crispr",
            "c1",
            8001,
            8400,
            strand,
            product="CRISPR array with 6 repeats of length 29, consensus sequence ...",
            repeat_consensus="GTGTTCCCCGCGCCAGCGGGGATAAACCG",
            spacers=[{"start": 8030, "stop": 8061}],
            repeats=[{"start": 8001, "stop": 8029}],
        ),
    )
    assert f.type == "crispr"
    assert f.interval == Interval("chr1", 8000, 8400, ".")
    assert f.product.startswith("CRISPR array")


@pytest.mark.parametrize("strand", ["?", None])
def test_gap_without_strand_is_unstranded(tmp_path: Path, strand: str | None) -> None:
    f = _one_feature(tmp_path, _nc("gap", "g1", 9001, 9100, strand, length=100))
    assert f.type == "gap"
    assert f.interval == Interval("chr1", 9000, 9100, ".")


@pytest.mark.parametrize(
    ("bakta_type", "contig", "expected"),
    [("oriC", "chr1", "oric"), ("oriV", "p1", "oriv")],
)
def test_origins(tmp_path: Path, bakta_type: str, contig: str, expected: str) -> None:
    f = _one_feature(tmp_path, _nc(bakta_type, "o1", 101, 500, "?", contig=contig))
    assert f.type == expected
    assert f.interval == Interval(contig, 100, 500, ".")


def test_edge_nc_feature_unwrapped(tmp_path: Path) -> None:
    f = _one_feature(tmp_path, _nc("oriC", "o1", 99_801, 200, "?", edge=True))
    assert f.interval == Interval("chr1", 99_800, 100_200, ".")


# ----------------------------------------------------------------------------------------------
# document-level behaviour


def test_unknown_types_counted_not_fatal(tmp_path: Path) -> None:
    doc = _doc(
        [
            _cds("GSIM_00005", 101, 400, "+"),
            _nc("mystery", "m1", 1001, 1100, "+"),
            _nc("mystery", "m2", 2001, 2100, "-"),
            _nc("tRNA", "t1", 3001, 3076, "+", anti_codon="cag", amino_acid="Leu"),
        ]
    )
    ann = _parse(tmp_path, doc)
    assert dict(ann.skipped_types) == {"mystery": 2}
    assert [g.locus_tag for g in ann.model.genes] == ["GSIM_00005"]
    assert [f.id for f in ann.model.features] == ["t1"]


def test_known_types_not_counted_as_skipped(tmp_path: Path) -> None:
    ann = _parse(tmp_path, _doc([_cds("GSIM_00005", 101, 400, "+")]))
    assert dict(ann.skipped_types) == {}


def test_empty_features_valid(tmp_path: Path) -> None:
    ann = _parse(tmp_path, _doc([]))
    assert ann.model.genes == ()
    assert ann.model.features == ()
    assert len(ann.model.replicons) == 3


def test_top_level_must_be_object(tmp_path: Path) -> None:
    mod = _bakta()
    with pytest.raises(mod.BaktaParseError):
        _parse(tmp_path, [_cds("GSIM_00005", 101, 400, "+")])


def test_genes_and_features_sorted(tmp_path: Path) -> None:
    doc = _doc(
        [
            _cds("L_40", 101, 400, "+", contig="p1"),
            _cds("L_30", 2001, 2300, "-"),
            _cds("L_20", 101, 400, "-"),
            _cds("L_10", 101, 400, "+"),
            _nc("tRNA", "t2", 4001, 4076, "+", anti_codon="cag", amino_acid="Leu"),
            _nc("tRNA", "t1", 1001, 1076, "+", anti_codon="cag", amino_acid="Leu"),
        ]
    )
    ann = _parse(tmp_path, doc)
    assert [g.locus_tag for g in ann.model.genes] == ["L_10", "L_20", "L_30", "L_40"]
    assert [f.id for f in ann.model.features] == ["t1", "t2"]


def test_full_document_round_trip(tmp_path: Path) -> None:
    """One document with every feature class Bakta 1.12 writes."""
    doc = _doc(
        [
            _cds("GSIM_00005", 101, 400, "+", db_xrefs=["COG:COG0001", "COG:H", "EC:5.4.3.8"]),
            _cds("GSIM_00010", 99901, 150, "-", edge=True),
            _cds("GSIM_00015", 3001, 3090, "+", type="sorf"),
            _cds("GSIM_00020", 4001, 4300, "+", pseudo=True),
            _cds("GSIM_00025", 101, 400, "+", contig="p1"),
            _cds("GSIM_00030", 1, 300, "-", contig="lin1"),
            _nc("rRNA", "r1", 5001, 6542, "+", gene="rrs", product="16S ribosomal RNA"),
            _nc("tRNA", "t1", 7001, 7076, "-", anti_codon="ugc", amino_acid="Ala"),
            _nc("tmRNA", "tm1", 7201, 7563, "+", gene="ssrA"),
            _nc("ncRNA", "n1", 7601, 7700, "+", **{"class": "riboswitch"}),
            _nc("ncRNA-region", "nr1", 7801, 7900, "+", **{"class": "RNA_thermometer"}),
            _nc("crispr", "c1", 8001, 8400, "?"),
            _nc("gap", "g1", 9001, 9100, None),
            _nc("oriC", "o1", 99_801, 99_900, "?"),
            _nc("oriV", "o2", 1001, 1100, "?", contig="p1"),
        ]
    )
    ann = _parse(tmp_path, doc)
    genes = _genes(ann)
    assert {k: g.feature_type for k, g in genes.items()} == {
        "GSIM_00005": "cds",
        "GSIM_00010": "cds",
        "GSIM_00015": "sorf",
        "GSIM_00020": "pseudogene",
        "GSIM_00025": "cds",
        "GSIM_00030": "cds",
    }
    reps = ann.model.replicons
    assert [k for k, g in genes.items() if g.interval.crosses_origin(reps[g.interval.contig])] == [
        "GSIM_00010"
    ]
    assert genes["GSIM_00030"].partial is True
    assert genes["GSIM_00030"].nt == _revcomp(LIN[0:300])
    assert genes["GSIM_00025"].nt == P1[100:400]
    assert sorted(f.type for f in ann.model.features) == sorted(
        ["rrna", "trna", "tmrna", "ncrna", "ncrna_region", "crispr", "gap", "oric", "oriv"]
    )
    assert dict(ann.skipped_types) == {}
    hits = ann.model.overlapping(Interval("chr1", 0, 10, "+"))
    assert [h.locus_tag for h in hits] == ["GSIM_00010"]
