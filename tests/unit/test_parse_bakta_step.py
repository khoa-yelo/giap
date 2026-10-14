"""The ``parse_bakta`` step (M2-04; GIAP_SPEC.md §6.2, §7.1.2, §9.3, §9.4, §6.5, §6.6).

Runs the registered step through ``giap.steps.registry.run_step`` (M2-02) on a hand-written Bakta
JSON and on giapsim's ``bakta_like.json`` (M3-18). Work-directory outputs are internal, so Parquet
coordinates stay 0-based half-open (Bakta 1-based closed ``start`` → ``start - 1``, ``stop`` →
``end``); only user-facing writers convert to 1-based (§6.5).

Decisions pinned here (ambiguities in the spec):
- Unknown Bakta feature types (M2-03 ``skipped_types``) are reported: ``annotation_meta.json`` has
  ``skipped_types`` (type → count, ``{}`` when none) and the step logs a WARNING naming each type.
- ncRNA ``class`` is not needed (M2-03 decision 7): ``features.parquet`` has exactly the
  spec columns.
- ``rrna16s.fa`` sequences come from the feature's ``nt`` or, when absent (giapsim, recorded JSON
  without ``nt``), from the genome, strand-aware.
- Outputs do not depend on ``cfg.seed`` (the step draws no randomness).
"""

from __future__ import annotations

import csv
import importlib
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow.parquet as pq
import pytest

from giap.config import GiapConfig

OUT_FILES = {
    "genes": "genes.parquet",
    "features": "features.parquet",
    "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa",
    "trna_counts": "trna_counts.tsv",
    "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}
GENE_COLUMNS = [
    "locus_tag",
    "gene",
    "product",
    "feature_type",
    "contig",
    "start",
    "end",
    "strand",
    "partial",
    "crosses_origin",
    "hypothetical",
    "aa",
    "nt",
    "start_codon",
    "ribosomal_protein",
    "amr",
]
FEATURE_COLUMNS = [
    "feature_id",
    "type",
    "contig",
    "start",
    "end",
    "strand",
    "gene",
    "product",
    "amino_acid",
    "anticodon",
]
META_KEYS = {
    "bakta_version",
    "db_version",
    "db_type",
    "annotation_source",
    "anti_sd",
    "flags",
    "n_genes",
    "n_features",
    "skipped_types",
}

_COMP = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


# ------------------------------------------------------------------ registry access


@pytest.fixture
def reg() -> Any:
    """The real registry with ``parse_bakta`` (robust to earlier temporary registries)."""
    from giap.steps import registry

    mod = importlib.import_module("giap.steps.annotation")
    if "parse_bakta" not in registry._REGISTRY:
        importlib.reload(mod)
    return registry


def _run(
    reg: Any, bakta_json: Path, genome: Path, outdir: Path, cfg: GiapConfig | None = None
) -> dict[str, Path]:
    outs = {k: outdir / v for k, v in OUT_FILES.items()}
    reg.run_step(
        "parse_bakta",
        cfg or GiapConfig(seed=1),
        {"bakta_json": bakta_json, "genome": genome},
        outs,
    )
    return outs


def test_step_declaration(reg: Any) -> None:
    spec = reg.get("parse_bakta")
    assert spec.kind == "py"
    assert spec.module == "core"
    assert spec.inputs == ("bakta_json", "genome")
    assert spec.outputs == tuple(OUT_FILES)
    assert spec.optional_inputs == ()
    assert spec.precomputed is None  # precomputed.bakta_json replaces `bakta`, never parse_bakta


# ------------------------------------------------------------------ hand-written input

CHR_LEN, PL_LEN = 3_000, 1_200
_rng = np.random.default_rng(20261014)
_CHR = list("".join(_rng.choice(list("ACGT"), CHR_LEN)))
# 16S rRNA on the minus strand at Bakta 1001..1600; its 5'→3' sequence ends in the E. coli tail
_CHR[1000:1600] = list(_rc("".join(_CHR[1000:1587]) + "GATCACCTCCTTA"))
CHR = "".join(_CHR)
PL = "".join(_rng.choice(list("ACGT"), PL_LEN))
S16 = _rc(CHR[1000:1600])
assert S16.endswith("GATCACCTCCTTA")


def _cds(locus: str, start: int, stop: int, strand: str, **kw: Any) -> dict[str, Any]:
    return {
        "type": "cds",
        "sequence": "chr1",
        "start": start,
        "stop": stop,
        "strand": strand,
        "locus": locus,
        "aa": "MKV",
        **kw,
    }


def _bakta_dict(*, db_type: str = "light", with_nt: bool = True) -> dict[str, Any]:
    seqs = [
        # JSON lists the plasmid first; replicons.parquet follows the FASTA order
        {"id": "p1", "type": "plasmid", "topology": "circular", "length": PL_LEN},
        {"id": "chr1", "type": "chromosome", "topology": "circular", "length": CHR_LEN},
    ]
    if with_nt:
        seqs[0]["nt"], seqs[1]["nt"] = PL, CHR
    feats = [
        _cds(
            "T_00010",
            101,
            400,
            "+",
            gene="rplB",
            product="50S ribosomal protein L2",
            db_xrefs=["COG:COG0090", "GO:0003735", "KEGG:K02886"],
        ),
        _cds(
            "T_00020",
            1701,
            2000,
            "-",
            product="hypothetical protein",
            expert=[{"type": "amrfinder", "rank": 1}],
        ),
        _cds(
            "T_00030",
            2901,
            99,
            "+",
            edge=True,
            product="DNA gyrase subunit B",
            gene="gyrB",
            db_xrefs=["EC:5.6.2.2"],
        ),
        {
            "type": "cds",
            "sequence": "p1",
            "start": 11,
            "stop": 310,
            "strand": "+",
            "locus": "T_00040",
            "aa": "MKV",
            "product": "replication protein",
            "pseudo": True,
        },
        {
            "type": "rRNA",
            "sequence": "chr1",
            "start": 1001,
            "stop": 1600,
            "strand": "-",
            "id": "T_rrs1",
            "gene": "rrs",
            "product": "16S ribosomal RNA",
        },
        {
            "type": "tRNA",
            "sequence": "chr1",
            "start": 2101,
            "stop": 2176,
            "strand": "+",
            "id": "T_t1",
            "product": "tRNA-Leu",
            "amino_acid": "Leu",
            "anti_codon": "cag",
        },
        {
            "type": "tRNA",
            "sequence": "chr1",
            "start": 2201,
            "stop": 2276,
            "strand": "-",
            "id": "T_t2",
            "product": "tRNA-Leu",
            "amino_acid": "Leu",
            "anti_codon": "CAG",
        },
        {
            "type": "tRNA",
            "sequence": "p1",
            "start": 501,
            "stop": 576,
            "strand": "+",
            "id": "T_t3",
            "product": "tRNA-Gly",
            "amino_acid": "Gly",
            "anti_codon": "gcc",
        },
        {
            "type": "crispr",
            "sequence": "chr1",
            "start": 2401,
            "stop": 2600,
            "strand": "?",
            "id": "T_cr1",
            "product": "CRISPR array",
        },
        {
            "type": "oriC",
            "sequence": "chr1",
            "start": 2951,
            "stop": 2990,
            "strand": "+",
            "id": "T_oriC",
        },
        {"type": "future_type", "sequence": "chr1", "start": 5, "stop": 50, "strand": "+"},
        {"type": "future_type", "sequence": "p1", "start": 5, "stop": 50, "strand": "+"},
    ]
    data: dict[str, Any] = {"sequences": seqs, "features": feats}
    if db_type != "none":
        data["version"] = {"bakta": "1.12.2", "db": {"version": "6.0", "type": db_type}}
    return data


def _write_inputs(d: Path, **kw: Any) -> tuple[Path, Path]:
    d.mkdir(parents=True, exist_ok=True)
    js, fa = d / "bakta.json", d / "genome.norm.fa"
    js.write_text(json.dumps(_bakta_dict(**kw)))
    fa.write_text(f">chr1\n{CHR}\n>p1\n{PL}\n")
    return js, fa


@pytest.fixture
def ran(reg: Any, tmp_path: Path) -> dict[str, Path]:
    js, fa = _write_inputs(tmp_path / "in")
    return _run(reg, js, fa, tmp_path / "work")


def _rows(path: Path) -> list[dict[str, Any]]:
    return pq.read_table(path).to_pylist()


def _meta(outs: dict[str, Path]) -> dict[str, Any]:
    meta: dict[str, Any] = json.loads(outs["annotation_meta"].read_text())
    return meta


def _fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, list[str]]] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            recs.append((line[1:].strip(), []))
        elif line:
            recs[-1][1].append(line.strip())
    return [(h, "".join(s)) for h, s in recs]


def test_every_output_written(ran: dict[str, Path]) -> None:
    for key, path in ran.items():
        assert path.is_file(), key
    assert not list(path.parent.glob(".tmp.*"))


def test_genes_parquet(ran: dict[str, Path]) -> None:
    table = pq.read_table(ran["genes"])
    assert table.column_names == GENE_COLUMNS
    rows = {r["locus_tag"]: r for r in table.to_pylist()}
    assert list(rows) == ["T_00010", "T_00020", "T_00030", "T_00040"]  # FASTA order, then start
    g1 = rows["T_00010"]
    assert (g1["contig"], g1["start"], g1["end"], g1["strand"]) == ("chr1", 100, 400, "+")
    assert g1["gene"] == "rplB" and g1["ribosomal_protein"] is True and g1["amr"] is False
    assert g1["feature_type"] == "cds" and g1["hypothetical"] is False
    assert g1["nt"] == CHR[100:400] and g1["start_codon"] == CHR[100:103] and g1["aa"] == "MKV"
    g2 = rows["T_00020"]
    assert (g2["start"], g2["end"], g2["strand"]) == (1700, 2000, "-")
    assert g2["nt"] == _rc(CHR[1700:2000])
    assert g2["amr"] is True and g2["hypothetical"] is True and g2["ribosomal_protein"] is False
    assert g2["gene"] is None
    g3 = rows["T_00030"]
    assert (g3["start"], g3["end"], g3["crosses_origin"]) == (2900, 3099, True)
    assert g3["nt"] == CHR[2900:] + CHR[:99]
    assert rows["T_00010"]["crosses_origin"] is False
    g4 = rows["T_00040"]
    assert (g4["contig"], g4["start"], g4["end"], g4["feature_type"]) == (
        "p1",
        10,
        310,
        "pseudogene",
    )
    assert all(r["partial"] is False for r in rows.values())


def test_features_parquet(ran: dict[str, Path]) -> None:
    table = pq.read_table(ran["features"])
    assert table.column_names == FEATURE_COLUMNS  # no ncRNA `class` (M2-03 decision 7)
    rows = {r["feature_id"]: r for r in table.to_pylist()}
    assert list(rows) == ["T_rrs1", "T_t1", "T_t2", "T_cr1", "T_oriC", "T_t3"]
    assert rows["T_rrs1"] == {
        "feature_id": "T_rrs1",
        "type": "rrna",
        "contig": "chr1",
        "start": 1000,
        "end": 1600,
        "strand": "-",
        "gene": "rrs",
        "product": "16S ribosomal RNA",
        "amino_acid": None,
        "anticodon": None,
    }
    assert (rows["T_t1"]["anticodon"], rows["T_t1"]["amino_acid"]) == ("CAG", "Leu")
    assert rows["T_t1"]["type"] == "trna"
    assert (rows["T_cr1"]["type"], rows["T_cr1"]["strand"]) == ("crispr", ".")
    assert (rows["T_oriC"]["type"], rows["T_oriC"]["start"], rows["T_oriC"]["end"]) == (
        "oric",
        2950,
        2990,
    )


def test_replicons_parquet_in_fasta_order(ran: dict[str, Path]) -> None:
    table = pq.read_table(ran["replicons"])
    assert table.column_names == ["contig", "length", "topology", "kind"]
    assert table.to_pylist() == [
        {"contig": "chr1", "length": CHR_LEN, "topology": "circular", "kind": "chromosome"},
        {"contig": "p1", "length": PL_LEN, "topology": "circular", "kind": "plasmid"},
    ]


def test_rrna16s_fasta(ran: dict[str, Path]) -> None:
    assert _fasta(ran["rrna16s"]) == [("T_rrs1", S16)]


def test_trna_counts_tsv(ran: dict[str, Path]) -> None:
    from giap.io.tsv import read_tsv

    header, rows = read_tsv(ran["trna_counts"])
    assert header == ["anticodon", "amino_acid", "n_copies"]
    assert rows == [["CAG", "Leu", "2"], ["GCC", "Gly", "1"]]


def test_xrefs_parquet(ran: dict[str, Path]) -> None:
    table = pq.read_table(ran["xrefs"])
    assert table.column_names == ["locus_tag", "db", "id"]
    assert sorted(tuple(r.values()) for r in table.to_pylist()) == [
        ("T_00010", "COG", "COG0090"),
        ("T_00010", "GO", "GO:0003735"),
        ("T_00010", "KEGG", "K02886"),
        ("T_00030", "EC", "5.6.2.2"),
    ]


def test_annotation_meta_light(ran: dict[str, Path]) -> None:
    meta = _meta(ran)
    assert set(meta) >= META_KEYS
    assert meta["bakta_version"] == "1.12.2"
    assert meta["db_version"] == "6.0"
    assert meta["db_type"] == "light"
    assert meta["annotation_source"] == "bakta_light"
    assert meta["anti_sd"] == "ACCTCCTTA"
    assert meta["flags"] == ["bakta_db_light"]
    assert meta["n_genes"] == 4
    assert meta["n_features"] == 6


def test_unknown_feature_types_reported(
    reg: Any, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    js, fa = _write_inputs(tmp_path / "in")
    with caplog.at_level(logging.WARNING):
        outs = _run(reg, js, fa, tmp_path / "work")
    assert _meta(outs)["skipped_types"] == {"future_type": 2}
    warned = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("future_type" in r.getMessage() for r in warned)


@pytest.mark.parametrize(
    ("db_type", "source", "flags"),
    [
        ("full", "bakta_full", []),
        ("light", "bakta_light", ["bakta_db_light"]),
        ("none", "user", []),
    ],
)
def test_annotation_source_from_db_type(
    reg: Any, tmp_path: Path, db_type: str, source: str, flags: list[str]
) -> None:
    js, fa = _write_inputs(tmp_path / "in", db_type=db_type)
    meta = _meta(_run(reg, js, fa, tmp_path / "work"))
    assert meta["annotation_source"] == source
    assert meta["flags"] == flags
    assert meta["db_type"] == ("unknown" if db_type == "none" else db_type)


def test_json_without_nt_uses_genome(reg: Any, tmp_path: Path, ran: dict[str, Path]) -> None:
    """Recorded JSON drops ``sequences[].nt`` (M2 lead decision 4); the genome input fills it."""
    js, fa = _write_inputs(tmp_path / "in2", with_nt=False)
    outs = _run(reg, js, fa, tmp_path / "work2")
    for key in ("genes", "features", "replicons", "rrna16s", "trna_counts", "xrefs"):
        assert outs[key].read_bytes() == ran[key].read_bytes(), key


def test_anti_sd_default_flag_without_16s(reg: Any, tmp_path: Path) -> None:
    data = _bakta_dict(db_type="full")
    data["features"] = [f for f in data["features"] if f.get("gene") != "rrs"]
    (tmp_path / "in").mkdir()
    js, fa = tmp_path / "in" / "bakta.json", tmp_path / "in" / "genome.norm.fa"
    js.write_text(json.dumps(data))
    fa.write_text(f">chr1\n{CHR}\n>p1\n{PL}\n")
    outs = _run(reg, js, fa, tmp_path / "work")
    meta = _meta(outs)
    assert meta["anti_sd"] == "ACCTCCTTA"
    assert meta["flags"] == ["anti_sd_default"]
    assert _fasta(outs["rrna16s"]) == []


def test_flags_sorted_with_light_db_and_default(reg: Any, tmp_path: Path) -> None:
    data = _bakta_dict(db_type="light")
    data["features"] = [f for f in data["features"] if f.get("gene") != "rrs"]
    (tmp_path / "in").mkdir()
    js, fa = tmp_path / "in" / "bakta.json", tmp_path / "in" / "genome.norm.fa"
    js.write_text(json.dumps(data))
    fa.write_text(f">chr1\n{CHR}\n>p1\n{PL}\n")
    meta = _meta(_run(reg, js, fa, tmp_path / "work"))
    assert meta["flags"] == ["anti_sd_default", "bakta_db_light"]


def test_configured_anti_sd_override(reg: Any, tmp_path: Path) -> None:
    js, fa = _write_inputs(tmp_path / "in")
    cfg = GiapConfig.model_validate({"seed": 1, "annotation": {"ostir": {"anti_sd": "TCCTCCACT"}}})
    meta = _meta(_run(reg, js, fa, tmp_path / "work", cfg))
    assert meta["anti_sd"] == "TCCTCCACT"
    assert "anti_sd_default" not in meta["flags"]


# ------------------------------------------------------------------ determinism (TE7) and resume


def test_byte_identical_rerun(reg: Any, tmp_path: Path, ran: dict[str, Path]) -> None:
    js, fa = _write_inputs(tmp_path / "in2")
    again = _run(reg, js, fa, tmp_path / "work2", GiapConfig(seed=99))
    for key in OUT_FILES:
        assert again[key].read_bytes() == ran[key].read_bytes(), key


def test_rerun_in_place_replaces_stale_outputs(reg: Any, tmp_path: Path) -> None:
    js, fa = _write_inputs(tmp_path / "in")
    work = tmp_path / "work"
    first = {k: p.read_bytes() for k, p in _run(reg, js, fa, work).items()}
    for name in OUT_FILES.values():
        (work / name).write_bytes(b"stale")
    second = _run(reg, js, fa, work)
    assert {k: p.read_bytes() for k, p in second.items()} == first


def test_failure_leaves_no_outputs(reg: Any, tmp_path: Path) -> None:
    """A malformed JSON fails the step (StepError) and no declared output survives (M2-02)."""
    js, fa = _write_inputs(tmp_path / "in")
    data = _bakta_dict()
    data["features"].append(_cds("T_00010", 2501, 2700, "+"))  # duplicate locus tag
    js.write_text(json.dumps(data))
    work = tmp_path / "work"
    work.mkdir()
    (work / "genes.parquet").write_bytes(b"stale")
    with pytest.raises(reg.StepError, match="parse_bakta"):
        _run(reg, js, fa, work)
    assert not any(work.iterdir())


def test_genome_mismatch_fails(reg: Any, tmp_path: Path) -> None:
    js, fa = _write_inputs(tmp_path / "in", with_nt=False)
    fa.write_text(f">chr1\n{CHR[:-10]}\n>p1\n{PL}\n")
    with pytest.raises(reg.StepError):
        _run(reg, js, fa, tmp_path / "work")


# ------------------------------------------------------------------ giapsim bakta_like.json


def _tsv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_simulator_m0_gene_count(reg: Any, tmp_path: Path) -> None:
    from giapsim import SimConfig, simulate

    res = simulate(SimConfig(genome_length=50_000), seed=1, outdir=tmp_path / "sim")
    outs = _run(reg, res.bakta_json, res.genome_fasta, tmp_path / "work")
    for key, path in outs.items():
        assert path.is_file(), key
    assert pq.read_metadata(outs["genes"]).num_rows == res.n_genes
    meta = _meta(outs)
    assert meta["annotation_source"] == "user"  # simulator JSON has no DB type
    assert meta["db_type"] == "unknown"
    assert meta["n_genes"] == res.n_genes
    assert meta["skipped_types"] == {}
    assert meta["flags"] == ["anti_sd_default"]  # M0 simulator has no rRNA


@pytest.fixture
def sim_rich(tmp_path: Path) -> Any:
    from giapsim import SimConfig, simulate

    cfg = SimConfig(
        genome_length=60_000,
        n_rrna_operons=2,
        n_trna=5,
        n_plasmids=1,
        plasmid_length=(5_000, 8_000),
        origin_crossing_gene=True,
        n_pseudogenes=2,
        n_sorfs=2,
        linear_contig_length=4_000,
    )
    return simulate(cfg, seed=3, outdir=tmp_path / "sim")


def test_simulator_rich(reg: Any, tmp_path: Path, sim_rich: Any) -> None:
    res = sim_rich
    outs = _run(reg, res.bakta_json, res.genome_fasta, tmp_path / "work")
    truth = _tsv_dicts(res.features_truth)
    gene_truth = {r["feature_id"]: r for r in truth if r["type"] in ("cds", "sorf")}
    genes = {r["locus_tag"]: r for r in _rows(outs["genes"])}
    assert set(genes) == set(gene_truth)
    for locus, t in gene_truth.items():
        g = genes[locus]
        # truth is 1-based closed (end unwrapped); Parquet is 0-based half-open
        assert (g["contig"], g["start"] + 1, g["end"], g["strand"]) == (
            t["contig"],
            int(t["start"]),
            int(t["end"]),
            t["strand"],
        ), locus
        want_type = "pseudogene" if t["pseudo_reason"] else t["type"]
        assert g["feature_type"] == want_type, locus
    assert sum(g["crosses_origin"] for g in genes.values()) == 1

    nc_truth = [r for r in truth if r["type"] not in ("cds", "sorf")]
    assert pq.read_metadata(outs["features"]).num_rows == len(nc_truth)

    reps = _rows(outs["replicons"])
    assert [(r["contig"], r["topology"], r["kind"]) for r in reps] == [
        ("chr1", "circular", "chromosome"),
        ("plasmid_1", "circular", "plasmid"),
        ("linear_1", "linear", "unknown"),
    ]

    s16 = _fasta(outs["rrna16s"])
    assert len(s16) == 2
    assert all(seq.endswith("GATCACCTCCTTA") and len(seq) == 1542 for _, seq in s16)
    meta = _meta(outs)
    assert meta["anti_sd"] == "ACCTCCTTA"
    assert meta["flags"] == []
    assert meta["annotation_source"] == "user"
    assert meta["n_genes"] == len(gene_truth)
    assert meta["n_features"] == len(nc_truth)

    from giap.io.tsv import read_tsv

    _, rows = read_tsv(outs["trna_counts"])
    assert sum(int(r[2] or 0) for r in rows) == 5
    assert [r[0] for r in rows] == sorted(r[0] for r in rows)


def test_simulator_rerun_byte_identical(reg: Any, tmp_path: Path, sim_rich: Any) -> None:
    a = _run(reg, sim_rich.bakta_json, sim_rich.genome_fasta, tmp_path / "w1")
    b = _run(reg, sim_rich.bakta_json, sim_rich.genome_fasta, tmp_path / "w2", GiapConfig(seed=7))
    for key in OUT_FILES:
        assert a[key].read_bytes() == b[key].read_bytes(), key
