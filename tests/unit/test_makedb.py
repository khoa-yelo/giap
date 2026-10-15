"""`giap makedb`: strain search database (issue M3-23; GIAP_SPEC.md §7.9.2, §5.1, §6.6).

Examples 1-9, the edge cases and the acceptance criteria of docs/specs/M3-23-makedb.md on

* the giapsim 50-kb genome (example 1 and a full six-frame comparison),
* a hand-built genome (``tests/helpers/makedb_genomes.py``: circular ``chr1`` of 3,300 nt and
  linear ``lin1`` of 660 nt) that holds examples 4-7 and 9 and the edge cases,
* a real Bakta 1.12.2 excerpt (``tests/data/recorded/makedb``, see its PROVENANCE.md).

Expected ORFs come from the brute-force oracle ``tests/helpers/makedb_oracle.py``.

Conventions these tests pin down where the spec is silent (resolved by the test author):

* A Bakta feature without ``product`` is written as ``>{locus_tag} hypothetical protein``
  (example 1: the giapsim JSON has no ``product`` and the expected header says so).
* Patched proteins are translated from the genome with table 11 and, like Bakta proteins, start
  with ``M`` when the first codon is an alternative start (GTG/TTG); six-frame ORFs and the
  upstream part of ``nterm_extension`` entries are translated codon by codon (GTG -> V).
* ``nterm_extension`` header/id_map interval: the upstream part plus the 30 junction codons
  (pinned by lead decision 9). Its ``locus_tag`` is the annotated gene.
* ``db_id_map.tsv`` (lead decision 8) follows the TSV conventions (§6.6): ``NA`` for missing
  values, strand ``+``/``-``; one row per FASTA record, in FASTA order; ``db_id`` is the first
  header token and unique.
* An ORF containing ``X`` is never written (whether it is dropped or split is not asserted).
* Patched records follow all CDS/sORF records and precede all ORFs.
* The CLI finds the contaminant library at ``DB_DIR/contaminants/contaminants.fasta`` (``giap
  setup`` layout, ``--db-dir``); when it is absent and ``--no-contaminants`` is not given, the
  exit code is 2 (§5.1: a database is missing).
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import sys
import time
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest

import giap
from giap.annotate.bakta import parse_bakta_json
from giap.config.model import InputConfig
from giap.exitcodes import ExitCode
from giap.io.protein_fasta import read_protein_fasta
from giap.model.coords import from_one_based_closed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

from makedb_checks import (
    assert_orfs_match_oracle,
    by_id,
    db_id,
    fasta_records,
    id_map,
    orf_ids,
    wrapped_exactly,
)
from makedb_genomes import (
    G1,
    G1_PRODUCT,
    G2,
    G3,
    G4,
    G5,
    P1,
    S1,
    ExampleGenome,
    build_example,
    write_bakta_json,
    write_fasta,
    write_patch_gff3,
)
from makedb_oracle import in_frame_overlap, strand_nt, translate

DATA = Path(__file__).resolve().parents[1] / "data"
CONTAMINANTS = DATA / "contaminants_mini.faa"
RECORDED = DATA / "recorded" / "makedb"
KINDS = {"cds", "sorf", "patch", "orf", "contaminant"}
ORF_CLASSES = {"intergenic", "cds_outframe", "antisense", "nterm_extension"}


# ---------------------------------------------------------------- helpers


def _api() -> Any:
    from giap.makedb import build

    return build


def _opts(annotation: Path, outdir: Path, **kw: Any) -> Any:
    return _api().MakedbOptions(
        annotation=annotation, genome=kw.pop("genome", None), outdir=outdir, **kw
    )


def build_db(annotation: Path, outdir: Path, **kw: Any) -> Path:
    _api().build_strain_db(_opts(annotation, outdir, **kw))
    return outdir


@pytest.fixture(scope="module")
def example() -> ExampleGenome:
    return build_example()


@pytest.fixture
def example_json(tmp_path: Path, example: ExampleGenome) -> Path:
    return write_bakta_json(tmp_path / "bakta.json", example.contigs, example.features)


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import giapsim
    from giapsim import SimConfig

    cfg = SimConfig(
        n_sorfs=3,
        n_pseudogenes=2,
        n_trna=2,
        n_rrna_operons=1,
        n_nterm_extensions=3,
        origin_crossing_gene=True,
        n_novel_orfs=(("intergenic", 2), ("cds_outframe", 1), ("antisense", 1)),
    )
    return giapsim.simulate(cfg, 11, tmp_path_factory.mktemp("sim"))


@pytest.fixture(scope="module")
def sim_db(sim: Any, tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("simdb") / "DB"
    return build_db(sim.bakta_json, out, six_frame_orfs=True)


# ---------------------------------------------------------------- example 1 (sim, defaults)


def test_example1_one_record_per_cds_and_sorf(sim: Any, tmp_path: Path) -> None:
    out = build_db(sim.bakta_json, tmp_path / "DB")
    gm = parse_bakta_json(sim.bakta_json).model
    recs = fasta_records(out / "strain_db.faa")
    coding = [g for g in gm.genes if g.feature_type in ("cds", "sorf")]
    assert {g.feature_type for g in coding} == {"cds", "sorf"}  # precondition
    assert len(recs) == len(coding)  # defaults: no ORFs, no patch, contaminants=None
    for (header, seq), g in zip(recs, coding, strict=True):
        assert header == f"{g.locus_tag} {g.product or 'hypothetical protein'}"
        assert seq == g.aa
    assert ("GSIM_00005 hypothetical protein", gm.by_locus("GSIM_00005").aa) in recs


def test_example1_pseudogenes_are_not_written(sim: Any, tmp_path: Path) -> None:
    out = build_db(sim.bakta_json, tmp_path / "DB")
    gm = parse_bakta_json(sim.bakta_json).model
    pseudo = {g.locus_tag for g in gm.genes if g.feature_type == "pseudogene"}
    assert len(pseudo) == 2  # precondition
    assert not pseudo & set(by_id(out))


def test_example1_id_map_rows_for_annotated_proteins(sim: Any, tmp_path: Path) -> None:
    out = build_db(sim.bakta_json, tmp_path / "DB")
    gm = parse_bakta_json(sim.bakta_json).model
    rows = id_map(out)
    coding = [g for g in gm.genes if g.feature_type in ("cds", "sorf")]
    assert [r["db_id"] for r in rows] == [g.locus_tag for g in coding]
    for r, g in zip(rows, coding, strict=True):
        assert r["kind"] == g.feature_type
        assert r["locus_tag"] == g.locus_tag
        assert r["contig"] == g.interval.contig
        assert (int(r["start"] or 0), int(r["end"] or 0)) == (g.interval.start + 1, g.interval.end)
        assert r["strand"] == g.interval.strand


# ---------------------------------------------------------------- contaminants (examples 2, 3)


def test_example2_no_contaminants(example_json: Path, tmp_path: Path) -> None:
    out = build_db(example_json, tmp_path / "DB", contaminants=None)
    assert not [i for i in by_id(out) if i.startswith("CONT_")]
    assert not [r for r in id_map(out) if r["kind"] == "contaminant"]


def test_example3_contaminants_as_cont_records_in_library_order(
    example_json: Path, tmp_path: Path
) -> None:
    out = build_db(example_json, tmp_path / "DB", contaminants=CONTAMINANTS)
    lib = read_protein_fasta(CONTAMINANTS)
    recs = fasta_records(out / "strain_db.faa")
    cont = [(db_id(h), s) for h, s in recs if h.startswith("CONT_")]
    assert cont == [(f"CONT_{r.id}", r.sequence) for r in lib]
    assert "CONT_P00761" in dict(cont)
    # contaminants come last
    n = len(cont)
    assert all(h.startswith("CONT_") for h, _ in recs[-n:])
    rows = [r for r in id_map(out) if r["kind"] == "contaminant"]
    assert [r["db_id"] for r in rows] == [c for c, _ in cont]
    for r in rows:
        assert all(r[c] is None for c in ("locus_tag", "class", "contig", "start", "end", "strand"))


def test_contaminant_ids_carry_a_recognised_contaminant_prefix(
    example_json: Path, tmp_path: Path
) -> None:
    out = build_db(example_json, tmp_path / "DB", contaminants=CONTAMINANTS)
    prefixes = tuple(InputConfig().contaminant_prefixes)
    cont = [i for i in by_id(out) if i.startswith("CONT_")]
    assert cont and all(i.startswith(prefixes) for i in cont)


# ---------------------------------------------------------------- six-frame ORFs (examples 4-7)


@pytest.fixture
def orf_db(example_json: Path, tmp_path: Path) -> Path:
    return build_db(example_json, tmp_path / "DB", six_frame_orfs=True)


def test_example4_intergenic_orf(orf_db: Path, example: ExampleGenome) -> None:
    recs = by_id(orf_db)
    assert recs[example.expect["intergenic"]] == example.expect["intergenic_aa"]


def test_example5_antisense_orf_on_the_minus_strand(orf_db: Path, example: ExampleGenome) -> None:
    recs = by_id(orf_db)
    assert recs[example.expect["antisense"]] == example.expect["antisense_aa"]


def test_example6_orf_of_19_aa_is_absent_by_default(orf_db: Path, example: ExampleGenome) -> None:
    ids = orf_ids(orf_db)
    assert not [i for i in ids if i.endswith("_chr1_1501_1557_p")]
    assert example.expect["orf20"] in ids  # 20 aa: kept (>= min_orf_aa)


def test_min_orf_aa_option_lowers_the_threshold(
    example_json: Path, tmp_path: Path, example: ExampleGenome
) -> None:
    out = build_db(example_json, tmp_path / "DB", six_frame_orfs=True, min_orf_aa=19)
    assert example.expect["orf19"] in orf_ids(out)
    out21 = build_db(example_json, tmp_path / "DB21", six_frame_orfs=True, min_orf_aa=21)
    assert example.expect["orf20"] not in orf_ids(out21)
    assert example.expect["intergenic"] in orf_ids(out21)


def test_example7_nterm_extension_plus_strand(orf_db: Path, example: ExampleGenome) -> None:
    recs = by_id(orf_db)
    hits = {i: s for i, s in recs.items() if i.startswith("GIAPORF_nterm_extension_chr1_565_")}
    assert set(hits) == {"GIAPORF_nterm_extension_chr1_565_690_p"}  # upstream + 30 codons
    assert list(hits.values()) == [example.expect["nterm_plus_aa"]]
    # 12 upstream residues, translated codon by codon (the GTG at 564 stays V), + 30 aa of G1
    seq = next(iter(hits.values()))
    assert seq[:12] == example.expect["nterm_plus_up"] and seq[0] == "V"
    assert seq[12:] == example.aa[G1][:30]


def test_example7_nterm_extension_minus_strand(orf_db: Path, example: ExampleGenome) -> None:
    recs = by_id(orf_db)
    hits = [
        s
        for i, s in recs.items()
        if i == "GIAPORF_nterm_extension_chr1_2211_2330_m"  # upstream + 30 codons
    ]
    assert hits == [example.expect["nterm_minus_aa"]]


def test_nterm_extension_id_map_row(orf_db: Path) -> None:
    rows = [r for r in id_map(orf_db) if r["class"] == "nterm_extension"]
    assert sorted(r["locus_tag"] or "" for r in rows) == sorted([G1, G3])
    assert {r["kind"] for r in rows} == {"orf"}


def test_upstream_part_shorter_than_7_aa_gives_no_nterm_entry(orf_db: Path) -> None:
    rows = [r for r in id_map(orf_db) if r["class"] == "nterm_extension"]
    assert not {G2, G4, G5, S1} & {r["locus_tag"] for r in rows}  # 1-aa upstream parts


def test_orf_in_frame_with_a_cds_or_sorf_is_dropped(
    orf_db: Path, example: ExampleGenome, example_json: Path
) -> None:
    gm = parse_bakta_json(example_json).model
    coding = [g for g in gm.genes if g.feature_type in ("cds", "sorf")]
    for row in id_map(orf_db):
        if row["kind"] != "orf" or row["class"] == "nterm_extension":
            continue
        rep = gm.replicons[row["contig"] or ""]
        iv = from_one_based_closed(
            row["contig"] or "",
            int(row["start"] or 0),
            int(row["end"] or 0),
            row["strand"],
            rep,  # type: ignore[arg-type]
        )
        same = [g for g in coding if g.interval.contig == iv.contig]
        assert not any(in_frame_overlap(iv.start, iv.end, iv.strand, g, rep) for g in same), row
    assert not [i for i in orf_ids(orf_db) if i.endswith(example.expect["sorf_frame"])]


def test_rna_overlap_and_pseudogene_orfs_are_omitted(orf_db: Path, example: ExampleGenome) -> None:
    ids = orf_ids(orf_db)
    assert not [i for i in ids if i.endswith(example.expect["rna"])]
    assert not [i for i in ids if i.endswith(example.expect["pseudo"])]
    assert not [i for i in ids if "rna_overlap" in i or "_pseudogene_" in i]


def test_orf_classes_are_the_four_spec_classes(orf_db: Path) -> None:
    for row in id_map(orf_db):
        if row["kind"] == "orf":
            assert row["class"] in ORF_CLASSES
            cls = row["class"]
            sfx = "p" if row["strand"] == "+" else "m"
            assert row["db_id"] == (
                f"GIAPORF_{cls}_{row['contig']}_{row['start']}_{row['end']}_{sfx}"
            )


def test_orfs_have_no_stop_and_no_x(orf_db: Path) -> None:
    recs = by_id(orf_db)
    for i, s in recs.items():
        assert "*" not in s, i
        if i.startswith("GIAPORF_"):
            assert "X" not in s, i


def test_orf_sequences_are_genome_translations_without_i_to_l(
    orf_db: Path, example: ExampleGenome
) -> None:
    genome = example.genome
    for i, s in by_id(orf_db).items():
        if not i.startswith("GIAPORF_") or "_nterm_extension_" in i:
            continue
        *_, contig, start1, end1, sfx = i.split("_")
        nt = strand_nt(genome[contig], int(start1) - 1, int(end1), "+" if sfx == "p" else "-")
        assert s == translate(nt), i
    assert any("I" in s for i, s in by_id(orf_db).items() if i.startswith("GIAPORF_"))


def test_origin_crossing_orf_header_is_unwrapped(orf_db: Path, example: ExampleGenome) -> None:
    recs = by_id(orf_db)
    assert recs[example.expect["origin"]] == example.expect["origin_aa"]
    for i in recs:
        if i.startswith("GIAPORF_"):
            *_, start1, end1, _sfx = i.split("_")
            assert int(start1) <= int(end1), f"wrapped coordinates in {i}"


def test_linear_contig_truncated_orfs(
    orf_db: Path, example: ExampleGenome, example_json: Path, tmp_path: Path
) -> None:
    recs = by_id(orf_db)
    assert recs[example.expect["lin_end"]] == example.expect["lin_end_aa"]  # 25 aa, no stop
    assert not [i for i in recs if i.endswith(example.expect["lin_start19"])]  # 19 aa
    out = build_db(example_json, tmp_path / "DB19", six_frame_orfs=True, min_orf_aa=19)
    assert [i for i in orf_ids(out) if i.endswith(example.expect["lin_start19"])]
    for row in id_map(orf_db):
        if row["contig"] == "lin1" and row["kind"] == "orf":
            assert int(row["end"] or 0) <= 660


def test_orfs_sorted_by_contig_start_end_strand(orf_db: Path) -> None:
    rows = [r for r in id_map(orf_db) if r["kind"] == "orf"]
    order = {"chr1": 0, "lin1": 1}
    keys = [
        (order[r["contig"] or ""], int(r["start"] or 0), int(r["end"] or 0), r["strand"] == "-")
        for r in rows
    ]
    assert keys == sorted(keys)
    assert {r["contig"] for r in rows} == {"chr1", "lin1"}


def test_six_frame_orfs_equal_the_oracle(
    orf_db: Path, example: ExampleGenome, example_json: Path
) -> None:
    assert_orfs_match_oracle(orf_db, example_json, example.genome, 20, skip=("lin1", 200, 386))


def test_six_frame_orfs_equal_the_oracle_min_10(
    example_json: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    out = build_db(example_json, tmp_path / "DB", six_frame_orfs=True, min_orf_aa=10)
    assert_orfs_match_oracle(out, example_json, example.genome, 10, skip=("lin1", 200, 386))


def test_six_frame_orfs_function_matches_the_written_orfs(
    orf_db: Path, example_json: Path, example: ExampleGenome
) -> None:
    from giap.mapping.sixframe import build_sixframe

    api = _api()
    ann = parse_bakta_json(example_json)
    index = build_sixframe(dict(ann.sequences), ann.model.replicons)
    orfs = api.six_frame_orfs(index, ann.model, dict(ann.sequences), 20)
    assert isinstance(orfs, list)
    assert len(orfs) == len(orf_ids(orf_db))
    assert len(api.six_frame_orfs(index, ann.model, dict(ann.sequences), 10)) > len(orfs)


# ---------------------------------------------------------------- Bakta proteins (edge cases)


def test_terminal_stop_is_stripped_and_duplicates_are_all_written(
    example_json: Path, tmp_path: Path, example: ExampleGenome
) -> None:
    out = build_db(example_json, tmp_path / "DB")
    recs = by_id(out)
    assert recs[G4] == example.aa[G4] and not recs[G4].endswith("*")
    assert recs[G5] == recs[G4]  # same protein under two IDs: both written
    rows = {r["db_id"]: r for r in id_map(out)}
    assert (rows[G4]["strand"], rows[G5]["strand"]) == ("+", "-")


def test_headers_carry_the_product(example_json: Path, tmp_path: Path) -> None:
    out = build_db(example_json, tmp_path / "DB")
    headers = dict((db_id(h), h) for h, _ in fasta_records(out / "strain_db.faa"))
    assert headers[G1] == f"{G1} {G1_PRODUCT}"
    assert headers[G2] == f"{G2} hypothetical protein"
    assert headers[S1] == f"{S1} sORF"
    assert P1 not in headers


def test_annotated_records_in_gene_model_order(example_json: Path, tmp_path: Path) -> None:
    out = build_db(example_json, tmp_path / "DB")
    gm = parse_bakta_json(example_json).model
    want = [g.locus_tag for g in gm.genes if g.feature_type in ("cds", "sorf")]
    assert list(by_id(out)) == want  # the JSON lists features in reverse locus order


def test_genome_fasta_option_gives_the_same_database(
    example_json: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    fa = write_fasta(tmp_path / "genome.fasta", example.contigs)
    a = build_db(example_json, tmp_path / "A", six_frame_orfs=True)
    b = build_db(example_json, tmp_path / "B", six_frame_orfs=True, genome=fa)
    assert (a / "strain_db.faa").read_bytes() == (b / "strain_db.faa").read_bytes()
    assert (a / "db_id_map.tsv").read_bytes() == (b / "db_id_map.tsv").read_bytes()


# ---------------------------------------------------------------- example 8: determinism


def test_example8_byte_identical_reruns(example_json: Path, tmp_path: Path) -> None:
    kw = {"six_frame_orfs": True, "contaminants": CONTAMINANTS}
    a = build_db(example_json, tmp_path / "A", **kw)
    b = build_db(example_json, tmp_path / "B", **kw)
    for name in ("strain_db.faa", "db_id_map.tsv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    # rerun into the same directory overwrites deterministically
    first = (a / "strain_db.faa").read_bytes()
    build_db(example_json, tmp_path / "A", **kw)
    assert (a / "strain_db.faa").read_bytes() == first


def test_fasta_has_no_timestamp(example_json: Path, tmp_path: Path) -> None:
    out = build_db(example_json, tmp_path / "DB", six_frame_orfs=True)
    text = (out / "strain_db.faa").read_text()
    assert not re.search(r"\d{4}-\d{2}-\d{2}|\d{2}:\d{2}:\d{2}", text)


# ---------------------------------------------------------------- example 9: patch


@pytest.fixture
def patch_db(example_json: Path, example: ExampleGenome, tmp_path: Path) -> Path:
    patch = write_patch_gff3(tmp_path / "annotation_patch.gff3", example)
    return build_db(
        example_json, tmp_path / "DB", patch=patch, six_frame_orfs=True, contaminants=CONTAMINANTS
    )


def test_example9_patch_novel_orf(patch_db: Path, example: ExampleGenome) -> None:
    recs = by_id(patch_db)
    assert recs["NOVEL_chr1_1201_p"] == example.expect["intergenic_aa"]  # ATG...; stop dropped
    row = next(r for r in id_map(patch_db) if r["db_id"] == "NOVEL_chr1_1201_p")
    assert row["kind"] == "patch"
    assert (row["contig"], row["start"], row["strand"]) == ("chr1", "1201", "+")
    assert row["end"] in {"1290", "1293"}


def test_patch_nterm_extension_alternative_start_is_methionine(
    patch_db: Path, example: ExampleGenome
) -> None:
    rows = [r for r in id_map(patch_db) if r["kind"] == "patch" and r["locus_tag"] == G1]
    assert len(rows) == 1
    seq = by_id(patch_db)[rows[0]["db_id"] or ""]
    assert seq == example.expect["patch_nterm_aa"]  # GTG start -> M, then the genome
    assert by_id(patch_db)[G1] == example.aa[G1]  # the Bakta protein is still written


def test_record_order_cds_patch_orfs_contaminants(patch_db: Path) -> None:
    kinds = [r["kind"] for r in id_map(patch_db)]
    rank = {"cds": 0, "sorf": 0, "patch": 1, "orf": 2, "contaminant": 3}
    assert [rank[k or ""] for k in kinds] == sorted(rank[k or ""] for k in kinds)
    assert {"cds", "sorf", "patch", "orf", "contaminant"} == set(kinds)


def test_db_ids_unique_and_id_map_follows_fasta(patch_db: Path) -> None:
    ids = [db_id(h) for h, _ in fasta_records(patch_db / "strain_db.faa")]
    assert len(ids) == len(set(ids))
    assert [r["db_id"] for r in id_map(patch_db)] == ids
    assert {r["kind"] for r in id_map(patch_db)} <= KINDS


def test_fasta_wrapped_at_60(patch_db: Path) -> None:
    wrapped_exactly(patch_db / "strain_db.faa")
    for _, seq in fasta_records(patch_db / "strain_db.faa"):
        assert seq and seq.isalpha() and seq.isupper()


def test_no_decoys(patch_db: Path) -> None:
    prefixes = tuple(InputConfig().decoy_prefixes)
    assert not [i for i in by_id(patch_db) if i.startswith(prefixes) or "decoy" in i.lower()]
    assert not [p for p in patch_db.iterdir() if "decoy" in p.name.lower()]


# ---------------------------------------------------------------- manifest


def _walk(x: Any) -> Iterator[Any]:
    yield x
    if isinstance(x, dict):
        for v in x.values():
            yield from _walk(v)
    elif isinstance(x, list):
        for v in x:
            yield from _walk(v)


def test_manifest_version_inputs_options_counts(patch_db: Path, example_json: Path) -> None:
    m = json.loads((patch_db / "makedb_manifest.json").read_text())
    assert isinstance(m, dict)
    values = list(_walk(m))
    assert giap.__version__ in values
    sha = hashlib.sha256(example_json.read_bytes()).hexdigest()
    assert sha in values
    assert hashlib.sha256(CONTAMINANTS.read_bytes()).hexdigest() in values
    counts: dict[str, int] = {}
    for r in id_map(patch_db):
        counts[r["kind"] or ""] = counts.get(r["kind"] or "", 0) + 1
    assert any(
        isinstance(d, dict) and all(d.get(k) == n for k, n in counts.items()) for d in values
    ), f"no per-kind counts {counts} in the manifest"
    assert 20 in values and True in values  # options: min_orf_aa, six_frame_orfs


def test_manifest_written_without_six_frame_orfs(example_json: Path, tmp_path: Path) -> None:
    out = build_db(example_json, tmp_path / "DB")
    assert sorted(p.name for p in out.iterdir()) == [
        "db_id_map.tsv",
        "makedb_manifest.json",
        "strain_db.faa",
    ]


# ---------------------------------------------------------------- real Bakta excerpt


def test_recorded_bakta_proteins_equal_baktas_own_faa(tmp_path: Path) -> None:
    out = build_db(RECORDED / "sim100k_plasmid_1.bakta.json", tmp_path / "DB")
    want = [
        (r.header, r.sequence) for r in read_protein_fasta(RECORDED / "sim100k_plasmid_1.bakta.faa")
    ]
    assert fasta_records(out / "strain_db.faa") == want
    # KABHIH_0098 starts with TTG: Bakta's protein (and the database) starts with M
    assert dict(want)["KABHIH_0098 hypothetical protein"][0] == "M"


def test_recorded_bakta_six_frame_orfs_equal_the_oracle(tmp_path: Path) -> None:
    ann_path = RECORDED / "sim100k_plasmid_1.bakta.json"
    out = build_db(ann_path, tmp_path / "DB", six_frame_orfs=True)
    ann = parse_bakta_json(ann_path)
    assert_orfs_match_oracle(out, ann_path, dict(ann.sequences), 20)


# ---------------------------------------------------------------- sim: whole-genome oracle


def test_sim_six_frame_orfs_equal_the_oracle(sim: Any, sim_db: Path) -> None:
    ann = parse_bakta_json(sim.bakta_json)
    assert_orfs_match_oracle(sim_db, sim.bakta_json, dict(ann.sequences), 20)


def test_sim_planted_nterm_extensions_are_written(sim: Any, sim_db: Path) -> None:
    rows = [r for r in id_map(sim_db) if r["class"] == "nterm_extension"]
    assert len(rows) >= 3  # three planted extensions with >= 6 stop-free codons upstream + start


def test_sim_classes_present(sim_db: Path) -> None:
    classes = {r["class"] for r in id_map(sim_db) if r["kind"] == "orf"}
    assert {"intergenic", "cds_outframe", "antisense", "nterm_extension"} <= classes


# ---------------------------------------------------------------- end-to-end readability


def test_database_is_readable_by_the_m1_and_m3_readers(patch_db: Path, example_json: Path) -> None:
    recs = read_protein_fasta(patch_db / "strain_db.faa")
    rows = id_map(patch_db)
    assert [r.id for r in recs] == [r["db_id"] for r in rows]
    for r in recs:
        assert r.internal_stops == 0
        assert not dict(r.other)
        assert r.id == r.header.split(" ", 1)[0]
    ann = parse_bakta_json(example_json)
    for row in rows:
        if row["kind"] == "contaminant":
            continue
        rep = ann.model.replicons[row["contig"] or ""]
        iv = from_one_based_closed(
            row["contig"] or "",
            int(row["start"] or 0),
            int(row["end"] or 0),
            row["strand"],
            rep,  # type: ignore[arg-type]
        )
        assert iv.strand in ("+", "-")


def test_database_peptides_map_with_the_exact_matcher(orf_db: Path, example: ExampleGenome) -> None:
    """Every ORF's tryptic peptides map exactly (M3-04) inside its id_map locus."""
    from giap.genetable.covariates import digest
    from giap.mapping.decoys import build_decoys
    from giap.mapping.exact import match_exact
    from giap.mapping.sixframe import build_sixframe
    from giap.model.coords import Replicon
    from giap.rng import for_step

    reps = {
        "chr1": Replicon("chr1", len(example.genome["chr1"]), "circular"),
        "lin1": Replicon("lin1", len(example.genome["lin1"]), "linear"),
    }
    idx = build_sixframe(example.genome, reps)
    decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", 0))}
    rows = {r["db_id"]: r for r in id_map(orf_db)}
    checked = 0
    for rid, seq in by_id(orf_db).items():
        row = rows[rid]
        if row["kind"] != "orf" or row["class"] == "nterm_extension":
            continue
        peps = [p for p in digest(seq, "trypsin") if 7 <= len(p) <= 40]
        if not peps:
            continue
        hits = match_exact(peps, idx, decoys, {p: "reverse" for p in peps})
        hits = hits[hits["search"] == "target"]
        rep = reps[row["contig"] or ""]
        s0, e0 = int(row["start"] or 0) - 1, int(row["end"] or 0)
        shifts = (0, rep.length, -rep.length) if rep.topology == "circular" else (0,)
        for p in peps:
            h = hits[(hits["peptide"] == p) & (hits["strand"] == row["strand"])]
            assert any(
                s0 <= hs + d and he + d <= e0
                for hs, he in zip(h["start"], h["end"], strict=True)
                for d in shifts
            ), (rid, p)
            checked += 1
    assert checked > 0


# ---------------------------------------------------------------- CLI


def _main(argv: list[str]) -> int:
    from giap.cli import main

    with pytest.raises(SystemExit) as info:
        main(argv)
    code = info.value.code
    return 0 if code is None else int(code)


def test_cli_makedb_exit_0(example_json: Path, tmp_path: Path) -> None:
    out = tmp_path / "DB"
    argv = ["makedb", "--annotation", str(example_json), "--outdir", str(out), "--no-contaminants"]
    assert _main(argv) == ExitCode.SUCCESS
    assert {p.name for p in out.iterdir()} >= {
        "strain_db.faa",
        "db_id_map.tsv",
        "makedb_manifest.json",
    }
    assert G1 in by_id(out) and not orf_ids(out)


def test_cli_matches_the_python_api(
    example_json: Path, example: ExampleGenome, tmp_path: Path
) -> None:
    fa = write_fasta(tmp_path / "g.fasta", example.contigs)
    patch = write_patch_gff3(tmp_path / "annotation_patch.gff3", example)
    argv = [
        "makedb",
        "--annotation",
        str(example_json),
        "--genome",
        str(fa),
        "--outdir",
        str(tmp_path / "CLI"),
        "--six-frame-orfs",
        "--min-orf-aa",
        "19",
        "--patch",
        str(patch),
        "--no-contaminants",
    ]
    assert _main(argv) == ExitCode.SUCCESS
    api = build_db(
        example_json, tmp_path / "API", genome=fa, six_frame_orfs=True, min_orf_aa=19, patch=patch
    )
    for name in ("strain_db.faa", "db_id_map.tsv"):
        assert (tmp_path / "CLI" / name).read_bytes() == (api / name).read_bytes(), name
    assert example.expect["orf19"] in orf_ids(tmp_path / "CLI")


def test_cli_without_annotation_exits_2(
    example: ExampleGenome, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    fa = write_fasta(tmp_path / "g.fasta", example.contigs)
    argv = ["makedb", "--genome", str(fa), "--outdir", str(tmp_path / "DB"), "--no-contaminants"]
    assert _main(argv) == ExitCode.ENVIRONMENT
    assert "Bakta wrapper required (M2)" in capsys.readouterr().err
    assert not (tmp_path / "DB" / "strain_db.faa").exists()


def test_cli_genome_is_optional(example_json: Path, tmp_path: Path) -> None:
    argv = [
        "makedb",
        "--annotation",
        str(example_json),
        "--outdir",
        str(tmp_path / "DB"),
        "--no-contaminants",
        "--six-frame-orfs",
    ]
    assert _main(argv) == ExitCode.SUCCESS
    assert orf_ids(tmp_path / "DB")


def test_cli_contaminants_from_db_dir(example_json: Path, tmp_path: Path) -> None:
    db = tmp_path / "giapdb"
    (db / "contaminants").mkdir(parents=True)
    shutil.copy(CONTAMINANTS, db / "contaminants" / "contaminants.fasta")
    argv = [
        "makedb",
        "--annotation",
        str(example_json),
        "--outdir",
        str(tmp_path / "DB"),
        "--db-dir",
        str(db),
    ]
    assert _main(argv) == ExitCode.SUCCESS
    assert "CONT_P00761" in by_id(tmp_path / "DB")


def test_cli_missing_contaminant_library_exits_2(
    example_json: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "empty_db"
    empty.mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(empty))
    argv = [
        "makedb",
        "--annotation",
        str(example_json),
        "--outdir",
        str(tmp_path / "DB"),
        "--db-dir",
        str(empty),
    ]
    assert _main(argv) == ExitCode.ENVIRONMENT


def test_cli_rejects_min_orf_aa_below_1(example_json: Path, tmp_path: Path) -> None:
    argv = [
        "makedb",
        "--annotation",
        str(example_json),
        "--outdir",
        str(tmp_path / "DB"),
        "--no-contaminants",
        "--min-orf-aa",
        "0",
    ]
    assert _main(argv) == ExitCode.INVALID_INPUT


# ---------------------------------------------------------------- benchmark (acceptance)


@pytest.mark.slow
def test_six_frame_orfs_on_5_mb_within_30_s(tmp_path: Path) -> None:
    import numpy as np
    from makedb_genomes import Contig

    rng = np.random.default_rng(23)
    seq = "".join(np.array(list("ACGT"))[rng.integers(0, 4, 5_000_000)])
    contig = Contig("chr1", len(seq), "circular", seq=list(seq))
    ann = write_bakta_json(tmp_path / "big.json", [contig], [])
    t0 = time.perf_counter()
    build_db(ann, tmp_path / "DB", six_frame_orfs=True)
    assert time.perf_counter() - t0 <= 30.0
