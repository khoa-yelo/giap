"""``annotation_patch.gff3`` -> ``giap makedb --patch`` round trip (issue M7-07 example 7 and
rule 8; GIAP_SPEC.md §7.6.6 "patch -> makedb round trip"; M3-23's collision rule ``<id>_patch``).

Two arms: the hand-built genome of examples 1 and 6 (both strands, and an origin-crossing ORF), and
a giapsim genome with planted novel ORFs and N-terminal extensions whose ``Orf`` records are built
from the simulator truth. For every new ORF the database protein equals its ``|best_ostir`` (or
``|longest``) record of ``novel_orfs.faa``; every N-terminal extension gives ``<locus>_patch``,
the extended protein starting with ``M``; ``db_id_map.tsv`` rows have kind ``patch``.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import pytest
from helpers.novel_outputs_fixture import (
    Case,
    ex1_ex6,
    fasta_records,
    mirror_case,
    orfs_of,
    protein_from,
    scored,
    table,
    wrap_case,
    write_all,
)

from giap.annotate.bakta import parse_bakta_json
from giap.io.protein_fasta import read_protein_fasta
from giap.io.tsv import read_tsv
from giap.model.coords import Replicon
from giap.novel.candidates import Orf


def _bakta_json(path: Path, case: Case) -> Path:
    seqs = [
        {"id": r.id, "description": "", "nt": case.genome[r.id], "length": r.length,
         "complete": True, "type": r.kind, "topology": r.topology}
        for r in case.replicons.values()
    ]  # fmt: skip
    feats = []
    for g in case.genes:
        iv = g.interval
        seq = case.genome[iv.contig]
        feats.append({
            "type": "cds", "sequence": iv.contig, "start": iv.start + 1, "stop": iv.end,
            "strand": iv.strand, "locus": g.locus_tag, "id": g.locus_tag,
            "aa": protein_from(seq, iv.start, iv.end, iv.strand), "product": "gene product",
        })  # fmt: skip
    data = {"version": {"bakta": "1.12.2", "db": {"version": "6.0", "type": "full"}},
            "sequences": seqs, "features": feats}  # fmt: skip
    path.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    return path


def _makedb(annotation: Path, patch: Path, outdir: Path) -> tuple[dict[str, str], list[dict]]:
    from giap.makedb.build import MakedbOptions, build_strain_db

    build_strain_db(MakedbOptions(annotation=annotation, genome=None, outdir=outdir, patch=patch))
    db = {r.id: r.sequence for r in read_protein_fasta(outdir / "strain_db.faa")}
    cols, rows = read_tsv(outdir / "db_id_map.tsv")
    return db, [dict(zip(cols, r, strict=True)) for r in rows]


def _faa(paths: dict[str, Path]) -> dict[str, str]:
    out = {}
    for header, seq in fasta_records(paths["novel_orfs.faa"]):
        oid, _, variants = header.partition("|")
        if "best_ostir" in variants.split(",") or oid not in out:
            out[oid] = seq
    return out


# ---------------------------------------------------------------- example 7 (hand-built)


@pytest.mark.parametrize("minus", [False, True])
def test_example7_makedb_patch_round_trip(tmp_path: Path, minus: bool) -> None:
    case = mirror_case(ex1_ex6()) if minus else ex1_ex6()
    orfs = orfs_of(case)
    by_cls = {o.cls: o for o in orfs}
    new, ext = by_cls["intergenic"], by_cls["nterm_extension"]
    gtg = 937 if minus else 60
    sc = scored(orfs, best={new.orf_id: gtg}, called={ext.orf_id: True})
    paths = write_all(tmp_path / "tables", case, table(case, orfs, sc))
    db, rows = _makedb(_bakta_json(tmp_path / "bakta.json", case), paths["annotation_patch.gff3"],
                       tmp_path / "DB")  # fmt: skip
    faa = dict(fasta_records(paths["novel_orfs.faa"]))
    assert db[new.orf_id] == faa[f"{new.orf_id}|best_ostir"]
    assert len(db[new.orf_id]) == 20
    seq = case.genome["c1"]
    want_g = protein_from(seq, 100, 751, "-") if minus else protein_from(seq, 249, 900, "+")
    assert db["G_patch"] == want_g
    assert want_g[0] == "M" and want_g.endswith(db["G"][1:])  # the TTG initiator read as M
    assert db["G_patch"] == _faa(paths)[ext.orf_id]  # no OSTIR score: the |longest record
    patch = {r["db_id"]: r for r in rows if r["kind"] == "patch"}
    assert set(patch) == {new.orf_id, "G_patch"}
    assert patch["G_patch"]["locus_tag"] == "G"
    if minus:
        assert (patch[new.orf_id]["start"], patch[new.orf_id]["end"]) == ("878", "940")
        assert (patch["G_patch"]["start"], patch["G_patch"]["end"]) == ("101", "751")
    else:
        assert (patch[new.orf_id]["start"], patch[new.orf_id]["end"]) == ("61", "123")
        assert (patch["G_patch"]["start"], patch["G_patch"]["end"]) == ("250", "900")


def test_origin_crossing_orf_round_trip(tmp_path: Path) -> None:
    case = wrap_case()
    orfs = orfs_of(case)
    (orf,) = orfs
    paths = write_all(tmp_path / "tables", case, table(case, orfs, scored(orfs)))
    db, rows = _makedb(_bakta_json(tmp_path / "bakta.json", case), paths["annotation_patch.gff3"],
                       tmp_path / "DB")  # fmt: skip
    assert db[orf.orf_id] == _faa(paths)[orf.orf_id] == protein_from(
        case.genome["c1"], 972, 1050, "+"
    )  # fmt: skip
    (row,) = [r for r in rows if r["kind"] == "patch"]
    assert (row["start"], row["end"]) == ("973", "1053")


def test_patch_without_called_orfs_adds_nothing(tmp_path: Path) -> None:
    case = ex1_ex6()
    orfs = orfs_of(case)
    sc = scored(orfs, called={o.orf_id: False for o in orfs})
    paths = write_all(tmp_path / "tables", case, table(case, orfs, sc))
    db, rows = _makedb(_bakta_json(tmp_path / "bakta.json", case), paths["annotation_patch.gff3"],
                       tmp_path / "DB")  # fmt: skip
    assert [r["kind"] for r in rows] == ["cds"]
    assert set(db) == {"G"}


# ---------------------------------------------------------------- simulated genome


def _sim(seed: int, out: Path) -> Any:
    import giapsim
    from giapsim import SimConfig

    cfg = SimConfig(
        n_nterm_extensions=3,
        n_novel_orfs=(("intergenic", 2), ("cds_outframe", 1), ("antisense", 1)),
    )
    return giapsim.simulate(cfg, seed, out)


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _orf(cls: str, contig: str, lo: int, hi: int, strand: str, start: int, stop: int,
         gene: str | None) -> Orf:  # fmt: skip
    """An ``Orf`` for codon-aligned ``[lo, hi)`` (start codon at ``start``, stop at ``stop``)."""
    return Orf(
        orf_id=f"NOVEL_{contig}_{lo + 1}_{'p' if strand == '+' else 'm'}",
        cls=cls,  # type: ignore[arg-type]
        contig=contig,
        strand=strand,  # type: ignore[arg-type]
        frame=1 + lo % 3 if strand == "+" else -1,
        stop=stop,
        starts=(start,),
        start_longest=start,
        end_5p_open=False,
        end_3p_open=False,
        peptides=("PEPTLDEA", "PEPTLDEK"),
        best_q=0.001,
        overlapping_gene=gene,
    )


@pytest.mark.parametrize("seed", [11, 12])
def test_simulated_planted_orfs_and_extensions_round_trip(tmp_path: Path, seed: int) -> None:
    sim = _sim(seed, tmp_path / "sim")
    ann = parse_bakta_json(sim.bakta_json)
    genome = dict(ann.sequences)
    reps: dict[str, Replicon] = dict(ann.model.replicons)
    orfs, want = [], {}
    for r in _tsv(sim.novel_orfs_truth):
        lo, hi = int(r["start"]) - 1, int(r["end"])
        plus = r["strand"] == "+"
        o = _orf(r["class"], r["contig"], lo, hi, r["strand"], lo if plus else hi - 3,
                 hi - 3 if plus else lo, None)  # fmt: skip
        orfs.append(o)
        want[o.orf_id] = r["aa"]
    n_ext = 0
    for r in _tsv(sim.features_truth):
        if r["type"] != "cds" or not r["true_start"]:
            continue
        s0, e0, ts = int(r["start"]) - 1, int(r["end"]), int(r["true_start"])
        if r["strand"] == "+":
            lo, hi, start, stop = ts - 1, e0, ts - 1, e0 - 3
        else:
            lo, hi, start, stop = s0, ts, ts - 3, s0
        orfs.append(_orf("nterm_extension", r["contig"], lo, hi, r["strand"], start, stop,
                         r["feature_id"]))  # fmt: skip
        want[f"{r['feature_id']}_patch"] = protein_from(genome[r["contig"]], lo, hi, r["strand"])
        n_ext += 1
    assert n_ext == 3 and len(want) == 7
    case = Case(genome, reps, [], genes=list(ann.model.genes))
    paths = write_all(tmp_path / "tables", case, table(case, orfs, scored(orfs)))
    db, rows = _makedb(sim.bakta_json, paths["annotation_patch.gff3"], tmp_path / "DB")
    faa = _faa(paths)
    for key, aa in want.items():
        assert db[key] == aa, key
        assert aa[0] == "M" and "*" not in aa
        if not key.endswith("_patch"):
            assert faa[key] == aa, key
        else:
            tag = key.removesuffix("_patch")
            assert aa.endswith(db[tag][1:]) and len(aa) > len(db[tag])
    kinds = {r["db_id"]: r["kind"] for r in rows}
    assert {k for k, v in kinds.items() if v == "patch"} == set(want)
