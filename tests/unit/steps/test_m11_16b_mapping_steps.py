"""M11-16b ``steps/mapping`` wiring (docs/specs/M11-16-translation-tables-4-25.md item 5.3, AC 8,
examples 3, 5 and 6 through ``sixframe_index -> map_exact -> map_tolerant -> placement -> loci``).

Every step takes ``cfg.annotation.translation_table``: ``sixframe_index`` builds and records the
table, the readers refuse an index built for another table ("rebuild the six-frame index"), the
loci walk and the stop sets follow the table, and the ``mapping.stop_sites`` event ``data`` gains
``"translation_table": t``. As in ``test_loci_readthrough_step``, the toy has no decoy power, so the
test accepts every placed peptide before ``loci``.
"""

from __future__ import annotations

import json
import logging
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.readthrough_toys import gene, model, rc, write_annotation, write_genome

from gippa.config import resolve_config
from gippa.mapping.sixframe import load_stop_cells
from gippa.model import Replicon
from gippa.steps import load_step_modules
from gippa.steps import registry as reg

EX3 = "ATGAAAGCCTGAGCCGAAGCCTGGAAATAA"  # MK A (TGA) A E A W K *
EVENT = "mapping.stop_sites"
PLACEMENT_OUT = {
    "placement": "placement.parquet",
    "mapping_fdr": "mapping_fdr.parquet",
    "stratum0_classes": "s0.parquet",
    "placement_costs": "costs.parquet",
}


@pytest.fixture(autouse=True)
def _steps() -> None:
    load_step_modules()


def _cfg(table: int) -> Any:
    return resolve_config(overrides=["seed=7", f"annotation.translation_table={table}"])


def _peptides(path: Path, peps: list[str]) -> Path:
    n = len(peps)
    pd.DataFrame({
        "peptide": peps, "source": ["database"] * n, "q_value": [0.001] * n,
        "pep": [0.01] * n, "score": [np.nan] * n, "is_decoy": [False] * n,
        "is_contaminant": [False] * n, "protein_ids": [["REF"]] * n,
        "modified_peptides": [[p] for p in peps], "n_rows": [1] * n, "n_runs": [1] * n,
        "intensity": [1.0e6] * n, "search_db": ["other"] * n,
    }).to_parquet(path, index=False)  # fmt: skip
    return path


class Chain:
    """The mapping steps on one toy world under one table."""

    def __init__(self, tmp: Path, genome: dict[str, str], gm: Any, peps: list[str], table: int):
        self.tmp, self.cfg = tmp, _cfg(table)
        ann = write_annotation(gm, tmp / "ann")
        self.paths: dict[str, Path] = {
            "genome": write_genome(genome, tmp / "genome.fa"),
            "replicons": ann["replicons"], "genes": ann["genes"], "features": ann["features"],
            "peptides": _peptides(tmp / "peptides.parquet", peps),
        }  # fmt: skip

    def run(self, name: str, outs: dict[str, str], cfg: Any = None) -> None:
        spec = reg.get(name)
        given = {k: self.paths[k] for k in spec.inputs}
        given |= {k: self.paths[k] for k in spec.optional_inputs if k in self.paths}
        out = {k: self.tmp / v for k, v in outs.items()}
        reg.run_step(name, cfg or self.cfg, given, out)
        self.paths.update(out)

    def index(self) -> Chain:
        self.run("sixframe_index", {"sixframe": "sixframe.npz"})
        return self

    def exact(self) -> pd.DataFrame:
        self.run("map_exact", {"hits_exact": "hits_exact.parquet"})
        hx = pd.read_parquet(self.paths["hits_exact"])
        return hx[hx["search"] == "target"]

    def loci(self) -> pd.DataFrame:
        if "hits_exact" not in self.paths:
            self.exact()
        self.run("map_tolerant", {"hits_tolerant": "hits_tolerant.parquet"})
        self.run("placement", PLACEMENT_OUT)
        per = pd.read_parquet(self.paths["placement"])
        per["accepted"] = True
        per["placement_q"] = 0.001
        per.to_parquet(self.paths["placement"], index=False)
        self.run("loci", {"mapping": "mapping"})
        return pd.read_csv(
            self.paths["mapping"] / "peptide_hits.tsv", sep="\t", dtype=str, keep_default_na=False
        )


def _ex3(circular: bool = False, minus: bool = False) -> tuple[dict[str, str], Any]:
    if circular:  # the gene across the origin, starting at 20
        nt = EX3[10:] + EX3[:10]
        reps = {"g": Replicon("g", 30, "circular")}
        return {"g": nt}, model(reps, [gene("G", "g", 20, 50, "+")])
    reps = {"g": Replicon("g", 30, "linear")}
    if minus:
        return {"g": rc(EX3)}, model(reps, [gene("G", "g", 0, 30, "-")])
    return {"g": EX3}, model(reps, [gene("G", "g", 0, 30, "+")])


def _rows(t: pd.DataFrame, pep: str) -> list[tuple[str, int, int, str]]:
    r = t[t["peptide"] == pep]
    return sorted(
        zip(r["contig"], r["start"].astype(int), r["end"].astype(int), r["strand"], strict=True)
    )


# ---------------------------------------------------------------- example 3 through the steps


@pytest.mark.parametrize(
    ("table", "want_w", "want_g"),
    [
        (4, [("g", 6, 27, "+")], []),
        (11, [], []),
        (25, [], [("g", 6, 27, "+")]),
    ],
)
def test_example_3_map_exact_follows_the_table(
    tmp_path: Path, table: int, want_w: list[Any], want_g: list[Any]
) -> None:
    genome, gm = _ex3()
    ch = Chain(tmp_path, genome, gm, ["AWAEAWK", "AGAEAWK"], table).index()
    with np.load(ch.paths["sixframe"]) as npz:
        assert json.loads(str(npz["meta"]))["translation_table"] == table
    t = ch.exact()
    assert _rows(t, "AWAEAWK") == want_w
    assert _rows(t, "AGAEAWK") == want_g


def test_example_3_minus_strand_under_4(tmp_path: Path) -> None:
    genome, gm = _ex3(minus=True)
    t = Chain(tmp_path, genome, gm, ["AWAEAWK"], 4).index().exact()
    assert _rows(t, "AWAEAWK") == [("g", 3, 24, "-")]


def test_example_3_circular_origin_spanning_under_4(tmp_path: Path) -> None:
    genome, gm = _ex3(circular=True)
    ch = Chain(tmp_path, genome, gm, ["AWAEAWK"], 4).index()
    assert _rows(ch.exact(), "AWAEAWK") == [("g", 26, 47, "+")]
    hits = ch.loci()
    (row,) = hits[hits["peptide"] == "AWAEAWK"].to_dict("records")
    assert "origin_spanning" in row["flags"].split(";")
    assert (row["locus_class"], row["gene"]) == ("cds_inframe", "G")


def test_example_3_loci_cds_inframe_under_4(tmp_path: Path) -> None:
    genome, gm = _ex3()
    hits = Chain(tmp_path, genome, gm, ["AWAEAWK"], 4).index().loci()
    (row,) = hits[hits["peptide"] == "AWAEAWK"].to_dict("records")
    assert (row["locus_class"], row["gene"]) == ("cds_inframe", "G")
    assert row["flags"] == "NA"  # empty flags are written as NA (test_mapping_outputs)


# ---------------------------------------------------------------- the index and the event


@pytest.mark.parametrize("table", (4, 11, 25))
def test_sixframe_index_records_the_table_and_the_event(
    tmp_path: Path, table: int, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.INFO)
    genome, gm = _ex3()
    ch = Chain(tmp_path, genome, gm, ["AWAEAWK"], table).index()
    with np.load(ch.paths["sixframe"]) as npz:
        meta = json.loads(str(npz["meta"]))
    assert (meta["format"], meta["translation_table"]) == (4, table)  # 4: M11-20
    (ev,) = [r for r in caplog.records if getattr(r, "event", None) == EVENT]
    data = ev.data  # type: ignore[attr-defined]
    assert data["translation_table"] == table
    assert (data["terminal"], data["not_a_stop"]) == (1, 0)


def _sec_world(codon: str = "TGA", letter: str = "U") -> tuple[dict[str, str], Any]:
    """Linear ``p`` 90 nt: G1 cds ``+`` ``[3, 51)`` = M K D E F H (codon) P Q R S T V Y W, TAA;
    the protein has ``letter`` at residue 6 (codon ``[21, 24)``)."""
    code = {"M": "ATG", "K": "AAA", "D": "GAT", "E": "GAA", "F": "TTT", "H": "CAT",
            "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT", "T": "ACC", "V": "GTT",
            "Y": "TAT", "W": "TGG"}  # fmt: skip
    body = "MKDEFH?PQRSTVYW"
    g = "".join(codon if a == "?" else code[a] for a in body) + "TAA"
    nt = "GCT" + g + "GCT" * 13
    assert len(nt) == 90 and nt[21:24] == codon
    reps = {"p": Replicon("p", 90, "linear")}
    gm = model(reps, [replace(gene("G1", "p", 3, 51, "+"), aa=body.replace("?", letter))])
    return {"p": nt}, gm


def test_example_5_under_4_through_the_steps(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """A ``U`` over TGA under 4: no recoded site, ``codon_mismatch`` 1; the peptide reads W at the
    TGA, maps exactly and is ``cds_inframe`` with no ``sec_site`` / ``readthrough``."""
    caplog.set_level(logging.INFO)
    genome, gm = _sec_world()
    pep = "KDEFHWPQRSTV"
    ch = Chain(tmp_path, genome, gm, [pep], 4).index()
    (ev,) = [r for r in caplog.records if getattr(r, "event", None) == EVENT]
    assert ev.data == {  # type: ignore[attr-defined]
        "terminal": 1, "recoded": 0, "not_a_stop": 0, "cells_target": 1,
        "cells_decoy": {"reverse": 1},
        "recoded_skipped": {"length_mismatch": 0, "codon_mismatch": 1},
        "translation_table": 4,
    }  # fmt: skip
    assert _rows(ch.exact(), pep) == [("p", 6, 42, "+")]
    (row,) = ch.loci().query("peptide == @pep").to_dict("records")
    assert (row["locus_class"], row["gene"]) == ("cds_inframe", "G1")
    assert "sec_site" not in row["flags"] and "readthrough" not in row["flags"]


def test_example_5_under_11_through_the_steps_keeps_sec_site(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.INFO)
    genome, gm = _sec_world()
    pep = "KDEFHWPQRSTV"
    ch = Chain(tmp_path, genome, gm, [pep], 11).index()
    (ev,) = [r for r in caplog.records if getattr(r, "event", None) == EVENT]
    data = ev.data  # type: ignore[attr-defined]
    assert (data["recoded"], data["translation_table"]) == (1, 11)
    assert data["recoded_skipped"] == {"length_mismatch": 0, "codon_mismatch": 0}
    assert _rows(ch.exact(), pep) == []
    (row,) = ch.loci().query("peptide == @pep").to_dict("records")
    assert "sec_site" in row["flags"].split(";")


@pytest.mark.parametrize("table", (4, 11, 25))
def test_example_6_pyl_site_through_the_steps(tmp_path: Path, table: int) -> None:
    genome, gm = _sec_world("TAG", "O")
    pep = "KDEFHWPQRSTV"
    ch = Chain(tmp_path, genome, gm, [pep], table).index()
    cells = load_stop_cells(ch.paths["sixframe"], table=table)
    assert [(s.codon_start, s.codon, s.kind) for s in cells.sites] == [
        (21, "TAG", "recoded"), (48, "TAA", "terminal"),
    ]  # fmt: skip
    (row,) = ch.loci().query("peptide == @pep").to_dict("records")
    assert "pyl_site" in row["flags"].split(";")
    assert "sec_site" not in row["flags"]


# ---------------------------------------------------------------- refusal across tables


@pytest.mark.parametrize(
    ("step", "outs"),
    [
        ("map_exact", {"hits_exact": "hx2.parquet"}),
        ("map_tolerant", {"hits_tolerant": "ht2.parquet"}),
        ("placement", {k: f"x_{v}" for k, v in PLACEMENT_OUT.items()}),
        ("loci", {"mapping": "mapping2"}),
    ],
)
def test_a_reader_refuses_an_index_built_for_another_table(
    tmp_path: Path, step: str, outs: dict[str, str]
) -> None:
    """The index is built under 4; every reader given it under 11 stops with the rebuild message
    (item 4: the format check covers files passed by hand)."""
    genome, gm = _sec_world()
    ch = Chain(tmp_path, genome, gm, ["KDEFHWPQRSTV"], 4).index()
    ch.exact()
    ch.run("map_tolerant", {"hits_tolerant": "hits_tolerant.parquet"})
    ch.run("placement", PLACEMENT_OUT)
    msg = "built for translation table 4; this run uses 11; rebuild the six-frame index"
    with pytest.raises(reg.StepError, match=msg):
        ch.run(step, outs, cfg=_cfg(11))
