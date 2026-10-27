"""The ``loci`` step and the ``readthrough`` flag (M7-16b rule 16; spec Q14 lead decision: an
**optional** ``sixframe`` input, always wired in ``mapping.nf``; acceptance "A format-1 index given
to ``loci`` raises. Without ``sixframe``, ``loci`` emits ``mapping.readthrough_unchecked`` and flags
nothing. ``mapping.nf``'s ``loci`` ``stepIn`` includes the six-frame index").

World: ``helpers.readthrough_index.world()`` with two database peptides, one over ``A1``'s TAA
(``a:[297, 300) +``, frame +0, p = 99) read as ``W`` and one exact interior peptide. The chain runs
``sixframe_index -> map_exact -> map_tolerant -> placement``; the test then accepts both peptides in
the ``placement`` table (a two-peptide toy has no decoy power, spec methods item 7), so ``loci``
sees one tolerant row that consumes ``A1``'s cell.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from helpers.readthrough_index import FIXTURE, SEED, inputs, world
from helpers.readthrough_toys import write_annotation

from gippa.config import resolve_config
from gippa.mapping.sixframe import build_sixframe, decode
from gippa.steps import load_step_modules
from gippa.steps import registry as reg

EVENT = "mapping.readthrough_unchecked"
REPO = Path(__file__).resolve().parents[2]
MAPPING_NF = REPO / "workflow" / "subworkflows" / "local" / "mapping.nf"


@pytest.fixture(autouse=True)
def _steps() -> None:
    load_step_modules()


def _cfg(relation: str = "same"):  # type: ignore[no-untyped-def]
    return resolve_config(overrides=[f"seed={SEED}", f"genome_relation={relation}"])


def _peptides(tmp: Path) -> tuple[Path, str, str, tuple[int, int]]:
    genome, reps, _ = world()
    f0 = decode(build_sixframe(genome, reps).frames[0].aa)
    assert f0[99] == "*"
    a = 92
    assert f0[a : a + 15].count("*") == 1 and "X" not in f0[a : a + 15]
    rt = f0[a:99] + "W" + f0[100 : a + 15]
    inner = next(f0[s : s + 12] for s in range(0, 80) if "*" not in f0[s : s + 12])
    df = pd.DataFrame(
        {
            "peptide": [rt, inner],
            "source": ["database"] * 2,
            "q_value": [0.001] * 2,
            "pep": [0.01] * 2,
            "score": [np.nan] * 2,
            "is_decoy": [False] * 2,
            "is_contaminant": [False] * 2,
            "protein_ids": [["REF_1"], ["REF_2"]],
            "modified_peptides": [[rt], [inner]],
            "n_rows": [1, 1],
            "n_runs": [1, 1],
            "intensity": [1.0e6] * 2,
            "search_db": ["other"] * 2,
        }
    )
    path = tmp / "peptides_rt.parquet"
    df.to_parquet(path, index=False)
    return path, rt, inner, (3 * a, 3 * (a + 15))


def _chain(tmp: Path, relation: str = "same") -> tuple[dict[str, Path], str, str, tuple[int, int]]:
    """Every ``loci`` input (``sixframe`` included) after the chain, with both peptides accepted."""
    cfg = _cfg(relation)
    ins = inputs(tmp)
    _g, _r, gm = world()
    ins["features"] = write_annotation(gm, tmp / "ann2")["features"]
    ins["peptides"], rt, inner, iv = _peptides(tmp)
    paths = dict(ins)

    def run(name: str, outs: dict[str, str]) -> None:
        spec = reg.get(name)
        given = {k: paths[k] for k in spec.inputs}
        given |= {k: paths[k] for k in spec.optional_inputs if k in paths}
        out_paths = {k: tmp / v for k, v in outs.items()}
        reg.run_step(name, cfg, given, out_paths)
        paths.update(out_paths)

    run("sixframe_index", {"sixframe": "sixframe.npz"})
    run("map_exact", {"hits_exact": "hits_exact.parquet"})
    run("map_tolerant", {"hits_tolerant": "hits_tolerant.parquet"})
    run(
        "placement",
        {
            "placement": "placement.parquet",
            "mapping_fdr": "mapping_fdr.parquet",
            "stratum0_classes": "s0.parquet",
            "placement_costs": "costs.parquet",
        },
    )
    per = pd.read_parquet(paths["placement"])
    assert set(per["peptide"]) == {rt, inner}, per
    per["accepted"] = True
    per["placement_q"] = 0.001
    per.to_parquet(paths["placement"], index=False)
    return paths, rt, inner, iv


def _loci(
    tmp: Path,
    paths: dict[str, Path],
    relation: str = "same",
    *,
    six: bool = True,
    sixframe: Path | None = None,
) -> Path:
    spec = reg.get("loci")
    given = {k: paths[k] for k in spec.inputs}
    if six:
        given["sixframe"] = sixframe or paths["sixframe"]
    out = tmp / "mapping"
    reg.run_step("loci", _cfg(relation), given, {"mapping": out})
    return out


def _hits(out: Path) -> pd.DataFrame:
    return pd.read_csv(out / "peptide_hits.tsv", sep="\t", dtype=str, keep_default_na=False)


def _gff_attrs(out: Path, name: str) -> dict[str, str]:
    for ln in (out / "mappings.gff3").read_text().splitlines():
        f = ln.split("\t")
        if (
            len(f) == 9
            and f[2] == "nucleotide_to_protein_match"
            and f";Name={name};" in f";{f[8]};"
        ):
            return dict(kv.split("=", 1) for kv in f[8].split(";"))
    raise AssertionError(f"no GFF3 feature named {name}")


def _events(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    return [r for r in caplog.records if getattr(r, "event", None) == EVENT]


# ---------------------------------------------------------------- registration and wiring


def test_loci_declares_sixframe_as_an_optional_input() -> None:
    spec = reg.get("loci")
    assert "sixframe" in spec.optional_inputs
    assert "sixframe" not in spec.inputs


def test_mapping_nf_wires_the_six_frame_index_into_loci() -> None:
    text = re.sub(r"(?m)//.*$", "", MAPPING_NF.read_text(encoding="utf-8"))
    m = re.search(r"stepIn\(\s*'loci'\s*,(.*)", text)
    assert m, "no stepIn('loci', ...) in mapping.nf"
    six_var = re.search(r"def\s+(\w+)\s*=\s*kv\(\s*'sixframe'", text)
    assert six_var, "no kv('sixframe', ...) binding in mapping.nf"
    tokens = set(re.findall(r"\b\w+\b", m.group(1)))
    assert six_var.group(1) in tokens or "sixframe" in tokens, m.group(0)


# ---------------------------------------------------------------- the flag through the step


@pytest.mark.parametrize(
    ("relation", "want"), [("same", "readthrough"), ("related", "readthrough;strain_difference")]
)
def test_loci_with_sixframe_flags_the_consuming_row(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, relation: str, want: str
) -> None:
    caplog.set_level(logging.INFO)
    paths, rt, inner, (s, e) = _chain(tmp_path, relation)
    out = _loci(tmp_path, paths, relation)
    t = _hits(out).set_index("peptide")
    assert (t.loc[rt, "start"], t.loc[rt, "end"], t.loc[rt, "route"]) == (
        str(s + 1),
        str(e),
        "tolerant",
    )
    assert t.loc[rt, "flags"] == want
    assert "readthrough" not in t.loc[inner, "flags"]
    assert "readthrough_codon" not in (out / "peptide_hits.tsv").read_text()
    a = _gff_attrs(out, rt)
    assert a["readthrough_codon"] == "298"  # A1's TAA [297, 300) -> 1-based 298
    assert "readthrough_codon" not in _gff_attrs(out, inner)
    assert _events(caplog) == []


def test_loci_without_sixframe_warns_once_and_flags_nothing(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.INFO)
    paths, rt, _inner, _iv = _chain(tmp_path)
    out = _loci(tmp_path, paths, six=False)
    assert not any("readthrough" in f for f in _hits(out)["flags"])
    assert "readthrough_codon" not in (out / "mappings.gff3").read_text()
    ev = _events(caplog)
    assert len(ev) == 1 and ev[0].levelno == logging.WARNING
    assert rt in set(_hits(out)["peptide"])  # the row is still written, only unflagged


def test_loci_with_a_format_1_index_raises(tmp_path: Path) -> None:
    paths, _rt, _inner, _iv = _chain(tmp_path)
    with pytest.raises(reg.StepError, match="rebuild the six-frame index"):
        _loci(tmp_path, paths, sixframe=FIXTURE)
