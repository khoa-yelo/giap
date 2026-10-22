"""Contract test: ``parse_ostir_starts`` on a recorded multi-start OSTIR 1.1.3 run (issue M7-08
rule 9, example 5; GIAP_SPEC.md §5.4, §7.1.1, §7.6.2).

Rule 9: one real OSTIR run with a multi-start ``start``/``end`` window on three sim100k
candidates, recorded under ``tests/data/recorded/ostir/1.1.3/sim100k_novel/`` with
``PROVENANCE.md`` (``ostir_input.csv`` as ``novel_candidates`` writes it, ``ostir.csv`` as OSTIR
wrote it). The existing recordings only hold single-start rows.

Requested starts of an input row = the 1-based window positions ``p`` in ``[start, end]``, in frame
with ``start`` (``p - start`` divisible by 3), whose codon ``seq[p-1:p+2]`` is ATG, GTG or TTG
(``annotation.start_codons``) -- exactly the ORF's candidate starts (M7-06 rule 8).

Example 5: one feature row per requested start that OSTIR scored, keyed by the recorded
1-based ``start_position``, with exactly the recorded values; positions OSTIR reported outside
the requested set are dropped.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

REC = (
    Path(__file__).resolve().parents[1] / "data" / "recorded" / "ostir" / "1.1.3" / "sim100k_novel"
)
STARTS = ("ATG", "GTG", "TTG")


def _dicts(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


@pytest.fixture(scope="module")
def rec() -> Path:
    for name in ("ostir_input.csv", "ostir.csv", "PROVENANCE.md"):
        assert (REC / name).is_file(), f"rule 9 recording missing: {REC / name}"
    return REC


def _requested(inp: Path) -> dict[str, set[int]]:
    out: dict[str, set[int]] = {}
    for r in _dicts(inp):
        s, e, seq = int(r["start"]), int(r["end"]), r["seq"]
        out[r["id"]] = {p for p in range(s, e + 1, 3) if seq[p - 1 : p + 2] in STARTS}
    return out


def test_recording_is_three_novel_candidates_with_a_multi_start_window(rec: Path) -> None:
    rows = _dicts(rec / "ostir_input.csv")
    assert list(rows[0]) == ["id", "seq", "start", "end", "anti-Shine-Dalgarno"]
    assert len(rows) == 3
    assert all(r["id"].startswith("NOVEL_") for r in rows)
    assert any(int(r["end"]) > int(r["start"]) for r in rows)
    prov = (rec / "PROVENANCE.md").read_text(encoding="utf-8")
    assert "1.1.3" in prov and "sim100k" in prov


def test_recorded_output_has_several_starts_for_one_id(rec: Path) -> None:
    req = _requested(rec / "ostir_input.csv")
    seen: dict[str, set[int]] = {}
    for r in _dicts(rec / "ostir.csv"):
        if int(r["start_position"]) in req.get(r["name"], set()):
            seen.setdefault(r["name"], set()).add(int(r["start_position"]))
    assert any(len(v) >= 2 for v in seen.values()), seen


def test_example_5_parse_ostir_starts_on_the_recording(rec: Path) -> None:
    from giap.annotate.ostir import RbsFeatures, parse_ostir_starts

    req = _requested(rec / "ostir_input.csv")
    got = parse_ostir_starts(rec / "ostir.csv", req)
    want: dict[str, dict[int, dict[str, str]]] = {}
    for r in _dicts(rec / "ostir.csv"):
        pos = int(r["start_position"])
        if pos in req.get(r["name"], set()):
            want.setdefault(r["name"], {}).setdefault(pos, r)
    assert {k: set(v) for k, v in got.items()} == {k: set(v) for k, v in want.items()}
    for name, by_pos in want.items():
        assert set(by_pos) <= req[name]
        for pos, r in by_pos.items():
            f = got[name][pos]
            assert isinstance(f, RbsFeatures)
            assert f.tir_ostir == float(r["expression"])
            assert f.sd_dg == float(r["dG_rRNA:mRNA"])
            assert f.mrna_dg == float(r["dG_mRNA"])
            assert f.spacing_dg == float(r["dG_spacing"])
            assert f.standby_dg == float(r["dG_standby"])
            assert f.spacer_bp == int(r["RBS_distance_bp"])
            assert f.start_codon == r["start_codon"]
