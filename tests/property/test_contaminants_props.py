"""M3-15 properties: contaminant matching vs brute-force oracles, unmapped-bucket invariants.

Input-frame contracts as in tests/unit/test_contaminants.py (module docstring).
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
from contaminant_oracle import (
    OUTFMT_COLUMNS,
    blast_line,
    il,
    oracle_best,
    oracle_passing,
    oracle_peptide_matches,
)

from giap.io.protein_fasta import read_protein_fasta
from giap.mapping.placement_fdr import PER_PEPTIDE_COLUMNS

PREFIXES = ["CON__", "Cont_"]
AA = "ACDEFGHIKLMNPQRSTVWY"
SMALL = "AILK"  # small alphabet: many I/L collisions and substring hits
SETTINGS = settings(
    max_examples=150, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)


def _c() -> ModuleType:
    return importlib.import_module("giap.mapping.contaminants")


def _library(tmp: Path, entries: list[tuple[str, str]]) -> list[Any]:
    fa = tmp / "lib.faa"
    fa.write_text("".join(f">{i} x\n{s}\n" for i, s in entries))
    return read_protein_fasta(fa)


lib_entries = st.dictionaries(
    st.from_regex(r"[A-Z][0-9]{1,3}", fullmatch=True),
    st.text(SMALL, min_size=7, max_size=40),
    max_size=5,
).map(lambda d: sorted(d.items()))
peptide_lists = st.lists(st.text(SMALL, min_size=7, max_size=12), max_size=15)


@SETTINGS
@given(entries=lib_entries, peps=peptide_lists)
def test_peptide_matches_equal_oracle(
    tmp_path_factory: Any, entries: list[tuple[str, str]], peps: list[str]
) -> None:
    lib = _library(tmp_path_factory.mktemp("lib"), entries)
    got = _c().match_contaminant_peptides(peps, lib)
    assert dict(got) == oracle_peptide_matches(peps, entries)


@SETTINGS
@given(entries=lib_entries, peps=peptide_lists, data=st.data())
def test_peptide_matches_i_l_invariant_and_order_free(
    tmp_path_factory: Any, entries: list[tuple[str, str]], peps: list[str], data: st.DataObject
) -> None:
    c = _c()
    lib = _library(tmp_path_factory.mktemp("lib"), entries)
    got = c.match_contaminant_peptides(peps, lib)
    flipped = [p.replace("L", "#").replace("I", "L").replace("#", "I") for p in peps]
    got2 = c.match_contaminant_peptides(flipped, lib)
    for p, f in zip(peps, flipped, strict=True):
        assert got.get(p) == got2.get(f)
    perm = data.draw(st.permutations(lib))
    assert dict(c.match_contaminant_peptides(list(reversed(peps)), perm)) == dict(got)
    for p, lid in got.items():
        seq = next(s for i, s in entries if i == lid)
        assert il(p) in il(seq)


hsp = st.fixed_dictionaries(
    {
        "qseqid": st.sampled_from(["Q1", "Q2", "Q3"]),
        "sseqid": st.sampled_from(["A1", "B2", "C3", "D4"]),
        "pident": st.sampled_from([50.0, 85.0, 89.999, 90.0, 92.5, 100.0]),
        "bitscore": st.sampled_from([150, 300]),
        "qlen": st.integers(10, 200),
        "a": st.floats(0, 1),
        "b": st.floats(0, 1),
    }
)


def _line(h: dict[str, Any]) -> str:
    qlen = h["qlen"]
    lo, hi = sorted((1 + int(h["a"] * (qlen - 1)), 1 + int(h["b"] * (qlen - 1))))
    return blast_line(
        qseqid=h["qseqid"],
        sseqid=h["sseqid"],
        pident=h["pident"],
        bitscore=h["bitscore"],
        qstart=lo,
        qend=hi,
        qlen=qlen,
    )


@SETTINGS
@given(hsps=st.lists(hsp, max_size=12), data=st.data())
def test_protein_matches_agree_with_oracle(
    tmp_path_factory: Any, hsps: list[dict[str, Any]], data: st.DataObject
) -> None:
    lines = [_line(h) for h in hsps]
    tmp = tmp_path_factory.mktemp("bp")
    path = tmp / "blastp.tsv"
    path.write_text("".join(x + "\n" for x in lines))
    c = _c()
    got = dict(c.match_contaminant_proteins(path))
    rows = [dict(zip(OUTFMT_COLUMNS, x.split("\t"), strict=True)) for x in lines]
    passing = oracle_passing(rows)
    assert set(got) == set(passing)
    for q, s in got.items():
        assert s in passing[q]
    assert got == oracle_best(rows)  # M3-15 lead decision 2 tie rule
    # deterministic under row order
    path.write_text("".join(x + "\n" for x in data.draw(st.permutations(lines))))
    assert dict(c.match_contaminant_proteins(path)) == got


# --------------------------------------------------------------------- unmapped_table


@st.composite
def scenario(draw: st.DrawFn) -> dict[str, Any]:
    peps = draw(st.lists(st.text(AA, min_size=7, max_size=20), unique=True, max_size=12))
    state: dict[str, tuple[bool, bool]] = {}  # peptide -> (has_hit, accepted)
    ids: dict[str, str | None] = {}
    for p in peps:
        has = draw(st.booleans())
        state[p] = (has, has and draw(st.booleans()))
        ids[p] = draw(
            st.sampled_from([None, "GSIM_1", "CON__X1", "CON__A;CON__B", "CON__A;GSIM_2", "Cont_Z"])
        )
    pep_matches = {p: draw(st.sampled_from(["L1", "L2"])) for p in peps if draw(st.booleans())}
    prots = draw(
        st.lists(
            st.tuples(
                st.sampled_from(["P1", "P2", "P3", "CON__P4", "P5"]),
                st.booleans(),
                st.booleans(),
            ),
            unique_by=lambda t: t[0],
            max_size=5,
        )
    )
    prot_matches = {p: "ALB" for p, _, _ in prots if draw(st.booleans())}
    return {
        "peps": peps,
        "state": state,
        "ids": ids,
        "pep_matches": pep_matches,
        "prots": prots,
        "prot_matches": prot_matches,
        "min_len": draw(st.integers(7, 12)),
        "seed": draw(st.integers(0, 2**16)),
    }


def _frames(s: dict[str, Any], seed: int | None = None) -> tuple[pd.DataFrame, ...]:
    peps = pd.DataFrame(
        {
            "peptide": s["peps"],
            "protein_ids": pd.Series(
                [s["ids"][p].split(";") if s["ids"][p] else [] for p in s["peps"]], dtype=object
            ),
            "is_decoy": [False] * len(s["peps"]),
            "q_value": [0.0] * len(s["peps"]),
            "is_contaminant": [False] * len(s["peps"]),
        }
    )
    placed = [p for p in s["peps"] if s["state"][p][0]]
    plc = pd.DataFrame(
        {
            "peptide": placed,
            "source": ["database"] * len(placed),
            "decoy_type": ["reverse"] * len(placed),
            "length": [len(p) for p in placed],
            "m_target": [0] * len(placed),
            "e_target": [1e-3] * len(placed),
            "placement_q": [0.001 if s["state"][p][1] else 0.2 for p in placed],
            "accepted": pd.Series([s["state"][p][1] for p in placed], dtype=bool),
            "flags": [""] * len(placed),
        },
        columns=list(PER_PEPTIDE_COLUMNS),
    )
    prots = pd.DataFrame(
        {
            "protein_id": [p for p, _, _ in s["prots"]],
            "sequence": ["M" + "A" * (50 + i) for i in range(len(s["prots"]))],
            "assigned": pd.Series([a for _, a, _ in s["prots"]], dtype=bool),
            "has_alignment": pd.Series([h for _, _, h in s["prots"]], dtype=bool),
        }
    )
    if seed is not None:
        peps, plc, prots = (
            df.sample(frac=1, random_state=seed).reset_index(drop=True) for df in (peps, plc, prots)
        )
    return peps, plc, prots


def _run(s: dict[str, Any], seed: int | None = None) -> list[Any]:
    peps, plc, prots = _frames(s, seed)
    return list(
        _c().unmapped_table(
            peps,
            plc,
            prots,
            s["pep_matches"],
            s["prot_matches"],
            PREFIXES,
            tolerant_min_len=s["min_len"],
        )
    )


def _all_prefixed(ids: str | None) -> str | None:
    if ids is None:
        return None
    parts = ids.split(";")
    return parts[0] if all(x.startswith(tuple(PREFIXES)) for x in parts) else None


@settings(max_examples=200, deadline=None)
@given(s=scenario())
def test_unmapped_table_matches_rules(s: dict[str, Any]) -> None:
    c = _c()
    rows = _run(s)
    keys = [(r.kind, r.sequence_id) for r in rows]
    assert keys == sorted(keys) and len(set(keys)) == len(keys)
    got = {(r.kind, r.sequence_id): r for r in rows}

    expected: dict[tuple[str, str], tuple[str, str, str | None]] = {}
    for p in s["peps"]:
        has, acc = s["state"][p]
        reason = c.peptide_reason(has, acc, len(p), s["min_len"])
        if acc:
            assert reason is None
            continue
        assert reason == (
            "rejected_placement_fdr"
            if has
            else ("too_short" if len(p) < s["min_len"] else "no_hit")
        )
        if p in s["pep_matches"]:
            expected[("peptide", p)] = ("contaminant", reason, s["pep_matches"][p])
        elif (first := _all_prefixed(s["ids"][p])) is not None:
            expected[("peptide", p)] = ("contaminant", "contaminant_id", first)
        else:
            expected[("peptide", p)] = ("unmapped", reason, None)
    for pid, assigned, aligned in s["prots"]:
        if assigned:
            continue
        reason = "below_identity" if aligned else "no_hit"
        if pid in s["prot_matches"]:
            expected[("protein", pid)] = ("contaminant", reason, s["prot_matches"][pid])
        elif pid.startswith(tuple(PREFIXES)):
            expected[("protein", pid)] = ("contaminant", "contaminant_id", pid)
        else:
            expected[("protein", pid)] = ("unmapped", reason, None)

    assert {k: (r.label, r.reason, r.contaminant_match) for k, r in got.items()} == expected
    for (kind, sid), r in got.items():
        if kind == "peptide":
            assert r.sequence == sid


@settings(max_examples=100, deadline=None)
@given(s=scenario())
def test_unmapped_table_order_independent(s: dict[str, Any]) -> None:
    assert _run(s, seed=s["seed"]) == _run(s)
