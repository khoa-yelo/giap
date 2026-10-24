"""Contract test: DIA-NN 2.3.1 pilot reports through ``convert_diann`` (M10-03; GIAP_SPEC.md §5.4).

Spec: docs/specs/M10-03-diann-convert.md, facts F1-F10, "The contract test", lead decisions 1-4.

Data live on FarmShare only (Khoa, 2026-10-23: nothing derived from the DIA-NN reports is
committed). The tests read

* ``$GIAP_ROOT/pilot/recorded/diann/2.3.1/<arm>/report.sample.parquet`` (M10-02's stratified
  sample, 900 targets + 100 decoys), and, for the ``slow`` ``full_report`` variant,
  ``$GIAP_ROOT/pilot/<arm>/<arm>_report.parquet``;
* the arm's search FASTAs: A1 ``$GIAP_DATA_DIR/proteomes/UP000000625.fasta``, A2 the committed
  ``validation/pilot/db_a2/strain_db.faa.gz``; both arms
  ``$GIAP_DB_DIR/contaminants/contaminants.fasta``;
* ``validation/pilot/inputs/acquisition.json`` (M10-01) for the run labels.

They **skip, never fail**, when ``GIAP_ROOT`` is unset (CI) or a file is missing, and read no
environment at import time. No byte-golden file is read or written. Recorded values (decoy rows
without q, contaminant rows, target/decoy sequence collisions) are printed and attached with
``record_property``; run with ``-s`` or ``-rA`` to see them.
"""

from __future__ import annotations

import json
import math
import os
import re
import subprocess
import sys
from collections import Counter
from functools import cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[2]
ARMS = ("A1", "A2")
SEARCH_DB = {"A1": "other", "A2": "same"}
STRAIN_DB = REPO / "validation" / "pilot" / "db_a2" / "strain_db.faa.gz"
ACQUISITION = REPO / "validation" / "pilot" / "inputs" / "acquisition.json"
CONT = "Cont_"
NO_ROOT = "GIAP_ROOT is not set: the DIA-NN samples live on FarmShare only (Khoa 2026-10-23)"
FIRST_CHOICE = (
    "Stripped.Sequence",
    "Modified.Sequence",
    "Protein.Ids",
    "Global.Q.Value",
    "PEP",
    "Decoy",
    "Precursor.Normalised",
    "Run",
    "Precursor.Charge",
    "Protein.Group",
    "PG.MaxLFQ",
)


# ------------------------------------------------------------------ data (skip when absent)


def _root() -> Path:
    root = os.environ.get("GIAP_ROOT")
    if not root:
        pytest.skip(NO_ROOT)
    return Path(root)


def _need(path: Path, why: str) -> Path:
    if not path.is_file():
        pytest.skip(f"{path} is missing: {why}")
    return path


def _outside_checkout(path: Path) -> Path:
    resolved = path.resolve()
    assert REPO.resolve() not in resolved.parents, f"{path} lies inside the checkout (Khoa's rule)"
    return path


def sample_path(arm: str) -> Path:
    p = _root() / "pilot" / "recorded" / "diann" / "2.3.1" / arm / "report.sample.parquet"
    return _outside_checkout(_need(p, "M10-02 report job not run"))


def full_report_path(arm: str) -> Path:
    p = _root() / "pilot" / arm / f"{arm}_report.parquet"
    return _outside_checkout(_need(p, "M10-02 report job not run"))


def fasta_paths(arm: str) -> list[Path]:
    _root()
    env = {k: os.environ.get(k) for k in ("GIAP_DATA_DIR", "GIAP_DB_DIR")}
    for k, v in env.items():
        if not v:
            pytest.skip(f"{k} is not set (source scripts/farmshare_env.sh)")
    cont = Path(str(env["GIAP_DB_DIR"])) / "contaminants" / "contaminants.fasta"
    arm_db = (
        Path(str(env["GIAP_DATA_DIR"])) / "proteomes" / "UP000000625.fasta"
        if arm == "A1"
        else STRAIN_DB
    )
    return [_need(arm_db, "search FASTA absent"), _need(cont, "giap setup --contaminants not run")]


@cache
def _raw(path: Path) -> pd.DataFrame:
    return pd.read_parquet(path, engine="pyarrow")


@cache
def _converted(path: Path, search_db: str) -> Any:
    from giap.config.model import InputConfig
    from giap.io.converters.diann import convert_diann

    return convert_diann([path], search_db=search_db, cfg=InputConfig())


@cache
def _fasta_ids(paths: tuple[Path, ...]) -> frozenset[str]:
    from giap.io.protein_fasta import read_protein_fasta

    return frozenset(m for p in paths for r in read_protein_fasta(p) for m in r.members)


def _record(record_property: Any, name: str, value: Any) -> None:
    record_property(name, value)
    print(f"[M10-03 recorded] {name} = {value}")


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def _decoy_mask(raw: pd.DataFrame) -> np.ndarray:
    return pd.to_numeric(raw["Decoy"].astype("float64")).to_numpy() == 1


# ------------------------------------------------------------------ shared checks (sample + full)


def check_f1_columns(path: Path) -> None:
    import pyarrow.parquet as pq

    names = set(pq.read_schema(path).names)
    missing = [c for c in FIRST_CHOICE if c not in names]
    assert not missing, (
        f"first-choice columns missing from {path.name}: {missing} "
        "(finding for the lead, not a silent fallback; spec F1)"
    )


def check_f2_types(path: Path, arm: str) -> None:
    res = _converted(path, SEARCH_DB[arm])
    bad = [
        f"{i.severity}:{i.code}: {i.message}" for i in res.issues if i.severity in ("error", "warn")
    ]
    assert not bad, bad
    df = res.peptides
    assert str(df["charge"].dtype) == "Int64" and not df["charge"].isna().any()
    assert str(df["is_decoy"].dtype) == "boolean" and not df["is_decoy"].isna().any()


def check_f6_runs(path: Path, arm: str) -> None:
    raw, df = _raw(path), _converted(path, SEARCH_DB[arm]).peptides
    assert Counter(df["run"].astype(str)) == Counter(raw["Run"].astype(str)), "run is verbatim"
    if not ACQUISITION.is_file():
        pytest.skip(f"{ACQUISITION.relative_to(REPO)} is missing (M10-01)")
    doc = json.loads(ACQUISITION.read_text())
    labels = {r["run"] for r in doc["runs"]}
    stems = {Path(str(r.get("source_folder", ""))).name.removesuffix(".d") for r in doc["runs"]}
    seen = {Path(str(v)).name.removesuffix(".d") for v in raw["Run"].astype(str).unique()}
    assert seen in (labels, stems), (sorted(seen), sorted(labels), sorted(stems))


def check_f7_protein_ids(path: Path, arm: str) -> None:
    from giap.io.protein_fasta import resolve_protein_id

    ids = _fasta_ids(tuple(fasta_paths(arm)))
    res = _converted(path, SEARCH_DB[arm])
    items: set[str] = set()
    for v in res.peptides["protein_ids"].dropna().astype(str):
        items.update(resolve_protein_id(v)[1])
    if res.abundance is not None:
        for v in res.abundance["protein_id"].astype(str):
            items.update(resolve_protein_id(v)[1])
    unknown = sorted(items - ids)
    dropped = [i for i in unknown if f"{CONT}{i}" in ids]
    assert not dropped, (
        f"DIA-NN dropped the {CONT} prefix ({dropped[:5]}): the contaminant FASTA headers must "
        "change and DIA-NN re-run; finding for M10-02 (M10-03 lead decision 3)"
    )
    assert not unknown, f"{len(unknown)} protein IDs not in the search FASTAs: {unknown[:10]}"


def check_f10_rows(path: Path, arm: str) -> None:
    raw, df = _raw(path), _converted(path, SEARCH_DB[arm]).peptides
    assert len(df) == len(raw)
    want = Counter(
        zip(
            raw["Run"].astype(str),
            raw["Stripped.Sequence"].astype(str),
            raw["Modified.Sequence"].astype(str),
            raw["Precursor.Charge"].astype(int),
            _decoy_mask(raw),
            strict=True,
        )
    )
    got = Counter(
        zip(
            df["run"].astype(str),
            df["peptide"].astype(str),
            df["modified_peptide"].astype(str),
            df["charge"].astype(int),
            df["is_decoy"].astype(bool),
            strict=True,
        )
    )
    assert got == want


# ------------------------------------------------------------------ F1-F10 on the sample


@pytest.mark.parametrize("arm", ARMS)
def test_f1_first_choice_columns(arm: str) -> None:
    check_f1_columns(sample_path(arm))


@pytest.mark.parametrize("arm", ARMS)
def test_f2_types_convert_cleanly(arm: str) -> None:
    check_f2_types(sample_path(arm), arm)


@pytest.mark.parametrize("arm", ARMS)
def test_f3_q_value_is_global_q(arm: str, record_property: Any) -> None:
    path = sample_path(arm)
    raw, df = _raw(path), _converted(path, SEARCH_DB[arm]).peptides
    q = pd.to_numeric(raw["Global.Q.Value"]).astype("float64")
    dec = _decoy_mask(raw)

    def k(q: Any) -> float | None:
        return None if _na(q) else float(q)

    want = Counter(
        zip(
            raw["Run"].astype(str),
            raw["Modified.Sequence"].astype(str),
            raw["Precursor.Charge"].astype(int),
            dec,
            map(k, q),
            strict=True,
        )
    )
    got = Counter(
        zip(
            df["run"].astype(str),
            df["modified_peptide"].astype(str),
            df["charge"].astype(int),
            df["is_decoy"].astype(bool),
            map(k, df["q_value"]),
            strict=True,
        )
    )
    assert got == want, "q_value equals float64(Global.Q.Value) row by row"
    assert not q[~dec].isna().any(), "Global.Q.Value is populated on every target row"
    present = q.dropna()
    assert ((present >= 0) & (present <= 1)).all()
    _record(record_property, f"{arm}_decoy_rows_q_na", int(q[dec].isna().sum()))


@pytest.mark.parametrize("arm", ARMS)
def test_f4_decoys(arm: str, record_property: Any) -> None:
    path = sample_path(arm)
    raw, df = _raw(path), _converted(path, SEARCH_DB[arm]).peptides
    n = int(_decoy_mask(raw).sum())
    assert int(df["is_decoy"].sum()) == n
    assert n == 100, f"the sample holds 100 decoys (M10-02), found {n}"
    ids = df.loc[df["is_decoy"].astype(bool), "protein_ids"]
    assert not ids.isna().any() and (ids.astype(str).str.len() > 0).all(), (
        "every decoy row has protein_ids (ADR-0009 pairing, ADR-0047 picked tier)"
    )
    tgt = set(df.loc[~df["is_decoy"].astype(bool), "peptide"].astype(str))
    collide = len(tgt & set(df.loc[df["is_decoy"].astype(bool), "peptide"].astype(str)))
    _record(record_property, f"{arm}_target_decoy_sequence_collisions", collide)


@pytest.mark.parametrize("arm", ARMS)
def test_f5_modifications(arm: str) -> None:
    from giap.io.converters.strip import strip_modifications

    path = sample_path(arm)
    df = _converted(path, SEARCH_DB[arm]).peptides
    tokens: set[str] = set()
    for pep, mod in zip(df["peptide"].astype(str), df["modified_peptide"].astype(str), strict=True):
        assert strip_modifications(mod) == pep, (mod, pep)
        assert re.fullmatch(r"[A-Z]{7,30}", pep), pep
        tokens.update(re.findall(r"UniMod:(\d+)", mod))
    assert tokens <= {"4", "35"}, f"modifications outside UniMod:4/35 (cfg finding): {tokens}"


@pytest.mark.parametrize("arm", ARMS)
def test_f6_run_names(arm: str) -> None:
    check_f6_runs(sample_path(arm), arm)


@pytest.mark.parametrize("arm", ARMS)
def test_f7_protein_ids_are_search_fasta_ids(arm: str) -> None:
    check_f7_protein_ids(sample_path(arm), arm)


@pytest.mark.parametrize("arm", ARMS)
def test_f8_contaminant_prefix(arm: str, record_property: Any) -> None:
    """The ``Cont_`` tag survives (F7 checks every ID); rows with one are recorded, not asserted."""
    path = sample_path(arm)
    df = _converted(path, SEARCH_DB[arm]).peptides
    ids = df["protein_ids"].fillna("").astype(str)
    n = int(ids.map(lambda s: any(x.startswith(CONT) for x in s.split(";"))).sum())
    _record(record_property, f"{arm}_contaminant_rows", n)
    check_f7_protein_ids(path, arm)


@pytest.mark.parametrize("arm", ARMS)
def test_f9_maxlfq_matches_an_independent_oracle(arm: str) -> None:
    path = sample_path(arm)
    raw = _raw(path)
    ab = _converted(path, SEARCH_DB[arm]).abundance
    assert ab is not None
    assert not ab["protein_id"].duplicated().any()
    t = raw.loc[~_decoy_mask(raw), ["Protein.Group", "Run", "PG.MaxLFQ"]].copy()
    t["Protein.Group"] = t["Protein.Group"].astype(str)
    t["Run"] = t["Run"].astype(str)
    t["v"] = pd.to_numeric(t["PG.MaxLFQ"]).astype("float64")
    t = t[t["v"].notna() & (t["v"] != 0)]
    per_run = t.groupby(["Protein.Group", "Run"])["v"].agg(["min", "max"])
    assert (per_run["min"] == per_run["max"]).all(), "one PG.MaxLFQ per (group, run)"
    want = per_run["min"].groupby(level=0).median()
    want = want[want > 0]
    got = dict(zip(ab["protein_id"].astype(str), ab["maxlfq"].astype(float), strict=True))
    assert set(got) == set(want.index)
    for g, v in want.items():
        assert got[str(g)] == pytest.approx(float(v), rel=1e-12), g


@pytest.mark.parametrize("arm", ARMS)
def test_f10_one_row_per_report_row(arm: str) -> None:
    check_f10_rows(sample_path(arm), arm)


# ------------------------------------------------------------------ the CLI and the round trip


def _cli(report: Path, out: Path, ab: Path, search_db: str) -> subprocess.CompletedProcess[str]:
    env = {**os.environ, "PYTHONPATH": str(REPO / "src")}
    cmd = [
        sys.executable,
        "-c",
        "from giap.cli import cli; cli()",
        "convert",
        "--from",
        "diann",
        str(report),
        "--out",
        str(out),
        "--abundance-out",
        str(ab),
        "--search-db",
        search_db,
    ]
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600, check=False)


@pytest.mark.parametrize("arm", ARMS)
def test_cli_is_deterministic_and_round_trips(arm: str, tmp_path: Path) -> None:
    from giap.io.abundance import read_abundance
    from giap.io.peptides import read_peptide_table

    path = sample_path(arm)
    outs = []
    for k in (1, 2):
        d = tmp_path / str(k)
        d.mkdir()
        res = _cli(path, d / "a.tsv.gz", d / "ab.tsv.gz", SEARCH_DB[arm])
        assert res.returncode == 0, res.stderr
        outs.append(((d / "a.tsv.gz").read_bytes(), (d / "ab.tsv.gz").read_bytes()))
    assert outs[0] == outs[1], "a second call gives byte-identical files"
    df, issues = read_peptide_table(tmp_path / "1" / "a.tsv.gz")
    assert not [i for i in issues if i.severity == "error"], issues
    assert len(df) == len(_raw(path))
    _, ab_issues = read_abundance(tmp_path / "1" / "ab.tsv.gz")
    assert not [i for i in ab_issues if i.severity == "error"], ab_issues


# ------------------------------------------------------------------ the full report (slow)


@pytest.mark.slow
@pytest.mark.parametrize("check", ["f1", "f2", "f6", "f7", "f10"])
@pytest.mark.parametrize("arm", ARMS)
def test_full_report(arm: str, check: str) -> None:
    """F1, F2, F6, F7 and F10 on ``$GIAP_ROOT/pilot/<arm>/<arm>_report.parquet`` (job step 7)."""
    path = full_report_path(arm)
    if check == "f1":
        check_f1_columns(path)
    elif check == "f2":
        check_f2_types(path, arm)
    elif check == "f6":
        check_f6_runs(path, arm)
    elif check == "f7":
        check_f7_protein_ids(path, arm)
    else:
        check_f10_rows(path, arm)
