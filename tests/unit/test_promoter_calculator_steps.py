"""Steps ``promoter_calculator`` (ext) and ``tx_layer`` (M2-14; M2-index step contract;
GIPPA_SPEC.md §6.2, §7.1.1, §7.1.4).

Contract pinned here (test author; the implementer keeps ``workflow/tests/data/annotate/
step_registry.json`` and ``annotate.nf`` consistent with it):

* ``promoter_calculator`` (ext): ``tu_upstream`` → ``pc``; ``precomputed="tx_rates"``. One call per
  ``tu_upstream/<tu_id>.fa``; the CSV is ``<tu_id>.csv`` (the merged PROMOTER_CALCULATOR module).
* ``tx_layer``: inputs ``pc`` (one CSV or a list, one per TU, the file stem is the TU id),
  ``tus`` and ``replicons``; output ``tx_rates`` (``tx_rates.parquet``); ``precomputed="tx_rates"``.
  The M2-index row lists only ``pc/*.csv`` and ``tus.parquet``, but genome TSSs of origin-crossing
  regions must be reduced modulo the replicon length (``tx_layer(tus, hits, replicons)``), so
  ``annotate.nf`` must pass ``reps`` to ``stepIn('tx_layer', ...)`` like ``ostir_input``.
* ``tx_rates.parquet`` columns exactly ``contig, strand, tss, tx_rate, tu_id`` (string, string,
  int64 0-based, float64, string); the first four are the precomputed ``tx_rates`` layer's
  columns. Deterministic: byte-identical on rerun and independent of the order of the ``pc`` files.
"""

from __future__ import annotations

import importlib
import json
import re
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from helpers.promoter_calculator_fixtures import header_only, write_pc_csv

from gippa.config import GippaConfig

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "annotate.nf"
FIXTURE = REPO / "workflow" / "tests" / "data" / "annotate" / "step_registry.json"
_TU_COLS = [
    ("tu_id", pa.string()), ("contig", pa.string()), ("start", pa.int64()), ("end", pa.int64()),
    ("strand", pa.string()), ("genes", pa.list_(pa.string())), ("upstream_start", pa.int64()),
    ("upstream_end", pa.int64()),
]  # fmt: skip
_REP_COLS = [
    ("contig", pa.string()), ("length", pa.int64()), ("topology", pa.string()),
    ("kind", pa.string()),
]  # fmt: skip


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


def _table(path: Path, cols: list[tuple[str, pa.DataType]], rows: list[tuple[object, ...]]) -> Path:
    schema = pa.schema([pa.field(n, t) for n, t in cols])
    arrays = [pa.array([r[i] for r in rows], type=t) for i, (_, t) in enumerate(cols)]
    pq.write_table(pa.Table.from_arrays(arrays, schema=schema), path)
    return path


@pytest.fixture
def work(tmp_path: Path) -> dict[str, object]:
    reps = _table(
        tmp_path / "replicons.parquet",
        _REP_COLS,
        [("chr1", 10_000, "circular", "chromosome"), ("p1", 2_000, "linear", "plasmid")],
    )
    tus = _table(
        tmp_path / "tus.parquet",
        _TU_COLS,
        [
            ("TU_1", "chr1", 1000, 1600, "+", ["g1"], 850, 1050),
            ("TU_2", "chr1", 3000, 3300, "+", ["g2"], 2900, 3050),
            ("TU_3", "chr1", 5000, 5600, "-", ["g3", "g4"], 5550, 5800),
            ("TU_4", "chr1", 100, 700, "+", ["g5"], 9800, 10150),  # region across the origin
            ("TU_5", "p1", 200, 500, "-", ["g6"], 450, 700),
        ],
    )
    pc = tmp_path / "pc"
    csvs = [
        write_pc_csv(pc / "TU_1.csv", [(101, "+", 1234.5), (60, "-", 999.0), (121, "+", 30.0)]),
        header_only(pc / "TU_2.csv"),
        write_pc_csv(pc / "TU_3.csv", [(51, "+", 80.0)], crlf=False),
        write_pc_csv(pc / "TU_4.csv", [(251, "+", 2.0), (400, "+", 77.0)]),  # 400 > 350: dropped
        write_pc_csv(pc / "TU_5.csv", [(1, "+", 4.0)]),
    ]
    return {"replicons": reps, "tus": tus, "pc": csvs}


WANT = {
    ("chr1", "+", 950, 1234.5, "TU_1"),
    ("chr1", "+", 970, 30.0, "TU_1"),
    ("chr1", "-", 5749, 80.0, "TU_3"),
    ("chr1", "+", 50, 2.0, "TU_4"),
    ("p1", "-", 699, 4.0, "TU_5"),
}


def _run(reg: ModuleType, work: dict[str, object], out: Path, **over: object) -> Path:
    ins = {**work, **over}
    reg.run_step("tx_layer", GippaConfig(seed=1), ins, {"tx_rates": out})
    return out


# ---------------------------------------------------------------- registry


def test_promoter_calculator_is_a_registered_external_step(steps: ModuleType) -> None:
    s = steps.get("promoter_calculator")
    assert s.kind == "ext" and s.func is None
    assert s.inputs == ("tu_upstream",) and s.outputs == ("pc",)
    assert s.precomputed == "tx_rates" and s.module == "core"


def test_tx_layer_is_a_registered_py_step(steps: ModuleType) -> None:
    s = steps.get("tx_layer")
    assert s.kind == "py"
    assert set(s.inputs) == {"pc", "tus", "replicons"} and s.optional_inputs == ()
    assert s.outputs == ("tx_rates",)
    assert s.precomputed == "tx_rates" and s.module == "core"


def test_tx_layer_runs_after_operon_pairs(steps: ModuleType) -> None:
    order = steps.topo_order(["tx_layer", "promoter_calculator", "operon_pairs", "parse_bakta"])
    assert (
        order.index("operon_pairs") < order.index("promoter_calculator") < order.index("tx_layer")
    )


def test_step_registry_fixture_has_tx_layer() -> None:
    got = json.loads(FIXTURE.read_text(encoding="utf-8"))
    assert got.get("tx_layer") == {
        "inputs": ["pc", "replicons", "tus"],
        "optional_inputs": [],
        "outputs": ["tx_rates"],
    }


def test_annotate_nf_passes_replicons_to_tx_layer() -> None:
    text = re.sub(r"(?m)//.*$", "", SUBWF.read_text(encoding="utf-8"))
    m = re.search(r"stepIn\(\s*'tx_layer'\s*,\s*\[(.*?)\]\s*\)\s*\n", text, flags=re.S)
    assert m, "no stepIn('tx_layer', [...]) in annotate.nf"
    assert re.search(r"\breps\b", m.group(1)), f"tx_layer inputs lack reps: {m.group(1)}"
    assert re.search(r"\btus\b", m.group(1)) and "'pc'" in m.group(1)


# ---------------------------------------------------------------- run_step


def test_tx_layer_step_writes_the_layer(
    steps: ModuleType, work: dict[str, object], tmp_path: Path
) -> None:
    out = _run(steps, work, tmp_path / "out" / "tx_rates.parquet")
    t = pq.read_table(out)
    assert t.column_names == ["contig", "strand", "tss", "tx_rate", "tu_id"]
    assert [f.type for f in t.schema] == [
        pa.string(),
        pa.string(),
        pa.int64(),
        pa.float64(),
        pa.string(),
    ]
    rows = {tuple(r.values()) for r in t.to_pylist()}
    assert rows == WANT
    assert len(t) == len(WANT)


def test_tx_layer_step_is_deterministic_and_order_independent(
    steps: ModuleType, work: dict[str, object], tmp_path: Path
) -> None:
    a = _run(steps, work, tmp_path / "a" / "tx_rates.parquet")
    b = _run(steps, work, tmp_path / "b" / "tx_rates.parquet")
    c = _run(steps, work, tmp_path / "c" / "tx_rates.parquet", pc=list(reversed(work["pc"])))  # type: ignore[call-overload]
    assert a.read_bytes() == b.read_bytes() == c.read_bytes()


def test_tx_layer_step_accepts_a_single_csv(
    steps: ModuleType, work: dict[str, object], tmp_path: Path
) -> None:
    one = work["pc"][0]  # type: ignore[index]
    out = _run(steps, work, tmp_path / "tx_rates.parquet", pc=one)
    rows = {tuple(r.values()) for r in pq.read_table(out).to_pylist()}
    assert rows == {r for r in WANT if r[4] == "TU_1"}


def test_tx_layer_step_with_only_header_only_csvs_writes_an_empty_layer(
    steps: ModuleType, work: dict[str, object], tmp_path: Path
) -> None:
    pc = [header_only(tmp_path / "pc0" / f"TU_{i}.csv") for i in (1, 2, 3)]
    out = _run(steps, work, tmp_path / "tx_rates.parquet", pc=pc)
    t = pq.read_table(out)
    assert len(t) == 0
    assert t.column_names == ["contig", "strand", "tss", "tx_rate", "tu_id"]


def test_tx_layer_tss_is_a_valid_precomputed_position(
    steps: ModuleType, work: dict[str, object], tmp_path: Path
) -> None:
    """0 <= tss < length: the layer converts 1:1 into a valid precomputed ``tx_rates`` TSV."""
    from gippa.annotate.precomputed import PRECOMPUTED_COLUMNS

    out = _run(steps, work, tmp_path / "tx_rates.parquet")
    lengths = {"chr1": 10_000, "p1": 2_000}
    for r in pq.read_table(out).to_pylist():
        assert 0 <= r["tss"] < lengths[r["contig"]]
    assert list(PRECOMPUTED_COLUMNS["tx_rates"]) == ["contig", "strand", "tss", "tx_rate"]
