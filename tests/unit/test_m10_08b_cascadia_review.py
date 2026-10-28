"""PR #334 review (M10-08b): cascadia_shards.py chunk grid, merge, markers and tmp names.

Review items 1 and 3 and the nits:

* item 1: the chunk grid must not depend on which ``reserved_<a>_<b>`` markers exist. Chunks are cut
  on the global grid ``[k*size, (k+1)*size)``; a reserved block only trims the cells it touches.
  Releasing a marker (a failed ``predict`` task) or finding its shard file must never make ``work``
  recompute a batch, and ``merge`` must then accept the result.
* item 3: ``merge`` refuses gaps, overlaps, duplicate starts and a missing tail, mixes shard and
  chunk files sorted by block start, ignores ``.tmp`` files and passes the threshold and
  ``width * cycle_time`` to ``write_results``; ``_reserved`` parses half-open markers and ignores
  other files; ``chunk_plan`` edge cases.
* nits: per-record checks in ``merge``; reserved blocks clamped to ``n_batches``; a unique ``.tmp``
  name per job; ``prep`` takes no MODEL argument (CLI and sbatch).

Cascadia, torch and lightning are never imported: ``work``, ``predict`` and ``merge`` import them
inside the function, so the tests put small fakes into ``sys.modules`` (pickle stands in for
``torch.save``/``torch.load``; a prediction is the tuple ``("pred", batch_index)``). No network, no
$GIPPA_ROOT. No statistical claim, so no calibration test.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import os
import pickle
import sys
import types
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "pilot" / "cascadia_shards.py"
SBATCH = REPO / "scripts" / "slurm" / "cascadia_array.sbatch"


def _load():  # type: ignore[no-untyped-def]
    sys.path.insert(0, str(SCRIPT.parent))
    spec = importlib.util.spec_from_file_location("cascadia_shards_review", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _load()


def _pairs(plan: list[range]) -> list[tuple[int, int]]:
    return [(b.start, b.stop) for b in plan]


def _grid_reference(n: int, reserved: list[tuple[int, int]], size: int) -> list[tuple[int, int]]:
    """Global grid cells of ``range(n)``, minus the reserved batches, as maximal runs per cell."""
    held = {b for a, c in reserved for b in range(a, c)}
    out: list[tuple[int, int]] = []
    for s in range(0, n, size):
        run: int | None = None
        for b in range(s, min(s + size, n)):
            if b in held:
                if run is not None:
                    out.append((run, b))
                run = None
            elif run is None:
                run = b
        if run is not None:
            out.append((run, min(s + size, n)))
    return out


# --- item 1: the chunk grid does not shift with the reserved set ----------------------------------


def test_review_example_grid_is_the_same_with_and_without_the_reserved_block() -> None:
    """The review's case: chunk_plan(4471, [(0, 746)], 20) must not start a grid at 746."""
    held = _pairs(S.chunk_plan(4471, [(0, 746)], 20))
    free = _pairs(S.chunk_plan(4471, [], 20))
    assert held[0] == (746, 760)
    assert held[1:] == [c for c in free if c[0] >= 760]


@given(st.integers(0, 3000), st.integers(1, 64), st.data())
def test_chunks_lie_on_the_global_grid_whatever_the_markers(
    n_batches: int, size: int, data: st.DataObject
) -> None:
    hi = n_batches + 2 * size
    blocks = data.draw(
        st.lists(st.tuples(st.integers(0, hi), st.integers(0, hi)), max_size=4).map(
            lambda xs: [(min(a, b), max(a, b)) for a, b in xs if a != b]
        )
    )
    plan = S.chunk_plan(n_batches, blocks, size)
    assert _pairs(plan) == _grid_reference(n_batches, blocks, size)
    # a grid cell that no reserved block touches is the same chunk as with no markers at all
    free = set(_pairs(S.chunk_plan(n_batches, [], size)))
    held = set(_pairs(plan))
    for s, e in free:
        if not any(a < e and s < c for a, c in blocks):
            assert (s, e) in held


# --- fakes for torch / lightning / cascadia -------------------------------------------------------


class _Recorder:
    def __init__(self) -> None:
        self.predicted: list[int] = []  # batch indices, every time one is predicted
        self.saves: list[str] = []  # paths torch.save wrote to
        self.write_results: list[tuple[Any, ...]] = []


def _install_fakes(monkeypatch: pytest.MonkeyPatch, rec: _Recorder, bs: int) -> None:
    torch = types.ModuleType("torch")

    def save(obj: Any, path: Any) -> None:
        rec.saves.append(str(path))
        with open(path, "wb") as fh:
            pickle.dump(obj, fh)

    def load(path: Any, weights_only: bool = True) -> Any:
        with open(path, "rb") as fh:
            return pickle.load(fh)

    torch.save = save  # type: ignore[attr-defined]
    torch.load = load  # type: ignore[attr-defined]
    torch.set_num_threads = lambda n: None  # type: ignore[attr-defined]
    utils = types.ModuleType("torch.utils")
    data = types.ModuleType("torch.utils.data")

    def subset(ds: Any, rows: Any) -> list[int]:
        return list(rows)

    def loader(rows: list[int], batch_size: int, **_: Any) -> list[list[int]]:
        return [rows[i : i + batch_size] for i in range(0, len(rows), batch_size)]

    data.Subset = subset  # type: ignore[attr-defined]
    data.DataLoader = loader  # type: ignore[attr-defined]
    torch.utils = utils  # type: ignore[attr-defined]
    utils.data = data  # type: ignore[attr-defined]

    pl = types.ModuleType("pytorch_lightning")

    class Trainer:
        def __init__(self, **_: Any) -> None:
            pass

        def predict(self, model: Any, dataloaders: list[list[int]]) -> list[tuple[str, int]]:
            out = []
            for batch in dataloaders:
                assert batch[0] % bs == 0, "a batch must start on a batch boundary"
                rec.predicted.append(batch[0] // bs)
                out.append(("pred", batch[0] // bs))
            return out

    pl.Trainer = Trainer  # type: ignore[attr-defined]

    cascadia = types.ModuleType("cascadia")
    cmodel = types.ModuleType("cascadia.model")

    class AugmentedSpec2Pep:
        @staticmethod
        def load_from_checkpoint(*_: Any, **__: Any) -> object:
            return object()

    cmodel.AugmentedSpec2Pep = AugmentedSpec2Pep  # type: ignore[attr-defined]
    cutils = types.ModuleType("cascadia.utils")
    cutils.write_results = lambda *a: rec.write_results.append(a)  # type: ignore[attr-defined]
    cascadia.model = cmodel  # type: ignore[attr-defined]
    cascadia.utils = cutils  # type: ignore[attr-defined]

    for name, mod in {
        "torch": torch, "torch.utils": utils, "torch.utils.data": data,
        "pytorch_lightning": pl, "cascadia": cascadia, "cascadia.model": cmodel,
        "cascadia.utils": cutils,
    }.items():  # fmt: skip
        monkeypatch.setitem(sys.modules, name, mod)

    class _DS:
        collate_fn = None

    monkeypatch.setattr(S, "_tokenizer", lambda: None)
    monkeypatch.setattr(S, "_dataset", lambda tok, asf, index: _DS())
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)


BS = 4


def _workdir(tmp_path: Path, n_batches: int, name: str = "shards") -> Path:
    wd = tmp_path / "WT_X" / name
    wd.mkdir(parents=True)
    meta = {"mzml": "/fake/WT_X.mzML", "n_spectra": (n_batches - 1) * BS + 2 if n_batches else 0,
            "batch_size": BS, "n_batches": n_batches, "isolation_window_size": 1.5,
            "cycle_time": 3.0, "width": 2, "max_charge": 4}  # fmt: skip
    (wd / S.META).write_text(json.dumps(meta))
    return wd


def _block(wd: Path, name: str, a: int, b: int, preds: list[Any] | None = None) -> None:
    n = json.loads((wd / S.META).read_text())["n_spectra"]
    rec = {"block": [a, b], "rows": max(0, min(b * BS, n) - a * BS),
           "preds": [("pred", k) for k in range(a, b)] if preds is None else preds}  # fmt: skip
    with open(wd / name, "wb") as fh:
        pickle.dump(rec, fh)


def _chunk(wd: Path, a: int, b: int, preds: list[Any] | None = None) -> None:
    _block(wd, f"chunk_{a:05d}_{b:05d}.pt", a, b, preds)


# --- item 1, behaviourally: releasing a marker never recomputes a batch ---------------------------


def test_released_marker_recomputes_nothing_and_merge_accepts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rec = _Recorder()
    _install_fakes(monkeypatch, rec, BS)
    wd = _workdir(tmp_path, 23)
    (wd / "reserved_0_7").touch()
    assert S.work("model.ckpt", 1, 1e9, 5, [wd]) == 0
    (wd / "reserved_0_7").unlink()  # the predict task for [0, 7) failed; the watcher releases it
    assert S.work("model.ckpt", 2, 1e9, 5, [wd]) == 0
    assert sorted(rec.predicted) == list(range(23)), "every batch predicted exactly once"
    assert S.merge(wd, str(tmp_path / "out"), 0.8) == 0
    assert rec.write_results[0][0] == [("pred", k) for k in range(23)]


def test_existing_shard_file_counts_as_covered_without_its_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rec = _Recorder()
    _install_fakes(monkeypatch, rec, BS)
    wd = _workdir(tmp_path, 23)
    _block(wd, "shard_000_of_003.pt", 0, 8)  # written by predict; marker already removed
    assert S.work("model.ckpt", 1, 1e9, 5, [wd]) == 0
    assert sorted(rec.predicted) == list(range(8, 23))
    assert S.merge(wd, str(tmp_path / "out"), 0.8) == 0


def test_work_skips_reserved_batches_and_tiles_the_rest(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rec = _Recorder()
    _install_fakes(monkeypatch, rec, BS)
    wd = _workdir(tmp_path, 23)
    (wd / "reserved_0_7").touch()
    assert S.work("model.ckpt", 1, 1e9, 5, [wd]) == 0
    assert sorted(rec.predicted) == list(range(7, 23))
    assert not list(wd.glob("*.tmp*")), "no tmp file left behind"


# --- nit 3: a unique tmp name per job -------------------------------------------------------------


def _tmp_saves(rec: _Recorder, wd: Path, final: str) -> list[str]:
    return [p for p in rec.saves if Path(p).name.startswith(final) and Path(p).name != final]


def test_work_tmp_name_is_unique_per_job(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rec = _Recorder()
    _install_fakes(monkeypatch, rec, BS)
    wd = _workdir(tmp_path, 3)
    final = "chunk_00000_00003.pt"
    for gen, (job, pid) in enumerate([("111", 4001), ("222", 4002)], start=1):
        monkeypatch.setenv("SLURM_JOB_ID", job)
        monkeypatch.setenv("SLURM_ARRAY_JOB_ID", job)
        monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "0")
        monkeypatch.setattr(os, "getpid", lambda pid=pid: pid)
        assert S.work("model.ckpt", gen, 1e9, 5, [wd]) == 0
        assert (wd / final).exists()
        if gen == 1:
            (wd / final).unlink()
    tmps = _tmp_saves(rec, wd, final)
    assert len(tmps) == 2, rec.saves
    assert tmps[0] != tmps[1], "two jobs on one chunk must not share a tmp path"
    assert all(Path(p).parent == wd for p in tmps), (
        "tmp must sit next to the target (atomic rename)"
    )
    assert all(str(wd / final) + ".tmp" != p for p in tmps)
    assert not list(wd.glob("*.tmp*"))


def test_predict_tmp_name_is_unique_per_job(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rec = _Recorder()
    _install_fakes(monkeypatch, rec, BS)
    wd = _workdir(tmp_path, 6)
    final = "shard_000_of_002.pt"
    for job, pid in [("111", 4001), ("222", 4002)]:
        monkeypatch.setenv("SLURM_JOB_ID", job)
        monkeypatch.setenv("SLURM_ARRAY_JOB_ID", job)
        monkeypatch.setattr(os, "getpid", lambda pid=pid: pid)
        assert S.predict("model.ckpt", wd, 0, 2, None) == 0
        assert (wd / final).exists()
    tmps = _tmp_saves(rec, wd, final)
    assert len(tmps) == 2 and tmps[0] != tmps[1], rec.saves
    assert all(Path(p).parent == wd for p in tmps)


# --- item 3: merge --------------------------------------------------------------------------------


@pytest.fixture
def merger(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):  # type: ignore[no-untyped-def]
    rec = _Recorder()
    _install_fakes(monkeypatch, rec, BS)
    wd = _workdir(tmp_path, 10)

    def run() -> int:
        return int(S.merge(wd, str(tmp_path / "WT_X"), 0.65))

    return wd, rec, run


def test_merge_single_file(merger) -> None:  # type: ignore[no-untyped-def]
    wd, rec, run = merger
    _block(wd, "shard_000_of_001.pt", 0, 10)
    assert run() == 0
    assert rec.write_results[0][0] == [("pred", k) for k in range(10)]


def test_merge_mixes_shards_and_chunks_sorted_by_block_start(merger) -> None:  # type: ignore[no-untyped-def]
    wd, rec, run = merger
    _chunk(wd, 8, 10)
    _chunk(wd, 4, 6)
    _block(wd, "shard_001_of_002.pt", 6, 8)
    _block(wd, "shard_000_of_002.pt", 0, 4)
    (wd / "chunk_00004_00006.pt.4001.tmp").write_bytes(b"torn")
    (wd / "chunk_00004_00006.pt.tmp").write_bytes(b"torn")
    assert run() == 0
    assert rec.write_results[0][0] == [("pred", k) for k in range(10)]


def test_merge_passes_paths_threshold_and_time_width(merger, tmp_path: Path) -> None:  # type: ignore[no-untyped-def]
    wd, rec, run = merger
    _chunk(wd, 0, 10)
    assert run() == 0
    (args,) = rec.write_results
    assert args[1:] == (str(tmp_path / "WT_X"), "/fake/WT_X.mzML", 1.5, 0.65, 2 * 3.0)


@pytest.mark.parametrize(
    "blocks",
    [
        pytest.param([(0, 5), (4, 10)], id="overlap"),
        pytest.param([(0, 5), (6, 10)], id="gap"),
        pytest.param([(0, 8)], id="missing-tail"),
        pytest.param([(2, 10)], id="missing-head"),
        pytest.param([(0, 10), (10, 12)], id="past-the-end"),
    ],
)
def test_merge_refuses_gaps_and_overlaps(merger, blocks) -> None:  # type: ignore[no-untyped-def]
    wd, rec, run = merger
    for a, b in blocks:
        _chunk(wd, a, b)
    assert run() == 2
    assert rec.write_results == []


def test_merge_refuses_a_duplicate_block_start(merger) -> None:  # type: ignore[no-untyped-def]
    wd, rec, run = merger
    _block(wd, "shard_000_of_002.pt", 0, 5)
    _chunk(wd, 0, 5)
    _chunk(wd, 5, 10)
    assert run() == 2
    assert rec.write_results == []


def test_merge_refuses_no_files(merger) -> None:  # type: ignore[no-untyped-def]
    _, rec, run = merger
    assert run() == 2
    assert rec.write_results == []


def test_merge_refuses_a_record_whose_preds_do_not_match_its_block(merger) -> None:  # type: ignore[no-untyped-def]
    """Nit 1: the total can balance while one record is short and another long."""
    wd, rec, run = merger
    _chunk(wd, 0, 5, preds=[("pred", k) for k in range(4)])
    _chunk(wd, 5, 10, preds=[("pred", k) for k in range(5, 11)])
    assert run() == 2
    assert rec.write_results == []


def test_merge_refuses_inconsistent_rows(merger) -> None:  # type: ignore[no-untyped-def]
    """Nit 1: ``rows`` must match the block (BS per batch, the last batch short)."""
    wd, _rec, run = merger
    _chunk(wd, 0, 10)
    with open(wd / "chunk_00000_00010.pt", "rb") as fh:
        r = pickle.load(fh)
    r["rows"] = 3
    with open(wd / "chunk_00000_00010.pt", "wb") as fh:
        pickle.dump(r, fh)
    assert run() == 2


# --- item 3: _reserved ----------------------------------------------------------------------------


def test_reserved_parses_markers_half_open_and_sorted(tmp_path: Path) -> None:
    (tmp_path / "reserved_746_900").touch()
    (tmp_path / "reserved_0_746").touch()
    assert S._reserved(tmp_path) == [(0, 746), (746, 900)]
    assert _pairs(S.chunk_plan(760, S._reserved(tmp_path), 20)) == []
    assert _pairs(S.chunk_plan(910, S._reserved(tmp_path), 20)) == [(900, 910)]


@pytest.mark.parametrize(
    "stray",
    ["reserved_0_746.bak", "reserved_0_746.tmp", "reserved_a_b", "reserved_1_2_3", "reserved__5",
     "chunk_00000_00020.pt", "shard_000_of_006.pt", "meta.json"],
)  # fmt: skip
def test_reserved_ignores_other_files(tmp_path: Path, stray: str) -> None:
    (tmp_path / "reserved_0_746").touch()
    (tmp_path / stray).touch()
    assert S._reserved(tmp_path) == [(0, 746)]


def test_reserved_empty_dir(tmp_path: Path) -> None:
    assert S._reserved(tmp_path) == []


# --- item 3 and nit 2: chunk_plan edges -----------------------------------------------------------


@pytest.mark.parametrize(
    ("n", "reserved", "size", "want"),
    [
        pytest.param(0, [], 20, [], id="no-batches"),
        pytest.param(0, [(0, 5)], 20, [], id="no-batches-reserved"),
        pytest.param(7, [], 7, [(0, 7)], id="size-equals-n"),
        pytest.param(7, [], 50, [(0, 7)], id="size-above-n"),
        pytest.param(100, [(0, 46)], 20, [(46, 60), (60, 80), (80, 100)], id="reserved-start"),
        pytest.param(30, [(25, 30)], 10, [(0, 10), (10, 20), (20, 25)], id="reserved-end"),
        pytest.param(
            50,
            [(12, 27)],
            10,
            [(0, 10), (10, 12), (27, 30), (30, 40), (40, 50)],
            id="reserved-middle",
        ),
        pytest.param(10, [(0, 3), (3, 7)], 4, [(7, 8), (8, 10)], id="adjacent"),
        pytest.param(10, [(3, 7), (0, 3)], 4, [(7, 8), (8, 10)], id="unsorted"),
        pytest.param(20, [(2, 9), (4, 6)], 5, [(0, 2), (9, 10), (10, 15), (15, 20)], id="nested"),
        pytest.param(
            20, [(2, 6), (4, 9)], 5, [(0, 2), (9, 10), (10, 15), (15, 20)], id="overlapping"
        ),
        pytest.param(10, [(0, 10)], 4, [], id="all-reserved"),
        pytest.param(10, [(15, 20)], 4, [(0, 4), (4, 8), (8, 10)], id="marker-past-end"),
        pytest.param(10, [(8, 20)], 4, [(0, 4), (4, 8)], id="marker-straddles-end"),
    ],
)
def test_chunk_plan_edges(
    n: int, reserved: list[tuple[int, int]], size: int, want: list[Any]
) -> None:
    assert _pairs(S.chunk_plan(n, reserved, size)) == want


@given(st.integers(0, 500), st.integers(1, 40), st.data())
def test_chunk_plan_stays_inside_n_batches(n: int, size: int, data: st.DataObject) -> None:
    """Nit 2: markers past the end (or unsorted, overlapping) never create chunks beyond n."""
    raw = data.draw(st.lists(st.tuples(st.integers(0, 3 * n + 50), st.integers(1, 60)), max_size=5))
    blocks = [(a, a + w) for a, w in raw]
    plan = S.chunk_plan(n, blocks, size)
    got = [b for blk in plan for b in blk]
    assert got == [b for b in range(n) if not any(a <= b < c for a, c in blocks)]
    assert all(0 < len(blk) <= size and blk.stop <= n for blk in plan)


# --- nit 6: prep takes no MODEL argument ----------------------------------------------------------


def test_prep_cli_takes_mzml_and_workdir_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "unimod_tables.xml"
    table.write_text("<unimod/>")
    calls: list[tuple[Any, ...]] = []
    monkeypatch.setattr(S.cascadia_offline, "install", lambda path: None)
    monkeypatch.setattr(S, "prep", lambda *a: calls.append(a) or 0)
    assert S.main([str(table), "prep", "/data/WT_1.mzML", str(tmp_path / "wd")]) == 0
    assert calls == [("/data/WT_1.mzML", tmp_path / "wd", 32, 2, 4)]
    with pytest.raises(SystemExit):
        S.main([str(table), "prep", "/data/WT_1.mzML", "model.ckpt", str(tmp_path / "wd")])


def test_sbatch_prep_step_passes_no_model() -> None:
    lines = [ln for ln in SBATCH.read_text().splitlines() if ln.strip().startswith("prep)")]
    assert len(lines) == 1, lines
    assert "$MODEL" not in lines[0] and "${MODEL" not in lines[0], lines[0]
    assert '"$WORK"' in lines[0]


# --- nit 4: torch, lightning and cascadia imports are function-local ------------------------------


def test_heavy_imports_in_the_shard_driver_are_function_local() -> None:
    tree = ast.parse(SCRIPT.read_text())
    heavy = {"cascadia", "torch", "pytorch_lightning", "lightning"}
    top = []
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for sub in ast.walk(node):
            if isinstance(sub, ast.Import):
                top += [a.name for a in sub.names if a.name.split(".")[0] in heavy]
            elif (
                isinstance(sub, ast.ImportFrom) and sub.module and sub.module.split(".")[0] in heavy
            ):
                top.append(sub.module)
    assert not top, top
