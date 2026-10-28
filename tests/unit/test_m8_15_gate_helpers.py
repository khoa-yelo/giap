"""Fast unit tests of the M8-15 gate glue (PR #266 code review B1, B2, nits 1 and 3).

- ``prefetch_indices`` (B2): unsharded -> block 1 ``0..R-1``; ``--shard`` block 1 -> the shard's
  range; ``--shard --block 2`` -> the shard's range offset by R (the replicates ``calibrate``
  runs); confirm run -> none.
- ``sharded`` (B1): true under ``GIPPA_CALIBRATION_REPLICATES`` or ``GIPPA_CALIBRATION_CONFIRM``.
- ``cpu_share`` / ``workers`` (nit 1): the job's CPUs divided by the xdist worker count.
- ``pick_one_per_collection`` (nit 3): one distinct term per collection in ``ORA_ORDER``,
  ``go_slim`` excluding taken IDs, an empty collection's slot substituted by an untaken term (a
  substitute loop without the exclusion fails ``test_substitute_never_repeats_a_term``).
- ``by_index``, ``slim_obo`` and ``prefetch`` (equal to a sequential run; ``done=``).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

CAL = Path(__file__).resolve().parents[1] / "calibration"
if str(CAL) not in sys.path:
    sys.path.insert(0, str(CAL))

import _criteria as C  # noqa: E402
from helpers import function_gate_chain as F  # noqa: E402

ENV = (C.SHARD_ENV, C.BLOCK_ENV, C.CONFIRM_ENV, "SLURM_CPUS_PER_TASK",
       "PYTEST_XDIST_WORKER_COUNT", "GIPPA_M815_WORKERS")  # fmt: skip


@pytest.fixture(autouse=True)
def _clean_env(monkeypatch: pytest.MonkeyPatch) -> None:
    for v in ENV:
        monkeypatch.delenv(v, raising=False)


def _pick() -> Any:
    from test_function_calibration import pick_one_per_collection

    return pick_one_per_collection


# --------------------------------------------------------------------------- B2, B1


def test_prefetch_indices_unsharded_is_block1() -> None:
    assert F.prefetch_indices(10) == list(range(10))


def test_prefetch_indices_block1_shard(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.SHARD_ENV, "5:10")
    assert F.prefetch_indices(10) == [5, 6, 7, 8, 9]


def test_prefetch_indices_block2_shard_is_offset_by_R(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.SHARD_ENV, "0:5")
    monkeypatch.setenv(C.BLOCK_ENV, "2")
    assert F.prefetch_indices(10) == [10, 11, 12, 13, 14]


def test_prefetch_indices_match_what_calibrate_runs(monkeypatch: pytest.MonkeyPatch,
                                                     tmp_path: Path) -> None:  # fmt: skip
    monkeypatch.setenv(C.SHARD_ENV, "2:4")
    monkeypatch.setenv(C.BLOCK_ENV, "2")
    R, stream, ran = 6, "ora.m815.unit.prefetch", []
    fn = F.by_index(stream, R, lambda i: ran.append(i) or 0.5)
    C.calibrate("ora", "level", 0.05, fn, R, stream=stream, log_dir=tmp_path)
    assert ran == F.prefetch_indices(R) == [8, 9]


def test_prefetch_indices_none_in_a_confirm_run(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(C.CONFIRM_ENV, "[]")
    assert F.prefetch_indices(10) == []


@pytest.mark.parametrize("var", [C.SHARD_ENV, C.CONFIRM_ENV])
def test_sharded(monkeypatch: pytest.MonkeyPatch, var: str) -> None:
    assert not F.sharded()
    monkeypatch.setenv(var, "0:1" if var == C.SHARD_ENV else "[]")
    assert F.sharded()


# --------------------------------------------------------------------------- nit 1


@pytest.mark.parametrize(("cpus", "xdist", "share"), [(None, None, 1), ("16", None, 16),
                                                      ("16", "16", 1), ("16", "4", 4),
                                                      ("4", "16", 1)])  # fmt: skip
def test_cpu_share(monkeypatch: pytest.MonkeyPatch, cpus: str | None, xdist: str | None,
                   share: int) -> None:  # fmt: skip
    if cpus is not None:
        monkeypatch.setenv("SLURM_CPUS_PER_TASK", cpus)
    if xdist is not None:
        monkeypatch.setenv("PYTEST_XDIST_WORKER_COUNT", xdist)
    assert F.cpu_share() == share
    assert F.workers() == min(4, share)


# --------------------------------------------------------------------------- nit 3: picks


def _table(terms: dict[str, list[str]]) -> pd.DataFrame:
    rows = [(c, t, 0.5) for c, ts in terms.items() for t in ts]
    rows.append(("go_bp", "GO:9000999", float("nan")))  # untested: never picked
    return pd.DataFrame(rows, columns=["collection", "term", "p_value"])


FULL = {"cog": ["A", "B"], "kegg": ["map90000", "map90001"], "go_bp": ["GO:9000010"],
        "go_mf": ["GO:9000020"], "go_cc": ["GO:9000030"],
        "go_slim": ["GO:9000010", "GO:9000020", "GO:9000040"]}  # fmt: skip


@pytest.mark.parametrize("seed", range(20))
def test_one_term_per_collection_in_order(seed: int) -> None:
    from test_function_calibration import ORA_ORDER

    picks, short = _pick()(_table(FULL), np.random.default_rng(seed))
    assert [c for c, _ in picks] == list(ORA_ORDER) and short == 0
    assert picks[-1] == ("go_slim", "GO:9000040")  # slim excludes the picked GO IDs
    assert len({t for _, t in picks}) == 6


def test_pick_is_reproducible() -> None:
    a = _pick()(_table(FULL), np.random.default_rng(3))
    b = _pick()(_table(FULL), np.random.default_rng(3))
    assert a == b


@pytest.mark.parametrize("seed", range(40))
def test_substitute_never_repeats_a_term(seed: int) -> None:
    """KEGG empty: its slot is filled after the loop, when exactly one term ID is still untaken,
    so a substitute loop that ignored the taken IDs would repeat a pick in most seeds."""
    terms = {**FULL, "kegg": [], "go_bp": ["GO:9000010", "GO:9000011"],
             "go_slim": ["GO:9000010", "GO:9000020", "GO:9000030", "GO:9000040"],
             "cog": ["A"]}  # fmt: skip
    picks, short = _pick()(_table(terms), np.random.default_rng(seed))
    ids = [t for _, t in picks]
    assert short == 1 and len(picks) == 6
    assert len(set(ids)) == 6, picks
    assert "GO:9000999" not in ids


def test_exhausted_pool_fails_closed() -> None:
    terms = {"cog": ["A"], "kegg": [], "go_bp": ["GO:9000010"], "go_mf": ["GO:9000020"],
             "go_cc": ["GO:9000030"], "go_slim": ["GO:9000040"]}  # fmt: skip
    with pytest.raises(AssertionError):
        _pick()(_table(terms), np.random.default_rng(0))


# --------------------------------------------------------------------------- nit 3: by_index


def test_by_index_maps_both_blocks_and_fails_closed() -> None:
    R, stream = 4, "ora.m815.unit.by_index"
    fn = F.by_index(stream, R, lambda i: i)
    assert [fn(C.replicate_rng(stream, i)) for i in range(2 * R)] == list(range(2 * R))
    with pytest.raises(KeyError):
        fn(C.replicate_rng(stream, 2 * R))
    with pytest.raises(KeyError):
        fn(C.replicate_rng(stream + ".other", 0))


# --------------------------------------------------------------------------- nit 3: slim_obo

OBO = """format-version: 1.2
ontology: go

[Term]
id: GO:9000000
name: root bp
namespace: biological_process

""" + "\n".join(
    f"[Term]\nid: GO:{9000001 + k}\nname: t{k}\nnamespace: biological_process\n"
    f"is_a: GO:9000000\n{'relationship: part_of GO:9000000' if k % 2 else ''}\n"
    for k in range(9)
)


def test_slim_obo_every_third_non_root_term_bare(tmp_path: Path) -> None:
    from gippa.io.function_resources import read_go

    slim = F.slim_obo(OBO)
    ids = re.findall(r"^id: (\S+)", slim, re.M)
    assert ids == ["GO:9000001", "GO:9000004", "GO:9000007"]  # roots excluded, every 3rd by ID
    assert "is_a:" not in slim and "relationship:" not in slim
    path = tmp_path / "slim.obo"
    path.write_text(slim)
    dag = read_go(path)
    assert sorted(dag.names) == ids
    assert set(dag.namespace.values()) == {"biological_process"}


def test_slim_obo_of_a_world_parses_and_stays_reserved(tmp_path: Path) -> None:
    import gippasim.function_sim as fs
    import gippasim.scenarios as sc
    from gippa.io.function_resources import read_go

    w = fs.draw(sc.ora_null(400), np.random.default_rng(0))
    path = tmp_path / "slim.obo"
    path.write_text(F.slim_obo(w.obo))
    names = read_go(path).names
    assert names and all(9000000 <= int(t.removeprefix("GO:")) <= 9999998 for t in names)


# --------------------------------------------------------------------------- prefetch


def square(k: int) -> dict[str, int]:
    return {"v": k * k, "full": 1}


def test_prefetch_equals_sequential_and_recomputes_incomplete(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "2")
    cache: dict[Any, Any] = {3: {"v": -1}}  # incomplete entry: no "full"
    F.prefetch(cache, square, list(range(6)), done=lambda v: "full" in v)
    assert cache == {k: square(k) for k in range(6)}
    cache[7] = {"v": -7}
    F.prefetch(cache, square, [7])  # no done=: a cached key is kept
    assert cache[7] == {"v": -7}
