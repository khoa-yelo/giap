"""M8-17b: the #204 review nits 1 and 2 on ``gippa.io.fgsea_r`` (docs/backlog.md "M8-17a review
nits"; docs/specs/M8-17-r-fgsea-engine.md §3 "Provenance", §4 step 6; A4).

NIT 1. ``read_results`` requires each job's rows to be **exactly** its in-range sets
(``min_size <= #members <= max_size``) in the sent order. Python sends only in-range sets for a
main job and a ``[1, n]`` single set for a re-test (§8), and fgsea 1.36.2 returns one row per
in-range set (NA statistics included), so a missing row means the output was truncated or the
names did not match: ``FgseaEngineError``. A set outside its job's range still adds no row.

NIT 2. ``versions.txt`` must carry the ``seed`` and ``jobs`` lines of §3 and the pinned
``data.table`` 1.18.6.1: ``read_versions`` rejects a file without them, and the engine rejects a
``seed`` other than the call's or a ``jobs`` count other than the number of sent jobs.

No R: the recorded outputs and the fake R front end of ``helpers.fgsea_r_data``.
"""

from __future__ import annotations

import json
import stat
import sys
from pathlib import Path

import numpy as np
import pytest
from helpers.fgsea_r_data import R_FRONT_END, RECORDED, RESULTS_HEADER, fake_prefix, load_call

BY_HAND = RECORDED / "by_hand"
SIM_MAIN = RECORDED / "sim_main"


def _lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").split("\n")[:-1]


def _write(path: Path, lines: list[str]) -> Path:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------- NIT 1


def test_all_in_range_rows_present_is_fine() -> None:
    from gippa.io.fgsea_r import read_results

    res = read_results(SIM_MAIN / "results.tsv", load_call(SIM_MAIN))
    assert sum(len(v) for v in res.values()) == 32


@pytest.mark.parametrize("drop", [1, 10, 32])
def test_a_dropped_row_is_an_error(tmp_path: Path, drop: int) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    lines = _lines(SIM_MAIN / "results.tsv")
    path = _write(tmp_path / "results.tsv", lines[:drop] + lines[drop + 1 :])
    with pytest.raises(FgseaEngineError):
        read_results(path, load_call(SIM_MAIN))


def test_a_job_without_rows_is_an_error(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    lines = _lines(SIM_MAIN / "results.tsv")
    path = _write(tmp_path / "results.tsv", [x for x in lines if not x.startswith("main:kegg\t")])
    with pytest.raises(FgseaEngineError):
        read_results(path, load_call(SIM_MAIN))


def test_header_only_with_in_range_sets_is_an_error(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    path = _write(tmp_path / "results.tsv", [RESULTS_HEADER])
    with pytest.raises(FgseaEngineError):
        read_results(path, load_call(BY_HAND))


def test_a_dropped_retest_row_is_an_error(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    d = RECORDED / "sim_retest"
    lines = _lines(d / "results.tsv")
    path = _write(tmp_path / "results.tsv", lines[:-1])
    with pytest.raises(FgseaEngineError):
        read_results(path, load_call(d))


def _call_with_out_of_range_set():
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import EngineCall, GseaJob

    ids = tuple(f"TU_{i}" for i in range(1, 9))
    sets = {"IN": ("TU_1", "TU_2", "TU_3"), "SMALL": ("TU_4",), "IN2": ("TU_5", "TU_6", "TU_7")}
    return EngineCall(
        ids, np.linspace(2.0, -2.0, 8), (GseaJob("main:cog", sets, 2, 7),), GseaConfig(), 99
    )


def _row(term: str, size: int, edge: str) -> str:
    return f"main:cog\t{term}\t{size}\t0.5\t1.2\t0.03\t0.2\t{edge}"


def test_out_of_range_set_without_row_is_fine(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import read_results

    path = _write(
        tmp_path / "results.tsv", [RESULTS_HEADER, _row("IN", 3, "TU_1"), _row("IN2", 3, "TU_5")]
    )
    res = read_results(path, _call_with_out_of_range_set())
    assert [r.term for r in res["main:cog"]] == ["IN", "IN2"]


def test_in_range_set_missing_next_to_out_of_range_set_is_an_error(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    path = _write(tmp_path / "results.tsv", [RESULTS_HEADER, _row("IN", 3, "TU_1")])
    with pytest.raises(FgseaEngineError):
        read_results(path, _call_with_out_of_range_set())


# ---------------------------------------------------------------- NIT 2: read_versions


def _versions(tmp_path: Path, drop: str | None = None, edit: tuple[str, str] | None = None) -> Path:
    lines = _lines(BY_HAND / "versions.txt")
    if drop is not None:
        lines = [x for x in lines if not x.startswith(f"{drop}\t")]
    if edit is not None:
        lines = [x.replace(*edit) for x in lines]
    return _write(tmp_path / "versions.txt", lines)


@pytest.mark.parametrize("key", ["seed", "jobs", "data.table"])
def test_read_versions_requires_the_line(tmp_path: Path, key: str) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_versions

    with pytest.raises(FgseaEngineError):
        read_versions(_versions(tmp_path, drop=key))


@pytest.mark.parametrize(
    "edit",
    [
        ("seed\t1234567", "seed\tabc"),
        ("jobs\t2", "jobs\ttwo"),
        ("seed\t1234567", "seed\t"),
        ("jobs\t2", "jobs\t-1"),
    ],
)
def test_read_versions_rejects_non_integer_seed_or_jobs(
    tmp_path: Path, edit: tuple[str, str]
) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_versions

    with pytest.raises(FgseaEngineError):
        read_versions(_versions(tmp_path, edit=edit))


@pytest.mark.parametrize("version", ["1.17.0", "1.18.6", ""])
def test_read_versions_rejects_other_data_table(tmp_path: Path, version: str) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_versions

    with pytest.raises(FgseaEngineError, match=r"data\.table"):
        read_versions(_versions(tmp_path, edit=("data.table\t1.18.6.1", f"data.table\t{version}")))


def test_read_versions_accepts_every_recording() -> None:
    from gippa.io.fgsea_r import read_versions

    for d in sorted(p for p in RECORDED.iterdir() if p.is_dir()):
        assert read_versions(d / "versions.txt").fgsea == "1.36.2", d.name


# ---------------------------------------------------------------- NIT 2: the engine


_EDITOR = r'''#!{python} -I
"""Runs the fake R front end, then rewrites versions.txt (tests only)."""
import json, subprocess, sys
from pathlib import Path

prefix = Path({prefix!r})
code = subprocess.call([{real!r}, *sys.argv[1:]])
edit = json.loads((prefix / "edit.json").read_text())
argv = sys.argv[1:]
opts = dict(zip(argv[2::2], argv[3::2]))
v = Path.cwd() / opts["--versions"]
if v.exists():
    lines = v.read_text().split("\n")[:-1]
    out = []
    for line in lines:
        key = line.split("\t", 1)[0]
        if key in edit:
            if edit[key] is None:
                continue
            line = key + "\t" + str(edit[key])
        out.append(line)
    v.write_text("\n".join(out) + "\n")
sys.exit(code)
'''


def _editing_prefix(root: Path, edit: dict[str, object]) -> Path:
    prefix = fake_prefix(root)
    exe = prefix / "bin" / R_FRONT_END
    real = prefix / "bin" / "fake_r_real"
    exe.rename(real)
    exe.write_text(_EDITOR.format(python=sys.executable, prefix=str(prefix), real=str(real)))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    (prefix / "edit.json").write_text(json.dumps(edit))
    return prefix


def _call():
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import EngineCall, GseaJob

    ids = tuple(f"TU_{i}" for i in (5, 2, 7, 1, 4, 6, 3, 8))
    return EngineCall(
        tu_ids=ids,
        scores=np.array([3.0, 2.5, 1.0, 0.5, -0.25, -1.0, -2.0, -3.5]),
        jobs=(
            GseaJob("main:cog", {"C1": ("TU_1", "TU_2", "TU_3"), "C2": ("TU_4", "TU_5")}, 2, 7),
            GseaJob("main:kegg", {"K1": ("TU_6", "TU_7", "TU_8")}, 2, 7),
        ),
        params=GseaConfig(),
        seed=424242,
    )


def _run(prefix: Path, work: Path):
    from gippa.io.fgsea_r import r_fgsea_engine, resolve_launcher

    eng = r_fgsea_engine(
        resolve_launcher({"GIPPA_FGSEA_PREFIX": str(prefix)}), timeout_s=60, tmp_root=work
    )
    return eng(_call())


@pytest.fixture
def work(tmp_path: Path) -> Path:
    d = tmp_path / "work"
    d.mkdir()
    return d


def test_editor_fake_passes_through_unchanged(tmp_path: Path, work: Path) -> None:
    res = _run(_editing_prefix(tmp_path / "env", {}), work)
    assert set(res) == {"main:cog", "main:kegg"}


@pytest.mark.parametrize(
    "edit",
    [
        {"seed": 424243},
        {"seed": None},
        {"jobs": 3},
        {"jobs": 1},
        {"jobs": None},
        {"data.table": "1.17.0"},
    ],
    ids=["seed_other", "seed_missing", "jobs_more", "jobs_fewer", "jobs_missing", "data_table"],
)
def test_engine_rejects_versions_not_matching_the_call(
    tmp_path: Path, work: Path, edit: dict[str, object]
) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, FgseaEngineUnavailable

    with pytest.raises(FgseaEngineError) as info:
        _run(_editing_prefix(tmp_path / "env", edit), work)
    assert not isinstance(info.value, FgseaEngineUnavailable)
    assert list(work.iterdir()) == []
