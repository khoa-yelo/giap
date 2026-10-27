"""T2: ``gippa.io.fgsea_r`` file I/O and the recorded-output parse (M8-17a §3, §4, §6; A4).

- ``write_call`` writes ``stats.tsv``, ``jobs.tsv`` and ``sets.tsv`` byte for byte (UTF-8, tab,
  header, ``\\n``; scores with ``repr(float)``; jobs and sets in call order, members long form).
- ``read_results`` parses ``results.tsv`` (``NA`` → NaN, empty leading edge → ()) and rejects
  every malformed output with ``FgseaEngineError``.
- ``read_versions`` reads ``versions.txt`` and rejects any fgsea other than 1.36.2.
- Every recorded directory parses, and ``write_call`` of the rebuilt call reproduces its inputs.
"""

from __future__ import annotations

import math
import struct
from pathlib import Path

import numpy as np
import pytest
from helpers.fgsea_r_data import RECORDED, RESULTS_HEADER, load_call
from hypothesis import given, settings
from hypothesis import strategies as st

BY_HAND = RECORDED / "by_hand"
REQUIRED_INPUTS = ("sim_main", "sim_retest", "ties", "unbalanced", "tiny_p")
RECORDING_FILES = (
    "stats.tsv",
    "jobs.tsv",
    "sets.tsv",
    "command.txt",
    "results.tsv",
    "versions.txt",
    "PROVENANCE.md",
)
SETTINGS = settings(max_examples=150, deadline=None)


def _call(
    tu_ids: tuple[str, ...],
    scores: list[float],
    jobs: list[tuple[str, dict[str, tuple[str, ...]], int, int]],
    seed: int = 12345,
):
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import EngineCall, GseaJob

    return EngineCall(
        tu_ids=tu_ids,
        scores=np.array(scores, dtype=np.float64),
        jobs=tuple(GseaJob(job_id=j, sets=s, min_size=lo, max_size=hi) for j, s, lo, hi in jobs),
        params=GseaConfig(),
        seed=seed,
    )


def _small_call():
    return _call(
        ("TU_3", "TU_1", "TU_2"),
        [0.1 + 0.2, -1.5, 2.5e-300],
        [
            ("main:cog", {"C1": ("TU_1", "TU_2"), "C2": ("TU_3",)}, 1, 3),
            ("retest:go_bp:GO:0008150", {"GO:0008150": ("TU_2", "TU_3")}, 1, 3),
        ],
    )


# ---------------------------------------------------------------- write_call


def test_write_call_byte_exact(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import write_call

    write_call(_small_call(), tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["jobs.tsv", "sets.tsv", "stats.tsv"]
    assert (tmp_path / "stats.tsv").read_bytes() == (
        b"tu_id\tscore\nTU_3\t0.30000000000000004\nTU_1\t-1.5\nTU_2\t2.5e-300\n"
    )
    assert (tmp_path / "jobs.tsv").read_bytes() == (
        b"job_id\tmin_size\tmax_size\nmain:cog\t1\t3\nretest:go_bp:GO:0008150\t1\t3\n"
    )
    assert (tmp_path / "sets.tsv").read_bytes() == (
        b"job_id\tset_id\tmember\n"
        b"main:cog\tC1\tTU_1\n"
        b"main:cog\tC1\tTU_2\n"
        b"main:cog\tC2\tTU_3\n"
        b"retest:go_bp:GO:0008150\tGO:0008150\tTU_2\n"
        b"retest:go_bp:GO:0008150\tGO:0008150\tTU_3\n"
    )


def _bits(x: float) -> bytes:
    return struct.pack("<d", x)


@SETTINGS
@given(
    st.lists(
        st.floats(allow_nan=False, allow_infinity=False, allow_subnormal=True),
        min_size=1,
        max_size=40,
    )
)
def test_write_call_scores_round_trip(
    tmp_path_factory: pytest.TempPathFactory, scores: list[float]
) -> None:
    """Every finite double (incl. -0.0 and subnormals) is read back bit for bit, in call order."""
    from gippa.io.fgsea_r import write_call

    d = tmp_path_factory.mktemp("rt")
    ids = tuple(f"TU_{i:04d}" for i in reversed(range(len(scores))))
    write_call(_call(ids, scores, [("main:x", {"S": (ids[0],)}, 1, len(ids))]), d)
    lines = (d / "stats.tsv").read_text(encoding="utf-8").split("\n")
    assert lines[0] == "tu_id\tscore" and lines[-1] == ""
    rows = [line.split("\t") for line in lines[1:-1]]
    assert tuple(r[0] for r in rows) == ids
    assert [_bits(float(r[1])) for r in rows] == [_bits(s) for s in scores]
    assert [r[1] for r in rows] == [repr(float(s)) for s in scores]


def test_write_call_rejects_non_finite_scores(tmp_path: Path) -> None:
    """NA scores never reach R (spec §8 edge cases): ValueError before any file or process."""
    from gippa.io.fgsea_r import write_call

    with pytest.raises(ValueError):
        call = _call(("A", "B"), [1.0, math.nan], [("main:x", {"S": ("A",)}, 1, 2)])
        write_call(call, tmp_path)
    assert not (tmp_path / "stats.tsv").exists()


# ---------------------------------------------------------------- read_results: by_hand


def test_read_results_by_hand_values() -> None:
    from gippa.io.fgsea_r import read_results

    result = read_results(BY_HAND / "results.tsv", load_call(BY_HAND))
    assert set(result) == {"main:cog", "main:kegg"}
    cog = result["main:cog"]
    assert [r.term for r in cog] == ["C0001", "C0002", "C0003"]
    c1 = cog[0]
    assert c1.size == 4
    assert c1.es == float("0.86956521739130432")
    assert c1.nes == float("1.6340718861470217")
    assert c1.p == float("0.0089820359281437123")
    assert c1.log2err == float("0.38073040072279097")
    assert c1.leading_edge == ("TU_07", "TU_03", "TU_01")
    assert cog[1].es < 0 and cog[1].leading_edge == ("TU_06", "TU_04")

    kegg = result["main:kegg"]
    assert [r.term for r in kegg] == ["K00001", "K00002", "K00003"]
    unbalanced = kegg[0]  # example 3: pval, NES, log2err NA; ES kept
    assert unbalanced.es == float("0.81000000000000005")
    assert math.isnan(unbalanced.p) and math.isnan(unbalanced.nes)
    assert math.isnan(unbalanced.log2err)
    no_err = kegg[1]  # example 4: finite tiny p, log2err NA
    assert no_err.p == float("3.2000000000000001e-15")
    assert math.isnan(no_err.log2err)
    assert no_err.leading_edge == ("TU_06", "TU_10", "TU_04")
    empty = kegg[2]
    assert empty.es == 0.0 and empty.p == 1.0 and empty.log2err == 0.0
    assert empty.leading_edge == ()


def test_read_versions_by_hand() -> None:
    from gippa.io.fgsea_r import read_versions

    prov = read_versions(BY_HAND / "versions.txt")
    assert prov.fgsea == "1.36.2"
    assert prov.r == "R version 4.5.3 (2026-03-11)"
    assert len(prov.warnings) == 3
    flat = [" ".join(map(str, w)) if not isinstance(w, str) else w for w in prov.warnings]
    assert sum("main:kegg" in w for w in flat) == 2
    assert any("unbalanced" in w for w in flat)
    assert any("ties" in w for w in flat)


@pytest.mark.parametrize("version", ["1.38.0", "1.36.1", ""])
def test_read_versions_rejects_other_fgsea(tmp_path: Path, version: str) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_versions

    text = (BY_HAND / "versions.txt").read_text(encoding="utf-8")
    path = tmp_path / "versions.txt"
    path.write_text(text.replace("fgsea\t1.36.2", f"fgsea\t{version}"), encoding="utf-8")
    with pytest.raises(FgseaEngineError, match=r"fgsea"):
        read_versions(path)


def test_read_versions_rejects_missing_fgsea_line(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_versions

    text = (BY_HAND / "versions.txt").read_text(encoding="utf-8")
    path = tmp_path / "versions.txt"
    path.write_text("\n".join(x for x in text.split("\n") if not x.startswith("fgsea\t")))
    with pytest.raises(FgseaEngineError):
        read_versions(path)


# ---------------------------------------------------------------- read_results: mutations


def _mutate(tmp_path: Path, old: str, new: str, count: int = 1) -> Path:
    text = (BY_HAND / "results.tsv").read_text(encoding="utf-8")
    assert old in text, old
    path = tmp_path / "results.tsv"
    path.write_text(text.replace(old, new, count), encoding="utf-8")
    return path


ROW_C1 = (
    "main:cog\tC0001\t4\t0.86956521739130432\t1.6340718861470217\t0.0089820359281437123\t"
    "0.38073040072279097\tTU_07;TU_03;TU_01"
)
ROW_C2 = (
    "main:cog\tC0002\t3\t-0.81818181818181823\t-1.5079365079365079\t0.032258064516129031\t"
    "0.32177591807536657\tTU_06;TU_04"
)
MUTATIONS = {
    "wrong header": (RESULTS_HEADER, RESULTS_HEADER.replace("ES", "es", 1)),
    "reordered header": (RESULTS_HEADER, RESULTS_HEADER.replace("ES\tNES", "NES\tES")),
    "unknown job": ("main:kegg\tK00002", "main:bogus\tK00002"),
    "unknown set": ("main:cog\tC0003", "main:cog\tC9999"),
    "set from another job": ("main:cog\tC0003", "main:cog\tK00002"),
    "size mismatch": ("main:cog\tC0001\t4", "main:cog\tC0001\t5"),
    "non-integer size": ("main:cog\tC0001\t4", "main:cog\tC0001\t4.0"),
    "edge outside set": ("TU_07;TU_03;TU_01", "TU_07;TU_03;TU_02"),
    "edge not ranked": ("TU_06;TU_04", "TU_06;TU_99"),
    "non-numeric ES": ("\t0.86956521739130432\t", "\tabc\t"),
    "non-numeric pval": ("\t0.032258064516129031\t", "\tnull\t"),
    "short row": ("\t0.38073040072279097\tTU_07;TU_03;TU_01", "\t0.38073040072279097"),
    "out of order": (ROW_C1 + "\n" + ROW_C2, ROW_C2 + "\n" + ROW_C1),
    "duplicate row": (ROW_C1, ROW_C1 + "\n" + ROW_C1),
}


@pytest.mark.parametrize("name", sorted(MUTATIONS))
def test_read_results_rejects_mutation(tmp_path: Path, name: str) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    old, new = MUTATIONS[name]
    path = _mutate(tmp_path, old, new)
    with pytest.raises(FgseaEngineError):
        read_results(path, load_call(BY_HAND))


def test_read_results_rejects_missing_file(tmp_path: Path) -> None:
    from gippa.io.fgsea_r import FgseaEngineError, read_results

    with pytest.raises(FgseaEngineError):
        read_results(tmp_path / "results.tsv", load_call(BY_HAND))


# ---------------------------------------------------------------- read_results: property


def _r_num(x: float) -> str:
    """R's ``format(x, digits = 17)`` up to spelling: 17 significant digits, or ``NA``."""
    return "NA" if math.isnan(x) else f"{x:.17g}"


_num = st.one_of(st.just(math.nan), st.floats(allow_nan=False, allow_infinity=False, width=64))


@SETTINGS
@given(
    st.lists(st.tuples(_num, _num, _num, _num, st.integers(0, 4)), min_size=1, max_size=6),
    st.data(),
)
def test_read_results_round_trip(
    tmp_path_factory: pytest.TempPathFactory, rows: list, data: st.DataObject
) -> None:
    from gippa.io.fgsea_r import read_results

    ids = tuple(f"TU_{i}" for i in range(8))
    sets = {f"S{k}": ids[k : k + 4] for k in range(len(rows))}
    call = _call(ids, [float(8 - i) for i in range(8)], [("main:x", sets, 1, 8)])
    lines = [RESULTS_HEADER]
    edges = []
    for (term, members), (es, nes, p, err, n_edge) in zip(sets.items(), rows, strict=True):
        edge = tuple(data.draw(st.permutations(members)))[:n_edge]
        edges.append(edge)
        lines.append(
            "\t".join(
                [
                    "main:x",
                    term,
                    str(len(members)),
                    _r_num(es),
                    _r_num(nes),
                    _r_num(p),
                    _r_num(err),
                    ";".join(edge),
                ]
            )
        )
    d = tmp_path_factory.mktemp("rr")
    (d / "results.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    got = read_results(d / "results.tsv", call)["main:x"]
    assert [r.term for r in got] == list(sets)
    for r, (es, nes, p, err, _), edge in zip(got, rows, edges, strict=True):
        for have, want in ((r.es, es), (r.nes, nes), (r.p, p), (r.log2err, err)):
            assert (math.isnan(have) and math.isnan(want)) or _bits(have) == _bits(want)
        assert r.leading_edge == edge
        assert r.size == len(sets[r.term])


# ---------------------------------------------------------------- recorded directories (§11.4)


def _recorded_dirs() -> list[Path]:
    return sorted(p for p in RECORDED.iterdir() if p.is_dir())


def test_recorded_readme_marks_by_hand() -> None:
    text = (RECORDED / "README.md").read_text(encoding="utf-8")
    assert "by_hand" in text and "recorded-by-hand" in text


@pytest.mark.parametrize("name", REQUIRED_INPUTS)
def test_required_recording_present(name: str) -> None:
    d = RECORDED / name
    for f in RECORDING_FILES:
        assert (d / f).is_file(), f"{name}/{f}"
        assert (d / f).stat().st_size < 1024**2, f"{name}/{f} is over 1 MB"
    prov = (d / "PROVENANCE.md").read_text(encoding="utf-8")
    assert "recorded-by-hand" not in prov
    assert "sha256" in prov
    # recorded before ADR-0064 under the wrapper's old name: the command as PROVENANCE.md has it
    cmd = (d / "command.txt").read_text(encoding="utf-8").strip()
    assert cmd.split()[0].endswith("_fgsea.R") and f"`{cmd}`" in prov


@pytest.mark.parametrize("d", _recorded_dirs(), ids=lambda p: p.name)
def test_recorded_dir_parses_and_round_trips(tmp_path: Path, d: Path) -> None:
    """T7 precondition: ``write_call(load_call(d))`` gives the recorded inputs byte for byte."""
    from gippa.io.fgsea_r import FGSEA_VERSION, read_results, read_versions, write_call

    call = load_call(d)
    write_call(call, tmp_path)
    for f in ("stats.tsv", "jobs.tsv", "sets.tsv"):
        assert (tmp_path / f).read_bytes() == (d / f).read_bytes(), f
    result = read_results(d / "results.tsv", call)
    assert set(result) <= {j.job_id for j in call.jobs}
    assert read_versions(d / "versions.txt").fgsea == FGSEA_VERSION


def _rows_of(name: str):
    from gippa.io.fgsea_r import read_results

    d = RECORDED / name
    result = read_results(d / "results.tsv", load_call(d))
    return [r for rows in result.values() for r in rows]


def test_recorded_unbalanced_has_na_pval() -> None:
    rows = _rows_of("unbalanced")
    assert any(math.isnan(r.p) and math.isnan(r.nes) and not math.isnan(r.es) for r in rows)


def test_recorded_tiny_p_has_finite_log2err() -> None:
    rows = _rows_of("tiny_p")
    assert any(r.p < 1e-20 and math.isfinite(r.log2err) for r in rows)


def test_recorded_sim_main_and_retest_job_ids() -> None:
    main = load_call(RECORDED / "sim_main")
    ids = [j.job_id for j in main.jobs]
    assert ids and all(i.startswith("main:") for i in ids)
    assert {i.split(":", 1)[1] for i in ids} <= {"cog", "kegg", "go_bp", "go_slim"}
    retest = load_call(RECORDED / "sim_retest")
    assert retest.jobs and all(j.job_id.startswith("retest:") for j in retest.jobs)
    assert all(len(j.sets) == 1 and j.min_size == 1 for j in retest.jobs)
    assert all(j.max_size == len(retest.tu_ids) for j in retest.jobs)
    assert retest.tu_ids == main.tu_ids and retest.seed == main.seed
