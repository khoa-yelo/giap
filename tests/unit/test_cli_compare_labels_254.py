"""``giap compare`` run labels that would break V16 (PR #254 review B4; M9-23 lead decision 5).

The ``compare_detection`` intersection names join labels with ``&`` and V16's ``runs`` column
joins them with ``;``; ``none`` and ``other`` are the reserved intersection rows. A label holding
``&`` or ``;``, or equal to ``none`` or ``other``, is rejected by ``giap.io.compare_io._labels``
with a ``CompareInputError`` naming the label, so ``giap compare`` exits 1 before any table is
written (the M9-22 rule "no table is written when the command exits 1").
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from click.testing import CliRunner, Result
from helpers.compare_data import gene, write_run

from giap.exitcodes import ExitCode

TABLES = ("compare_detection", "compare_abundance", "compare_stability")
N = 12
BAD = ["a&b", "a;b", "none", "other"]


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "db_empty"))


def _genes(shift: float) -> list[dict[str, object]]:
    vals = np.sin(np.arange(N) * 1.7) * 3 + np.arange(N) * 0.4
    return [gene(f"g{i:02d}", start=1000 * i, log2=float(v + shift * (i % 3)), tu=f"TU_{i // 2}")
            for i, v in enumerate(vals)]  # fmt: skip


def _compare(*args: str) -> Result:
    from giap.cli import cli

    res = CliRunner().invoke(cli, ["compare", *args, "--set", "stats.bootstrap=50"])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _nothing_written(cmp: Path) -> bool:
    return (
        not any((cmp / "tables" / f"{t}.tsv").exists() for t in TABLES)
        and not (cmp / "report" / "compare.html").exists()
    )


@pytest.mark.parametrize("bad", BAD)
def test_labels_function_rejects_reserved_and_separator_labels(bad: str) -> None:
    from giap.io.compare_io import CompareInputError, _labels

    with pytest.raises(CompareInputError) as err:
        _labels([Path("runA"), Path("runB")], f"{bad},ok")
    assert bad in str(err.value), str(err.value)


@pytest.mark.parametrize("bad", BAD)
def test_labels_function_rejects_them_in_any_position(bad: str) -> None:
    from giap.io.compare_io import CompareInputError, _labels

    with pytest.raises(CompareInputError):
        _labels([Path("runA"), Path("runB"), Path("runC")], f"x,y,{bad}")


@pytest.mark.parametrize("good", ["A1", "run-1", "a.b", "a_b", "None2", "others"])
def test_ordinary_labels_still_pass(good: str) -> None:
    from giap.io.compare_io import _labels

    assert _labels([Path("runA"), Path("runB")], f"{good},z") == [good, "z"]


@pytest.mark.parametrize("bad", BAD)
def test_cli_exits_1_before_any_table_is_written(tmp_path: Path, bad: str) -> None:
    a = write_run(tmp_path / "runA", _genes(0.0))
    b = write_run(tmp_path / "runB", _genes(0.5))
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp), "--labels", f"{bad},c")
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert bad in res.output, res.output
    assert _nothing_written(cmp)


def test_cli_rejects_a_default_label_with_an_ampersand(tmp_path: Path) -> None:
    """Default labels are the directory basenames; ``_labels`` returns them, so they are checked
    too (a run directory named ``wt&ko`` would crash V16 just the same)."""
    a = write_run(tmp_path / "wt&ko", _genes(0.0))
    b = write_run(tmp_path / "runB", _genes(0.5))
    cmp = tmp_path / "CMP"
    res = _compare(str(a), str(b), "--outdir", str(cmp))
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert "wt&ko" in res.output
    assert _nothing_written(cmp)
