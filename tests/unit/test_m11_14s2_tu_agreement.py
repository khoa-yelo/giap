"""``scripts/tu_agreement.py``: TU pair agreement with RegulonDB (M11-14 14s2, weekly arm).

Synthetic truth and ``operons.tsv`` with hand-computed counts (#400 review B1): ambiguous and
conflicting pairs are dropped, pair orientation does not matter, a pair with one gene outside
every TU is not called, a pair across the origin is an ordinary pair, nothing called gives NaN
precision, and the printed line carries the counts.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from types import ModuleType

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "tu_agreement.py"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("tu_agreement", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ta = _load()

# g1..g6 on a circular replicon; g6-g1 spans the origin.
TRUTH = [
    {"gene_a": "g1", "gene_b": "g2", "label": "positive"},  # in TU A: TP
    {"gene_a": "g3", "gene_b": "g2", "label": "negative"},  # reversed, g2 in A, g3 in B: TN
    {"gene_a": "g3", "gene_b": "g4", "label": "positive"},  # both in B: TP
    {"gene_a": "g4", "gene_b": "g5", "label": "positive"},  # g5 in no TU: FN
    {"gene_a": "g6", "gene_b": "g1", "label": "negative"},  # across the origin, both in A: FP
    {"gene_a": "g5", "gene_b": "g6", "label": "ambiguous"},  # dropped
    {"gene_a": "g7", "gene_b": "g8", "label": "positive"},  # labelled twice, both ways: dropped
    {"gene_a": "g8", "gene_b": "g7", "label": "negative"},
]
OPERONS = [
    {"tu_id": "A", "genes": "g6;g1;g2"},
    {"tu_id": "B", "genes": "g3,g4"},
    {"tu_id": "C", "genes": "g7;g8"},
]


def test_truth_drops_ambiguous_and_conflicting_pairs() -> None:
    truth = ta.truth_labels(TRUTH)
    assert len(truth) == 5
    assert frozenset(("g5", "g6")) not in truth
    assert frozenset(("g7", "g8")) not in truth
    assert truth[frozenset(("g2", "g3"))] is False


def test_hand_computed_counts() -> None:
    a = ta.agreement(ta.truth_labels(TRUTH), ta.tu_of(OPERONS))
    assert (a["tp"], a["fp"], a["fn"]) == (2, 1, 1)
    assert a["precision"] == pytest.approx(2 / 3)
    assert a["recall"] == pytest.approx(2 / 3)


def test_nothing_called_gives_nan_precision() -> None:
    a = ta.agreement(ta.truth_labels(TRUTH), {})
    assert (a["tp"], a["fp"], a["fn"]) == (0, 0, 3)
    assert math.isnan(a["precision"])
    assert a["recall"] == 0


def test_na_genes_and_missing_tu_id() -> None:
    assert "NA" not in ta.tu_of([{"tu_id": "A", "genes": "g1;NA"}])
    with pytest.raises(ValueError):
        ta.tu_of([{"genes": "g1;g2"}])


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    path.write_text("\t".join(cols) + "\n" + "".join("\t".join(r[c] for c in cols) + "\n"
                                                      for r in rows))  # fmt: skip


def test_main_prints_the_counts(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    reg = tmp_path / "regulondb" / "14.5.0"
    _write(reg / "truth" / "pair_truth.tsv", TRUTH)
    _write(tmp_path / "fb" / "tables" / "operons.tsv", OPERONS)
    bad = tmp_path / "bad" / "tables" / "operons.tsv"
    _write(bad, [{"genes": "g1;g2", "strand": "+"}])
    arms = {"fallback": "fb", "uniop": "none", "broken": "bad"}
    argv = ["--regulondb-dir", str(reg)]
    for name, sub in arms.items():
        argv += ["--arm", f"{name}={tmp_path / sub}"]
    code = ta.main(argv)
    out = capsys.readouterr().out.strip()
    assert code == 0
    assert out.startswith("TU pairs vs RegulonDB 14.5.0 (5 labelled): "), out
    assert "fallback P=0.667 R=0.667 (TP 2 FP 1 FN 1)" in out, out
    assert "uniop: not available" in out, out
    assert "broken: not available" in out, out


def test_arm_without_equals_is_refused(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc:
        ta.main(["--regulondb-dir", str(tmp_path), "--arm", "fallback"])
    assert exc.value.code == 2
