"""``run_evidence.py`` for A3 (rules A10, A11), rule set C for N labels (C7), and the 4-argument
``pilot_compare.sbatch`` (M10-10a criteria 5-7).

Spec: docs/specs/M10-10-a3-run.md, "The A3 run and its evidence", "The three-way compare",
examples 7 and 8. Run trees, ``sacct`` rows and the compare-job sandbox come from M10-04's kit
(``helpers.pilot_runs_m10_04``, PR #296); the compare tables are the output of the real
``giap compare`` on three tiny arms. Interface conventions (``expect["a2"]`` / ``--a2``,
``pair_cmp`` / ``--pair-cmp``, ``expected_warnings``, ``warnings``):
``tests/helpers/a3_run_m10_10.py``.
"""

from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_run_m10_10 import (
    REFERENCE_ROLE,
    a3_cli_args,
    a3_expect,
    giap_compare_n,
    m10_04,
    three_arms,
    write_a3_arm,
    write_mapping_fdr,
)

LABELS3 = ("A1", "A2", "A3")
PAIRS3 = [("A1", "A2"), ("A1", "A3"), ("A2", "A3")]


def _arm(out: Path, out_a2: Path, sacct: list[dict[str, str]] | None = None,
         **exp: Any) -> dict[str, Any]:  # fmt: skip
    kit = m10_04()
    mod = kit.load_run_evidence()
    rows = kit.sacct_rows() if sacct is None else sacct
    doc = mod.arm_evidence(out, rows, a3_expect(out_a2, **exp))
    assert isinstance(doc, dict)
    json.dumps(doc)
    return doc


@pytest.fixture
def a2(tmp_path: Path) -> Path:
    return m10_04().write_arm(tmp_path / "out_a2")


# ------------------------------------------------------------------ criterion 5: example 7


def test_example7_a3_arm_passes(tmp_path: Path, a2: Path) -> None:
    out = write_a3_arm(tmp_path / "out_a3", reference=m10_04().SHA_PEPTIDES)
    doc = _arm(out, a2)
    assert doc["failures"] == [], doc["failures"]
    assert doc["arm"] == "A3"
    assert doc["expected_warnings"] == ["no_spectrum_pairing"]
    assert doc["inputs"][REFERENCE_ROLE] == m10_04().SHA_PEPTIDES


def test_example7_cli_exit_0(tmp_path: Path, a2: Path) -> None:
    kit = m10_04()
    out = write_a3_arm(tmp_path / "out_a3", reference=kit.SHA_PEPTIDES)
    dest = tmp_path / "pass2" / "evidence" / "a3.json"
    res = kit.run_evidence_cli(tmp_path, *a3_cli_args(out, dest, a2))
    assert res.returncode == 0, res.stdout + res.stderr
    doc = json.loads(dest.read_text())
    assert doc["failures"] == [] and doc["arm"] == "A3"
    assert doc["expected_warnings"] == ["no_spectrum_pairing"]


def test_no_spectrum_pairing_absent_is_not_a_failure(tmp_path: Path, a2: Path) -> None:
    out = write_a3_arm(tmp_path / "out_a3", reference=m10_04().SHA_PEPTIDES,
                       spectrum_warning=False)  # fmt: skip
    doc = _arm(out, a2)
    assert doc["failures"] == []
    assert "no_spectrum_pairing" not in doc.get("expected_warnings", [])


# ------------------------------------------------------------------ A10


@pytest.mark.parametrize("case", ["no_file", "placement_only"])
def test_a10_no_denovo_precision(tmp_path: Path, a2: Path, case: str) -> None:
    out = write_a3_arm(tmp_path / "out_a3", reference=m10_04().SHA_PEPTIDES, fdr=False)
    if case == "placement_only":
        write_mapping_fdr(out, n_bins=0)
    assert _arm(out, a2)["failures"] == ["no_denovo_precision"]


@pytest.mark.parametrize("flags", ["c_assumed", "c_assumed;unreported", "unreported;c_assumed"])
def test_a10_c_assumed_on_any_row(tmp_path: Path, a2: Path, flags: str) -> None:
    out = write_a3_arm(tmp_path / "out_a3", reference=m10_04().SHA_PEPTIDES)
    path = out / "mapping" / "mapping_fdr.tsv"
    with path.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    header = list(rows[0])
    dn = [i for i, r in enumerate(rows) if r["section"] == "denovo_precision"]
    rows[dn[4]]["flags"] = flags  # one bin only
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    assert _arm(out, a2)["failures"] == ["c_assumed"]


def test_a10_unreported_alone_is_not_c_assumed(tmp_path: Path, a2: Path) -> None:
    out = write_a3_arm(tmp_path / "out_a3", reference=m10_04().SHA_PEPTIDES,
                       fdr_flags="unreported")  # fmt: skip
    assert _arm(out, a2)["failures"] == []


# ------------------------------------------------------------------ A11 and example 8


def test_a11_reference_is_not_a2s_table(tmp_path: Path, a2: Path) -> None:
    kit = m10_04()
    wrong = "5" * 64
    out = write_a3_arm(tmp_path / "out_a3", reference=wrong)
    inputs = {**kit.EXPECT_INPUTS, "peptides": "4" * 64, REFERENCE_ROLE: wrong}
    doc = _arm(out, a2, inputs=inputs)  # A5 agrees with the manifest; only A11 can fail
    assert doc["failures"] == ["reference_not_a2"]


def test_a11_follows_the_accepted_out_a2(tmp_path: Path) -> None:
    """The reference must equal *this* out_a2's peptides sha, not a constant."""
    kit = m10_04()
    other = "6" * 64
    a2 = kit.write_arm(tmp_path / "out_a2", inputs={"peptides": other})
    out = write_a3_arm(tmp_path / "out_a3", reference=other)
    inputs = {**kit.EXPECT_INPUTS, "peptides": "4" * 64, REFERENCE_ROLE: other}
    assert _arm(out, a2, inputs=inputs)["failures"] == []
    out2 = write_a3_arm(tmp_path / "out_a3b", reference=kit.SHA_PEPTIDES)
    inputs[REFERENCE_ROLE] = kit.SHA_PEPTIDES
    assert _arm(out2, a2, inputs=inputs)["failures"] == ["reference_not_a2"]


def test_example8_submitted_without_reference(tmp_path: Path, a2: Path) -> None:
    kit = m10_04()
    out = write_a3_arm(tmp_path / "out_a3", reference=None, fdr_flags="c_assumed",
                       spectrum_warning=False)  # fmt: skip
    failures = _arm(out, a2)["failures"]
    assert {"c_assumed", "reference_not_a2"} <= set(failures), failures
    assert f"input_missing:{REFERENCE_ROLE}" in failures, "A5 covers the four roles"
    dest = tmp_path / "a3.json"
    res = kit.run_evidence_cli(tmp_path, *a3_cli_args(out, dest, a2))
    assert res.returncode == 1, res.stdout + res.stderr


# ------------------------------------------------------------------ A1-A9 unchanged for --arm A3


@pytest.mark.parametrize(
    ("case", "expected"),
    [("no_manifest", "no_manifest"), ("dirty", "dirty_checkout"), ("cpus", "cpu_cap:16"),
     ("peptides", "input_mismatch:peptides"), ("report", "no_report")],
)  # fmt: skip
def test_rule_set_a_applies_to_a3(tmp_path: Path, a2: Path, case: str, expected: str) -> None:
    kit = m10_04()
    ref = kit.SHA_PEPTIDES
    kw: dict[str, Any] = {}
    sacct = None
    if case == "no_manifest":
        kw["manifest"] = False
    elif case == "dirty":
        kw["dirty"] = True
    elif case == "report":
        kw["report"] = None
    elif case == "cpus":
        sacct = kit.sacct_rows(tasks=dict.fromkeys(kit.TASK_JOBS, 16))
    out = write_a3_arm(tmp_path / "out_a3", reference=ref, **kw)
    exp: dict[str, Any] = {}
    if case == "peptides":
        exp["inputs"] = {**kit.EXPECT_INPUTS, "peptides": "0" * 64, REFERENCE_ROLE: ref}
    failures = _arm(out, a2, sacct=sacct, **exp)["failures"]
    assert expected in failures, failures
    if case == "no_manifest":
        assert failures == ["no_manifest"]


def test_reference_sha_is_checked_by_a5(tmp_path: Path, a2: Path) -> None:
    kit = m10_04()
    out = write_a3_arm(tmp_path / "out_a3", reference=kit.SHA_PEPTIDES)
    inputs = {**kit.EXPECT_INPUTS, "peptides": "4" * 64, REFERENCE_ROLE: "0" * 64}
    assert f"input_mismatch:{REFERENCE_ROLE}" in _arm(out, a2, inputs=inputs)["failures"]


def test_a1_a2_arms_ignore_a10_a11(tmp_path: Path) -> None:
    """A10/A11 are A3 rules: an A2 arm without mapping_fdr.tsv or a reference still passes."""
    kit = m10_04()
    out = kit.write_arm(tmp_path / "out_a2")
    assert kit.arm_doc(out)["failures"] == []


# ------------------------------------------------------------------ criterion 6: rule set C, N = 3


@pytest.fixture(scope="module")
def compare3(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """``CMP``, ``CMP2`` (same seed) and ``cmp_a1_a2`` from the real ``giap compare``."""
    tmp = tmp_path_factory.mktemp("m10_10_cmp")
    a1, a2, a3 = three_arms(tmp)
    for name in ("cmp_a1_a2_a3", "cmp2"):
        res = giap_compare_n([a1, a2, a3], LABELS3, tmp / name, env_root=tmp / "env")
        assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
    res = giap_compare_n([a1, a2], LABELS3[:2], tmp / "cmp_a1_a2", env_root=tmp / "env")
    assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
    return tmp


@pytest.fixture
def cmp3(compare3: Path, tmp_path: Path) -> Path:
    for name in ("cmp_a1_a2_a3", "cmp2", "cmp_a1_a2", "out_a1"):
        shutil.copytree(compare3 / name, tmp_path / name)
    return tmp_path


def _compare(root: Path, labels: tuple[str, ...] = LABELS3, *, cmp: str = "cmp_a1_a2_a3",
             pair_cmp: Path | None = None) -> dict[str, Any]:  # fmt: skip
    mod = m10_04().load_run_evidence()
    kw: dict[str, Any] = {"genes": root / "out_a1" / "tables" / "genes.tsv"}
    if pair_cmp is not None:
        kw["pair_cmp"] = pair_cmp
    doc = mod.compare_evidence(root / cmp, root / "cmp2", labels, **kw)
    assert isinstance(doc, dict)
    json.dumps(doc)
    return doc


def _read(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return list(r.fieldnames or []), list(r)


def _write(path: Path, header: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def _both(root: Path, table: str) -> list[Path]:
    """The table in CMP and CMP2 (edited alike, so C6 stays quiet)."""
    return [root / d / "tables" / f"{table}.tsv" for d in ("cmp_a1_a2_a3", "cmp2")]


def test_three_way_compare_passes(cmp3: Path) -> None:
    doc = _compare(cmp3)
    assert doc["failures"] == [], doc["failures"]
    assert doc["labels"] == list(LABELS3)
    _, ab = _read(cmp3 / "cmp_a1_a2_a3" / "tables" / "compare_abundance.tsv")
    assert [(r["run_a"], r["run_b"]) for r in ab] == PAIRS3
    expected = [w for w in doc["expected_warnings"] if "missing_in:A3" in w and "abundance" in w]
    assert len(expected) == 2, doc["expected_warnings"]  # (A1, A3) and (A2, A3)


def test_c2_abundance_pairs_generalised(cmp3: Path) -> None:
    for path in _both(cmp3, "compare_abundance"):
        header, rows = _read(path)
        _write(path, header, [r for r in rows if (r["run_a"], r["run_b"]) != ("A1", "A3")])
    assert "abundance_rows:2" in _compare(cmp3)["failures"]


def test_c2_abundance_pairs_in_label_order(cmp3: Path) -> None:
    for path in _both(cmp3, "compare_abundance"):
        header, rows = _read(path)
        _write(path, header, [rows[0], rows[2], rows[1]])
    assert "abundance_rows:3" in _compare(cmp3)["failures"]


def test_c3_matrix_n_times_g(cmp3: Path) -> None:
    n_matrix = 0
    for path in _both(cmp3, "compare_detection"):
        header, rows = _read(path)
        matrix = [i for i, r in enumerate(rows) if r["section"] == "matrix"]
        n_matrix = len(matrix)
        drop = next(i for i in matrix if rows[i]["run"] == "A3")
        _write(path, header, [r for i, r in enumerate(rows) if i != drop])
    assert n_matrix == 3 * m10_04().N_GENES, "N x G matrix rows from giap compare"
    assert f"matrix_rows:{n_matrix - 1}" in _compare(cmp3)["failures"]


def test_c4_one_row_per_family_and_pair(cmp3: Path) -> None:
    for path in _both(cmp3, "compare_stability"):
        header, rows = _read(path)
        keep = [
            r
            for r in rows
            if not (r["family"] == "novel_orfs" and (r["run_a"], r["run_b"]) == ("A2", "A3"))
        ]
        _write(path, header, keep)  # fmt: skip
    assert "stability_family:novel_orfs" in _compare(cmp3)["failures"]


def test_c6_cmp2_differs(cmp3: Path) -> None:
    path = cmp3 / "cmp2" / "tables" / "compare_stability.tsv"
    path.write_bytes(path.read_bytes() + b"\n")
    assert "nondeterministic:compare_stability" in _compare(cmp3)["failures"]


def test_c7_equal_pair_rows_give_no_warning(cmp3: Path) -> None:
    """``pass1/cmp_a1_a2`` made of the (A1, A2) rows of the three-way tables."""
    pair = cmp3 / "pass1_cmp"
    shutil.copytree(cmp3 / "cmp_a1_a2", pair)
    for t in ("compare_abundance", "compare_stability"):
        header, rows = _read(cmp3 / "cmp_a1_a2_a3" / "tables" / f"{t}.tsv")
        _write(pair / "tables" / f"{t}.tsv", header,
               [r for r in rows if (r["run_a"], r["run_b"]) == ("A1", "A2")])  # fmt: skip
    doc = _compare(cmp3, pair_cmp=pair)
    assert doc["failures"] == []
    assert not [w for w in doc.get("warnings", []) if w.startswith("a1_a2_rows_differ")]


@pytest.mark.parametrize("table", ["compare_abundance", "compare_stability"])
def test_c7_differing_pair_rows_warn_never_fail(cmp3: Path, table: str) -> None:
    pair = cmp3 / "pass1_cmp"
    shutil.copytree(cmp3 / "cmp_a1_a2", pair)
    for t in ("compare_abundance", "compare_stability"):
        header, rows = _read(cmp3 / "cmp_a1_a2_a3" / "tables" / f"{t}.tsv")
        keep = [dict(r) for r in rows if (r["run_a"], r["run_b"]) == ("A1", "A2")]
        if t == table:  # one cell of one (A1, A2) row differs from the three-way table
            col = next(c for c in header if c not in ("run_a", "run_b", "family", "method")
                       and keep[0][c] not in ("", "NA"))  # fmt: skip
            keep[0][col] = keep[0][col] + "9"
        _write(pair / "tables" / f"{t}.tsv", header, keep)
    doc = _compare(cmp3, pair_cmp=pair)
    assert doc["failures"] == [], "C7 is a warning, never a failure"
    assert f"a1_a2_rows_differ:{table}" in doc["warnings"], doc.get("warnings")


def test_two_labels_keep_m10_04_behaviour(cmp3: Path) -> None:
    root = cmp3
    shutil.copytree(root / "cmp_a1_a2", root / "cmp2", dirs_exist_ok=True)
    doc = _compare(root, ("A1", "A2"), cmp="cmp_a1_a2")
    assert doc["failures"] == [], doc["failures"]
    assert doc["labels"] == ["A1", "A2"]


def test_compare_cli_accepts_three_labels(cmp3: Path) -> None:
    kit = m10_04()
    dest = cmp3 / "evidence" / "compare3.json"
    args = ["compare", str(cmp3 / "cmp_a1_a2_a3"), "--second", str(cmp3 / "cmp2"),
            "--labels", "A1,A2,A3", "--genes", str(cmp3 / "out_a1" / "tables" / "genes.tsv"),
            "--out", str(dest)]  # fmt: skip
    res = kit.run_evidence_cli(cmp3, *args)
    assert res.returncode == 0, res.stdout + res.stderr
    assert json.loads(dest.read_text())["failures"] == []


# ------------------------------------------------------------- criterion 7: pilot_compare.sbatch


@pytest.fixture
def job(tmp_path: Path):  # type: ignore[no-untyped-def]
    kit = m10_04()
    box = kit.compare_job_sandbox(tmp_path)
    pass1, pass2 = box.runs / "pass1", box.runs / "pass2"
    a1, a2 = kit.two_arms(pass1)
    a3 = write_a3_arm(pass2 / "out_a3", reference=kit.SHA_PEPTIDES)
    return box, a1, a2, a3, pass2


def _compare_calls(box) -> list[list[str]]:  # type: ignore[no-untyped-def]
    return [c for c in box.python_calls() if "compare" in c and any("giap.cli" in a for a in c)]


def test_four_arguments_run_the_three_way_compare(job) -> None:  # type: ignore[no-untyped-def]
    box, a1, a2, a3, pass2 = job
    cmp = pass2 / "cmp_a1_a2_a3"
    res = box.run(
        str(a1), str(a2), str(a3), str(cmp), env={"GIAP_PYTHON": str(box.recording_python())}
    )
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    ev = pass2 / "evidence" / "compare3.json"
    assert ev.is_file(), "the evidence goes to $(dirname CMP)/evidence/compare3.json"
    doc = json.loads(ev.read_text())
    assert doc["failures"] == [] and doc["labels"] == list(LABELS3)
    calls = _compare_calls(box)
    assert len(calls) == 2, calls
    for c in calls:
        assert c[c.index("--labels") + 1] == "A1,A2,A3"
        assert [str(a1), str(a2), str(a3)] == [a for a in c if a in (str(a1), str(a2), str(a3))]
    assert box.repo_is_clean()


def test_three_arguments_keep_m10_04_behaviour(job) -> None:  # type: ignore[no-untyped-def]
    box, a1, a2, _, _ = job
    cmp = box.runs / "pass1" / "cmp_a1_a2"
    res = box.run(str(a1), str(a2), str(cmp), env={"GIAP_PYTHON": str(box.recording_python())})
    assert res.returncode == 0, res.stdout[-3000:] + res.stderr[-3000:]
    assert (box.runs / "pass1" / "evidence" / "compare.json").is_file()
    assert _compare_calls(box), "the 3-argument form runs giap compare"
    assert not (box.runs / "pass1" / "evidence" / "compare3.json").exists()
    for c in _compare_calls(box):
        assert c[c.index("--labels") + 1] == "A1,A2"


@pytest.mark.parametrize("n", [2, 5])
def test_two_or_five_arguments_exit_2(job, n: int) -> None:  # type: ignore[no-untyped-def]
    box, a1, a2, a3, pass2 = job
    cmp = pass2 / "cmp_x"
    five = [str(a1), str(a2), str(a3), str(pass2 / "x"), str(cmp)]
    args = [str(a1), str(cmp)] if n == 2 else five
    res = box.run(*args)
    assert res.returncode == 2, res.stdout + res.stderr
    assert not cmp.exists()
    assert _compare_calls(box) == []
