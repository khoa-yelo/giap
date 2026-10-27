"""``scripts/dev/rename_compare.py`` (ADR-0064 item 8; dev tooling) on small temporary trees."""

from __future__ import annotations

import gzip
import hashlib
import importlib
import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import pytest

DEV = Path(__file__).resolve().parents[2] / "scripts/dev"
PNG = b"\x89PNG\r\n\x1a\n\0\0giap"


@pytest.fixture(scope="module")
def rc() -> ModuleType:
    spec = importlib.util.spec_from_file_location("rename_compare", DEV / "rename_compare.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["rename_compare"] = mod
    spec.loader.exec_module(mod)
    return mod


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _run(root: Path, name: str, commit: str, table: str, fasta: str, png: bytes = PNG) -> Path:
    files = {
        "tables/genes.tsv": table.encode(),
        "annotation/annotation.gff3": f"c1\t{name.upper()}\tCDS\t1\t9\t.\t+\t0\t.\n".encode(),
        "novel/proteins.faa": fasta.encode(),
        "plots/V01.png": png,
        "mapping/hits.tsv.gz": gzip.compress(f"id\tpeptide\n{name.upper()}ORF_1\tPEPK\n".encode()),
        "logs/events.jsonl": f'{{"ts":"{commit}","step":"{name}","level":"INFO"}}\n'.encode(),
        "logs/trace.tsv": f"task_id\thash\n1\t{commit}\n".encode(),
    }
    for rel, body in files.items():
        (root / rel).parent.mkdir(parents=True, exist_ok=True)
        (root / rel).write_bytes(body)
    manifest = {
        f"{name}_version": "0.0.1",
        "git_commit": commit,
        "git_dirty": False,
        "provenance_digest": commit * 2,
        "started_at": commit,
        "command_line": [name, "run"],
        "outputs": [{"path": "tables/genes.tsv", "sha256": _sha(files["tables/genes.tsv"])}],
        "inputs": [{"path": "/elsewhere/genome.fna", "sha256": "ab" * 32}],
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest))
    return root


def _pair(tmp_path: Path, **a_over: str) -> tuple[Path, Path]:
    b = _run(tmp_path / "B", "giap", "a" * 40, "gene\tsource\tn\ng1\tgiapsim\t3\n",
             ">p1 giap\nMKGIAPLL\n")  # fmt: skip
    args = {"table": "gene\tsource\tn\ng1\tgippasim\t3\n", "fasta": ">p1 gippa\nMKGIAPLL\n"}
    args |= a_over
    a = _run(tmp_path / "A", "gippa", "b" * 40, args["table"], args["fasta"])
    return b, a


def test_compare_passes_on_a_name_only_difference(rc: ModuleType, tmp_path: Path) -> None:
    b, a = _pair(tmp_path)
    assert rc.compare(b, a, rc.DEFAULT_EXCLUDES) == []
    assert rc.main([str(b), str(a)]) == 0


def test_compare_fails_on_a_one_byte_data_difference(rc: ModuleType, tmp_path: Path) -> None:
    b, a = _pair(tmp_path, table="gene\tsource\tn\ng1\tgippasim\t4\n")
    diffs = rc.compare(b, a, rc.DEFAULT_EXCLUDES)
    assert any(d.startswith("content differs: tables/genes.tsv") for d in diffs), diffs
    assert rc.main([str(b), str(a)]) == 1


def test_compare_fails_on_a_renamed_sequence(rc: ModuleType, tmp_path: Path) -> None:
    b, a = _pair(tmp_path, fasta=">p1 gippa\nMKGIPPALL\n")
    diffs = rc.compare(b, a, rc.DEFAULT_EXCLUDES)
    assert [d.split("\n")[0] for d in diffs] == [
        "content differs: novel/proteins.faa | novel/proteins.faa"
    ]


def test_compare_peptide_column_is_exact(rc: ModuleType) -> None:
    b = rc.tsv_lines("id\tpeptide\nGIAPORF_1\tAGIAPK\n")
    assert b == rc.tsv_lines("id\tpeptide\nGIPPAORF_1\tAGIAPK\n")
    assert b != rc.tsv_lines("id\tpeptide\nGIPPAORF_1\tAGIPPAK\n")


def test_compare_png_must_be_byte_equal(rc: ModuleType, tmp_path: Path) -> None:
    b, a = _pair(tmp_path)
    (a / "plots/V01.png").write_bytes(PNG.replace(b"giap", b"gippa"))
    assert rc.compare(b, a, rc.DEFAULT_EXCLUDES) == ["bytes differ: plots/V01.png | plots/V01.png"]


def test_compare_recomputes_recorded_digests(rc: ModuleType, tmp_path: Path) -> None:
    b, a = _pair(tmp_path)
    man = json.loads((a / "run_manifest.json").read_text())
    man["outputs"][0]["sha256"] = "0" * 64
    (a / "run_manifest.json").write_text(json.dumps(man))
    diffs = rc.compare(b, a, rc.DEFAULT_EXCLUDES)
    assert any(d.startswith("A: run_manifest.json/outputs/0/sha256: recorded") for d in diffs)


def test_compare_pairs_files_by_normalised_path(rc: ModuleType, tmp_path: Path) -> None:
    b, a = _pair(tmp_path)
    (a / "extra.txt").write_text("x")
    assert rc.compare(b, a, rc.DEFAULT_EXCLUDES) == ["only in A: extra.txt"]
    assert rc.compare(b, a, ["extra.txt"]) != []  # logs/trace.tsv is compared once not excluded


# ---------------------------------------------------------------- the #429 review probe cases


def _trees(tmp_path: Path, b: dict[str, bytes], a: dict[str, bytes]) -> tuple[Path, Path]:
    for side, fs in (("B", b), ("A", a)):
        for rel, body in fs.items():
            p = tmp_path / side / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_bytes(body)
    return tmp_path / "B", tmp_path / "A"


def _differs(rc: ModuleType, tmp_path: Path, rel: str, b: bytes, a: bytes) -> list[str]:
    diffs: list[str] = rc.compare(*_trees(tmp_path, {rel: b}, {rel: a}), rc.DEFAULT_EXCLUDES)
    return diffs


def test_probe_csv_quoting_and_column_count(rc: ModuleType, tmp_path: Path) -> None:
    assert _differs(rc, tmp_path, "t.csv", b'id,product\nx,"a,b"\n', b"id,product\nx,a,b\n")


def test_probe_csv_name_only_difference_passes(rc: ModuleType, tmp_path: Path) -> None:
    b, a = b'id,source\n"x,1",giap\n', b'id,source\n"x,1",gippa\n'
    assert _differs(rc, tmp_path, "t.csv", b, a) == []


def test_probe_diann_sequence_columns(rc: ModuleType, tmp_path: Path) -> None:
    b = b"Precursor.Id\tStripped.Sequence\nGIAP2\tGIAP\n"
    a = b"Precursor.Id\tStripped.Sequence\nGIPPA2\tGIPPA\n"
    assert _differs(rc, tmp_path, "r.tsv", b, a)
    only_id = b"Precursor.Id\tStripped.Sequence\nGIPPA2\tGIAP\n"
    assert rc.tsv_lines(b.decode()) != rc.tsv_lines(only_id.decode())


def test_probe_json_peptides_list(rc: ModuleType, tmp_path: Path) -> None:
    b, a = json.dumps({"peptides": ["GIAP"]}), json.dumps({"peptides": ["GIPPA"]})
    assert _differs(rc, tmp_path, "x.json", b.encode(), a.encode())


def test_probe_giapsim_peptide_in_a_text_column(rc: ModuleType, tmp_path: Path) -> None:
    b, a = b"id\tnote\n1\tGIAPSIMK\n", b"id\tnote\n1\tGIPPASIMK\n"
    assert _differs(rc, tmp_path, "p.tsv", b, a)


def test_probe_line_endings_and_final_newline(rc: ModuleType, tmp_path: Path) -> None:
    assert _differs(rc, tmp_path, "a.txt", b"x\r\ny\n", b"x\ny\n")
    assert _differs(rc, tmp_path / "2", "a.txt", b"x\ny\n", b"x\ny")


def test_probe_json_key_order_and_whitespace(rc: ModuleType, tmp_path: Path) -> None:
    assert _differs(rc, tmp_path, "a.json", b'{"a":1,"b":2}', b'{\n "b": 2, "a": 1}\n')
    assert _differs(rc, tmp_path / "2", "a.json", b'{"a":1,"b":2}', b'{"a": 1, "b": 2}')


def test_probe_run_manifest_key_order_counts(rc: ModuleType, tmp_path: Path) -> None:
    b = b'{"git_commit":"x","a":1,"b":2}'
    assert _differs(rc, tmp_path, "run_manifest.json", b, b'{"git_commit":"y","a":1,"b":2}') == []
    assert _differs(rc, tmp_path / "2", "run_manifest.json", b, b'{"git_commit":"x","b":2,"a":1}')


def test_probe_nested_git_commit_is_compared(rc: ModuleType, tmp_path: Path) -> None:
    b, a = b'{"tool":{"git_commit":"x","v":1}}', b'{"tool":{"git_commit":"y","v":1}}'
    assert _differs(rc, tmp_path, "a.json", b, a)
    assert _differs(rc, tmp_path / "2", "run_manifest.json", b, a)


def test_probe_peptide_column_exact(rc: ModuleType, tmp_path: Path) -> None:
    b, a = b"id\tpeptide\n1\tGIAP\n", b"id\tpeptide\n1\tGIPPA\n"
    assert _differs(rc, tmp_path, "p.tsv", b, a)


def test_probe_logs_are_compared_but_ts_and_trace_are_not(rc: ModuleType, tmp_path: Path) -> None:
    ev = '{{"ts":"{}","step":"s","msg":"{}"}}\n'
    b = {
        "logs/events.jsonl": ev.format("t1", "giap ok").encode(),
        "logs/steps/s.events.jsonl": ev.format("t1", "x").encode(),
        "logs/internal_error.log": b"Traceback giap\n",
        "logs/trace.tsv": b"a\n",
    }
    a = {
        "logs/events.jsonl": ev.format("t2", "gippa ok").encode(),
        "logs/steps/s.events.jsonl": ev.format("t2", "x").encode(),
        "logs/internal_error.log": b"Traceback gippa\n",
        "logs/trace.tsv": b"b\n",
    }
    roots = _trees(tmp_path, b, a)
    assert rc.compare(*roots, rc.DEFAULT_EXCLUDES) == []
    assert rc.excluded_files(roots[0], rc.DEFAULT_EXCLUDES) == ["logs/trace.tsv"]
    for rel, body in (
        ("logs/events.jsonl", ev.format("t2", "gippa changed")),
        ("logs/steps/s.events.jsonl", ev.format("t2", "y")),
        ("logs/internal_error.log", "Traceback gippa 2\n"),
    ):
        (roots[1] / rel).write_text(body)
        assert any(d.startswith(f"content differs: {rel}") for d in rc.compare(
            *roots, rc.DEFAULT_EXCLUDES))  # fmt: skip
    assert set(rc.DEFAULT_EXCLUDES) == set(rc.EXCLUDE_REASONS)
    assert not any(g in ("logs/*", "logs/") for g in rc.DEFAULT_EXCLUDES)


def test_main_prints_excluded_files(
    rc: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    b, a = _pair(tmp_path)
    assert rc.main([str(b), str(a)]) == 0
    out = capsys.readouterr().out
    assert "excluded file B: logs/trace.tsv" in out
    assert "excluded logs/trace.tsv: Nextflow trace" in out


def test_volatile_matches_the_manifest_module(rc: ModuleType) -> None:
    try:
        manifest = importlib.import_module("gippa.manifest")
    except ImportError:
        manifest = importlib.import_module("giap.manifest")
    assert (*rc.VOLATILE, rc.NAME_ONLY) == tuple(manifest.VOLATILE)
    assert rc.NAME_ONLY == manifest._NAME_ONLY
