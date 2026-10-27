"""``scripts/dev/rename_gippa.py`` and ``scripts/dev/rename_compare.py`` (ADR-0064 items 2 and 8;
dev tooling) on small temporary trees."""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

DEV = Path(__file__).resolve().parents[2] / "scripts/dev"


def _load(name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, DEV / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rg() -> ModuleType:
    return _load("rename_gippa")


@pytest.fixture(scope="module")
def rc() -> ModuleType:
    return _load("rename_compare")


GZ = gzip.compress(b">p giap\nMKGIAPLL\n", mtime=0)
PNG = b"\x89PNG\r\n\x1a\n\0\0giap"
FARMSHARE = """\
# GIAP environment on FarmShare. Source it; never execute it.
export GIAP_REPO="${GIAP_REPO:-/home/classes/bios/270/khoa/giap}"
export GIAP_ROOT="${GIAP_ROOT:-$GIAP_REPO/local}"
export GIAP_DB_DIR="$GIAP_ROOT/db"
export MAMBA_ROOT_PREFIX="$GIAP_ROOT/micromamba"
export PATH="$GIAP_ROOT/bin:$MAMBA_ROOT_PREFIX/envs/giap-wf/bin:$PATH"
"""
TREE: dict[str, bytes | str] = {
    "src/giap/__init__.py": 'from giap.config import GiapConfig\nNAME = "GIAP"\n',
    "src/giap/io/plot_files.py": '#: fixed SVG id salt\nHASHSALT = "giap"\n',
    "src/giap/novel/ids.py": (
        'PREFIX = "GIAPORF_"\nSEQ = "MKGIAPLL"\nTAB = "\\tgiap"\nB = b"giap"\n'
    ),
    "tests/sim/giapsim/__init__.py": 'SOURCE = "giapsim"\n',
    "tests/unit/test_giapsim_x.py": "import giapsim\n",
    "workflow/tests/bin_giap/giap": '#!/bin/sh\nexec giap "$@"\n',
    "README.md": "GIAP: see GIAP_SPEC.md; run `giap run` from /home/classes/bios/270/khoa/giap.\n",
    "GIAP_SPEC.md": "# GIAP spec\n",
    "docs/adr/ADR-0005-x.md": "GIAP decided\n",
    "docs/adr/ADR-0065-later.md": "GIAP later\n",
    "docs/specs/M1-01-giap-x.md": "giap spec\n",
    "PROGRESS.md": "# GIAP progress\n",
    "tests/data/recorded/uniop/x/uniop_cds.gff": "chr1\tGIAP\tCDS\t1\t9\t.\t+\t0\tID=a\n",
    "tests/data/snapshot/PROVENANCE.md": "made with giap run\n",
    "tests/data/snapshot/OUT/run_manifest.json": '{"giap_version": "0.0.1"}\n',
    "tests/data/proteins.faa": ">p1 GIAP fixture\nMGIAPK\n",
    "tests/data/hits.tsv": "protein_id\tpeptide\nGIAPORF_1\tGIAP\ngiap_2\tMKGIAPL\n",
    "tests/data/blob.gz": GZ,
    "tests/data/plot.png": PNG,
    "validation/reports/run.md": "the GIAP run\n",
    ".github/workflows/ci.yml": "run: mypy src/giap\n",
    ".gitignore": ".giap/\n",
    "scripts/farmshare_env.sh": FARMSHARE,
    "scripts/session_start.sh": "squeue --me 2>/dev/null | grep ' giap-' || true\n",
    "scripts/slurm/calibration_sharded.sh": 'case "$name" in giap-*) ;; *) continue ;; esac\n',
    "scripts/dev/record.sh": 'touch "$d/.giap-complete"\n',
}
KEPT = [
    "GIAP_SPEC.md", "docs/adr/ADR-0005-x.md", "docs/specs/M1-01-giap-x.md", "PROGRESS.md",
    "tests/data/recorded/uniop/x/uniop_cds.gff", "tests/data/snapshot/PROVENANCE.md",
    "tests/data/proteins.faa", "tests/data/blob.gz", "tests/data/plot.png",
    "validation/reports/run.md", ".github/workflows/ci.yml",
]  # fmt: skip


def _git(repo: Path, *args: str) -> str:
    cmd = ["git", "-C", str(repo), "-c", "user.email=t@t", "-c", "user.name=t", *args]
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def _repo(tmp_path: Path, tree: dict[str, bytes | str] = TREE) -> Path:
    repo = tmp_path / "repo"
    for rel, body in tree.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(body if isinstance(body, bytes) else body.encode("utf-8"))
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    return repo


def _snapshot(repo: Path) -> dict[str, bytes]:
    files = _git(repo, "ls-files", "-z").split("\0")
    return {f: (repo / f).read_bytes() for f in files if f}


def _body(rel: str) -> bytes:
    body = TREE[rel]
    return body if isinstance(body, bytes) else body.encode("utf-8")


def test_apply_respects_the_keep_list_and_moves_paths(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path)
    rg.main(["--apply", "--repo", str(repo)])
    snap = _snapshot(repo)
    for rel in KEPT:
        assert snap[rel] == _body(rel), rel
    assert "src/giap/__init__.py" not in snap
    assert (
        snap["src/gippa/__init__.py"] == b'from gippa.config import GippaConfig\nNAME = "GIPPA"\n'
    )
    assert snap["tests/sim/gippasim/__init__.py"] == b'SOURCE = "gippasim"\n'
    assert "tests/unit/test_gippasim_x.py" in snap
    assert snap["workflow/tests/bin_gippa/gippa"] == b'#!/bin/sh\nexec gippa "$@"\n'
    assert snap["docs/adr/ADR-0065-later.md"] == b"GIPPA later\n"  # not a historical ADR
    ids = snap["src/gippa/novel/ids.py"].decode()
    assert 'PREFIX = "GIPPAORF_"' in ids and 'SEQ = "MKGIAPLL"' in ids  # identifier vs peptide
    assert 'TAB = "\\tgippa"' in ids and 'B = b"gippa"' in ids  # no word boundaries
    readme = snap["README.md"].decode()
    assert readme.startswith("GIPPA: see GIAP_SPEC.md; run `gippa run`")
    assert "/home/classes/bios/270/khoa/giap." in readme  # the checkout does not move
    assert (
        snap["tests/data/hits.tsv"] == b"protein_id\tpeptide\nGIPPAORF_1\tGIAP\ngippa_2\tMKGIAPL\n"
    )
    assert snap["tests/data/snapshot/OUT/run_manifest.json"] == b'{"gippa_version": "0.0.1"}\n'
    assert b'HASHSALT = "giap"' in snap["src/gippa/io/plot_files.py"]
    assert snap["scripts/dev/record.sh"] == _body("scripts/dev/record.sh")


def test_transitions(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path)
    rg.main(["--apply", "--repo", str(repo)])
    env = (repo / "scripts/farmshare_env.sh").read_text()
    assert (
        'export GIPPA_REPO="${GIPPA_REPO:-${GIAP_REPO:-/home/classes/bios/270/khoa/giap}}"' in env
    )
    assert "envs/gippa-wf/bin:$MAMBA_ROOT_PREFIX/envs/giap-wf/bin:" in env
    out = subprocess.run(
        ["bash", "-c", f'. "{repo}/scripts/farmshare_env.sh"; echo "$GIPPA_ROOT $GIAP_DB_DIR"'],
        env={"GIAP_ROOT": "/old/root", "PATH": "/usr/bin:/bin"}, capture_output=True, text=True,
        check=True,
    ).stdout.split()  # fmt: skip
    assert out == ["/old/root", "/old/root/db"]  # an old override is honoured; old names exported
    assert "grep -E ' (giap|gippa)-'" in (repo / "scripts/session_start.sh").read_text()
    cal = (repo / "scripts/slurm/calibration_sharded.sh").read_text()
    assert 'case "$name" in giap-*|gippa-*)' in cal
    assert (repo / ".gitignore").read_text() == ".gippa/\n.giap/\n"
    assert "ADR-0064 item 7" in (repo / "src/gippa/io/plot_files.py").read_text()


def test_apply_twice_equals_once(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path)
    rg.main(["--apply", "--repo", str(repo)])
    once = _snapshot(repo)
    rg.main(["--apply", "--repo", str(repo)])
    assert _snapshot(repo) == once


def test_check_catches_leftovers(
    rg: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    repo = _repo(tmp_path)
    assert rg.main(["--check", "--repo", str(repo)]) == 1  # nothing renamed yet
    rg.main(["--apply", "--repo", str(repo)])
    capsys.readouterr()
    assert rg.main(["--check", "--repo", str(repo)]) == 0
    out = capsys.readouterr().out
    assert "GIAP_SPEC token" in out and "uncategorised: 0" in out
    assert "REVIEW possible peptide (review): src/gippa/novel/ids.py:2: MKGIAPLL" in out
    (repo / "src/gippa/late.py").write_text("import giap\n")
    (repo / "scripts/giap_new.sh").write_text("echo\n")
    _git(repo, "add", "-A")
    assert rg.main(["--check", "--repo", str(repo)]) == 1
    out = capsys.readouterr().out
    assert "LEFTOVER src/gippa/late.py:1: import giap" in out
    assert "LEFTOVER path scripts/giap_new.sh (should be scripts/gippa_new.sh)" in out


def test_plan_lists_github_lines_without_writing(
    rg: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    repo = _repo(tmp_path)
    before = _snapshot(repo)
    assert rg.main(["--plan", "--repo", str(repo)]) == 0
    out = capsys.readouterr().out
    assert "+ run: mypy src/gippa" in out
    assert "git mv src/giap src/gippa" in out
    assert _snapshot(repo) == before


def test_the_tooling_never_renames_itself(rg: ModuleType) -> None:
    for rel in ("scripts/dev/rename_gippa.py", "scripts/dev/rename_compare.py",
                "tests/unit/test_rename_gippa_tooling.py"):  # fmt: skip
        assert rg.keep_category(rel) == "rename-tooling" and rg.path_stays(rel)


def test_binary_and_gz_untouched(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path)
    rg.main(["--apply", "--repo", str(repo)])
    assert (repo / "tests/data/blob.gz").read_bytes() == GZ
    assert (repo / "tests/data/plot.png").read_bytes() == PNG


# ---------------------------------------------------------------- the comparator


def _sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def _run(root: Path, name: str, commit: str, table: str, fasta: str, png: bytes = PNG) -> Path:
    files = {
        "tables/genes.tsv": table.encode(),
        "annotation/annotation.gff3": f"c1\t{name.upper()}\tCDS\t1\t9\t.\t+\t0\t.\n".encode(),
        "novel/proteins.faa": fasta.encode(),
        "plots/V01.png": png,
        "mapping/hits.tsv.gz": gzip.compress(f"id\tpeptide\n{name.upper()}ORF_1\tPEPK\n".encode()),
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
    (root / f"logs/{name}-x.log").parent.mkdir(parents=True, exist_ok=True)
    (root / f"logs/{name}-x.log").write_text(commit)
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
    assert rc.compare(b, a, ["extra.txt"]) != []  # logs/ is compared once not excluded
