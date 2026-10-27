"""``scripts/dev/rename_gippa.py`` (ADR-0064 item 2; dev tooling) on small temporary trees."""

from __future__ import annotations

import gzip
import importlib.util
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
    "tests/data/recorded/uniop/x/uniop_cds.gff", "tests/data/snapshot/PROVENANCE.md",
    "tests/data/proteins.faa", "tests/data/blob.gz", "tests/data/plot.png",
    ".github/workflows/ci.yml",
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
    assert readme.startswith("GIPPA: see GIPPA_SPEC.md; run `gippa run`")
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
    assert "checkout path (does not move)" in out and "uncategorised: 0" in out
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


def test_a_transition_whose_anchor_is_missing_raises(rg: ModuleType) -> None:
    text = "squeue --me | grep -c ' giap-' || true\n"  # a later main changed the R14 line
    with pytest.raises(rg.TransitionError, match=r"\(giap\|gippa\)-"):
        rg.transform("scripts/session_start.sh", text)
    with pytest.raises(rg.TransitionError):
        rg.transform("scripts/probe_env.py", "giap_jobs = [ln for ln in q if 'giap-' in ln]\n")
    with pytest.raises(rg.TransitionError, match="GIAP_X fallback"):
        rg.transform("scripts/farmshare_env.sh", "export GIAP_X=${GIAP_X:-1}\n")


def test_check_fails_when_a_transition_end_state_is_absent(
    rg: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    repo = _repo(tmp_path)
    rg.main(["--apply", "--repo", str(repo)])
    assert rg.main(["--check", "--repo", str(repo)]) == 0
    # gippa-only R14 counting with a changed anchor: no giap left, no transition fires
    (repo / "scripts/session_start.sh").write_text("squeue --me | grep -c ' gippa-' || true\n")
    capsys.readouterr()
    assert rg.main(["--check", "--repo", str(repo)]) == 1
    assert "LEFTOVER transition scripts/session_start.sh" in capsys.readouterr().out


def test_apply_refuses_a_broken_transition_before_writing(rg: ModuleType, tmp_path: Path) -> None:
    tree = {**TREE, "scripts/session_start.sh": "squeue --me | grep -c ' giap-' || true\n"}
    repo = _repo(tmp_path, tree)
    before = _snapshot(repo)
    with pytest.raises(SystemExit, match=r"session_start\.sh"):
        rg.main(["--apply", "--repo", str(repo)])
    assert _snapshot(repo) == before


def test_indented_export_is_aliased(rg: ModuleType) -> None:
    env = FARMSHARE + 'if [ -d /scratch ]; then\n  export GIAP_TMP="${GIAP_TMP:-/scratch}"\nfi\n'
    out = rg.transform("scripts/farmshare_env.sh", env)
    assert '  export GIPPA_TMP="${GIPPA_TMP:-${GIAP_TMP:-/scratch}}"' in out
    assert 'export GIAP_TMP="$GIPPA_TMP"' in out.split(rg.ALIAS_MARK)[1]
    assert rg.transform("scripts/farmshare_env.sh", out) == out  # idempotent


def test_probe_env_counts_both_prefixes(rg: ModuleType) -> None:
    line = (
        "    giap_jobs = [line for line in queue.splitlines() "
        """if " giap-" in f" {line.split(' ', 1)[-1]}"]\n"""
    )
    out = rg.transform("scripts/probe_env.py", line)
    assert 'for p in (" giap-", " gippa-"))' in out and "gippa_jobs = [" in out
    assert rg.transform("scripts/probe_env.py", out) == out


def test_setup_lock_is_never_renamed(rg: ModuleType) -> None:
    text = 'LOCK = ".giap_setup.lock"\n'
    assert rg.transform("src/gippa/setupdb/run.py", text) == text


def test_a_move_onto_an_existing_file_raises(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path, {"d/x_giap.txt": "a\n", "d/x_gippa.txt": "b\n"})
    with pytest.raises(SystemExit, match=r"d/x_gippa\.txt exists, cannot move d/x_giap\.txt"):
        rg.main(["--apply", "--repo", str(repo)])
    with pytest.raises(SystemExit, match="exists"):
        rg.move_paths(repo)  # the guard inside the loop too (no git mv -k, no endless loop)
    assert _snapshot(repo) == {"d/x_giap.txt": b"a\n", "d/x_gippa.txt": b"b\n"}


# ---------------------------------------------------------------- Q-K32: every document rewritten

ADR64 = "docs/adr/ADR-0064-rename-giap-to-gippa.md"
README_430 = """\
- the first release: PyPI and the hosted docs;
- the final pilot report, and the code rename: until it lands, the package, the command and the
  environment variables are `giap` (`giap run`, `GIAP_DB_DIR`);
- genetic code 4 is planned after 1.0.
"""
DOCS_TREE: dict[str, bytes | str] = {
    **TREE,
    "README.md": "# GIAP\n\n**GIAP** (genome-integrated analysis of proteomics data) takes\n"
    + README_430,
    "GIAP_SPEC.md": "# GIAP spec\n\n**GIAP** is short for *genome-integrated analysis of\n"
    "proteomics data*.\n",
    "mkdocs.yml": "site_name: GIAP\nsite_description: Genome-integrated analysis of proteomics "
    f"data for bacterial genomes\nnav:\n  - adr/{ADR64.split('/')[-1]}\n",
    "docs/adr/ADR-0012-x.md": "GIAP decided; see GIAP_SPEC.md §6\n",
    ADR64: "# ADR-0064: Rename GIAP to GIPPA\n`giap` becomes `gippa`\n",
    "docs/adr/index.md": "- [ADR-0012: GIAP x](ADR-0012-x.md): accepted\n"
    f"- [ADR-0064: Rename GIAP to GIPPA]({ADR64.split('/')[-1]}): accepted-by-khoa\n",
    "docs/specs/M2-08-giap-setup-core.md": "`giap setup` (GIAP_SPEC.md §7)\n",
    "docs/backlog.md": "# GIAP backlog\n",
    "docs/validation/run.md": "giap run on FarmShare\n",
    "docs/marker.md": "touch `.giap-complete`; lock `.giap_setup.lock`\n",
}


def test_documents_are_rewritten_and_the_spec_moves(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path, DOCS_TREE)
    rg.main(["--apply", "--repo", str(repo)])
    snap = _snapshot(repo)
    assert "GIAP_SPEC.md" not in snap
    assert snap["GIPPA_SPEC.md"].decode().startswith("# GIPPA spec\n\n**GIPPA** is short for")
    assert snap["docs/adr/ADR-0012-x.md"] == b"GIPPA decided; see GIPPA_SPEC.md \xc2\xa76\n"
    spec = snap["docs/specs/M2-08-gippa-setup-core.md"]
    assert spec == b"`gippa setup` (GIPPA_SPEC.md \xc2\xa77)\n"
    assert snap["docs/backlog.md"] == b"# GIPPA backlog\n"
    assert snap["PROGRESS.md"] == b"# GIPPA progress\n"
    assert snap["validation/reports/run.md"] == b"the GIPPA run\n"
    assert snap["docs/validation/run.md"] == b"gippa run on FarmShare\n"
    assert snap[ADR64].decode() == DOCS_TREE[ADR64]  # the rename record keeps the old name
    index = snap["docs/adr/index.md"].decode()
    assert index.startswith("- [ADR-0012: GIPPA x](ADR-0012-x.md)")
    assert f"- [ADR-0064: Rename GIAP to GIPPA]({ADR64.split('/')[-1]})" in index
    assert f"  - adr/{ADR64.split('/')[-1]}\n" in snap["mkdocs.yml"].decode()  # link still valid


def test_the_expansion_is_replaced_in_every_variant(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path, DOCS_TREE)
    rg.main(["--apply", "--repo", str(repo)])
    new = "genome integrated prokaryotic proteomics analysis"
    readme = (repo / "README.md").read_text()
    assert f"**GIPPA** ({new}) takes" in readme
    spec = (repo / "GIPPA_SPEC.md").read_text()
    assert f"short for *{new}*." in spec  # wrapped over two lines before
    mk = (repo / "mkdocs.yml").read_text()
    assert f"site_description: {new.capitalize()}\n" in mk  # capitalised; "for bacterial..." gone
    assert rg.transform("x/meta.yml", "Genome-integrated annotation of proteomics (this)") == (
        f"{new.capitalize()} (this)"
    )
    with pytest.raises(rg.TransitionError, match="old expansion"):
        rg.transform("docs/x.md", "genome-integrated analytics of proteomics data\n")


def test_readme_drops_the_rename_clause(rg: ModuleType, tmp_path: Path) -> None:
    repo = _repo(tmp_path, DOCS_TREE)
    rg.main(["--apply", "--repo", str(repo)])
    readme = (repo / "README.md").read_text()
    assert "- the final pilot report;\n- genetic code 4" in readme
    assert "until it lands" not in readme and "code rename" not in readme
    rg.main(["--apply", "--repo", str(repo)])
    assert (repo / "README.md").read_text() == readme
    with pytest.raises(rg.TransitionError, match="README rename clause"):
        rg.transform("README.md", "- the code rename: until it lands, `giap` is the name\n")


def test_check_flags_giap_left_in_docs_but_not_kept_data_tokens(
    rg: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    repo = _repo(tmp_path, DOCS_TREE)
    rg.main(["--apply", "--repo", str(repo)])
    capsys.readouterr()
    assert rg.main(["--check", "--repo", str(repo)]) == 0
    out = capsys.readouterr().out
    assert "docs leftovers: 0" in out and f"kept in docs (rename-record): {ADR64}" in out
    assert (repo / "docs/marker.md").read_text() == DOCS_TREE["docs/marker.md"]
    # a peptide-like span the content map keeps is still a docs leftover
    (repo / "docs/late.md").write_text("ok `.giap-complete`\nthe MKGIAPLL story\n")
    _git(repo, "add", "-A")
    assert rg.main(["--check", "--repo", str(repo)]) == 1
    out = capsys.readouterr().out
    assert "LEFTOVER doc docs/late.md:2: the MKGIAPLL story" in out
    assert "docs/late.md:1" not in out and "docs leftovers: 1" in out


def test_the_old_to_new_arrow_does_not_read_gippa_to_gippa(rg: ModuleType) -> None:
    text = "after the GIAP → GIPPA rename PR; rename GIAP → GIPPA (repo)\n"
    out = rg.transform("docs/release.md", text)
    assert out == "after the GIPPA rename PR; rename to GIPPA (repo)\n"
    with pytest.raises(rg.TransitionError, match="GIPPA → GIPPA"):
        rg.transform("docs/x.md", "GIAP->GIPPA\n")
