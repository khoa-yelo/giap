"""``scripts/dev/rename_gippa.py`` documentation layer (Q-K32 review fixes of #434).

Self-arrows in any case or code span, blockquote pages, every tracked ``*.md`` checked, the loose
expansion guard and protected spans that never pass the expansion rewrite.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest

DEV = Path(__file__).resolve().parents[2] / "scripts/dev"


@pytest.fixture(scope="module")
def rg() -> ModuleType:
    spec = importlib.util.spec_from_file_location("rename_gippa", DEV / "rename_gippa.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["rename_gippa"] = mod
    spec.loader.exec_module(mod)
    return mod


def _git(repo: Path, *args: str) -> str:
    cmd = ["git", "-C", str(repo), "-c", "user.email=t@t", "-c", "user.name=t", *args]
    return subprocess.run(cmd, check=True, capture_output=True, text=True).stdout


def _repo(tmp_path: Path, tree: dict[str, str]) -> Path:
    repo = tmp_path / "repo"
    for rel, body in tree.items():
        p = repo / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(body)
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    return repo


def _add(repo: Path, rel: str, body: str) -> None:
    (repo / rel).parent.mkdir(parents=True, exist_ok=True)
    (repo / rel).write_text(body)
    _git(repo, "add", "-A")


M11_14 = "docs/specs/M11-14-container-install.md"


def test_a_code_span_self_arrow_is_rewritten(rg: ModuleType) -> None:
    line = "  published. The rename script also renames the image prefix (`giap-` → `gippa-`) in\n"
    assert rg.transform(M11_14, line) == (
        "  published. The rename script also renames the image prefix to `gippa-` in\n"
    )
    assert rg.transform("docs/x.md", "the env prefix (giap_ -> gippa_).\n") == (
        "the env prefix to gippa_.\n"
    )


@pytest.mark.parametrize(
    "text",
    [
        "`giap-wf` → `gippa-wf`\n",
        "Giap → gippa\n",
        "GiapConfig -> GippaConfig\n",
        "rename giap to gippa\n",
        "giapsim → `gippa run` → `giap-x`, `giap-x` → `gippa-x`\n",
    ],
)
def test_any_self_arrow_left_raises(rg: ModuleType, text: str) -> None:
    with pytest.raises(rg.TransitionError, match="GIPPA → GIPPA"):
        rg.transform("docs/x.md", text)


def test_different_tokens_around_an_arrow_are_fine(rg: ModuleType) -> None:
    text = '"""Example 4 (TE7): giapsim -> ``giap run`` -> M4, twice."""\n'
    assert rg.transform("tests/integration/t.py", text) == (
        '"""Example 4 (TE7): gippasim -> ``gippa run`` -> M4, twice."""\n'
    )


def test_a_blockquote_page_is_text_not_sequence(rg: ModuleType) -> None:
    assert rg.file_kind("docs/x.md", b"> GIAP note\n") == "text"
    assert rg.file_kind("notes.txt", b"> GIAP note\n") == "text"
    assert rg.file_kind("tests/data/x.dat", b">p1\nMGIAPK\n") == "sequence"
    assert rg.file_kind("docs/tutorials/data/x.faa", b">p1\nMGIAPK\n") == "sequence"


def test_a_blockquote_page_is_renamed_and_a_leftover_is_flagged(
    rg: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    repo = _repo(tmp_path, {"docs/x.md": "> GIAP note\n"})
    rg.main(["--apply", "--repo", str(repo)])
    assert (repo / "docs/x.md").read_text() == "> GIPPA note\n"
    capsys.readouterr()
    assert rg.main(["--check", "--repo", str(repo)]) == 0
    _add(repo, "docs/x.md", "> the MKGIAPLL note\n")  # a protected span the map keeps
    assert rg.main(["--check", "--repo", str(repo)]) == 1
    out = capsys.readouterr().out
    assert "LEFTOVER doc docs/x.md:1: > the MKGIAPLL note" in out


DOCS = {
    "NOTICE.md": "GIAP notice\n",
    "CHANGELOG.md": "# GIAP changelog\n",
    "changelog.d/+x.feature.md": "giap run\n",
    ".claude/agents/lead.md": "GIAP lead\n",
    "validation/reports/r.md": "GIAP report\n",
    "validation/pilot/README.md": "GIAP pilot\n",
    "validation/pilot/reproduce.md": "giap run\n",
    "validation/pilot/REPORT.md": "GIAP run (pinned)\n",
    "tests/data/snapshot/PROVENANCE.md": "made with giap\n",
    "src/x/notes.rst": "GIAP notes\n",
}


def test_every_markdown_file_is_in_the_docs_layer(
    rg: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    for path in DOCS:
        assert rg.is_doc(path), path
    assert not rg.is_doc("src/x/y.py")
    repo = _repo(tmp_path, DOCS)
    rg.main(["--apply", "--repo", str(repo)])
    assert (repo / "NOTICE.md").read_text() == "GIPPA notice\n"
    assert (repo / "validation/pilot/REPORT.md").read_text() == "GIAP run (pinned)\n"
    capsys.readouterr()
    assert rg.main(["--check", "--repo", str(repo)]) == 0
    out = capsys.readouterr().out
    assert "kept in docs (pilot-artefact): validation/pilot/REPORT.md" in out
    assert "kept in docs (provenance): tests/data/snapshot/PROVENANCE.md" in out
    leftovers = [p for p in DOCS if not rg.keep_category(p)]
    for path in leftovers:
        _add(repo, path, "the MKGIAPLL story\n")
    assert rg.main(["--check", "--repo", str(repo)]) == 1
    out = capsys.readouterr().out
    for path in leftovers:
        assert f"LEFTOVER doc {path}:1: the MKGIAPLL story" in out, path
    assert f"docs leftovers: {len(leftovers)}" in out


def test_proteomics_prose_is_not_an_old_expansion(rg: ModuleType) -> None:
    text = "a statistical analysis of proteomics data\n"
    assert rg.transform("docs/x.md", text) == text
    with pytest.raises(rg.TransitionError, match="old expansion"):
        rg.transform("docs/x.md", "genome-integrated proteomics\n")


def test_protected_spans_never_pass_the_expansion(rg: ModuleType) -> None:
    line = "- [ADR-0064: genome-integrated analysis of proteomics data](a.md): accepted\n"
    text = line + "GIAP: genome-integrated analysis of proteomics data\n"
    out = rg.transform("docs/adr/index.md", text)
    assert out == line + "GIPPA: genome integrated prokaryotic proteomics analysis\n"


def test_this_file_is_rename_tooling(rg: ModuleType) -> None:
    assert rg.keep_category("tests/unit/test_rename_gippa_docs.py") == "rename-tooling"
    assert rg.path_stays("tests/unit/test_rename_gippa_docs.py")
