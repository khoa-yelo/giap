"""M11-12 acceptance 2 (the CI step in ``guards``) and 3 plus definition item 2 (the 1.0 README)."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import yaml

REPO = Path(__file__).resolve().parents[2]
CI: dict[str, Any] = yaml.safe_load((REPO / ".github/workflows/ci.yml").read_text())
README = (REPO / "README.md").read_text()
LICENCE_BODY = (
    "MIT. Third-party tools keep their own licences and run only as separate programs; see\n"
    "[`NOTICE.md`](NOTICE.md)."
)

# ---- CI step (acceptance 2) ---------------------------------------------------------------------


def _guard_steps() -> list[dict[str, Any]]:
    steps: list[dict[str, Any]] = CI["jobs"]["guards"]["steps"]
    return [s for s in steps if "check_changelog.py" in str(s.get("run", ""))]


def test_guards_runs_check_changelog_against_the_base_ref() -> None:
    steps = _guard_steps()
    assert len(steps) == 1, "guards has no (or more than one) check_changelog.py step"
    run = str(steps[0]["run"])
    assert "scripts/ci/check_changelog.py origin/${{ github.base_ref }}" in run


def test_check_changelog_runs_only_on_pull_request() -> None:
    assert CI["jobs"]["guards"]["if"].strip() == "github.event_name == 'pull_request'"
    step = _guard_steps()[0]
    assert "if" not in step or "pull_request" in str(step["if"])


def test_check_changelog_step_is_not_soft() -> None:
    """Severity fail (Q-L9): the step must not swallow the exit code."""
    step = _guard_steps()[0]
    assert not step.get("continue-on-error", False)
    assert "|| true" not in str(step["run"])


def test_check_changelog_runs_only_in_guards() -> None:
    for name, job in CI["jobs"].items():
        if name == "guards":
            continue
        assert "check_changelog.py" not in yaml.safe_dump(job), (
            f"job {name} also runs check_changelog.py"
        )


# ---- README (definition item 2, acceptance 3) ---------------------------------------------------

_LINK = re.compile(r"!?\[[^\]]*\]\(\s*<?([^)\s>]+)>?(?:\s+\"[^\"]*\")?\s*\)")
_REF = re.compile(r"^\s*\[[^\]]+\]:\s*(\S+)", re.M)
_HTML = re.compile(r"""(?:href|src)=["']([^"']+)["']""")


def _targets() -> list[str]:
    body = re.sub(r"```.*?```", "", README, flags=re.S)
    return _LINK.findall(body) + _REF.findall(body) + _HTML.findall(body)


def _relative(t: str) -> bool:
    return not re.match(r"^[a-z][a-z0-9+.-]*:", t, re.I) and not t.startswith(("#", "//"))


def test_every_relative_readme_link_resolves() -> None:
    rel = [t for t in _targets() if _relative(t)]
    assert rel, "README has no relative links at all"
    broken = []
    for t in rel:
        path = t.split("#", 1)[0].split("?", 1)[0]
        if path and not (REPO / path).exists():
            broken.append(t)
    assert not broken, f"broken relative links in README.md: {broken}"


def _headings() -> list[tuple[int, str, int]]:
    """(level, title, offset) for every ATX heading outside fenced code."""
    out = []
    in_fence = False
    offset = 0
    for line in README.splitlines(keepends=True):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
        elif not in_fence:
            m = re.match(r"^(#{1,6})\s+(.*?)\s*#*\s*$", line)
            if m:
                out.append((len(m.group(1)), m.group(2), offset))
        offset += len(line)
    return out


def _section(pattern: str) -> tuple[int, str]:
    hs = _headings()
    for i, (lvl, title, off) in enumerate(hs):
        if re.search(pattern, title, re.I):
            end = len(README)
            for lvl2, _, off2 in hs[i + 1 :]:
                if lvl2 <= lvl:
                    end = off2
                    break
            return lvl, README[off:end]
    raise AssertionError(
        f"README has no heading matching {pattern!r}; headings: {[h[1] for h in hs]}"
    )


def test_status_line_is_no_longer_m0() -> None:
    assert "under development" not in README
    assert "milestone M0" not in README


def test_pypi_version_badge() -> None:
    images = re.findall(r"!\[[^\]]*\]\(([^)\s]+)", README)
    assert any("pypi" in u.lower() and "gippa" in u for u in images), (
        f"no PyPI version badge: {images}"
    )


def test_docs_site_link() -> None:
    urls = [t for t in _targets() if t.startswith("https://")]
    assert any("github.io" in u or "readthedocs" in u for u in urls), f"no docs-site link: {urls}"


def test_install_section_has_pip_and_bioconda() -> None:
    _, sec = _section(r"^install")
    assert "pip install gippa" in sec
    assert re.search(r"(conda|mamba|micromamba) install .*bioconda.*gippa|bioconda::gippa", sec), sec


def test_quickstart_is_four_lines_and_points_at_the_site() -> None:
    _, sec = _section(r"quick ?start")
    blocks = re.findall(r"```[^\n]*\n(.*?)```", sec, flags=re.S)
    assert blocks, "quickstart has no code block"
    lines = [ln for ln in blocks[0].splitlines() if ln.strip()]
    assert 1 <= len(lines) <= 4, lines
    assert re.search(r"https://\S+", sec), "quickstart does not point at the docs site"


def test_citing_section_points_at_citation_cff() -> None:
    _, sec = _section(r"^citing gippa")
    assert "CITATION.cff" in sec
    assert "(CITATION.cff)" in sec, "CITATION.cff should be a relative link"


def test_licence_section_unchanged() -> None:
    _, sec = _section(r"^licen[cs]e$")
    assert LICENCE_BODY in sec


def test_farmshare_development_moved_under_developers() -> None:
    lvl_dev, dev = _section(r"^developers$")
    lvl_fs, fs = _section(r"development on farmshare")
    assert lvl_fs > lvl_dev, "`Development on FarmShare` must be a sub-heading of `Developers`"
    assert fs in dev
    assert "scripts/farmshare_env.sh" in fs


def test_readme_keeps_its_first_paragraph() -> None:
    """M11-02: the CITATION.cff abstract is the README's first paragraph."""
    assert "**GIPPA** (genome integrated prokaryotic proteomics analysis)" in README
