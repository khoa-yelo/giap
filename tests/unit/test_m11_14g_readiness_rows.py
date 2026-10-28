"""Readiness rows under 14g: calibration evidence reused by ``reusable(C, R)``, and H4 image keys
(M11-14 item 4 "Evidence reuse", examples 12, 13, 13b, 13d, edge cases "An image with no process"
and "A tool fails to build", AC 13; ADR-0060 §1, §4).

Offline. The calibration rows run ``readiness.py compute --offline`` on the M11-03 fixture
repository (``tests/helpers/release_readiness.py``), whose ``nightly.tsv`` rows are logged for the
fixture commit R; each test adds commits after R, so R becomes the evidence commit C.

Interface fixed by these tests:

* R3 and R4 accept a nightly row at commit C iff ``reusable(C, HEAD)``. A reused row (C != HEAD)
  says ``reused from <first 7 hex of C> (reusable: diff ⊆ NON_STATS_PATHS)`` (value or evidence);
  a row for HEAD itself does not say "reused from".
* H4 additionally recomputes ``gippa.io.image_key.key(tool, <repo>)`` for the GHCR tools: every
  directory under ``workflow/containers/`` except ``uniop`` (unpublished), as in
  ``scripts/release/image_key.py``. H4 fails when
  - an ``images.json`` entry's key differs from the recomputed one (``stale image: <tool>``),
  - a GHCR tool has no ``images.json`` entry (the 14a ``--check`` passes such a tool; H4 must not),
  - an entry has ``"processes": []`` and the tool is not ``fgsea``,
  - an entry names a tool with no image inputs, or ``images.json`` is not valid JSON;
  each reason names the tool, and H4 never raises. ``fgsea`` is pinned in ``images.json`` only
  and needs no ``withName`` line. With no GHCR directory and no ``images.json`` H4 is unchanged.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pytest
from helpers.release_readiness import (
    COPIED,
    REPO,
    VERSION,
    Fixture,
    git,
    load_module,
    passing_nightly,
    row,
    text_of,
)

from gippa.io.image_key import key

REUSED = "(reusable: diff ⊆ NON_STATS_PATHS)"


@pytest.fixture
def fx(tmp_path: Path) -> Fixture:
    return Fixture.build(tmp_path)


def _calib(fx: Fixture) -> tuple[dict[str, Any], dict[str, Any]]:
    _, doc = fx.compute()
    return row(doc, "R3"), row(doc, "R4")


# ------------------------------------------------------------------ calibration rows (R3, R4)


def test_rows_for_head_itself_are_not_marked_reused(fx: Fixture) -> None:
    r3, r4 = _calib(fx)
    assert (r3["status"], r4["status"]) == ("pass", "pass")
    assert "reused from" not in text_of(r3) + text_of(r4)


def test_example12_viz_docs_version_line_reuses_the_rows(fx: Fixture) -> None:
    c = fx.r
    fx.write_text("src/gippa/viz/locus.py", '"""locus plot (fixture)."""\n')
    fx.write_text("docs/x.md", "# x\n")
    fx.write_text("src/gippa/__init__.py", '"""GIPPA."""\n\n__version__ = "1.0.1"\n')
    fx.commit("1.0.1: viz, docs, version")
    for r in _calib(fx):
        assert r["status"] == "pass", r
        assert f"reused from {c[:7]}" in text_of(r), r
        assert REUSED in text_of(r), r


@pytest.mark.parametrize("rel", ["src/gippa/cli.py", "src/gippa/stats/fdr.py"])
def test_example13_cli_or_stats_change_leaves_rows_pending(fx: Fixture, rel: str) -> None:
    fx.write_text("docs/x.md", "# x\n")
    fx.write_text(rel, '"""changed (fixture)."""\n')
    fx.commit(f"touch {rel}")
    for r in _calib(fx):
        assert r["status"] == "pending", r
        assert "reused from" not in text_of(r)


def test_example13b_readiness_change_leaves_rows_pending(fx: Fixture) -> None:
    fx.write_text("scripts/release/readiness.py", "# widened (fixture)\n")
    fx.write_text("docs/x.md", "# x\n")
    fx.commit("edit readiness")
    assert [r["status"] for r in _calib(fx)] == ["pending", "pending"]


def test_example13d_report_stats_is_not_under_report(fx: Fixture) -> None:
    fx.write_text("src/gippa/report_stats.py", '"""not report/ (fixture)."""\n')
    fx.commit("report_stats")
    assert [r["status"] for r in _calib(fx)] == ["pending", "pending"]


@pytest.mark.parametrize(
    "rel", ["workflow/conf/base.config", "benchmarks/baseline_extra.json", "uv.lock"]
)
def test_paths_outside_src_but_not_non_stats_leave_rows_pending(fx: Fixture, rel: str) -> None:
    """The M11-03 rule reused rows whenever CALIB_PATHS were unchanged; 14g reuses only diffs
    inside NON_STATS_PATHS, so a workflow/ (never list) or other change pays the rerun."""
    fx.write_text(rel, "changed\n")
    fx.commit(f"touch {rel}")
    assert [r["status"] for r in _calib(fx)] == ["pending", "pending"]


def test_a_move_out_of_stats_into_docs_leaves_rows_pending(fx: Fixture) -> None:
    git(fx.repo, "mv", "src/gippa/stats/bh.py", "docs/bh.py")
    git(fx.repo, "commit", "-q", "-m", "move bh into docs")
    assert [r["status"] for r in _calib(fx)] == ["pending", "pending"]


def test_report_change_reuses_and_cites_the_evidence_commit(fx: Fixture) -> None:
    c = fx.r
    fx.write_text("src/gippa/report/card.py", '"""card (fixture)."""\n')
    fx.commit("report card")
    r3, _ = _calib(fx)
    assert r3["status"] == "pass", r3
    assert any(fx.short(c) in e for e in r3["evidence"]), r3["evidence"]
    assert f"reused from {c[:7]}" in text_of(r3)


def test_evidence_from_a_non_ancestor_commit_is_not_reused(fx: Fixture) -> None:
    """Review #385: a nightly row logged on a side commit whose diff to HEAD is docs-only still
    leaves R3/R4 pending, because C must be an ancestor of R."""
    base = fx.r
    git(fx.repo, "checkout", "-q", "-b", "side")
    fx.write_text("docs/side.md", "# side\n")
    side = fx.commit("side: docs")
    git(fx.repo, "checkout", "-q", base)
    git(fx.repo, "checkout", "-q", "-B", "release")
    fx.write_text("docs/side.md", "# side\n")
    fx.write_text("docs/x.md", "# x\n")
    fx.commit("release: docs")
    fx.set_nightly(passing_nightly(fx.short(side)))
    for r in _calib(fx):
        assert r["status"] == "pending", r
        assert "reused from" not in text_of(r)


# ------------------------------------------------------------------ H4 image keys

OWNER = "ghcr.io/khoa-yelo/gippa"
D1, D2 = "1" * 64, "2" * 64


def _tree(tmp_path: Path) -> Path:
    """A directory with the H4 inputs of the checkout (containers.config, tools.py)."""
    root = tmp_path / "tree"
    for rel in COPIED:
        dst = root / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / rel, dst)
    return root


def _write(root: Path, rel: str, text: str) -> None:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text, encoding="utf-8")


def _ostir(root: Path, *, pin: bool = True) -> None:
    _write(
        root,
        "workflow/containers/ostir/Dockerfile",
        "FROM scratch\nCOPY workflow/envs/ostir.yml /x\n",
    )
    _write(root, "workflow/envs/ostir.yml", "name: gippa-ostir\ndependencies: [ostir=1.1.3]\n")
    if pin:
        cfg = root / "workflow/conf/containers.config"
        cfg.write_text(
            cfg.read_text(encoding="utf-8")
            + "\nprocess {\n"
            + f"  withName: 'OSTIR' {{ container = '{OWNER}-ostir@sha256:{D1}' }}\n"
            + f"  withName: 'OSTIR_NOVEL' {{ container = '{OWNER}-ostir@sha256:{D1}' }}\n}}\n",
            encoding="utf-8",
        )


def _fgsea(root: Path) -> None:
    _write(
        root,
        "workflow/containers/fgsea/Dockerfile",
        "FROM scratch\nCOPY workflow/bin/gippa_fgsea.R /opt/gippa/\n",
    )
    _write(root, "workflow/bin/gippa_fgsea.R", "#!/usr/bin/env Rscript\ncat('fgsea')\n")


def _entry(root: Path, tool: str, processes: list[str], k: str | None = None) -> dict[str, Any]:
    k = key(tool, root) if k is None else k
    digest = D1 if tool == "ostir" else D2
    return {"key": k, "tag": f"1.0-{k}", "ref": f"{OWNER}-{tool}@sha256:{digest}",
            "processes": processes, "built_from": "0" * 40, "licence": "MIT"}  # fmt: skip


def _images(root: Path, doc: dict[str, Any]) -> None:
    _write(root, "workflow/conf/images.json", json.dumps(doc, indent=2) + "\n")


def _h4(root: Path) -> dict[str, Any]:
    rd = load_module()
    out: dict[str, Any] = rd.h4(rd.Ctx(root, VERSION))
    return out


def test_h4_unchanged_without_ghcr_tools(tmp_path: Path) -> None:
    assert _h4(_tree(tmp_path))["status"] == "pass"


def test_h4_passes_with_current_keys(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _ostir(root)
    _fgsea(root)
    _images(root, {"fgsea": _entry(root, "fgsea", []),
                   "ostir": _entry(root, "ostir", ["OSTIR", "OSTIR_NOVEL"])})  # fmt: skip
    h4 = _h4(root)
    assert h4["status"] == "pass", h4


def test_h4_fails_on_a_stale_key(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _ostir(root)
    _images(root, {"ostir": _entry(root, "ostir", ["OSTIR", "OSTIR_NOVEL"])})
    _write(root, "workflow/envs/ostir.yml", "name: gippa-ostir\ndependencies: [ostir=1.1.4]\n")
    h4 = _h4(root)
    assert h4["status"] == "fail", h4
    assert "stale image: ostir" in h4["value"], h4


def test_h4_fails_on_a_stale_fgsea_script(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _fgsea(root)
    _images(root, {"fgsea": _entry(root, "fgsea", [])})
    _write(root, "workflow/bin/gippa_fgsea.R", "#!/usr/bin/env Rscript\ncat('fgsea v2')\n")
    h4 = _h4(root)
    assert h4["status"] == "fail" and "stale image: fgsea" in h4["value"], h4


@pytest.mark.parametrize("pin", [True, False], ids=["withName-pin", "conda-only"])
def test_h4_fails_for_a_ghcr_tool_without_an_images_entry(tmp_path: Path, pin: bool) -> None:
    """The 14a review: ``image_key.py --check`` passes an unpinned GHCR tool; H4 must not."""
    root = _tree(tmp_path)
    _ostir(root, pin=pin)
    _images(root, {})
    h4 = _h4(root)
    assert h4["status"] == "fail", h4
    assert "ostir" in h4["value"], h4


def test_h4_fails_for_a_ghcr_tool_when_images_json_is_absent(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _ostir(root)
    h4 = _h4(root)
    assert h4["status"] == "fail" and "ostir" in h4["value"], h4


def test_h4_accepts_a_process_less_entry_for_fgsea_only(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _ostir(root, pin=False)
    _images(root, {"ostir": _entry(root, "ostir", [])})
    h4 = _h4(root)
    assert h4["status"] == "fail", h4
    assert "ostir" in h4["value"], h4


def test_h4_fgsea_needs_no_withname_line(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _fgsea(root)
    _images(root, {"fgsea": _entry(root, "fgsea", [])})
    cfg = (root / "workflow/conf/containers.config").read_text(encoding="utf-8")
    assert "FGSEA" not in cfg
    h4 = _h4(root)
    assert h4["status"] == "pass", h4


def test_h4_fails_on_an_entry_without_image_inputs(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _images(root, {"ghost": _entry(root, "ghost", ["GHOST"], k="0" * 12)})
    h4 = _h4(root)
    assert h4["status"] == "fail" and "ghost" in h4["value"], h4


def test_h4_fails_on_unreadable_images_json(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _ostir(root)
    _write(root, "workflow/conf/images.json", "{not json\n")
    h4 = _h4(root)
    assert h4["status"] == "fail", h4
    assert "images.json" in text_of(h4), h4


def test_h4_on_the_checkout_names_each_ghcr_tool_not_current(tmp_path: Path) -> None:
    """Holds before and after 14j: every GHCR tool lacking a current images.json entry is named;
    uniop (never published) is never asked for one."""
    rd = load_module()
    h4 = rd.h4(rd.Ctx(REPO, VERSION))
    images_path = REPO / "workflow/conf/images.json"
    images = json.loads(images_path.read_text(encoding="utf-8")) if images_path.is_file() else {}
    tools = sorted(p.name for p in (REPO / "workflow/containers").iterdir()
                   if p.is_dir() and p.name not in {"uniop", "step", "step-r"})  # fmt: skip
    behind = [t for t in tools if (images.get(t) or {}).get("key") != key(t, REPO)]
    for t in behind:
        assert t in h4["value"], (t, h4)
    if behind:
        assert h4["status"] == "fail", h4
    assert not any("uniop" in part for part in h4["value"].split("; ")), h4


def test_h4_evidence_names_images_json(tmp_path: Path) -> None:
    root = _tree(tmp_path)
    _ostir(root)
    _images(root, {"ostir": _entry(root, "ostir", ["OSTIR", "OSTIR_NOVEL"])})
    assert "workflow/conf/images.json" in _h4(root)["evidence"]
