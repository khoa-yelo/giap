"""``plots`` and ``report`` in ``gippa run`` under both executors (M9-02 behaviour 8; acceptance
"the registry key test for ``plots`` / ``report`` (as ``test_nf_annotate_registry.py``)").

Local executor: ``PUBLISHED`` gains ``plots`` -> ``plots/`` and ``report`` -> ``report/``; both
steps run after every producer of a published key. Nextflow: ``workflow/subworkflows/local/
report.nf`` defines ``REPORT`` with ``GIPPA_STEP as PLOTS`` and ``GIPPA_STEP as REPORT_HTML``, whose
out keys equal the registry outputs; ``main.nf`` includes ``REPORT`` and calls it last. The
nf-test fixture ``workflow/tests/data/annotate/step_registry.json`` lists both steps (it is read
by workflow/tests/report.nf.test). No Nextflow needed.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from gippa.runner.local import PUBLISHED, PUBLISHED_FILES, publish
from gippa.steps import load_step_modules
from gippa.steps import registry as reg

REPO = Path(__file__).resolve().parents[2]
REPORT_NF = REPO / "workflow" / "subworkflows" / "local" / "report.nf"
MAIN_NF = REPO / "workflow" / "main.nf"
FIXTURE = REPO / "workflow" / "tests" / "data" / "annotate" / "step_registry.json"
ALIASES = {"PLOTS": "plots", "REPORT_HTML": "report"}


def _strip(text: str) -> str:
    return re.sub(r"(?m)//.*$", "", text)


# ---------------------------------------------------------------- local executor


def test_published_gains_plots_and_report_directories() -> None:
    assert PUBLISHED["plots"] == "plots"
    assert PUBLISHED["report"] == "report"
    assert "plots" not in PUBLISHED_FILES and "report" not in PUBLISHED_FILES


def test_publish_copies_the_plots_and_report_directories(tmp_path: Path) -> None:
    work, out = tmp_path / "work", tmp_path / "OUT"
    (work / "plots" / "data").mkdir(parents=True)
    (work / "plots" / "V07_obs_vs_exp.svg").write_text("<svg/>")
    (work / "plots" / "data" / "V07_obs_vs_exp.tsv").write_text("x\n")
    (work / "report").mkdir()
    (work / "report" / "report.html").write_text("<!doctype html>")
    publish(work, out)
    assert (out / "plots" / "V07_obs_vs_exp.svg").read_text() == "<svg/>"
    assert (out / "plots" / "data" / "V07_obs_vs_exp.tsv").read_text() == "x\n"
    assert (out / "report" / "report.html").read_text() == "<!doctype html>"


def test_plots_and_report_run_last() -> None:
    from gippa.runner.local import MEMBER_KEYS

    load_step_modules()
    py = [s.name for s in reg.all_steps() if s.kind == "py"]
    order = reg.topo_order(py, {k: parent for k, (parent, _) in MEMBER_KEYS.items()})
    pos = {n: i for i, n in enumerate(order)}
    published = (set(PUBLISHED) | set(PUBLISHED_FILES)) - {"plots", "report"}
    producers = {s.name for s in reg.all_steps() if s.kind == "py" and set(s.outputs) & published}
    assert producers, "no producer of a published key"
    for name in producers:
        assert pos[name] < pos["plots"], name
    assert pos["plots"] < pos["report"]


# ---------------------------------------------------------------- Nextflow


def _calls() -> dict[str, list[str]]:
    """alias -> the ``'key=...'`` out keys in its call ``ALIAS(files, keys, <out spec>)``."""
    text = _strip(REPORT_NF.read_text(encoding="utf-8"))
    out: dict[str, list[str]] = {}
    for m in re.finditer(r"(?<![\w.])(PLOTS|REPORT_HTML)\(", text):
        depth, i = 1, m.end()
        while depth:
            depth += {"(": 1, ")": -1}.get(text[i], 0)
            i += 1
        out[m.group(1)] = re.findall(r"['\"](\w+)=", text[m.end() : i - 1])
    return out


def test_report_nf_defines_report_with_two_gippa_step_aliases() -> None:
    assert REPORT_NF.is_file(), "workflow/subworkflows/local/report.nf does not exist"
    text = _strip(REPORT_NF.read_text(encoding="utf-8"))
    assert re.search(r"\bworkflow\s+REPORT\s*\{", text)
    for alias in ALIASES:
        assert re.search(rf"GIPPA_STEP\s+as\s+{alias}\b", text), alias


def test_report_nf_out_keys_equal_the_registry_outputs() -> None:
    load_step_modules()
    calls = _calls()
    assert set(calls) >= set(ALIASES), f"report.nf calls {sorted(calls)}"
    for alias, step in ALIASES.items():
        assert set(calls[alias]) == set(reg.get(step).outputs), (alias, calls[alias])


def test_main_nf_calls_report_last() -> None:
    text = _strip(MAIN_NF.read_text(encoding="utf-8"))
    assert re.search(r"include\s*\{\s*REPORT\s*\}\s*from\s*'\./subworkflows/local/report'", text)
    calls = [m.start() for m in re.finditer(r"\bREPORT\(", text)]
    assert len(calls) == 1, f"REPORT called {len(calls)}x"
    for other in ("ANNOTATE(", "MAPPING(", "GENETABLE(", "CONTEXT(", "NOVEL(", "OPERONS("):
        assert text.index(other) < calls[0], other


def test_nf_registry_fixture_lists_plots_and_report() -> None:
    load_step_modules()
    got = json.loads(FIXTURE.read_text(encoding="utf-8"))
    for name in ("plots", "report"):
        spec = reg.get(name)
        assert got.get(name) == {
            "inputs": sorted(spec.inputs),
            "optional_inputs": sorted(spec.optional_inputs),
            "outputs": sorted(spec.outputs),
        }, f"regenerate {FIXTURE.relative_to(REPO)} ({name})"
