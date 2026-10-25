"""M9-25 work item 5: ``scripts/dev/approval_packet.py`` writes the M9 approval packet for Khoa
(docs/specs/M9-25-snapshot-baselines.md; M9-index "Khoa's checkpoint"; GIAP_SPEC.md §14.8 item 3).

Interface (pinned by the test-author; the implementer may add files and options, not rename
these). ``python scripts/dev/approval_packet.py [--dest DIR] [--baselines DIR] [--k12 OUT]``,
with ``main(argv: list[str]) -> int``:

* ``--dest`` defaults to ``local/review/m9-approval``; ``--baselines`` to
  ``tests/snapshot/baseline``; the fixtures are ``tests/data/snapshot/{sim_plots,sim_empty}``.
* ``DEST/baseline/<name>.png``: a byte copy of every baseline PNG;
  ``DEST/deuteranopia/<name>.png``: ``giap.viz.checklist.simulate_deuteranopia`` of its RGB
  (alpha kept), same shape.
* with ``--k12 OUT`` (the M9-26 K-12 run): ``DEST/k12/plots/<stem>.png`` for every instance of
  every registered plot that does not skip on ``OUT`` (100 dpi) and
  ``DEST/k12/deuteranopia/<stem>.png`` for each. Without it the index says the K-12 run was not
  supplied; a ``--k12`` path that is not a directory exits 2 and writes no ``index.html``.
* ``DEST/checklist.tsv``: header ``fixture plot item passed detail``, the M9-24 checklist of every
  instance on each fixture (and ``k12``); ``DEST/score.json``: ``{fixture: score}`` with
  ``checklist.score`` over that fixture's results (the score at approval).
* a ``report.html`` (rich fixture) and a ``compare.html`` somewhere under ``DEST``.
* ``DEST/index.html``: relative links (``href`` or ``src``) to every PNG under ``baseline/``,
  ``deuteranopia/`` and ``k12/``, to ``checklist.tsv``, the report and compare pages; it states
  each fixture's score. No absolute path anywhere in it.
"""

from __future__ import annotations

import csv
import importlib.util
import json
import re
import shutil
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from helpers import snapshot_fixture as sf

SCRIPT = sf.REPO / "scripts" / "dev" / "approval_packet.py"


def _module() -> ModuleType:
    assert SCRIPT.is_file(), "scripts/dev/approval_packet.py (M9-25 work item 5) is missing"
    spec = importlib.util.spec_from_file_location("approval_packet_m9_25", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _baselines(tmp: Path) -> Path:
    """Two small stand-in baseline PNGs (the real ones are Khoa's to approve)."""
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    d = tmp / "baselines"
    d.mkdir()
    for name, colour in (("sim_plots-V07_obs_vs_exp", "#0072B2"), ("sim_empty-V07_obs_vs_exp",
                                                                   "#E69F00")):  # fmt: skip
        fig = Figure(figsize=(1.0, 0.8))
        ax = fig.add_subplot()
        ax.bar([0, 1], [1, 2], color=[colour, "#009E73"])
        fig.savefig(d / f"{name}.png", dpi=100)
    return d


@pytest.fixture(scope="module")
def packet(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    for fixture in sf.FIXTURES:
        why = sf.missing(fixture)
        if why is not None:
            pytest.fail(why)
    mod = _module()
    tmp = tmp_path_factory.mktemp("m9_approval")
    baselines = _baselines(tmp)
    k12 = tmp / "k12" / "OUT"
    shutil.copytree(sf.fixture_out(sf.RICH), k12)  # a stand-in for the M9-26 K-12 run
    dest = tmp / "packet"
    code = mod.main(["--dest", str(dest), "--baselines", str(baselines), "--k12", str(k12)])
    assert code == 0
    return dest, baselines


def _links(index: Path) -> set[str]:
    text = index.read_text(encoding="utf-8")
    return set(re.findall(r'(?:href|src)="([^"#]+)"', text))


def test_baselines_copied_byte_for_byte(packet: tuple[Path, Path]) -> None:
    dest, baselines = packet
    for png in sorted(baselines.iterdir()):
        assert (dest / "baseline" / png.name).read_bytes() == png.read_bytes()


def test_deuteranopia_simulation_of_every_baseline(packet: tuple[Path, Path]) -> None:
    from giap.io.plot_files import read_png
    from giap.viz.checklist import simulate_deuteranopia

    dest, baselines = packet
    for png in sorted(baselines.iterdir()):
        src, sim = read_png(png), read_png(dest / "deuteranopia" / png.name)
        assert sim.shape == src.shape
        want = simulate_deuteranopia(np.asarray(src[..., :3], dtype=float))
        assert np.allclose(np.asarray(sim[..., :3], dtype=float), want, atol=2.5 / 255)
        if src.shape[-1] == 4:
            assert np.allclose(sim[..., 3], src[..., 3], atol=1 / 255)


def test_k12_plots_and_their_simulations(packet: tuple[Path, Path]) -> None:
    from giap.io.plot_files import stem
    from giap.io.run_view import RunView
    from giap.steps.plots import load_plot_modules
    from giap.viz import base

    dest, _ = packet
    load_plot_modules()
    cfg = sf.fixture_cfg(sf.RICH)
    view = RunView.from_outdir(sf.fixture_out(sf.RICH))

    def stems() -> set[str]:
        out = set()
        for pid, spec in base.registered().items():
            res = spec.fn(view, cfg)
            if not isinstance(res, base.Skip):
                out |= {stem(pid, spec.name, i.suffix) for i in res}
        return out

    want = sf.styled(stems)
    plots = {p.stem for p in (dest / "k12" / "plots").glob("*.png")}
    assert want <= plots
    sims = {p.stem for p in (dest / "k12" / "deuteranopia").glob("*.png")}
    assert plots <= sims


def test_checklist_table_and_score(packet: tuple[Path, Path]) -> None:
    dest, _ = packet
    with (dest / "checklist.tsv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert list(rows[0]) == ["fixture", "plot", "item", "passed", "detail"]
    assert {r["fixture"] for r in rows} == {*sf.FIXTURES, "k12"}
    score = json.loads((dest / "score.json").read_text(encoding="utf-8"))
    assert set(score) >= set(sf.FIXTURES)
    for fixture in sf.FIXTURES:
        gated = [r["passed"] for r in rows if r["fixture"] == fixture and r["passed"] != "NA"]
        gated = [g for g in gated if g in ("true", "false")]
        expect = sum(g == "true" for g in gated) / len(gated) if gated else 1.0
        assert score[fixture] == pytest.approx(expect)


def test_index_links_everything_relatively(packet: tuple[Path, Path]) -> None:
    dest, _ = packet
    index = dest / "index.html"
    links = _links(index)
    for link in links:
        assert not link.startswith(("/", "file:", "http")), link
    resolved = {(dest / link).resolve() for link in links}
    for sub in ("baseline", "deuteranopia", "k12"):
        for png in sorted((dest / sub).rglob("*.png")):
            assert png.resolve() in resolved, png.relative_to(dest)
    assert (dest / "checklist.tsv").resolve() in resolved
    pages = {p.name for p in resolved if p.is_file() and p.suffix == ".html"}
    assert {"report.html", "compare.html"} <= pages
    text = index.read_text(encoding="utf-8")
    assert str(dest) not in text and str(sf.REPO) not in text
    score = json.loads((dest / "score.json").read_text(encoding="utf-8"))
    for fixture in sf.FIXTURES:
        assert fixture in text
        assert f"{score[fixture]:.2f}" in text


def test_k12_path_not_a_directory_exits_2(tmp_path: Path) -> None:
    mod = _module()
    dest = tmp_path / "packet"
    code = mod.main(["--dest", str(dest), "--baselines", str(_baselines(tmp_path)),
                     "--k12", str(tmp_path / "absent")])  # fmt: skip
    assert code == 2
    assert not (dest / "index.html").exists()


def test_default_destination_is_local_review() -> None:
    mod = _module()
    assert Path(mod.DEFAULT_DEST) == sf.REPO / "local" / "review" / "m9-approval"
    assert Path(mod.DEFAULT_BASELINES) == sf.BASELINE_DIR
