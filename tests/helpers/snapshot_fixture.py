"""The M9-25 snapshot fixtures, plot cases and report-structure snapshots
(docs/specs/M9-25-snapshot-baselines.md work items 1-3; GIAP_SPEC.md §10.1, §10.3).

Layout pinned by the test-author (the implementer's ``scripts/dev/make_snapshot_fixture.py``
writes it; the tests read it):

* ``tests/data/snapshot/<fixture>/OUT/``: the ``OUT/`` subset ``RunView`` reads, with its
  ``resolved_config.yaml`` (set ``report.max_locus_plots: 2`` there to keep renders short).
  ``FIXTURES = ("sim_plots", "sim_empty")``: the rich giapsim run and the header-only tree.
* ``tests/data/snapshot/<fixture>/CMP/``: the V16 / ``compare.html`` input, a frozen
  ``giap compare`` subset: ``CMP/tables/compare_{detection,abundance,stability}.tsv`` and
  ``CMP/runs/<label>/`` (two or more runs; each with ``tables/genes.tsv`` and
  ``preflight/preflight.json``, the same genome SHA-256). Labels are the sorted directory names.
* ``tests/snapshot/baseline/<fixture>-<stem>.png``: one pytest-mpl baseline per ``Case``.
* ``tests/snapshot/report_structure.json``: ``{"<fixture>/report.html": [[tag, id], ...],
  "<fixture>/compare.html": [...]}`` for both fixtures; regenerate with
  ``python -m helpers.snapshot_fixture --write-structure`` (from ``tests/``, after Khoa's
  approval only).

Plot cases: every instance of every registered plot on ``OUT/`` (V16 from ``CMP/`` through
``giap.viz.compare.compare_plot``), except that V03 and V13 keep only their first two locus
instances (V03: suffix not ``None``; V13: suffix ``locus_<id>``); the index and panel (a)
instances stay. Figures are drawn under ``matplotlib.style.context("default")`` +
``giap.viz.style.apply()``, which is what the ``plots`` step sees and what the tests' pytest-mpl
marker (``style="default"``) restores at save time.
"""

from __future__ import annotations

import contextlib
import functools
import json
import os
import shutil
import sys
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from html.parser import HTMLParser
from importlib import resources
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
SNAPSHOT = REPO / "tests" / "data" / "snapshot"
FIXTURES = ("sim_plots", "sim_empty")
RICH, EMPTY = FIXTURES
BASELINE_DIR = REPO / "tests" / "snapshot" / "baseline"
STRUCTURE_JSON = REPO / "tests" / "snapshot" / "report_structure.json"
SNAPSHOT_DPI = 100
TOLERANCE = 2
LOCUS_KEEP = 2
COMPARE_TABLES = ("compare_detection", "compare_abundance", "compare_stability")
MAX_FIXTURE_BYTES = 2_000_000  # spec work item 1: "< 2 MB" (both fixtures together)
MAX_FILE_BYTES = 5_000_000  # R11
MAX_BASELINE_BYTES = 300_000  # spec work item 2
#: A directory shared by every xdist worker of one pytest session, set by
#: ``tests/snapshot/conftest.py`` (``None`` elsewhere: scripts and other test layers render
#: afresh on every call). With it, ``raw_results`` is computed once per process and
#: ``render_plot_files`` once per session (file-locked; later calls copy the rendered tree).
SESSION_CACHE: Path | None = None
#: environment variable naming the directory of the session's case lists (see ``_case_keys``)
CASES_ENV = "GIAP_SNAPSHOT_CASES"


@dataclass(frozen=True)
class Case:
    """One snapshot: ``fixture`` x one instance ``stem`` of plot ``pid`` (or a setup error)."""

    fixture: str
    pid: str
    stem: str
    error: str | None = None

    @property
    def id(self) -> str:
        return f"{self.fixture}-{self.stem}"

    @property
    def filename(self) -> str:
        return f"{self.fixture}-{self.stem}.png"


def fixture_out(fixture: str) -> Path:
    return SNAPSHOT / fixture / "OUT"


def fixture_cmp(fixture: str) -> Path:
    return SNAPSHOT / fixture / "CMP"


def missing(fixture: str) -> str | None:
    """Why ``fixture`` cannot be used, or ``None``."""
    out, cmp = fixture_out(fixture), fixture_cmp(fixture)
    if not ((out / "tables").is_dir() or (out / "preflight").is_dir()):
        return f"M9-25 snapshot fixture {out.relative_to(REPO)} is not committed"
    if not (cmp / "tables").is_dir() or not (cmp / "runs").is_dir():
        return f"M9-25 compare fixture {cmp.relative_to(REPO)} (tables/, runs/) is not committed"
    return None


def fixture_cfg(fixture: str) -> Any:
    from giap.config import resolve_config

    out = fixture_out(fixture)
    runs = (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml")
    return resolve_config(next((p for p in runs if p.is_file()), None))


def styled(fn: Callable[[], Any]) -> Any:
    """``fn()`` under matplotlib's default style plus ``giap.viz.style.apply()``."""
    import matplotlib.style

    from giap.viz import style

    with matplotlib.style.context("default", after_reset=True):
        style.apply()
        return fn()


def compare_inputs(fixture: str) -> tuple[list[Path], Any, dict[str, Any]]:
    """``(run paths, LoadedRuns, compare tables)`` of ``CMP/`` (typed by the table schemas)."""
    from giap.io.compare_io import load_runs
    from giap.io.run_view import read_schema_tsv

    cmp = fixture_cmp(fixture)
    paths = sorted(p for p in (cmp / "runs").iterdir() if p.is_dir())
    loaded = load_runs(paths, None)
    tables = {}
    for name in COMPARE_TABLES:
        res = resources.files("giap.schemas").joinpath("tables", f"{name}.schema.json")
        schema = json.loads(res.read_text(encoding="utf-8"))
        tables[name] = read_schema_tsv(cmp / "tables" / f"{name}.tsv", name, schema)
    return paths, loaded, tables


def is_locus(pid: str, suffix: str | None) -> bool:
    if pid == "V03":
        return suffix is not None
    if pid == "V13":
        return (suffix or "").startswith("locus_")
    return False


def keep_instances(pid: str, insts: list[Any]) -> list[Any]:
    """All non-locus instances and the first ``LOCUS_KEEP`` locus instances, in order."""
    out, n = [], 0
    for inst in insts:
        if is_locus(pid, inst.suffix):
            n += 1
            if n > LOCUS_KEEP:
                continue
        out.append(inst)
    return out


def raw_results(fixture: str) -> dict[str, Any]:
    """``pid -> list[Instance] | Skip`` of every registered plot (V16 from ``CMP/``); call
    inside ``styled``. Under ``SESSION_CACHE`` the result is computed once per process (the
    callers only read the instances)."""
    if SESSION_CACHE is not None:
        return _raw_results_cached(fixture)
    return _raw_results(fixture)


@functools.cache
def _raw_results_cached(fixture: str) -> dict[str, Any]:
    return _raw_results(fixture)


def _raw_results(fixture: str, pids: tuple[str, ...] = ()) -> dict[str, Any]:
    """``raw_results`` of every registered plot, or only of ``pids``."""
    from giap.io.run_view import RunView
    from giap.steps.plots import load_plot_modules
    from giap.viz import base
    from giap.viz.compare import compare_plot

    load_plot_modules()
    cfg = fixture_cfg(fixture)
    view = RunView.from_outdir(fixture_out(fixture))
    out: dict[str, Any] = {}
    for pid, spec in base.registered().items():
        if pids and pid not in pids:
            continue
        if pid == "V16":
            _, loaded, tables = compare_inputs(fixture)
            det, ab = tables["compare_detection"], tables["compare_abundance"]
            out[pid] = [compare_plot(det, ab, loaded.genes)]
        else:
            out[pid] = spec.fn(view, cfg)
    return out


def _snapshot_instances(results: dict[str, Any]) -> dict[tuple[str, str], Any]:
    from giap.io.plot_files import stem
    from giap.viz import base

    found: dict[tuple[str, str], Any] = {}
    for pid, result in results.items():
        if isinstance(result, base.Skip):
            continue
        name = base.registered()[pid].name
        for inst in keep_instances(pid, list(result)):
            found[(pid, stem(pid, name, inst.suffix))] = inst
    return found


#: fixture -> ``rendered(fixture)``, once computed in this process
_RENDERED: dict[str, dict[tuple[str, str], Any]] = {}


def rendered(fixture: str) -> dict[tuple[str, str], Any]:
    """``(pid, stem) -> Instance`` of the snapshot cases of ``fixture`` (cached per process)."""
    if fixture not in _RENDERED:
        _RENDERED[fixture] = styled(lambda: _snapshot_instances(raw_results(fixture)))
    return _RENDERED[fixture]


@functools.cache
def _rendered_pid(fixture: str, pid: str) -> dict[tuple[str, str], Any]:
    """``rendered(fixture)`` restricted to plot ``pid``, drawing only that plot."""
    out: dict[tuple[str, str], Any] = styled(
        lambda: _snapshot_instances(_raw_results(fixture, (pid,)))
    )
    return out


def _case_keys(fixture: str) -> list[tuple[str, str]]:
    """The ``(pid, stem)`` keys of ``rendered(fixture)``. Under ``$GIAP_SNAPSHOT_CASES`` (set by
    ``tests/snapshot/conftest.py`` for the session and inherited by the subprocesses the tests
    start) the first process to render writes them there and later processes read them."""
    where = os.environ.get(CASES_ENV)
    path = Path(where) / f"cases-{fixture}.json" if where else None
    if path is not None and path.is_file() and fixture not in _RENDERED:
        return [(pid, s) for pid, s in json.loads(path.read_text(encoding="utf-8"))]
    keys = list(rendered(fixture))
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(f"{path.name}.{os.getpid()}.part")
        part.write_text(json.dumps(keys), encoding="utf-8")
        os.replace(part, path)
    return keys


def cases() -> list[Case]:
    """Every snapshot case of both fixtures; a fixture that is absent or fails to render gives
    one error case, so the parametrisation never silently shrinks."""
    out: list[Case] = []
    for fixture in FIXTURES:
        why = missing(fixture)
        if why is None:
            try:
                keys = _case_keys(fixture)
            except Exception as exc:  # reported by the error case, never at collection
                why = f"rendering {fixture} raised {type(exc).__name__}: {exc}"
        if why is not None:
            out.append(Case(fixture, "-", "setup", why))
            continue
        out += [Case(fixture, pid, s) for pid, s in keys]
    return out


def figure(case: Case) -> Any:
    """The figure of ``case``. A process that took its cases from ``$GIAP_SNAPSHOT_CASES``
    without rendering, outside the snapshot session (a test's subprocess), draws only
    ``case.pid``."""
    key = (case.pid, case.stem)
    if case.fixture in _RENDERED or SESSION_CACHE is not None:
        return rendered(case.fixture)[key].figure
    return _rendered_pid(case.fixture, case.pid)[key].figure


# ---------------------------------------------------------------- report and compare pages


@contextlib.contextmanager
def _locked(path: Path) -> Iterator[None]:
    """An exclusive ``flock`` on ``path`` (no lock where ``fcntl`` is absent)."""
    try:
        import fcntl
    except ImportError:  # pragma: no cover - not POSIX
        yield
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)


def render_plot_files(fixture: str, plots_dir: Path, formats: tuple[str, ...]) -> None:
    """Every instance (no locus cap) of every registered plot of ``OUT/`` into ``plots_dir``.

    Under ``SESSION_CACHE`` the first call of the session renders into the cache (one worker at
    a time, under a lock) and every call copies that tree: the renders are byte-identical
    (``test_snapshot_byte_identity.py``), so the copy is what a fresh render would write."""
    if SESSION_CACHE is None:
        _render_plot_files(fixture, plots_dir, formats)
        return
    done = SESSION_CACHE / f"plots-{fixture}-{'-'.join(formats)}"
    with _locked(SESSION_CACHE / f"{done.name}.lock"):
        if not done.is_dir():
            part = SESSION_CACHE / f"{done.name}.part"
            shutil.rmtree(part, ignore_errors=True)
            _render_plot_files(fixture, part, formats)
            part.rename(done)
    shutil.copytree(done, plots_dir, dirs_exist_ok=True)


def _render_plot_files(fixture: str, plots_dir: Path, formats: tuple[str, ...]) -> None:
    from giap.io.plot_files import write_instance
    from giap.io.run_view import RunView
    from giap.steps.plots import load_plot_modules
    from giap.viz import base

    def run() -> None:
        load_plot_modules()
        cfg = fixture_cfg(fixture)
        view = RunView.from_outdir(fixture_out(fixture))
        for pid, spec in base.registered().items():
            result = spec.fn(view, cfg)
            if isinstance(result, base.Skip):
                continue
            for inst in result:
                write_instance(inst, pid, spec.name, plots_dir, list(formats), SNAPSHOT_DPI)

    styled(run)


def build_report(fixture: str, tmp: Path) -> Path:
    """``report.html`` of ``fixture`` via ``giap.steps.plots.write_report`` (plots under
    ``tmp/plots``); returns the page path."""
    from giap.io.run_view import RunView
    from giap.steps.plots import write_report

    plots = tmp / "plots"
    render_plot_files(fixture, plots, ("svg", "png"))
    view = RunView.from_outdir(fixture_out(fixture))
    return write_report(view, plots, tmp / "report", fixture_cfg(fixture))


def build_compare(fixture: str, tmp: Path) -> Path:
    """``compare.html`` of ``fixture``'s ``CMP/`` through the ``giap compare`` report tail
    (``giap.cli._compare_report``) on a copy under ``tmp``; returns the page path."""
    from giap.cli import _compare_report

    cmp = tmp / "CMP"
    shutil.copytree(fixture_cmp(fixture), cmp)
    src_paths, _, tables = compare_inputs(fixture)
    paths = [cmp / "runs" / p.name for p in src_paths]
    from giap.io.compare_io import load_runs

    loaded = load_runs(paths, None)
    styled(lambda: _compare_report(cmp, paths, loaded, tables, fixture_cfg(fixture)))
    return cmp / "report" / "compare.html"


# ---------------------------------------------------------------- HTML structure


class _Structure(HTMLParser):
    """Collects the §10.3 structure of a GIAP page, in document order:

    * ``["section", id]`` for every ``<section>``;
    * ``["h2", id-or-text]`` and ``["h3", id-or-text]`` (the ``id`` attribute when present, else
      the heading text with whitespace collapsed);
    * ``["div", "plot:<data-plot>"]`` for every ``<div class="plot">``, ``["div", id]`` for
      every other ``<div>`` with an ``id`` (the Plotly twins);
    * ``["svg", id-or-null]`` for every outermost ``<svg>`` (an inlined figure);
    * ``["figure", stem]`` for every ``<figcaption>`` (its first word, the instance stem).
    """

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.items: list[list[str | None]] = []
        self._svg = 0
        self._text: list[str] | None = None
        self._open: tuple[str, str | None] | None = None

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        a = dict(attrs)
        if tag == "svg":
            if self._svg == 0:
                self.items.append(["svg", a.get("id")])
            self._svg += 1
            return
        if self._svg:
            return
        if tag == "section":
            self.items.append(["section", a.get("id")])
        elif tag == "div" and "plot" in (a.get("class") or "").split():
            self.items.append(["div", f"plot:{a.get('data-plot')}"])
        elif tag == "div" and a.get("id"):
            self.items.append(["div", a.get("id")])
        elif tag in ("h2", "h3", "figcaption"):
            self._open, self._text = (tag, a.get("id")), []

    def handle_startendtag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "svg":
            if self._svg == 0:
                self.items.append(["svg", dict(attrs).get("id")])
            return
        self.handle_starttag(tag, attrs)

    def handle_endtag(self, tag: str) -> None:
        if tag == "svg":
            self._svg = max(0, self._svg - 1)
            return
        if self._open is not None and tag == self._open[0] and self._text is not None:
            text = " ".join("".join(self._text).split())
            name, ident = self._open
            if name == "figcaption":
                self.items.append(["figure", text.split(" ")[0] if text else ""])
            else:
                self.items.append([name, ident if ident else text])
            self._open, self._text = None, None

    def handle_data(self, data: str) -> None:
        if self._text is not None and not self._svg:
            self._text.append(data)


def structure(html: str) -> list[list[str | None]]:
    """The ordered ``[tag, id]`` list of ``html`` (see ``_Structure``)."""
    parser = _Structure()
    parser.feed(html)
    parser.close()
    return parser.items


def structures(tmp: Path) -> dict[str, list[list[str | None]]]:
    """The four snapshot entries: ``<fixture>/report.html`` and ``<fixture>/compare.html``."""
    out: dict[str, list[list[str | None]]] = {}
    for fixture in FIXTURES:
        page = build_report(fixture, tmp / fixture / "run")
        out[f"{fixture}/report.html"] = structure(page.read_text(encoding="utf-8"))
        page = build_compare(fixture, tmp / fixture / "cmp")
        out[f"{fixture}/compare.html"] = structure(page.read_text(encoding="utf-8"))
    return out


def _main(argv: list[str]) -> int:
    import tempfile

    if argv != ["--write-structure"]:
        print("usage: python -m helpers.snapshot_fixture --write-structure", file=sys.stderr)
        return 2
    with tempfile.TemporaryDirectory() as tmp:
        data = structures(Path(tmp))
    text = json.dumps(data, indent=1, sort_keys=True) + "\n"
    STRUCTURE_JSON.write_text(text, encoding="utf-8")
    print(f"wrote {STRUCTURE_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
