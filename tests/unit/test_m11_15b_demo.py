"""M11-15b ``gippa demo``: bundled SIM100K data, ``make_demo_layers.py``, the CLI (fast tests).

docs/specs/M11-15-portable-clusters.md item 2 (examples 2.1-2.5, AC 1, 2 and 4). Nothing here
runs the pipeline: ``gippa.cli._run_local`` (the ``--executor local`` path) and
``gippa.cli._run_preflight`` (the first step of the ``gippa run`` Nextflow path, used by
``--full``) are replaced by recorders. The full local-executor demo run and its byte-identity
check against an independent ``gippa run --executor local`` are in
tests/integration/test_m11_15b_demo_run.py (``slow``).

Layout pinned here (spec: "under ``gippa/data/demo/``: ``sim100k.fna.gz``,
``sim100k_peptides.tsv.gz``, the precomputed layers (gzip) and ``PROVENANCE.md``"): the layers
are ``layers/bakta.json.gz`` and ``layers/<key>.tsv.gz`` for every other
``PrecomputedConfig`` field, each the gzip of the file gippasim writes for SIM100K (seed
20261014), written with ``mtime=0``.

No network, no real ``$GIPPA_DB_DIR``, no Nextflow.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import inspect
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
DEMO = REPO / "src" / "gippa" / "data" / "demo"
LAYERS = DEMO / "layers"
GENOME = REPO / "tests" / "data" / "genomes" / "sim100k.fna.gz"
PEPTIDES = REPO / "docs" / "tutorials" / "data" / "sim100k_peptides.tsv.gz"
SCRIPT = REPO / "scripts" / "docs" / "make_demo_layers.py"
CHECK_DIST = REPO / "scripts" / "release" / "check_dist.py"
MIB = 1024 * 1024
SEED = 20261014
BANNER = "Synthetic demo data (SIM100K, gippasim seed 20261014)"
REFUSAL = "exists; use --outdir or remove it"
DONE = re.compile(r"^gippa demo: done in (\d+) s; report: (\S.*)$", re.M)
DEMO_MODULES = ["core", "mobile", "bgc"]


def _layer_keys() -> list[str]:
    from gippa.config.model import PrecomputedConfig

    return [k for k in PrecomputedConfig.model_fields if k != "bakta_json"]


LAYER_KEYS = _layer_keys()


def _layer_file(key: str) -> Path:
    return LAYERS / ("bakta.json.gz" if key == "bakta_json" else f"{key}.tsv.gz")


def _plain(path: Path) -> bytes:
    """File bytes, gunzipped when the file is gzip."""
    data = path.read_bytes()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


# ------------------------------------------------------------------ bundled data (AC 1)


def test_demo_directory_is_package_data() -> None:
    from importlib.resources import files

    root = files("gippa") / "data" / "demo"
    assert root.is_dir(), "gippa/data/demo/ is missing from the package (M11-15 item 2)"
    for name in ("sim100k.fna.gz", "sim100k_peptides.tsv.gz", "PROVENANCE.md"):
        assert (root / name).is_file(), name


def test_wheel_build_includes_the_demo_directory() -> None:
    """The wheel's ``only-include`` covers ``src/gippa``, and nothing excludes ``data/demo``."""
    import tomllib

    hatch = tomllib.loads((REPO / "pyproject.toml").read_text())["tool"]["hatch"]["build"]
    wheel = hatch["targets"]["wheel"]
    assert "src/gippa" in wheel["only-include"]
    excluded = [*wheel.get("exclude", []), *hatch.get("exclude", [])]
    assert not [e for e in excluded if "demo" in e or e.endswith(".gz")], excluded
    assert DEMO.is_dir()


@pytest.mark.parametrize("key", ["bakta_json", *LAYER_KEYS])
def test_every_precomputed_layer_is_bundled(key: str) -> None:
    """Every ``PrecomputedConfig`` field (incl. ``rbs_novel``, ADR-0024) has a gzip layer."""
    path = _layer_file(key)
    assert path.is_file(), f"{path.relative_to(REPO)} missing"
    assert path.read_bytes()[:2] == b"\x1f\x8b", f"{path.name} is not gzip"


def test_demo_directory_at_most_1_mib() -> None:
    assert DEMO.is_dir(), "gippa/data/demo/ missing"
    total = sum(p.stat().st_size for p in DEMO.rglob("*") if p.is_file())
    assert 0 < total <= MIB, f"gippa/data/demo/ is {total} bytes (cap {MIB})"


def test_no_cache_or_stray_files_in_demo_directory() -> None:
    assert DEMO.is_dir()
    allowed = {"sim100k.fna.gz", "sim100k_peptides.tsv.gz", "PROVENANCE.md"}
    allowed |= {str(_layer_file(k).relative_to(DEMO)) for k in ["bakta_json", *LAYER_KEYS]}
    present = {str(p.relative_to(DEMO)) for p in DEMO.rglob("*") if p.is_file()}
    assert present - allowed == set()


def test_wheel_genome_is_the_sim100k_test_genome_byte_for_byte() -> None:
    assert (DEMO / "sim100k.fna.gz").read_bytes() == GENOME.read_bytes()


def test_wheel_peptides_are_the_tutorial_table_byte_for_byte() -> None:
    assert (DEMO / "sim100k_peptides.tsv.gz").read_bytes() == PEPTIDES.read_bytes()


@pytest.mark.parametrize("key", ["bakta_json", *LAYER_KEYS])
def test_bundled_gzip_has_no_mtime(key: str) -> None:
    """Reproducible bytes, as make_demo_peptides.py writes (``mtime=0``)."""
    head = _layer_file(key).read_bytes()[:10]
    assert head[:2] == b"\x1f\x8b"
    assert head[4:8] == b"\x00\x00\x00\x00", "gzip header carries an mtime"


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    import gippasim

    return gippasim.simulate(
        gippasim.SIM100K, gippasim.SIM100K_SEED, tmp_path_factory.mktemp("sim100k") / "sim"
    )


def test_sim100k_seed_is_the_banner_seed() -> None:
    import gippasim

    assert gippasim.SIM100K_SEED == SEED


@pytest.mark.parametrize("key", [k for k in LAYER_KEYS if k != "rbs_novel"])
def test_bundled_layer_equals_gippasim_output(sim: Any, key: str) -> None:
    assert _layer_file(key).is_file(), f"{key} layer missing"
    assert _plain(_layer_file(key)) == Path(sim.precomputed[key]).read_bytes()


def test_bundled_rbs_novel_equals_gippasim_output(sim: Any) -> None:
    """``rbs_novel`` is not in ``SimResult.precomputed``; gippasim writes it beside the others."""
    path = _layer_file("rbs_novel")
    assert path.is_file(), "rbs_novel layer missing (ADR-0024)"
    written = Path(sim.bakta_json).parent / "precomputed" / "rbs_novel.tsv"
    assert _plain(path) == written.read_bytes()


def test_bundled_bakta_json_equals_gippasim_output(sim: Any) -> None:
    assert _plain(_layer_file("bakta_json")) == Path(sim.bakta_json).read_bytes()


def test_bundled_layers_carry_no_absolute_paths() -> None:
    for key in ["bakta_json", *LAYER_KEYS]:
        assert _layer_file(key).is_file(), key
        text = _plain(_layer_file(key)).decode()
        assert "/tmp/" not in text and "/home/" not in text, key


def test_provenance_names_seed_generator_and_files() -> None:
    text = (DEMO / "PROVENANCE.md").read_text()
    low = text.lower()
    assert "synthetic" in low
    assert str(SEED) in text
    assert "make_demo_layers.py" in text
    assert "make_demo_peptides.py" in text
    for name in ("sim100k.fna.gz", "sim100k_peptides.tsv.gz", "rbs_novel"):
        assert name in text, name


# ------------------------------------------------------------------ make_demo_layers.py (AC 2)


def _script(*args: str, script: Path = SCRIPT) -> subprocess.CompletedProcess[str]:
    assert script.is_file(), f"{script} does not exist (M11-15 item 2)"
    return subprocess.run(
        [sys.executable, str(script), *args], capture_output=True, text=True, timeout=300
    )


def _tree(root: Path) -> dict[str, tuple[int, str]]:
    return {
        str(p.relative_to(root)): (p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).hexdigest())
        for p in sorted(root.rglob("*"))
        if p.is_file()
    }


def test_make_demo_layers_check_passes_and_writes_nothing() -> None:
    assert DEMO.is_dir()
    before = _tree(DEMO)
    res = _script("--check")
    assert res.returncode == 0, res.stdout + res.stderr
    assert _tree(DEMO) == before


def test_make_demo_layers_has_check_option_and_docstring() -> None:
    assert SCRIPT.is_file()
    text = SCRIPT.read_text()
    assert "--check" in text
    assert "SIM100K" in text and "rbs_novel" in text


@pytest.fixture
def fake_checkout(tmp_path: Path) -> Iterator[Path]:
    """A checkout-shaped tree with a copy of the script and of ``gippa/data/demo``; ``src`` holds
    only the demo data (the real ``gippa`` stays importable through PYTHONPATH), so the script's
    checkout-relative paths resolve inside ``tmp_path``."""
    assert SCRIPT.is_file() and DEMO.is_dir()
    root = tmp_path / "checkout"
    (root / "scripts" / "docs").mkdir(parents=True)
    shutil.copy2(SCRIPT, root / "scripts" / "docs" / SCRIPT.name)
    shutil.copytree(DEMO, root / "src" / "gippa" / "data" / "demo")
    (root / "tests" / "data" / "genomes").mkdir(parents=True)
    shutil.copy2(GENOME, root / "tests" / "data" / "genomes" / GENOME.name)
    (root / "tests" / "sim").symlink_to(REPO / "tests" / "sim")
    (root / "docs" / "tutorials" / "data").mkdir(parents=True)
    shutil.copy2(PEPTIDES, root / "docs" / "tutorials" / "data" / PEPTIDES.name)
    yield root


def test_make_demo_layers_check_fails_on_a_changed_layer(fake_checkout: Path) -> None:
    layer = fake_checkout / "src" / "gippa" / "data" / "demo" / "layers" / "rbs_novel.tsv.gz"
    rows = gzip.decompress(layer.read_bytes()).splitlines(keepends=True)
    layer.write_bytes(gzip.compress(b"".join(rows[:-1]), mtime=0))
    res = _script("--check", script=fake_checkout / "scripts" / "docs" / SCRIPT.name)
    assert res.returncode == 1, res.stdout + res.stderr
    assert "rbs_novel" in res.stdout + res.stderr


def test_make_demo_layers_rewrites_a_removed_layer(fake_checkout: Path) -> None:
    demo = fake_checkout / "src" / "gippa" / "data" / "demo"
    want = (demo / "layers" / "promoters.tsv.gz").read_bytes()
    (demo / "layers" / "promoters.tsv.gz").unlink()
    script = fake_checkout / "scripts" / "docs" / SCRIPT.name
    res = _script(script=script)
    assert res.returncode == 0, res.stdout + res.stderr
    assert (demo / "layers" / "promoters.tsv.gz").read_bytes() == want
    again = _script("--check", script=script)
    assert again.returncode == 0, again.stdout + again.stderr


# ------------------------------------------------------------------ check_dist rule (AC 1)


def _check_dist() -> Any:
    spec = importlib.util.spec_from_file_location("_check_dist_m11_15b", CHECK_DIST)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses look the module up by name
    spec.loader.exec_module(mod)
    return mod


def _demo_members(size: int = 1000) -> dict[str, bytes]:
    names = ["sim100k.fna.gz", "sim100k_peptides.tsv.gz", "PROVENANCE.md"]
    names += [str(_layer_file(k).relative_to(DEMO)) for k in ["bakta_json", *LAYER_KEYS]]
    return {f"gippa/data/demo/{n}": b"x" * size for n in names}


def _fake_wheel(dist: Path, cd: Any, demo: dict[str, bytes]) -> Path:
    version = cd.checkout_version()
    path = dist / f"gippa-{version}-py3-none-any.whl"
    di = f"gippa-{version}.dist-info"
    files = {p: b"x\n" for p in cd.wheel_required()}
    files[f"{di}/METADATA"] = (
        f"Metadata-Version: 2.4\nName: gippa\nVersion: {version}\nLicense-Expression: MIT\n"
    ).encode()
    files[f"{di}/entry_points.txt"] = b"[console_scripts]\ngippa = gippa.cli:main\n"
    files[f"{di}/licenses/LICENSE"] = b"MIT License\n"
    files.update(demo)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as zf:
        for arc, data in files.items():
            info = zipfile.ZipInfo(arc, date_time=(2026, 10, 27, 0, 0, 0))
            info.external_attr = (0o100000 | 0o644) << 16
            zf.writestr(info, data)
    return path


def test_check_dist_accepts_a_wheel_with_the_demo(tmp_path: Path) -> None:
    cd = _check_dist()
    whl = _fake_wheel(tmp_path, cd, _demo_members())
    assert cd.check_wheel(whl, cd.checkout_version()) == []


def test_check_dist_fails_a_wheel_without_the_demo(tmp_path: Path) -> None:
    cd = _check_dist()
    whl = _fake_wheel(tmp_path, cd, {})
    fails = cd.check_wheel(whl, cd.checkout_version())
    assert any("gippa/data/demo" in where for _rule, where in fails), fails


@pytest.mark.parametrize("name", ["sim100k.fna.gz", "sim100k_peptides.tsv.gz", "PROVENANCE.md"])
def test_check_dist_fails_a_wheel_missing_a_demo_file(tmp_path: Path, name: str) -> None:
    cd = _check_dist()
    demo = _demo_members()
    del demo[f"gippa/data/demo/{name}"]
    fails = cd.check_wheel(_fake_wheel(tmp_path, cd, demo), cd.checkout_version())
    assert any(name in where for _rule, where in fails), fails


def test_check_dist_fails_a_demo_directory_over_1_mib(tmp_path: Path) -> None:
    cd = _check_dist()
    demo = _demo_members(size=MIB // 8)  # 11 files x 128 KiB: each small, together > 1 MiB
    assert sum(map(len, demo.values())) > MIB
    fails = cd.check_wheel(_fake_wheel(tmp_path, cd, demo), cd.checkout_version())
    assert any("gippa/data/demo" in where for _rule, where in fails), fails


# ------------------------------------------------------------------ CLI


@pytest.fixture
def clean_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """No DB directory, no container runtime, no numba cache override; cwd is a fresh dir."""
    for name in ("GIPPA_DB_DIR", "GIPPA_CONTAINER_RUNTIME", "NUMBA_CACHE_DIR"):
        monkeypatch.delenv(name, raising=False)
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


class Recorder:
    """Replaces a ``gippa.cli`` function; records its bound arguments and the environment."""

    def __init__(self, original: Callable[..., Any], then: Callable[[dict[str, Any]], Any]):
        self.sig = inspect.signature(original)
        self.then = then
        self.calls: list[dict[str, Any]] = []
        self.env: list[dict[str, str]] = []

    def __call__(self, *args: Any, **kwargs: Any) -> Any:
        bound = self.sig.bind(*args, **kwargs)
        bound.apply_defaults()
        self.calls.append(dict(bound.arguments))
        self.env.append(dict(os.environ))
        return self.then(dict(bound.arguments))


def _fake_local_success(arguments: dict[str, Any]) -> None:
    report = Path(arguments["outdir"]) / "report" / "report.html"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(f"<html><body>{BANNER}</body></html>\n")
    sys.exit(0)


@pytest.fixture
def local(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    import gippa.cli as gc

    rec = Recorder(gc._run_local, _fake_local_success)
    monkeypatch.setattr(gc, "_run_local", rec)
    return rec


@pytest.fixture
def no_nextflow(monkeypatch: pytest.MonkeyPatch) -> None:
    import gippa.cli as gc

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("gippa demo (without --full) must not launch Nextflow")

    monkeypatch.setattr(gc, "_run_nextflow", boom)
    monkeypatch.setattr(gc, "_run_preflight", boom)


def _invoke(*args: str) -> Any:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["demo", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def test_demo_is_a_visible_command() -> None:
    from gippa.cli import cli

    assert "demo" in cli.commands
    assert not cli.commands["demo"].hidden


def test_demo_help_lists_options_and_default_outdir() -> None:
    res = _invoke("--help")
    assert res.exit_code == 0, res.output
    for opt in ("--outdir", "--full", "--open"):
        assert opt in res.output
    from gippa.cli import cli

    params = {p.name: p for p in cli.commands["demo"].params}
    assert str(params["outdir"].default) in ("gippa-demo", "./gippa-demo")
    assert params["full"].default is False
    assert params["open"].default is False


def test_example_2_1_prints_done_line_with_absolute_report_path(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    res = _invoke()
    assert res.exit_code == 0, res.output
    m = DONE.search(res.output)
    assert m is not None, res.output
    report = Path(m.group(2).strip())
    assert report.is_absolute()
    assert report == (clean_env / "gippa-demo" / "report" / "report.html").resolve()
    assert report.is_file()
    assert len(local.calls) == 1


def test_relative_outdir_is_printed_absolute(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    res = _invoke("--outdir", "sub/mydemo")
    assert res.exit_code == 0, res.output
    m = DONE.search(res.output)
    assert m is not None, res.output
    assert Path(m.group(2).strip()) == (clean_env / "sub/mydemo/report/report.html").resolve()


def test_demo_runs_the_local_executor_on_bundled_inputs(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    res = _invoke()
    assert res.exit_code == 0, res.output
    (call,) = local.calls
    cfg = call["cfg"]
    assert _plain(Path(call["genome"])) == _plain(DEMO / "sim100k.fna.gz")
    peptides = call["inputs"].get("peptides")
    assert peptides is not None
    assert _plain(Path(peptides)) == _plain(DEMO / "sim100k_peptides.tsv.gz")
    assert sorted(cfg.modules) == sorted(DEMO_MODULES)


def test_demo_sets_every_precomputed_layer_from_the_bundle(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    """Every ``PrecomputedConfig`` field the demo modules read, ``rbs_novel`` included."""
    res = _invoke()
    assert res.exit_code == 0, res.output
    (call,) = local.calls
    pre = call["cfg"].annotation.precomputed
    for key in LAYER_KEYS:
        path = getattr(pre, key)
        assert path is not None, f"annotation.precomputed.{key} not set"
        assert _plain(Path(path)) == _plain(_layer_file(key)), key
    bakta = pre.bakta_json or call["inputs"].get("annotation")
    assert bakta is not None
    assert _plain(Path(bakta)) == _plain(_layer_file("bakta_json"))


def test_demo_uses_no_database_directory(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    res = _invoke()
    assert res.exit_code == 0, res.output
    (call,) = local.calls
    assert call["inputs"].get("db_dir") is None
    assert call["info"].get("db_dir") is None


def test_demo_needs_no_nextflow_on_path(
    clean_env: Path, local: Recorder, no_nextflow: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = clean_env.parent / "emptybin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    res = _invoke()
    assert res.exit_code == 0, res.output
    assert DONE.search(res.output) is not None


def test_numba_cache_dir_points_inside_outdir_during_the_run(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    """numba's ``cache=True`` must not write into a read-only site-packages (pip/conda installs)."""
    res = _invoke("--outdir", "d")
    assert res.exit_code == 0, res.output
    (env,) = local.env
    value = env.get("NUMBA_CACHE_DIR")
    assert value, "NUMBA_CACHE_DIR not set for the demo run"
    assert Path(value).resolve().is_relative_to((clean_env / "d").resolve())
    assert "NUMBA_CACHE_DIR" not in os.environ, "NUMBA_CACHE_DIR leaked after the demo"


def test_numba_cache_dir_set_by_the_user_is_kept(
    clean_env: Path, local: Recorder, no_nextflow: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    mine = str(clean_env.parent / "my-numba-cache")
    monkeypatch.setenv("NUMBA_CACHE_DIR", mine)
    res = _invoke()
    assert res.exit_code == 0, res.output
    (env,) = local.env
    assert env.get("NUMBA_CACHE_DIR") == mine
    assert os.environ.get("NUMBA_CACHE_DIR") == mine


# --- example 2.2: an existing outdir is refused ---------------------------------------------


def test_example_2_2_existing_default_outdir_refused_untouched(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    out = clean_env / "gippa-demo"
    out.mkdir()
    (out / "keep.txt").write_text("mine\n")
    res = _invoke()
    assert res.exit_code == 1, res.output
    assert f"gippa-demo {REFUSAL}" in res.output
    assert local.calls == []
    assert sorted(p.name for p in out.iterdir()) == ["keep.txt"]
    assert (out / "keep.txt").read_text() == "mine\n"
    assert DONE.search(res.output) is None


def test_example_2_2_second_demo_run_refused(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    first = _invoke()
    assert first.exit_code == 0, first.output
    report = clean_env / "gippa-demo" / "report" / "report.html"
    before = report.read_bytes()
    second = _invoke()
    assert second.exit_code == 1, second.output
    assert REFUSAL in second.output
    assert len(local.calls) == 1
    assert report.read_bytes() == before


def test_existing_empty_outdir_is_refused_too(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    (clean_env / "gippa-demo").mkdir()
    res = _invoke()
    assert res.exit_code == 1, res.output
    assert REFUSAL in res.output
    assert local.calls == []


@settings(
    max_examples=20, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    name=st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_-", min_size=1, max_size=12),
    as_file=st.booleans(),
    payload=st.binary(max_size=64),
)
def test_property_any_existing_outdir_is_refused_and_left_byte_identical(
    clean_env: Path, local: Recorder, no_nextflow: None, name: str, as_file: bool, payload: bytes
) -> None:
    if name.startswith("-"):
        name = "d" + name
    with tempfile.TemporaryDirectory(dir=clean_env) as tmp:
        out = Path(tmp) / name
        if as_file:
            out.write_bytes(payload)
        else:
            out.mkdir()
            (out / "f.bin").write_bytes(payload)
        before = _tree(Path(tmp))
        n_calls = len(local.calls)
        res = _invoke("--outdir", str(out))
        assert res.exit_code == 1, res.output
        assert REFUSAL in res.output
        assert name in res.output
        assert len(local.calls) == n_calls
        assert _tree(Path(tmp)) == before


@pytest.mark.skipif(os.name != "posix" or os.geteuid() == 0, reason="needs a non-root POSIX user")
def test_read_only_cwd_error_names_outdir(
    clean_env: Path, local: Recorder, no_nextflow: None
) -> None:
    clean_env.chmod(0o555)
    try:
        res = _invoke()
    finally:
        clean_env.chmod(0o755)
    assert res.exit_code == 1, res.output
    assert "--outdir" in res.output
    assert "Traceback" not in res.output
    assert local.calls == []
    assert not (clean_env / "gippa-demo").exists()


# --- example 2.5: --open ---------------------------------------------------------------------


@pytest.fixture
def browser(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    import webbrowser

    opened: list[str] = []

    def fake_open(url: str, *_a: Any, **_k: Any) -> bool:
        opened.append(url)
        return False  # no browser available (as on a headless host)

    monkeypatch.setattr(webbrowser, "open", fake_open)
    return opened


def _is_report(url: str, report: Path) -> bool:
    return url in (str(report), report.as_uri())


def test_example_2_5_open_without_display_prints_path_and_exits_0(
    clean_env: Path, local: Recorder, no_nextflow: None, browser: list[str]
) -> None:
    res = _invoke("--open")
    assert res.exit_code == 0, res.output
    m = DONE.search(res.output)
    assert m is not None, res.output
    report = Path(m.group(2).strip())
    assert len(browser) == 1
    assert _is_report(browser[0], report), browser


def test_open_survives_a_browser_error(
    clean_env: Path, local: Recorder, no_nextflow: None, monkeypatch: pytest.MonkeyPatch
) -> None:
    import webbrowser

    def broken(*_a: Any, **_k: Any) -> bool:
        raise webbrowser.Error("could not locate runnable browser")

    monkeypatch.setattr(webbrowser, "open", broken)
    res = _invoke("--open")
    assert res.exit_code == 0, res.output
    assert DONE.search(res.output) is not None


def test_no_open_flag_never_starts_a_browser(
    clean_env: Path, local: Recorder, no_nextflow: None, browser: list[str]
) -> None:
    res = _invoke()
    assert res.exit_code == 0, res.output
    assert browser == []


def test_failed_demo_does_not_open_a_browser(
    clean_env: Path, no_nextflow: None, browser: list[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    import gippa.cli as gc

    rec = Recorder(gc._run_local, lambda _a: sys.exit(4))
    monkeypatch.setattr(gc, "_run_local", rec)
    res = _invoke("--open")
    assert res.exit_code == 4, res.output
    assert browser == []
    assert DONE.search(res.output) is None


# --- --full (examples 2.3, 2.4) ------------------------------------------------------------


MISSING_DBS = "environment.tools: missing: bakta_db, genomad_db; run: gippa setup"


def _fake_preflight_failure(_arguments: dict[str, Any]) -> None:
    import click

    click.echo(MISSING_DBS, err=True)
    sys.exit(2)


@pytest.fixture
def preflight(monkeypatch: pytest.MonkeyPatch) -> Recorder:
    import gippa.cli as gc

    rec = Recorder(gc._run_preflight, _fake_preflight_failure)
    monkeypatch.setattr(gc, "_run_preflight", rec)

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("Nextflow launched although preflight failed")

    monkeypatch.setattr(gc, "_run_nextflow", boom)

    def no_local(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("--full must go through gippa run's Nextflow path")

    monkeypatch.setattr(gc, "_run_local", no_local)
    return rec


def test_example_2_3_full_before_setup_exits_2_with_preflight_message(
    clean_env: Path, preflight: Recorder
) -> None:
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    assert "run: gippa setup" in res.output
    assert "missing: bakta_db" in res.output
    assert len(preflight.calls) == 1
    assert DONE.search(res.output) is None


def test_full_runs_the_bundled_inputs_with_real_tools(clean_env: Path, preflight: Recorder) -> None:
    """Same genome and peptides; no precomputed layers (the real tools run: an install check)."""
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    (call,) = preflight.calls
    inputs = call["inputs"]
    assert _plain(Path(inputs["genome"])) == _plain(DEMO / "sim100k.fna.gz")
    assert _plain(Path(inputs["peptides"])) == _plain(DEMO / "sim100k_peptides.tsv.gz")
    pre = call["cfg"].annotation.precomputed
    assert all(getattr(pre, k) is None for k in ["bakta_json", *LAYER_KEYS]), pre


def test_full_without_function_resources_leaves_function_out(
    clean_env: Path, preflight: Recorder, monkeypatch: pytest.MonkeyPatch
) -> None:
    db = clean_env.parent / "db"
    db.mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    (call,) = preflight.calls
    assert sorted(call["cfg"].modules) == sorted(DEMO_MODULES)


def test_full_with_function_resources_adds_function(
    clean_env: Path, preflight: Recorder, monkeypatch: pytest.MonkeyPatch
) -> None:
    from helpers.function_db import write_function_db

    db = clean_env.parent / "db"
    write_function_db(db, pd.DataFrame({"locus_tag": [], "db": [], "id": []}))
    monkeypatch.setenv("GIPPA_DB_DIR", str(db))
    res = _invoke("--full")
    assert res.exit_code == 2, res.output
    (call,) = preflight.calls
    assert sorted(call["cfg"].modules) == sorted([*DEMO_MODULES, "function"])


def test_example_2_4_full_uses_the_profile_gippa_run_would_use(
    clean_env: Path, preflight: Recorder, tmp_path: Path
) -> None:
    """``--full`` passes no profile of its own: whatever ``gippa run`` resolves without
    ``--profile`` (the detected profile once 15a lands) reaches preflight."""
    from gippa.cli import cli

    direct = CliRunner().invoke(
        cli,
        ["run", "--genome", str(GENOME), "--outdir", str(tmp_path / "direct")],
    )
    assert direct.exit_code == 2, direct.output
    assert len(preflight.calls) == 1, direct.output
    want = preflight.calls[0]["profile"]
    res = _invoke("--full", "--outdir", str(tmp_path / "viademo"))
    assert res.exit_code == 2, res.output
    assert len(preflight.calls) == 2, res.output
    assert preflight.calls[1]["profile"] == want


def test_full_refuses_an_existing_outdir(clean_env: Path, preflight: Recorder) -> None:
    (clean_env / "gippa-demo").mkdir()
    res = _invoke("--full")
    assert res.exit_code == 1, res.output
    assert REFUSAL in res.output
    assert preflight.calls == []


# ------------------------------------------------------------------ manifest schema


def test_manifest_schema_accepts_demo_true_and_keeps_it_optional() -> None:
    """``run_manifest.schema.json`` has ``additionalProperties: false``: the demo's
    ``demo: true`` must be a declared boolean property, not required (ordinary runs omit it)."""
    from helpers.manifest_data import build, schema, schema_errors

    assert "demo" in schema()["properties"]
    assert "demo" not in schema()["required"]
    m = build()
    m.pop("demo", None)
    assert schema_errors(m) == []
    m["demo"] = True
    assert schema_errors(m) == []


@pytest.mark.parametrize("bad", ["yes", 1, None])
def test_manifest_schema_rejects_a_non_boolean_demo(bad: object) -> None:
    from helpers.manifest_data import build, schema_errors

    m = build()
    m["demo"] = bad
    assert schema_errors(m), bad
