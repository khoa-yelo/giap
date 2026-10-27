"""``scripts/pilot/cascadia_offline.py``: the offline Unimod launcher (M10-08a, deliverable 4).

Spec: docs/specs/M10-08-a3-feasibility.md, deliverable 4 and acceptance criterion 6. Cascadia 0.0.7
downloads ``http://www.unimod.org/xml/unimod_tables.xml`` at the end of every run (pyteomics
``mass.unimod.Unimod()`` with no path); R6 forbids it, so the launcher

(a) refuses a ``UNIMOD_XML`` that is not an existing local file (exit 2, before importing
    ``cascadia``);
(b) sets ``pyteomics.mass.unimod._unimod_xml_download_url`` to its absolute path;
(c) hands the rest of ``argv`` to ``cascadia.cascadia.main()`` unchanged.

Review B1 of PR #301 (lead decision): Cascadia 0.0.7 also reaches Unimod at the **start** of every
run: ``PeptideTokenizer.from_massivekb()`` -> ``pyteomics.proforma.UnimodResolver.load_database()``
-> psims ``obo_cache.resolve("http://www.unimod.org/obo/unimod.obo")`` -> ``Unimod()`` ->
``lxml.etree.parse("http://www.unimod.org/xml/unimod_tables.xml")``. So the launcher also

(d) registers a psims resolver for the Unimod OBO URL built on the local table, before Cascadia is
    imported, so no http URL is ever parsed;
(e) installs a Python socket guard: any connection attempt under the launcher is refused.

Stub ``cascadia``, ``pyteomics``, ``psims`` and ``lxml`` packages
(``helpers.a3_kit.write_stub_packages``) record what the launcher did; nothing real is imported.
The production call is ``python -I`` inside the ``gippa-cascadia`` env, where ``-I`` ignores
``PYTHONPATH``: most tests use ``PYTHONPATH`` without ``-I`` and two tests run ``-I`` against a
throw-away venv (no pip) holding the stubs.
No statistical claim here, so no calibration test (spec §8.4).
"""

from __future__ import annotations

import ast
import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest
from helpers.a3_kit import (
    DEFAULT_UNIMOD_URL,
    LAUNCHER,
    STUB_UNIMOD_XML,
    make_stub_venv,
    require,
    write_stub_packages,
)

PARSE_LOG = "stub_parse.jsonl"
ARGS = ["sequence", "/data/mzml/WT_1.mzML", "/db/cascadia/cascadia_0.0.7.ckpt", "-o",
        "/root/pilot/A3_cascadia/WT_1/WT_1", "-b", "32", "-t", "0.8"]  # fmt: skip


@pytest.fixture
def stubs(tmp_path: Path) -> Path:
    return write_stub_packages(tmp_path / "site")


@pytest.fixture
def unimod(tmp_path: Path) -> Path:
    p = tmp_path / "db" / "unimod_tables.xml"
    p.parent.mkdir(parents=True)
    p.write_bytes(STUB_UNIMOD_XML)
    return p


def _launch(
    tmp: Path,
    site: Path | None,
    *argv: str,
    python: str = sys.executable,
    isolated: bool = False,
    rc: int = 0,
    cwd: Path | None = None,
    extra_env: dict[str, str] | None = None,
) -> tuple[subprocess.CompletedProcess[str], dict[str, object] | None, bool]:
    require(LAUNCHER, "deliverable 4")
    log, imported = tmp / "stub_main.json", tmp / "stub_import.log"
    for p in (log, imported, tmp / PARSE_LOG):
        p.unlink(missing_ok=True)
    env = {"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "STUB_LOG": str(log),
           "STUB_IMPORT_LOG": str(imported), "STUB_RC": str(rc),
           "STUB_PARSE_LOG": str(tmp / PARSE_LOG), **(extra_env or {})}  # fmt: skip
    if site is not None:
        env["PYTHONPATH"] = str(site)
    cmd = [python, *(["-I"] if isolated else []), str(LAUNCHER), *argv]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60,
                         cwd=str(cwd or tmp), check=False)  # fmt: skip
    rec = json.loads(log.read_text()) if log.is_file() else None
    return res, rec, imported.is_file()


def _same_path(a: object, b: Path) -> bool:
    return isinstance(a, str) and os.path.isabs(a) and os.path.realpath(a) == os.path.realpath(b)


def test_sets_the_unimod_url_before_cascadia_runs(
    tmp_path: Path, stubs: Path, unimod: Path
) -> None:
    res, rec, imported = _launch(tmp_path, stubs, str(unimod), *ARGS)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None, "cascadia.cascadia.main() was not called"
    assert imported
    assert rec["url_in_main"] != DEFAULT_UNIMOD_URL
    assert _same_path(rec["url_in_main"], unimod), rec
    # set before cascadia is even imported (deliverable 4: (b) before (c))
    assert _same_path(rec["url_at_import"], unimod), rec


def test_argv_is_passed_through_unchanged(tmp_path: Path, stubs: Path, unimod: Path) -> None:
    odd = [*ARGS, "--path with spaces", "", "-w", "2"]
    res, rec, _ = _launch(tmp_path, stubs, str(unimod), *odd)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and rec["argv"] == odd


def test_relative_unimod_path_becomes_absolute(tmp_path: Path, stubs: Path, unimod: Path) -> None:
    rel = os.path.relpath(unimod, tmp_path)
    res, rec, _ = _launch(tmp_path, stubs, rel, *ARGS, cwd=tmp_path)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and _same_path(rec["url_in_main"], unimod), rec


def test_cascadia_exit_code_is_propagated(tmp_path: Path, stubs: Path, unimod: Path) -> None:
    res, rec, _ = _launch(tmp_path, stubs, str(unimod), *ARGS, rc=5)
    assert rec is not None
    assert res.returncode == 5


@pytest.mark.parametrize(
    "bad",
    [
        "http://www.unimod.org/xml/unimod_tables.xml",
        "https://www.unimod.org/xml/unimod_tables.xml",
        "file:///etc/unimod_tables.xml",
        "missing",
        "dir",
    ],
)
def test_refuses_a_non_local_or_missing_table_without_importing_cascadia(
    tmp_path: Path, stubs: Path, bad: str
) -> None:
    if bad == "missing":
        bad = str(tmp_path / "absent" / "unimod_tables.xml")
    elif bad == "dir":
        (tmp_path / "a_dir").mkdir()
        bad = str(tmp_path / "a_dir")
    res, rec, imported = _launch(tmp_path, stubs, bad, *ARGS)
    assert res.returncode == 2, res.stdout + res.stderr
    assert rec is None, "cascadia main() ran"
    assert not imported, "cascadia was imported before the table was checked"


def test_no_arguments_exits_2(tmp_path: Path, stubs: Path) -> None:
    res, rec, imported = _launch(tmp_path, stubs)
    assert res.returncode == 2, res.stdout + res.stderr
    assert rec is None and not imported


def test_runs_isolated_inside_an_env_that_has_cascadia(tmp_path: Path, unimod: Path) -> None:
    # The production call: `micromamba run -n gippa-cascadia python -I cascadia_offline.py …`.
    # No gippa on the path, no PYTHONPATH: only the env's own site-packages.
    env = make_stub_venv(tmp_path / "gippa-cascadia")
    res, rec, _ = _launch(tmp_path, None, str(unimod), *ARGS, python=str(env / "bin" / "python"),
                          isolated=True)  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and rec["argv"] == ARGS
    assert _same_path(rec["url_in_main"], unimod)


def test_launcher_imports_no_gippa_and_no_network_modules() -> None:
    tree = ast.parse(require(LAUNCHER, "deliverable 4").read_text())
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            mods.add(node.module.split(".")[0])
    assert "gippa" not in mods  # it runs inside gippa-cascadia, where gippa is not installed
    # `socket` is allowed: the launcher imports it to install the connection guard (B1).
    assert not mods & {"http", "requests", "ftplib", "urllib3"}


# ------------------------------------------------------------- B1: the start-of-run Unimod path


def _parsed(tmp: Path) -> list[str]:
    log = tmp / PARSE_LOG
    if not log.is_file():
        return []
    return [json.loads(ln)["source"] for ln in log.read_text().splitlines() if ln]


@pytest.mark.parametrize("isolated", [False, True], ids=["pythonpath", "venv-I"])
def test_psims_unimod_resolver_uses_the_local_table(
    tmp_path: Path, unimod: Path, isolated: bool
) -> None:
    # The stub cascadia main() runs proforma.UnimodResolver().load_database(), as every Cascadia
    # `sequence` run does when it builds its tokenizer. Without a local resolver, psims' default
    # parses http://www.unimod.org/xml/unimod_tables.xml (a plain-HTTP request on a permissive
    # lxml build).
    if isolated:
        env = make_stub_venv(tmp_path / "gippa-cascadia")
        res, rec, _ = _launch(tmp_path, None, str(unimod), *ARGS,
                              python=str(env / "bin" / "python"), isolated=True)  # fmt: skip
    else:
        res, rec, _ = _launch(tmp_path, write_stub_packages(tmp_path / "site"), str(unimod), *ARGS)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None, "cascadia.cascadia.main() was not called"
    parsed = _parsed(tmp_path)
    remote = [s for s in parsed if "://" in s]
    assert not remote, f"a Unimod URL was parsed under the launcher: {remote}"
    assert any(_same_path(s, unimod) for s in parsed), f"the local table was not parsed: {parsed}"
    assert _same_path(rec["unimod_db_source"], unimod), rec


def test_psims_resolver_is_registered_before_cascadia_is_imported(
    tmp_path: Path, stubs: Path, unimod: Path
) -> None:
    res, rec, _ = _launch(tmp_path, stubs, str(unimod), *ARGS)
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None
    assert rec["resolver_default_at_import"] is False, (
        "psims' default Unimod resolver was still registered when cascadia was imported"
    )


def test_refused_table_registers_nothing_and_parses_nothing(tmp_path: Path, stubs: Path) -> None:
    res, rec, imported = _launch(tmp_path, stubs, DEFAULT_UNIMOD_URL, *ARGS)
    assert res.returncode == 2, res.stdout + res.stderr
    assert rec is None and not imported
    assert _parsed(tmp_path) == []


# ------------------------------------------------------------------ B1: the socket guard


@pytest.fixture
def listener():  # type: ignore[no-untyped-def]
    import socket

    srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    srv.bind(("127.0.0.1", 0))
    srv.listen(8)
    try:
        yield srv
    finally:
        srv.close()


def test_python_sockets_are_refused_under_the_launcher(  # type: ignore[no-untyped-def]
    tmp_path: Path, stubs: Path, unimod: Path, listener
) -> None:
    import select

    port = listener.getsockname()[1]
    res, rec, _ = _launch(tmp_path, stubs, str(unimod), *ARGS,
                          extra_env={"STUB_NET_PORT": str(port)})  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    assert rec is not None and "net" in rec, rec
    net = rec["net"]
    assert set(net) == {"create_connection", "socket.connect", "socket.connect_ex"}, net
    leaked = {k: v for k, v in net.items() if v == "connected"}
    assert not leaked, f"a Python-level connection went through under the launcher: {net}"
    readable, _, _ = select.select([listener], [], [], 0.2)
    assert not readable, "the local listener received a connection"


def test_network_attempts_succeed_without_the_launcher(  # type: ignore[no-untyped-def]
    tmp_path: Path, stubs: Path, listener
) -> None:
    # Control: the stub really connects when the launcher is not in the way, so the refusal
    # above comes from the launcher's guard, not from the sandbox.
    log = tmp_path / "direct.json"
    code = ("import sys; from cascadia import cascadia; sys.argv = ['cascadia', 'sequence']; "
            "cascadia.main()")  # fmt: skip
    env = {"PATH": "/usr/bin:/bin", "LANG": "C.UTF-8", "PYTHONPATH": str(stubs),
           "STUB_LOG": str(log), "STUB_NET_PORT": str(listener.getsockname()[1])}  # fmt: skip
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env,
                         timeout=60, cwd=str(tmp_path), check=False)  # fmt: skip
    assert res.returncode == 0, res.stdout + res.stderr
    net = json.loads(log.read_text())["net"]
    assert net["create_connection"] == "connected", net
    assert net["socket.connect"] == "connected", net


@pytest.mark.skipif(not hasattr(__import__("socket"), "AF_UNIX"), reason="no AF_UNIX")
def test_unix_sockets_stay_open_ip_sockets_stay_refused(tmp_path: Path) -> None:
    """torch's DataLoader workers pass tensors over AF_UNIX sockets; a guard that refuses them
    kills every multi-worker run (job 1889042, lead fix 2026-10-24). IP stays refused."""
    code = textwrap.dedent(
        f"""
        import socket, sys
        sys.path.insert(0, {str(LAUNCHER.parent)!r})
        import cascadia_offline as L
        srv = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        srv.bind({str(tmp_path / "s.sock")!r}); srv.listen(1)
        L._guard_sockets()
        c = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        c.connect({str(tmp_path / "s.sock")!r})
        assert c.connect_ex({str(tmp_path / "s.sock")!r}) in (0, 106, 56), "connect_ex"
        for fam, addr in ((socket.AF_INET, ("127.0.0.1", 9)),):
            for how in ("connect", "connect_ex"):
                try:
                    getattr(socket.socket(fam, socket.SOCK_STREAM), how)(addr)
                except OSError as e:
                    assert "network access is disabled" in str(e), e
                else:
                    raise SystemExit(f"{{how}} on {{fam}} went through")
        try:
            socket.create_connection(("127.0.0.1", 9))
        except OSError as e:
            assert "network access is disabled" in str(e), e
        else:
            raise SystemExit("create_connection went through")
        print("ok")
        """
    )
    res = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True,
                         timeout=60, check=False)  # fmt: skip
    assert res.returncode == 0 and res.stdout.strip() == "ok", res.stdout + res.stderr
