"""ANTISMASH uses the image's databases on a fresh setup (#515).

``workflow/modules/local/antismash/main.nf`` always passed ``--databases ${db_dir}/antismash``, but
``gippa setup`` under a container profile records ``antismash: ok (databases in the antiSMASH
image)`` (``setupdb.run._antismash_image``, url ``image``) and never creates that directory, so the
task died: "No matching database in location .../db/antismash/nrps_pks/transATor".

Interface for the implementer:

``script:`` of ``workflow/modules/local/antismash/main.nf``
    Passes ``--databases ${db_dir}/antismash`` only when that path is a directory with content
    (``items._check_antismash``'s rule: ``setup`` adopts exactly such a directory, so existing
    installs such as FarmShare's ``local/db/antismash`` keep the flag), otherwise omits the flag.
    **The test runs in the task's shell, not in Groovy:** the task must look at the path it will
    hand to the tool, as seen inside the container (``-B db_dir:db_dir:ro``), not at what the head
    job sees. So Groovy ``def``s stay as they are (``gz``, ``input``), and the only Groovy
    interpolations are ``db_dir``, ``task.cpus``, ``task.process``, ``gbff``, ``input`` and the
    existing ``gz`` ternary; shell variables are written ``\\$x``. Nextflow runs the script under
    ``bash -ue``, so an unset variable or a failed ``ls`` must not abort a task. ``stub:`` is
    unchanged.

``gippa.annotate.antismash_cmd.antismash_command(gbff, outdir, threads, databases: Path | None)``
    The argument stays required and positional; ``None`` omits ``--databases``. The module's two
    command lines are exactly this builder's two forms, in this order (pinned below).

``meta.yml``, the builder's docstring and M2-18 lead decisions 2 and 4 ("always passes") describe
the old rule and need updating. GIPPA_SPEC.md §7.1.1's command has no ``--databases`` at all.

Real Nextflow run of the same cases: tests/integration/test_fresh_run_nextflow.py.
"""

from __future__ import annotations

import gzip
import json
import os
import re
import stat
import subprocess
import textwrap
from pathlib import Path

import pytest
import yaml

from gippa.annotate.antismash_cmd import antismash_command

MODULE = Path(__file__).resolve().parents[2] / "workflow" / "modules" / "local" / "antismash"
PROCESS = "GIPPA:ANNOTATE:ANTISMASH"
CPUS = 4
# antiSMASH's own failure, as in the issue, when --databases names a directory without them.
FAKE_ANTISMASH = """#!/bin/bash
if [ "${1:-}" = --version ]; then echo 'antiSMASH 8.0.4'; exit 0; fi
printf '%s\\n' "$@" > "$FAKE_ARGV"
prev=
for a in "$@"; do
  if [ "$prev" = --databases ] && [ ! -d "$a/nrps_pks" ]; then
    echo "ERROR antismash.detection.nrps_pks_domains: No matching database in location $a" >&2
    exit 1
  fi
  prev="$a"
done
mkdir -p "$PWD/antismash"
"""


def _blocks() -> tuple[str, str]:
    """(``script:`` section, ``stub:`` section) of the module, Groovy source."""
    text = (MODULE / "main.nf").read_text(encoding="utf-8")
    m = re.search(r"^\s*script:\s*$(.*?)^\s*stub:\s*$(.*)", text, re.M | re.S)
    assert m, "module has no script:/stub: sections"
    return m[1], m[2]


def _interpolate(src: str, values: dict[str, str]) -> str:
    """Resolve ``${expr}`` by exact expression (nested ``${}`` in a ternary included) and the
    GString escapes; a shell ``\\$x`` stays ``$x``. An unknown expression fails: see the module
    docstring for why the decision belongs to the shell."""
    out: list[str] = []
    i = 0
    while i < len(src):
        ch = src[i]
        if ch == "\\" and i + 1 < len(src):
            nxt = src[i + 1]
            if nxt != "\n":  # backslash-newline is a Groovy line continuation: both go
                out.append({"n": "\n", "t": "\t"}.get(nxt, nxt))
            i += 2
        elif src.startswith("${", i):
            depth, j = 1, i + 2
            while depth:
                depth += {"{": 1, "}": -1}.get(src[j], 0)
                j += 1
            expr = src[i + 2 : j - 1]
            assert expr in values, f"Groovy expression ${{{expr}}} in the module: {sorted(values)}"
            out.append(values[expr])
            i = j
        else:
            out.append(ch)
            i += 1
    return "".join(out)


def _triple(section: str) -> str:
    m = re.search(r'"""(.*?)"""', section, re.S)
    assert m, 'no """ block'
    return m[1]


def _render_script(db_dir: Path, gbff: str) -> str:
    """The ``script:`` block as Nextflow hands it to bash, for ``${db_dir}`` and a staged GBFF."""
    section = _blocks()[0]
    head = section.split('"""', 1)[0]
    defs = re.findall(r"^\s*def\s+(\w+)", head, re.M)
    assert defs == ["gz", "input"], f"Groovy defs in script: (decide in the shell): {defs}"
    gz = gbff.lower().endswith(".gz")
    meta = "genome"
    plain = f"{meta}.gbff" if gz else gbff
    values = {
        "db_dir": str(db_dir),
        "task.cpus": str(CPUS),
        "task.process": PROCESS,
        "gbff": gbff,
        "input": plain,
        "gz ? \"gzip -dc ${gbff} > ${input}\" : ''": f"gzip -dc {gbff} > {plain}" if gz else "",
    }
    return textwrap.dedent(_interpolate(_triple(section), values)).strip() + "\n"


def _bash_ue(script: str, work: Path, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    """Nextflow's default ``process.shell``: ``/bin/bash -ue``."""
    return subprocess.run(["/bin/bash", "-ue", "-c", script], cwd=work, env=env, text=True,
                          capture_output=True, check=False, timeout=60)  # fmt: skip


def _run_task(
    tmp_path: Path, db_dir: Path, gbff: str = "genome.gbff"
) -> tuple[subprocess.CompletedProcess[str], list[str], Path]:
    """Run the rendered script with a fake ``antismash``: (result, recorded argv, task dir)."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "antismash"
    exe.write_text(FAKE_ANTISMASH, encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    work = tmp_path / "task"
    work.mkdir()
    if gbff.endswith(".gz"):
        (work / gbff).write_bytes(gzip.compress(b"LOCUS x\n"))
    else:
        (work / gbff).write_text("LOCUS x\n", encoding="utf-8")
    argv_file = tmp_path / "argv.txt"
    env = {"PATH": f"{bindir}{os.pathsep}/usr/bin:/bin", "LC_ALL": "C", "FAKE_ARGV": str(argv_file)}
    res = _bash_ue(_render_script(db_dir, gbff), work, env)
    argv = argv_file.read_text(encoding="utf-8").splitlines() if argv_file.is_file() else []
    return res, argv, work


def _with_databases(db_dir: Path, gbff: str = "genome.gbff") -> list[str]:
    return antismash_command(Path(gbff), Path("antismash"), CPUS, db_dir / "antismash")[1:]


def _image_databases(gbff: str = "genome.gbff") -> list[str]:
    return antismash_command(Path(gbff), Path("antismash"), CPUS, None)[1:]


def _existing_install(db_dir: Path) -> Path:
    """The layout of FarmShare's ``local/db/antismash`` (trimmed): a marker and the data dirs."""
    root = db_dir / "antismash"
    for sub in ("as-js", "clusterblast", "nrps_pks/transATor", "pfam"):
        (root / sub).mkdir(parents=True)
    (root / ".giap-complete").touch()
    (root / "pfam" / "Pfam-A.hmm").write_text("HMMER3/f\n", encoding="utf-8")
    return root


# ---------------------------------------------------------------- the builder


def test_builder_none_omits_the_flag() -> None:
    got = antismash_command(Path("g.gbff"), Path("o"), 8, None)
    assert got == ["antismash", "--genefinding-tool", "none", "--minimal", "--cpus", "8",
                   "--output-dir", "o", "g.gbff"]  # fmt: skip


def test_builder_with_a_path_keeps_the_flag() -> None:
    got = antismash_command(Path("g.gbff"), Path("o"), 8, Path("d/antismash"))
    assert got[got.index("--databases") + 1] == "d/antismash"


# ---------------------------------------------------------------- fresh setup: the image's


def test_fresh_setup_runs_without_the_flag(tmp_path: Path) -> None:
    """Example of #515: the DB directory exists, ``antismash/`` does not."""
    db = tmp_path / "db"
    db.mkdir()
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases()


def test_empty_directory_runs_without_the_flag(tmp_path: Path) -> None:
    db = tmp_path / "db"
    (db / "antismash").mkdir(parents=True)
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases()


def test_a_file_at_the_path_runs_without_the_flag(tmp_path: Path) -> None:
    db = tmp_path / "db"
    db.mkdir()
    (db / "antismash").write_text("not a directory\n", encoding="utf-8")
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases()


def test_a_dangling_symlink_runs_without_the_flag(tmp_path: Path) -> None:
    db = tmp_path / "db"
    db.mkdir()
    (db / "antismash").symlink_to(tmp_path / "gone")
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases()


def test_no_db_directory_at_all_runs_without_the_flag(tmp_path: Path) -> None:
    """``annotate.nf`` falls back to ``${launchDir}/db``, which need not exist."""
    res, argv, _ = _run_task(tmp_path, tmp_path / "no_such_db")
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases()


def test_gzipped_gbff_without_databases(tmp_path: Path) -> None:
    """The decompress line and the conditional live in one script."""
    db = tmp_path / "db"
    db.mkdir()
    res, argv, work = _run_task(tmp_path, db, "bakta.gbff.gz")
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases("genome.gbff")
    assert (work / "genome.gbff").read_text(encoding="utf-8") == "LOCUS x\n"


def test_setup_from_the_image_then_the_task(tmp_path: Path) -> None:
    """The chain of the issue: ``install`` records the image and creates nothing; the task then
    runs on the image's databases."""
    from gippa.setupdb.items import SOURCES
    from gippa.setupdb.run import install

    class Runner:
        def run(self, _tool: str, _argv: list[str]) -> int:
            return 0

    class NoDownload:
        def fetch(self, _url: str, _dest: Path) -> None:
            raise AssertionError("the image's databases need no download")

    db = tmp_path / "db"
    db.mkdir()
    rec = install(SOURCES["antismash"], db, Runner(), NoDownload())
    assert (rec.status, rec.url) == ("ok", "image")
    assert not (db / "antismash").exists()
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _image_databases()


# ---------------------------------------------------------------- existing installs keep the flag


def test_existing_install_keeps_the_flag(tmp_path: Path) -> None:
    db = tmp_path / "db"
    _existing_install(db)
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _with_databases(db)


def test_symlinked_install_keeps_the_flag(tmp_path: Path) -> None:
    """A DB directory that links to a shared copy."""
    shared = tmp_path / "shared"
    _existing_install(shared)
    db = tmp_path / "db"
    db.mkdir()
    (db / "antismash").symlink_to(shared / "antismash")
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert argv == _with_databases(db)


@pytest.mark.parametrize("state", ["absent", "empty", "file", "installed"])
def test_flag_agrees_with_what_setup_adopts(tmp_path: Path, state: str) -> None:
    """The flag is passed exactly when ``gippa setup`` verifies the directory as installed."""
    from gippa.setupdb.items import SOURCES, verify

    db = tmp_path / "db"
    db.mkdir()
    if state == "empty":
        (db / "antismash").mkdir()
    elif state == "file":
        (db / "antismash").write_text("x\n", encoding="utf-8")
    elif state == "installed":
        _existing_install(db)
    adopted = verify(SOURCES["antismash"], db).status == "ok"
    assert adopted is (state == "installed")
    res, argv, _ = _run_task(tmp_path, db)
    assert res.returncode == 0, res.stderr
    assert ("--databases" in argv) is adopted, argv


# ---------------------------------------------------------------- the stub path


def test_stub_writes_its_outputs_without_the_tool_or_the_databases(tmp_path: Path) -> None:
    """``-stub-run`` needs neither: no ``antismash`` on ``PATH``, no database directory."""
    src = _triple(_blocks()[1])
    assert "db_dir" not in src and "--databases" not in src, src
    script = _interpolate(src, {"meta_id": "genome", "task.process": PROCESS})
    work = tmp_path / "task"
    work.mkdir()
    res = _bash_ue(textwrap.dedent(script).strip() + "\n", work, {"PATH": "/usr/bin:/bin"})
    assert res.returncode == 0, res.stderr
    doc = json.loads((work / "antismash" / "genome.json").read_text(encoding="utf-8"))
    assert doc == {"records": []}
    versions = yaml.safe_load((work / "versions.yml").read_text(encoding="utf-8"))
    assert versions == {PROCESS: {"antismash": "8.0.4"}}
