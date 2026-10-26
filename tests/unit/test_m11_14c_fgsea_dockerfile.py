"""M11-14c: the fgsea image with ``giap_fgsea.R`` baked in (static), the NOTICE line, and R5.

Spec docs/specs/M11-14-container-install.md item 1 ("fgsea (14c)"), item 2 (key inputs; "fgsea
adds ``workflow/bin/giap_fgsea.R``"; 14b's locked, checked builds), AC 7 (the script is a key
input of ``fgsea``); ADR-0060 §1; R5 (R is GPL-2 and fgsea MIT: they run only as a separate
process, never imported). Static: the Dockerfile, locks, NOTICE.md and ``src/`` are only read.

Readings chosen where the spec leaves room (consistent with 14b's
``tests/unit/test_m11_14b_dockerfile_locks.py``):

* The Dockerfile installs ``workflow/envs/locks/fgsea.explicit.txt`` with ``micromamba
  install|create -f <its COPY destination>`` and never solves ``fgsea.yml``; it fetches no URL
  (no ``ADD http…``; any URL ``ADD`` would need ``--checksum=sha256:<64 hex>`` anyway) and installs
  no R package at build time (no ``install.packages``/``BiocManager``).
* ``workflow/bin/giap_fgsea.R`` is copied to exactly ``/opt/<pkg>/giap_fgsea.R`` (``<pkg>`` =
  ``[project] name`` of pyproject.toml), the path the image argv names.
* Every ``FROM`` is pinned by digest. ``ENV PATH`` puts ``/opt/conda/bin`` (the micromamba root
  prefix) first, because ``apptainer exec`` does not run the image entrypoint that would activate
  the base env, and the argv calls ``Rscript`` by name.
* Labels: ``org.opencontainers.image.title`` is ``<pkg>-fgsea`` and
  ``org.opencontainers.image.licenses`` names both ``MIT`` (fgsea) and ``GPL-2`` (R).
* NOTICE.md has a paragraph naming the ``<pkg>-fgsea`` image with ``GPL-2`` and ``MIT``.
"""

from __future__ import annotations

import ast
import re
import sys
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
DOCKERFILE = REPO / "workflow" / "containers" / "fgsea" / "Dockerfile"
LOCK = "workflow/envs/locks/fgsea.explicit.txt"
SCRIPT = "workflow/bin/giap_fgsea.R"
BAKED = f"/opt/{PKG}/giap_fgsea.R"
FGSEA_R = REPO / "src" / "giap" / "io" / "fgsea_r.py"


# ---------------------------------------------------------------- helpers (as in 14b)


def _text() -> str:
    assert DOCKERFILE.is_file(), "workflow/containers/fgsea/Dockerfile is missing (M11-14c)"
    return DOCKERFILE.read_text(encoding="utf-8")


def _instructions(text: str) -> list[tuple[str, str]]:
    """``(INSTRUCTION, rest)`` per instruction; comments dropped, continuations joined."""
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    joined = re.sub(r"\\[ \t]*\n", " ", "\n".join(lines))
    out = []
    for line in joined.splitlines():
        words = line.strip().split(None, 1)
        if len(words) == 2:
            out.append((words[0].upper(), words[1].strip()))
    return out


def _split_flags(rest: str) -> tuple[list[str], list[str]]:
    words = rest.split()
    return [w for w in words if w.startswith("--")], [w for w in words if not w.startswith("--")]


def _copies() -> dict[str, str]:
    """COPY source -> destination (build-context copies only)."""
    out = {}
    for op, rest in _instructions(_text()):
        if op == "COPY":
            flags, paths = _split_flags(rest)
            if any(f.startswith("--from") for f in flags):
                continue
            for src in paths[:-1]:
                out[src.removeprefix("./")] = paths[-1]
    return out


def _runs() -> list[str]:
    return [rest for op, rest in _instructions(_text()) if op == "RUN"]


def _labels() -> str:
    return " ".join(rest for op, rest in _instructions(_text()) if op == "LABEL")


# ---------------------------------------------------------------- the Dockerfile


def test_fgsea_dockerfile_exists() -> None:
    assert DOCKERFILE.is_file(), "workflow/containers/fgsea/Dockerfile is missing (M11-14c)"


def test_every_base_image_is_pinned_by_digest() -> None:
    froms = [rest for op, rest in _instructions(_text()) if op == "FROM"]
    assert froms, "no FROM"
    for rest in froms:
        image = next(w for w in rest.split() if not w.startswith("--"))
        assert re.search(r"@sha256:[0-9a-f]{64}$", image), f"base image not pinned: {image}"


def test_installs_the_explicit_lock() -> None:
    copies = _copies()
    assert LOCK in copies, f"COPY {LOCK} missing"
    mamba = [r for r in _runs() if re.search(r"micromamba\s+(install|create)", r)]
    assert mamba, "no micromamba install/create"
    for r in mamba:
        assert re.search(rf"-f\s+{re.escape(copies[LOCK])}(\s|$)", r), r


def test_never_solves_the_env_file() -> None:
    for r in _runs():
        if re.search(r"micromamba\s+(install|create)", r):
            assert ".yml" not in r, r
    assert "workflow/envs/fgsea.yml" not in _copies() or all(
        _copies()["workflow/envs/fgsea.yml"] not in r for r in _runs()
    )


def test_bakes_the_script_at_the_argv_path() -> None:
    copies = _copies()
    assert SCRIPT in copies, f"COPY {SCRIPT} missing"
    dest = copies[SCRIPT]
    assert dest in {BAKED, f"/opt/{PKG}/"}, dest


def test_fetches_no_url_and_installs_no_r_package_at_build_time() -> None:
    for op, rest in _instructions(_text()):
        if op == "ADD":
            _, paths = _split_flags(rest)
            urls = [p for p in paths[:-1] if re.match(r"https?://", p)]
            assert not urls, f"ADD of a URL in the fgsea image: {urls}"
    for r in _runs():
        assert "install.packages" not in r, r
        assert "BiocManager" not in r, r
        assert not re.search(r"\b(curl|wget)\b", r), r


def test_rscript_is_on_path_without_the_entrypoint() -> None:
    """``apptainer exec`` skips the micromamba entrypoint; ``Rscript`` must still resolve."""
    envs = [rest for op, rest in _instructions(_text()) if op == "ENV"]
    path = [e for e in envs if re.match(r"PATH[= ]", e)]
    assert path, "no ENV PATH"
    assert any("/opt/conda/bin" in e for e in path), path


def test_labels_name_the_image_and_both_licences() -> None:
    labels = _labels()
    assert re.search(rf'org\.opencontainers\.image\.title="?{re.escape(PKG)}-fgsea"?', labels)
    m = re.search(r'org\.opencontainers\.image\.licenses="([^"]+)"', labels)
    assert m, labels
    assert "MIT" in m.group(1) and "GPL-2" in m.group(1), m.group(1)


def test_no_pip_lock_is_needed() -> None:
    """fgsea.yml has no pip section; the image installs nothing with pip."""
    assert not (REPO / "workflow/envs/locks/fgsea.pip.txt").exists()
    assert not [r for r in _runs() if "pip install" in r]


def test_the_lock_carries_the_registry_pins() -> None:
    from giap.annotate.tools import TOOLS

    text = (REPO / LOCK).read_text(encoding="utf-8")
    assert f"/bioconductor-fgsea-{TOOLS['fgsea'].version}-" in text
    assert "/r-base-4.5.3-" in text
    assert "/r-data.table-1.18.6.1-" in text


# ---------------------------------------------------------------- the key (item 2, AC 7)


def test_dockerfile_and_script_are_key_inputs() -> None:
    from giap.io.image_key import key_inputs

    got = key_inputs("fgsea", REPO)
    for rel in ("workflow/containers/fgsea/Dockerfile", SCRIPT, LOCK, "workflow/envs/fgsea.yml"):
        assert rel in got, rel


def test_every_copy_source_is_a_key_input() -> None:
    from giap.io.image_key import uncovered_copy_sources

    _text()
    assert uncovered_copy_sources("fgsea", REPO) == []


def test_script_is_executable_and_keyed_as_such() -> None:
    """``COPY`` keeps the exec bit, so the key records 100755 for the script."""
    from giap.io.image_key import key_text

    _text()
    line = next(ln for ln in key_text("fgsea", REPO).splitlines() if ln.startswith(SCRIPT + "\t"))
    assert line.split("\t")[1] == "100755"


# ---------------------------------------------------------------- NOTICE (item 1)


def test_notice_names_the_fgsea_image_with_r_and_fgsea_licences() -> None:
    notice = (REPO / "NOTICE.md").read_text(encoding="utf-8")
    paras = [p for p in re.split(r"\n\s*\n", notice) if f"{PKG}-fgsea" in p]
    assert paras, f"NOTICE.md does not name the {PKG}-fgsea image"
    assert any("GPL-2" in p and "MIT" in p for p in paras), paras


# ---------------------------------------------------------------- R5: only a separate process


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.add(node.module.split(".")[0])
    return names


def test_fgsea_r_imports_only_the_stdlib_and_giap() -> None:
    """The image path shells out to ``docker``/``apptainer``: no SDK, no R bridge."""
    extra = {n for n in _imports(FGSEA_R) if n != "giap" and n not in sys.stdlib_module_names}
    assert extra == set(), extra
    assert {"subprocess"} <= _imports(FGSEA_R)


def test_check_imports_still_passes_on_src() -> None:
    sys.path.insert(0, str(REPO / "scripts" / "ci"))
    try:
        import check_imports
    finally:
        sys.path.pop(0)
    assert check_imports.violations(REPO / "src") == []
    assert "rpy2" in check_imports.FORBIDDEN


@pytest.mark.parametrize("name", ["docker", "rpy2", "spython"])
def test_no_runtime_sdk_anywhere_in_src(name: str) -> None:
    for path in (REPO / "src" / "giap").rglob("*.py"):
        assert name not in _imports(path), f"{path} imports {name}"
