"""Content-addressed image keys and tags (M11-14 item 2, examples 1-4 and 2b, AC 1; ADR-0060 §1).

Spec docs/specs/M11-14-container-install.md. Everything runs on temporary trees; nothing touches a
registry or the network.

Interface fixed by these tests. The key must work from an installed wheel, where ``scripts/`` is
not shipped, so ``key()`` lives in ``src/gippa/io/image_key.py`` (spec amendment from #370);
``scripts/release/image_key.py`` and readiness H4 are thin callers of it.

``gippa.io.image_key``
    ``root`` is the directory holding ``workflow/``; it is optional everywhere and defaults to
    ``gippa.runner.nextflow.workflow_main().parent.parent`` looked up at call time (the checkout,
    or the installed package whose bundled ``workflow/`` M11-01 ships).
    ``key_inputs(tool, root) -> list[str]``
        repo-relative POSIX paths (relative to ``root``, the directory holding ``workflow/``) of
        every file under ``workflow/containers/<tool>/`` (recursive), ``workflow/envs/<tool>.yml``,
        ``workflow/envs/locks/<tool>.explicit.txt`` and ``<tool>.pip.txt`` when present, and for
        ``fgsea`` also ``workflow/bin/gippa_fgsea.R``; sorted by code point (``sorted`` on ``str``).
    ``key_text(tool, root) -> str``
        one ``<path>\\t<mode>\\t<sha256 hex of the bytes>\\n`` line per input, in that order;
        ``<mode>`` is ``100755`` when any exec bit of the file is set, else ``100644`` (git's two
        modes, without git; Docker ``COPY`` preserves the bit).
    ``key(tool, root) -> str``
        the first 12 hex of sha256 over ``key_text`` encoded as UTF-8.
    ``tag(version, key) -> str``
        ``<version>-<key>``; ``git:<40 hex>`` is written ``git-<first 8 hex>``, every character
        outside ``[A-Za-z0-9_.-]`` becomes ``_``; a result that does not fully match
        ``[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}`` raises ``ValueError``.
    ``uncovered_copy_sources(tool, root) -> list[str]``
        the build-context sources of each ``COPY`` in ``workflow/containers/<tool>/Dockerfile``
        (flags such as ``--chown`` skipped, ``COPY --from=…`` ignored, shell and JSON forms) that
        are not key inputs, in file order.
    ``image_cache(environ, db_dir) -> Path``
        the image-cache helper (tests in ``test_m11_14_image_cache.py``).

``scripts/release/image_key.py [--root DIR] [--images PATH] (TOOL ... | --plan | --check)``
    ``--root`` defaults to the repository, ``--images`` to ``<root>/workflow/conf/images.json``
    (an object ``{tool: {"key", "tag", "ref", "processes", "built_from", "licence"}}``).
    ``TOOL ...`` prints one ``<tool>\\t<key>\\t<tag>`` line per tool (version from
    ``gippa.annotate.tools.TOOLS``). The GHCR tools are the directories under
    ``workflow/containers/`` other than ``uniop`` (never published).
    ``--plan`` prints ``will rebuild: <t> (<old key> → <new key>)`` per GHCR tool whose key
    differs from its ``images.json`` key (or has no entry), or the single line
    ``will rebuild: none``; exit 0.
    ``--check`` recomputes the key of every ``images.json`` entry and exits 1 with
    ``stale image: <tool> (inputs changed since <key>)`` (one line per stale tool) when it
    differs; entries that match, and GHCR tools without an entry, pass (exit 0).
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "image_key.py"
TAG_RE = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}")
PROMOTECH_SHA = "56251ad9b883ef831b4753fc623d5ec970fe65e0"


def _ik():  # type: ignore[no-untyped-def]
    import gippa.io.image_key as c

    return c


def _write(root: Path, rel: str, data: bytes | str, mode: int = 0o644) -> Path:
    p = root / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_bytes(data.encode() if isinstance(data, str) else data)
    p.chmod(mode)
    return p


OSTIR_FILES = {
    "workflow/containers/ostir/Dockerfile": "FROM scratch\nCOPY workflow/envs/ostir.yml /tmp/e\n",
    "workflow/envs/ostir.yml": "name: gippa-ostir\ndependencies:\n  - ostir==1.1.3\n",
    "workflow/envs/locks/ostir.explicit.txt": "@EXPLICIT\nhttps://x/a.conda#00\n",
    "workflow/envs/locks/ostir.pip.txt": "ostir==1.1.3 --hash=sha256:00\n",
}
DISTRACTORS = {
    "workflow/containers/promotech/Dockerfile": "FROM scratch\n",
    "workflow/envs/promotech.yml": "name: gippa-promotech\n",
    "workflow/envs/ostirx.yml": "not ostir\n",
    "workflow/envs/locks/ostirx.explicit.txt": "not ostir\n",
    "workflow/envs/bakta.yml": "name: gippa-bakta\n",
    "workflow/bin/gippa_fgsea.R": "#!/usr/bin/env Rscript\n",
}


def _ostir_tree(root: Path, files: dict[str, str] | None = None) -> Path:
    for rel, text in {**OSTIR_FILES, **DISTRACTORS, **(files or {})}.items():
        _write(root, rel, text)
    return root


def _expected_text(root: Path, rels: list[str]) -> str:
    out = []
    for rel in sorted(rels):
        p = root / rel
        mode = "100755" if p.stat().st_mode & 0o111 else "100644"
        out.append(f"{rel}\t{mode}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n")
    return "".join(out)


def _expected_key(root: Path, rels: list[str]) -> str:
    return hashlib.sha256(_expected_text(root, rels).encode("utf-8")).hexdigest()[:12]


def _cli(*args: str) -> subprocess.CompletedProcess[str]:
    assert SCRIPT.is_file(), "scripts/release/image_key.py is missing (M11-14 14a)"
    env = {
        **os.environ,
        "PYTHONPATH": f"{REPO / 'src'}{os.pathsep}{os.environ.get('PYTHONPATH', '')}",
    }
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
        timeout=120,
        env=env,
    )


# --- example 1 ------------------------------------------------------------------------------


def test_example1_ostir_inputs_are_the_four_files(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path)
    assert _ik().key_inputs("ostir", root) == sorted(OSTIR_FILES)


def test_example1_key_text_is_four_sorted_path_mode_sha_lines(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path)
    text = _ik().key_text("ostir", root)
    assert text == _expected_text(root, list(OSTIR_FILES))
    lines = text.splitlines()
    assert len(lines) == 4
    assert [ln.split("\t")[0] for ln in lines] == sorted(OSTIR_FILES)
    assert all(re.fullmatch(r"[^\t]+\t100644\t[0-9a-f]{64}", ln) for ln in lines)


def test_example1_key_is_12_hex_of_sha256_over_the_lines(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path)
    key = _ik().key("ostir", root)
    assert re.fullmatch(r"[0-9a-f]{12}", key)
    assert key == _expected_key(root, list(OSTIR_FILES))


def test_example1_tag_is_version_dash_key(tmp_path: Path) -> None:
    c = _ik()
    root = _ostir_tree(tmp_path)
    key = c.key("ostir", root)
    assert c.tag("1.1.3", key) == f"1.1.3-{key}"


def test_example1_cli_prints_tool_key_and_tag(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path)
    key = _expected_key(root, list(OSTIR_FILES))
    proc = _cli("--root", str(root), "ostir")
    assert proc.returncode == 0, proc.stderr
    assert f"ostir\t{key}\t1.1.3-{key}" in proc.stdout.splitlines()


def test_missing_locks_are_simply_not_inputs(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path)
    (root / "workflow/envs/locks/ostir.pip.txt").unlink()
    rels = [r for r in OSTIR_FILES if not r.endswith(".pip.txt")]
    assert _ik().key_inputs("ostir", root) == sorted(rels)
    assert _ik().key("ostir", root) == _expected_key(root, rels)


def test_nested_container_files_are_inputs(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path, {"workflow/containers/ostir/sub/helper.sh": "echo hi\n"})
    assert "workflow/containers/ostir/sub/helper.sh" in _ik().key_inputs("ostir", root)


def test_paths_sort_by_code_point_not_locale(tmp_path: Path) -> None:
    extra = {
        "workflow/containers/ostir/a.txt": "a\n",
        "workflow/containers/ostir/Z.txt": "Z\n",
        "workflow/containers/ostir/_u.txt": "_\n",
        "workflow/containers/ostir/é.txt": "e\n",
    }
    root = _ostir_tree(tmp_path, extra)
    got = _ik().key_inputs("ostir", root)
    assert got == sorted([*OSTIR_FILES, *extra])
    # 'D' < 'Z' < '_' < 'a' < 'é' by code point
    names = [p.rsplit("/", 1)[1] for p in got if p.startswith("workflow/containers/ostir/")]
    assert names == ["Dockerfile", "Z.txt", "_u.txt", "a.txt", "é.txt"]
    assert _ik().key("ostir", root) == _expected_key(root, [*OSTIR_FILES, *extra])


def test_fgsea_includes_the_baked_in_script(tmp_path: Path) -> None:
    root = _ostir_tree(
        tmp_path,
        {
            "workflow/containers/fgsea/Dockerfile": "FROM scratch\n",
            "workflow/envs/fgsea.yml": "name: gippa-fgsea\n",
            "workflow/envs/locks/fgsea.explicit.txt": "@EXPLICIT\n",
        },
    )
    c = _ik()
    assert "workflow/bin/gippa_fgsea.R" in c.key_inputs("fgsea", root)
    before = c.key("fgsea", root)
    _write(root, "workflow/bin/gippa_fgsea.R", "#!/usr/bin/env Rscript\n# changed\n")
    assert c.key("fgsea", root) != before
    assert "workflow/bin/gippa_fgsea.R" not in c.key_inputs("ostir", root)


# --- example 2 and 2b -----------------------------------------------------------------------


def test_example2_one_byte_in_env_file_changes_only_that_key(tmp_path: Path) -> None:
    c = _ik()
    root = _ostir_tree(tmp_path)
    ostir, promotech = c.key("ostir", root), c.key("promotech", root)
    env = root / "workflow/envs/ostir.yml"
    data = bytearray(env.read_bytes())
    data[0] ^= 0x01
    env.write_bytes(bytes(data))
    assert c.key("ostir", root) != ostir
    assert c.key("promotech", root) == promotech


def test_example2b_exec_bit_alone_changes_the_key(tmp_path: Path) -> None:
    c = _ik()
    root = _ostir_tree(tmp_path)
    before = c.key("ostir", root)
    (root / "workflow/containers/ostir/Dockerfile").chmod(0o755)
    text = c.key_text("ostir", root)
    assert "workflow/containers/ostir/Dockerfile\t100755\t" in text
    assert c.key("ostir", root) != before
    (root / "workflow/containers/ostir/Dockerfile").chmod(0o644)
    assert c.key("ostir", root) == before


@pytest.mark.parametrize("mode", [0o654, 0o645, 0o744])
def test_any_exec_bit_is_mode_100755(tmp_path: Path, mode: int) -> None:
    root = _ostir_tree(tmp_path)
    (root / "workflow/containers/ostir/Dockerfile").chmod(mode)
    assert "workflow/containers/ostir/Dockerfile\t100755\t" in _ik().key_text("ostir", root)


def test_key_ignores_mtime_and_listing_order(tmp_path: Path) -> None:
    c = _ik()
    a = _ostir_tree(tmp_path / "a")
    b = tmp_path / "b"
    for rel, text in reversed(list({**DISTRACTORS, **OSTIR_FILES}.items())):
        _write(b, rel, text)
    assert c.key("ostir", a) == c.key("ostir", b)
    before = c.key("ostir", a)
    for rel in OSTIR_FILES:
        os.utime(a / rel, (1_000_000, 1_000_000))
    assert c.key("ostir", a) == before
    assert c.key("ostir", a) == c.key("ostir", a)


@settings(max_examples=25, deadline=None)
@given(
    files=st.dictionaries(
        st.text("abcXYZ_-.0123", min_size=1, max_size=8).filter(lambda s: s not in {".", ".."}),
        st.tuples(st.binary(max_size=64), st.booleans()),
        min_size=1,
        max_size=6,
    )
)
def test_property_key_matches_the_definition(
    tmp_path_factory: pytest.TempPathFactory, files: dict[str, tuple[bytes, bool]]
) -> None:
    root = tmp_path_factory.mktemp("tree")
    rels = []
    for name, (data, exe) in files.items():
        rel = f"workflow/containers/ostir/{name}"
        _write(root, rel, data, 0o755 if exe else 0o644)
        rels.append(rel)
    _write(root, "workflow/envs/ostir.yml", "name: gippa-ostir\n")
    rels.append("workflow/envs/ostir.yml")
    c = _ik()
    assert c.key_inputs("ostir", root) == sorted(rels)
    assert c.key_text("ostir", root) == _expected_text(root, rels)
    assert c.key("ostir", root) == _expected_key(root, rels)


# --- example 4 and tag rules ----------------------------------------------------------------


def test_example4_git_version_is_shortened() -> None:
    key = "0b7e19c4aa21"
    assert _ik().tag(f"git:{PROMOTECH_SHA}", key) == f"git-56251ad9-{key}"


def test_example4_promotech_cli_tag_uses_the_registry_version() -> None:
    proc = _cli("promotech")
    assert proc.returncode == 0, proc.stderr
    line = next(ln for ln in proc.stdout.splitlines() if ln.startswith("promotech\t"))
    _, key, tag = line.split("\t")
    assert re.fullmatch(r"[0-9a-f]{12}", key)
    assert tag == f"git-56251ad9-{key}"
    assert key == _ik().key("promotech", REPO)


@pytest.mark.parametrize(
    ("version", "tag"),
    [
        ("1.1.3", "1.1.3-0123456789ab"),
        ("1.2.4+local", "1.2.4_local-0123456789ab"),
        ("v2 beta/1", "v2_beta_1-0123456789ab"),
    ],
)
def test_tag_replaces_characters_outside_the_tag_alphabet(version: str, tag: str) -> None:
    assert _ik().tag(version, "0123456789ab") == tag


@pytest.mark.parametrize("version", [".1", "-1", "x" * 120])
def test_tag_that_cannot_match_the_grammar_is_rejected(version: str) -> None:
    with pytest.raises(ValueError):
        _ik().tag(version, "0123456789ab")


@settings(max_examples=200, deadline=None)
@given(version=st.text(min_size=1, max_size=40))
def test_property_tag_is_sanitized_version_and_key_or_an_error(version: str) -> None:
    key = "abcdef012345"
    m = re.fullmatch(r"git:([0-9a-f]{40})", version)
    base = f"git-{m.group(1)[:8]}" if m else re.sub(r"[^A-Za-z0-9_.-]", "_", version)
    want = f"{base}-{key}"
    if TAG_RE.fullmatch(want):
        assert _ik().tag(version, key) == want
    else:
        with pytest.raises(ValueError):
            _ik().tag(version, key)


# --- static COPY check ----------------------------------------------------------------------


def _ghcr_dockerfiles() -> list[str]:
    return sorted(p.parent.name for p in (REPO / "workflow/containers").glob("*/Dockerfile"))


def test_repo_has_dockerfiles_to_check() -> None:
    assert {"ostir", "promotech", "uniop"} <= set(_ghcr_dockerfiles())


@pytest.mark.parametrize("tool", _ghcr_dockerfiles())
def test_every_repo_dockerfile_copy_source_is_a_key_input(tool: str) -> None:
    assert _ik().uncovered_copy_sources(tool, REPO) == []


def test_copy_of_a_file_outside_the_inputs_is_reported(tmp_path: Path) -> None:
    docker = (
        "FROM scratch AS base\n"
        "COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/ostir.yml /tmp/env.yml\n"
        "COPY --from=base /opt/x /opt/y\n"
        "COPY workflow/bin/helper.py workflow/envs/locks/ostir.pip.txt /opt/\n"
        'COPY ["workflow/envs/promotech.yml", "/tmp/p.yml"]\n'
    )
    root = _ostir_tree(
        tmp_path,
        {"workflow/containers/ostir/Dockerfile": docker, "workflow/bin/helper.py": "x\n"},
    )
    got = _ik().uncovered_copy_sources("ostir", root)
    assert got == ["workflow/bin/helper.py", "workflow/envs/promotech.yml"]


# --- --plan (examples 2 and 3) and --check (example 2) --------------------------------------


def _repo_tree(root: Path) -> Path:
    """Copy the key inputs of every tool (containers, envs, locks, the fgsea script) to root."""
    for rel in ("workflow/containers", "workflow/envs"):
        shutil.copytree(REPO / rel, root / rel)
    (root / "workflow/bin").mkdir(parents=True)
    shutil.copy2(REPO / "workflow/bin/gippa_fgsea.R", root / "workflow/bin/gippa_fgsea.R")
    (root / "workflow/conf").mkdir(parents=True)
    return root


def _ghcr_tools(root: Path) -> list[str]:
    return sorted(
        p.name for p in (root / "workflow/containers").iterdir() if p.is_dir() and p.name != "uniop"
    )


def _pin_all(root: Path) -> dict[str, str]:
    from gippa.annotate.tools import TOOLS

    c = _ik()
    images: dict[str, dict[str, object]] = {}
    keys = {}
    for i, tool in enumerate(_ghcr_tools(root)):
        key = c.key(tool, root)
        keys[tool] = key
        images[tool] = {
            "key": key,
            "tag": c.tag(TOOLS[tool].version, key),
            "ref": f"ghcr.io/khoa-yelo/gippa-{tool}@sha256:{format(i + 1, 'x') * 64}",
            "processes": [tool.upper()],
            "built_from": "0123456789abcdef0123456789abcdef01234567",
            "licence": TOOLS[tool].licence,
        }
    (root / "workflow/conf/images.json").write_text(json.dumps(images, indent=2) + "\n")
    return keys


def _rebuild_lines(out: str) -> list[str]:
    return [ln for ln in out.splitlines() if ln.startswith("will rebuild:")]


def test_plan_with_every_key_current_rebuilds_none(tmp_path: Path) -> None:
    root = _repo_tree(tmp_path)
    _pin_all(root)
    proc = _cli("--root", str(root), "--plan")
    assert proc.returncode == 0, proc.stderr
    assert _rebuild_lines(proc.stdout) == ["will rebuild: none"]


def test_example3_comment_in_upstream_tool_env_changes_no_ghcr_key(tmp_path: Path) -> None:
    root = _repo_tree(tmp_path)
    keys = _pin_all(root)
    bakta = root / "workflow/envs/bakta.yml"
    bakta.write_text(bakta.read_text() + "# a comment\n")
    c = _ik()
    assert {t: c.key(t, root) for t in keys} == keys
    proc = _cli("--root", str(root), "--plan")
    assert proc.returncode == 0, proc.stderr
    assert _rebuild_lines(proc.stdout) == ["will rebuild: none"]
    assert _cli("--root", str(root), "--check").returncode == 0


def test_example2_plan_names_only_the_changed_tool(tmp_path: Path) -> None:
    root = _repo_tree(tmp_path)
    keys = _pin_all(root)
    env = root / "workflow/envs/ostir.yml"
    env.write_bytes(env.read_bytes() + b" ")
    new = _ik().key("ostir", root)
    assert new != keys["ostir"]
    proc = _cli("--root", str(root), "--plan")
    assert proc.returncode == 0, proc.stderr
    assert _rebuild_lines(proc.stdout) == [f"will rebuild: ostir ({keys['ostir']} → {new})"]


def test_plan_lists_a_ghcr_tool_with_no_entry(tmp_path: Path) -> None:
    root = _repo_tree(tmp_path)
    _pin_all(root)
    images = json.loads((root / "workflow/conf/images.json").read_text())
    del images["terminatornet"]
    (root / "workflow/conf/images.json").write_text(json.dumps(images))
    proc = _cli("--root", str(root), "--plan")
    assert proc.returncode == 0, proc.stderr
    lines = _rebuild_lines(proc.stdout)
    assert len(lines) == 1
    assert lines[0].startswith("will rebuild: terminatornet (")
    assert _ik().key("terminatornet", root) in lines[0]


def test_plan_never_lists_uniop(tmp_path: Path) -> None:
    root = _repo_tree(tmp_path)
    _pin_all(root)
    uniop = root / "workflow/envs/uniop.yml"
    uniop.write_bytes(uniop.read_bytes() + b"# changed\n")
    proc = _cli("--root", str(root), "--plan")
    assert proc.returncode == 0, proc.stderr
    assert _rebuild_lines(proc.stdout) == ["will rebuild: none"]


def test_example2_check_reports_the_stale_tool_only(tmp_path: Path) -> None:
    root = _repo_tree(tmp_path)
    keys = _pin_all(root)
    assert _cli("--root", str(root), "--check").returncode == 0
    env = root / "workflow/envs/ostir.yml"
    data = bytearray(env.read_bytes())
    data[-1] ^= 0x01
    env.write_bytes(bytes(data))
    proc = _cli("--root", str(root), "--check")
    assert proc.returncode == 1
    out = proc.stdout + proc.stderr
    assert f"stale image: ostir (inputs changed since {keys['ostir']})" in out
    assert out.count("stale image:") == 1


def test_check_on_the_repository_passes() -> None:
    proc = _cli("--check")
    assert proc.returncode == 0, proc.stdout + proc.stderr


# --- the installed package computes the same keys ------------------------------------------


def _installed(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Pretend to be a wheel: the bundled workflow under site/gippa/workflow (modes kept)."""
    pkg = tmp_path / "site" / "gippa"
    for rel in ("workflow/containers", "workflow/envs"):
        shutil.copytree(REPO / rel, pkg / rel)
    (pkg / "workflow/bin").mkdir(parents=True)
    shutil.copy2(REPO / "workflow/bin/gippa_fgsea.R", pkg / "workflow/bin/gippa_fgsea.R")
    (pkg / "workflow/main.nf").write_text("// bundled\n")
    monkeypatch.setattr("gippa.runner.nextflow.workflow_main", lambda: pkg / "workflow" / "main.nf")
    return pkg


def test_default_root_is_the_checkout() -> None:
    ik = _ik()
    for tool in [*_ghcr_dockerfiles(), "fgsea"]:
        assert ik.key(tool) == ik.key(tool, REPO), tool


def test_installed_package_key_equals_the_repo_key(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    ik = _ik()
    want = {t: ik.key(t, REPO) for t in [*_ghcr_dockerfiles(), "fgsea"]}
    pkg = _installed(monkeypatch, tmp_path)
    assert {t: ik.key(t) for t in want} == want
    # the default really reads the bundled copy: a change there changes the key
    env = pkg / "workflow/envs/ostir.yml"
    env.write_text(env.read_text() + "# bundled edit\n")
    assert ik.key("ostir") != want["ostir"]
    assert ik.key("ostir", REPO) == want["ostir"]


# --- review round 1 (lead): directory sources, continuation lines, unknown tools -------------


def test_directory_copy_source_with_a_non_input_file_is_reported(tmp_path: Path) -> None:
    docker = "FROM scratch\nCOPY workflow/containers/ostir/ /opt/\nCOPY workflow/bin/ /opt/\n"
    root = _ostir_tree(
        tmp_path,
        {"workflow/containers/ostir/Dockerfile": docker, "workflow/bin/other.py": "x\n"},
    )
    assert _ik().uncovered_copy_sources("ostir", root) == ["workflow/bin/"]


def test_copy_continuation_lines_are_joined(tmp_path: Path) -> None:
    docker = "FROM scratch\nCOPY workflow/envs/ostir.yml \\\n     workflow/bin/helper.py /opt/\n"
    root = _ostir_tree(
        tmp_path,
        {"workflow/containers/ostir/Dockerfile": docker, "workflow/bin/helper.py": "x\n"},
    )
    assert _ik().uncovered_copy_sources("ostir", root) == ["workflow/bin/helper.py"]


def test_copy_source_sharing_a_prefix_with_an_input_is_reported(tmp_path: Path) -> None:
    docker = "FROM scratch\nCOPY workflow/envs/ostir /opt/\n"
    root = _ostir_tree(tmp_path, {"workflow/containers/ostir/Dockerfile": docker})
    assert _ik().uncovered_copy_sources("ostir", root) == ["workflow/envs/ostir"]


def test_key_of_a_tool_with_no_inputs_is_an_error(tmp_path: Path) -> None:
    root = _ostir_tree(tmp_path)
    with pytest.raises(ValueError, match="typo"):
        _ik().key("typo", root)


def test_cli_unknown_tool_says_so() -> None:
    proc = _cli("typo")
    assert proc.returncode == 2
    assert "unknown tool typo" in proc.stderr
    assert "Traceback" not in proc.stderr
