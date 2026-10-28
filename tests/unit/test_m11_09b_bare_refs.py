"""Bare image references in ``workflow/conf/containers.config`` (M11-09b; ADR-0059).

Spec docs/specs/M11-09b-bare-image-refs.md: rules 1-5, examples 1 and 4-9, acceptance criteria 1,
2, 5, 6, 7 and 9, and the trace-row note of the test plan. Defaults of the spec questions: H4 is
strict (Q1), Docker Hub short names (Q2), readers accept an optional ``docker://`` (Q3).

No network and no real ``$GIPPA_ROOT``: ``update_digests.py`` edits a copy in ``tmp_path``,
``resolve_digest.resolve`` runs against a monkeypatched ``urllib.request.urlopen``, readiness runs
on the ``helpers.release_readiness`` fixture repository, and ``ToolRunner`` reads a throw-away
``GIPPA_REPO``.

Where the spec leaves the wording open: H4's reason for a non-``withName`` ``container =`` value is
checked only for ``docker:// prefix:`` plus the image name (the spec writes ``<ref>``, which may or
may not keep the prefix).
"""

from __future__ import annotations

import importlib.util
import io
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers.fake_env import image_file_name as fake_image_file_name
from helpers.manifest_data import trace_row
from helpers.release_readiness import Fixture, row, text_of
from hypothesis import given
from hypothesis import strategies as st

from gippa.annotate.tools import read_containers_config
from gippa.preflight.environment import _image_file_name

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "workflow" / "conf" / "containers.config"
UPDATE_DIGESTS = REPO / "scripts" / "release" / "update_digests.py"
RESOLVE_DIGEST = REPO / "scripts" / "resolve_digest.py"

BARE = re.compile(r"[a-z0-9][a-z0-9._/-]*@sha256:[0-9a-f]{64}")
PIN_VALUE = re.compile(r"(container\s*=\s*')(?:docker://)?")

# The pins on origin/main at 613d653c (AC 1: process set and digests unchanged).
PINNED = {
    "BAKTA": "oschwengers/bakta@sha256:"
    "1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18",
    "TRANSTERMHP": "quay.io/biocontainers/transtermhp@sha256:"
    "fc6fe0404be9340a883ec7a7b90415dfc3826c2bca7333a0314335c7bfee2d03",
    "GENOMAD": "quay.io/biocontainers/genomad@sha256:"
    "d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c",
    "ANTISMASH": "antismash/standalone@sha256:"
    "d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d",
    "MINIPROT": "quay.io/biocontainers/miniprot@sha256:"
    "2eb53fea53743b56b1c10ac4b202469b92db959ffbb8bafd83156fdbde80e22f",
    "BLAST": "quay.io/biocontainers/blast@sha256:"
    "9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c",
    "TBLASTN": "quay.io/biocontainers/blast@sha256:"
    "9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c",
}
TAG_LINES = [
    "// Tags the digests were resolved from (2026-10-14):",
    "//   oschwengers/bakta:v1.12.2",
    "//   quay.io/biocontainers/transtermhp:2.09--h9948957_2",
    "//   quay.io/biocontainers/genomad:1.12.0--pyhdfd78af_0",
    "//   antismash/standalone:8.0.4",
    "//   quay.io/biocontainers/miniprot:0.18--h577a1d6_0",
    "//   quay.io/biocontainers/blast:2.16.0--h66d330f_5",
]
# The six Apptainer cache files already in $NXF_APPTAINER_CACHEDIR (evidence 4; AC 2).
CACHE_NAMES = {
    "oschwengers-bakta@sha256-1e24ae2f435c46f876a322298ab96d994d15aadc2824e1f16a77fa9a1bc07f18.img",
    "quay.io-biocontainers-transtermhp@sha256-"
    "fc6fe0404be9340a883ec7a7b90415dfc3826c2bca7333a0314335c7bfee2d03.img",
    "quay.io-biocontainers-genomad@sha256-"
    "d54652ff868af72604ea0e54c541c5c92168c7ac821a422de05dcbad65509e5c.img",
    "antismash-standalone@sha256-"
    "d3f390e687cfe02f1387be40256dddbf678e91bc56b96c36da07b56ec8bda29d.img",
    "quay.io-biocontainers-miniprot@sha256-"
    "2eb53fea53743b56b1c10ac4b202469b92db959ffbb8bafd83156fdbde80e22f.img",
    "quay.io-biocontainers-blast@sha256-"
    "9dfc69f990c0aeb936276ee591ed32919a79f46dfa34060055e05a050a17959c.img",
}
BAKTA_D = PINNED["BAKTA"].split(":")[1]


def _text() -> str:
    return CONFIG.read_text(encoding="utf-8")


def _code_lines(text: str) -> list[str]:
    return [ln for ln in text.splitlines() if not ln.lstrip().startswith("//")]


def prefixed(text: str) -> str:
    """The config with ``docker://`` on every ``container = '…'`` value (the pre-M11-09b form)."""
    return PIN_VALUE.sub(r"\1docker://", text)


def bare(text: str) -> str:
    """The config with no ``docker://`` on any ``container = '…'`` value."""
    return PIN_VALUE.sub(r"\1", text)


def _load(path: Path, name: str) -> ModuleType:
    """A script as a module, by path only (nothing from its directory goes on ``sys.path``)."""
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None, f"{path} missing"
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ AC 1, example 1: the config


def test_every_pin_is_a_bare_digest_reference() -> None:
    pins = read_containers_config(CONFIG)
    assert pins, "no pins read from containers.config"
    not_bare = {p: ref for p, ref in pins.items() if not BARE.fullmatch(ref)}
    assert not not_bare, f"pins that are not bare digest references: {not_bare}"


def test_no_code_line_of_the_config_contains_docker_scheme() -> None:
    hits = [ln for ln in _code_lines(_text()) if "docker://" in ln]
    assert hits == [], hits


def test_bakta_line_is_example_1() -> None:
    line = f"  withName: 'BAKTA' {{ container = '{PINNED['BAKTA']}' }}"
    assert line in _text().splitlines()


def test_process_set_digests_and_tag_comments_are_unchanged() -> None:
    pins = read_containers_config(CONFIG)
    assert set(pins) == set(PINNED)
    have = {ref.removeprefix("docker://") for ref in pins.values()}
    assert have == set(PINNED.values())
    lines = _text().splitlines()
    start = lines.index(TAG_LINES[0])
    assert lines[start : start + len(TAG_LINES)] == TAG_LINES


def test_header_comment_states_the_bare_form_and_cites_the_adr() -> None:
    header = "\n".join(_text().splitlines()[:5])
    assert "ADR-0059" in header, header
    assert "docker://<image>" not in header, "the header still says images run as docker://…"


# ------------------------------------------------------------------ AC 2, example 4: cache names


@pytest.mark.parametrize("name_of", [_image_file_name, fake_image_file_name])
def test_example4_both_forms_give_the_same_cache_name(name_of: Any) -> None:
    ref = PINNED["BAKTA"]
    want = f"oschwengers-bakta@sha256-{BAKTA_D}.img"
    assert name_of(ref) == want
    assert name_of(f"docker://{ref}") == want


def test_cache_names_from_new_config_equal_those_from_the_old(tmp_path: Path) -> None:
    old = tmp_path / "old.config"
    old.write_text(prefixed(_text()), encoding="utf-8")
    old_refs = read_containers_config(old).values()
    assert all(r.startswith("docker://") for r in old_refs), "old config built wrongly"
    new_names = {_image_file_name(r) for r in read_containers_config(CONFIG).values()}
    old_names = {_image_file_name(r) for r in old_refs}
    assert new_names == old_names == CACHE_NAMES


@given(
    image=st.from_regex(
        r"[a-z0-9][a-z0-9._-]{0,8}(/[a-z0-9][a-z0-9._-]{0,8}){0,2}", fullmatch=True
    ),
    digest=st.from_regex(r"[0-9a-f]{64}", fullmatch=True),
)
def test_cache_name_ignores_the_prefix_property(image: str, digest: str) -> None:
    ref = f"{image}@sha256:{digest}"
    name = _image_file_name(ref)
    assert name == _image_file_name(f"docker://{ref}") == fake_image_file_name(ref)
    assert name.endswith(f"@sha256-{digest}.img") and "/" not in name and ":" not in name


# ------------------------------------------------------------------ AC 9: no reader needs docker://

# `container = 'docker://` (optionally escaped quote) with the scheme as a required literal, not in
# an optional group such as `(?:docker://)?`.
_REQUIRES_PREFIX = re.compile(r"container\s*=\s*\\?['\"]docker://")


def test_no_code_requires_the_docker_prefix_to_find_a_pin() -> None:
    hits = []
    for top in ("scripts", "src", "workflow"):
        for path in sorted((REPO / top).rglob("*")):
            if not path.is_file() or path.suffix in {".pyc", ".img", ".gz"}:
                continue
            try:
                text = path.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                continue
            for n, ln in enumerate(text.splitlines(), 1):
                if ln.lstrip().startswith(("//", "#")):
                    continue
                if _REQUIRES_PREFIX.search(ln):
                    hits.append(f"{path.relative_to(REPO)}:{n}: {ln.strip()}")
    assert hits == [], "\n".join(hits)


# ------------------------------------------------------------------ example 5: update_digests.py

OSTIR_E = "e" * 64
OSTIR_REF = f"ghcr.io/khoa-yelo/gippa-ostir@sha256:{OSTIR_E}"
COMMIT = "0123456789abcdef0123456789abcdef01234567"


def _update(tmp_path: Path, config_text: str, digests: dict[str, str]) -> tuple[Any, Path]:
    cfg = tmp_path / "containers.config"
    cfg.write_text(config_text, encoding="utf-8")
    dj = tmp_path / "digests.json"
    dj.write_text(json.dumps(digests), encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, "-I", str(UPDATE_DIGESTS), "--digests", str(dj), "--config", str(cfg),
         "--date", "2026-10-25", "--source-commit", COMMIT],
        capture_output=True, text=True, timeout=60, check=False,
    )  # fmt: skip
    return proc, cfg


def test_example5_update_digests_writes_a_bare_line(tmp_path: Path) -> None:
    proc, cfg = _update(tmp_path, _text(), {"OSTIR": OSTIR_REF})
    assert proc.returncode == 0, proc.stderr
    lines = cfg.read_text(encoding="utf-8").splitlines()
    assert f"  withName: 'OSTIR' {{ container = '{OSTIR_REF}' }}" in lines
    assert not any("docker://" in ln for ln in _code_lines("\n".join(lines)))


def test_update_digests_replaces_a_legacy_prefixed_line_in_place(tmp_path: Path) -> None:
    legacy = f"docker://ghcr.io/khoa-yelo/gippa-ostir@sha256:{'a' * 64}"
    old = f"  withName: 'OSTIR' {{ container = '{legacy}' }}"
    lines = _text().splitlines()
    end = max(i for i, ln in enumerate(lines) if ln == "}")
    lines.insert(end, old)
    before = "\n".join(lines) + "\n"
    proc, cfg = _update(tmp_path, before, {"OSTIR": OSTIR_REF})
    assert proc.returncode == 0, proc.stderr
    after = cfg.read_text(encoding="utf-8").splitlines()
    new = f"  withName: 'OSTIR' {{ container = '{OSTIR_REF}' }}"
    assert after.count(new) == 1 and old not in after
    # in place: still the last line of the process block (the tag comment above it grows by one)
    close = max(i for i, ln in enumerate(after) if ln == "}")
    assert after.index(new) == close - 1, "the legacy line is replaced in place"
    # every pin line other than OSTIR's is byte-identical
    pins_before = [ln for ln in before.splitlines() if "withName" in ln and "OSTIR" not in ln]
    pins_after = [ln for ln in after if "withName" in ln and "OSTIR" not in ln]
    assert pins_before == pins_after


def test_update_digests_pin_line_is_bare_and_readable(tmp_path: Path) -> None:
    mod = _load(UPDATE_DIGESTS, "gippa_m11_09b_update_digests")
    line = mod.pin_line("OSTIR", OSTIR_REF)
    assert line == f"  withName: 'OSTIR' {{ container = '{OSTIR_REF}' }}"
    cfg = tmp_path / "c.config"
    cfg.write_text(f"process {{\n{line}\n}}\n", encoding="utf-8")
    assert read_containers_config(cfg) == {"OSTIR": OSTIR_REF}


def test_update_digests_still_rejects_a_prefixed_input(tmp_path: Path) -> None:
    proc, cfg = _update(tmp_path, _text(), {"OSTIR": f"docker://{OSTIR_REF}"})
    assert proc.returncode == 1 and "OSTIR" in proc.stderr
    assert cfg.read_text(encoding="utf-8") == _text(), "config left untouched"


# ------------------------------------------------------------------ example 6: resolve_digest.py

BLAST_D = PINNED["BLAST"].split(":")[1]


class _Resp(io.BytesIO):
    def __init__(self, headers: dict[str, str]) -> None:
        super().__init__(b"")
        self.headers = headers

    def __enter__(self) -> _Resp:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()


@pytest.fixture
def resolve_mod(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    mod = _load(RESOLVE_DIGEST, "gippa_m11_09b_resolve_digest")
    calls: list[str] = []

    def urlopen(req: Any, timeout: float = 0) -> _Resp:
        url = req if isinstance(req, str) else req.full_url
        calls.append(url)
        assert url == "https://quay.io/v2/biocontainers/blast/manifests/2.16.0--h66d330f_5", url
        return _Resp({"Docker-Content-Digest": f"sha256:{BLAST_D}"})

    monkeypatch.setattr(mod.urllib.request, "urlopen", urlopen)
    mod.calls = calls  # type: ignore[attr-defined]
    return mod


@pytest.mark.parametrize("prefix", ["", "docker://"])
def test_example6_resolve_prints_a_bare_reference(resolve_mod: ModuleType, prefix: str) -> None:
    out = resolve_mod.resolve(f"{prefix}quay.io/biocontainers/blast:2.16.0--h66d330f_5")
    assert out == f"quay.io/biocontainers/blast@sha256:{BLAST_D}"
    assert resolve_mod.calls, "the fake registry was never asked"  # type: ignore[attr-defined]


def test_resolve_main_prints_the_bare_reference(
    resolve_mod: ModuleType, capsys: pytest.CaptureFixture[str]
) -> None:
    assert resolve_mod.main(["quay.io/biocontainers/blast:2.16.0--h66d330f_5"]) == 0
    assert capsys.readouterr().out.split() == [f"quay.io/biocontainers/blast@sha256:{BLAST_D}"]


@pytest.mark.parametrize(
    "ref",
    ["quay.io/biocontainers/blast:2.16.0--h66d330f_5", "oschwengers/bakta:v1.12.2", "ubuntu:24.04"],
)
def test_split_ref_accepts_a_prefixed_input(ref: str) -> None:
    mod = _load(RESOLVE_DIGEST, "gippa_m11_09b_resolve_digest_split")
    assert mod.split_ref(f"docker://{ref}") == mod.split_ref(ref)


# ------------------------------------------------------------------ examples 7-9: readiness H4


def _h4(tmp_path: Path, edit: Any = None) -> dict[str, Any]:
    fx = Fixture.build(tmp_path)
    p = fx.repo / "workflow/conf/containers.config"
    text = bare(p.read_text(encoding="utf-8"))  # the M11-09b config, whatever the checkout holds
    if edit is not None:
        text = edit(text)
    p.write_text(text, encoding="utf-8")
    _, doc = fx.compute()
    return row(doc, "H4")


def _bakta(value: str) -> Any:
    def edit(text: str) -> str:
        new = re.sub(
            r"(withName: 'BAKTA' \{ container = ')[^']+'", lambda m: f"{m.group(1)}{value}'", text
        )
        assert value in new
        return new

    return edit


def test_example7_h4_passes_on_the_bare_config(tmp_path: Path) -> None:
    h4 = _h4(tmp_path)
    assert h4["status"] == "pass", h4
    assert h4["value"].startswith("7 digests; "), h4


def test_h4_expected_text_names_the_bare_form(tmp_path: Path) -> None:
    h4 = _h4(tmp_path)
    assert "@sha256:<64 hex>" in h4["threshold"] and "no docker:// prefix" in h4["threshold"], h4
    assert not h4["threshold"].startswith("docker://"), h4


def test_example8_h4_fails_on_a_prefixed_pin(tmp_path: Path) -> None:
    h4 = _h4(tmp_path, _bakta(f"docker://{PINNED['BAKTA']}"))
    assert h4["status"] == "fail", h4
    assert "docker:// prefix: BAKTA" in text_of(h4), h4
    assert "tag, not digest" not in text_of(h4), h4


def test_example9_h4_tag_still_fails_as_tag(tmp_path: Path) -> None:
    h4 = _h4(tmp_path, _bakta("oschwengers/bakta:v1.12.2"))
    assert h4["status"] == "fail" and "tag, not digest: BAKTA" in text_of(h4), h4
    # the threshold itself says "no docker:// prefix"; the reason must not
    assert "docker:// prefix" not in text_of({k: v for k, v in h4.items() if k != "threshold"}), h4


def test_h4_prefixed_tag_is_still_a_tag(tmp_path: Path) -> None:
    h4 = _h4(tmp_path, _bakta("docker://oschwengers/bakta:v1.12.2"))
    assert h4["status"] == "fail" and "tag, not digest: BAKTA" in text_of(h4), h4


def test_h4_mixed_config_names_each_prefixed_process(tmp_path: Path) -> None:
    def edit(text: str) -> str:
        for p in ("GENOMAD", "MINIPROT"):
            text = text.replace(
                f"withName: '{p}' {{ container = '", f"withName: '{p}' {{ container = 'docker://"
            )
        return text

    h4 = _h4(tmp_path, edit)
    t = text_of(h4)
    assert h4["status"] == "fail", h4
    assert "docker:// prefix: GENOMAD" in t and "docker:// prefix: MINIPROT" in t, h4
    assert "docker:// prefix: BAKTA" not in t, h4


def test_h4_other_container_value_with_prefix_fails(tmp_path: Path) -> None:
    d = f"sha256:{'c' * 64}"
    h4 = _h4(tmp_path, lambda t: t + f"\nprocess.container = 'docker://ubuntu@{d}'\n")
    assert h4["status"] == "fail", h4
    t = text_of(h4)
    assert "docker:// prefix: " in t and f"ubuntu@{d}" in t, h4


def test_h4_other_container_value_bare_passes(tmp_path: Path) -> None:
    d = f"sha256:{'c' * 64}"
    h4 = _h4(tmp_path, lambda t: t + f"\nprocess.container = 'ubuntu@{d}'\n")
    assert h4["status"] == "pass", h4


# ------------------------------------------------------------------ AC 7: ToolRunner.command


@pytest.mark.parametrize("form", ["bare", "prefixed"])
def test_toolrunner_finds_the_pinned_image_in_either_form(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, form: str
) -> None:
    from gippa.setupdb import ToolRunner

    repo = tmp_path / "repo"
    (repo / "workflow" / "conf").mkdir(parents=True)
    text = bare(_text()) if form == "bare" else prefixed(_text())
    (repo / "workflow" / "conf" / "containers.config").write_text(text, encoding="utf-8")
    cache = tmp_path / "cache"
    cache.mkdir()
    img = cache / f"oschwengers-bakta@sha256-{BAKTA_D}.img"
    img.write_bytes(b"")
    monkeypatch.setenv("GIPPA_REPO", str(repo))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "root"))
    monkeypatch.setenv("MAMBA_ROOT_PREFIX", str(tmp_path / "mamba"))
    monkeypatch.setenv("NXF_APPTAINER_CACHEDIR", str(cache))
    argv = ["bakta_db", "download", "--type", "light"]
    cmd = ToolRunner().command("bakta", argv)
    assert cmd[:2] == ["apptainer", "exec"], cmd
    i = cmd.index(str(img))
    assert cmd[i + 1 :] == argv, cmd


# ------------------------------------------------------------------ trace row with a bare value


def test_trace_row_with_bare_container_gives_the_same_digest(tmp_path: Path) -> None:
    from gippa.io import manifest_sources

    reps = []
    for value in (PINNED["MINIPROT"], f"docker://{PINNED['MINIPROT']}"):
        d = tmp_path / "work" / str(len(reps))
        d.mkdir(parents=True)
        (d / "versions.yml").write_text('"ANNOTATE:MINIPROT":\n    miniprot: "0.18"\n')
        r = trace_row("ANNOTATE:MINIPROT (genome)", workdir=str(d), container=value)
        reps.append(manifest_sources.tool_reports([r], tmp_path / "work", tmp_path / "locks"))
    want = f"sha256:{PINNED['MINIPROT'].split(':')[1]}"
    assert [rep[0]["container_digest"] for rep in reps] == [want, want], reps


def test_fixture_helpers_are_not_vacuous(tmp_path: Path) -> None:
    """``bare`` and ``prefixed`` really change the 7 pin values (guards the tests above)."""
    assert len(re.findall(r"container = 'docker://", prefixed(_text()))) == 7
    assert "docker://" not in "\n".join(_code_lines(bare(prefixed(_text()))))
    shutil.copyfile(CONFIG, tmp_path / "c.config")
    assert len(read_containers_config(tmp_path / "c.config")) == 7
