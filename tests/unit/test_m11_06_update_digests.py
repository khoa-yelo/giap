"""``scripts/release/update_digests.py`` on copies of the real ``containers.config`` (M11-06).

Spec docs/specs/M11-06-ghcr-containers.md item 3, examples 1-5 and the "two processes, one image"
edge case; acceptance criterion 1. No network: every run edits a copy in a temporary directory.

CLI under test (item 3, plus two flags this file fixes so runs are deterministic):

    update_digests.py --digests digests.json --config <containers.config>
                      --date YYYY-MM-DD --source-commit <40-hex git sha>

``--date`` and ``--source-commit`` are the build date and source commit that item 3 appends to the
"Tags the digests were resolved from" comment (in CI: today's UTC date and ``$GITHUB_SHA``).
Exit 0 on success; exit 1 with a one-line reason on stderr and the config left untouched otherwise.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "update_digests.py"
CONFIG = REPO / "workflow" / "conf" / "containers.config"
ORIGINAL = CONFIG.read_text()

DATE = "2026-10-25"
COMMIT = "0123456789abcdef0123456789abcdef01234567"
SHA7 = COMMIT[:7]

# process -> image tool (item 3; second-round aliases from workflow/subworkflows/local/novel.nf)
PROCESSES = {
    "PROMOTECH": "promotech",
    "PROMOTECH_NOVEL": "promotech",
    "PROMOTER_CALCULATOR": "promoter_calculator",
    "TERMINATORNET": "terminatornet",
    "BACTERMFINDER": "bactermfinder",
    "OSTIR": "ostir",
    "OSTIR_NOVEL": "ostir",
}
CONDA_ONLY_LINE = "//   uniop, promotech, promoter_calculator, terminatornet, bactermfinder, ostir"
TAGS_HEADER = "// Tags the digests were resolved from"


def _ref(tool: str, hexd: str) -> str:
    return f"ghcr.io/khoa-yelo/gippa-{tool}@sha256:{hexd}"


def _pin_line(process: str, ref: str) -> str:
    return f"  withName: '{process}' {{ container = '{ref}' }}"


def _run(
    workdir: Path,
    digests: dict[str, str],
    config_text: str = ORIGINAL,
    *,
    write_config: bool = True,
) -> tuple[subprocess.CompletedProcess[str], Path]:
    cfg = workdir / "containers.config"
    if write_config:
        cfg.write_text(config_text)
    dj = workdir / "digests.json"
    dj.write_text(json.dumps(digests))
    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--digests",
            str(dj),
            "--config",
            str(cfg),
            "--date",
            DATE,
            "--source-commit",
            COMMIT,
        ],
        capture_output=True,
        text=True,
        cwd=workdir,
        timeout=60,
    )
    return proc, cfg


def _is_subsequence(needle: list[str], hay: list[str]) -> bool:
    it = iter(hay)
    return all(any(h == n for h in it) for n in needle)


def _process_block(lines: list[str]) -> tuple[int, int]:
    start = lines.index("process {")
    end = max(i for i, ln in enumerate(lines) if ln == "}")
    return start, end


def _pins(text: str, process: str) -> list[str]:
    return [ln for ln in text.splitlines() if re.match(rf"\s*withName: '{process}' ", ln)]


def test_script_exists() -> None:
    assert SCRIPT.is_file(), "scripts/release/update_digests.py is missing (M11-06 item 3)"


def test_real_config_has_the_lines_these_tests_rely_on() -> None:
    lines = ORIGINAL.splitlines()
    assert CONDA_ONLY_LINE in lines
    assert any(ln.startswith(TAGS_HEADER) for ln in lines)
    assert "process {" in lines
    assert len(re.findall(r"withName:", ORIGINAL)) == 7


# --- example 1 ------------------------------------------------------------------------------


def test_example1_adds_one_ostir_line_inside_process_block(tmp_path: Path) -> None:
    ref = _ref("ostir", "a" * 64)
    proc, cfg = _run(tmp_path, {"OSTIR": ref})
    assert proc.returncode == 0, proc.stderr
    new = cfg.read_text()
    lines = new.splitlines()
    pins = _pins(new, "OSTIR")
    assert pins == [_pin_line("OSTIR", ref)]
    start, end = _process_block(lines)
    assert start < lines.index(pins[0]) < end
    assert len(re.findall(r"withName:", new)) == 8


def test_example1_removes_ostir_from_conda_only_list(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"OSTIR": _ref("ostir", "a" * 64)})
    assert proc.returncode == 0, proc.stderr
    lines = cfg.read_text().splitlines()
    assert CONDA_ONLY_LINE not in lines
    assert "//   uniop, promotech, promoter_calculator, terminatornet, bactermfinder" in lines


def test_example1_keeps_every_other_line_byte_identical(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"OSTIR": _ref("ostir", "a" * 64)})
    assert proc.returncode == 0, proc.stderr
    new_lines = cfg.read_text().splitlines()
    kept = [ln for ln in ORIGINAL.splitlines() if ln != CONDA_ONLY_LINE]
    assert _is_subsequence(kept, new_lines), "an existing line was changed, dropped or reordered"
    # the 7 existing pins are untouched
    for ln in ORIGINAL.splitlines():
        if "withName:" in ln:
            assert ln in new_lines
    # everything new is the OSTIR pin, the edited conda-only list, or a comment line
    orig = set(ORIGINAL.splitlines())
    added = [ln for ln in new_lines if ln not in orig]
    for ln in added:
        assert ln.startswith("//") or ln == _pin_line("OSTIR", _ref("ostir", "a" * 64)), ln
    assert cfg.read_text().endswith("\n")


def test_example1_records_build_date_and_source_commit_under_tags_comment(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"OSTIR": _ref("ostir", "a" * 64)})
    assert proc.returncode == 0, proc.stderr
    lines = cfg.read_text().splitlines()
    orig = set(ORIGINAL.splitlines())
    tags_at = next(i for i, ln in enumerate(lines) if ln.startswith(TAGS_HEADER))
    start, _ = _process_block(lines)
    new_comments = [(i, ln) for i, ln in enumerate(lines) if ln not in orig and ln.startswith("//")]
    new_comments = [(i, ln) for i, ln in new_comments if not ln.startswith("//   uniop")]
    assert new_comments, "no build date / source commit appended to the Tags comment"
    assert all(tags_at < i < start for i, _ in new_comments), "appended outside the Tags comment"
    blob = "\n".join(ln for _, ln in new_comments)
    assert DATE in blob
    assert SHA7 in blob
    assert "gippa-ostir" in blob


# --- example 2 ------------------------------------------------------------------------------


def test_example2_second_run_is_a_no_op(tmp_path: Path) -> None:
    digests = {"OSTIR": _ref("ostir", "a" * 64)}
    p1, cfg = _run(tmp_path, digests)
    assert p1.returncode == 0, p1.stderr
    once = cfg.read_bytes()
    p2, _ = _run(tmp_path, digests, write_config=False)
    assert p2.returncode == 0, p2.stderr
    assert cfg.read_bytes() == once


# --- example 3 ------------------------------------------------------------------------------


def test_example3_older_digest_is_replaced_in_place(tmp_path: Path) -> None:
    old = _pin_line("OSTIR", _ref("ostir", "0" * 64))
    lines = ORIGINAL.splitlines()
    g = next(i for i, ln in enumerate(lines) if "withName: 'GENOMAD'" in ln)
    lines.insert(g + 1, old)
    seeded = "\n".join(lines) + "\n"
    new_ref = _ref("ostir", "b" * 64)
    proc, cfg = _run(tmp_path, {"OSTIR": new_ref}, seeded)
    assert proc.returncode == 0, proc.stderr
    out = cfg.read_text()
    assert _pins(out, "OSTIR") == [_pin_line("OSTIR", new_ref)]
    assert "0" * 64 not in out
    out_lines = out.splitlines()
    i = out_lines.index(_pin_line("OSTIR", new_ref))
    # same position inside process { }: between the GENOMAD and ANTISMASH pins as before
    assert "withName: 'GENOMAD'" in out_lines[i - 1]
    assert "withName: 'ANTISMASH'" in out_lines[i + 1]
    start, _ = _process_block(out_lines)
    assert i - start == g + 1 - lines.index("process {")


# --- example 4 ------------------------------------------------------------------------------


def test_example4_bad_digest_exits_1_and_leaves_config(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"OSTIR": "ghcr.io/khoa-yelo/gippa-ostir@sha256:abc"})
    assert proc.returncode == 1
    assert "bad digest for OSTIR" in proc.stderr
    assert cfg.read_text() == ORIGINAL


_HEX = "0123456789abcdef"
_bad_hex = st.one_of(
    st.text(_HEX, min_size=0, max_size=63),
    st.text(_HEX, min_size=65, max_size=80),
    st.text("0123456789ABCDEF", min_size=64, max_size=64).filter(lambda s: s != s.lower()),
    st.tuples(st.text(_HEX, min_size=63, max_size=63), st.sampled_from("gzG-_ :/")).map(
        lambda t: t[0] + t[1]
    ),
)


@settings(max_examples=12, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(bad=_bad_hex)
def test_property_any_non_64_lowercase_hex_digest_is_rejected(bad: str) -> None:
    with tempfile.TemporaryDirectory() as d:
        proc, cfg = _run(Path(d), {"OSTIR": f"ghcr.io/khoa-yelo/gippa-ostir@sha256:{bad}"})
        assert proc.returncode == 1
        assert "bad digest for OSTIR" in proc.stderr
        assert cfg.read_text() == ORIGINAL


def test_tag_instead_of_digest_is_rejected(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"OSTIR": "ghcr.io/khoa-yelo/gippa-ostir:1.1.3"})
    assert proc.returncode == 1
    assert "bad digest for OSTIR" in proc.stderr
    assert cfg.read_text() == ORIGINAL


def test_one_bad_entry_rejects_the_whole_batch(tmp_path: Path) -> None:
    proc, cfg = _run(
        tmp_path,
        {
            "PROMOTECH": _ref("promotech", "c" * 64),
            "OSTIR": "ghcr.io/khoa-yelo/gippa-ostir@sha256:abc",
        },
    )
    assert proc.returncode == 1
    assert cfg.read_text() == ORIGINAL


# --- example 5 ------------------------------------------------------------------------------


def test_example5_uniop_is_never_pinned(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"UNIOP": _ref("uniop", "d" * 64)})
    assert proc.returncode == 1
    assert "UNIOP must not be pinned to a public image (unlicensed; R5)" in proc.stderr
    assert cfg.read_text() == ORIGINAL


# --- other input guards (chosen reading: process names and images must match item 3) ---------


def test_unknown_process_is_rejected(tmp_path: Path) -> None:
    proc, cfg = _run(tmp_path, {"NOT_A_PROCESS": _ref("ostir", "e" * 64)})
    assert proc.returncode == 1
    assert "NOT_A_PROCESS" in proc.stderr
    assert cfg.read_text() == ORIGINAL


@pytest.mark.parametrize(
    "ref",
    [
        "ghcr.io/khoa-yelo/gippa-promotech@sha256:" + "f" * 64,  # another tool's image
        "docker.io/khoa-yelo/gippa-ostir@sha256:" + "f" * 64,  # not GHCR
        "ghcr.io/someone-else/gippa-ostir@sha256:" + "f" * 64,  # not the project's owner
    ],
)
def test_image_must_be_the_process_own_ghcr_image(tmp_path: Path, ref: str) -> None:
    proc, cfg = _run(tmp_path, {"OSTIR": ref})
    assert proc.returncode == 1
    assert "OSTIR" in proc.stderr
    assert cfg.read_text() == ORIGINAL


# --- edge case: two processes, one image -----------------------------------------------------


def test_two_processes_one_image_get_the_same_digest(tmp_path: Path) -> None:
    ref = _ref("promotech", "1" * 64)
    proc, cfg = _run(tmp_path, {"PROMOTECH": ref, "PROMOTECH_NOVEL": ref})
    assert proc.returncode == 0, proc.stderr
    out = cfg.read_text()
    assert _pins(out, "PROMOTECH") == [_pin_line("PROMOTECH", ref)]
    assert _pins(out, "PROMOTECH_NOVEL") == [_pin_line("PROMOTECH_NOVEL", ref)]
    assert (
        "//   uniop, promoter_calculator, terminatornet, bactermfinder, ostir" in out.splitlines()
    )


def test_all_five_images_leave_only_uniop_conda_only(tmp_path: Path) -> None:
    digests = {p: _ref(t, format(i + 2, "x") * 64) for i, (p, t) in enumerate(PROCESSES.items())}
    digests["PROMOTECH_NOVEL"] = digests["PROMOTECH"]
    digests["OSTIR_NOVEL"] = digests["OSTIR"]
    proc, cfg = _run(tmp_path, digests)
    assert proc.returncode == 0, proc.stderr
    out = cfg.read_text()
    for p, ref in digests.items():
        assert _pins(out, p) == [_pin_line(p, ref)]
    assert len(re.findall(r"withName:", out)) == 17  # includes the 3 sharded TerminatorNet pins
    lines = out.splitlines()
    assert "//   uniop" in lines
    # fgsea's own comment is untouched (conda-only under every profile, M8-17 decision 9)
    for ln in ORIGINAL.splitlines():
        if "fgsea" in ln:
            assert ln in lines


# --- property: any subset of pins ------------------------------------------------------------


@settings(max_examples=10, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    chosen=st.dictionaries(
        st.sampled_from(sorted(PROCESSES)),
        st.text(_HEX, min_size=64, max_size=64),
        min_size=1,
    )
)
def test_property_pins_are_exact_preserving_and_idempotent(chosen: dict[str, str]) -> None:
    hexes: dict[str, str] = {}  # one digest per image: OSTIR and OSTIR_NOVEL share one (#372)
    for p, h in chosen.items():
        hexes.setdefault(PROCESSES[p], h)
    digests = {p: _ref(PROCESSES[p], hexes[PROCESSES[p]]) for p in chosen}
    with tempfile.TemporaryDirectory() as d:
        proc, cfg = _run(Path(d), digests)
        assert proc.returncode == 0, proc.stderr
        once = cfg.read_text()
        for p, ref in digests.items():
            assert _pins(once, p) == [_pin_line(p, ref)]
        shards = 3 if "TERMINATORNET" in digests else 0  # the sharded processes share its image
        assert len(re.findall(r"withName:", once)) == 7 + len(digests) + shards
        kept = [ln for ln in ORIGINAL.splitlines() if ln != CONDA_ONLY_LINE]
        assert _is_subsequence(kept, once.splitlines())
        start, end = _process_block(once.splitlines())
        for p, ref in digests.items():
            assert start < once.splitlines().index(_pin_line(p, ref)) < end
        proc2, _ = _run(Path(d), digests, write_config=False)
        assert proc2.returncode == 0, proc2.stderr
        assert cfg.read_text() == once
