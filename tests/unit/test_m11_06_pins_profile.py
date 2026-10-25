"""Pins vs the apptainer profile, and the NOTICE sentence (M11-06 items 3-4, AC 4).

Invariants on the real files (they hold before and after the digest PR):

* every process pinned in ``containers.config`` gets ``conda = ''`` under ``-profile apptainer``
  (item 4: "once pinned"); a pin without it would build the conda env inside the container run;
* UniOP is never pinned and never loses conda (R5/R11); fgsea is never pinned (M8-17 decision 9);
* the processes ``update_digests.py`` may pin exist in the workflow (PROMOTECH/OSTIR plus their
  ``_NOVEL`` aliases).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CONTAINERS = (REPO / "workflow/conf/containers.config").read_text()
NEXTFLOW = (REPO / "workflow/nextflow.config").read_text()
NOTICE = (REPO / "NOTICE.md").read_text()

PINNABLE = {
    "PROMOTECH",
    "PROMOTECH_NOVEL",
    "PROMOTER_CALCULATOR",
    "TERMINATORNET",
    "BACTERMFINDER",
    "OSTIR",
    "OSTIR_NOVEL",
}


def _pinned() -> set[str]:
    return set(re.findall(r"withName:\s*'([A-Z0-9_]+)'\s*\{\s*container\s*=", CONTAINERS))


def _apptainer_block() -> str:
    start = NEXTFLOW.index("apptainer {")
    depth = 0
    for i in range(start, len(NEXTFLOW)):
        if NEXTFLOW[i] == "{":
            depth += 1
        elif NEXTFLOW[i] == "}":
            depth -= 1
            if depth == 0:
                return NEXTFLOW[start : i + 1]
    raise AssertionError("unbalanced apptainer profile")


def _conda_off() -> set[str]:
    names: set[str] = set()
    for pat in re.findall(r"withName:\s*'([^']+)'\s*\{\s*conda\s*=\s*''\s*\}", _apptainer_block()):
        names |= set(pat.split("|"))
    return names


def test_every_pin_has_conda_off_under_apptainer() -> None:
    missing = _pinned() - _conda_off()
    assert not missing, f"pinned but still conda under -profile apptainer: {sorted(missing)}"


def test_conda_off_only_for_pinned_processes() -> None:
    extra = _conda_off() - _pinned()
    assert not extra, f"conda = '' without a container pin: {sorted(extra)}"


def test_uniop_and_fgsea_keep_conda() -> None:
    assert "UNIOP" not in _pinned() and "UNIOP" not in _conda_off()
    assert not {p for p in _pinned() | _conda_off() if "FGSEA" in p}
    assert "giap-uniop" not in CONTAINERS


def test_pinnable_processes_exist_in_the_workflow() -> None:
    nf = "\n".join(p.read_text() for p in (REPO / "workflow").rglob("*.nf"))
    for name in PINNABLE:
        assert re.search(rf"(process {name}\b|as {name}\b)", nf), name


GHCR_PIN = re.compile(r"container = '(?:docker://)?(ghcr\.io/[^']+)'")


def test_ghcr_pins_use_the_project_namespace() -> None:
    for ref in GHCR_PIN.findall(CONTAINERS):
        assert re.fullmatch(r"ghcr\.io/khoa-yelo/giap-[a-z_]+@sha256:[0-9a-f]{64}", ref), ref


@pytest.mark.parametrize("prefix", ["", "docker://"])
def test_ghcr_pin_pattern_finds_a_written_pin_line(prefix: str) -> None:
    """The config holds no GHCR pin yet, so the loop above can run zero times (#354 review S3)."""
    ref = "ghcr.io/khoa-yelo/giap-ostir@sha256:" + "a" * 64
    line = f"  withName: 'OSTIR' {{ container = '{prefix}{ref}' }}"
    assert GHCR_PIN.findall(line) == [ref]


def test_notice_states_ghcr_images_contain_gpl_tools_with_source_links() -> None:
    """AC 4: the GHCR images contain GPL-3.0 tools fetched from upstream, with source links."""
    paras = re.split(r"\n\s*\n", NOTICE)
    hits = [p for p in paras if re.search(r"GHCR|ghcr\.io", p)]
    assert hits, "NOTICE.md has no paragraph about the GHCR images"
    para = " ".join(hits)
    assert "GPL-3.0" in para
    for tool in ("Promotech", "Promoter Calculator", "OSTIR", "BacTermFinder"):
        assert tool in para, tool
    assert len(re.findall(r"https://", para)) >= 4, "link each GPL tool's upstream source"
