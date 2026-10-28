"""The report's "command line" provenance row carries no absolute path in any form (M9-02
behaviour 6, GIPPA_SPEC.md §6.6; #323 review B1 and S1).

Rule pinned here (the review's fix, lead brief for #323):

* an absolute path (``/``-rooted) under the manifest's ``outdir`` shows as ``OUT/<relative>``
  (``OUT`` for the outdir itself, trailing ``/`` or not);
* any other absolute path shows as its file name (``PurePosixPath.name``);
* the path may stand alone, follow the last ``=`` of a token (``--outdir=/x``,
  ``--set k=/x``, ``--set=k=/x``), be one item of a comma list (``a.fq,/x/b.fq``, also after an
  ``=``), or be attached to a short option (``-o/x``, shown as ``-o<shown>``);
* relative arguments are unchanged.

The manifest keeps the full ``command_line``; only the report row is redacted.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

OUTDIR = "/abs/run1"

#: an absolute path in a rendered command line: one starting a token, following ``=`` or ``,``,
#: or attached to a short option (``-o/x``)
ABS = re.compile(r"(?:^|[\s=,])/|(?:^|\s)-[A-Za-z]/")


def _cmd(args: list[str], outdir: str | None = OUTDIR) -> str:
    from gippa.io.run_view import RunView
    from gippa.report.content import provenance

    doc: dict[str, Any] = {"command_line": ["gippa", "run", *args]}
    if outdir is not None:
        doc["outdir"] = outdir
    view = RunView.from_inputs({})  # no files: only the manifest doc is read
    lines = dict(provenance(view, doc, None)["lines"])
    row = lines["command line"]
    assert row.startswith("gippa run"), row
    return str(row)[len("gippa run") :].strip()


# ---------------------------------------------------------------- B1: the leaking forms


@pytest.mark.parametrize(
    ("args", "shown"),
    [
        # --set KEY=VALUE, two tokens (shared_options; config path keys such as bakta_json)
        (["--set", "db=/home/u/db"], "--set db=db"),
        (["--set", "annotation.precomputed.bakta_json=/home/u/x/sim.json"],
         "--set annotation.precomputed.bakta_json=sim.json"),
        (["--set", "k=/abs/run1/tables/x.tsv"], "--set k=OUT/tables/x.tsv"),
        (["--set", "k=/abs/run1"], "--set k=OUT"),
        (["--set", "k=/home/u/db/"], "--set k=db"),
        # --set=KEY=VALUE, one token: the value after the last "="
        (["--set=k=/home/u/c"], "--set=k=c"),
        (["--set=k=/abs/run1/sub/c"], "--set=k=OUT/sub/c"),
        # comma lists, alone and after "="
        (["--reads", "a.fq,/home/u/b.fq"], "--reads a.fq,b.fq"),
        (["--reads", "/home/u/a.fq,/abs/run1/b.fq"], "--reads a.fq,OUT/b.fq"),
        (["--reads=a.fq,/home/u/b.fq"], "--reads=a.fq,b.fq"),
        (["--set", "k=/home/u/a.fa,/abs/run1/b.fa"], "--set k=a.fa,OUT/b.fa"),
        # an attached short option
        (["-o/home/u/x"], "-ox"),
        (["-o/abs/run1/x"], "-oOUT/x"),
    ],
)  # fmt: skip
def test_b1_leaking_forms_are_redacted(args: list[str], shown: str) -> None:
    got = _cmd(args)
    assert got == shown
    assert not ABS.search(got), got
    assert "/home/u" not in got and OUTDIR not in got


# ---------------------------------------------------------------- forms that already work


@pytest.mark.parametrize(
    ("args", "shown"),
    [
        (["--outdir", "/abs/run1"], "--outdir OUT"),
        (["--outdir", "/abs/run1/"], "--outdir OUT"),
        (["--outdir=/abs/run1/x"], "--outdir=OUT/x"),
        (["--outdir=/abs/run1"], "--outdir=OUT"),
        # a sibling with the same string prefix is not under OUT
        (["--genome", "/abs/run1x/y"], "--genome y"),
        (["--genome", "/home/u/data/"], "--genome data"),
        (["--genome", "/home/u/my data/genome v2.fa"], "--genome genome v2.fa"),
        (["--genome", "/abs/run1/sub dir/f.fa"], "--genome OUT/sub dir/f.fa"),
        (["--genome", "/home/u/g.fa", "--outdir", "/abs/run1"], "--genome g.fa --outdir OUT"),
    ],
)
def test_forms_that_already_work(args: list[str], shown: str) -> None:
    got = _cmd(args)
    assert got == shown
    assert not ABS.search(got), got


@pytest.mark.parametrize(
    "args",
    [
        ["--genome", "genome.fa", "--peptides", "peptides.tsv", "--seed", "7"],
        ["--genome", "data/g.fa"],
        ["--set", "k=v"],
        ["--set", "annotation.bakta.db=light"],
        ["--set=k=rel/x.tsv"],
        ["--set", "modules=[core,function]"],
        ["--reads", "a.fq,b.fq"],
        ["-o", "rel/x"],
        ["-n2"],
        ["--outdir", "OUT"],
        ["--genome", "~/g.fa"],
        ["--config=~/x.yaml"],
    ],
)
def test_relative_arguments_are_unchanged(args: list[str]) -> None:
    assert _cmd(args) == " ".join(args)


def test_no_outdir_in_the_manifest_gives_file_names() -> None:
    assert _cmd(["--outdir", "/abs/run1", "--set", "k=/abs/run1/x"], outdir=None) == (
        "--outdir run1 --set k=x"
    )


# ---------------------------------------------------------------- the page itself


def test_report_html_command_row_has_no_absolute_path(tmp_path: Path) -> None:
    """Every B1 form at once, through ``build_report`` on a written ``run_manifest.json``."""
    from helpers.report_content import manifest, write_manifest
    from helpers.report_fixtures import cfg, parse, run_out

    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    out = run_out(tmp_path)
    secret = "/home/someone/private"
    line = ["gippa", "run", "--genome", f"{secret}/g.fa", "--outdir", str(out),
            "--set", f"db={secret}/db", f"--set=k={secret}/c.yaml",
            "--reads", f"a.fq,{secret}/b.fq", f"-o{secret}/x", f"--outdir={out}/sub"]  # fmt: skip
    write_manifest(out, manifest(command_line=line, outdir=str(out)))
    page = build_report(RunView.from_outdir(out), {}, {}, cfg())
    text = " ".join(parse(page).section("provenance").text().split())
    assert "gippa run --genome g.fa --outdir OUT --set db=db --set=k=c.yaml" in text, text
    assert "--reads a.fq,b.fq -ox --outdir=OUT/sub" in text, text
    assert secret not in page
    assert str(tmp_path) not in page


# ---------------------------------------------------------------- property


_SEG = st.text(alphabet="abcdefgxyz0123456789_.-", min_size=1, max_size=8).filter(
    lambda s: s not in (".", "..") and not s.startswith("-")
)


@st.composite
def _abs_arg(draw: st.DrawFn) -> tuple[str, str]:
    """``(token, the absolute path inside it)`` in one of the B1 forms."""
    under = draw(st.booleans())
    segs = draw(st.lists(_SEG, min_size=1, max_size=4))
    path = "/".join([OUTDIR if under else "/home/u", *segs])
    if draw(st.booleans()):
        path += "/"
    rel = draw(_SEG)
    key = draw(st.sampled_from(["k", "annotation.precomputed.bakta_json", "a.b"]))
    form = draw(st.sampled_from(["bare", "eq", "set_eq", "comma", "eq_comma", "short"]))
    token = {
        "bare": path,
        "eq": f"--genome={path}",
        "set_eq": f"--set={key}={path}",
        "comma": f"{rel},{path}",
        "eq_comma": f"{key}={rel},{path}",
        "short": f"-o{path}",
    }[form]
    return token, path


@settings(max_examples=150, deadline=None)
@given(st.lists(_abs_arg(), min_size=1, max_size=5))
def test_property_no_absolute_path_survives(parts: list[tuple[str, str]]) -> None:
    got = _cmd([t for t, _ in parts])
    assert not ABS.search(got), got
    assert "/home/u" not in got and OUTDIR not in got, got


@settings(max_examples=100, deadline=None)
@given(st.lists(_SEG, min_size=1, max_size=6))
def test_property_relative_tokens_are_unchanged(segs: list[str]) -> None:
    args = ["/".join(segs), "--set", f"k={'/'.join(segs)}", ",".join(segs)]
    assert _cmd(args) == " ".join(args)
