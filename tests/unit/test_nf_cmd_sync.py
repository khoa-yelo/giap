"""Python command builders vs. Nextflow ``script:`` blocks (M2-07 example 6; M2-index decision 2).

Each ``giap.annotate.<tool>_cmd`` builder and ``workflow/modules/local/<tool>/main.nf`` describe the
same command (tool images do not contain GIAP). Every flag the builder can emit must appear in the
process script. Builders landing later (M2-12…M2-18) add their pair to ``BUILDERS`` in their own PR;
``test_every_cmd_builder_has_a_sync_pair`` makes that mandatory.

Independently of the builders, the §7.1.1 essentials of every tool command are checked in the
script blocks, including TerminatorNet's substring-matching parser (§7.1.1 gotcha: only the
documented flags ``-f -s -n -o``, ``-s`` always).
"""

from __future__ import annotations

import pkgutil
import re
from collections.abc import Callable
from pathlib import Path

import pytest

import giap.annotate
from giap.annotate import bakta_cmd, terminatornet_cmd, transtermhp_cmd, uniop_cmd

REPO = Path(__file__).resolve().parents[2]
MODULES = REPO / "workflow" / "modules" / "local"


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", text)


def _script(tool: str) -> str:
    path = MODULES / tool / "main.nf"
    assert path.is_file(), f"workflow/modules/local/{tool}/main.nf missing"
    m = re.search(r"\bscript\s*:(.*?)\bstub\s*:", _strip_comments(path.read_text()), flags=re.S)
    assert m, f"{tool}: no script: block followed by stub:"
    return m.group(1)


def _has_flag(text: str, flag: str) -> bool:
    return re.search(rf"(?<![\w-]){re.escape(flag)}(?![\w-])", text) is not None


def _flags(argv: list[str]) -> set[str]:
    return {a for a in argv if re.fullmatch(r"--?[A-Za-z][\w-]*", a)}


# ---------------------------------------------------------------- builder <-> script pairs


def _bakta_variants() -> list[list[str]]:
    common = dict(
        genome=Path("genome.norm.fa"),
        outdir=Path("bakta"),
        prefix="P",
        db=Path("db"),
        threads=4,
        all_circular=False,
        replicons_tsv=Path("replicons.tsv"),
        extra_args=[],
    )
    return [
        bakta_cmd.bakta_command(**common, complete=True, regions_gff3=None),  # type: ignore[arg-type]
        bakta_cmd.bakta_command(**common, complete=False, regions_gff3=Path("regions.gff3")),  # type: ignore[arg-type]
        bakta_cmd.bakta_command(
            **{**common, "all_circular": True},  # type: ignore[arg-type]
            complete="auto",
            regions_gff3=None,
        ),
    ]


def test_bakta_script_passes_replicons_unconditionally_decision_11() -> None:
    """M2-07 decision 11: --replicons is never swapped out for --complete in the BAKTA script."""
    script = _script("bakta")
    assert not re.search(r"mode=--complete\b", script), (
        "BAKTA script replaces --replicons with --complete; decision 11 keeps --replicons always"
    )


def _uniop_variants() -> list[list[str]]:
    return [uniop_cmd.uniop_command(Path("uniop_cds.gff"), Path("out"), Path("/opt/uniop"))]


def _terminatornet_variants() -> list[list[str]]:
    home = Path("/opt/terminatornet")
    return [
        terminatornet_cmd.terminatornet_command(
            Path("genome.fa"), Path("terminatornet.tsv"), 1, home
        ),
        terminatornet_cmd.terminatornet_command(
            Path("genome.fa"), Path("terminatornet.tsv"), 8, home
        ),
    ]


def _transtermhp_variants() -> list[list[str]]:
    return [
        transtermhp_cmd.transterm_command(
            Path("genome.fa"), Path("genes.coords"), Path("/opt/data/expterm.dat")
        )
    ]


BUILDERS: dict[str, Callable[[], list[list[str]]]] = {
    "bakta": _bakta_variants,
    "uniop": _uniop_variants,  # M2-12 (merged before M2-07)
    "terminatornet": _terminatornet_variants,  # M2-15
    "transtermhp": _transtermhp_variants,  # M2-15
}


@pytest.mark.parametrize("tool", sorted(BUILDERS))
def test_every_builder_flag_appears_in_the_script(tool: str) -> None:
    script = _script(tool)
    flags = set().union(*(_flags(v) for v in BUILDERS[tool]()))
    assert flags, f"{tool}: builder produced no flags"
    missing = sorted(f for f in flags if not _has_flag(script, f))
    assert not missing, f"{tool}: builder flags missing from main.nf script: {missing}"


@pytest.mark.parametrize("tool", sorted(BUILDERS))
def test_script_starts_with_the_builders_executable(tool: str) -> None:
    exe = BUILDERS[tool]()[0][0]
    assert re.search(rf"(?m)^\s*{re.escape(exe)}\b", _script(tool)), f"{tool}: `{exe}` not run"


def test_every_cmd_builder_has_a_sync_pair() -> None:
    cmd_modules = {
        m.name[: -len("_cmd")]
        for m in pkgutil.iter_modules(giap.annotate.__path__)
        if m.name.endswith("_cmd")
    }
    assert "bakta" in cmd_modules
    assert cmd_modules <= set(BUILDERS), (
        f"builders without a main.nf sync pair: {sorted(cmd_modules - set(BUILDERS))}"
    )


# ---------------------------------------------------------------- §7.1.1 essentials per tool

ESSENTIALS: dict[str, tuple[str, list[str]]] = {
    # tool: (executable / entry point regex, flags that must appear)
    "bakta": (r"\bbakta\b", ["--db", "--output", "--prefix", "--threads", "--keep-contig-headers",
                             "--complete", "--replicons", "--regions"]),
    "uniop": (r"UniOP", ["-f", "--gff_source", "-t"]),
    "promotech": (r"promotech\.py", ["-s", "-m", "-f", "-o"]),
    "promoter_calculator": (r"\bpromoter-calculator\b", ["-i", "-o", "-j"]),
    "terminatornet": (r"TerminatorNet\.py", ["-f", "-s", "-n", "-o"]),
    "transtermhp": (r"\btransterm\b", ["-p"]),
    "ostir": (r"\bostir\b", ["-i", "-o", "-j"]),
    "genomad": (r"\bgenomad\s+end-to-end\b", ["--cleanup", "--splits", "--threads"]),
    "antismash": (r"\bantismash\b", ["--genefinding-tool", "--minimal", "--cpus", "--output-dir"]),
}  # fmt: skip

LITERALS: dict[str, list[str]] = {
    "uniop": ["prodigal_gff"],
    "promotech": ["RF-HOT"],
    "promoter_calculator": ["-j 1"],
    "terminatornet": ["genome.fa"],
    "transtermhp": ["genes.coords", "genome.fa"],
    "ostir": ["ostir_input.csv", "ostir.csv"],
    "genomad": ["--splits 8", "genomad_db"],
    "antismash": ["--genefinding-tool none"],
    "bakta": ["--keep-contig-headers"],
}


@pytest.mark.parametrize("tool", sorted(ESSENTIALS))
def test_script_has_the_spec_command(tool: str) -> None:
    exe, flags = ESSENTIALS[tool]
    script = _script(tool)
    assert re.search(exe, script), f"{tool}: entry point /{exe}/ not in script"
    missing = [f for f in flags if not _has_flag(script, f)]
    assert not missing, f"{tool}: flags missing from script: {missing}"
    for lit in LITERALS.get(tool, []):
        assert lit in script, f"{tool}: {lit!r} not in script"


@pytest.mark.parametrize(
    "tool", ["bakta", "promotech", "terminatornet", "ostir", "genomad", "antismash"]
)
def test_thread_count_comes_from_task_cpus(tool: str) -> None:
    assert "task.cpus" in _script(tool), f"{tool}: threads must be ${{task.cpus}}"


def _command_line(script: str, entry: str) -> str:
    """The (possibly backslash-continued) shell line that runs ``entry`` (not its version call)."""
    joined = re.sub(r"\\\n\s*", " ", script)
    lines = [ln for ln in joined.splitlines() if re.search(entry, ln) and "version" not in ln]
    assert len(lines) == 1, f"expected one line running {entry}, got {lines}"
    return lines[0]


def test_terminatornet_passes_only_documented_flags() -> None:
    """§7.1.1: the parser matches flags by substring; pass only -f -s -n -o, always -s."""
    line = _command_line(_script("terminatornet"), r"TerminatorNet\.py")
    flags = re.findall(r"(?<![\w$/{.-])-{1,2}[A-Za-z][\w-]*", line)
    assert sorted(set(flags)) == ["-f", "-n", "-o", "-s"], f"TerminatorNet flags: {flags}"
    assert re.search(r"(?<![\w-])-f\s+genome\.fa(?!\S)", line), "genome must be staged as genome.fa"


def test_promoter_calculator_runs_one_region_per_call() -> None:
    """§7.1.1 / G11: exactly one region per call, -j 1."""
    line = _command_line(_script("promoter_calculator"), r"promoter-calculator")
    assert re.search(r"(?<![\w-])-j\s+1(?!\d)", line)
    assert len(re.findall(r"(?<![\w-])-i(?![\w-])", line)) == 1
