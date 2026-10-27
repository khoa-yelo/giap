"""User guide and get-started pages (M11-08 acceptance 2: the ✓ checks, examples 1-4, option drift,
the local-executor rule and the row-level-data guard). Spec: docs/specs/M11-08-user-guide.md."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
GUIDE = [f"get-started/{p}.md" for p in ("install", "setup-databases", "manual-databases",
                                                 "quickstart")] + [
    f"user-guide/{p}.md"
    for p in ("inputs", "run", "outputs", "interpreting-results", "plots-and-report", "compare",
              "makedb", "convert", "licences", "runtimes")
]  # fmt: skip
TUTORIALS = ["tutorials/ecoli-k12.md", "tutorials/nissle-pilot.md"]  # guarded from the stubs on
CHECKED = GUIDE + TUTORIALS
PLOTS = {"V01": "circos", "V02": "abundance_landscape", "V03": "locus", "V05": "mapping_sankey",
         "V06": "mapping_fdr", "V07": "obs_vs_exp", "V08": "cog_abundance", "V10": "correlogram",
         "V11": "operon_pairs", "V12": "variants", "V13": "novel_orfs", "V14": "enrichment",
         "V15": "determinants", "V16": "compare", "V17": "qc",
         "V18": "stat_diagnostics"}  # fmt: skip
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")  # validation/pilot/metrics.py (M10-05)
FENCE = re.compile(r"^```.*?^```", re.S | re.M)


def _text(rel: str) -> str:
    return (DOCS / rel).read_text(encoding="utf-8")


def _cli_options() -> dict[str, set[str]]:
    opts: dict[str, set[str]] = {}
    for block in _text("reference/cli.md").split("\n## gippa ")[1:]:
        name, _, body = block.partition("\n")
        rows = [ln.split("|")[1] for ln in body.splitlines() if ln.startswith("| `")]
        found = {o for r in rows for o in re.findall(r"`(--[a-z0-9-]+)`", r)}
        opts[name.strip()] = found | {"--help"}
    return opts


def _flat_nav(node: object) -> list[str]:
    if isinstance(node, str):
        return [node]
    items = node.values() if isinstance(node, dict) else node
    return [p for item in items for p in _flat_nav(item)]  # type: ignore[union-attr]


class _Loader(yaml.SafeLoader):
    """mkdocs.yml with its ``!relative`` tags read as plain values."""


_Loader.add_multi_constructor("!", lambda loader, suffix, node: None)


def test_pages_exist_and_are_in_nav() -> None:
    cfg = yaml.load((REPO / "mkdocs.yml").read_text(encoding="utf-8"), Loader=_Loader)
    nav = _flat_nav(cfg["nav"])
    for rel in GUIDE:
        assert "Filled by M11-08" not in _text(rel), rel
        assert rel in nav, rel


@pytest.mark.parametrize("rel", CHECKED)
def test_options_exist_in_cli_reference(rel: str) -> None:
    """Option drift: every `--option` and every `gippa <cmd> ...` line names options of cli.md."""
    cli, text = _cli_options(), _text(rel)
    every = set().union(*cli.values())
    for raw in re.findall(r"`([^`]+)`", FENCE.sub("", text)):  # a span may wrap a line
        span = " ".join(raw.split())
        if span.startswith("--"):
            assert re.match(r"--[a-z0-9-]+", span)[0] in every, (rel, span)  # type: ignore[index]
    lines = [m for b in FENCE.findall(text) for m in b.replace("\\\n", " ").splitlines()]
    spans = [" ".join(s.split()) for s in re.findall(r"`(gippa\s[^`]+)`", FENCE.sub("", text))]
    for cmdline in [ln.split(" #")[0].strip() for ln in lines] + spans:
        words = cmdline.split()
        if len(words) > 1 and words[0] == "gippa" and words[1] in cli:
            for opt in re.findall(r"(?<![\w-])--[a-z0-9-]+", cmdline):
                assert opt in cli[words[1]], (rel, words[1], opt)


@pytest.mark.parametrize("rel", CHECKED)
def test_no_local_executor_or_unmeasured_runtime(rel: str) -> None:
    text = _text(rel)
    assert "--executor local" not in text, rel  # §6.4: developer guide only (example 3)
    for line in text.splitlines():  # example 4: §13.2's estimates only inside a quote of it
        if re.search(r"\b35\s*min", line):
            assert line.startswith(">") and "13.2" in line, (rel, line)


def test_tutorials_hold_no_row_level_data() -> None:
    """The M10-05 guard regex on docs/tutorials/ (no peptide-like string of 7+ residues)."""
    for page in sorted((DOCS / "tutorials").glob("*.md")):
        assert not PEPTIDE_LIKE.search(page.read_text(encoding="utf-8")), page


def test_licences_name_the_restrictions() -> None:
    lines = _text("user-guide/licences.md").splitlines()
    # example 1, anchored to the "What you need to know" bullets (not only the tables)
    assert any(ln.startswith("- **geNOMAD") and "non-commercial" in ln for ln in lines)
    assert any(ln.startswith("- **KEGG") and "academic" in ln for ln in lines)
    text = "\n".join(lines)
    for needle in ("UniOP has no licence", "never redistributes", "DIA-NN is proprietary",
                   "separate programs", "AGPL", '--8<-- "NOTICE.md"'):  # fmt: skip
        assert needle in text, needle


def test_install_setup_run_outputs_contents() -> None:
    install = _text("get-started/install.md")
    assert 'requires-python = ">=3.11"' in (REPO / "pyproject.toml").read_text(encoding="utf-8")
    assert "| Python | Python ≥ 3.11 |" in install
    setup, cli = _text("get-started/setup-databases.md"), _cli_options()
    for flag in ("--bakta", "--genomad", "--promotech-models", "--go", "--kegg", "--cog",
                 "--contaminants", "--antismash"):  # fmt: skip
        assert f"`{flag}`" in setup and flag in cli["setup"], flag
    assert "db_manifest.json" in setup and "Users never download RegulonDB" in setup
    run = _text("user-guide/run.md")
    for opt in ("--modules", "--publication", "--seed", "--resume", "--strict", "--threads"):
        assert f"`{opt}`" in run and opt in cli["run"], opt  # example 2 for --threads
    for code in range(6):
        assert re.search(rf"^\| {code} \| ", run, re.M), code
    out = _text("user-guide/outputs.md")
    for needle in ("../reference/tables.md", "1-based closed", "0-based half-open", "_part1",
                   "_part2", "Is_circular=true", "modules_off"):  # fmt: skip
        assert needle in out, needle
    assert "partial_excluded" in _text("user-guide/interpreting-results.md")


def test_plot_catalogue_matches_registry() -> None:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import base

    load_plot_modules()
    assert {k: s.name for k, s in base.registered().items()} == PLOTS
    text = _text("user-guide/plots-and-report.md")
    for pid, name in PLOTS.items():
        row = rf"^\| {pid} \| `{name}` \|.*plots/data/{pid}_{name}\.tsv"
        assert re.search(row, text, re.M), pid


def test_install_pins_and_profiles_match_the_code() -> None:
    from gippa.cli import _PROFILES

    install = _text("get-started/install.md")
    boot = (REPO / "scripts" / "bootstrap_farmshare.sh").read_text(encoding="utf-8")
    nf = re.search(r"NEXTFLOW_VERSION:-([0-9.]+)", boot)[1]  # type: ignore[index]
    jdk = re.search(r"openjdk=(\d+)", boot)[1]  # type: ignore[index]
    assert f"nextflow={nf} openjdk={jdk}" in install and f"| Nextflow | {nf}" in install
    assert "--override-channels" in install
    section = install.split("## Execution profiles")[1].split("\n## ")[0]
    assert set(re.findall(r"^\| `([a-z_]+)` \|", section, re.M)) == set(_PROFILES)


def test_exit_code_meanings_match_exitcodes() -> None:
    from gippa.exitcodes import ExitCode

    words = {"SUCCESS": "success", "INVALID_INPUT": "invalid input",
             "ENVIRONMENT": "environment error", "WARNINGS_STRICT": "`--strict`",
             "STEP_FAILED": "step failed", "INTERNAL": "internal error"}  # fmt: skip
    run = _text("user-guide/run.md")
    for code in ExitCode:
        row = re.search(rf"^\| {int(code)} \| (.*)$", run, re.M)
        assert row and words[code.name] in row[1], code.name


def test_coordinate_example_is_consistent() -> None:
    """The BED rows of outputs.md follow from its genes.tsv rows (§6.5: 0-based half-open, split
    at the origin into _part1/_part2)."""
    out = _text("user-guide/outputs.md")
    length = int(re.search(r"##sequence-region chr1 1 (\d+)", out)[1])  # type: ignore[index]
    rows = re.findall(r"^(G_\d+)\s+chr1\s+(\d+)\s+(\d+)\s+([+-])\s+(true|false)$", out, re.M)
    assert {r[3] for r in rows} == {"+", "-"} and "true" in {r[4] for r in rows}
    for tag, start, end, strand, crosses in rows:
        s0, e = int(start) - 1, int(end)
        assert (e > length) == (crosses == "true")
        parts = [(s0, length, f"{tag}_part1"), (0, e - length, f"{tag}_part2")]
        for a, b, name in parts if e > length else [(s0, e, tag)]:
            assert re.search(rf"^chr1\s+{a}\s+{b}\s+{name}\s+\d+\s+\{strand}$", out, re.M), name
