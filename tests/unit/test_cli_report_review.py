"""PR #205 review round, ``gippa report OUT`` (M9-02 behaviour 7; GIPPA_SPEC.md §5.1, §7.9.4;
review B3 and B5 and the lead decisions of that round).

* B3: without ``--config``, ``gippa report OUT`` resolves ``OUT/resolved_config.yaml`` when it
  exists, so every alpha a card or plot title uses, ``report.formats`` and ``report.dpi`` follow the
  run; ``--set`` still applies on top; defaults only when the file is absent; an explicit
  ``--config`` wins over the run's file.
* B5: an unreadable or schema-invalid table never crashes ``gippa report``: exit 0, the report
  is written, the affected card carries a note naming the file.

Fake plots (tests/helpers/report_fixtures.py); never touches the real ``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

import re
import shutil
import struct
from pathlib import Path
from types import ModuleType

import pytest
import yaml
from click.testing import CliRunner, Result
from helpers.report_fixtures import parse, register_fake, run_out, swap_registry
from helpers.viz_out_tree import write_schema_table

_PAIR = {"gene_a": "g1", "gene_b": "g2", "uniop_class": "OP", "promoter_between_yes": False,
         "terminator_between_yes": True, "abundance_pattern": "discordant",
         "flag": "internal_terminator_supported", "lfdr": 0.5, "evidence_note": "x"}  # fmt: skip


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _report(out: Path, *args: str) -> Result:
    from gippa.cli import cli

    return CliRunner().invoke(cli, ["report", str(out), *args])


def _ok(res: Result) -> None:
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output


def _write_config(path: Path, *overrides: str) -> Path:
    """A full resolved config, as ``gippa run`` writes ``OUT/resolved_config.yaml``."""
    from gippa.config import resolve_config

    cfg = resolve_config(overrides=list(overrides))
    path.write_text(yaml.safe_dump(cfg.model_dump(mode="json"), sort_keys=True), encoding="utf-8")
    return path


def _flags(out: Path) -> None:
    rows = [{**_PAIR, "pair_id": "PAIR_1", "q_value": 0.03},
            {**_PAIR, "pair_id": "PAIR_2", "q_value": 0.08}]  # fmt: skip
    write_schema_table(out / "tables" / "operon_flags.tsv", "operon_flags", rows)


def _operon_card(out: Path) -> tuple[str, str]:
    page = parse((out / "report" / "report.html").read_text(encoding="utf-8"))
    card = page.cards()["operon_flags"]
    value = " ".join(card.by_class("value")[0].text().split())
    return value, " ".join(card.text().split())


def _png_width(path: Path) -> int:
    head = path.read_bytes()[:24]
    assert head[:8] == b"\x89PNG\r\n\x1a\n", path
    return int(struct.unpack(">I", head[16:20])[0])


# ---------------------------------------------------------------- B3 the run's config


def test_b3_alpha_follows_the_runs_resolved_config(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    """The reviewer's probe: the run used ``stats.alpha.operon_flags=0.05``; ``gippa report OUT``
    with no ``--config`` counts and labels at 0.05, not at the default 0.1."""
    out = run_out(tmp_path)
    _write_config(out / "resolved_config.yaml", "stats.alpha.operon_flags=0.05",
                  "report.formats=[svg]", "report.dpi=50")  # fmt: skip
    _flags(out)
    _ok(_report(out))
    value, text = _operon_card(out)
    assert re.fullmatch(r"\D*1\D*", value), text
    assert "0.05" in text and "0.1" not in text, text


def test_b3_set_still_applies_over_the_runs_config(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    _write_config(out / "resolved_config.yaml", "stats.alpha.operon_flags=0.05",
                  "report.formats=[svg]", "report.dpi=50")  # fmt: skip
    _flags(out)
    _ok(_report(out, "--set", "stats.alpha.operon_flags=0.01"))
    value, text = _operon_card(out)
    assert re.fullmatch(r"\D*0\D*", value), text
    assert "0.01" in text, text


def test_b3_explicit_config_wins_over_the_runs_file(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    _write_config(out / "resolved_config.yaml", "stats.alpha.operon_flags=0.05",
                  "report.formats=[svg]", "report.dpi=50")  # fmt: skip
    _flags(out)
    other = _write_config(tmp_path / "other.yaml", "stats.alpha.operon_flags=0.2",
                          "report.formats=[svg]", "report.dpi=50")  # fmt: skip
    _ok(_report(out, "--config", str(other)))
    value, text = _operon_card(out)
    assert re.fullmatch(r"\D*2\D*", value), text
    assert "0.2" in text and "0.05" not in text, text


def test_b3_defaults_when_the_run_has_no_resolved_config(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    out = run_out(tmp_path)
    _flags(out)
    _ok(_report(out, "--set", "report.formats=[svg]", "--set", "report.dpi=50"))
    value, text = _operon_card(out)
    assert re.fullmatch(r"\D*2\D*", value), text
    assert "0.1" in text, text


def test_b3_formats_and_dpi_follow_the_runs_config(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    register_fake("V07")
    run = run_out(tmp_path / "run")
    _write_config(run / "resolved_config.yaml", "report.formats=[png]", "report.dpi=30")
    _ok(_report(run))
    assert list((run / "plots").glob("V07_*.png")), sorted(p.name for p in run.rglob("*"))
    assert not list((run / "plots").glob("V07_*.svg")), "formats=[png] wrote an SVG"
    assert not list((run / "plots").glob("V07_*.pdf"))
    explicit = run_out(tmp_path / "explicit")
    _ok(_report(explicit, "--set", "report.formats=[png]", "--set", "report.dpi=30"))
    default_dpi = run_out(tmp_path / "default")
    _ok(_report(default_dpi, "--set", "report.formats=[png]"))
    w_run = _png_width(next((run / "plots").glob("V07_*.png")))
    w_explicit = _png_width(next((explicit / "plots").glob("V07_*.png")))
    w_default = _png_width(next((default_dpi / "plots").glob("V07_*.png")))
    assert w_run == w_explicit, (w_run, w_explicit)
    assert w_run < w_default, (w_run, w_default)


def test_b3_no_config_report_equals_report_with_the_runs_config(
    fake_registry: ModuleType, tmp_path: Path
) -> None:
    """``gippa report OUT`` reproduces what ``gippa report OUT --config OUT/resolved_config.yaml``
    writes, byte for byte (the rebuilt report agrees with the run's thresholds)."""
    register_fake("V07")
    a = run_out(tmp_path / "a")
    _write_config(a / "resolved_config.yaml", "stats.alpha.operon_flags=0.05",
                  "stats.alpha.context=0.01", "report.formats=[svg]", "report.dpi=50")  # fmt: skip
    _flags(a)
    b = tmp_path / "b" / "OUT"
    shutil.copytree(a, b)
    _ok(_report(a))
    _ok(_report(b, "--config", str(b / "resolved_config.yaml")))
    html_a, html_b = (o / "report" / "report.html" for o in (a, b))
    assert html_a.read_bytes() == html_b.read_bytes()


# ---------------------------------------------------------------- B5 broken tables


@pytest.mark.parametrize(
    ("rel", "body"),
    [
        ("tables/genes.tsv", b"locus_tag\tdetected\ng1\ttrue\n"),
        ("tables/islands.tsv", b"\xff\xfe\x00\x81binary\x00\n"),
        ("tables/novel_orfs.tsv", b"orf_id\tclass\tcalled\nN1\tsorf\ttrue\n"),
    ],
    ids=["genes_missing_columns", "islands_binary", "novel_orfs_missing_columns"],
)
def test_b5_gippa_report_exits_0_and_writes_the_report_on_a_broken_table(
    fake_registry: ModuleType, tmp_path: Path, rel: str, body: bytes
) -> None:
    out = run_out(tmp_path)
    (out / rel).write_bytes(body)
    res = _report(out, "--set", "report.formats=[svg]", "--set", "report.dpi=50")
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        pytest.fail(f"gippa report raised {type(res.exception).__name__}: {res.exception}")
    assert res.exit_code == 0, res.output
    html = out / "report" / "report.html"
    assert html.is_file()
    page = parse(html.read_text(encoding="utf-8"))
    name = rel.rsplit("/", 1)[1]
    assert name in " ".join(page.section("summary").text().split())
