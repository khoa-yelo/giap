"""Genome-wide Promotech is rejected at config load (M2-13 Q2, lead decision 4).

The config key is ``annotation.promoters.mode`` (GIAP_SPEC.md §5.2; ``giap run --promoters``);
the lead decision names it ``promotech.mode``. ``mode: genome`` must fail model validation with a
message containing "not implemented", so ``giap run`` and ``giap preflight`` exit 1.
"""

from __future__ import annotations

import contextlib
from pathlib import Path

import pytest
from click.testing import CliRunner, Result

GENOME_CFG = "annotation:\n  promoters:\n    mode: genome\n"


def _cfg(tmp_path: Path) -> Path:
    p = tmp_path / "giap.yaml"
    p.write_text(GENOME_CFG)
    return p


def _fasta(tmp_path: Path) -> Path:
    p = tmp_path / "g.fa"
    p.write_text(">c1\n" + "ACGT" * 50 + "\n")
    return p


def _text(r: Result) -> str:
    out = r.output
    with contextlib.suppress(ValueError, AttributeError):  # stderr may be mixed into output
        out += r.stderr
    return out


def _says_genome_not_implemented(r: Result) -> bool:
    """One message line naming the mode, not the generic ``giap run`` stub."""
    return any("genome" in ln and "not implemented" in ln for ln in _text(r).splitlines())


def test_default_mode_targeted_is_accepted() -> None:
    from giap.config import GiapConfig

    assert GiapConfig().annotation.promoters.mode == "targeted"


def test_model_rejects_genome_mode() -> None:
    from giap.config import GiapConfig

    with pytest.raises(ValueError, match="not implemented"):
        GiapConfig.model_validate({"annotation": {"promoters": {"mode": "genome"}}})


def test_resolve_config_rejects_genome_mode_from_file(tmp_path: Path) -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError, match=r"genome.*not implemented|not implemented.*genome"):
        resolve_config(_cfg(tmp_path))


def test_resolve_config_rejects_genome_mode_from_set() -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError, match="not implemented"):
        resolve_config(None, ["annotation.promoters.mode=genome"])


def test_preflight_exits_1_on_genome_mode(tmp_path: Path) -> None:
    from giap.cli import cli

    r = CliRunner().invoke(
        cli, ["preflight", "--config", str(_cfg(tmp_path)), "--genome", str(_fasta(tmp_path))]
    )
    assert r.exit_code == 1, (r.exit_code, r.output, r.exception)
    assert _says_genome_not_implemented(r), _text(r)


def test_run_exits_1_on_promoters_genome_flag(tmp_path: Path) -> None:
    from giap.cli import cli

    args = ["run", "--promoters", "genome", "--genome", str(_fasta(tmp_path))]
    r = CliRunner().invoke(cli, [*args, "--outdir", str(tmp_path / "out")])
    assert r.exit_code == 1, (r.exit_code, r.output, r.exception)
    assert _says_genome_not_implemented(r), _text(r)
    assert not (tmp_path / "out").exists()
