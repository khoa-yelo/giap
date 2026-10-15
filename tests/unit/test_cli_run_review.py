"""Code review of PR #78 (real ``giap run``, Slurm 1848591): lead decisions 3a-f.

docs/specs/M2-07-nextflow-tool-processes.md, "Lead decisions for the run-wiring follow-up", item 3:

a. ``--resume`` must cache: preflight rewrites ``work/preflight/*`` only when the content changed
   (unchanged files keep their mtime). The cached-task count of a second stub run is pinned in
   tests/integration/test_nf_run_review.py (needs Nextflow).
b. The apptainer profile binds ``params.db_dir`` read-only (``-B <db_dir>:<db_dir>:ro``) when set.
c. ``--profile`` accepts comma-joined known profiles; unknown names exit 1.
d. Preflight runs before the ``nextflow`` presence check.
e. ``db_dir`` is passed to preflight/``read_manifest`` explicitly; no ``os.environ`` mutation.
f. ``--strict`` with a warning exits 3 before launch; a missing manifest exits 2 with no launch.

A fake ``nextflow`` on ``PATH`` records its argv; nothing real is launched.
"""

from __future__ import annotations

import json
import os
import random
import re
import stat
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner

from giap.exitcodes import ExitCode
from giap.preflight.environment import MANIFEST, requirements

REPO = Path(__file__).resolve().parents[2]
NF_CONFIG = REPO / "workflow" / "nextflow.config"
FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""
ALL_MODULES = ["core", "mobile", "bgc", "function"]


def _write_db(path: Path, manifest: bool = True) -> Path:
    path.mkdir(parents=True)
    if manifest:
        data = {
            r.name: {"version": r.min_version or "1.0"}
            for mode in ("targeted", "genome")
            for r in requirements(ALL_MODULES, mode)
        }
        (path / MANIFEST).write_text(json.dumps(data))
    return path


def _fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            fh.writelines(seq[i : i + 60] + "\n" for i in range(0, len(seq), 60))
    return path


def _seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    rng = random.Random(3)
    return _fasta(tmp_path / "in" / "genome.fa", [("chr1 [topology=circular]", _seq(rng, 6000))])


@pytest.fixture
def proteins(tmp_path: Path) -> Path:
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    recs = [(f"prot{i}", "M" + "".join(rng.choice(aa) for _ in range(150))) for i in range(3)]
    return _fasta(tmp_path / "in" / "proteins.faa", recs)


@pytest.fixture
def db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    path = _write_db(tmp_path / "db_env")
    monkeypatch.setenv("GIAP_DB_DIR", str(path))
    return path


def _install_fake_nf(bindir: Path) -> None:
    bindir.mkdir(parents=True, exist_ok=True)
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "fakebin"
    _install_fake_nf(bindir)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    return args


def _run(argv: list[str]) -> tuple[int, str]:
    from giap.cli import cli

    result = CliRunner().invoke(cli, ["run", *argv])
    out = result.output + (result.stderr if result.stderr_bytes is not None else "")
    if result.exception is not None and not isinstance(result.exception, SystemExit):
        raise result.exception
    return result.exit_code, out


def _main(argv: list[str]) -> int:
    """Through the console-script entry point (usage errors → exit 1, §5.1)."""
    from giap.cli import main

    with pytest.raises(SystemExit) as exc:
        main(["run", *argv])
    code = exc.value.code
    return int(code) if code is not None else 0


def _profile_arg(fake_nf: Path) -> str:
    args = fake_nf.read_text().splitlines()
    return args[args.index("-profile") + 1]


# ------------------------------------------------------ a. --resume keeps preflight mtimes


def test_a_rerun_keeps_unchanged_preflight_files_mtime(
    genome: Path, proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    """Nextflow's resume hashes input files by path, size and mtime: a rewrite of identical
    ``work/preflight/*`` content invalidates every downstream cached task."""
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out), "--resume"]
    code, text = _run(argv)
    assert code == ExitCode.SUCCESS, text
    work = out / "work" / "preflight"
    files = sorted(p for p in work.rglob("*") if p.is_file())
    assert files, "first run wrote no work/preflight files"
    old = 1_000_000_000  # 2001-09-09: any rewrite moves the mtime away from it
    for p in files:
        os.utime(p, (old, old))
    fake_nf.unlink()
    code, text = _run(argv)
    assert code == ExitCode.SUCCESS, text
    assert fake_nf.exists(), "second run did not launch"
    moved = {str(p.relative_to(out)): p.stat().st_mtime for p in files if p.stat().st_mtime != old}
    assert not moved, f"unchanged preflight files rewritten on --resume: {moved}"


def test_a_changed_input_rewrites_the_preflight_file(
    proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    """Content-addressed, not frozen: a changed genome yields a new normalised genome."""
    rng = random.Random(11)
    g = _fasta(tmp_path / "in" / "g.fa", [("chr1 [topology=circular]", _seq(rng, 6000))])
    out = tmp_path / "OUT"
    argv = ["--genome", str(g), "--proteins", str(proteins), "--outdir", str(out), "--resume"]
    code, text = _run(argv)
    assert code == ExitCode.SUCCESS, text
    norm = out / "work" / "preflight" / "genome.norm.fa"
    before = norm.read_text()
    _fasta(g, [("chr1 [topology=circular]", _seq(rng, 6000))])
    code, text = _run(argv)
    assert code == ExitCode.SUCCESS, text
    assert norm.read_text() != before, "changed genome not re-normalised"


# ------------------------------------------------------ b. apptainer binds params.db_dir


def _apptainer_profile_text() -> str:
    text = NF_CONFIG.read_text()
    m = re.search(r"^\s*apptainer\s*\{", text, re.M)
    assert m is not None, "no apptainer profile in workflow/nextflow.config"
    depth, i = 0, m.end() - 1
    for j in range(i, len(text)):
        depth += {"{": 1, "}": -1}.get(text[j], 0)
        if depth == 0:
            return text[i : j + 1]
    raise AssertionError("unbalanced apptainer profile")


def test_b_apptainer_profile_binds_db_dir_read_only() -> None:
    """Config text: the profile builds a ``-B <db_dir>:<db_dir>:ro`` bind from ``params.db_dir``
    (``apptainer.runOptions`` or ``containerOptions``). Evaluated end to end in
    tests/integration/test_nf_run_review.py."""
    block = _apptainer_profile_text()
    assert "params.db_dir" in block, "apptainer profile does not reference params.db_dir"
    assert re.search(r"runOptions|containerOptions", block)
    assert re.search(r"-B\s*\$\{?params\.db_dir\}?:\$\{?params\.db_dir\}?:ro", block), (
        "no read-only `-B ${params.db_dir}:${params.db_dir}:ro` bind in the apptainer profile"
    )


# ------------------------------------------------------ c. comma-joined profiles


@pytest.mark.parametrize(
    "profile", ["farmshare_local,apptainer", "farmshare,apptainer", "farmshare_local"]
)
def test_c_comma_joined_known_profiles_are_accepted(
    profile: str, genome: Path, proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out)]
    code, text = _run([*argv, "--profile", profile])
    assert code == ExitCode.SUCCESS, text
    assert _profile_arg(fake_nf) == profile


@pytest.mark.parametrize("profile", ["bogus", "farmshare,bogus", "farmshare,,apptainer"])
def test_c_unknown_profile_exits_1_without_launch(
    profile: str, genome: Path, proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out)]
    assert _main([*argv, "--profile", profile]) == ExitCode.INVALID_INPUT
    assert not fake_nf.exists(), "nextflow launched with an unknown profile"


# ------------------------------------------------------ d. preflight before nextflow check


def test_d_bad_genome_without_nextflow_exits_1_with_preflight_report(
    proteins: Path, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "emptybin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    bad = tmp_path / "in" / "bad.fa"
    bad.parent.mkdir(exist_ok=True)
    bad.write_text("this is not a FASTA file\n")
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(bad), "--proteins", str(proteins), "--outdir", str(out)])
    assert code == ExitCode.INVALID_INPUT, text
    assert (out / "preflight" / "preflight.json").is_file(), "no preflight report written"
    report = json.loads((out / "preflight" / "preflight.json").read_text())
    assert "error" in json.dumps(report), report


def test_d_good_input_without_nextflow_still_exits_2(
    genome: Path, proteins: Path, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "emptybin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    out = tmp_path / "OUT"
    code, text = _run(["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out)])
    assert code == ExitCode.ENVIRONMENT, text
    assert "nextflow" in text.lower(), text


# ------------------------------------------------------ e. explicit db_dir, no env mutation


def _spy_read_manifest(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    import giap.preflight.environment as env_mod
    import giap.preflight.run as run_mod

    real = env_mod.read_manifest
    calls: list[dict[str, Any]] = []

    def spy(db_dir: str | None) -> Any:
        calls.append({"arg": db_dir, "env": os.environ.get("GIAP_DB_DIR")})
        return real(db_dir)

    monkeypatch.setattr(env_mod, "read_manifest", spy)
    if hasattr(run_mod, "read_manifest"):
        monkeypatch.setattr(run_mod, "read_manifest", spy)
    return calls


@pytest.mark.parametrize("env_set", [True, False], ids=["env-set", "env-unset"])
def test_e_db_dir_reaches_preflight_without_touching_the_environment(
    env_set: bool,
    genome: Path,
    proteins: Path,
    fake_nf: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    env_db = tmp_path / "db_env_no_manifest"
    env_db.mkdir()
    if env_set:
        monkeypatch.setenv("GIAP_DB_DIR", str(env_db))
    else:
        monkeypatch.delenv("GIAP_DB_DIR", raising=False)
    before = os.environ.get("GIAP_DB_DIR")
    given = _write_db(tmp_path / "db_given")
    calls = _spy_read_manifest(monkeypatch)
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out)]
    code, text = _run([*argv, "--db-dir", str(given)])
    assert code == ExitCode.SUCCESS, text
    assert os.environ.get("GIAP_DB_DIR") == before
    assert calls, "preflight never read the manifest"
    for call in calls:
        assert call["arg"] is not None and Path(call["arg"]).resolve() == given.resolve(), call
        assert call["env"] == before, f"$GIAP_DB_DIR mutated during preflight: {call}"


# ------------------------------------------------------ f. strict / missing manifest


def test_f_strict_warning_exits_3_before_launch(
    proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    rng = random.Random(13)
    # `chr|1` is sanitised to `chr_1`: a genome.ids warning (and 6 kb is below min_total_bp_warn)
    g = _fasta(tmp_path / "in" / "w.fa", [("chr|1 [topology=circular]", _seq(rng, 6000))])
    out = tmp_path / "OUT"
    argv = ["--genome", str(g), "--proteins", str(proteins), "--outdir", str(out)]
    code, text = _run([*argv, "--strict"])
    assert code == ExitCode.WARNINGS_STRICT, text
    assert not fake_nf.exists(), "nextflow launched after a --strict warning"
    fake_nf.unlink(missing_ok=True)
    code, text = _run(argv)
    assert code == ExitCode.SUCCESS, text  # without --strict the same input launches


def test_f_missing_manifest_exits_2_without_launch(
    genome: Path, proteins: Path, fake_nf: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    nodb = _write_db(tmp_path / "db_empty", manifest=False)
    monkeypatch.setenv("GIAP_DB_DIR", str(nodb))
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out)]
    code, text = _run(argv)
    assert code == ExitCode.ENVIRONMENT, text
    assert MANIFEST in text, text
    assert not fake_nf.exists(), "nextflow launched without a database manifest"


def test_f_missing_manifest_in_given_db_dir_exits_2_without_launch(
    genome: Path, proteins: Path, db: Path, fake_nf: Path, tmp_path: Path
) -> None:
    """``$GIAP_DB_DIR`` has a manifest but ``--db-dir`` (which wins) does not."""
    nodb = _write_db(tmp_path / "db_empty", manifest=False)
    out = tmp_path / "OUT"
    argv = ["--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out)]
    code, text = _run([*argv, "--db-dir", str(nodb)])
    assert code == ExitCode.ENVIRONMENT, text
    assert not fake_nf.exists(), "nextflow launched without a database manifest"
