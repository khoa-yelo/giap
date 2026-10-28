"""``diann_pilot.sbatch`` and ``fetch_pxd078459.sbatch`` run against stubs (M3-24, decisions 3-6).

Every run happens inside ``helpers.pilot_kit.pilot_sandbox``: a throw-away checkout, fake
``apptainer``/``singularity``/``wget``/``curl``/``sbatch`` first on a minimal ``PATH``, an
unreachable HTTP proxy and a 12-byte stand-in for the image. **No container is executed, no
DIA-NN runs on real data and nothing reaches the network** (decision 6).

Conventions settled by the test author:

* ``sbatch scripts/slurm/diann_pilot.sbatch A1|A2`` - the arm is the only positional argument;
  no arm or an unknown arm is a usage error (exit 2, as in deliverable 3).
* ``DIANN_SIF``, ``RAW_DIR`` and ``OUT_DIR`` are the environment overrides; ``RAW_DIR`` holds the
  timsTOF ``.d`` runs of PXD078459 and defaults under ``$GIPPA_DATA_DIR``.
* A missing image exits 2 before any container call, naming the path and ``DIANN_SIF``.
* Step 1 writes ``--out-lib``; the stub creates both that file and its ``.predicted.speclib``
  sibling (DIA-NN's own convention), and step 2 must pass one of them to ``--lib``.
* An ``--help``/``--version`` call is the smoke check decision 6 allows and is not counted.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any

import pytest
from helpers.pilot_kit import (
    BINDS,
    CFG,
    DIANN_BIN,
    PilotSandbox,
    pilot_sandbox,
    read_cfg,
)

ARMS = ("A1", "A2")


@pytest.fixture
def pilot(tmp_path: Path) -> PilotSandbox:
    return pilot_sandbox(tmp_path)


def _argv(call: dict[str, Any]) -> list[str]:
    inner = list(call["inner"])
    assert inner and inner[0] == DIANN_BIN, f"the container runs {inner[:1]}, not {DIANN_BIN}"
    return inner[1:]


def _flag_values(argv: list[str], flag: str) -> list[str]:
    return [argv[i + 1] for i, a in enumerate(argv) if a == flag and i + 1 < len(argv)]


def _steps(pilot: PilotSandbox) -> list[dict[str, Any]]:
    calls = pilot.diann_calls()
    assert calls, "no DIA-NN call reached the stubbed container; stdout/stderr of the run is above"
    return calls


def _run(pilot: PilotSandbox, arm: str = "A1", **env: str):  # type: ignore[no-untyped-def]
    res = pilot.run("diann_pilot", arm, **env)
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    return res


# ---------------------------------------------------------------- decision 3: the image


def test_missing_image_exits_2_before_any_container_call(tmp_path: Path) -> None:
    pilot = pilot_sandbox(tmp_path, sif=False)
    res = pilot.run("diann_pilot", "A1")
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    message = res.stdout + res.stderr
    assert str(pilot.sif) in message, message
    assert "DIANN_SIF" in message, message
    assert not pilot.diann_calls(), "the container was called although the image is missing"


def test_diann_sif_override_is_honoured(tmp_path: Path) -> None:
    pilot = pilot_sandbox(tmp_path)
    other = pilot.sb.tmp / "images" / "elsewhere.sif"
    other.write_text("stub sif\n")
    other.chmod(0o444)  # M10-02: the pilot runs only a read-only copy whose sha256 is the pin
    _run(pilot, "A1", DIANN_SIF=str(other))
    images = {c["image"] for c in _steps(pilot)}
    assert images == {str(other)}, images


def test_the_default_image_is_the_path_of_decision_1(tmp_path: Path) -> None:
    """With ``DIANN_SIF`` unset the pilot runs the verified copy
    ``$GIPPA_ROOT/images/diann_2.3.1.sif`` (M10-02 lead decision 1; diann_image.sbatch makes it
    from Khoa's image); without that copy it exits 2 naming the path and diann_image.sbatch."""
    pilot = pilot_sandbox(tmp_path)
    pilot.sb.env.pop("DIANN_SIF")
    res = pilot.run("diann_pilot", "A1")
    copy = pilot.sb.root / "images" / "diann_2.3.1.sif"
    assert res.returncode == 2, res.stdout + res.stderr
    assert str(copy) in res.stdout + res.stderr and "diann_image.sbatch" in res.stdout + res.stderr
    copy.parent.mkdir(parents=True, exist_ok=True)
    copy.write_text("stub sif\n")
    copy.chmod(0o444)
    res = pilot.run("diann_pilot", "A1")
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    assert {c["image"] for c in _steps(pilot)} == {str(copy)}


def test_the_container_call_binds_the_cluster_filesystems(pilot: PilotSandbox) -> None:
    _run(pilot, "A1")
    for call in _steps(pilot):
        binds = {b.split(":")[0] for b in call["binds"]}
        missing = [b for b in BINDS if b not in binds]
        assert not missing, f"missing binds {missing} in {call['argv']}"
        assert call["tool"] in {"apptainer", "singularity"}, call["tool"]


def test_nothing_is_written_into_the_repository(pilot: PilotSandbox) -> None:
    """Decision 3: the repo keeps only configs, README and scripts; outputs go to $GIPPA_ROOT."""
    _run(pilot, "A2")
    assert pilot.repo_is_clean() == "", pilot.repo_is_clean()
    for call in _steps(pilot):
        inside = [a for a in call["argv"] if a.startswith(str(pilot.sb.repo))
                  and Path(a).suffix in {".parquet", ".speclib", ".tsv"}]  # fmt: skip
        assert not inside, f"outputs inside the checkout: {inside}"


# ---------------------------------------------------------------- decision 4: two steps


@pytest.mark.parametrize("arm", ARMS)
def test_step1_predicts_then_step2_searches(pilot: PilotSandbox, arm: str) -> None:
    _run(pilot, arm)
    calls = _steps(pilot)
    assert len(calls) == 2, f"expected exactly two DIA-NN calls, got {len(calls)}"
    step1, step2 = (_argv(c) for c in calls)
    assert "--gen-spec-lib" in step1 and "--predictor" in step1, step1
    assert _flag_values(step1, "--out-lib"), step1
    assert not _flag_values(step1, "--f"), f"step 1 must not take raw files: {step1}"
    assert not _flag_values(step1, "--lib"), step1
    assert _flag_values(step2, "--lib"), step2


@pytest.mark.parametrize("arm", ARMS)
def test_step2_searches_the_library_step1_wrote(pilot: PilotSandbox, arm: str) -> None:
    _run(pilot, arm)
    first, second = _steps(pilot)
    (lib,) = set(_flag_values(_argv(second), "--lib"))
    assert lib in set(first["created"]), (
        f"step 2's --lib {lib} is not a file step 1 produced: {first['created']}"
    )


def test_step2_passes_every_raw_run(pilot: PilotSandbox) -> None:
    _run(pilot, "A1")
    step2 = _argv(_steps(pilot)[1])
    raws = _flag_values(step2, "--f")
    expected = sorted(str(p) for p in pilot.raw_dir.glob("*.d"))
    assert sorted(raws) == expected, (raws, expected)


def test_step1_is_skipped_when_the_predicted_library_exists(pilot: PilotSandbox) -> None:
    _run(pilot, "A1")
    first_run = len(_steps(pilot))
    (pilot.sb.state / "diann.jsonl").unlink()
    _run(pilot, "A1")
    again = _steps(pilot)
    assert first_run == 2 and len(again) == 1, (
        f"the second run made {len(again)} DIA-NN calls; step 1 must be skipped (decision 4)"
    )
    assert _flag_values(_argv(again[0]), "--lib"), _argv(again[0])


def test_the_search_reports_decoys_and_writes_under_out_dir(pilot: PilotSandbox) -> None:
    _run(pilot, "A1")
    step2 = _argv(_steps(pilot)[1])
    assert "--report-decoys" in step2, step2
    (out,) = set(_flag_values(step2, "--out"))
    assert out.startswith(str(pilot.out_dir)), (out, str(pilot.out_dir))
    assert Path(out).is_file(), out


def test_threads_follow_the_slurm_allocation(pilot: PilotSandbox) -> None:
    _run(pilot, "A1", SLURM_CPUS_PER_TASK="16")
    for call in _steps(pilot):
        assert _flag_values(_argv(call), "--threads") == ["16"], _argv(call)


@pytest.mark.parametrize("arm", ARMS)
def test_every_path_given_to_the_container_is_absolute(pilot: PilotSandbox, arm: str) -> None:
    """The container's working directory is not the submitter's; DIA-NN gets absolute paths."""
    _run(pilot, arm)
    for call in _steps(pilot):
        argv = _argv(call)
        for flag in ("--fasta", "--f", "--lib", "--out", "--out-lib", "--cfg"):
            for value in _flag_values(argv, flag):
                assert value.startswith("/"), f"{flag} {value} is not absolute"


@pytest.mark.parametrize(("arm", "token"), [("A1", "UP000000625"), ("A2", "strain_db")])
def test_the_arm_decides_the_database(pilot: PilotSandbox, arm: str, token: str) -> None:
    _run(pilot, arm)
    fastas = _flag_values(_argv(_steps(pilot)[0]), "--fasta")
    assert any(token in f for f in fastas), (arm, fastas)
    assert all(Path(f).is_file() for f in fastas), fastas
    other = "strain_db" if arm == "A1" else "UP000000625"
    assert not any(other in f for f in fastas), (arm, fastas)


@pytest.mark.parametrize("arm", ARMS)
def test_the_cfg_flags_reach_dia_nn(pilot: PilotSandbox, arm: str) -> None:
    """The cfg file is the source of truth for the search settings, not the script."""
    _run(pilot, arm)
    for step, call in zip(("step1", "step2"), _steps(pilot), strict=True):
        argv = _argv(call)
        for flag, args in read_cfg(CFG[arm])[step]:
            if flag == "--fasta":  # expanded by the script; checked above
                continue
            assert flag in argv, f"{arm} [{step}]: {flag} never reached DIA-NN"
            if args:
                i = argv.index(flag)
                assert list(args) == argv[i + 1 : i + 1 + len(args)], (flag, args, argv)


# ---------------------------------------------------------------- usage errors


@pytest.mark.parametrize("args", [(), ("A3",), ("a1", "A2")], ids=["none", "unknown", "two"])
def test_a_bad_arm_is_a_usage_error(pilot: PilotSandbox, args: tuple[str, ...]) -> None:
    res = pilot.run("diann_pilot", *args)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert re.search(r"A1|A2|arm|usage", res.stdout + res.stderr, re.I), res.stdout + res.stderr
    assert not pilot.diann_calls(), "DIA-NN ran despite the usage error"


def test_an_empty_raw_directory_fails_before_the_search(tmp_path: Path) -> None:
    pilot = pilot_sandbox(tmp_path, raws=0)
    res = pilot.run("diann_pilot", "A1")
    assert res.returncode != 0, f"{res.stdout}\n{res.stderr}"
    assert str(pilot.raw_dir) in res.stdout + res.stderr, res.stdout + res.stderr
    searches = [c for c in pilot.diann_calls() if _flag_values(_argv(c), "--f")]
    assert not searches, "the search ran without raw files"


# ---------------------------------------------------------------- fetch_pxd078459.sbatch


# M10-01 decision 1: the fetch reads a committed manifest (or this override) and unpacks archives;
# the stub curl returns text, so the override names files inside the three ``<run>.d`` folders.
FETCH_OVERRIDE = {
    "PXD078459_BASE_URL": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/10/PXD078459",
    "PXD078459_FILES": " ".join(f"WT_{r}.d/{f}" for r in (1, 2, 3)
                                for f in ("analysis.tdf", "analysis.tdf_bin")),
}  # fmt: skip


def test_fetch_downloads_three_runs_and_records_their_checksums(pilot: PilotSandbox) -> None:
    """§15.2 step 4 and deliverable 5: three wild-type runs under $GIPPA_DATA_DIR with
    SHA256SUMS."""
    res = pilot.run("fetch_pxd078459", **FETCH_OVERRIDE)
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    calls = pilot.fetch_calls()
    assert calls, "nothing was downloaded (the stubbed wget/curl was never called)"
    urls = [u for c in calls for u in c["urls"]]
    assert len(urls) >= 3, urls
    assert all(u.startswith(("https://", "ftps://")) for u in urls), urls
    assert not any(c["insecure"] for c in calls), "R6: TLS verification was disabled"
    dest = pilot.sb.root / "data" / "pxd078459"
    sums = dest / "SHA256SUMS"
    assert sums.is_file(), f"{sums} was not written"
    listed = {}
    for line in sums.read_text().splitlines():
        m = re.match(r"^([0-9a-f]{64})\s+\*?(\S+)$", line.strip())
        assert m, f"unparsable SHA256SUMS line: {line!r}"
        listed[Path(m.group(2)).as_posix()] = m.group(1)
    created = [Path(p) for c in calls for p in c["created"]]
    assert created, calls
    names = {(p.relative_to(dest) if p.is_relative_to(dest) else p).as_posix() for p in created}
    assert names <= set(listed), (sorted(names - set(listed)), sorted(listed))
    for rel, want in listed.items():
        assert (dest / rel).is_file(), f"{rel} is listed but not under {dest}"
        assert hashlib.sha256((dest / rel).read_bytes()).hexdigest() == want, rel


def test_fetch_writes_only_under_gippa_data_dir(pilot: PilotSandbox) -> None:
    res = pilot.run("fetch_pxd078459", **FETCH_OVERRIDE)
    assert res.returncode == 0, f"{res.stdout}\n{res.stderr}"
    assert pilot.repo_is_clean() == "", pilot.repo_is_clean()
    dest = str(pilot.sb.root / "data" / "pxd078459")
    for call in pilot.fetch_calls():
        for path in call["created"]:
            assert str(path).startswith(dest), f"{path} is outside {dest} (R11)"


def test_fetch_is_the_only_script_that_downloads(pilot: PilotSandbox) -> None:
    """Decision 6: preparing the kit never touches PXD078459 or the network from the other jobs."""
    _run(pilot, "A1")
    assert not pilot.fetch_calls(), pilot.fetch_calls()
