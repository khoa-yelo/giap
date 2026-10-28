"""M11-19: under docker every pinned image the run needs must already be local (runs never pull).

docs/specs/M11-19-nextflow-docker-pull-never.md item 2 (the docker arm of
``missing_tool_reasons``, rules 1-5 in order, probing cost, which tools), item 3 (where it does
not apply), item 4 (podman-docker), examples 8-17, acceptance criteria 2-5.

Interface pinned here (implementer, ``src/gippa/preflight/environment.py``):

- runtime docker = ``container_runtime(profile) == "docker"`` (``docker`` and no ``apptainer``),
  with or without ``conda``: a tool pinned in the bundled ``containers.config`` is present iff
  ``docker image inspect <ref>`` succeeds; otherwise one reason per tool, first match wins:
  1. ``pull_never_problem(docker_client_version())`` (Docker < 20.10 / podman < 3.0);
  3. docker not on ``PATH``;
  4. the daemon is not reachable (Docker's stderr markers, or the inspect timed out);
  5. the image is not present locally.
- ``docker_image_reason(ref) -> str``: ``""`` when local, else reason 3, 4 or 5; never pulls.
- ``docker --version`` at most once per call and only when a tool needs a docker image;
  ``image inspect`` once per present distinct ref (BLAST and TBLASTN share one), sorted tool
  order; after a no-daemon result no other ref is inspected.
- the fgsea row (``images.json`` ref) and a local UniOP docker build get rule 1 too; fgsea's
  absent text is rule 5's.

Fakes only: a fake workflow tree (``workflow_main`` patched) and a logging fake ``docker`` first
on ``PATH`` (after the autouse container guard): ``--version`` prints ``FAKE_DOCKER_VERSION``,
``image inspect <ref>`` exits 0 iff ``<ref>`` is a line of ``FAKE_DOCKER_TAGS``, or prints
``FAKE_DOCKER_INSPECT_ERR`` and exits 1. Every call is logged. No network, no real docker, never
the real ``$GIPPA_ROOT``/``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

import json
import os
import random
import stat
import subprocess
import sys
import tomllib
from pathlib import Path
from typing import Any

import pytest

from gippa.annotate.tools import TOOLS, ToolSpec, read_containers_config, required_tools
from gippa.config import resolve_config

REPO = Path(__file__).resolve().parents[2]
WORKFLOW = REPO / "workflow"
CONTAINERS_CONFIG = WORKFLOW / "conf" / "containers.config"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
SETUP = f"{PKG} setup --images --profile docker"
FGSEA_REF = f"ghcr.io/khoa-yelo/{PKG}-fgsea@sha256:" + "c3" * 32
LOCAL_TAG = f"localhost/{PKG}-uniop:{TOOLS['uniop'].version}-0123456789ab"
NO_DAEMON = (
    "Cannot connect to the Docker daemon at unix:///var/run/docker.sock. "
    "Is the docker daemon running?"
)
CONDA_REASON = "conda-only tool under a container profile without conda: {}"

DOCKER_ROWS = ["docker", "standard,docker", "docker,conda", "slurm,docker"]


def absent(ref: str) -> str:
    """Rule 5."""
    return (
        f"image {ref} is not present locally and runs never pull it (docker image inspect {ref} "
        f"failed); run: {SETUP}"
    )


def no_docker(ref: str) -> str:
    """Rule 3."""
    return (
        f"image {ref} cannot be checked: docker is not on PATH (runs never pull images); "
        f"install Docker >= 20.10, then run: {SETUP}"
    )


def no_daemon(ref: str, why: str) -> str:
    """Rule 4."""
    return (
        f"image {ref} cannot be checked: the docker daemon is not reachable ({why}); start "
        f"Docker, then run: {SETUP} if the image is missing"
    )


def too_old(family: str, xy: str, floor: str) -> str:
    """Rule 1."""
    return (
        f"{family} {xy} is older than {floor}, the first version whose docker run has "
        "--pull=never, and runs never pull images; upgrade Docker (Docker >= 20.10, or podman "
        ">= 3.0 with podman-docker)"
    )


FAKE_DOCKER = r'''#!{python} -I
"""Fake docker (tests only, M11-19): logs every call; never pulls, never runs."""
import json, os, sys
argv = sys.argv[1:]
with open(os.environ["FAKE_DOCKER_LOG"], "a") as fh:
    fh.write(json.dumps(argv) + "\n")
if argv == ["--version"]:
    line = os.environ.get("FAKE_DOCKER_VERSION", "Docker version 27.3.1, build ce12230")
    if line.startswith("podman"):
        sys.stderr.write("Emulate Docker CLI using podman. Create /etc/containers/nodocker to "
                         "quiet msg.\n")
    sys.stdout.write(line + "\n")
    sys.exit(0)
if argv[:2] == ["image", "inspect"]:
    err = os.environ.get("FAKE_DOCKER_INSPECT_ERR", "")
    if err:
        sys.stderr.write(err + "\n")
        sys.exit(1)
    refs = [a for a in argv[2:] if not a.startswith("-") and not a.startswith("{{")]
    try:
        with open(os.environ["FAKE_DOCKER_TAGS"]) as fh:
            tags = {{l.strip() for l in fh if l.strip()}}
    except OSError:
        tags = set()
    if refs and all(r in tags for r in refs):
        sys.stdout.write(json.dumps([{{"Id": "sha256:" + "ab" * 32}}]) + "\n")
        sys.exit(0)
    sys.stdout.write("[]\n")
    if os.environ.get("FAKE_DOCKER_VERSION", "").startswith("podman"):
        sys.stderr.write("Error: " + refs[-1] + ": image not known\n")
    else:
        sys.stderr.write("Error: No such image: " + refs[-1] + "\n")
    sys.exit(1)
sys.stderr.write("fake docker: unexpected call\n")
sys.exit(3)
'''


# ---------------------------------------------------------------- helpers


def _env() -> Any:
    from gippa.preflight import environment

    return environment


def _fn(name: str) -> Any:
    fn = getattr(_env(), name, None)
    assert fn is not None, f"gippa.preflight.environment.{name} is not implemented (M11-19)"
    return fn


def pins() -> dict[str, str]:
    """Registry tool name → pinned ref (TBLASTN shares blast's)."""
    return {
        p.lower(): ref
        for p, ref in read_containers_config(CONTAINERS_CONFIG).items()
        if p.lower() in TOOLS
    }


PINNED = pins()
ALL_REFS = sorted(set(PINNED.values()))


@pytest.fixture(autouse=True)
def fake(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> dict[str, Path]:
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "NXF_APPTAINER_CACHEDIR", "GIPPA_FGSEA_PREFIX",
                "GIPPA_CONTAINER_RUNTIME", "FAKE_DOCKER_VERSION",
                "FAKE_DOCKER_INSPECT_ERR"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "docker"
    exe.write_text(FAKE_DOCKER.format(python=sys.executable))
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    log, tags = tmp_path / "docker.log", tmp_path / "docker_tags.txt"
    monkeypatch.setenv("FAKE_DOCKER_LOG", str(log))
    monkeypatch.setenv("FAKE_DOCKER_TAGS", str(tags))
    wf = tmp_path / "wf"
    (wf / "conf").mkdir(parents=True)
    (wf / "main.nf").write_text("// fake\n")
    (wf / "envs").symlink_to(WORKFLOW / "envs")
    (wf / "bin").symlink_to(WORKFLOW / "bin")
    (wf / "conf" / "containers.config").write_text(CONTAINERS_CONFIG.read_text(encoding="utf-8"))
    images = json.loads((WORKFLOW / "conf" / "images.json").read_text())
    images["fgsea"] = {"key": "0123456789ab", "tag": "1.36.2-0123456789ab", "ref": FGSEA_REF,
                       "processes": [], "built_from": "0" * 40, "licence": "MIT"}  # fmt: skip
    (wf / "conf" / "images.json").write_text(json.dumps(images, indent=2) + "\n")
    monkeypatch.setattr(_env(), "workflow_main", lambda: wf / "main.nf")
    cache = tmp_path / "cache"
    cache.mkdir()
    envs = tmp_path / "envs"
    envs.mkdir()
    return {"log": log, "tags": tags, "cache": cache, "envs": envs, "bin": bindir}


def set_tags(fake: dict[str, Path], refs: list[str]) -> None:
    fake["tags"].write_text("".join(r + "\n" for r in refs))


def calls(fake: dict[str, Path]) -> list[list[str]]:
    log = fake["log"]
    return [json.loads(x) for x in log.read_text().splitlines()] if log.exists() else []


def inspected(fake: dict[str, Path]) -> list[str]:
    """The refs of every ``image inspect`` call, in order."""
    out = []
    for c in calls(fake):
        if c[:2] == ["image", "inspect"]:
            out += [a for a in c[2:] if not a.startswith("-") and not a.startswith("{{")]
    return out


def versions(fake: dict[str, Path]) -> int:
    return sum(1 for c in calls(fake) if c == ["--version"])


def assert_never_pulls(fake: dict[str, Path]) -> None:
    """AC 4: the log holds only ``--version`` and ``image inspect`` lines."""
    for c in calls(fake):
        assert c == ["--version"] or c[:2] == ["image", "inspect"], c


def default_tools(*overrides: str) -> list[ToolSpec]:
    return list(required_tools(resolve_config(overrides=list(overrides))))


def sorted_tools(tools: list[ToolSpec]) -> list[ToolSpec]:
    return sorted(tools, key=lambda t: t.name)


def reasons(
    fake: dict[str, Path],
    tools: list[ToolSpec],
    profile: str = "docker",
    local_tools: Path | None = None,
) -> dict[str, str]:
    got = _fn("missing_tool_reasons")(tools, profile, fake["cache"], None, local_tools)
    assert sorted(got) == _fn("missing_tools")(tools, profile, fake["cache"], None, local_tools)
    return dict(got)


def pinned_of(tools: list[ToolSpec]) -> dict[str, str]:
    return {t.name: PINNED[t.name] for t in tools if t.name in PINNED}


# ---------------------------------------------------------------- example 8


def test_example_8_all_refs_local_every_pinned_tool_present(fake: dict[str, Path]) -> None:
    set_tags(fake, ALL_REFS)
    tools = default_tools()
    need = pinned_of(tools)
    assert set(need) == {"bakta", "transtermhp", "genomad", "antismash", "miniprot", "blast"}
    got = _fn("missing_tool_reasons")(tools, "docker", fake["cache"], None, None)
    assert not set(got) & set(need), got
    assert versions(fake) == 1, calls(fake)
    seen = inspected(fake)
    assert sorted(seen) == ALL_REFS, f"one inspect per distinct ref (blast's once): {seen}"
    assert_never_pulls(fake)


# ---------------------------------------------------------------- example 9


def test_example_9_absent_antismash_ref(fake: dict[str, Path]) -> None:
    ref = PINNED["antismash"]
    set_tags(fake, [r for r in ALL_REFS if r != ref])
    got = reasons(fake, default_tools())
    assert got.get("antismash") == absent(ref), got
    assert not set(got) & (set(PINNED) - {"antismash"}), got
    assert_never_pulls(fake)


PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


def _tools_check(tmp_path: Path, profile: str, *overrides: str) -> Any:
    from helpers.fake_env import write_setup_manifest

    from gippa.preflight.run import PreflightInputs, run_preflight

    rng = random.Random(0)
    g = tmp_path / "in" / "g.fa"
    g.parent.mkdir(parents=True, exist_ok=True)
    g.write_text(">chr1 [topology=circular]\n" + "".join(rng.choice("ACGT") for _ in range(120_000))
                 + "\n")  # fmt: skip
    p = tmp_path / "in" / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    cfg = resolve_config(overrides=["preflight.min_n50_warn=1",
                                    "preflight.min_largest_contig_warn=1", *overrides])  # fmt: skip
    db = write_setup_manifest(tmp_path / "db")
    report, _ = run_preflight(PreflightInputs(genome=g, proteins=p), cfg, None, "run",
                              tmp_path / "OUT", db, profile=profile)  # fmt: skip
    checks = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert len(checks) == 1, checks
    return checks[0]


def test_example_9_run_preflight_reports_antismash(tmp_path: Path, fake: dict[str, Path]) -> None:
    ref = PINNED["antismash"]
    set_tags(fake, [r for r in ALL_REFS if r != ref])
    check = _tools_check(tmp_path, "docker")
    assert check.level == "error" and check.exit_code == 2, check
    assert "antismash" in check.value["missing"], check.value
    assert f"antismash: {absent(ref)}" in check.message, check.message
    assert_never_pulls(fake)


# ---------------------------------------------------------------- examples 10, 11: which tools


@pytest.mark.parametrize(
    "override",
    ["modules=[core,function]", "annotation.precomputed.regions=/x/regions.tsv"],
    ids=["example-10-modules", "example-11-precomputed-regions"],
)
def test_examples_10_11_unneeded_images_are_not_checked(
    fake: dict[str, Path], override: str
) -> None:
    ref = PINNED["antismash"]
    set_tags(fake, [r for r in ALL_REFS if r != ref])
    tools = default_tools(override)
    got = reasons(fake, tools)
    assert "antismash" not in got and "genomad" not in got, got
    seen = set(inspected(fake))
    assert ref not in seen and PINNED["genomad"] not in seen, seen
    assert seen == {PINNED[n] for n in ("bakta", "transtermhp", "miniprot", "blast")}, seen


def test_precomputed_bakta_json_skips_the_bakta_image(fake: dict[str, Path]) -> None:
    set_tags(fake, [])
    got = reasons(fake, default_tools("annotation.precomputed.bakta_json=/x/b.json"))
    assert "bakta" not in got, got
    assert PINNED["bakta"] not in inspected(fake)
    assert got.get("antismash") == absent(PINNED["antismash"]), got


# ---------------------------------------------------------------- example 12: rule 1 first


@pytest.mark.parametrize(
    ("line", "family", "xy", "floor"),
    [
        ("Docker version 19.03.15, build 99e3ed8919", "docker", "19.3", "20.10"),
        ("Docker version 1.13.1, build 7d71120/1.13.1", "docker", "1.13", "20.10"),
        ("podman version 2.2.1", "podman", "2.2", "3.0"),
    ],
)
def test_example_12_old_docker_blocks_every_pinned_tool_before_inspect(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch, line: str, family: str, xy: str,
    floor: str
) -> None:  # fmt: skip
    monkeypatch.setenv("FAKE_DOCKER_VERSION", line)
    set_tags(fake, ALL_REFS)
    tools = default_tools()
    got = reasons(fake, tools)
    for name in pinned_of(tools):
        assert got.get(name) == too_old(family, xy, floor), (name, got)
    assert inspected(fake) == [], "rule 1 comes before any image inspect"
    assert_never_pulls(fake)


def test_rule_1_reason_does_not_send_the_user_to_setup(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "Docker version 19.03.15, build 99e3ed8919")
    got = reasons(fake, [TOOLS["bakta"]])
    assert SETUP not in got["bakta"], got


def test_unparseable_version_never_blocks(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Q2: nerdctl (or any odd CLI) passes rule 1; presence still decides."""
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "nerdctl version 1.7.6")
    set_tags(fake, [PINNED["bakta"]])
    assert reasons(fake, [TOOLS["bakta"], TOOLS["miniprot"]]) == {
        "miniprot": absent(PINNED["miniprot"])
    }


# ---------------------------------------------------------------- example 13: no docker


def test_example_13_no_docker_on_path(
    tmp_path: Path, fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "nobin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    tools = default_tools()
    got = reasons(fake, tools)
    for name, ref in pinned_of(tools).items():
        assert got.get(name) == no_docker(ref), (name, got)
        assert "docker is not on PATH" in got[name] and SETUP in got[name]
    assert calls(fake) == []


# ---------------------------------------------------------------- example 14: no daemon


@pytest.mark.parametrize(
    "stderr",
    [
        NO_DAEMON,
        "permission denied while trying to connect to the Docker daemon socket at "
        'unix:///var/run/docker.sock: Get "http://%2Fvar%2Frun%2Fdocker.sock/v1.24/images/json"',
        'error during connect: Get "http://docker.example:2375/v1.47/images/x/json": dial tcp: '
        "lookup docker.example: no such host",
        # Docker 29+ (moby client v0.1.0) wording
        "permission denied while trying to connect to the docker API at "
        "unix:///var/run/docker.sock",
        "failed to connect to the docker API at unix:///var/run/docker.sock; check if the path "
        "is correct and if the daemon is running: dial unix /var/run/docker.sock: connect: no "
        "such file or directory",
        "failed to connect to the docker API at tcp://docker.example:2375: lookup docker.example: "
        "no such host",
    ],
    ids=[
        "cannot-connect",
        "permission-denied",
        "error-during-connect",
        "v29-permission-denied",
        "v29-no-socket",
        "v29-dns",
    ],
)
def test_example_14_daemon_unreachable(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch, stderr: str
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_INSPECT_ERR", stderr + "\nsecond line, not quoted")
    tools = sorted_tools(default_tools())
    got = _fn("missing_tool_reasons")(tools, "docker", fake["cache"], None, None)
    need = pinned_of(tools)
    first = sorted(need)[0]
    assert got.get(first) == no_daemon(need[first], stderr), got
    prefix = "the docker daemon is not reachable ("
    for name, ref in need.items():
        assert got.get(name, "").startswith(f"image {ref} cannot be checked: {prefix}"), got
        assert got[name].endswith(f"start Docker, then run: {SETUP} if the image is missing")
    seen = inspected(fake)
    assert seen and set(seen) == {need[first]}, f"no inspect after a no-daemon result: {seen}"
    assert len(seen) <= 2, seen
    assert_never_pulls(fake)


def test_inspect_timeout_counts_as_no_daemon(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    mod = _env()
    real = mod.subprocess.run

    def run(argv: Any, *a: Any, **k: Any) -> Any:
        if list(argv)[1:3] == ["image", "inspect"]:
            raise subprocess.TimeoutExpired(argv, k.get("timeout") or 1)
        return real(argv, *a, **k)

    monkeypatch.setattr(mod.subprocess, "run", run)
    ref = PINNED["bakta"]
    assert _fn("docker_image_reason")(ref) == no_daemon(ref, "timed out")


# ---------------------------------------------------------------- docker_image_reason


def test_docker_image_reason_present_and_absent(fake: dict[str, Path]) -> None:
    ref = PINNED["bakta"]
    set_tags(fake, [ref])
    assert _fn("docker_image_reason")(ref) == ""
    other = PINNED["genomad"]
    assert _fn("docker_image_reason")(other) == absent(other)
    assert_never_pulls(fake)


def test_docker_image_reason_without_docker(
    tmp_path: Path, fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "nobin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    ref = PINNED["bakta"]
    assert _fn("docker_image_reason")(ref) == no_docker(ref)


def test_docker_image_reason_no_daemon(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_INSPECT_ERR", NO_DAEMON)
    ref = PINNED["bakta"]
    assert _fn("docker_image_reason")(ref) == no_daemon(ref, NO_DAEMON)


def test_wrong_local_digest_is_absent(fake: dict[str, Path]) -> None:
    """Edge case: an older bakta pulled by tag does not satisfy the pin."""
    ref = PINNED["bakta"]
    set_tags(fake, [ref.split("@")[0] + ":v1.9.4", ref.split("@")[0] + "@sha256:" + "f" * 64])
    assert reasons(fake, [TOOLS["bakta"]]) == {"bakta": absent(ref)}


# ---------------------------------------------------------------- example 15 and the docker rows


@pytest.mark.parametrize("profile", DOCKER_ROWS)
def test_docker_rows_absent_pinned_image_is_missing(fake: dict[str, Path], profile: str) -> None:
    """AC 3 (example 15 for ``docker,conda``): every profile of the docker row checks images."""
    set_tags(fake, [])
    tools = [TOOLS[n] for n in sorted(PINNED)]
    got = reasons(fake, tools, profile)
    assert got == {n: absent(r) for n, r in PINNED.items()}, got
    assert_never_pulls(fake)


def test_example_15_docker_conda_conda_only_rows_unchanged(fake: dict[str, Path]) -> None:
    set_tags(fake, [])
    tools = default_tools()
    got = reasons(fake, tools, "docker,conda")
    for name, ref in pinned_of(tools).items():
        assert got.get(name) == absent(ref), got
    for t in tools:
        if t.name not in PINNED and not t.in_process:
            assert got.get(t.name) == f"prebuilt env {t.env_name}: $GIPPA_ROOT is not set", got


def test_present_in_every_docker_row(fake: dict[str, Path]) -> None:
    set_tags(fake, ALL_REFS)
    tools = [TOOLS[n] for n in sorted(PINNED)]
    for profile in DOCKER_ROWS:
        assert reasons(fake, tools, profile) == {}, profile
    assert sorted(set(inspected(fake))) == ALL_REFS


# ---------------------------------------------------------------- example 16: not docker


def test_example_16_apptainer_is_unchanged_and_never_asks_docker(fake: dict[str, Path]) -> None:
    """Regression guard (AC 3, Q12): apptainer keeps its cache rule; ``apptainer,docker`` too."""
    tools = [TOOLS[n] for n in sorted(PINNED)]
    for profile in ("apptainer", "farmshare,apptainer", "slurm,apptainer", "apptainer,docker"):
        got = reasons(fake, tools, profile)
        assert set(got) == set(PINNED), (profile, got)
        for why in got.values():
            assert "setup --images --profile apptainer" in why and SETUP not in why, why
    assert calls(fake) == []


@pytest.mark.parametrize("profile", ["test,docker", "test,docker,stub", "stub,docker", "test,stub"])
def test_example_16_test_and_stub_need_nothing(fake: dict[str, Path], profile: str) -> None:
    """Regression guard: ``test``/``stub`` need no tools, so no docker call at all."""
    assert reasons(fake, default_tools(), profile) == {}
    assert calls(fake) == []


# ---------------------------------------------------------------- example 17: podman-docker


def test_example_17_podman_with_all_refs(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "podman version 4.9.3")
    set_tags(fake, ALL_REFS)
    tools = [TOOLS[n] for n in sorted(PINNED)]
    assert reasons(fake, tools) == {}
    assert sorted(set(inspected(fake))) == ALL_REFS


def test_podman_missing_image_is_rule_5(
    fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Podman has no daemon: ``image not known`` falls to rule 5."""
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "podman version 4.9.3")
    set_tags(fake, [])
    assert reasons(fake, [TOOLS["bakta"]]) == {"bakta": absent(PINNED["bakta"])}


# ---------------------------------------------------------------- probing cost


def test_no_docker_call_without_a_docker_image(fake: dict[str, Path]) -> None:
    """Regression guard: an all-conda ``docker,conda`` selection and an empty one ask nothing."""
    conda_only = [t for t in default_tools() if t.name not in PINNED and not t.in_process]
    assert conda_only
    reasons(fake, conda_only, "docker,conda")
    reasons(fake, [], "docker")
    assert calls(fake) == []


def test_blast_and_tblastn_share_one_inspect(fake: dict[str, Path]) -> None:
    set_tags(fake, ALL_REFS)
    _fn("missing_tool_reasons")([TOOLS["blast"]], "docker", fake["cache"], None, None)
    assert inspected(fake) == [PINNED["blast"]]
    assert versions(fake) == 1


def test_version_probed_once_for_many_images(fake: dict[str, Path]) -> None:
    set_tags(fake, [])
    tools = [TOOLS[n] for n in sorted(PINNED)]
    _fn("missing_tool_reasons")(tools, "docker", fake["cache"], None, None)
    assert versions(fake) == 1, calls(fake)
    assert set(inspected(fake)) == set(ALL_REFS)


# ---------------------------------------------------------------- fgsea and local UniOP (shared)


def test_fgsea_absent_reason_is_rule_5_with_the_images_json_ref(fake: dict[str, Path]) -> None:
    """AC 5: the fgsea row's absent text is rule 5's, with the ``images.json`` ref."""
    set_tags(fake, [])
    assert reasons(fake, [TOOLS["fgsea"]]) == {"fgsea": absent(FGSEA_REF)}


def test_fgsea_without_docker_is_rule_3(
    tmp_path: Path, fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    empty = tmp_path / "nobin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    assert reasons(fake, [TOOLS["fgsea"]]) == {"fgsea": no_docker(FGSEA_REF)}


def test_fgsea_old_docker_is_rule_1(fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "Docker version 19.03.15, build 99e3ed8919")
    set_tags(fake, [FGSEA_REF])
    assert reasons(fake, [TOOLS["fgsea"]]) == {"fgsea": too_old("docker", "19.3", "20.10")}


def _local_tools(db: Path, tag: str) -> Path:
    db.mkdir(parents=True, exist_ok=True)
    path = db / "local_tools.config"
    path.write_text(
        "profiles {\n"
        f"    docker {{ process {{ withName: 'UNIOP' {{ container = '{tag}'; conda = '' }} }} }}\n"
        "}\n",
        encoding="utf-8",
    )
    return path


def test_local_uniop_docker_build_old_docker_is_rule_1(
    tmp_path: Path, fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The local UniOP task now runs with ``--pull=never`` too: rule 1 applies when the build is
    otherwise usable."""
    monkeypatch.setenv("FAKE_DOCKER_VERSION", "Docker version 19.03.15, build 99e3ed8919")
    set_tags(fake, [LOCAL_TAG])
    local = _local_tools(tmp_path / "db", LOCAL_TAG)
    got = reasons(fake, [TOOLS["uniop"]], "docker", local)
    assert got == {"uniop": too_old("docker", "19.3", "20.10")}, got


def test_local_uniop_docker_build_new_docker_unchanged(
    tmp_path: Path, fake: dict[str, Path]
) -> None:
    """Regression guard: a supported Docker and a present local tag still pass."""
    set_tags(fake, [LOCAL_TAG])
    local = _local_tools(tmp_path / "db", LOCAL_TAG)
    assert reasons(fake, [TOOLS["uniop"]], "docker", local) == {}
    assert_never_pulls(fake)


# ---------------------------------------------------------------- AC 2: every reason names setup


@pytest.mark.parametrize("case", ["absent", "no-docker", "no-daemon"])
def test_ac2_every_non_rule_1_reason_names_setup_and_the_full_ref(
    tmp_path: Path, fake: dict[str, Path], monkeypatch: pytest.MonkeyPatch, case: str
) -> None:
    set_tags(fake, [])
    if case == "no-docker":
        empty = tmp_path / "nobin"
        empty.mkdir()
        monkeypatch.setenv("PATH", str(empty))
    elif case == "no-daemon":
        monkeypatch.setenv("FAKE_DOCKER_INSPECT_ERR", NO_DAEMON)
    tools = [TOOLS[n] for n in sorted(PINNED)]
    got = reasons(fake, tools)
    assert set(got) == set(PINNED), got
    for name, why in got.items():
        assert SETUP in why and PINNED[name] in why, why
        assert "@sha256:" + PINNED[name].split("@sha256:")[1] in why
