"""M11-15a: ``detect_profile`` (executor D1-D4, runtime R0-R5, combination, hints), host probe.

Spec: docs/specs/M11-15-portable-clusters.md item 1 ("Profile auto-detection (15a)"), the executor
and runtime tables, examples 1.1-1.13 (1.14 is ``gippa setup``: test_m11_15a_cli_auto.py), edge
cases (a), (b), the 10 s timeouts and "ties ... resolve by the fixed rule order, never by PATH
order" (cross-cutting edge cases); acceptance criteria 1 and 5 (the pinned constant).

Interface pinned here (readings chosen where the spec leaves room):

- New module ``gippa.runner.profiles``.
- ``Probe`` (frozen dataclass, every field defaulted, ``Probe()`` = "nothing found"):
  ``on_path: frozenset[str]`` (the names ``shutil.which`` found among ``sbatch``, ``apptainer``,
  ``singularity``, ``docker``, ``conda``, ``mamba``, ``micromamba``), ``docker_ok: bool``
  (``docker info`` exited 0 within 10 s), ``slurm_job_id: str | None``, ``cluster_name: str |
  None``, ``setup_runtime: str | None`` (``setup.runtime`` of ``<db_dir>/db_manifest.json``).
- ``detect_profile(probe, sites=()) -> Detected``: pure (no PATH lookup, no subprocess).
  ``Detected.profile`` is the comma-joined profile, executor first (``"slurm,apptainer"``), or
  ``None`` when no runtime is found (R5). ``Detected.reason`` is the text printed after
  ``auto:`` (contains ``FarmShare detected`` / ``site <s> detected`` / ``sbatch on PATH`` for D1-D3
  and ``as recorded by gippa setup`` when R0 chose the runtime). ``Detected.hints`` holds the
  hints of edge cases (a) and (b); the R0 warning of example 1.10 is in ``Detected.warnings``
  (or, accepted too, in ``hints``).
- ``Site(name, cluster_name, runtime)`` and ``read_sites(sites_dir) -> tuple[Site, ...]`` parse
  the ``// site:``, ``// runtime:`` and ``// detect: cluster_name=`` headers of
  ``<sites_dir>/*.config`` (``_``-prefixed files are not sites), for D2 (fixture directory only;
  the first real site lands with 15j).
- ``FARMSHARE_CLUSTER == "farmshare"`` (D1's pinned constant).
- ``probe_host(environ, db_dir) -> Probe``: the real probe, through ``shutil.which`` and
  ``subprocess.run`` (looked up as module attributes, so the tests patch them), each subprocess
  with ``timeout=10``; a timeout or ``OSError`` counts as not found. ``SLURM_CLUSTER_NAME`` wins
  over ``scontrol show config``'s ``ClusterName``.
- R0 on a Slurm host (``sbatch`` on PATH) never selects docker (item 1 "Combination"): a recorded
  ``docker`` is passed over for R1/R2/R4 there; with docker the only runtime, example 1.12 applies.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

HINT_A = (
    "sbatch on PATH but only docker found: tasks run on this host; install Apptainer or pass "
    "--profile"
)
HINT_B = "the Nextflow head process runs on this host for the whole run; submit gippa run as a job"
NAMES = ("sbatch", "apptainer", "singularity", "docker", "conda", "mamba", "micromamba")
EXECUTORS = {"farmshare", "slurm"}
RUNTIMES = {"apptainer", "singularity", "docker", "conda"}


def _mod() -> Any:
    try:
        from gippa.runner import profiles
    except ImportError as exc:
        pytest.fail(f"gippa.runner.profiles is not implemented (M11-15a): {exc}")
    return profiles


def _probe(*names: str, **kw: Any) -> Any:
    kw.setdefault("docker_ok", "docker" in names)
    return _mod().Probe(on_path=frozenset(names), **kw)


def _detect(probe: Any, sites: Any = ()) -> Any:
    return _mod().detect_profile(probe, sites=sites)


def _messages(d: Any) -> list[str]:
    return [*d.hints, *getattr(d, "warnings", ())]


@pytest.fixture(autouse=True)
def _no_host(monkeypatch: pytest.MonkeyPatch) -> None:
    """``detect_profile`` is pure: any PATH lookup or subprocess fails the test."""

    def boom(*_a: Any, **_k: Any) -> Any:
        raise AssertionError("detect_profile touched the host")

    monkeypatch.setattr(shutil, "which", boom)
    monkeypatch.setattr(subprocess, "run", boom)
    monkeypatch.setattr(subprocess, "Popen", boom)


# ---------------------------------------------------------------- examples 1.1-1.13


EXAMPLES: list[tuple[str, dict[str, Any], str | None]] = [
    ("1.1 laptop docker", {"names": ("docker",), "docker_ok": True}, "docker"),
    ("1.2 workstation apptainer", {"names": ("apptainer",)}, "apptainer"),
    ("1.3 cluster", {"names": ("sbatch", "apptainer"), "slurm_job_id": "42"}, "slurm,apptainer"),
    ("1.4 singularityce", {"names": ("sbatch", "singularity"), "slurm_job_id": "42"},
     "slurm,singularity"),
    ("1.5 farmshare", {"names": ("sbatch", "apptainer"), "cluster_name": "farmshare"},
     "farmshare,apptainer"),
    ("1.6 docker+apptainer, daemon up", {"names": ("docker", "apptainer"), "docker_ok": True},
     "docker"),
    ("1.6 docker+apptainer, daemon down", {"names": ("docker", "apptainer"), "docker_ok": False},
     "apptainer"),
    ("1.7 slurm conda", {"names": ("sbatch", "micromamba"), "slurm_job_id": "42"}, "slurm,conda"),
    ("1.8 nothing", {"names": ()}, None),
    ("1.9 setup docker", {"names": ("docker", "apptainer"), "docker_ok": True,
                          "setup_runtime": "docker"}, "docker"),
    ("1.10 setup apptainer gone", {"names": ("docker",), "docker_ok": True,
                                   "setup_runtime": "apptainer"}, "docker"),
    ("1.11 farmshare setup apptainer", {"names": ("sbatch", "apptainer"),
                                        "cluster_name": "farmshare",
                                        "setup_runtime": "apptainer"}, "farmshare,apptainer"),
    ("1.12 slurm docker only", {"names": ("sbatch", "docker"), "docker_ok": True,
                                "slurm_job_id": "42"}, "docker"),
    ("1.13 login node", {"names": ("sbatch", "apptainer")}, "slurm,apptainer"),
]  # fmt: skip


@pytest.mark.parametrize(("label", "spec", "want"), EXAMPLES, ids=[e[0] for e in EXAMPLES])
def test_examples(label: str, spec: dict[str, Any], want: str | None) -> None:
    spec = dict(spec)
    names = spec.pop("names")
    assert _detect(_probe(*names, **spec)).profile == want, label


def test_1_5_farmshare_is_todays_default_and_ignores_other_runtimes() -> None:
    for names in (("sbatch", "apptainer"), ("sbatch", "docker", "singularity", "micromamba")):
        d = _detect(_probe(*names, cluster_name="farmshare", docker_ok=True))
        assert d.profile == "farmshare,apptainer"
        assert "FarmShare detected" in d.reason


def test_1_9_r0_reason_names_setup() -> None:
    d = _detect(_probe("docker", "apptainer", docker_ok=True, setup_runtime="docker"))
    assert "as recorded by gippa setup" in d.reason


def test_r0_beats_the_host_order() -> None:
    """R0 first: a recorded apptainer wins over R3 docker on a host without Slurm."""
    d = _detect(_probe("docker", "apptainer", docker_ok=True, setup_runtime="apptainer"))
    assert d.profile == "apptainer"
    assert "as recorded by gippa setup" in d.reason


def test_r0_beats_r1_on_a_slurm_host() -> None:
    d = _detect(_probe("sbatch", "apptainer", "singularity", slurm_job_id="1",
                       setup_runtime="singularity"))  # fmt: skip
    assert d.profile == "slurm,singularity"
    assert "sbatch on PATH" in d.reason and "as recorded by gippa setup" in d.reason


def test_1_10_unusable_recorded_runtime_warns() -> None:
    d = _detect(_probe("docker", docker_ok=True, setup_runtime="apptainer"))
    assert d.profile == "docker"
    assert any("setup used apptainer, not available here" in m for m in _messages(d))
    assert "as recorded by gippa setup" not in d.reason


def test_r0_docker_needs_docker_info() -> None:
    """R0 'still usable here': for docker, ``docker info`` exits 0."""
    d = _detect(_probe("docker", "apptainer", docker_ok=False, setup_runtime="docker"))
    assert d.profile == "apptainer"
    assert any("setup used docker, not available here" in m for m in _messages(d))


def test_r0_never_selects_docker_for_a_slurm_executor() -> None:
    d = _detect(_probe("sbatch", "docker", "apptainer", docker_ok=True, slurm_job_id="1",
                       setup_runtime="docker"))  # fmt: skip
    assert d.profile == "slurm,apptainer"


def test_1_11_farmshare_never_local_apptainer() -> None:
    d = _detect(_probe("sbatch", "apptainer", cluster_name="farmshare", setup_runtime="docker",
                       docker_ok=True))  # fmt: skip
    assert d.profile == "farmshare,apptainer"


def test_1_12_hint_a() -> None:
    d = _detect(_probe("sbatch", "docker", docker_ok=True, slurm_job_id="1"))
    assert d.profile == "docker"
    assert any(HINT_A in h for h in d.hints)


def test_1_13_hint_b_on_a_login_node() -> None:
    d = _detect(_probe("sbatch", "apptainer"))
    assert d.profile == "slurm,apptainer"
    assert "sbatch on PATH" in d.reason
    assert any(HINT_B in h for h in d.hints)


def test_1_3_no_hint_inside_an_allocation() -> None:
    d = _detect(_probe("sbatch", "apptainer", slurm_job_id="42"))
    assert not any(HINT_B in h for h in d.hints)
    assert not any(HINT_A in h for h in d.hints)


def test_no_hint_on_a_host_without_slurm() -> None:
    assert _detect(_probe("apptainer")).hints == ()


# ---------------------------------------------------------------- rules beyond the examples


@pytest.mark.parametrize("conda", ["conda", "mamba", "micromamba"])
def test_r4_any_conda_frontend(conda: str) -> None:
    assert _detect(_probe(conda)).profile == "conda"
    assert _detect(_probe("sbatch", conda, slurm_job_id="1")).profile == "slurm,conda"


def test_ties_resolve_by_rule_order() -> None:
    assert _detect(_probe("singularity", "apptainer")).profile == "apptainer"
    assert _detect(_probe("micromamba", "singularity")).profile == "singularity"
    assert _detect(_probe("micromamba", "docker", docker_ok=True)).profile == "docker"
    every = _probe("sbatch", "micromamba", "singularity", "apptainer", "docker", docker_ok=True,
                   slurm_job_id="1")  # fmt: skip
    assert _detect(every).profile == "slurm,apptainer"


def test_docker_without_daemon_is_not_a_runtime() -> None:
    assert _detect(_probe("docker", docker_ok=False)).profile is None


def test_sbatch_without_any_runtime_is_r5() -> None:
    assert _detect(_probe("sbatch", slurm_job_id="1")).profile is None


def test_slurm_job_id_alone_is_not_slurm() -> None:
    """D3 needs ``sbatch`` on PATH (the launching host submits)."""
    assert _detect(_probe("apptainer", slurm_job_id="1")).profile == "apptainer"


def test_other_cluster_name_is_not_farmshare() -> None:
    d = _detect(_probe("sbatch", "apptainer", cluster_name="sherlock", slurm_job_id="1"))
    assert d.profile == "slurm,apptainer"


def test_pinned_cluster_constant() -> None:
    """AC 5: D1 compares against the pinned constant ``farmshare``."""
    assert _mod().FARMSHARE_CLUSTER == "farmshare"


# ---------------------------------------------------------------- D2: a fixture site directory


def _site_dir(tmp_path: Path) -> Path:
    d = tmp_path / "sites"
    d.mkdir()
    (d / "acme.config").write_text(
        "// site: acme\n// scheduler: slurm\n// runtime: singularity\n"
        "// maintainer: @acme-hpc\n// detect: cluster_name=acme\n// tested: no\n"
        "process { executor = 'slurm' }\n",
        encoding="utf-8",
    )
    (d / "nodetect.config").write_text(
        "// site: nodetect\n// scheduler: slurm\n// runtime: apptainer\n"
        "// maintainer: @someone\n// tested: no\nprocess { executor = 'slurm' }\n",
        encoding="utf-8",
    )
    (d / "_template.config").write_text(
        "// site: template\n// runtime: apptainer\n// detect: cluster_name=acme\n",
        encoding="utf-8",
    )
    return d


def test_read_sites_parses_headers(tmp_path: Path) -> None:
    sites = {s.name: s for s in _mod().read_sites(_site_dir(tmp_path))}
    assert set(sites) == {"acme", "nodetect"}, "files starting with _ are not sites"
    assert (sites["acme"].cluster_name, sites["acme"].runtime) == ("acme", "singularity")
    assert sites["nodetect"].cluster_name is None


def test_read_sites_missing_dir_is_empty(tmp_path: Path) -> None:
    assert tuple(_mod().read_sites(tmp_path / "absent")) == ()


def test_d2_site_detected(tmp_path: Path) -> None:
    """Example 10.6: ``detect: cluster_name=acme`` and the cluster name ``acme``."""
    sites = _mod().read_sites(_site_dir(tmp_path))
    d = _detect(_probe("sbatch", "apptainer", "docker", docker_ok=True, cluster_name="acme"),
                sites)  # fmt: skip
    assert d.profile == "acme,singularity"
    assert "site acme detected" in d.reason


def test_d2_ignores_recorded_runtime(tmp_path: Path) -> None:
    """The runtime step is skipped for D1/D2: the site's ``// runtime:`` header decides."""
    sites = _mod().read_sites(_site_dir(tmp_path))
    d = _detect(_probe("sbatch", "apptainer", cluster_name="acme", setup_runtime="apptainer"),
                sites)  # fmt: skip
    assert d.profile == "acme,singularity"


def test_d2_other_cluster_falls_through(tmp_path: Path) -> None:
    sites = _mod().read_sites(_site_dir(tmp_path))
    d = _detect(_probe("sbatch", "apptainer", cluster_name="other", slurm_job_id="1"), sites)
    assert d.profile == "slurm,apptainer"


def test_d1_before_d2(tmp_path: Path) -> None:
    d = tmp_path / "sites"
    d.mkdir()
    (d / "fs2.config").write_text(
        "// site: fs2\n// runtime: singularity\n// detect: cluster_name=farmshare\n",
        encoding="utf-8",
    )
    got = _detect(_probe("sbatch", "apptainer", cluster_name="farmshare"), _mod().read_sites(d))
    assert got.profile == "farmshare,apptainer"


# ---------------------------------------------------------------- properties


@st.composite
def probes(draw: st.DrawFn) -> Any:
    names = draw(st.frozensets(st.sampled_from(NAMES)))
    return _mod().Probe(
        on_path=names,
        docker_ok=draw(st.booleans()),
        slurm_job_id=draw(st.sampled_from([None, "7"])),
        cluster_name=draw(st.sampled_from([None, "farmshare", "sherlock"])),
        setup_runtime=draw(st.sampled_from([None, "apptainer", "singularity", "docker"])),
    )


def _usable(p: Any) -> set[str]:
    out = {n for n in ("apptainer", "singularity") if n in p.on_path}
    if "docker" in p.on_path and p.docker_ok:
        out.add("docker")
    if p.on_path & {"conda", "mamba", "micromamba"}:
        out.add("conda")
    return out


@settings(max_examples=300, deadline=None)
@given(probes())
def test_property_combination_rules(p: Any) -> None:
    d = _detect(p)
    assert d == _detect(p), "deterministic"
    if p.cluster_name == "farmshare":
        assert d.profile == "farmshare,apptainer"
        return
    usable = _usable(p)
    if not usable:
        assert d.profile is None
        return
    assert d.profile is not None
    parts = d.profile.split(",")
    assert len([x for x in parts if x in EXECUTORS]) <= 1, "at most one executor profile"
    assert len([x for x in parts if x in RUNTIMES]) == 1, "exactly one runtime"
    assert set(parts) <= EXECUTORS | RUNTIMES
    runtime = next(x for x in parts if x in RUNTIMES)
    assert runtime in usable, "the chosen runtime is usable on this host"
    if "docker" in parts:
        assert parts == ["docker"], "docker never runs under a Slurm executor"
    slurm_ok = "sbatch" in p.on_path and bool(usable - {"docker"})
    assert ("slurm" in parts) == slurm_ok
    if p.setup_runtime in usable and not ("sbatch" in p.on_path and p.setup_runtime == "docker"):
        assert runtime == p.setup_runtime, "R0 wins when usable"


@settings(max_examples=200, deadline=None)
@given(probes())
def test_property_hints(p: Any) -> None:
    d = _detect(p)
    b = any(HINT_B in h for h in d.hints)
    a = any(HINT_A in h for h in d.hints)
    if p.cluster_name == "farmshare":
        return  # D1: the hints are not specified for FarmShare
    if d.profile is not None and d.profile.startswith("slurm,"):
        assert b == (p.slurm_job_id is None)
    else:
        assert not b
    assert a == (d.profile == "docker" and "sbatch" in p.on_path)


# ---------------------------------------------------------------- probe_host (mocked host)


class _Host:
    """A fake host: ``which`` answers from ``names``; ``run`` answers ``docker info`` and
    ``scontrol show config`` and records every argv with its timeout."""

    def __init__(self, names: set[str], docker: Any = 0, scontrol: Any = "") -> None:
        self.names, self.docker, self.scontrol = names, docker, scontrol
        self.calls: list[tuple[list[str], Any]] = []

    def which(self, name: str, *_a: Any, **_k: Any) -> str | None:
        return f"/fake/bin/{name}" if name in self.names else None

    def run(self, argv: Any, *_a: Any, **kw: Any) -> Any:
        argv = [str(a) for a in argv]
        self.calls.append((argv, kw.get("timeout")))
        prog = Path(argv[0]).name
        answer = self.docker if prog == "docker" else self.scontrol
        if isinstance(answer, BaseException):
            raise answer
        rc, out = (answer, "") if isinstance(answer, int) else (0, answer)
        return subprocess.CompletedProcess(argv, rc, stdout=out, stderr="")


@pytest.fixture
def host(monkeypatch: pytest.MonkeyPatch) -> Any:
    def install(**kw: Any) -> _Host:
        h = _Host(**kw)
        monkeypatch.setattr(shutil, "which", h.which)
        monkeypatch.setattr(subprocess, "run", h.run)
        return h

    return install


def _probe_host(environ: dict[str, str], db_dir: Path | None = None) -> Any:
    return _mod().probe_host(environ, db_dir)


@pytest.mark.real_host_probe
def test_probe_host_path_lookups(host: Any) -> None:
    host(names={"sbatch", "apptainer", "micromamba", "ls"})
    p = _probe_host({})
    assert p.on_path == frozenset({"sbatch", "apptainer", "micromamba"})
    assert p.docker_ok is False
    assert p.slurm_job_id is None and p.setup_runtime is None


@pytest.mark.real_host_probe
def test_probe_host_docker_info_with_timeout(host: Any) -> None:
    h = host(names={"docker"}, docker=0)
    assert _probe_host({}).docker_ok is True
    docker_calls = [(a, t) for a, t in h.calls if Path(a[0]).name == "docker"]
    assert docker_calls and docker_calls[0][0][1:] == ["info"]
    assert docker_calls[0][1] == 10


@pytest.mark.real_host_probe
@pytest.mark.parametrize(
    "answer", [1, subprocess.TimeoutExpired(["docker", "info"], 10), OSError("broken")]
)
def test_probe_host_docker_unreachable(host: Any, answer: Any) -> None:
    host(names={"docker"}, docker=answer)
    p = _probe_host({})
    assert "docker" in p.on_path and p.docker_ok is False


@pytest.mark.real_host_probe
def test_probe_host_no_docker_no_subprocess(host: Any) -> None:
    h = host(names={"apptainer"})
    _probe_host({})
    assert not h.calls


@pytest.mark.real_host_probe
def test_probe_host_slurm_env(host: Any) -> None:
    h = host(names={"sbatch", "scontrol"}, scontrol="ClusterName = other\n")
    p = _probe_host({"SLURM_JOB_ID": "1900311", "SLURM_CLUSTER_NAME": "farmshare"})
    assert (p.slurm_job_id, p.cluster_name) == ("1900311", "farmshare")
    assert not [a for a, _t in h.calls if Path(a[0]).name == "scontrol"], (
        "SLURM_CLUSTER_NAME wins; scontrol is not run"
    )


@pytest.mark.real_host_probe
def test_probe_host_scontrol_cluster_name(host: Any) -> None:
    text = "Configuration data as of 2026-10-27\nAccountingStorageType   = x\n" \
           "ClusterName             = farmshare\nCompleteWait = 0 sec\n"  # fmt: skip
    h = host(names={"sbatch", "scontrol"}, scontrol=text)
    assert _probe_host({}).cluster_name == "farmshare"
    sc = [(a, t) for a, t in h.calls if Path(a[0]).name == "scontrol"]
    assert sc and sc[0][0][1:] == ["show", "config"] and sc[0][1] == 10


@pytest.mark.real_host_probe
@pytest.mark.parametrize("answer", [subprocess.TimeoutExpired(["scontrol"], 10), 1])
def test_probe_host_scontrol_slow_or_failing(host: Any, answer: Any) -> None:
    host(names={"sbatch", "scontrol"}, scontrol=answer)
    assert _probe_host({}).cluster_name is None


@pytest.mark.real_host_probe
def test_probe_host_scontrol_absent(host: Any) -> None:
    h = host(names={"sbatch"})
    assert _probe_host({}).cluster_name is None
    assert not h.calls


def _manifest(db: Path, setup: dict[str, Any] | None) -> Path:
    db.mkdir(parents=True, exist_ok=True)
    doc: dict[str, Any] = {"schema_version": "1.1", "items": []}
    if setup is not None:
        doc["setup"] = setup
    (db / "db_manifest.json").write_text(json.dumps(doc), encoding="utf-8")
    return db


@pytest.mark.real_host_probe
def test_probe_host_reads_setup_runtime(host: Any, tmp_path: Path) -> None:
    host(names=set())
    db = _manifest(tmp_path / "db", {"runtime": "singularity"})
    assert _probe_host({}, db).setup_runtime == "singularity"


@pytest.mark.real_host_probe
@pytest.mark.parametrize("setup", [None, {}])
def test_probe_host_no_setup_runtime(host: Any, tmp_path: Path, setup: Any) -> None:
    host(names=set())
    db = _manifest(tmp_path / "db", setup)
    assert _probe_host({}, db).setup_runtime is None
    assert _probe_host({}, tmp_path / "missing").setup_runtime is None
    assert _probe_host({}, None).setup_runtime is None


@pytest.mark.real_host_probe
def test_probe_host_corrupt_manifest(host: Any, tmp_path: Path) -> None:
    host(names=set())
    db = tmp_path / "db"
    db.mkdir()
    (db / "db_manifest.json").write_text("{not json", encoding="utf-8")
    assert _probe_host({}, db).setup_runtime is None


@pytest.mark.real_host_probe
def test_probe_host_reads_no_real_environment(host: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Only the ``environ`` passed in counts, not ``os.environ``."""
    host(names=set())
    monkeypatch.setenv("SLURM_JOB_ID", "999")
    monkeypatch.setenv("SLURM_CLUSTER_NAME", "farmshare")
    p = _probe_host({})
    assert p.slurm_job_id is None and p.cluster_name is None
