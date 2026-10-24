"""``scripts/slurm/benchmark_standard.sbatch`` (docs/specs/M10-12-benchmark-standard.md).

Acceptance criteria 1-5 and examples 1, 2, 5-9: the static header and text checks, the dry run,
the preconditions, a stubbed full run and the validity rules. Nothing is submitted and nothing
reaches the network: the script runs in the sandbox of ``tests/helpers/bench_m10_12.py`` (fake
``srun``/``sacct``/``squeue``/``scontrol``/``giap``/``timeout``/``sleep`` first on ``PATH``, a
throw-away checkout with a local ``origin``, a temporary ``GIAP_ROOT``). The conventions the
spec leaves open are listed in that module's docstring.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import re
import subprocess
import sys
from datetime import date, timedelta
from pathlib import Path
from typing import Any

import pytest
from helpers.bench_m10_12 import (
    ARRAY_JOB,
    COTENANCY_KEYS,
    JOB_BASE,
    K12_GENOME,
    K12_GENOME_REL,
    K12_GENOME_SHA256,
    K12_JSON,
    K12_JSON_REL,
    K12_JSON_SHA256,
    NODE,
    ROTATION,
    RUN_KEYS,
    SACCT_KEYS,
    SCRIPT,
    SEED,
    TRACE_KEYS,
    Bench,
    cell_of,
    commit_files,
    flag,
    giap_argv,
    gz_variant,
    make_bench,
    planned_lines,
    run_dir_of,
    sets,
)
from helpers.slurm_sandbox import REPO, assert_honours_giap_repo

BENCH_ID = "2026-10-30_1f56790"


def _text() -> str:
    assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} does not exist (M10-12)"
    return SCRIPT.read_text()


def _code() -> str:
    """The script without comment lines (``#SBATCH`` directives kept)."""
    keep = [ln for ln in _text().splitlines()
            if not ln.lstrip().startswith("#") or ln.startswith("#SBATCH")]  # fmt: skip
    return "\n".join(keep)


def _sbatch(key: str, text: str | None = None) -> str | None:
    m = re.search(rf"^#SBATCH\s+--{key}(?:[= ](\S+))?\s*$", text or _text(), re.M)
    return (m.group(1) or "") if m else None


def mem_mb(text: str) -> float:
    """Slurm ``--mem`` value → MB (Slurm's M is MiB; no suffix means M)."""
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([KMGT]?)", text.upper())
    assert m, text
    return float(m.group(1)) * {"K": 1 / 1024, "M": 1, "": 1, "G": 1024, "T": 1024**2}[m.group(2)]


def seconds(text: str) -> float:
    """``timeout`` duration (``5h``, ``30m``, ``90``) or Slurm ``--time`` → seconds."""
    m = re.fullmatch(r"(\d+(?:\.\d+)?)([smhd]?)", text)
    if m:
        return float(m.group(1)) * {"": 1, "s": 1, "m": 60, "h": 3600, "d": 86400}[m.group(2)]
    days, _, rest = text.rpartition("-")
    parts = [int(p) for p in rest.split(":")]
    h, mnt, s = {1: (0, parts[0], 0), 2: (0, *parts), 3: tuple(parts)}[len(parts)]
    return (int(days) if days else 0) * 86400 + h * 3600 + mnt * 60 + s


def mem_per_cpu_ok(header: str) -> bool:
    """Example 1: ``--mem / --cpus-per-task`` within ``MaxMemPerCPU`` (4,000 MB on normal)."""
    mem, cpus = _sbatch("mem", header), _sbatch("cpus-per-task", header)
    assert mem and cpus, header
    return mem_mb(mem) / int(cpus) <= 4000


# ================================================================== criterion 1: static checks


def test_script_parses() -> None:
    _text()
    proc = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


def test_bash32_compatible_and_not_exempt() -> None:
    assert "bash32: exempt" not in _text(), "bash 3.2-compatible like the other scripts"
    proc = subprocess.run([sys.executable, str(REPO / "scripts" / "ci" / "check_bash32.py"),
                           str(SCRIPT)], capture_output=True, text=True, cwd=REPO)  # fmt: skip
    assert proc.returncode == 0, proc.stdout + proc.stderr


def test_honours_giap_repo() -> None:
    assert_honours_giap_repo(_text(), "benchmark_standard.sbatch")
    assert re.search(r"^\s*(source|\.)\s+\S*scripts/farmshare_env\.sh", _code(), re.M)


def test_header() -> None:
    assert _sbatch("job-name") == "giap-bench"
    assert _sbatch("partition") == "normal"
    assert _sbatch("qos") == "normal"
    assert _sbatch("cpus-per-task") == "16"
    assert _sbatch("mem") == "64000M", "cluster fact 1 / L1: --mem=64000M, not 64G"
    assert _sbatch("constraint") == "CPU_SKU:6148", "one pinned node class (wheat)"
    assert _sbatch("array") == "0-2%1", "3 repeats, one at a time (L3)"
    assert _sbatch("exclusive") is None, "not exclusive (L4)"
    assert _sbatch("gres") is None
    out = _sbatch("output")
    assert out is not None and re.fullmatch(r"\S*/logs/%x\.%A_%a\.out", out), out


def test_header_mem_per_cpu_within_max() -> None:
    assert mem_per_cpu_ok(_text())


@pytest.mark.parametrize(
    ("mem", "ok"), [("64G", False), ("65536M", False), ("64000M", True), ("60G", True)]
)
def test_example1_mem_rule(mem: str, ok: bool) -> None:
    """Example 1: 65,536 / 16 = 4,096 MB > MaxMemPerCPU, so Slurm would allocate 18 CPUs."""
    header = re.sub(r"^#SBATCH\s+--mem=\S+", f"#SBATCH --mem={mem}", _text(), flags=re.M)
    assert _sbatch("mem", header) == mem
    assert mem_per_cpu_ok(header) is ok


def test_time_covers_four_runs_at_the_default_timeout() -> None:
    m = re.search(r"BENCH_RUN_TIMEOUT:?-([0-9]+[smhd]?)\}", _code())
    assert m and m.group(1) == "5h", "BENCH_RUN_TIMEOUT defaults to 5h"
    wall = _sbatch("time")
    assert wall is not None and seconds(wall) >= 4 * seconds(m.group(1)), wall


def test_resume_occurs_nowhere() -> None:
    text = _text()
    assert "--resume" not in text and not re.search(r"(?<![\w-])-resume\b", text)


def test_profile_default_is_farmshare_local() -> None:
    assert re.search(r"BENCH_PROFILE:?-farmshare_local,apptainer\}", _code()), _code()


def test_run_is_wrapped_in_timeout_with_term_and_kill_after() -> None:
    code = _code()
    assert re.search(r"\btimeout\s+--signal=TERM\s+--kill-after=5m\b", code), code


def test_pins_the_committed_k12_inputs() -> None:
    """Convention: the sha256 of the committed genome and JSON are pinned under scripts/."""
    texts = [_text()] + [p.read_text(errors="ignore") for p in (REPO / "scripts").rglob("*")
                         if p.is_file() and p.suffix in {".sh", ".sbatch", ".txt", ".sha256",
                                                         ".tsv", ".json", ""}]  # fmt: skip
    for sha in (K12_GENOME_SHA256, K12_JSON_SHA256):
        assert any(sha in t for t in texts), f"sha256 {sha} pinned nowhere under scripts/"


def test_never_syncs_the_shared_venv_or_pushes() -> None:
    code = _code()
    assert not re.search(r"\buv\s+sync\b|\bpip\s+install\b|\bgit\s+push\b|\bgit\s+commit\b", code)
    for m in re.finditer(r"\buv\s+run\b[^\n]*", code):
        assert "--no-sync" in m.group(0), m.group(0)


# ================================================================== the sandbox


@pytest.fixture
def bench(tmp_path: Path) -> Bench:
    _text()
    return make_bench(tmp_path)


def _out(res: subprocess.CompletedProcess[str]) -> str:
    return f"exit {res.returncode}\n--- stdout\n{res.stdout}\n--- stderr\n{res.stderr}"


# ================================================================== criterion 2: dry run


@pytest.mark.parametrize("task", [0, 1, 2])
def test_dry_run_order_and_flags(bench: Bench, task: int) -> None:
    """Example 2 (task 1) and the rotation table (tasks 0 and 2)."""
    res = bench.run(task, BENCH_DRY_RUN="1", BENCH_ID=BENCH_ID)
    assert res.returncode == 0, _out(res)
    lines = planned_lines(res.stdout)
    cells = [cell_of(giap_argv(t)) for t in lines]
    assert cells == ROTATION[task], _out(res)
    root = bench.bench_dir / BENCH_ID / f"r{task + 1}"
    db = str(bench.root / "db")
    for toks, cell in zip(lines, cells, strict=True):
        assert cell is not None
        genome, variant = cell.split("/")
        srun = toks[: toks.index("timeout")] if "timeout" in toks else toks
        for opt in ("--ntasks=1", "--cpus-per-task=16", "--mem=64000M", "--export=ALL"):
            assert opt in srun, (opt, toks)
        tmo = toks[toks.index("timeout") :]
        assert "--signal=TERM" in tmo and "--kill-after=5m" in tmo and "5h" in tmo, toks
        g = giap_argv(toks)
        assert run_dir_of(g) == root / f"{genome}_{variant}", g
        assert flag(g, "--threads") == "16" and flag(g, "--memory") == "56.GB", g
        assert flag(g, "--profile") == "farmshare_local,apptainer", g
        assert flag(g, "--seed") == str(SEED), g
        assert flag(g, "--db-dir") == db, g
        assert not any(a in ("--resume", "-resume") for a in toks), toks
        bakta_set = [s for s in sets(g) if s.startswith("annotation.bakta.db")]
        if variant == "standard":
            assert bakta_set == ["annotation.bakta.db=full"], g
            assert flag(g, "--annotation") is None and flag(g, "--modules") is None, g
        else:
            assert bakta_set == [], g
            assert flag(g, "--modules") == "core", g
            ann = flag(g, "--annotation")
            assert ann is not None and ann.endswith(".json") and not ann.endswith(".gz"), g
            if genome == "k12":
                assert ann.endswith("/inputs/k12_mg1655.json"), ann
            else:
                assert ann.endswith("bakta_like.json"), ann
        pep = flag(g, "--peptides")
        assert pep is not None, g
        if genome == "k12":
            assert pep == str(bench.root / "mvp04" / "20261023" / "peptides.tsv"), pep
            assert flag(g, "--config") == str(bench.root / "mvp04" / "20261023" / "config.yaml")
            assert Path(flag(g, "--genome") or "").name.startswith("k12_mg1655"), g
        else:
            assert pep.endswith("peptides.tsv"), pep
            assert (flag(g, "--genome") or "").endswith("genome.fasta"), g


def test_dry_run_runs_and_creates_nothing(bench: Bench) -> None:
    res = bench.run(0, BENCH_DRY_RUN="1")
    assert res.returncode == 0, _out(res)
    assert len(planned_lines(res.stdout)) == 4
    assert not bench.jsonl("srun.jsonl"), "a dry run calls no srun"
    assert not bench.giap_runs(), "a dry run calls no giap run"
    assert not list(bench.bench_dir.glob("*/r*/*")), "a dry run creates no RUN_DIR"
    again = bench.run(0)  # a real run afterwards is not refused as "RUN_DIR exists"
    assert again.returncode == 0, _out(again)


def test_dry_run_default_bench_id(bench: Bench) -> None:
    res = bench.run(0, unset=("BENCH_ID",), BENCH_DRY_RUN="1")
    assert res.returncode == 0, _out(res)
    days = {(date.today() - timedelta(days=d)).isoformat() for d in (0, 1)}
    for toks in planned_lines(res.stdout):
        rd = run_dir_of(giap_argv(toks))
        assert rd is not None
        bench_id = rd.parent.parent.name
        day, _, sha7 = bench_id.partition("_")
        assert day in days and sha7 == bench.commit[:7], bench_id


def test_dry_run_bench_dir_and_subsets(bench: Bench) -> None:
    other = bench.root / "elsewhere"
    res = bench.run(1, BENCH_DRY_RUN="1", BENCH_DIR=str(other), BENCH_GENOMES="k12",
                    BENCH_RUN_TIMEOUT="7m")  # fmt: skip
    assert res.returncode == 0, _out(res)
    lines = planned_lines(res.stdout)
    assert [cell_of(giap_argv(t)) for t in lines] == ["k12/core", "k12/standard"]
    for toks in lines:
        assert str(run_dir_of(giap_argv(toks))).startswith(str(other)), toks
        assert "7m" in toks[toks.index("timeout") :], toks
    res = bench.run(2, BENCH_DRY_RUN="1", BENCH_VARIANTS="core")
    assert [cell_of(giap_argv(t)) for t in planned_lines(res.stdout)] == ["sim5m/core", "k12/core"]


def test_dry_run_replacement_attempt(bench: Bench) -> None:
    res = bench.run(1, BENCH_DRY_RUN="1", BENCH_ID=BENCH_ID, BENCH_ONLY="k12/core",
                    BENCH_ATTEMPT="b")  # fmt: skip
    assert res.returncode == 0, _out(res)
    lines = planned_lines(res.stdout)
    assert len(lines) == 1, _out(res)
    assert run_dir_of(giap_argv(lines[0])) == bench.bench_dir / BENCH_ID / "r2.b" / "k12_core"


# ================================================================== criterion 3: preconditions


def _assert_refused(bench: Bench, res: subprocess.CompletedProcess[str]) -> None:
    assert res.returncode == 2, _out(res)
    assert not bench.jsonl("srun.jsonl"), "no srun before the preconditions pass"
    assert not bench.giap_runs()
    assert not list(bench.root.rglob("bench_run.json"))
    assert (res.stdout + res.stderr).strip(), "a one-line reason"


@pytest.mark.parametrize("profile", ["farmshare,apptainer", "apptainer,farmshare", "farmshare"])
def test_slurm_executor_profile_refused(bench: Bench, profile: str) -> None:
    res = bench.run(0, BENCH_PROFILE=profile)
    _assert_refused(bench, res)
    assert re.search(r"R14|farmshare_local", res.stdout + res.stderr), _out(res)


def test_bench_commit_unset(bench: Bench) -> None:
    _assert_refused(bench, bench.run(0, unset=("BENCH_COMMIT",)))


@pytest.mark.parametrize("kind", ["short", "unknown", "not_hex"])
def test_bench_commit_short_or_unknown(bench: Bench, kind: str) -> None:
    value = {"short": bench.commit[:7], "unknown": "0123456789abcdef0123456789abcdef01234567",
             "not_hex": "z" * 40}[kind]  # fmt: skip
    _assert_refused(bench, bench.run(0, BENCH_COMMIT=value))


def test_existing_run_dir(bench: Bench) -> None:
    rd = bench.bench_dir / BENCH_ID / "r1" / "k12_standard"
    rd.mkdir(parents=True)
    (rd / "keep.txt").write_text("old\n")
    _assert_refused(bench, bench.run(0, BENCH_ID=BENCH_ID))
    assert sorted(p.name for p in rd.iterdir()) == ["keep.txt"]


@pytest.mark.parametrize(
    "env",
    [{"BENCH_GENOMES": ""}, {"BENCH_GENOMES": "k12 ecoli"}, {"BENCH_VARIANTS": ""},
     {"BENCH_VARIANTS": "standard fast"}],
    ids=["genomes_empty", "genome_unknown", "variants_empty", "variant_unknown"],
)  # fmt: skip
def test_empty_or_unknown_subsets(bench: Bench, env: dict[str, str]) -> None:
    _assert_refused(bench, bench.run(0, **env))


@pytest.mark.parametrize("name", ["peptides.tsv", "config.yaml"])
def test_missing_k12_inputs(bench: Bench, name: str) -> None:
    (bench.root / "mvp04" / "20261023" / name).unlink()
    _assert_refused(bench, bench.run(0))


@pytest.mark.parametrize("rel", [K12_GENOME_REL, K12_JSON_REL], ids=["genome", "json"])
def test_k12_sha256_differs_from_the_committed_file(bench: Bench, rel: str) -> None:
    data = (K12_GENOME if rel == K12_GENOME_REL else K12_JSON).read_bytes()
    sha = commit_files(bench.sb, {rel: gz_variant(data)}, "drifted input")
    bench.configure(commit=sha)
    _assert_refused(bench, bench.run(0, BENCH_COMMIT=sha))


def test_failing_setup_check(bench: Bench) -> None:
    bench.configure(setup_rc=1)
    _assert_refused(bench, bench.run(0))
    setup = [r["argv"] for r in bench.jsonl("giap.jsonl") if r["argv"][:1] == ["setup"]]
    assert setup and "--check" in setup[0], setup


# ================================================================== criterion 4: stubbed full run


@pytest.fixture(scope="module")
def full(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Task 0 with all four cells, then the same task into another BENCH_DIR (criterion 4)."""
    _text()
    b = make_bench(tmp_path_factory.mktemp("full"))
    first = b.run(0, BENCH_ID=BENCH_ID)
    runs1 = len(b.giap_runs())
    second = b.run(0, BENCH_ID=BENCH_ID, BENCH_DIR=str(b.root / "bench2"),
                   SLURM_JOB_ID=str(JOB_BASE + 100))  # fmt: skip
    return {"b": b, "first": first, "second": second, "runs1": runs1}


def _records(b: Bench, bench_dir: Path) -> dict[str, dict[str, Any]]:
    out = {}
    for p in sorted((bench_dir / BENCH_ID).glob("r*/*/bench_run.json")):
        out[p.parent.name] = json.loads(p.read_text())
    return out


def test_full_run_exit_and_four_records(full: dict[str, Any]) -> None:
    b, res = full["b"], full["first"]
    assert res.returncode == 0, _out(res)
    recs = _records(b, b.bench_dir)
    assert sorted(recs) == ["k12_core", "k12_standard", "sim5m_core", "sim5m_standard"], recs
    for name, r in recs.items():
        assert r["status"] == "ok", (name, r["status"], r["invalid_reasons"], r["warnings"])


def test_full_run_records_validate(full: dict[str, Any]) -> None:
    b = full["b"]
    for name, r in _records(b, b.bench_dir).items():
        missing = RUN_KEYS - set(r)
        assert not missing, (name, missing)
        assert r["schema"] == "giap-bench-run/1"
        assert set(r["sacct"]) >= SACCT_KEYS, r["sacct"]
        assert set(r["trace"]) >= TRACE_KEYS, r["trace"]
        assert set(r["cotenancy"]) >= COTENANCY_KEYS, r["cotenancy"]
        assert isinstance(r["argv"], list) and isinstance(r["inputs"], dict)
        assert all({"path", "sha256"} <= set(v) for v in r["inputs"].values()), r["inputs"]
        assert isinstance(r["invalid_reasons"], list) and isinstance(r["warnings"], list)


def test_full_run_record_values(full: dict[str, Any]) -> None:
    b = full["b"]
    recs = _records(b, b.bench_dir)
    for i, cell in enumerate(ROTATION[0]):
        genome, variant = cell.split("/")
        r = recs[f"{genome}_{variant}"]
        assert (r["genome"], r["variant"], r["repeat"], r["attempt"]) == (genome, variant, 1, "a")
        assert r["order_index"] in (i, i + 1), "position in the rotation (0- or 1-based)"
        assert r["bench_id"] == BENCH_ID and r["commit"] == b.commit
        assert r["node"] == NODE and r["cpu_sku"] in ("CPU_SKU:6148", "6148"), r["cpu_sku"]
        assert str(r["array_task"]) == "0"
        assert ARRAY_JOB in str(r["job_id"]) or str(JOB_BASE) in str(r["job_id"]), r["job_id"]
        assert r["wall_s"] == 3912, "wall_s is the srun step's sacct Elapsed (01:05:12)"
        s = r["sacct"]
        assert s["elapsed_s"] == 3912 and s["alloc_cpus"] == 16 and s["state"] == "COMPLETED"
        assert s["max_rss_bytes"] == 21_600_000 * 1024, s
        assert s["total_cpu_s"] == pytest.approx(13 * 3600 + 2 * 60 + 24), s
        assert r["exit_code"] == 0 and r["invalid_reasons"] == []
        t = r["trace"]
        assert t["n_cached"] == 0 and t["n_failed"] == 0 and t["n_retried"] == 0
        assert t["n_tasks"] == (5 if variant == "standard" else 3), t
        assert t["hosts"] == [NODE]
        assert t["task_peak_rss_max_bytes"] == pytest.approx(
            (20.6 if variant == "standard" else 3) * 1024**3, rel=1e-6
        )
        assert ("bakta" in r["steps"]) is (variant == "standard"), r["steps"]
        assert r["steps"]["operons"]["realtime_s"] == 180
        assert r["steps"]["operons"]["peak_rss_bytes"] == pytest.approx(2 * 1024**3, rel=1e-6)
    steps = [r["step_id"] for r in recs.values()]
    assert len(set(map(str, steps))) == 4, steps


def test_full_run_record_format(full: dict[str, Any]) -> None:
    b = full["b"]
    for p in (b.bench_dir / BENCH_ID).glob("r*/*/bench_run.json"):
        text = p.read_text()
        obj = json.loads(text)
        assert text in (
            json.dumps(obj, indent=2, sort_keys=True) + "\n",
            json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        ), p


def test_full_run_files_kept_and_work_deleted(full: dict[str, Any]) -> None:
    b = full["b"]
    for rd in (b.bench_dir / BENCH_ID).glob("r1/*"):
        assert (rd / "runtimes.md").is_file(), rd
        assert "| step |" in (rd / "runtimes.md").read_text()
        assert (rd / "cotenancy.tsv").is_file(), rd
        rows = [ln for ln in (rd / "cotenancy.tsv").read_text().splitlines() if ln.strip()]
        assert len(rows) >= 2, "a sample at the start and one at the end"
        assert not (rd / "out" / "work").exists(), "work/ is deleted after the run"
        for keep in ("logs/trace.tsv", "run_manifest.json", "report"):
            assert (rd / "out" / keep).exists(), (rd, keep)


def test_full_run_srun_and_timeout_shape(full: dict[str, Any]) -> None:
    b = full["b"]
    steps = [s for s in b.jsonl("srun.jsonl") if s["cell"]][:4]
    assert [s["cell"] for s in steps] == ROTATION[0]
    for s in steps:
        for opt in ("--ntasks=1", "--cpus-per-task=16", "--mem=64000M", "--export=ALL"):
            assert opt in s["opts"], s
        assert Path(s["cmd"][0]).name == "timeout", s
    tmo = [t for t in b.jsonl("timeout.jsonl") if t["cell"]][:4]
    for t in tmo:
        assert "--signal=TERM" in t["opts"] and "--kill-after=5m" in t["opts"], t
        assert t["duration"] == "5h", t


def test_full_run_giap_argv(full: dict[str, Any]) -> None:
    b = full["b"]
    for r in b.giap_runs()[:4]:
        g = r["argv"]
        assert "--resume" not in g and "-resume" not in g
        assert flag(g, "--profile") == "farmshare_local,apptainer"
        assert flag(g, "--threads") == "16" and flag(g, "--memory") == "56.GB"
        assert r["env"]["SLURM_CPUS_PER_TASK"] == "16", "giap runs inside the srun step"


def test_full_run_k12_core_json_is_decompressed(full: dict[str, Any]) -> None:
    b = full["b"]
    core = [r for r in b.giap_runs()[:4] if r["cell"] == "k12/core"]
    assert len(core) == 1
    ann = flag(core[0]["argv"], "--annotation")
    assert ann and ann.endswith(".json"), ann
    assert Path(ann).read_bytes() == gzip.decompress(K12_JSON.read_bytes())


def test_full_run_inputs_sha256(full: dict[str, Any]) -> None:
    b = full["b"]
    recs = _records(b, b.bench_dir)
    raw = gzip.decompress(K12_GENOME.read_bytes())
    ok = {K12_GENOME_SHA256, hashlib.sha256(raw).hexdigest()}
    for name in ("k12_standard", "k12_core"):
        shas = {v["sha256"] for v in recs[name]["inputs"].values()}
        assert shas & ok, (name, recs[name]["inputs"])
    sims = [{v["sha256"] for v in recs[n]["inputs"].values()} for n in ("sim5m_standard",
                                                                       "sim5m_core")]  # fmt: skip
    assert all(s for s in sims)


def test_full_run_sim5m_from_the_pinned_worktree(full: dict[str, Any]) -> None:
    b = full["b"]
    calls = b.jsonl("giapsim.jsonl")
    assert calls, "SIM5M generated by the worktree's giapsim"
    assert {(c["config"], c["seed"]) for c in calls} == {("SIM5M", SEED)}, calls
    for r in b.giap_runs()[:4]:
        if r["cell"].startswith("sim5m/"):
            gen = Path(flag(r["argv"], "--genome") or "")
            assert gen.is_file() and gen.name == "genome.fasta", r["argv"]


def test_full_run_pinned_worktree_on_pythonpath(full: dict[str, Any]) -> None:
    b = full["b"]
    wt = b.root / "ci" / f"bench-{b.commit[:7]}"
    head = subprocess.run(["git", "-C", str(wt), "rev-parse", "HEAD"], capture_output=True,
                          text=True, env=b.sb.env)  # fmt: skip
    assert head.stdout.strip() == b.commit, head.stderr
    for r in b.giap_runs()[:4]:
        pp = (r["env"]["PYTHONPATH"] or "").split(":")
        assert str(wt / "src") in pp and str(wt / "tests" / "sim") in pp, pp


def test_full_run_fresh_numba_cache_and_launch_dir(full: dict[str, Any]) -> None:
    b = full["b"]
    runs = b.giap_runs()[:4]
    numba = [r["env"]["NUMBA_CACHE_DIR"] for r in runs]
    assert len(set(numba)) == 4, numba
    assert all(r["numba_state"] == "empty" for r in runs), [r["numba_state"] for r in runs]
    assert all(str(n).startswith(str(b.nodetmp)) for n in numba), (numba, b.nodetmp)
    cwds = [r["cwd"] for r in runs]
    assert len(set(cwds)) == 4, "a fresh launch directory per run"
    assert all(".nextflow" not in r["cwd_entries"] for r in runs)
    assert all(not c.startswith(str(b.root / "ci")) for c in cwds), "not the worktree"


def test_full_run_setup_check(full: dict[str, Any]) -> None:
    b = full["b"]
    setup = [r["argv"] for r in b.jsonl("giap.jsonl") if r["argv"][:1] == ["setup"]]
    assert setup, "giap setup --check runs first"
    assert "--check" in setup[0] and flag(setup[0], "--db-dir") == str(b.root / "db"), setup


def test_cotenancy_looks_only_at_our_own_jobs(full: dict[str, Any]) -> None:
    """R14: ``squeue -u $USER``; another user's same-node giap-* job never invalidates."""
    b = full["b"]
    calls = b.jsonl("squeue.jsonl")
    assert calls, "the sampler queries squeue"
    for argv in calls:
        assert any(a in ("-u", "--me") or a.startswith(("--user=", "-u")) for a in argv), argv
    for argv in b.jsonl("scontrol.jsonl"):
        if argv[:2] in (["show", "node"], ["show", "nodes"]):
            assert argv[2:3] == [NODE], argv
        elif argv[:2] == ["show", "job"]:
            assert argv[2:3] and argv[2].split("_")[0] in (ARRAY_JOB, str(JOB_BASE)), argv
        else:
            raise AssertionError(f"unexpected scontrol call {argv}")
    for r in _records(b, b.bench_dir).values():
        assert r["cotenancy"]["giap_other"] == [], "the benchmark's own job is excluded"
        assert r["cotenancy"]["samples"] >= 2


def test_second_task_identical_but_times_and_ids(full: dict[str, Any]) -> None:
    b, res = full["b"], full["second"]
    assert res.returncode == 0, _out(res)
    one, two = _records(b, b.bench_dir), _records(b, b.root / "bench2")
    assert sorted(one) == sorted(two)

    def norm(r: dict[str, Any], d: Path) -> str:
        r = {k: v for k, v in r.items()
             if k not in ("job_id", "step_id", "script_wall_s", "warnings")}  # fmt: skip
        r["cotenancy"] = {k: v for k, v in r["cotenancy"].items() if k != "samples"}
        text = json.dumps(r, sort_keys=True).replace(str(d), "<BENCH_DIR>")
        return re.sub(str(JOB_BASE + 100) + "|" + str(JOB_BASE), "<JOB>", text)

    for name in one:
        assert norm(one[name], b.bench_dir) == norm(two[name], b.root / "bench2"), name


# ================================================================== criterion 5: validity rules


def _k12(bench: Bench, task: int = 0, **env: str) -> subprocess.CompletedProcess[str]:
    return bench.run(task, BENCH_ID=BENCH_ID, BENCH_GENOMES="k12", **env)


def _pair(bench: Bench) -> tuple[dict[str, Any], dict[str, Any]]:
    return bench.record("r1/k12_standard", BENCH_ID), bench.record("r1/k12_core", BENCH_ID)


INVALID = [
    ("cached_tasks", {"giap": {"k12/core": {"cached": 3}}}, "cached_tasks:3"),
    ("manifest_resume", {"giap": {"k12/core": {"resume": True}}}, "manifest_resume"),
    ("commit_other", {"giap": {"k12/core": {"commit": "f" * 40}}}, "commit_mismatch"),
    ("commit_dirty", {"giap": {"k12/core": {"dirty": True}}}, "commit_mismatch"),
    ("host_mismatch", {"giap": {"k12/core": {"host": "wheat-07"}}}, "host_mismatch"),
    ("bakta_ran_in_core", {"giap": {"k12/core": {"bakta_in_core": True}}}, "bakta_ran_in_core"),
    ("image_created", {"giap": {"k12/core": {"image": True}}}, "image_or_env_created"),
    ("env_created", {"giap": {"k12/core": {"conda": True}}}, "image_or_env_created"),
    ("alloc_cpus", {"sacct": {"k12/core": {"alloc_cpus": 18}}}, "alloc_cpus:18"),
]


@pytest.mark.parametrize(("cfg", "reason"), [c[1:] for c in INVALID], ids=[c[0] for c in INVALID])
def test_invalid_run_and_the_other_runs_still_run(
    bench: Bench, cfg: dict[str, Any], reason: str
) -> None:
    """Examples 6 and 8 and criterion 5: one stub per reason; the task exits 1."""
    bench.configure(**cfg)
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    std, core = _pair(bench)
    assert core["status"] == "invalid", core
    assert reason in core["invalid_reasons"], core["invalid_reasons"]
    assert std["status"] == "ok", (std["status"], std["invalid_reasons"])


def test_cached_step_in_manifest_is_invalid(bench: Bench) -> None:
    bench.configure(giap={"k12/core": {"cached_step": True}})
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    core = _pair(bench)[1]
    assert core["status"] == "invalid", core
    assert any(r.startswith("cached") for r in core["invalid_reasons"]), core["invalid_reasons"]


def test_example7a_same_node_giap_job(bench: Bench) -> None:
    bench.configure(squeue=["777|giap-run-pilot-a1|2|" + NODE])
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    for r in _pair(bench):
        assert r["status"] == "invalid" and "same_node_giap:777" in r["invalid_reasons"], r
        others = r["cotenancy"]["giap_other"]
        assert [str(o["job"]) for o in others] == ["777"], others
        assert {"job", "name", "cpus", "node"} <= set(others[0])


def test_example7b_contended_other_nodes_is_valid(bench: Bench) -> None:
    bench.configure(squeue=["778|giap-diann-a2|16|wheat-07"])
    res = _k12(bench)
    assert res.returncode == 0, _out(res)
    for r in _pair(bench):
        assert r["status"] == "ok" and r["invalid_reasons"] == [], r
        assert "contended_other_nodes" in json.dumps(r["warnings"]), r["warnings"]
        o = r["cotenancy"]["giap_other"]
        assert len(o) == 1 and str(o[0]["job"]) == "778" and int(o[0]["cpus"]) == 16, o
        assert o[0]["node"] == "wheat-07" and o[0]["name"] == "giap-diann-a2", o


def test_small_giap_job_on_another_node_adds_no_flag(bench: Bench) -> None:
    bench.configure(squeue=["780|giap-tiny|2|wheat-07", "779|jupyter|4|" + NODE])
    res = _k12(bench)
    assert res.returncode == 0, _out(res)
    for r in _pair(bench):
        assert r["status"] == "ok", r
        assert "contended_other_nodes" not in json.dumps(r["warnings"]), r["warnings"]
        assert [str(o["job"]) for o in r["cotenancy"]["giap_other"]] == ["780"], "giap-* only"


def test_example7c_node_load_reported_never_invalidates(bench: Bench) -> None:
    bench.configure(scontrol={"cpualloc": [64, 78, 70]})
    res = _k12(bench)
    assert res.returncode == 0, _out(res)
    for r in _pair(bench):
        assert r["status"] == "ok" and r["invalid_reasons"] == [], r
        c = r["cotenancy"]
        assert c["node_cpualloc_max"] in (64, 70, 78), c
        assert 64 <= c["node_cpualloc_median"] <= c["node_cpualloc_max"], c


def test_node_class_outside_the_pinned_sku_is_invalid(bench: Bench) -> None:
    bench.configure(scontrol={"features": "CPU_GEN:GNR,CPU_SKU:9384X,CPU_FRQ:3.10GHz"})
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    for r in _pair(bench):
        assert r["status"] == "invalid", r
        assert any(x.startswith("node_class:") and "9384X" in x for x in r["invalid_reasons"])


@pytest.mark.parametrize(
    ("cfg", "status"),
    [
        ({"giap": {"k12/core": {"failed_in": "OPERONS", "cached": 3}}}, "failed"),
        ({"giap": {"k12/core": {"cached": 3}}, "timeout": {"k12/core": {"rc": 124, "run": True}}},
         "timeout"),
        ({"giap": {"k12/core": {"cached": 3, "resume": True}}}, "invalid"),
    ],
    ids=["failed_over_invalid", "timeout_over_invalid", "two_invalid_reasons"],
)  # fmt: skip
def test_status_precedence(bench: Bench, cfg: dict[str, Any], status: str) -> None:
    bench.configure(**cfg)
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    core = _pair(bench)[1]
    assert core["status"] == status, core
    if status == "invalid":
        assert {"cached_tasks:3", "manifest_resume"} <= set(core["invalid_reasons"]), core


def test_nextflow_retry_is_valid(bench: Bench) -> None:
    bench.configure(giap={"k12/core": {"retried": True}})
    res = _k12(bench)
    assert res.returncode == 0, _out(res)
    core = _pair(bench)[1]
    assert core["status"] == "ok" and core["trace"]["n_retried"] == 1, core["trace"]


def test_example5_failed_run_is_kept(bench: Bench) -> None:
    bench.configure(giap={"k12/core": {"failed_in": "OPERONS"}},
                    sacct={"k12/core": {"state": "FAILED", "exit": "4:0"}})  # fmt: skip
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    std, core = _pair(bench)
    assert core["status"] == "failed" and core["exit_code"] == 4, core
    assert core["trace"]["n_failed"] == 1 and core["sacct"]["state"] == "FAILED", core
    assert std["status"] == "ok"


def test_failure_before_nextflow_still_writes_a_record(bench: Bench) -> None:
    bench.configure(giap={"k12/standard": {"early_fail": True}})
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    std, core = _pair(bench)
    assert std["status"] == "failed" and std["exit_code"] == 1, std
    assert std["trace"]["n_tasks"] == 0 and std["steps"] == {}, std
    assert core["status"] == "ok", "the next cell still runs"


def test_out_of_memory_is_failed_with_its_state(bench: Bench) -> None:
    bench.configure(giap={"k12/core": {"rc": 1}},
                    sacct={"k12/core": {"state": "OUT_OF_MEMORY", "exit": "0:125"}})  # fmt: skip
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    core = _pair(bench)[1]
    assert core["status"] == "failed" and core["sacct"]["state"] == "OUT_OF_MEMORY", core
    assert len([r for r in bench.giap_runs() if r["cell"] == "k12/core"]) == 1, "never retried"


@pytest.mark.parametrize("rc", [124, 137])
def test_timeout(bench: Bench, rc: int) -> None:
    bench.configure(timeout={"k12/standard": {"rc": rc}})
    res = _k12(bench)
    assert res.returncode == 1, _out(res)
    std, core = _pair(bench)
    assert std["status"] == "timeout" and std["wall_s"] is not None, std
    assert core["status"] == "ok", "the remaining runs still run"


def test_run_timeout_override_reaches_timeout(bench: Bench) -> None:
    res = _k12(bench, BENCH_RUN_TIMEOUT="7m")
    assert res.returncode == 0, _out(res)
    assert {t["duration"] for t in bench.jsonl("timeout.jsonl") if t["cell"]} == {"7m"}


def test_sigterm_marks_the_current_run_timeout(bench: Bench) -> None:
    bench.configure(giap={"k12/standard": {"term": True}})
    res = _k12(bench)
    assert res.returncode == 143, _out(res)
    std = bench.record("r1/k12_standard", BENCH_ID)
    assert std["status"] == "timeout", std
    root = bench.bench_dir / BENCH_ID / "r1"
    assert not (root / "k12_core" / "bench_run.json").exists(), "not-started runs are absent"


def test_empty_max_rss_is_null(bench: Bench) -> None:
    bench.configure(sacct={"default": {"max_rss": ""}})
    res = _k12(bench)
    assert res.returncode == 0, _out(res)
    for r in _pair(bench):
        assert r["sacct"]["max_rss_bytes"] is None, r["sacct"]
        assert r["trace"]["task_peak_rss_max_bytes"] > 0


def test_sacct_unreachable_falls_back_to_script_wall(bench: Bench) -> None:
    bench.configure(sacct={"rc": 1})
    res = _k12(bench, BENCH_ONLY="k12/core")
    core = bench.record("r1/k12_core", BENCH_ID)
    assert len(bench.jsonl("sacct.jsonl")) >= 3, "retried 3 times"
    assert core["wall_s"] == core["script_wall_s"] and core["wall_s"] is not None, core
    assert "wall_from_script" in json.dumps(core["warnings"]), (core["warnings"], _out(res))


def test_wall_cross_check_warning(bench: Bench) -> None:
    """sacct Elapsed and the script's own clock must agree within 5 s, else a warning."""
    bench.configure(sacct={"k12/standard": {"elapsed_s": 0, "total_cpu": "00:00:00"}})
    res = _k12(bench)
    assert res.returncode == 0, _out(res)
    std, core = _pair(bench)
    assert std["wall_s"] == 0 and std["warnings"] == [], std["warnings"]
    assert core["wall_s"] == 3912 and core["warnings"] != [], "3,912 s vs a sub-second run"


def test_replacement_attempt_writes_its_own_dir(bench: Bench) -> None:
    res = bench.run(1, BENCH_ID=BENCH_ID, BENCH_ONLY="k12/core", BENCH_ATTEMPT="b")
    assert res.returncode == 0, _out(res)
    r = bench.record("r2.b/k12_core", BENCH_ID)
    assert (r["repeat"], r["attempt"], r["genome"], r["variant"]) == (2, "b", "k12", "core")
    assert [x["cell"] for x in bench.giap_runs()] == ["k12/core"]
    assert not (bench.bench_dir / BENCH_ID / "r2").exists()


def test_existing_attempt_dir_is_refused(bench: Bench) -> None:
    (bench.bench_dir / BENCH_ID / "r2.b" / "k12_core").mkdir(parents=True)
    _assert_refused(bench, bench.run(1, BENCH_ID=BENCH_ID, BENCH_ONLY="k12/core",
                                     BENCH_ATTEMPT="b"))  # fmt: skip
