"""M2-15b: ``plan`` / ``chunk`` / ``merge`` of ``workflow/bin/gippa_terminatornet_shard.py``
against a fake upstream (docs/specs/M2-15b-terminatornet-sharding.md items 2-4; examples 4, 5, 6,
7, 8, 9, 10).

Interface pinned here (test author, from the spec): every subcommand works in the cwd with the
staged names ``genome.fa``, ``plan.tsv``, ``model.pkl``, ``calls_<chunk_id>.pkl`` and
``terminatornet.tsv``; ``plan [--core N]`` (default CORE), ``chunk <chunk_id> --cpus N`` and
``merge``. Upstream ``TerminatorNet`` / ``Features`` are imported from ``$TERMINATORNET_HOME``. A
calls pickle (protocol 4) is a dict ``{chunk_id, core_start, core_end, calls}``.

The fake upstream (written into tmp_path, no real tool, no network) has upstream's names and
signatures: ``read_in_sequences`` drops records with no sequence, ``clean_sequences`` blanks
sequences under 4 nt, upper-cases and drops non-ACGT; ``model.pkl`` holds ``(MODEL, SCALER)``.
Its window scorer reads the module global ``Features.SEQUENCE`` (also past the chunk's core and
its length) and logs every evaluated window ``j`` to ``$FAKE_TNET_LOG``; its ``merge_overlapping``
logs the calls it receives and is a greedy sweep in which ``pred >= best_pred`` replaces the best.
"""

from __future__ import annotations

import csv
import json
import os
import pickle
import subprocess
import sys
from operator import itemgetter
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "bin" / "gippa_terminatornet_shard.py"

FEATURES_PY = """
import os
import numpy as np
SEQUENCE, GC_PERCENT = "", 0.5
HAIRPIN_LENGTH, TAIL_LENGTH, WINDOW, STEP = 52, 8, 60, 1
FEATURES = ["f_pred", "f_len"]

def reverse_complement(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def _log(j):
    p = os.environ.get("FAKE_TNET_LOG")
    if p:
        with open(p, "a") as fh:
            fh.write("%d %d\\n" % (len(SEQUENCE), j))

def scan_oneprocess_1(j):
    _log(j)
    assert GC_PERCENT == 0.5
    if SEQUENCE[j:j + 1] == "G":
        return (SEQUENCE[j:j + HAIRPIN_LENGTH], j, 40.0)
    return "", -1, -1

def fake_pred(seq, j):
    far = seq[min(j + 55, len(seq) - 1)]
    return 0.2 if seq[j + 1:j + 2] == "T" else 0.5 + ((7 * j + len(seq) + ord(far)) % 40) / 100

def scan_oneprocess_2(candidate):
    hrpn, j, tt = candidate
    feats = {"f_pred": fake_pred(SEQUENCE, j), "f_len": float(len(SEQUENCE)), "hairpin_start": 3,
             "hairpin_end": 30, "structure": ".((((....)))).", "GH": -12.5}
    return (feats, hrpn, "TTTTTTTT", "AAAAAAAA", tt, j)

class FakeScaler:
    def transform(self, x):
        return np.asarray(x, dtype=float)

class FakeModel:
    def predict_proba(self, x):
        p = float(x[0, 0])
        return np.array([[1.0 - p, p]])
"""

TERMINATORNET_PY = """
import json, os, sys
from Features import *
OUT = sys.stdout
BOTH_STRANDS = False

def read_in_sequences(filename):
    seqs, names, s, name = [], [], "", None
    for line in open(filename):
        if line.startswith(">"):
            if s:
                seqs.append(s); names.append(name)
            s, name = "", line[1:].strip()
        else:
            s += line.strip()
    if s:
        seqs.append(s); names.append(name)
    return seqs, names

def clean_sequences(sequences):
    long_seq = False
    for i, s in enumerate(sequences):
        if len(s) < 4:
            sequences[i] = ""
            continue
        s = s.upper().replace("U", "T")
        sequences[i] = "".join(c for c in s if c in "ACGT")
        long_seq = long_seq or len(sequences[i]) >= 501
    return sequences, 0.5, long_seq

def merge_overlapping(terms, seq_length=0, strand="+"):
    with open(os.environ["FAKE_TNET_MERGE_LOG"], "a") as fh:
        fh.write(json.dumps({"both": BOTH_STRANDS, "L": seq_length, "strand": strand,
                             "calls": [[t[0], t[1], float(t[7])] for t in terms]}) + "\\n")
    out, i = [], 0
    while i < len(terms):
        best, j = terms[i], i + 1
        while j < len(terms) and terms[j][0] <= best[1]:
            if terms[j][7] >= best[7]:
                best = terms[j]
            j += 1
        out.append([best[0] + 1, best[1] + 1, best[7], strand])
        i = j
    return out

def output_terms(terms, name):
    if not terms:
        OUT.write(name + "\\tNO TERMINATORS FOUND\\n")
    for i, t in enumerate(terms):
        OUT.write(name + "\\tTERM_" + str(i + 1) + "\\t" + "\\t".join(map(str, t)) + "\\n")
"""


@pytest.fixture
def home(tmp_path: Path) -> Path:
    h = tmp_path / "tnet_home"
    h.mkdir()
    (h / "Features.py").write_text(FEATURES_PY)
    (h / "TerminatorNet.py").write_text(TERMINATORNET_PY)
    return h


def _env(home: Path, **extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "TERMINATORNET_HOME")}
    env.update(TERMINATORNET_HOME=str(home), PYTHONNOUSERSITE="1", **extra)
    return env


def _run(home: Path, cwd: Path, *args: str, **extra: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=cwd,
        env=_env(home, **extra),
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )


def _workdir(home: Path, tmp: Path, fasta: str, name: str = "wd") -> Path:
    wd = tmp / name
    wd.mkdir()
    (wd / "genome.fa").write_text(fasta)
    code = (
        "import pickle, Features; "
        "pickle.dump((Features.FakeModel(), Features.FakeScaler()), open('model.pkl', 'wb'))"
    )
    subprocess.run(
        [sys.executable, "-c", code], cwd=wd, env=_env(home, PYTHONPATH=str(home)), check=True
    )
    return wd


def _plan(wd: Path) -> list[dict[str, str]]:
    with (wd / "plan.tsv").open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _seq(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=n))


def _fasta(records: list[tuple[str, str]]) -> str:
    return "".join(
        f">{n}\n" + "".join(s[i : i + 70] + "\n" for i in range(0, len(s), 70)) for n, s in records
    )


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _reference(seq: str) -> list[tuple]:
    """The unsharded fake scan: every window of the whole strand sequence, in window order."""
    out = []
    for j in range(max(len(seq) - 59, 1)):
        if seq[j : j + 1] != "G":
            continue
        far = seq[min(j + 55, len(seq) - 1)]
        p = 0.2 if seq[j + 1 : j + 2] == "T" else 0.5 + ((7 * j + len(seq) + ord(far)) % 40) / 100
        if p >= 0.5:
            out.append(
                (j + 3, j + 30, "AAAAAAAA", seq[j : j + 52], "TTTTTTTT", ".((((....)))).", -12.5, p)
            )
    return out


# ---------------------------------------------------------------- plan


def test_plan_tsv_examples_4_5_8_empty_record_two_nt_and_ns(home: Path, tmp_path: Path) -> None:
    a = _seq(700, 1)
    fasta = f">a\n{a}\n>e\n>t\nAC\n>c\nACGTNNNNacgt{'gattaca' * 10}\n"
    wd = _workdir(home, tmp_path, fasta)
    r = _run(home, wd, "plan", "--core", "250")
    assert r.returncode == 0, r.stderr
    header = (wd / "plan.tsv").read_text().splitlines()[0]
    assert header == "chunk_id\tseq_index\tseq_name\tstrand\tseq_length\tcore_start\tcore_end\tmode"
    rows = _plan(wd)
    got = [
        (x["chunk_id"], x["seq_name"], x["strand"], x["seq_length"], x["core_start"], x["core_end"])
        for x in rows
    ]
    a_cores = [("0", "250"), ("250", "500"), ("500", "641")]
    want = [(f"c00000p{k:04d}", "a", "+", "700", *c) for k, c in enumerate(a_cores)]
    want += [(f"c00000m{k:04d}", "a", "-", "700", *c) for k, c in enumerate(a_cores)]
    want += [("c00001p0000", "t", "+", "0", "0", "1"), ("c00001m0000", "t", "-", "0", "0", "1")]
    want += [("c00002p0000", "c", "+", "78", "0", "19"), ("c00002m0000", "c", "-", "78", "0", "19")]
    assert got == want
    assert {x["mode"] for x in rows} == {"long"}


def test_plan_example_10_all_records_short(home: Path, tmp_path: Path) -> None:
    wd = _workdir(home, tmp_path, _fasta([("x", _seq(500, 2)), ("y", _seq(40, 3))]))
    assert _run(home, wd, "plan").returncode == 0
    assert {x["mode"] for x in _plan(wd)} == {"short"}


def test_plan_no_sequence_exits_1_with_upstream_message(home: Path, tmp_path: Path) -> None:
    wd = _workdir(home, tmp_path, ">e1\n>e2\n")
    r = _run(home, wd, "plan")
    assert r.returncode == 1
    assert "could not read in genomic sequences" in r.stderr


# ---------------------------------------------------------------- chunk


def _chunk(home: Path, wd: Path, cid: str, cpus: int = 2) -> tuple[dict, list[int]]:
    log = wd / f"windows_{cid}_{cpus}.log"
    r = _run(home, wd, "chunk", cid, "--cpus", str(cpus), FAKE_TNET_LOG=str(log))
    assert r.returncode == 0, r.stderr
    raw = (wd / f"calls_{cid}.pkl").read_bytes()
    assert raw[:2] == b"\x80\x04", "pickle protocol 4"
    js = [int(line.split()[1]) for line in log.read_text().splitlines()] if log.exists() else []
    return pickle.loads(raw), js


@pytest.mark.parametrize(("length", "core"), [(2_000, 1_000), (1_500, 333)])
def test_chunks_own_their_windows_and_read_the_whole_strand(
    home: Path, tmp_path: Path, length: int, core: int
) -> None:
    s = _seq(length, length)
    wd = _workdir(home, tmp_path, _fasta([("chr", s)]))
    assert _run(home, wd, "plan", "--core", str(core)).returncode == 0
    for strand, seq in (("+", s), ("-", _rc(s))):
        rows = [x for x in _plan(wd) if x["strand"] == strand]
        concat: list[tuple] = []
        for x in rows:
            payload, js = _chunk(home, wd, x["chunk_id"])
            lo, hi = int(x["core_start"]), int(x["core_end"])
            assert sorted(js) == list(range(lo, hi)), (
                f"{x['chunk_id']} must scan exactly its own windows, once"
            )
            assert (payload["chunk_id"], payload["core_start"], payload["core_end"]) == (
                x["chunk_id"],
                lo,
                hi,
            )
            assert all(lo <= c[0] - 3 < hi for c in payload["calls"])
            assert all(type(c[7]) is np.float64 for c in payload["calls"]), "upstream pred type"
            concat += [tuple(c) for c in payload["calls"]]
        assert concat == _reference(seq), (
            f"{strand}: sharded calls differ from the single-chunk scan"
        )


def test_example_9_chunk_bytes_do_not_depend_on_cpus(home: Path, tmp_path: Path) -> None:
    # Pool results come back as fresh objects whose sharing depends on the pool's batching; pickle's
    # memo is identity-based, so equal strings must be pickled the same way at any CPU count.
    wd = _workdir(home, tmp_path, _fasta([("chr", _seq(1_200, 9))]))
    assert _run(home, wd, "plan", "--core", "400").returncode == 0
    blobs = set()
    for cpus in (1, 3):
        _chunk(home, wd, "c00000m0001", cpus)
        blobs.add((wd / "calls_c00000m0001.pkl").read_bytes())
    assert len(blobs) == 1


def test_example_5_chunk_of_a_cleaned_empty_record(home: Path, tmp_path: Path) -> None:
    wd = _workdir(home, tmp_path, ">t\nAC\n")
    assert _run(home, wd, "plan").returncode == 0
    payload, js = _chunk(home, wd, "c00000p0000")
    assert (
        js == [0]
        and payload["calls"] == []
        and (payload["core_start"], payload["core_end"]) == (0, 1)
    )


# ---------------------------------------------------------------- merge


def _call(start: int, stop: int, pred: float) -> tuple:
    return (start, stop, "AAAAAAAA", "GCGC", "TTTTTTTT", ".((..)).", -10.0, np.float64(pred))


def _merge_dir(
    home: Path, tmp: Path, records: list[tuple[str, str]], core: int, calls: dict[str, list]
) -> Path:
    wd = _workdir(home, tmp, _fasta(records))
    assert _run(home, wd, "plan", "--core", str(core)).returncode == 0
    for x in _plan(wd):
        payload = {
            "chunk_id": x["chunk_id"],
            "core_start": int(x["core_start"]),
            "core_end": int(x["core_end"]),
            "calls": calls.get(x["chunk_id"], []),
        }
        (wd / f"calls_{x['chunk_id']}.pkl").write_bytes(pickle.dumps(payload, protocol=4))
    return wd


def _merge(home: Path, wd: Path) -> tuple[subprocess.CompletedProcess[str], list[dict]]:
    log = wd / "merge.log"
    r = _run(home, wd, "merge", FAKE_TNET_MERGE_LOG=str(log))
    return r, [json.loads(line) for line in log.read_text().splitlines()] if log.exists() else []


@pytest.mark.parametrize(("pred_b", "winner"), [(0.93, (1001, 1031)), (0.92, (981, 1011))])
def test_example_7_tie_across_a_chunk_edge(
    home: Path, tmp_path: Path, pred_b: float, winner: tuple
) -> None:
    calls = {"c00000p0000": [_call(980, 1010, 0.93)], "c00000p0001": [_call(1000, 1030, pred_b)]}
    wd = _merge_dir(home, tmp_path, [("chr", _seq(2_000, 7))], 1_000, calls)
    r, log = _merge(home, wd)
    assert r.returncode == 0, r.stderr
    assert log[0] == {
        "both": True,
        "L": 0,
        "strand": "+",
        "calls": [[980, 1010, 0.93], [1000, 1030, pred_b]],
    }
    lines = (wd / "terminatornet.tsv").read_text().splitlines()
    assert lines[0].split("\t")[:4] == ["chr", "TERM_1", str(winner[0]), str(winner[1])]
    assert lines[1] == "chr\tNO TERMINATORS FOUND" and len(lines) == 2


@pytest.mark.parametrize("long_mode", [True, False])
def test_merge_block_order_and_mode_dependent_stop_sort(
    home: Path, tmp_path: Path, long_mode: bool
) -> None:
    n = 900 if long_mode else 400
    c0 = [_call(10, 80, 0.7), _call(12, 40, 0.8), _call(14, 80, 0.6)]
    c1 = [_call(150, 170, 0.9), _call(152, 160, 0.55)]
    calls = {"c00000m0000": c0, "c00000m0001": c1, "c00001p0000": [_call(5, 30, 0.6)]}
    wd = _merge_dir(home, tmp_path, [("a", _seq(n, 4)), ("b", _seq(300, 5))], 100, calls)
    r, log = _merge(home, wd)
    assert r.returncode == 0, r.stderr
    assert [(e["strand"], e["L"]) for e in log] == [("+", 0), ("+", 0), ("-", n), ("-", 300)]
    assert all(e["both"] for e in log), "TerminatorNet.BOTH_STRANDS must be True"
    flat = [[c[0], c[1], float(c[7])] for c in c0 + c1]
    want = sorted(flat, key=itemgetter(1)) if long_mode else flat
    assert log[2]["calls"] == want
    names = [line.split("\t")[0] for line in (wd / "terminatornet.tsv").read_text().splitlines()]
    assert names[0] == "a" and names[1] == "b" and names[2] == "a" and names[-1] == "b"


def test_merge_example_5_two_nt_record_no_terminators_on_both_strands(
    home: Path, tmp_path: Path
) -> None:
    wd = _merge_dir(home, tmp_path, [("t", "AC")], 250_000, {})
    r, _ = _merge(home, wd)
    assert r.returncode == 0, r.stderr
    assert (
        wd / "terminatornet.tsv"
    ).read_text() == "t\tNO TERMINATORS FOUND\nt\tNO TERMINATORS FOUND\n"


@pytest.mark.parametrize("fault", ["missing", "extra", "core_end", "chunk_id"])
def test_merge_fails_on_missing_extra_or_mismatched_calls(
    home: Path, tmp_path: Path, fault: str
) -> None:
    wd = _merge_dir(home, tmp_path, [("chr", _seq(2_000, 8))], 1_000, {})
    target = wd / "calls_c00000m0001.pkl"
    payload = pickle.loads(target.read_bytes())
    if fault == "missing":
        target.unlink()
    elif fault == "extra":
        extra = {"chunk_id": "c00001p0000", "core_start": 0, "core_end": 1, "calls": []}
        (wd / "calls_c00001p0000.pkl").write_bytes(pickle.dumps(extra, protocol=4))
    else:
        payload[fault] = 1_000 if fault == "core_end" else "c00000m0000"
        target.write_bytes(pickle.dumps(payload, protocol=4))
    r, _ = _merge(home, wd)
    assert r.returncode == 1, f"{fault}: exit {r.returncode}\n{r.stderr}"
    assert not (wd / "terminatornet.tsv").exists(), (
        "nothing may be written before validation passes"
    )
