"""Paths, parsers and a container sandbox for the M3-24 pilot-kit tests (GIAP_SPEC.md §15).

Nothing here reaches the network, the real DIA-NN image or the real ``$GIAP_ROOT``:

* the cfg/README/db_a2 parsers only read files tracked in the repository;
* ``pilot_sandbox`` builds on ``helpers.slurm_sandbox.make_sandbox`` (a throw-away checkout of
  ``scripts/``, ``ci/`` and ``validation/`` with fake Slurm tools first on a minimal ``PATH``) and
  adds fake ``apptainer``/``singularity``/``wget``/``curl``, an unreachable HTTP proxy and a
  ``DIANN_SIF`` that points at a 12-byte stand-in, so no container and no download ever runs.

Conventions settled by the test author (lead decisions 1-6 of docs/specs/M3-24-pilot-kit.md):

* A cfg file is an INI-like list of DIA-NN flags with two sections, ``[step1]`` (predict the
  in-silico library) and ``[step2]`` (search the raw files). Section headers are matched
  case-insensitively and may be written ``[step 1]``/``[step-1]``/``[step_1]``. Blank lines and
  ``#`` comment lines are ignored; a line holds one flag and its arguments (several flags on one
  line are parsed too). ``$VAR``/``${VAR}`` in a value is expanded by the sbatch script.
* ``scripts/slurm/diann_pilot.sbatch`` takes the arm (``A1``/``A2``) as its only positional
  argument and reads ``DIANN_SIF``, ``RAW_DIR`` and ``OUT_DIR`` from the environment.
* ``RAW_DIR`` holds the timsTOF runs as ``*.d`` directories (PXD078459 is diaPASEF, §15.1).
* The sandbox creates every ``--fasta`` the cfg files name **except** the ``db_a2`` database:
  only ``strain_db.faa.gz`` is committed, so the script itself must unpack it outside the
  checkout (decision 3: the repository is never written to).
"""

from __future__ import annotations

import inspect
import re
import shlex
import shutil
import string
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path

from helpers.slurm_sandbox import REPO, Sandbox, make_sandbox, write_exe

# ------------------------------------------------------------------ paths and constants

PILOT = REPO / "validation" / "pilot"
README = PILOT / "README.md"
CFG = {"A1": PILOT / "diann_A1.cfg", "A2": PILOT / "diann_A2.cfg"}
DB_A2 = PILOT / "db_a2"
STRAIN_DB = DB_A2 / "strain_db.faa.gz"
MANIFEST = DB_A2 / "makedb_manifest.json"
ID_MAP = DB_A2 / "db_id_map.tsv"  # ``giap makedb`` writes this name; see ID_MAP_NAMES
ID_MAP_NAMES = ("db_id_map.tsv", "id_map.tsv")

SLURM = REPO / "scripts" / "slurm"
SCRIPTS = {
    "diann_pilot": SLURM / "diann_pilot.sbatch",
    "denovo_pilot": SLURM / "denovo_pilot.sbatch",
    "fetch_pxd078459": SLURM / "fetch_pxd078459.sbatch",
}

# Lead decision 1: the image Khoa supplied, the binary inside it and the bind mounts.
DIANN_SIF = "/farmshare/home/classes/bios/270/jordan/diann_2.3.1.sif"
DIANN_BIN = "/diann-2.3.1/diann-linux"
DIANN_VERSION = "DIA-NN 2.3.1 Academia"
BINDS = ("/home/classes", "/scratch", "/farmshare")
SPEC_VERSION = "2.7.0"  # what GIAP_SPEC.md §15 names; decision 2 records the difference

NISSLE_DIR = REPO / "tests" / "data" / "recorded" / "bakta" / "1.12.2-full" / "nissle"
RECORDED_NISSLE = NISSLE_DIR / "nissle1917.json.gz"

# Flags of the §15.2 template that the sbatch script supplies, not the cfg file.
SBATCH_SUPPLIED = frozenset({"--f", "--lib", "--fasta", "--out", "--out-lib", "--threads"})


# ------------------------------------------------------------------ cfg parsing

Flag = tuple[str, tuple[str, ...]]

_SECTION = re.compile(r"^\[\s*step[\s_-]*([12])\s*\]$", re.I)


def read_cfg(path: Path) -> dict[str, list[Flag]]:
    """``{"step1": [(flag, args), …], "step2": [...]}`` of a pilot cfg file."""
    assert path.is_file(), f"{path.relative_to(REPO)} does not exist (M3-24)"
    out: dict[str, list[Flag]] = {}
    section: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        m = _SECTION.match(line)
        if m:
            section = f"step{m.group(1)}"
            out.setdefault(section, [])
            continue
        tokens = shlex.split(line, comments=True)
        if not tokens:
            continue
        assert section is not None, f"{path.name}: flag before any [step1]/[step2] header: {line}"
        out[section].extend(_split_flags(tokens))
    return out


def _split_flags(tokens: list[str]) -> Iterator[Flag]:
    flag: str | None = None
    args: list[str] = []
    for tok in tokens:
        if tok.startswith("--"):
            if flag is not None:
                yield flag, tuple(args)
            flag, args = tok, []
        else:
            assert flag is not None, f"value before any flag: {tokens}"
            args.append(tok)
    if flag is not None:
        yield flag, tuple(args)


def values(flags: Iterable[Flag], name: str) -> list[str]:
    """Every argument given to ``name`` (one entry per occurrence, "" for a bare flag)."""
    return [" ".join(args) for flag, args in flags if flag == name]


def names(flags: Iterable[Flag]) -> set[str]:
    return {flag for flag, _ in flags}


# ------------------------------------------------------------------ GIAP's own digest


def giap_cut_spec(enzyme: str = "trypsin_p") -> str:
    """DIA-NN ``--cut`` specification of ``giap.mapping.digest`` (expected: ``K*,R*``)."""
    from giap.mapping.digest import cut_sites

    residues = [
        aa
        for aa in string.ascii_uppercase
        if len(cut_sites(f"AAA{aa}AAA", enzyme)) == 2  # type: ignore[arg-type]
    ]
    blocked = len(cut_sites("AAAKPAAA", enzyme)) == 1  # type: ignore[arg-type]
    suffix = "" if blocked else "*"
    return ",".join(f"{r}{suffix}" for r in residues)


def giap_digest_flags() -> dict[str, str]:
    """The DIA-NN digest flags implied by ``giap.mapping.digest`` (M3-22, lead decision 5).

    ``--unimod4`` (fixed carbamidomethyl on C) is not a digest parameter: GIAP's masses assume it,
    so it is required as a literal.
    """
    from giap.mapping.digest import digest

    p = inspect.signature(digest).parameters
    return {
        "--cut": giap_cut_spec(),
        "--missed-cleavages": str(p["missed"].default),
        "--min-pep-len": str(p["min_len"].default),
        "--max-pep-len": str(p["max_len"].default),
        "--unimod4": "",
    }


def giap_excises_initiator_met() -> bool:
    """Whether GIAP's digest produces Met-excised peptides (it does not: digest.py decision 1)."""
    from giap.mapping.digest import digest

    return any(start == 1 for start, _, _ in digest("MAAAAAAAAK" + "AAAAAAAK"))


# ------------------------------------------------------------------ the spec's DIA-NN template


def spec_template_flags() -> set[str]:
    """Flags of the ``DIA-NN template`` code block of GIAP_SPEC.md §15.2."""
    text = (REPO / "GIAP_SPEC.md").read_text()
    i = text.index("DIA-NN template")
    block = text[i:].split("```", 2)[1]
    return {t for t in block.split() if t.startswith("--")}


# ------------------------------------------------------------------ the sandbox

_APPTAINER = r"""#!{python}
'''Fake apptainer/singularity: logs the call, never executes anything.

``--version`` and an inner ``--help``/``--version`` are the smoke check lead decision 6 allows.
A real call creates the files DIA-NN would write (``--out``, ``--out-lib`` and the
``.predicted.speclib`` sibling of ``--out-lib``) and fails when ``--lib`` does not exist.
'''
import json, os, sys

state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
rec = {"tool": os.path.basename(sys.argv[0]), "argv": argv, "cwd": os.getcwd(),
       "binds": [], "image": None, "inner": [], "created": [], "smoke": False}

def log(rc):
    rec["rc"] = rc
    with open(os.path.join(state, "diann.jsonl"), "a") as fh:
        fh.write(json.dumps(rec) + "\n")
    sys.exit(rc)

i = 0
while i < len(argv):
    a = argv[i]
    if a in ("-B", "--bind"):
        rec["binds"] += argv[i + 1].split(","); i += 2; continue
    if a.startswith("--bind="):
        rec["binds"] += a.split("=", 1)[1].split(","); i += 1; continue
    if a.endswith((".sif", ".simg", ".img")):
        rec["image"] = a
        rec["inner"] = argv[i + 1:]
        break
    i += 1

if rec["image"] is None:
    rec["smoke"] = True
    print("apptainer version 1.4.4")
    log(0)

if not os.path.exists(rec["image"]):
    print("FATAL: image not found: " + rec["image"], file=sys.stderr)
    log(255)

inner = rec["inner"][1:]
if not inner or set(inner) <= {"--help", "--version", "--dry-run"}:
    rec["smoke"] = True
    print("DIA-NN 2.3.1 Academia (Data-Independent Acquisition by Neural Networks)")
    print("Compiled on Dec  5 2025")
    log(0)

def value(flag):
    return inner[inner.index(flag) + 1] if flag in inner and inner.index(flag) + 1 < len(inner) \
        else None

def touch(path):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("stub DIA-NN output\n")
    rec["created"].append(path)

lib = value("--lib")
if lib is not None and not os.path.exists(lib):
    print("FATAL: library not found: " + lib, file=sys.stderr)
    log(1)

out_lib = value("--out-lib")
if out_lib is not None:
    touch(out_lib)
    touch(os.path.splitext(out_lib)[0] + ".predicted.speclib")
    touch(out_lib + ".predicted.speclib")
out = value("--out")
if out is not None:
    touch(out)
log(int(os.environ.get("FAKE_DIANN_RC", "0")))
"""

_FETCH = r"""#!{python}
'''Fake wget/curl: logs the call and writes a stand-in payload. Nothing reaches the network.'''
import json, os, sys

state = os.environ["FAKE_SLURM_DIR"]
tool = os.path.basename(sys.argv[0])
argv = sys.argv[1:]
urls, target, outdir, insecure, remote_name = [], None, None, False, False
i = 0
while i < len(argv):
    a = argv[i]
    if a in ("-o", "--output", "-O", "--output-document") and i + 1 < len(argv) \
            and not (tool == "curl" and a == "-O"):
        target = argv[i + 1]; i += 2; continue
    if a in ("-P", "--directory-prefix", "--output-dir") and i + 1 < len(argv):
        outdir = argv[i + 1]; i += 2; continue
    if a in ("-O", "--remote-name"):
        remote_name = True; i += 1; continue
    if a in ("-k", "--insecure", "--no-check-certificate"):
        insecure = True; i += 1; continue
    if a.startswith(("http://", "https://", "ftp://")):
        urls.append(a)
    i += 1

created = []
for n, url in enumerate(urls):
    path = target if (target and n == 0) else os.path.join(outdir or os.getcwd(),
                                                           os.path.basename(url) or "index")
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(("stub payload for " + url + "\n").encode())
    created.append(os.path.abspath(path))

rec = {"tool": tool, "argv": argv, "urls": urls, "created": created, "insecure": insecure,
       "remote_name": remote_name, "cwd": os.getcwd()}
with open(os.path.join(state, "fetch.jsonl"), "a") as fh:
    fh.write(json.dumps(rec) + "\n")
rc = int(os.environ.get("FAKE_FETCH_RC", "0"))
if rc:
    print(tool + ": fake failure", file=sys.stderr)
sys.exit(rc)
"""


@dataclass
class PilotSandbox:
    sb: Sandbox
    sif: Path
    raw_dir: Path
    out_dir: Path

    def run(self, script: str, *args: str, **env: str):  # type: ignore[no-untyped-def]
        """Run ``scripts/slurm/<script>.sbatch`` in the sandbox (never through Slurm)."""
        path = self.sb.repo / "scripts" / "slurm" / f"{script}.sbatch"
        return self.sb.run(["bash", str(path), *args], env=env, cwd=self.sb.tmp)

    def diann_calls(self, *, smoke: bool = False) -> list[dict[str, object]]:
        calls = self.sb.jsonl("diann.jsonl")
        return [c for c in calls if bool(c["smoke"]) == smoke]

    def fetch_calls(self) -> list[dict[str, object]]:
        return self.sb.jsonl("fetch.jsonl")

    def repo_is_clean(self) -> str:
        import subprocess

        cmd = ["git", "-C", str(self.sb.repo), "status", "--porcelain"]
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()


def pilot_sandbox(tmp: Path, *, sif: bool = True, raws: int = 3) -> PilotSandbox:
    """A sandbox checkout with fake container and download tools (nothing real runs)."""
    sb = make_sandbox(tmp)
    shutil.copytree(REPO / "validation", sb.repo / "validation",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))  # fmt: skip
    env_sh = sb.repo / "scripts" / "farmshare_env.sh"
    env_sh.write_text(env_sh.read_text() + '\nexport GIAP_TOOLS_DIR="$GIAP_ROOT/tools"\n')
    import subprocess

    subprocess.run(["git", "-C", str(sb.repo), "add", "-A"], check=True, capture_output=True,
                   env=sb.env)  # fmt: skip
    subprocess.run(["git", "-C", str(sb.repo), "commit", "-q", "-m", "validation"], check=True,
                   capture_output=True, env=sb.env)  # fmt: skip

    py = __import__("sys").executable
    for name in ("apptainer", "singularity"):
        write_exe(sb.bin / name, _APPTAINER.replace("{python}", py))
    for name in ("wget", "curl"):
        write_exe(sb.bin / name, _FETCH.replace("{python}", py))

    image = sb.tmp / "images" / "diann_2.3.1.sif"
    image.parent.mkdir(parents=True, exist_ok=True)
    if sif:
        image.write_text("stub sif\n")
    raw_dir = sb.root / "data" / "pxd078459"
    raw_dir.mkdir(parents=True, exist_ok=True)
    for i in range(raws):
        d = raw_dir / f"wt_rep{i + 1}.d"
        d.mkdir(exist_ok=True)
        (d / "analysis.tdf").write_text("stub timsTOF data\n")
    out_dir = sb.root / "pilot"
    sb.env.update({
        "DIANN_SIF": str(image),
        "RAW_DIR": str(raw_dir),
        "OUT_DIR": str(out_dir),
        "SLURM_JOB_ID": "4242",
        "SLURM_CPUS_PER_TASK": "16",
        "http_proxy": "http://127.0.0.1:9",
        "https_proxy": "http://127.0.0.1:9",
        "ftp_proxy": "http://127.0.0.1:9",
        "no_proxy": "",
    })  # fmt: skip
    _materialise_fastas(sb)
    return PilotSandbox(sb=sb, sif=image, raw_dir=raw_dir, out_dir=out_dir)


def _materialise_fastas(sb: Sandbox) -> None:
    """Create the ``--fasta`` inputs the cfg files name, outside the sandbox repository."""
    env = {
        "GIAP_ROOT": str(sb.root),
        "GIAP_DATA_DIR": str(sb.root / "data"),
        "GIAP_DB_DIR": str(sb.root / "db"),
        "GIAP_REPO": str(sb.repo),
        "GIAP_TOOLS_DIR": str(sb.root / "tools"),
    }
    for path in CFG.values():
        if not path.is_file():
            continue
        for section in read_cfg(path).values():
            for raw in values(section, "--fasta"):
                expanded = raw
                for key, val in env.items():
                    expanded = expanded.replace(f"${{{key}}}", val).replace(f"${key}", val)
                # db_a2 is committed gzipped: the script itself must unpack it (decision 3).
                if "db_a2" in expanded or not expanded.startswith("/"):
                    continue
                if str(sb.repo) in expanded:
                    continue
                target = Path(expanded)
                if target.exists():
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(">stub_protein\nMAAAKAAAAAAAKAAAAAAR\n")
